import csv,json,subprocess,sys,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import gemmi,numpy as np
ROOT=Path(__file__).resolve().parents[3];B=Path(__file__).resolve().parents[1];OLD=ROOT/'outputs/docking_audit_2026-10-24'
def tab(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
PLAN=tab(B/'inputs/47_structure_and_search_plan.tsv')
INV={r['gene']:r for r in tab(B/'inputs/47_verified_input_inventory.tsv')}
JAVA=OLD/'software/jdk-17.0.20.1+1-jre/bin/java.exe'
PRANK=OLD/'software/p2rank_2.5.1'
(B/'pockets').mkdir(exist_ok=True)
def predict(p):
    g=p['gene'];dest=B/'pockets'/g;dest.mkdir(exist_ok=True)
    source=B/'prepared/receptors'/f'{g}_raw.pdb'
    predicted=dest/f'{source.name}_predictions.csv'
    if predicted.exists():return
    cmd=[str(JAVA),'-Xmx2048m','-cp',str(PRANK/'bin/p2rank.jar')+';'+str(PRANK/'bin/lib/*'),
         'cz.siret.prank.program.Main','predict','-f',str(source),'-o',str(dest),'-threads','2','-c','alphafold']
    start=time.monotonic()
    with (B/'logs'/f'p2rank_{g}.log').open('w',encoding='utf8') as f:
        proc=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    (B/'audit'/f'{g}_pocket_run.json').write_text(json.dumps({'cmd':cmd,'exit_code':proc.returncode,'elapsed_seconds':time.monotonic()-start},indent=2),encoding='utf8')
    assert proc.returncode==0 and predicted.exists(),g
    print(g,'P2Rank complete',flush=True)

def get_atoms(path):
    st=gemmi.read_structure(str(path));return [(c.name,r.seqid.num,a) for c in st[0] for r in c for a in r if not a.element.is_hydrogen]
def xyz_for(path):return np.array([[a.pos.x,a.pos.y,a.pos.z] for _,_,a in get_atoms(path)])
def ligand_padding():
    # Padding must allow the ligand to extend beyond the target surface in global searches.
    from rdkit import Chem
    mol=next(Chem.ForwardSDMolSupplier((OLD/'prepared/spiropachysine_pH74.sdf').open('rb'),removeHs=False))
    pos=Chem.RemoveHs(mol).GetConformer().GetPositions()
    return max(6.,float(np.linalg.norm(pos-pos.mean(0),axis=1).max())+2.)

def summarize():
    rows=[];errors=[];padding=ligand_padding()
    for p in PLAN:
        g=p['gene'];target=B/'prepared/targets'/f'{g}_target.pdb'
        t_atoms=get_atoms(target);allxyz=xyz_for(target)
        target_residues={(c,n) for c,n,a in t_atoms}
        rule=p['site_rule'];basis={};pocketrank='';pocketprob=''
        try:
            if rule=='global_search':
                points=allxyz;pad=padding;basis={'source':'all_actual_target_heavy_atoms','padding_each_side_A':padding}
            elif rule.startswith('experimental_ligand:'):
                meta=json.loads((B/'audit'/f'{g}_extraction.json').read_text(encoding='utf8'))
                refs=meta['reference_ligands'];assert refs,(g,'missing experimental reference')
                # For multiple copies use the closest recorded copy of each desired CCD to the primary target chain.
                chosen=[]
                for ccd in p['known_ligand_CCD'].split('+'):
                    chosen.append(min((r for r in refs if r['CCD']==ccd),key=lambda r:(r['distance_A'],r['resnum'])))
                points=np.vstack([xyz_for(ROOT/r['path']) for r in chosen]);pad=6.
                basis={'source':'experimental_ligand_extent','references':chosen}
            else:
                pp=B/'pockets'/g/f'{g}_raw.pdb_predictions.csv' if p['source']=='PDB' else ROOT/INV[g]['pocket_predictions_path']
                pockets=[{k.strip():v.strip() for k,v in x.items()} for x in csv.DictReader(pp.open())]
                candidates=[]
                for pocket in pockets:
                    ids=[(r.rsplit('_',1)[0],int(r.rsplit('_',1)[1])) for r in pocket['residue_ids'].split()]
                    overlap=set(ids)&target_residues
                    if len(overlap)>=3 and len(overlap)/max(len(ids),1)>=.25:candidates.append((pocket,overlap))
                assert candidates,(g,'no target-associated pocket in selected actual structure')
                pocket,overlap=candidates[0]
                points=np.array([[a.pos.x,a.pos.y,a.pos.z] for c,n,a in t_atoms if (c,n) in overlap]);pad=6.
                pocketrank=pocket['rank'];pocketprob=pocket['probability']
                basis={'source':str(pp.relative_to(ROOT)),'pocket_rank':pocketrank,'pocket_probability':pocketprob,
                    'target_residues':[f'{c}:{n}' for c,n in sorted(overlap)],'target_association_rule':'at least3 target residues and >=25% of pocket residues'}
            lo=points.min(0);hi=points.max(0);center=(hi+lo)/2;size=np.maximum(hi-lo+2*pad,24.)
            assert np.isfinite(center).all() and np.isfinite(size).all()
            data={'gene':g,'source':p['source'],'structure':p['structure'],'site_rule':rule,
                'center':center.tolist(),'size':size.tolist(),'volume_A3':float(np.prod(size)),
                'basis':basis,'target_chains':sorted({c for c,n in target_residues})}
            (B/'prepared'/f'{g}_box.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf8')
            rows.append({'gene':g,'structure':p['structure'],'site_rule':rule,'pocket_rank':pocketrank,'pocket_probability':pocketprob,
                **{f'center_{a}_A':round(float(v),4) for a,v in zip('xyz',center)},
                **{f'size_{a}_A':round(float(v),4) for a,v in zip('xyz',size)},'volume_A3':round(float(np.prod(size)),1),
                'greater_than_30A_cube':bool(np.prod(size)>27000),'status':'search_region_fixed'})
            print(g,'BOX',list(np.round(size,2)),flush=True)
        except Exception as e:errors.append({'gene':g,'error':str(e)});print('ERROR',g,str(e),flush=True)
    if rows:
        with (B/'inputs/47_actual_search_regions.tsv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (B/'audit/Search_region_errors.json').write_text(json.dumps(errors,ensure_ascii=False,indent=2),encoding='utf8')
    if errors:raise SystemExit(1)
    assert len(rows)==47
if __name__=='__main__':
    candidates=[p for p in PLAN if p['source']=='PDB' and p['site_rule']=='P2Rank_first_target_associated_pocket']
    with ThreadPoolExecutor(max_workers=2) as ex:
        for fut in as_completed([ex.submit(predict,p) for p in candidates]):fut.result()
    summarize()
