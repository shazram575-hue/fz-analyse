"""Incremental pose checks and source-evidence summaries. Never invent missing jobs."""
import csv,json,hashlib,statistics,itertools
from pathlib import Path
from runtime_io import read_json
import numpy as np,pandas as pd,gemmi
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import rdMolAlign,rdMolChemicalFeatures
from meeko import PDBQTMolecule,RDKitMolCreate
ROOT=Path(__file__).resolve().parents[3];B=Path(__file__).resolve().parents[1];OLD=ROOT/'outputs/docking_audit_2026-10-24'
# RDKit's native Windows filename reader cannot open this Unicode feature-file path.
# Read the identical bytes with Python, then construct the same factory from the string.
rdMolChemicalFeatures.BuildFeatureFactory=lambda p:rdMolChemicalFeatures.BuildFeatureFactoryFromString(Path(p).read_text(encoding='utf8'))
from posebusters import PoseBusters
config=json.loads((OLD/'audit/PoseBusters_basic_geometry_config.json').read_text(encoding='utf8'))
config['modules']=[m for m in config['modules'] if m['function']!='intermolecular_distance']
config['modules'][0]['chosen_binary_test_output']=['mol_pred_loaded']
config['max_workers']=0
(B/'audit/PoseBusters_geometry_config.json').write_text(json.dumps(config,indent=2),encoding='utf8')
buster=PoseBusters(config=config)
for d in ['poses','audit/pose_checks']:(B/d).mkdir(parents=True,exist_ok=True)
def tab(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
PLAN={r['gene']:r for r in tab(B/'inputs/47_structure_and_search_plan.tsv')}
INV={r['gene']:r for r in tab(B/'inputs/47_verified_input_inventory.tsv')}
MUTATIONS={}
for r in tab(B/'inputs/Experimental_residue_mapping.tsv'):
    if r['matches']=='False':MUTATIONS.setdefault(r['gene'],set()).add((r['chain'],int(r['canonical_resnum'])))
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def write(p,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(rows)
def sdf_load(p):return next(Chem.ForwardSDMolSupplier(p.open('rb'),removeHs=False))
def load_dock(p):
    mols=RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(str(p),skip_typing=True))
    assert len(mols)==1 and mols[0] is not None
    return mols[0]
def one_pose(m,idx=0):
    mm=Chem.Mol(m);mm.RemoveAllConformers();mm.AddConformer(Chem.Conformer(m.GetConformer(idx)));return mm
def rms(a,b,aid=0,bid=0):return float(rdMolAlign.CalcRMS(Chem.RemoveHs(a),Chem.RemoveHs(b),prbId=aid,refId=bid,maxMatches=100000))
def receptor_atoms(g):
    paths=[B/'prepared/receptors'/f'{g}_prepared.pdb']
    cf=B/'prepared/cofactors'/f'{g}_retained.pdb'
    if cf.exists():paths.append(cf)
    atoms=[]
    for p in paths:
        st=gemmi.read_structure(str(p))
        for c in st[0]:
            for r in c:
                for a in r:
                    if a.element.is_hydrogen:continue
                    atoms.append((c.name,r.seqid.num,r.name,a.name,a.element.atomic_number,np.array([a.pos.x,a.pos.y,a.pos.z]),p==cf))
    return atoms
REC={}
def contacts(g,m):
    if g not in REC:
        atoms=receptor_atoms(g);xyz=np.array([a[5] for a in atoms]);REC[g]=(atoms,cKDTree(xyz))
    atoms,tree=REC[g];lig=Chem.RemoveHs(m);xyz=lig.GetConformer().GetPositions();dist,ids=tree.query(xyz)
    extraction=json.loads((B/'audit'/f'{g}_extraction.json').read_text(encoding='utf8'))
    targetchains=set(extraction['target_chains']);near={};minimum_scaled=100.;target_lig=set();partner_lig=set();ptable=Chem.GetPeriodicTable()
    for i,close in enumerate(tree.query_ball_point(xyz,5.0)):
        for j in close:
            a=atoms[j];d=float(np.linalg.norm(xyz[i]-a[5]));key=(a[0],a[1],a[2])
            if d<near.get(key,100.):near[key]=d
            metal=a[4] in [11,12,19,20,25,26,27,28,29,30]
            radius=(ptable.GetRcovalent(lig.GetAtomWithIdx(i).GetAtomicNum())+ptable.GetRcovalent(a[4])) if metal else (ptable.GetRvdw(lig.GetAtomWithIdx(i).GetAtomicNum())+ptable.GetRvdw(a[4]))
            minimum_scaled=min(minimum_scaled,d/radius)
            if d<=4 and not a[6]:
                (target_lig if a[0] in targetchains else partner_lig).add(i)
    target_contacts=[(k,d) for k,d in near.items() if k[0] in targetchains]
    confs=[]
    if PLAN[g]['source']=='AlphaFold':
        st=gemmi.read_structure(str(OLD/'structures'/f'{PLAN[g]["accession"]}_AF.pdb'))
        confidence={r.seqid.num:r['CA'][0].b_iso for r in st[0][0] if 'CA' in r}
        confs=[confidence[k[1]] for k,d in target_contacts if d<=4 and k[1] in confidence]
    closest=sorted(near.items(),key=lambda x:x[1])
    mutated=[(k,d) for k,d in closest if (k[0],k[1]) in MUTATIONS.get(g,set())]
    return {'min_receptor_heavy_distance_A':float(min(dist)), 'min_scaled_receptor_distance':minimum_scaled,
      'no_major_receptor_clash':minimum_scaled>=.75,'ligand_heavy_atoms_within4A_target':len(target_lig),
      'ligand_heavy_atoms_within4A_partner':len(partner_lig),'contact_residue_AF_pLDDT_mean':statistics.mean(confs) if confs else None,
      'sequence_difference_residues_within5A':';'.join(f'{k[0]}:{k[2]}{k[1]}({d:.2f})' for k,d in mutated),
      'contact_residues_within4A':';'.join(f'{k[0]}:{k[2]}{k[1]}({d:.2f})' for k,d in closest if d<=4),
      'closest_residues_within5A':';'.join(f'{k[0]}:{k[2]}{k[1]}({d:.2f})' for k,d in closest[:12])}

def evaluate_main(meta):
    g=meta['gene'];tag=meta['tag'];p=B/'docking'/f'{tag}.pdbqt';out=B/'audit/pose_checks'/f'{tag}.json'
    if out.exists():
        r=json.loads(out.read_text(encoding='utf8'))
        if r.get('pose_SHA256')==sha(p) and r.get('QA_version')==4:return r
    m=load_dock(p);top=one_pose(m);reference=sdf_load(B/'prepared/ligands'/f'spiropachysine_{meta["state"]}.sdf')
    refkey=Chem.MolToInchiKey(reference);checked=Chem.Mol(top)
    Chem.AssignAtomChiralTagsFromStructure(checked,replaceExistingTags=True)
    identity=Chem.MolToInchiKey(checked)==refkey
    pb=buster.bust(mol_pred=top,mol_true=reference,full_report=True)
    binary=[];chosen={key for module in config['modules'] for key in module['chosen_binary_test_output']}
    chosen|={module['rename_outputs'][key].lower().replace(' ','_') for module in config['modules']
             for key in module['chosen_binary_test_output'] if key in module.get('rename_outputs',{})}
    for name in pb.columns:
        value=pb.iloc[0][name]
        if str(name) in chosen and isinstance(value,(bool,np.bool_)):binary.append((str(name),bool(value)))
    assert len(binary)>=8,('PoseBusters tests missing',list(pb.columns),binary)
    failures=[name for name,ok in binary if not ok]
    # Write every returned pose with score and seed; top-pose checks refer to rank1.
    sdfpath=B/'poses'/f'{tag}.sdf'
    with sdfpath.open('w') as f:
        w=Chem.SDWriter(f)
        for i in range(m.GetNumConformers()):
            m.SetProp('pose_rank',str(i+1));m.SetProp('vina_kcal_mol',str(meta['scores_kcal_mol'][i]));m.SetProp('seed',str(meta['seed']));w.write(m,confId=i)
        w.close()
    row={'gene':g,'state':meta['state'],'seed':meta['seed'],'tag':tag,'score_kcal_mol':meta['best_vina_kcal_mol'],
      'returned_poses':m.GetNumConformers(),'QA_version':4,'stereochemistry_and_identity_match':identity,'PoseBusters_geometry_pass':not failures,
      'PoseBusters_failed_checks':';'.join(failures),**contacts(g,top),'pose_SHA256':sha(p),'sdf_path':str(sdfpath.relative_to(ROOT))}
    out.write_text(json.dumps(row,ensure_ascii=False,indent=2),encoding='utf8');print('QA',tag,identity,not failures,flush=True)
    return row

def main():
    rows=[];errors=[];controls=[]
    for path in sorted((B/'docking').glob('*.json')):
        meta=read_json(path)
        if meta.get('exit_code')!=0 or not meta.get('poses_returned'):continue
        try:rows.append(evaluate_main(meta))
        except Exception as e:errors.append({'tag':path.stem,'error':str(e)});print('ERROR',path.stem,str(e),flush=True)
    for path in sorted((B/'controls').glob('*.json')):
        meta=read_json(path)
        if meta.get('exit_code')!=0 or not meta.get('poses_returned'):continue
        try:
            m=load_dock(path.with_suffix('.pdbqt'));ref=sdf_load(B/'prepared/ligands'/f'{meta["gene"]}_{meta["state"]}_experimental.sdf')
            values=[rms(m,ref,i) for i in range(m.GetNumConformers())]
            controls.append({'gene':meta['gene'],'CCD':meta['state'],'seed':meta['seed'],'top_pose_RMSD_A':values[0],
              'best_returned_pose_RMSD_A':min(values),'top_pose_within2A':values[0]<=2,'any_pose_within2A':min(values)<=2,
              'top_score_kcal_mol':meta['best_vina_kcal_mol'],'RMSD_method':'symmetry-corrected heavy atoms in receptor frame, no ligand alignment'})
        except Exception as e:errors.append({'tag':path.stem,'error':str(e)})
    write(B/'results/All_completed_top_pose_checks.tsv',rows);write(B/'results/Redocking_control_results.tsv',controls)
    summaries=[];rmsrows=[]
    for g in PLAN:
        for state in ['pH74','neutral']:
            rr=[r for r in rows if r['gene']==g and r['state']==state]
            if not rr:continue
            pairvals=[]
            for a,b in itertools.combinations(rr,2):
                value=rms(load_dock(B/'docking'/f'{a["tag"]}.pdbqt'),load_dock(B/'docking'/f'{b["tag"]}.pdbqt'))
                pairvals.append(value);rmsrows.append({'gene':g,'state':state,'seed1':a['seed'],'seed2':b['seed'],'RMSD_A':value})
            scores=[r['score_kcal_mol'] for r in rr]
            summaries.append({'gene':g,'state':state,'completed_seeds':len(rr),'expected_seeds':3,'score_median_kcal_mol':statistics.median(scores),
              'score_min_kcal_mol':min(scores),'score_max_kcal_mol':max(scores),'pairwise_pose_RMSD_median_A':statistics.median(pairvals) if pairvals else None,
              'pairwise_pose_RMSD_max_A':max(pairvals) if pairvals else None,
              'all_top_poses_geometry_pass':all(r['PoseBusters_geometry_pass'] for r in rr),'all_top_poses_identity_pass':all(r['stereochemistry_and_identity_match'] for r in rr),
              'all_top_poses_no_major_clash':all(r['no_major_receptor_clash'] for r in rr),
              'min_ligand_heavy_atoms_within4A_target':min(r['ligand_heavy_atoms_within4A_target'] for r in rr),
              'contact_residue_AF_pLDDT_mean':statistics.mean([r['contact_residue_AF_pLDDT_mean'] for r in rr if r['contact_residue_AF_pLDDT_mean'] is not None]) if any(r['contact_residue_AF_pLDDT_mean'] is not None for r in rr) else None})
    write(B/'results/Target_state_summary.tsv',summaries);write(B/'results/Seed_pose_pairwise_RMSD.tsv',rmsrows)
    (B/'audit/Evaluation_errors.json').write_text(json.dumps(errors,ensure_ascii=False,indent=2),encoding='utf8')
    (B/'audit/Evaluation_progress.json').write_text(json.dumps({'main_top_poses_checked':len(rows),'expected':282,'controls_checked':len(controls),'errors':len(errors)},indent=2),encoding='utf8')
    print('CHECKED',len(rows),'main top poses;',len(controls),'controls;',len(errors),'errors',flush=True)
    if errors:raise SystemExit(1)
if __name__=='__main__':main()
