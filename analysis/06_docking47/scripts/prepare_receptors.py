"""Prepare actual target coordinates, preserve native assembly context, no missing-loop building."""
import os
os.environ.setdefault('OPENMM_CPU_THREADS','2')
import csv,json,sys,subprocess,time,gzip,xml.etree.ElementTree as ET,traceback,faulthandler
os.environ['PYTHONUTF8']='1'
sys.stdout.reconfigure(encoding='utf8');sys.stderr.reconfigure(encoding='utf8')
faulthandler.dump_traceback_later(120,repeat=True)
from pathlib import Path
import gemmi,numpy as np
from scipy.spatial import cKDTree
from pdbfixer import PDBFixer
from openmm.app import PDBFile

ROOT=Path(__file__).resolve().parents[3]; B=Path(__file__).resolve().parents[1]
OLD=ROOT/'outputs/docking_audit_2026-10-24'
def tab(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
PLAN=tab(B/'inputs/47_structure_and_search_plan.tsv')
MAPPINGS=tab(B/'inputs/Experimental_residue_mapping.tsv')
EXP={r['gene']:r for r in tab(B/'inputs/25_experimental_structure_review.tsv')}
for d in ['receptors','targets','cofactors','references']:(B/'prepared'/d).mkdir(exist_ok=True,parents=True)

def human_residues(code):
    out=set();root=ET.fromstring(gzip.decompress((B/'sources'/f'{code}_SIFTS.xml.gz').read_bytes()))
    for item in root.iter():
        if item.tag.split('}')[-1]!='residue':continue
        crosses=[x.attrib for x in item if x.tag.split('}')[-1]=='crossRefDb']
        if not any(x.get('dbSource')=='NCBI' and x.get('dbAccessionId')=='9606' for x in crosses):continue
        if not any(x.get('dbSource')=='UniProt' for x in crosses):continue
        p=next((x for x in crosses if x.get('dbSource')=='PDB'),None)
        if p:out.add((p['dbChainId'],p['dbResNum']))
    return out

def clean_res(res,number=None):
    r=res.clone()
    if number is not None:r.seqid=gemmi.SeqId(number,' ')
    # Pick the highest-occupancy alternate for each atom, retaining one coordinate.
    byname={}
    for a in r:
        if a.element.is_hydrogen:continue
        if a.name not in byname or a.occ>byname[a.name].occ:byname[a.name]=a.clone()
    while len(r):del r[-1]
    for a in byname.values():a.altloc='\0';r.add_atom(a)
    return r

def write_structure(path,model):
    st=gemmi.Structure();st.name=path.stem;st.add_model(model);st.write_pdb(str(path))

def extract(plan):
    g=plan['gene'];out=B/'prepared/receptors'/f'{g}_raw.pdb';target=B/'prepared/targets'/f'{g}_target.pdb'
    if plan['source']=='AlphaFold':
        src=OLD/'structures'/f'{plan["accession"]}_AF.pdb';st=gemmi.read_structure(str(src))
        model=gemmi.Model('1');ch=gemmi.Chain('A')
        for r in st[0][0]:
            if r.het_flag=='A':ch.add_residue(clean_res(r))
        model.add_chain(ch);write_structure(out,model);write_structure(target,model)
        return {'gene':g,'source':str(src.relative_to(ROOT)),'assembly':'AF monomer','target_chains':['A'],
                'context_chains':[],'residues':sum(len(c) for c in model),'retained_cofactors':[],
                'removed_heteroatoms':[],'reference_ligands':[]}
    code=plan['structure'];st=gemmi.read_structure(str(B/'structures'/f'{code}.cif'))
    mapping={(r['chain'],r.get('author_resid',r['author_resnum'])):int(r['canonical_resnum']) for r in MAPPINGS if r['gene']==g}
    policy=plan['assembly_policy'];targetchains=set(plan['target_chains'].split(';'))
    assembly=st.assemblies[0] if st.assemblies and policy!='target_chain' else None
    model=gemmi.make_assembly(assembly,st[0],gemmi.HowToNameCopiedChain.Dup) if assembly else st[0]
    allowed=human_residues(code) if policy=='target_and_native_human_complex_partners' else set()
    targetmodel=gemmi.Model('1');newmodel=gemmi.Model('1');context=[];tg=[];removed=[];kept=[];refs=[]
    native_xyz=[]
    for c in model:
        tc=gemmi.Chain(c.name);nc=gemmi.Chain(c.name)
        for res in c:
            if res.het_flag!='A':continue
            key=(c.name,str(res.seqid))
            if key in mapping:
                rr=clean_res(res,mapping[key]);tc.add_residue(rr);nc.add_residue(rr.clone())
                native_xyz.extend([[a.pos.x,a.pos.y,a.pos.z] for a in rr])
            elif (c.name,str(res.seqid)) in allowed and c.name not in targetchains:
                nc.add_residue(clean_res(res))
        if len(tc):targetmodel.add_chain(tc);tg.append(c.name)
        if len(nc):newmodel.add_chain(nc)
        if len(nc) and c.name not in tg:context.append(c.name)
    # A single-chain calculation uses the longest canonical chain; retain homomers only as planned.
    if policy=='target_chain' and len(tg)>1:
        best=EXP[g]['best_chain'];targetmodel=gemmi.Model('1');newmodel=gemmi.Model('1')
        cc=next(c for c in model if c.name==best);chain=gemmi.Chain(best)
        for r in cc:
            key=(best,str(r.seqid))
            if r.het_flag=='A' and key in mapping:chain.add_residue(clean_res(r,mapping[key]))
        targetmodel.add_chain(chain);newmodel.add_chain(chain.clone());tg=[best]
        native_xyz=[[a.pos.x,a.pos.y,a.pos.z] for r in chain for a in r]
    assert tg and native_xyz,(g,'target extraction empty')
    assert len(set(c.name for c in newmodel))==len(newmodel),(g,'assembly chain naming collision')
    # Non-water reference ligands remain in separate files, not inside the docking cavity.
    tree=cKDTree(native_xyz);cofactor_model=gemmi.Model('1')
    metal_names={'CA12':{'ZN'},'CYBRD1':{'HEM','ZN'},'ITPR3':{'ATP','CA','ZN'},
                 'KCNK1':{'K'},'SCD':{'ZN','FE'}}.get(g,set())
    desired=set(plan['known_ligand_CCD'].split('+')) if plan['known_ligand_CCD'] else set()
    best=EXP[g]['best_chain']
    for c in model:
        cf=gemmi.Chain(c.name)
        for res in c:
            if res.het_flag=='A' or res.name in ['HOH','WAT','DOD']:continue
            atoms=[a for a in res if not a.element.is_hydrogen]
            if not atoms:continue
            xyz=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in atoms]);distance=float(tree.query(xyz)[0].min())
            if res.name in desired and c.name==best:
                refmodel=gemmi.Model('1');refchain=gemmi.Chain(c.name);refchain.add_residue(clean_res(res));refmodel.add_chain(refchain)
                refpath=B/'prepared/references'/f'{g}_{res.name}_{res.seqid.num}.pdb'
                write_structure(refpath,refmodel)
                refs.append({'CCD':res.name,'chain':c.name,'resnum':res.seqid.num,'path':str(refpath.relative_to(ROOT)),'distance_A':distance})
            if res.name in metal_names and distance<8:
                cf.add_residue(clean_res(res));kept.append({'CCD':res.name,'chain':c.name,'resnum':res.seqid.num,'atoms':len(atoms)})
            elif distance<12:removed.append({'CCD':res.name,'chain':c.name,'resnum':res.seqid.num})
        if len(cf):cofactor_model.add_chain(cf)
    if len(cofactor_model):write_structure(B/'prepared/cofactors'/f'{g}_retained.pdb',cofactor_model)
    write_structure(out,newmodel);write_structure(target,targetmodel)
    return {'gene':g,'source':code,'assembly':assembly.name+':'+assembly.oligomeric_details if assembly else 'target chain from deposited coordinates',
      'target_chains':tg,'context_chains':context,'residues':sum(len(c) for c in newmodel),
      'target_residues':sum(len(c) for c in targetmodel),'retained_cofactors':kept,'removed_heteroatoms':removed,'reference_ligands':refs}

def repair_and_parameterize(plan,record):
    g=plan['gene'];base=B/'prepared/receptors'/g;raw=base.with_name(g+'_raw.pdb');fixed=base.with_name(g+'_fixed.pdb')
    meta=B/'audit'/f'{g}_receptor_preparation.json'
    if base.with_suffix('.pdbqt').exists() and meta.exists():
        prior=json.loads(meta.read_text(encoding='utf8'))
        if prior.get('connectivity_builder') in ['standard_residue_templates_v1','standard_residue_templates_v2'] and prior.get('exit_code')==0:return prior
    start=time.monotonic()
    if not fixed.exists():
        fixer=PDBFixer(filename=str(raw));fixer.findMissingResidues()
        missing_segments={str(k):v for k,v in fixer.missingResidues.items()};fixer.missingResidues={}
        fixer.findNonstandardResidues()
        assert not fixer.nonstandardResidues,(g,'nonstandard residues require explicit treatment',fixer.nonstandardResidues)
        fixer.findMissingAtoms()
        missing_atoms={f'{r.chain.id}:{r.id}':[a.name for a in atoms] for r,atoms in fixer.missingAtoms.items()}
        terminal={f'{r.chain.id}:{r.id}':atoms for r,atoms in fixer.missingTerminals.items()}
        fixer.missingTerminals={}
        if fixer.missingAtoms:fixer.addMissingAtoms(seed=20261024)
        with fixed.open('w') as f:PDBFile.writeFile(fixer.topology,fixer.positions,f,keepIds=True)
        record.update({'missing_residue_segments_not_built':missing_segments,'missing_heavy_atoms_added':missing_atoms,
                       'missing_terminal_atoms_not_added':terminal,'missing_terminal_atoms_added':{}})
        (B/'audit'/f'{g}_repair.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf8')
    else:record=json.loads((B/'audit'/f'{g}_repair.json').read_text(encoding='utf8'))
    # Remove only artificial terminal atoms from the earlier pilot preparation attempt.
    earlier_terminals=record.get('missing_terminal_atoms_added',{})
    if earlier_terminals:
        lines=fixed.read_text().splitlines()
        lines=[line for line in lines if not (line.startswith(('ATOM  ','HETATM')) and
          line[12:16].strip() in earlier_terminals.get(f'{line[21:22]}:{line[22:26].strip()}',[]))]
        fixed.write_text('\n'.join(lines)+'\n',encoding='ascii')
        record['missing_terminal_atoms_not_added']=earlier_terminals;record['missing_terminal_atoms_added']={}
        (B/'audit'/f'{g}_repair.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf8')
    cmd=[sys.executable,str(B/'scripts/parameterize_receptor.py'),'--read_pdb',str(fixed),'-o',str(base),'-p','--write_pdb',str(base.with_name(g+'_prepared.pdb'))]
    with (B/'logs'/f'prepare_{g}.log').open('w',encoding='utf8') as f:proc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,cwd=ROOT)
    record.update({'preparation_command':cmd,'exit_code':proc.returncode,'elapsed_seconds':time.monotonic()-start,
                   'connectivity_builder':'standard_residue_templates_v2'})
    meta.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf8')
    if proc.returncode:raise RuntimeError(g+': '+(B/'logs'/f'prepare_{g}.log').read_text(encoding='utf8',errors='replace')[-2200:])
    return record

if __name__=='__main__':
    requested=set(sys.argv[1:]);plans=[r for r in PLAN if not requested or r['gene'] in requested]
    mode=os.environ.get('DOCK47_PREP_MODE','all')
    errors=[];summary=[]
    for plan in plans:
        g=plan['gene']
        try:
            print(g,'PREPARING',flush=True)
            record=extract(plan)
            (B/'audit'/f'{g}_extraction.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf8')
            if mode!='extract':record=repair_and_parameterize(plan,record)
            summary.append(record);print(g,'READY' if mode!='extract' else 'EXTRACTED',record['residues'],'residues',flush=True)
        except Exception as e:
            traceback.print_exc()
            errors.append({'gene':g,'error':str(e)});print('ERROR',g,str(e),flush=True)
        (B/'audit/Preparation_errors.json').write_text(json.dumps(errors,ensure_ascii=False,indent=2),encoding='utf8')
    (B/'audit/Receptor_preparation_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
    if errors:raise SystemExit(1)
