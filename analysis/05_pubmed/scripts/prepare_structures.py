import json,csv,hashlib
from pathlib import Path
import gemmi,numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
B=Path(__file__).resolve().parents[1]
ROOT=B.parents[1]
P=B/'prepared';P.mkdir(exist_ok=True)
audit=[]
for code in ['9D7W','9D7X','8IJR','8W9Y']:
    s=gemmi.read_structure(str(B/'structures'/f'{code}.cif'))
    out=gemmi.Structure();out.name=code;model=gemmi.Model('1')
    canonical=json.loads((B/'sources'/('uniprot_Q6NT16.json' if code.startswith('9') else 'uniprot_Q96LT4.json')).read_text(encoding='utf8'))['sequence']['value']
    rows=[];coords=[]
    for c in s[0]:
        chain=gemmi.Chain(c.name)
        for res in c:
            if res.het_flag!='A':continue
            if code.startswith('9'):
                if not (232<=res.seqid.num<=638):continue
                number=res.seqid.num-207
            else:number=res.seqid.num
            letter=gemmi.find_tabulated_residue(res.name).one_letter_code
            assert canonical[number-1]==letter,(code,number,letter,canonical[number-1])
            rr=res.clone();rr.seqid=gemmi.SeqId(number,' ')
            for i in reversed(range(len(rr))):
                if rr[i].element.is_hydrogen:del rr[i]
            chain.add_residue(rr)
            rows.append({'chain':c.name,'original_resnum':res.seqid.num,'canonical_resnum':number,'resname':res.name})
            coords += [[a.pos.x,a.pos.y,a.pos.z] for a in rr]
        if len(chain):model.add_chain(chain)
    out.add_model(model);out.write_pdb(str(P/f'{code}_protein.pdb'))
    with (P/f'{code}_residue_mapping.tsv').open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    audit.append({'code':code,'chain_residue_counts':{c.name:len(c) for c in out[0]},'sequence_matches':True,'residue_span':[min(x['canonical_resnum'] for x in rows),max(x['canonical_resnum'] for x in rows)],'all_heteroatoms_removed':True,'fusion_removed':code.startswith('9')})

def ligand_from_ccd(code,chain,resname):
    structure=gemmi.read_structure(str(B/'structures'/f'{code}.cif'))
    res=next(r for r in structure[0][chain] if r.name==resname)
    position={a.name:[a.pos.x,a.pos.y,a.pos.z] for a in res if not a.element.is_hydrogen}
    ccd=gemmi.cif.read_file(str(B/'sources'/f'CCD_{resname}.cif')).sole_block()
    tab=ccd.find(['_chem_comp_atom.atom_id','_chem_comp_atom.type_symbol','_chem_comp_atom.charge'])
    rw=Chem.RWMol();lookup={}
    for row in tab:
        if row[1]=='H':continue
        atom=Chem.Atom(row[1]);atom.SetFormalCharge(int(row[2]));atom.SetProp('atom_name',row[0])
        lookup[row[0]]=rw.AddAtom(atom)
    bt={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for row in ccd.find(['_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2','_chem_comp_bond.value_order']):
        if row[0] in lookup and row[1] in lookup:rw.AddBond(lookup[row[0]],lookup[row[1]],bt[row[2]])
    m=rw.GetMol();Chem.SanitizeMol(m)
    conf=Chem.Conformer(m.GetNumAtoms());conf.Set3D(True)
    assert set(lookup)==set(position),(resname,set(lookup)-set(position))
    for name,index in lookup.items():conf.SetAtomPosition(index,position[name])
    m.AddConformer(conf);Chem.AssignStereochemistryFrom3D(m)
    m.SetProp('_Name',f'{code}_{chain}_{resname}_experimental')
    with (P/f'{code}_{resname}_reference.sdf').open('w') as f:
        w=Chem.SDWriter(f);w.write(m);w.close()
    heavycoords=np.array(m.GetConformer().GetPositions())
    center=(heavycoords.max(0)+heavycoords.min(0))/2
    box=np.maximum(heavycoords.max(0)-heavycoords.min(0)+12,24)
    return m,center,box
tbz,center,box=ligand_from_ccd('9D7W','A','YHL')
dag,dagcenter,dagbox=ligand_from_ccd('8IJR','A','Z0P')
ope,_,_=ligand_from_ccd('8IJR','A','OPE')
# Ligand-informed reaction-chamber box, common target ligand search volume >=24 A.
points=np.vstack([dag.GetConformer().GetPositions(),ope.GetConformer().GetPositions()])
samcenter=(points.max(0)+points.min(0))/2
sambox=np.maximum(points.max(0)-points.min(0)+12,24)
boxes={'SLC18B1_9D7W':{'receptor':'9D7W','center':center.tolist(),'size':box.tolist(),'basis':'experimental tetrabenazine heavy-atom extent + 6 A each side; min 24 A'},
 'SAMD8_8IJR':{'receptor':'8IJR','center':samcenter.tolist(),'size':sambox.tolist(),'basis':'experimental chain-A DAG+phosphoethanolamine extent +6 A each side; min 24 A'}}
(P/'ligand_informed_boxes.json').write_text(json.dumps(boxes,indent=2),encoding='utf8')
(B/'audit/experimental_structure_preparation.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
for name,m in [('TBZ',tbz),('DAG',dag)]:
    m=Chem.AddHs(m);m.RemoveAllConformers();params=AllChem.ETKDGv3();params.randomSeed=20261024
    assert AllChem.EmbedMolecule(m,params)==0
    AllChem.MMFFOptimizeMolecule(m,maxIters=2000)
    with (P/f'{name}_regenerated.sdf').open('w') as f:
        w=Chem.SDWriter(f);w.write(m);w.close()
print(json.dumps(boxes,indent=2));print(json.dumps(audit,indent=2))
