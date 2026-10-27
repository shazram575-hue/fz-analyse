import csv,json,time,shutil
from pathlib import Path
import requests,gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation,PDBQTWriterLegacy
ROOT=Path(__file__).resolve().parents[3];B=Path(__file__).resolve().parents[1];OLD=ROOT/'outputs/docking_audit_2026-10-24'
P=B/'prepared/ligands';P.mkdir(exist_ok=True)
def get(ccd):
    p=B/'sources'/f'CCD_{ccd}.cif'
    if not p.exists():
        old=OLD/'sources'/f'CCD_{ccd}.cif'
        if old.exists():shutil.copy2(old,p)
        else:
            r=requests.get(f'https://files.rcsb.org/ligands/download/{ccd}.cif',timeout=50);r.raise_for_status();p.write_bytes(r.content)
    return gemmi.cif.read_file(str(p)).sole_block()
def from_ccd(ccd,refpath):
    block=get(ccd);st=gemmi.read_structure(str(refpath));res=st[0][0][0]
    positions={a.name:[a.pos.x,a.pos.y,a.pos.z] for a in res if not a.element.is_hydrogen}
    rw=Chem.RWMol();lookup={}
    for row in block.find(['_chem_comp_atom.atom_id','_chem_comp_atom.type_symbol','_chem_comp_atom.charge']):
        if row[1]=='H':continue
        a=Chem.Atom(row[1]);a.SetFormalCharge(int(row[2]));a.SetProp('atom_name',row[0]);lookup[row[0]]=rw.AddAtom(a)
    assert set(lookup)==set(positions),(ccd,'missing observed heavy atoms',set(lookup)-set(positions))
    order={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for row in block.find(['_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2','_chem_comp_bond.value_order']):
        if row[0] in lookup and row[1] in lookup:rw.AddBond(lookup[row[0]],lookup[row[1]],order[row[2]])
    m=rw.GetMol();Chem.SanitizeMol(m);conf=Chem.Conformer(m.GetNumAtoms());conf.Set3D(True)
    for name,i in lookup.items():conf.SetAtomPosition(i,positions[name])
    m.AddConformer(conf);Chem.AssignStereochemistryFrom3D(m);return m
def sdf(path,m):
    with path.open('w') as f:w=Chem.SDWriter(f);w.write(m);w.close()
def pdbqt(path,m):
    setup=MoleculePreparation().prepare(m)[0];text,ok,error=PDBQTWriterLegacy.write_string(setup)
    assert ok,error;path.write_text(text,encoding='ascii')

def main():
    for state in ['pH74','neutral']:
        for ext in ['sdf','pdbqt']:shutil.copy2(OLD/'prepared'/f'spiropachysine_{state}.{ext}',P/f'spiropachysine_{state}.{ext}')
    controls=[]
    for gene,ccd in [('PPARG','BRL'),('PNP','JTC'),('SLC18B1','YHL')]:
        meta=json.loads((B/'audit'/f'{gene}_extraction.json').read_text(encoding='utf8'))
        ref=min((x for x in meta['reference_ligands'] if x['CCD']==ccd),key=lambda x:(x['distance_A'],x['resnum']))
        m=from_ccd(ccd,ROOT/ref['path']);sdf(P/f'{gene}_{ccd}_experimental.sdf',m)
        m=Chem.AddHs(m);m.RemoveAllConformers();params=AllChem.ETKDGv3();params.randomSeed=20261024
        assert AllChem.EmbedMolecule(m,params)==0,(gene,'control embedding')
        AllChem.MMFFOptimizeMolecule(m,maxIters=2000)
        sdf(P/f'{gene}_{ccd}_regenerated.sdf',m);pdbqt(P/f'{gene}_{ccd}.pdbqt',m)
        controls.append({'gene':gene,'CCD':ccd,'reference':str((P/f'{gene}_{ccd}_experimental.sdf').relative_to(ROOT)),
                         'ligand':str((P/f'{gene}_{ccd}.pdbqt').relative_to(ROOT)),'state':'CCD formal charge; regenerated stereochemistry-preserving independent conformer'})
        print('CONTROL READY',gene,ccd,flush=True)
    (B/'inputs/Redocking_controls.json').write_text(json.dumps(controls,indent=2),encoding='utf8')
    # Fixed cofactors use observed coordinates and CCD aromatic flags. Vina ignores partial charges.
    # Heme coordinated N is non-acceptor N; oxygen atoms OA; aromatic carbons A; metal types retained.
    records=[]
    for p in sorted((B/'prepared/cofactors').glob('*_retained.pdb')):
        gene=p.stem.replace('_retained','');st=gemmi.read_structure(str(p));lines=[];serial=1
        cofactor_rows=[]
        for c in st[0]:
            for r in c:
                flags={}
                if len(r)>1:
                    block=get(r.name)
                    flags={x[0]:x[1] for x in block.find(['_chem_comp_atom.atom_id','_chem_comp_atom.pdbx_aromatic_flag'])}
                for a in r:
                    if a.element.is_hydrogen:continue
                    el=a.element.name
                    if el=='C':typ='A' if flags.get(a.name)=='Y' else 'C'
                    elif el=='N':typ='N' if r.name=='HEM' else ('NA' if flags.get(a.name)=='Y' else 'N')
                    elif el=='O':typ='OA'
                    elif el=='S':typ='SA'
                    elif el in ['Zn','Fe','Mg','Ca','Cu','Na','K','P']:typ=el
                    else:raise ValueError((gene,r.name,a.name,el))
                    line=f'HETATM{serial:5d} {a.name:>4s} {r.name:>3s} {c.name[:1]:1s}{r.seqid.num:4d}    {a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}{1.:6.2f}{0.:6.2f}    {0.:6.3f} {typ:>2s}'
                    lines.append(line);serial+=1
                cofactor_rows.append({'chain':c.name,'resnum':r.seqid.num,'CCD':r.name,'atoms':len(r)})
        out=p.with_suffix('.pdbqt');out.write_text('\n'.join(lines)+'\n',encoding='ascii')
        records.append({'gene':gene,'retained':cofactor_rows,'pdbqt':str(out.relative_to(ROOT)),
            'typing':'observed heavy atoms; CCD aromatic flags; coordinated heme N non-acceptor; partial charges0 ignored by Vina',
            'limitations':'fixed structural cofactors; no metal-directional or redox chemistry term'})
        print('COFACTOR READY',gene,len(lines),'atoms',flush=True)
    (B/'audit/Cofactor_preparation.json').write_text(json.dumps(records,indent=2),encoding='utf8')
if __name__=='__main__':main()
