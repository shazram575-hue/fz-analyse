"""Meeko receptor CLI with standard-residue template connectivity rather than distance-inferred extra bonds.

Only this process changes the input graph builder; installed packages stay unchanged.
Coordinates and Meeko atom typing, charge assignment and residue templates are unchanged.
"""
import sys,json,re,itertools
import numpy as np
from scipy.spatial import cKDTree
from pathlib import Path
from rdkit import Chem
import meeko.polymer as polymer
from meeko import ResidueChemTemplates
from meeko.utils.rdkitutils import build_one_rdkit_mol_per_altloc
from meeko.cli.mk_prepare_receptor import main

templates=ResidueChemTemplates.create_from_defaults().residue_templates
original=polymer._aux_altloc_mol_build
def standard_connectivity(fields,requested_altloc,default_altloc):
    name=fields[0].resname
    if name not in templates:return original(fields,requested_altloc,default_altloc)
    mols=build_one_rdkit_mol_per_altloc(fields)
    if None in mols:key=None
    elif requested_altloc in mols:key=requested_altloc
    elif default_altloc in mols:key=default_altloc
    else:return original(fields,requested_altloc,default_altloc)
    mol,idx=mols[key]
    names=[a.atomname.strip() for a in fields if not a.altloc or a.altloc==key]
    assert len(names)==mol.GetNumAtoms(),(name,names,mol.GetNumAtoms())
    template=templates.get('C'+name) if 'OXT' in names else templates[name]
    if template is None:template=templates[name]
    lookup={n:i for i,n in enumerate(names)}
    rw=Chem.RWMol(mol)
    for bond in template.mol.GetBonds():
        n1=template.atom_names[bond.GetBeginAtomIdx()];n2=template.atom_names[bond.GetEndAtomIdx()]
        if n1 in lookup and n2 in lookup:
            # Meeko's raw graph is connectivity only; residue templates assign bond orders downstream.
            rw.AddBond(lookup[n1],lookup[n2],Chem.BondType.SINGLE)
    out=rw.GetMol()
    for a in out.GetAtoms():
        if a.GetAtomicNum()==7 and len(a.GetNeighbors())==4:a.SetFormalCharge(1)
    Chem.SanitizeMol(out)
    return out,idx,False,False
polymer._aux_altloc_mol_build=standard_connectivity
original_inter=polymer.find_inter_mols_bonds
def chemical_inter_residue_bonds(mols):
    # Same peptide/disulfide rules as v2, without Meeko's quadratic all-residue-pair scan.
    # The original covalent radii and 1.2 distance allowance are retained exactly.
    kept={};keys=list(mols);keyorder={r:i for i,r in enumerate(keys)};lookup={};sulfurs=[]
    for r,(mol,_) in mols.items():
        ch,n=r.rsplit(':',1);n=int(re.sub('[^0-9-]','',n));positions=mol.GetConformer().GetPositions()
        atoms={a.GetPDBResidueInfo().GetName().strip():a.GetIdx() for a in mol.GetAtoms()}
        lookup[(ch,n)]=(r,atoms,positions)
        if 'SG' in atoms:sulfurs.append((r,atoms['SG'],positions[atoms['SG']]))
    def add(r1,i,r2,j):
        if keyorder[r1]<keyorder[r2]:kept.setdefault((r1,r2),[]).append((i,j))
        else:kept.setdefault((r2,r1),[]).append((j,i))
    cutoff=1.2*(polymer.covalent_radius[6]+polymer.covalent_radius[7])
    for (ch,n),(r1,a1,p1) in lookup.items():
        nxt=lookup.get((ch,n+1))
        if nxt is None:continue
        r2,a2,p2=nxt
        for left,right in [('C','N'),('N','C')]:
            if left in a1 and right in a2 and np.linalg.norm(p1[a1[left]]-p2[a2[right]])<cutoff:add(r1,a1[left],r2,a2[right])
    if sulfurs:
        tree=cKDTree([x[2] for x in sulfurs]);cutoff=1.2*2*polymer.covalent_radius[16]
        for a,b in tree.query_pairs(cutoff):
            r1,i,_=sulfurs[a];r2,j,_=sulfurs[b]
            if r1!=r2:add(r1,i,r2,j)
    return kept
polymer.find_inter_mols_bonds=chemical_inter_residue_bonds
if __name__=='__main__':sys.exit(main())
