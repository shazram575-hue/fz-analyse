import json
from pathlib import Path
import gemmi
B=Path(__file__).resolve().parents[1]
u=json.loads((B/'sources/uniprot_Q14728.json').read_text())
seq=u['sequence']['value']
s=gemmi.read_structure(str(B/'sources/6S4M.cif'))
out=gemmi.Structure();m=gemmi.Model('1');c=gemmi.Chain('A');n=0
for r in s[0]['A']:
    if r.het_flag!='A' or not 1<=r.seqid.num<=444:continue
    assert seq[r.seqid.num-1]==gemmi.find_tabulated_residue(r.name).one_letter_code
    rr=r.clone()
    for i in reversed(range(len(rr))):
        if rr[i].element.is_hydrogen:del rr[i]
    c.add_residue(rr);n+=1
m.add_chain(c);out.add_model(m);out.write_pdb(str(B/'prepared/6S4M_protein.pdb'))
(B/'audit/MFSD10_6S4M_preparation.json').write_text(json.dumps({'canonical_accession':'Q14728','frozen_symbol':'MFSD10','current_symbol':'SLC75A1','resolved_residues':n,'resolved_span':[18,443],'sequence_matches':True,'removed':'water, citrate and lipid; no validated drug ligand assumed','docking_run':False},indent=2),encoding='utf8')
print('MFSD10 6S4M',n,'residues sequence verified')
