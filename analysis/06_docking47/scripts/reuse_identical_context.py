"""Reuse parameterization only for byte-identical protein atom records, never docking scores."""
import json,hashlib,shutil
from pathlib import Path
B=Path(__file__).resolve().parents[1]
def fingerprint(p):
    lines=[l for l in p.read_text(encoding='ascii').splitlines() if l.startswith(('ATOM  ','HETATM','TER'))]
    return hashlib.sha256('\n'.join(lines).encode('ascii')).hexdigest()
genes=['MT-ND2','MT-ND4','MT-ND4L'];ready=[]
for g in genes:
    p=B/'audit'/f'{g}_receptor_preparation.json'
    if p.exists() and json.loads(p.read_text(encoding='utf8')).get('exit_code')==0:ready.append(g)
if ready:
    source=ready[0];source_raw=B/'prepared/receptors'/f'{source}_raw.pdb'
    for g in genes:
        if g in ready:continue
        dest_raw=B/'prepared/receptors'/f'{g}_raw.pdb'
        if not dest_raw.exists() or fingerprint(source_raw)!=fingerprint(dest_raw):continue
        original=json.loads((B/'audit'/f'{source}_receptor_preparation.json').read_text(encoding='utf8'))
        own=json.loads((B/'audit'/f'{g}_extraction.json').read_text(encoding='utf8'))
        for suffix in ['.pdbqt','_fixed.pdb','_prepared.pdb']:
            shutil.copy2(B/'prepared/receptors'/f'{source}{suffix}',B/'prepared/receptors'/f'{g}{suffix}')
        record={**original,**own,'identical_context_parameterization_reused_from':source,
                'raw_protein_atom_records_SHA256':fingerprint(source_raw)}
        (B/'audit'/f'{g}_receptor_preparation.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf8')
        print('Reused identical receptor context',source,'->',g,flush=True)
