from pathlib import Path
import sys,json,hashlib,shutil
import pandas as pd
root=Path(sys.argv[1]);old=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26';out=root/'outputs/MFSD10_稀疏与聚合探索_2026-10-26'
freeze=json.loads((old/'frozen/Analysis_freeze_v1.json').read_text())
for p,h in freeze['SHA256'].items():assert hashlib.sha256((old/p).read_bytes()).hexdigest()==h
rows=[]
for cohort,path in [('GSE161529',old/'discovery/MFSD10_patient_level_associations.tsv'),('GSE176078',old/'validation/MFSD10_external_all_frozen_programs.tsv')]:
    d=pd.read_csv(path,sep='\t');d['cohort']=cohort;d['reused_without_refitting']=True;d['source']=str(path)
    d['analysis_unit']=d.association.map({'within':'Single-cell partial rho per patient, then patient meta-analysis','between':'Patient pseudobulk MFSD10 vs patient mean fixed program score'})
    rows.append(d)
pd.concat(rows).to_csv(out/'results/Existing_models_reused_all20.tsv',sep='\t',index=False)
files=[old/'frozen/Programs_v1.json',old/'frozen/Analysis_freeze_v1.json',old/'scripts/program_analysis.py',old/'frozen/Analysis_code_v1.py',old/'discovery/MFSD10_patient_level_associations.tsv',old/'validation/MFSD10_external_all_frozen_programs.tsv']
pd.DataFrame([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]).to_csv(out/'audit/Frozen_and_reused_before_SHA256.tsv',sep='\t',index=False)
p=out/'audit/PROGRESS.json';d=json.loads(p.read_text());d.update(old_models_reused=True,stage='building_nonoverlapping_state_units');p.write_text(json.dumps(d,indent=2),encoding='utf-8')
print('REUSED_20_ESTIMATES_AND_VERIFIED_FROZEN_HASHES')
