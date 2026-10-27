"""Documentation/integrity closeout only: never refit or optimize an association."""
from pathlib import Path
import sys,json,re,hashlib,importlib.metadata as im
import numpy as np,pandas as pd
root=Path(sys.argv[1]);out=root/'outputs/MFSD10_稀疏与聚合探索_2026-10-26'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
done=json.loads((out/'audit/COMPLETED.json').read_text());tests=pd.read_csv(out/'audit/Final_validation.tsv',sep='\t')
assert done['status']=='complete' and done['checks_passed']==32 and len(tests)==32 and tests.passed.all()
prior=pd.read_csv(out/'audit/Prior_results_preservation.tsv',sep='\t');assert len(prior)==5473 and prior.unchanged.all()
before=pd.read_csv(out/'audit/Frozen_and_reused_before_SHA256.tsv',sep='\t')
assert all(sha(r.path)==r.sha256 for _,r in before.iterrows())
stop=json.loads((out/'audit/STOP_THIS_ROUTE.json').read_text());assert stop['stop_automatic_expression_program_association_optimization'] and not stop['NMF_rerun'] and not stop['signature_changed']
a=pd.read_csv(out/'results/Aggregation_main_all10.tsv',sep='\t');c=pd.read_csv(out/'results/Comparison_all30_effects_CI_P_FDR_N.tsv',sep='\t')
assert len(a)==10 and len(c)==30
z=c[c.analysis=='within_patient_aggregate'].merge(a,on=['cohort','program_id'],suffixes=('_comparison','_source'))
for col in ['n_patients','effect_rho','ci_low','ci_high','p']:assert np.allclose(z[col+'_comparison'],z[col+'_source'])
assert np.allclose(z.FDR,z.FDR_main10)
links=[]
for p in out.glob('*.md'):
 for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
  if target.startswith(('https:','http:')):continue
  links.append(dict(document=p.name,target=target,exists=(p.parent/target).exists()))
assert all(r['exists'] for r in links)
pd.DataFrame(links).to_csv(out/'audit/Closeout_document_links.tsv',sep='\t',index=False)
for p in (out/'scripts').glob('*.py'):compile(p.read_text(encoding='utf-8'),str(p),'exec')
versions={p:im.version(p) for p in ['numpy','scipy','pandas','scikit-learn','anndata','matplotlib','statsmodels']}
(out/'logs/Software_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
progress=json.loads((out/'audit/PROGRESS.json').read_text());progress.update(closeout_complete=True,stage='complete',status='complete',stop_optimization=True,background_tasks=False)
(out/'audit/PROGRESS.json').write_text(json.dumps(progress,indent=2),encoding='utf-8')
files=[]
for p in sorted(out.rglob('*')):
 if not p.is_file() or '__pycache__' in p.parts or p.name in ['Deliverable_SHA256.tsv','CLOSEOUT_COMPLETED.json','06_closeout.log']:continue
 files.append(dict(relative_path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
pd.DataFrame(files).to_csv(out/'audit/Deliverable_SHA256.tsv',sep='\t',index=False)
result=dict(status='complete',closeout_only=True,scientific_analysis_rerun=False,original_audit_checks_passed=32,original_prior_files_preserved=5473,critical_frozen_hashes_rechecked=len(before),comparison_rows_verified=30,local_links_verified=len(links),deliverables_hashed=len(files),handoff_section=38,stop_optimization=True,background_tasks=False)
(out/'audit/CLOSEOUT_COMPLETED.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
