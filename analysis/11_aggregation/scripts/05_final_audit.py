from pathlib import Path
import sys,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad
root=Path(sys.argv[1]);out=root/'outputs/MFSD10_稀疏与聚合探索_2026-10-26';old=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26';checks=[]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
def check(name,value):
 checks.append(dict(check=name,passed=bool(value)))
 assert value,name
before=pd.read_csv(out/'audit/Frozen_and_reused_before_SHA256.tsv',sep='\t')
check('Frozen_models_unchanged',all(sha(r.path)==r.sha256 for _,r in before.iterrows()))
check('Original_helper_identical',sha(out/'scripts/reused_original_models.py')==sha(old/'frozen/Analysis_code_v1.py'))
prior=pd.read_csv(old/'audit/Deliverable_SHA256.tsv',sep='\t');pres=[]
for _,r in prior.iterrows():pres.append(dict(path=r.relative_path,unchanged=sha(old/r.relative_path)==r.sha256))
check('All_prior_results_preserved',all(r['unchanged'] for r in pres));pd.DataFrame(pres).to_csv(out/'audit/Prior_results_preservation.tsv',sep='\t',index=False)
for coh,fn,n in [('discovery','GSE161529_confirmed_malignant_singlets_counts.h5ad',8722),('validation','GSE176078_malignant_candidates_author_counts_clinical_v2.h5ad',24489)]:
 m=pd.read_csv(out/'metadata'/f'{coh}_cell_unit_membership.tsv.gz',sep='\t',index_col=0);src=ad.read_h5ad(old/'objects'/fn);total=np.asarray(src.X.sum(axis=0))
 check(coh+'_all_cells_once',len(m)==n and m.index.is_unique and set(m.index)==set(src.obs_names))
 check(coh+'_no_cross_patient',m.groupby('state_unit').patient_id.nunique().eq(1).all())
 check(coh+'_no_cross_sample',m.groupby('state_unit').unit_sample.nunique().eq(1).all())
 g=m.groupby(['patient_id','unit_sample']).agg(n=('state_unit','size'),k=('state_unit','nunique'));check(coh+'_fixed_grain20',g.k.eq(np.maximum(1,g.n//20)).all())
 hv=pd.read_csv(out/'metadata'/f'{coh}_neighborhood_features.tsv.gz',sep='\t');check(coh+'_target_excluded',not hv.feature_id.isin(['ENSG00000109736','MFSD10','SLC75A1']).any())
 for kind in ['state','random']:
  a=ad.read_h5ad(out/'objects'/f'{coh}_{kind}_raw_count_units.h5ad')
  check(coh+'_'+kind+'_counts_integer',np.isfinite(a.X.data).all() and (a.X.data>=0).all() and (a.X.data==np.rint(a.X.data)).all())
  check(coh+'_'+kind+'_UMI_conserved',np.array_equal(total,np.asarray(a.X.sum(axis=0))))
  check(coh+'_'+kind+'_cells_conserved',int(a.obs.n_cells.sum())==n);del a
 a=ad.read_h5ad(out/'objects'/f'{coh}_patient_pseudobulk_raw_counts.h5ad');check(coh+'_PB_conserved',np.array_equal(total,np.asarray(a.X.sum(axis=0))));del src,a;gc.collect()
check('Thirty_main_comparison_estimates',len(pd.read_csv(out/'results/Comparison_all30_effects_CI_P_FDR_N.tsv',sep='\t'))==30)
check('Sixty_all_program_sensitivity_tests',len(pd.read_csv(out/'results/Aggregation_all60_fixed_specifications.tsv',sep='\t'))==60)
a=pd.read_csv(out/'results/Aggregation_main_all10.tsv',sep='\t');check('Primary_biological_n',len(a)==10 and a[a.cohort=='discovery'].n_patients.eq(13).all() and a[a.cohort=='validation'].n_patients.eq(15).all())
f=pd.read_csv(out/'results/Primary_patients_not_estimable.tsv',sep='\t');check('Insufficient_units_reported',len(f)==1 and f.iloc[0].patient_id=='GSE176078_Case__4465' and f.iloc[0].n_units==6)
check('Stop_applied',json.loads((out/'audit/STOP_THIS_ROUTE.json').read_text())['stop_automatic_expression_program_association_optimization'])
for p in (out/'scripts').glob('*.py'):compile(p.read_text(encoding='utf-8'),str(p),'exec')
pd.DataFrame(checks).to_csv(out/'audit/Final_validation.tsv',sep='\t',index=False)
p=out/'audit/PROGRESS.json';d=json.loads(p.read_text());d.update(status='complete',stage='complete',stop_optimization=True,background_tasks=False);p.write_text(json.dumps(d,indent=2),encoding='utf-8')
done=dict(status='complete',analysis_type='exploratory_reanalysis',checks_passed=len(checks),all_prior_results_unchanged=len(prior),NMF_rerun=False,signature_changed=False,primary_n_discovery=13,primary_n_validation_aggregate=15,stop_optimization=True,background_tasks=False)
(out/'audit/COMPLETED.json').write_text(json.dumps(done,indent=2),encoding='utf-8');print(json.dumps(done,indent=2))
