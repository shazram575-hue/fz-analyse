from pathlib import Path
import sys,json,hashlib,re,csv,importlib.metadata as im
import numpy as np,pandas as pd,anndata as ad
from cnmf import cNMF
from cnmf.cnmf import load_df_from_npz

root=Path(sys.argv[1]);out=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26';checks=[]
def add(name,passed,detail=''):
    checks.append(dict(check=name,pass_=bool(passed),detail=str(detail)))
    if not passed:raise RuntimeError(name+': '+str(detail))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

freeze=json.loads((out/'frozen/Analysis_freeze_v1.json').read_text());bank=json.loads((out/'frozen/Programs_v1.json').read_text())
for name,value in freeze['SHA256'].items():add('Frozen_SHA256_'+name,sha(out/name)==value)
add('Five_fixed_programs',len(bank['meta_programs'])==5)
add('Program_target_exclusion',all('ENSG00000109736' not in m['genes'] and 'MFSD10' not in m['symbols'] and 'SLC75A1' not in m['symbols'] for m in bank['meta_programs']))
add('At_least_3_independent_patients',all(len(set(m['support_patients']))>=3 for m in bank['meta_programs']))
add('Exactly_50_unique_genes_each',all(len(set(m['genes']))==50 for m in bank['meta_programs']))
fit=json.loads((out/'metadata/Discovery_patient_fit_list.json').read_text())['patients'];add('Thirteen_discovery_patients',len(fit)==13)
n=0;seeds=[]
for pid in fit:
    c=cNMF(output_dir=str(out/'discovery/cnmf'),name=pid);par=load_df_from_npz(c.paths['nmf_replicate_parameters'])
    inp=ad.read_h5ad(out/'cache/discovery_patient_inputs'/f'{pid}.h5ad',backed='r')
    add('Input_excludes_target_'+pid,'ENSG00000109736' not in inp.var_names and not inp.var.gene_symbol.isin(['MFSD10','SLC75A1']).any());inp.file.close()
    add('Complete_300_repeats_'+pid,len(par)==300 and set(par.n_components)=={4,6,8})
    for _,r in par.iterrows():
        p=c.paths['iter_spectra']%(int(r.n_components),int(r['iter']));z=load_df_from_npz(p)
        if z.shape!=(int(r.n_components),2000) or not np.isfinite(z.to_numpy()).all() or (z.to_numpy()<0).any() or 'ENSG00000109736' in z.columns:raise RuntimeError('Invalid primary spectrum '+str(p))
        n+=1;seeds.append(dict(patient_id=pid,K=int(r.n_components),iteration=int(r['iter']),seed=int(r.nmf_seed)))
add('All_3900_primary_spectra_valid',n==3900);pd.DataFrame(seeds).to_csv(out/'audit/Primary_NMF_seed_manifest.tsv',sep='\t',index=False)
sens=pd.read_csv(out/'discovery/Parameter_sensitivity_summary.tsv',sep='\t');add('Nineteen_sensitivity_scenarios',len(sens)==19 and sens.status.eq('complete').all())
raw=json.loads((out/'audit/GSE176078_import_audit.json').read_text());add('Validation_sparse_input_QC',all(raw[k] for k in ['matrix_nonnegative_integer','barcodes_unique','metadata_one_to_one','UMI_matches_author_metadata','nFeature_matches_author_metadata']))
add('Validation_author_QC_no_new_loss',raw['Cancer_Epithelial_before_QC']==raw['Cancer_Epithelial_after_QC']==24489 and raw['QC_inconsistent_all_cells']==0)
for name,shape in [('GSE176078_all_cells_author_counts_clinical_v2.h5ad',(100064,29733)),('GSE176078_malignant_candidates_author_counts_clinical_v2.h5ad',(24489,29733))]:
    obj=ad.read_h5ad(out/'objects'/name,backed='r');add('Object_shape_'+name,obj.shape==shape);add('Object_counts_sparse_'+name,getattr(obj.X,'format',None) in ['csr','csc']);add('Object_unique_barcodes_'+name,obj.obs_names.is_unique);obj.file.close()
patients=pd.read_csv(out/'validation/GSE176078_all26_patient_summary.tsv',sep='\t');add('All26_clinical_cases_preserved',len(patients)==26 and patients.patient_id.is_unique)
ext=json.loads((out/'audit/External_validation_summary.json').read_text());add('Primary16_naive_4_treated_separate',ext['primary_naive_patients']==16 and ext['treated_patients']==4)
add('No_false_confirmatory_success',ext['discovery_confirmatory_associations']==0 and ext['confirmed_supported_associations']==0)
add('External_gene_coverage',pd.read_csv(out/'validation/Program_gene_coverage.tsv',sep='\t').fraction_mapped.ge(.8).all())
for cohort in ['discovery','validation']:
    p=out/cohort/'Coherence_primary_only_reference_sensitivity.tsv';d=pd.read_csv(p,sep='\t');add('Primary_reference_coherence_sensitivity_'+cohort,len(d)==5 and (d.median_excess>0).all() and d.BH_q.lt(.05).all())

# Validate earlier delivered scientific files against their existing manifests, without rewriting them.
protection=[]
stage=root/'outputs/第一阶段_2026-10-23';base=pd.read_csv(stage/'audit/Final_artifact_SHA256_manifest.tsv',sep='\t')
for _,r in base.iterrows():
    if not str(r.relative_path).startswith('outputs'):continue
    p=root/str(r.relative_path);actual=sha(p);protection.append(dict(stage='stage1',path=str(p),unchanged=actual.lower()==str(r.sha256).lower()))
prev=root/'outputs/MFSD10_表达分布摸底_2026-10-26';base=pd.read_csv(prev/'audit/Deliverable_SHA256.tsv',sep='\t')
for _,r in base.iterrows():
    p=prev/str(r.relative_path);protection.append(dict(stage='expression_survey',path=str(p),unchanged=sha(p).lower()==str(r.sha256).lower()))
add('Earlier_results_preserved',len(protection)==140 and all(r['unchanged'] for r in protection),'86 stage1 + 54 expression-survey artifacts')
pd.DataFrame(protection).to_csv(out/'audit/Earlier_deliverables_preservation.tsv',sep='\t',index=False)

sources=[]
for p in (root/'原始数据/单细胞/GSE176078').iterdir():
    if p.is_file():sources.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
pd.DataFrame(sources).to_csv(out/'audit/Original_GSE176078_input_SHA256.tsv',sep='\t',index=False)
add('Three_original_validation_inputs_present',len(sources)==3)
local=[]
for p in out.glob('*.md'):
    for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
        if link.startswith(('http','codex:')):continue
        local.append(dict(document=p.name,target=link,exists=(p.parent/link).exists()))
add('Report_local_links',all(r['exists'] for r in local));pd.DataFrame(local).to_csv(out/'audit/Report_local_links.tsv',sep='\t',index=False)
add('Seven_figure_pairs',len(list((out/'figures').glob('*.png')))==7 and len(list((out/'figures').glob('*.svg')))==7)
for p in (out/'scripts').glob('*.py'):compile(p.read_text(encoding='utf-8'),str(p),'exec')
add('Python_scripts_parse',True)
versions={p:im.version(p) for p in ['cnmf','numpy','pandas','scipy','scikit-learn','anndata','scanpy','statsmodels','matplotlib','h5py']};(out/'logs/Software_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
pd.DataFrame(checks).rename(columns={'pass_':'pass'}).to_csv(out/'audit/Final_validation.tsv',sep='\t',index=False)
progress=json.loads((out/'audit/PROGRESS.json').read_text());progress.update(status='complete',stage='complete',final_audit_complete=True,background_task_remaining=False);(out/'audit/PROGRESS.json').write_text(json.dumps(progress,indent=2),encoding='utf-8')
manifest=[]
for p in out.rglob('*'):
    if not p.is_file() or any(part in ['env','cache','__pycache__'] for part in p.relative_to(out).parts):continue
    if p.name in ['Deliverable_SHA256.tsv','PIPELINE_COMPLETED.json','09_final_audit.log']:continue # The current redirected stdout is still being written.
    manifest.append(dict(relative_path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
pd.DataFrame(manifest).to_csv(out/'audit/Deliverable_SHA256.tsv',sep='\t',index=False)
done=dict(status='complete',date='2026-10-26',validation_checks=len(checks),all_passed=True,primary_NMF_runs=3900,sensitivity_scenarios=19,sensitivity_initial_repeats=380,frozen_meta_programs=5,discovery_patients=13,validation_primary_naive_patients=16,validation_treated_patients_with_candidates=4,discovery_significant_associations=0,confirmatory_successes=0,external_exploratory_MP05_within_rho=.057081,interpretation='Shared-program gene coherence supported; no discovery-FDR-supported MFSD10 association was confirmed. External MP05 is exploratory only.',background_task_remaining=False)
(out/'audit/PIPELINE_COMPLETED.json').write_text(json.dumps(done,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(done,ensure_ascii=False,indent=2))
