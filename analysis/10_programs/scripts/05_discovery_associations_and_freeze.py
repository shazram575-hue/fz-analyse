from pathlib import Path
import sys,json,hashlib,shutil
import numpy as np,pandas as pd,anndata as ad
from program_analysis import prepare,associations,coherence

root=Path(sys.argv[1]);out=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26'
assert not (out/'frozen/Analysis_freeze_v1.json').exists(),'Already frozen; do not overwrite discovery estimates'
bank=json.loads((out/'frozen/Programs_v1.json').read_text());fit=json.loads((out/'metadata/Discovery_patient_fit_list.json').read_text())['patients']
obj=ad.read_h5ad(out/'objects/GSE161529_confirmed_malignant_singlets_counts.h5ad')
obs,scores,coverage,lookup,denom=prepare(obj,bank,'discovery')
coverage.to_csv(out/'discovery/Program_gene_coverage.tsv',sep='\t',index=False)
scores.to_csv(out/'discovery/Fixed_program_scores_all_singlets.tsv.gz',sep='\t')
result,per,means=associations(obs,scores,fit)
result.to_csv(out/'discovery/MFSD10_patient_level_associations.tsv',sep='\t',index=False)
per.to_csv(out/'discovery/MFSD10_within_patient_effects.tsv',sep='\t',index=False)
means.to_csv(out/'discovery/MFSD10_between_patient_inputs.tsv',sep='\t',index=False)
cp,cs=coherence(obj,obs,bank,lookup,denom,fit)
cp.to_csv(out/'discovery/Gene_coherence_by_patient.tsv',sep='\t',index=False);cs.to_csv(out/'discovery/Gene_coherence_summary.tsv',sep='\t',index=False)
directions=[]
for r in result.to_dict('records'):
    directions.append(dict(program_id=r['program_id'],association=r['association'],discovery_effect=r['effect_rho'],direction=int(np.sign(r['effect_rho'])) if np.isfinite(r['effect_rho']) else 0,discovery_p=r['p'],discovery_q=r['BH_q'],confirmatory=bool(np.isfinite(r['BH_q']) and r['BH_q']<.05)))
frozen={'version':'v1','discovery_cohort':'GSE161529','external_cohort':'GSE176078','external_primary_treatment':'Naive','external_other_treatments':'Separate descriptive/secondary report','external_patient_min_cells':50,'min_patients_within_inference':6,'min_patients_between_inference':12,'min_between_residual_df':5,'min_program_gene_coverage':.8,'target_not_used_for_program_selection':True,'scoring':'Fixed 50 genes mean log1p CP10K with MFSD10 excluded from denominator; no cohort joint normalization','within_model':'Per-patient rank-residual rho adjusting target-excluded logUMI/logFeatures/mitochondrial percentage; equal-patient Fisher-z; patient sign-flips and bootstrap CI','between_model':'Patient-level partial rank rho adjusting explicit clinical IHC, median target-excluded UMI and mitochondrial percentage; t approximation and Fisher CI; descriptive if minimum sample/df not met','discovery_multiple_testing':'BH jointly across MP x {within,between}','external_confirmatory_multiple_testing':'Prespecified one-sided direction, BH across discovery-supported association list; require same effect direction, q<0.05 and independent coherence support','external_all_programs_reporting':'All frozen programs two-sided, joint BH, regardless of discovery significance','coherence':'Depth/feature/mitochondrial residual gene-gene mean offdiagonal correlation, 25 expression-matched null gene sets, patient-level one-sided signed-rank test with BH; coverage>=80%','directions':directions}
shutil.copyfile(out/'scripts/program_analysis.py',out/'frozen/Analysis_code_v1.py')
hashes={}
for p in [out/'frozen/Programs_v1.json',out/'frozen/Program_genes_v1.tsv',out/'frozen/Analysis_code_v1.py',out/'scripts/program_analysis.py',out/'预冻结分析规则_v1.md',out/'discovery/MFSD10_patient_level_associations.tsv']:
    hashes[p.relative_to(out).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
frozen['SHA256']=hashes
(out/'frozen/Analysis_freeze_v1.json').write_text(json.dumps(frozen,indent=2),encoding='utf-8')
p=out/'audit/PROGRESS.json';d=json.loads(p.read_text());d.update(association_rules_frozen=True,stage='external_validation_ready',n_discovery_supported_associations=sum(a['confirmatory'] for a in directions));p.write_text(json.dumps(d,indent=2),encoding='utf-8')
print(result.to_string(index=False));print('ANALYSIS_AND_DIRECTIONS_FROZEN; EXTERNAL_ASSOCIATION_NOT_YET_READ')
