from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd,anndata as ad
from program_analysis import prepare,associations,coherence,bh
root=Path(sys.argv[1]);out=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26'
freeze=json.loads((out/'frozen/Analysis_freeze_v1.json').read_text())
for name,digest in freeze['SHA256'].items():assert hashlib.sha256((out/name).read_bytes()).hexdigest()==digest,f'Frozen discovery file changed: {name}'
bank=json.loads((out/'frozen/Programs_v1.json').read_text());obj=ad.read_h5ad(out/'objects/GSE176078_malignant_candidates_author_counts_clinical_v2.h5ad')
obs,scores,coverage,lookup,denom=prepare(obj,bank,'validation')
scores.to_csv(out/'validation/Fixed_program_scores_all_author_candidates.tsv.gz',sep='\t');coverage.to_csv(out/'validation/Program_gene_coverage.tsv',sep='\t',index=False)
counts=obs.groupby('patient_id',observed=True).size()
primary=obs.loc[obs.treatment_status=='Naive','patient_id'].unique().tolist();primary=[p for p in primary if counts[p]>=50]
treated=obs.loc[obs.treatment_status=='Treated','patient_id'].unique().tolist()
results,per,means=associations(obs,scores,primary)
results['analysis_role']='Primary_Naive';results.to_csv(out/'validation/MFSD10_external_all_frozen_programs.tsv',sep='\t',index=False)
per.to_csv(out/'validation/MFSD10_within_patient_effects_Naive.tsv',sep='\t',index=False);means.to_csv(out/'validation/MFSD10_between_patient_inputs_Naive.tsv',sep='\t',index=False)
cp,cs=coherence(obj,obs,bank,lookup,denom,primary)
cp.to_csv(out/'validation/Gene_coherence_by_patient_Naive.tsv',sep='\t',index=False);cs.to_csv(out/'validation/Gene_coherence_summary_Naive.tsv',sep='\t',index=False)
directions={(a['program_id'],a['association']):a['direction'] for a in freeze['directions'] if a['confirmatory']}
if directions:
    confirm,_,_=associations(obs,scores,primary,directions=directions)
    confirm=confirm[[tuple(r) in directions for r in confirm[['program_id','association']].to_numpy()]].copy()
    confirm['BH_q_confirmatory']=bh(confirm.p)
else:confirm=pd.DataFrame(columns=['program_id','association','effect_rho','p','BH_q_confirmatory'])
compare=results.merge(pd.DataFrame(freeze['directions']),on=['program_id','association'],how='left').merge(cs[['program_id','median_excess','BH_q']].rename(columns={'BH_q':'coherence_q'}),on='program_id',how='left')
if len(confirm):compare=compare.merge(confirm[['program_id','association','p','BH_q_confirmatory']].rename(columns={'p':'directional_p'}),on=['program_id','association'],how='left')
else:compare['directional_p']=np.nan;compare['BH_q_confirmatory']=np.nan
compare['same_direction']=np.sign(compare.effect_rho)==compare.direction
compare['coherence_supported']=(compare.median_excess>0)&(compare.coherence_q<.05)
compare['confirmatory_supported']=compare.confirmatory&compare.same_direction&(compare.BH_q_confirmatory<.05)&compare.coherence_supported
compare['interpretation']=np.select([compare.confirmatory_supported,compare.confirmatory&~compare.coherence_supported,compare.confirmatory&~(compare.BH_q_confirmatory<.05),~compare.confirmatory],['Supported_direction_effect_and_coherence','Not_supported_program_coherence','Not_supported_association','Exploratory_only_no_discovery_significant_association'],default='Not_supported')
compare.to_csv(out/'validation/Discovery_external_comparison_all.tsv',sep='\t',index=False);confirm.to_csv(out/'validation/Confirmatory_directional_tests.tsv',sep='\t',index=False)
if treated:
    tr,tp,tm=associations(obs,scores,treated);tr['analysis_role']='Treated_secondary';tr.to_csv(out/'validation/MFSD10_Treated_separate.tsv',sep='\t',index=False);tp.to_csv(out/'validation/MFSD10_Treated_patient_effects.tsv',sep='\t',index=False)
    tm.to_csv(out/'validation/MFSD10_Treated_patient_inputs.tsv',sep='\t',index=False)
# Lineage-mixture flag is a sensitivity only, without changing original author classification or primary object.
sens_obs=obs.loc[~obs.lineage_mixture_flag.astype(bool)].copy();sens_scores=scores.loc[sens_obs.index]
sens,_,_=associations(sens_obs,sens_scores,primary);sens.to_csv(out/'validation/MFSD10_lineage_flag_excluded_sensitivity.tsv',sep='\t',index=False)
(out/'audit/External_validation_summary.json').write_text(json.dumps({'primary_naive_patients':len(primary),'primary_naive_patient_ids':primary,'treated_patients':len(treated),'frozen_programs':len(bank['meta_programs']),'discovery_confirmatory_associations':len(directions),'confirmed_supported_associations':int(compare.confirmatory_supported.sum()),'programs_with_external_coherence_support':int(((cs.median_excess>0)&(cs.BH_q<.05)).sum()),'frozen_hashes_verified':True},indent=2),encoding='utf-8')
p=out/'audit/PROGRESS.json';d=json.loads(p.read_text());d.update(external_validation_complete=True,stage='diagnostics_and_final_report');p.write_text(json.dumps(d,indent=2),encoding='utf-8')
print(compare[['program_id','association','discovery_effect','effect_rho','BH_q','BH_q_confirmatory','coherence_supported','interpretation']].to_string(index=False));print('EXTERNAL_VALIDATION_COMPLETE')
