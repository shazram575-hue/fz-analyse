"""Independent recomputation of all ranks, cutoffs, sets, input preservation and hashes."""
import hashlib
import json
import math
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=pathlib.Path(__file__).resolve().parents[1]
PRIOR=ROOT.parent/'DeepPurpose_预测_2026-10-24'
WORKSPACE=ROOT.parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    out=ROOT/'audit';out.mkdir(exist_ok=True);checks=[]
    def check(name,value,detail=''):checks.append({'check':name,'pass':bool(value),'detail':str(detail)})
    protocol=json.loads((ROOT/'background'/'Frozen_independent_drug_side_protocol.json').read_text(encoding='utf-8'))
    for name,key in [('All_reviewed_human_protein_background.tsv','background_sha256'),('Protein_to_full_reference_gene_mapping.tsv','mapping_sha256')]:check(name+'_frozen',sha(ROOT/'background'/name)==protocol[key])
    identity=json.loads((ROOT/'structure_comparison'/'Published_A_Figure1_reconstruction.json').read_text())
    check('published_A_CID_all9_centers_and_full_key',identity['all9_centers_match'] and identity['complete_InChIKey_match'])
    bg=pd.read_csv(ROOT/'background'/'All_reviewed_human_protein_background.tsv',sep='\t',keep_default_na=False)
    check('20431_independent_unique_accessions',len(bg)==20431 and bg.uniprot_accession.is_unique)
    check('all_protein_sequence_hashes',all(hashlib.sha256(s.encode()).hexdigest()==h for s,h in zip(bg.protein_sequence,bg.sequence_sha256)))
    raw=pd.read_csv(ROOT/'results'/'All_independent_protein_model_scores.tsv',sep='\t')
    check('163448_finite_unique_scores',len(raw)==163448 and np.isfinite(raw.predicted_p_value).all() and not raw.duplicated(['model','uniprot_accession']).any())
    check('eight_models_no_MPNN',raw.model.nunique()==8 and not raw.model.str.contains('MPNN').any())
    runs=json.loads((ROOT/'predictions'/'All_model_run_summary.json').read_text())
    check('all_model_weights_and_score_files_hashes',all(sha(PRIOR/'models'/r['model']/'model.pt')==r['model_sha256'] and sha(ROOT/'predictions'/(r['model']+'_scores.tsv'))==r['score_sha256'] for r in runs))
    check('native_optimized_forward_checks_all8',len(runs)==8 and all(r['native_validation_n']==16 and r['native_vs_cached_drug_forward_max_abs_difference']<1e-4 for r in runs))
    # Compare shared sequences against the previously completed native prediction pipeline.
    comparisons=[]
    for run in runs:
        old=pd.read_csv(PRIOR/'predictions'/(run['model']+'_protein_scores.tsv'),sep='\t').groupby('sequence_sha256').predicted_p_value.median()
        current=raw[raw.model.eq(run['model'])].set_index('sequence_sha256').predicted_p_value
        shared=current.to_frame('current').join(old.rename('prior_native'),how='inner')
        diff=float((shared.current-shared.prior_native).abs().max())
        comparisons.append({'model':run['model'],'shared_protein_entries':len(shared),'shared_unique_sequences':shared.index.nunique(),'max_absolute_difference':diff})
    pd.DataFrame(comparisons).to_csv(out/'Prior_native_prediction_equivalence.tsv',sep='\t',index=False)
    check('shared_sequences_agree_with_prior_native_inference',all(r['shared_unique_sequences']>1800 and r['max_absolute_difference']<1e-4 for r in comparisons),comparisons)
    mapping=pd.read_csv(ROOT/'background'/'Protein_to_full_reference_gene_mapping.tsv',sep='\t')
    sets={};cutoffchecks=[]
    for endpoint in ['Kd','IC50']:
        matrix=raw[raw.endpoint.eq(endpoint)].pivot(index='uniprot_accession',columns='model',values='predicted_p_value').sort_index()
        pct=rankdata(matrix.to_numpy(),axis=0,method='average')/len(matrix)
        median=np.median(pct,axis=1);ranks=rankdata(-median,method='min').astype(int)
        reported=pd.read_csv(ROOT/'results'/(endpoint+'_independent_proteome_ranked.tsv'),sep='\t').set_index('uniprot_accession').reindex(matrix.index)
        check(endpoint+'_consensus_independent_recalculation',matrix.shape==(20431,4) and np.allclose(median,reported.consensus_percentile,atol=1e-12,rtol=0) and np.array_equal(ranks,reported.consensus_rank))
        for p in [1,5,10]:
            expected=set(matrix.index[ranks<=math.ceil(len(matrix)*p/100)])
            gene_expected=set(mapping.loc[mapping.uniprot_accession.isin(expected),'ensembl_stable'])
            selected=pd.read_csv(ROOT/'results'/'independent_drug_sets'/f'{endpoint}_top{p}percent_proteins.tsv',sep='\t')
            check(f'{endpoint}_top{p}_tie_aware_independent_cutoff',set(selected.uniprot_accession)==expected)
            sets[(endpoint,p)]=gene_expected
    disease=pd.read_csv(PRIOR/'main_input'/'Disease_matching_candidates_1881.tsv',sep='\t');diseaseids=set(disease.canonical_id)
    summary=pd.read_csv(ROOT/'results'/'Threshold_comparison_1_5_10percent.tsv',sep='\t')
    for p in [1,5,10]:
        kd=sets[('Kd',p)];ic=sets[('IC50',p)]
        for label,genes in {'Kd':kd,'IC50':ic,'Kd_union_IC50':kd|ic,'Kd_intersection_IC50':kd&ic}.items():
            name=f'top{p}percent_{label}'
            file=ROOT/'results'/'disease_intersections'/(name+'_disease1881.tsv')
            reported=set(pd.read_csv(file,sep='\t').canonical_id)
            druggenes=set(pd.read_csv(ROOT/'results'/'independent_drug_sets'/(name+'_genes.tsv'),sep='\t').canonical_id)
            count=summary[(summary.threshold_percent==p)&summary.drug_set.eq(label)].iloc[0]
            check(name+'_gene_sets_and_disease_intersection',druggenes==genes and reported==genes&diseaseids and len(reported)==count.disease1881_intersection and len(genes)==count.drug_side_genes)
    allgenes=pd.read_csv(ROOT/'results'/'All1881_disease_genes_independent_prediction_evidence.tsv',sep='\t')
    check('1881_genes_retained_with_background_coverage_flag',len(allgenes)==1881 and set(allgenes.canonical_id)==diseaseids and allgenes.in_reviewed_human_background.sum()==len(diseaseids&set(mapping.ensembl_stable)))
    scope=json.loads((PRIOR/'main_input'/'Confirmed_scope_and_provenance.json').read_text(encoding='utf-8'))
    check('all_confirmed_disease_input_hashes_unchanged',all(sha(PRIOR/'main_input'/name)==expected for name,expected in scope['files_sha256'].items()))
    frozen=pd.read_csv(WORKSPACE/'outputs'/'单细胞基础补全_2026-10-23'/'audit'/'Frozen_stage1_before.tsv',sep='\t')
    preserved=[]
    for row in frozen.itertuples(index=False):
        actual=sha(WORKSPACE/'outputs'/'第一阶段_2026-10-23'/row.relative_path)
        preserved.append({'relative_path':row.relative_path,'pass':actual.lower()==row.sha256.lower(),'sha256':actual})
    pd.DataFrame(preserved).to_csv(out/'Frozen_stage1_preservation.tsv',sep='\t',index=False)
    check('original7_phase1_files_unchanged',len(preserved)==7 and all(r['pass'] for r in preserved))
    pd.DataFrame(checks).to_csv(out/'Independent_result_validation.tsv',sep='\t',index=False)
    report={'n_checks':len(checks),'n_pass':sum(r['pass'] for r in checks),'failed':[r for r in checks if not r['pass']]}
    (out/'Independent_validation_summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=True));assert not report['failed']
if __name__=='__main__':main()
