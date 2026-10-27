"""Freeze independent drug sets first; only then read and intersect disease candidates."""
import hashlib
import json
import math
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=pathlib.Path(__file__).resolve().parents[1]
PRIOR=ROOT.parent/'DeepPurpose_预测_2026-10-24'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')
def main():
    assert (ROOT/'predictions'/'INDEPENDENT_RAW_PREDICTIONS_COMPLETE.txt').exists()
    out=ROOT/'results';out.mkdir(parents=True,exist_ok=True)
    drug=out/'independent_drug_sets';drug.mkdir(exist_ok=True)
    matches=out/'disease_intersections';matches.mkdir(exist_ok=True)
    figs=ROOT/'figures';figs.mkdir(exist_ok=True)
    background=pd.read_csv(ROOT/'background'/'All_reviewed_human_protein_background.tsv',sep='\t',keep_default_na=False)
    mapping=pd.read_csv(ROOT/'background'/'Protein_to_full_reference_gene_mapping.tsv',sep='\t').drop_duplicates(['uniprot_accession','ensembl_stable'])
    runs=json.loads((ROOT/'predictions'/'All_model_run_summary.json').read_text())
    raw=[];tables={};selected={};stability=[];drug_manifest=[]
    for endpoint in ['Kd','IC50']:
        matrix=pd.DataFrame(index=background.uniprot_accession)
        for run in [r for r in runs if r['endpoint']==endpoint]:
            source=ROOT/'predictions'/(run['model']+'_scores.tsv');assert sha(source)==run['score_sha256']
            scores=pd.read_csv(source,sep='\t')
            matrix[run['model']]=scores.set_index('uniprot_accession').predicted_p_value
            scores['model']=run['model'];scores['endpoint']=endpoint;raw.append(scores)
        assert matrix.shape==(20431,4) and np.isfinite(matrix.to_numpy()).all()
        percentiles=matrix.rank(method='average',pct=True);consensus=percentiles.median(axis=1)
        table=background.drop(columns='protein_sequence').set_index('uniprot_accession')
        table['consensus_percentile']=consensus;table['consensus_rank']=consensus.rank(ascending=False,method='min').astype(int)
        table['model_percentile_min']=percentiles.min(axis=1);table['model_percentile_max']=percentiles.max(axis=1)
        table['model_percentile_IQR']=percentiles.quantile(.75,axis=1)-percentiles.quantile(.25,axis=1)
        for model in matrix.columns:
            table[model+'_predicted_p_value']=matrix[model];table[model+'_percentile']=percentiles[model]
        table=table.reset_index().sort_values(['consensus_rank','uniprot_accession']);tables[endpoint]=table
        table.to_csv(out/(endpoint+'_independent_proteome_ranked.tsv'),sep='\t',index=False)
        for threshold in [1,5,10]:
            limit=math.ceil(len(table)*threshold/100)
            proteins=table.loc[table.consensus_rank.le(limit)].copy()
            accessions=set(proteins.uniprot_accession)
            genes=set(mapping.loc[mapping.uniprot_accession.isin(accessions),'ensembl_stable'])
            selected[(endpoint,threshold)]=(accessions,genes)
            name=f'{endpoint}_top{threshold}percent'
            proteins.to_csv(drug/(name+'_proteins.tsv'),sep='\t',index=False)
            gene_details=mapping[mapping.uniprot_accession.isin(accessions)].merge(table[['uniprot_accession','consensus_rank','consensus_percentile']],on='uniprot_accession')
            gene_details.to_csv(drug/(name+'_gene_mappings.tsv'),sep='\t',index=False)
            drug_manifest.append({'endpoint':endpoint,'threshold_percent':threshold,'nominal_protein_cutoff':limit,'n_selected_proteins':len(proteins),'n_mapped_genes':len(genes),'boundary_ties_included':True,'protein_set_sha256':sha(drug/(name+'_proteins.tsv')),'gene_mapping_sha256':sha(drug/(name+'_gene_mappings.tsv'))})
            for dropped in matrix.columns:
                reduced=percentiles.drop(columns=dropped).median(axis=1)
                reduced_rank=reduced.rank(ascending=False,method='min')
                alternative=set(reduced_rank[reduced_rank.le(limit)].index)
                stability.append({'endpoint':endpoint,'threshold_percent':threshold,'dropped_model':dropped,'base_n':len(accessions),'alternative_n':len(alternative),'shared_proteins':len(accessions&alternative),'Jaccard':len(accessions&alternative)/len(accessions|alternative)})
    # These independent drug outputs are finished before any disease input is read.
    save_json(drug/'Drug_side_sets_frozen_before_disease_intersection.json',drug_manifest)
    pd.concat(raw,ignore_index=True).to_csv(out/'All_independent_protein_model_scores.tsv',sep='\t',index=False)
    pd.DataFrame(stability).to_csv(out/'Leave_one_encoding_out_threshold_stability.tsv',sep='\t',index=False)
    diseasepath=PRIOR/'main_input'/'Disease_matching_candidates_1881.tsv'
    disease=pd.read_csv(diseasepath,sep='\t');disease_ids=set(disease.canonical_id)
    assert len(disease)==1881 and disease.canonical_id.is_unique
    extra=pd.read_csv(PRIOR/'main_input'/'Additional_diagnostic_only_76genes.tsv',sep='\t')
    combined=disease.copy();covered=set(mapping.ensembl_stable)
    combined['in_reviewed_human_background']=combined.canonical_id.isin(covered)
    summary=[]
    for endpoint,table in tables.items():
        linked=mapping.merge(table[['uniprot_accession','consensus_rank','consensus_percentile','CNN_truncates_after1000','nonstandard_amino_acids']],on='uniprot_accession')
        best=linked.groupby('ensembl_stable').agg(best_protein_rank=('consensus_rank','min'),best_protein_consensus_percentile=('consensus_percentile','max'),reviewed_proteins=('uniprot_accession',lambda s:';'.join(sorted(set(s)))),any_CNN_truncation=('CNN_truncates_after1000','any'),nonstandard_amino_acids=('nonstandard_amino_acids',lambda s:';'.join(sorted(set(s)-{''}))))
        combined=combined.merge(best.add_prefix(endpoint+'_'),how='left',left_on='canonical_id',right_index=True,validate='one_to_one')
    for threshold in [1,5,10]:
        kdprot,kd=selected[('Kd',threshold)];icprot,ic=selected[('IC50',threshold)]
        sets={'Kd':kd,'IC50':ic,'Kd_union_IC50':kd|ic,'Kd_intersection_IC50':kd&ic}
        for label,genes in sets.items():
            col=f'top{threshold}percent_{label}'
            combined[col]=combined.canonical_id.isin(genes)
            gene_list=pd.DataFrame({'canonical_id':sorted(genes)})
            gene_list.to_csv(drug/(col+'_genes.tsv'),sep='\t',index=False)
            matched=combined[combined[col]].copy()
            matched.sort_values(['Kd_best_protein_rank','IC50_best_protein_rank','canonical_id'],na_position='last').to_csv(matches/(col+'_disease1881.tsv'),sep='\t',index=False)
            ex=extra[extra.canonical_id.isin(genes)];ex.to_csv(matches/(col+'_diagnostic_extra76.tsv'),sep='\t',index=False)
            summary.append({'threshold_percent':threshold,'drug_set':label,'drug_side_genes':len(genes),'disease1881_intersection':len(disease_ids&genes),'diagnostic_extra76_intersection':len(ex),'Kd_selected_proteins':len(kdprot),'IC50_selected_proteins':len(icprot),'Kd_IC50_shared_proteins':len(kdprot&icprot)})
    combined.to_csv(out/'All1881_disease_genes_independent_prediction_evidence.tsv',sep='\t',index=False,na_rep='NA')
    combined[~combined.in_reviewed_human_background].to_csv(out/'Disease_genes_outside_reviewed_background.tsv',sep='\t',index=False,na_rep='NA')
    summary=pd.DataFrame(summary);summary.to_csv(out/'Threshold_comparison_1_5_10percent.tsv',sep='\t',index=False)
    cross=tables['Kd'].set_index('uniprot_accession').consensus_percentile.to_frame('Kd').join(tables['IC50'].set_index('uniprot_accession').consensus_percentile.rename('IC50'))
    stats={'background_proteins':len(background),'unique_sequences':int(background.sequence_sha256.nunique()),'raw_scores':sum(r['n_scores'] for r in runs),'disease_main_n':1881,'disease_covered':len(disease_ids&covered),'disease_outside_reviewed_background':len(disease_ids-covered),'mapped_full_reference_genes':len(covered),'background_proteins_CNN_truncated':int(background.CNN_truncates_after1000.sum()),'background_proteins_nonstandard_amino_acids':int(background.nonstandard_amino_acids.ne('').sum()),'Kd_IC50_consensus_Spearman':cross.Kd.corr(cross.IC50,method='spearman'),'disease_input_sha256':sha(diseasepath),'single_primary_threshold':None,'selection_unit':'protein accessions; map to gene sets before gene-level endpoint union/intersection; disease intersection last'}
    save_json(out/'Summary_statistics.json',stats)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(9,5.2))
    labels=['Kd','IC50','Kd_union_IC50','Kd_intersection_IC50'];colors=['#4477AA','#EE6677','#228833','#AA3377']
    xpos=np.arange(3);width=.19
    for j,(label,color) in enumerate(zip(labels,colors)):
        values=summary.loc[summary.drug_set.eq(label),'disease1881_intersection'].to_numpy()
        bars=ax.bar(xpos+(j-1.5)*width,values,width,color=color,label={'Kd_union_IC50':'Union','Kd_intersection_IC50':'Intersection'}.get(label,label))
        ax.bar_label(bars,padding=3,fontsize=9)
    ax.set_xticks(xpos,['Top 1%','Top 5%','Top 10%']);ax.set_ylabel('Genes intersecting disease set (n=1,881)');ax.set_ylim(0,summary.disease1881_intersection.max()*1.22+1)
    ax.legend(ncol=4,loc='upper left');fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(figs/('Disease_intersections_by_threshold.'+ext),dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.3,5.5));density=ax.hexbin(cross.Kd,cross.IC50,gridsize=45,mincnt=1,cmap='viridis',bins='log')
    ax.set_xlabel('Kd consensus percentile');ax.set_ylabel('IC50 consensus percentile');ax.set_title(f"Independent human proteome (n=20,431)\nSpearman = {stats['Kd_IC50_consensus_Spearman']:.3f}")
    fig.colorbar(density,ax=ax,label='Protein count (log color scale)');fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(figs/('Independent_endpoint_agreement.'+ext),dpi=180)
    plt.close(fig)
    print(json.dumps(stats,ensure_ascii=True));print(summary.to_string(index=False))
if __name__=='__main__':main()
