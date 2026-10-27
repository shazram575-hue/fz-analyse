"""Summarize all scores without hard target selection; endpoint-specific rank consensus."""
import json
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=pathlib.Path(__file__).resolve().parents[1]

def rank_consensus(gene_scores, endpoint):
    selected=gene_scores[gene_scores.endpoint.eq(endpoint)].copy()
    primary=selected[selected.role.eq("primary_encoding_ensemble")]
    wide=primary.pivot(index="canonical_id",columns="model",values="predicted_log_affinity")
    assert wide.shape[1]==4 and not wide.isna().any().any()
    percentiles=wide.rank(pct=True,method="average")
    consensus=percentiles.median(axis=1)
    table=pd.DataFrame({"canonical_id":wide.index,"consensus_percentile":consensus,
                        "consensus_rank":consensus.rank(ascending=False,method="min").astype(int),
                        "model_percentile_IQR":percentiles.quantile(.75,axis=1)-percentiles.quantile(.25,axis=1),
                        "model_percentile_min":percentiles.min(axis=1),"model_percentile_max":percentiles.max(axis=1),
                        "median_predicted_log_affinity_across4models":wide.median(axis=1),
                        "n_primary_models":4}).reset_index(drop=True)
    shifts=[];stability=[]
    base_ranks=consensus.rank(ascending=False,method="min")
    base_order=consensus.sort_values(ascending=False,kind="stable").index
    for dropped in wide.columns:
        reduced=percentiles.drop(columns=dropped).median(axis=1)
        reduced_ranks=reduced.rank(ascending=False,method="min")
        shifts.append((reduced_ranks-base_ranks).abs())
        order=reduced.sort_values(ascending=False,kind="stable").index
        stability.append({"endpoint":endpoint,"dropped_model":dropped,"n_genes":len(wide),
                          "top20_overlap":len(set(base_order[:20]) & set(order[:20])),
                          "top50_overlap":len(set(base_order[:50]) & set(order[:50])),
                          "rank_Spearman":base_ranks.corr(reduced_ranks,method="spearman"),
                          "role":"encoding_omission_sensitivity_not_independent_validation"})
    shift=pd.concat(shifts,axis=1).max(axis=1)
    table["leave_one_model_out_max_abs_rank_shift"]=table.canonical_id.map(shift)
    for model in wide.columns:
        table[model+"_predicted_log_affinity"]=table.canonical_id.map(wide[model])
        table[model+"_percentile"]=table.canonical_id.map(percentiles[model])
    sensitivity=selected[selected.role.eq("sensitivity_MPNN")].set_index("canonical_id")
    table["MPNN_sensitivity_predicted_log_affinity"]=table.canonical_id.map(sensitivity.predicted_log_affinity)
    table["MPNN_role"]="single_pair_sensitivity_only_original_batch_dependency"
    table["endpoint"]=endpoint
    return table,percentiles,stability

def main():
    source=ROOT/"predictions"
    assert (source/"RAW_PREDICTIONS_COMPLETE.txt").exists()
    results=ROOT/"results";results.mkdir(parents=True,exist_ok=True)
    figures=ROOT/"figures";figures.mkdir(parents=True,exist_ok=True)
    catalog=pd.read_csv(source/"Protein_sequence_catalog.tsv",sep="\t",keep_default_na=False)
    panel=json.loads((source/"All_model_run_summary.json").read_text(encoding="utf-8"))
    raw=pd.concat([pd.read_csv(source/(item["model"]+"_protein_scores.tsv"),sep="\t") for item in panel],ignore_index=True)
    raw=raw.rename(columns={"predicted_p_value":"predicted_log_affinity"})
    assert len(raw)==19930 and raw.predicted_log_affinity.notna().all()
    assert not raw.duplicated(["model","prediction_row_id"]).any()
    raw.to_csv(results/"All_protein_model_predictions.tsv",sep="\t",index=False)
    groups=["input_scope","canonical_id","gene_symbol","model","endpoint","role"]
    scores=raw.groupby(groups,sort=True).agg(predicted_log_affinity=("predicted_log_affinity","median"),
                                           protein_score_min=("predicted_log_affinity","min"),
                                           protein_score_max=("predicted_log_affinity","max"),
                                           n_unique_protein_sequences=("sequence_sha256","nunique")).reset_index()
    scores["protein_score_span"]=scores.protein_score_max-scores.protein_score_min
    scores.to_csv(results/"All_gene_model_scores.tsv",sep="\t",index=False)
    main_genes=pd.read_csv(ROOT/"main_input"/"Disease_matching_candidates_1881.tsv",sep="\t")
    extra=pd.read_csv(ROOT/"main_input"/"Additional_diagnostic_only_76genes.tsv",sep="\t")
    gene_flags=catalog.groupby("canonical_id").agg(n_input_protein_sequences=("sequence_sha256","nunique"),
                                                   max_sequence_length=("sequence_length","max"),
                                                   any_CNN_truncation=("CNN_truncates_after1000",lambda x:x.astype(str).str.lower().eq("true").any()),
                                                   any_GENCODE_translation_only=("mapping_tier",lambda x:x.eq("GENCODEv36_translation_only").any()),
                                                   any_nonstandard_amino_acid=("nonstandard_amino_acids",lambda x:x.ne("").any()),
                                                   mapping_methods=("mapping_method",lambda x:";".join(sorted(set(x))))).reset_index()
    main_genes=main_genes.merge(gene_flags,on="canonical_id",validate="one_to_one")
    main_scores=scores[scores.input_scope.eq("main")]
    all_tables={};all_percentiles={};stability=[]
    combined=main_genes.copy()
    for endpoint in ["Kd","IC50"]:
        table,percentiles,checks=rank_consensus(main_scores,endpoint)
        table=main_genes.merge(table,on="canonical_id",validate="one_to_one").sort_values(["consensus_rank","canonical_id"],kind="stable")
        table["interpretation"]="computational_prioritization_not_verified_binding_or_mechanism"
        table.to_csv(results/("Main_"+endpoint+"_consensus_ranked.tsv"),sep="\t",index=False,na_rep="NA")
        columns=["canonical_id","consensus_rank","consensus_percentile","model_percentile_IQR","leave_one_model_out_max_abs_rank_shift","median_predicted_log_affinity_across4models","MPNN_sensitivity_predicted_log_affinity"]
        combined=combined.merge(table[columns].rename(columns={c:endpoint+"_"+c for c in columns if c!="canonical_id"}),on="canonical_id",validate="one_to_one")
        all_tables[endpoint]=table;all_percentiles[endpoint]=percentiles;stability.extend(checks)
        extended,_,_=rank_consensus(scores,endpoint)
        extended=pd.concat([main_genes,extra.merge(gene_flags,on="canonical_id",validate="one_to_one")],ignore_index=True).merge(extended,on="canonical_id",validate="one_to_one")
        extended["diagnostic_only_extra"]=extended.canonical_id.isin(extra.canonical_id)
        extended.sort_values(["consensus_rank","canonical_id"]).to_csv(results/("Diagnostic_extended1957_"+endpoint+"_ranked.tsv"),sep="\t",index=False,na_rep="NA")
    combined.to_csv(results/"Main1881_combined_disease_model_evidence.tsv",sep="\t",index=False,na_rep="NA")
    pd.DataFrame(stability).to_csv(results/"Leave_one_encoding_out_stability.tsv",sep="\t",index=False)
    wide=main_scores.pivot(index="canonical_id",columns="model",values="predicted_log_affinity")
    correlation=wide.corr(method="spearman")
    correlation.to_csv(results/"Model_gene_rank_Spearman.tsv",sep="\t")
    plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False})
    names=list(correlation.columns)
    labels=[name.replace("_BindingDB_IC50"," [IC50]").replace("_BindingDB"," [Kd]") for name in names]
    fig,ax=plt.subplots(figsize=(10,9))
    image=ax.imshow(correlation,vmin=-1,vmax=1,cmap="RdBu_r")
    ax.set_xticks(range(len(names)),labels,rotation=60,ha="right");ax.set_yticks(range(len(names)),labels)
    for i in range(len(names)):
        for j in range(len(names)):ax.text(j,i,f"{correlation.iloc[i,j]:.2f}",ha="center",va="center",fontsize=7)
    fig.colorbar(image,ax=ax,label="Spearman correlation of gene-level predicted ranks")
    ax.set_title("Model agreement across 1,881 disease candidates\nMPNN shown as sensitivity only; shared training sources remain")
    fig.tight_layout()
    for ext in ["png","svg"]:fig.savefig(figures/("Model_rank_agreement."+ext),dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,6))
    ax.scatter(combined.Kd_consensus_percentile,combined.IC50_consensus_percentile,s=11,alpha=.45,c="#4E79A7",edgecolors="none")
    ax.set(xlabel="Kd panel: median model percentile",ylabel="IC50 panel: median model percentile",title="Endpoint-specific consensus agreement (n=1,881)",xlim=(0,1.02),ylim=(0,1.02))
    fig.tight_layout()
    for ext in ["png","svg"]:fig.savefig(figures/("Kd_vs_IC50_consensus."+ext),dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,9))
    for ax,endpoint in zip(axes,["Kd","IC50"]):
        top=all_tables[endpoint].head(20)
        matrix=all_percentiles[endpoint].reindex(top.canonical_id)
        heat=ax.imshow(matrix,cmap="viridis",vmin=0,vmax=1,aspect="auto")
        ax.set_yticks(range(20),top.gene_symbol)
        ax.set_xticks(range(4),[name.split("_BindingDB")[0] for name in matrix.columns],rotation=45,ha="right")
        ax.set_title(endpoint+" panel: first20 by consensus\nDisplay only; no target-selection cutoff")
    fig.subplots_adjust(left=.1,right=.86,bottom=.18,top=.92,wspace=.5)
    color_axis=fig.add_axes([.90,.22,.018,.56])
    fig.colorbar(heat,cax=color_axis,label="Per-model gene percentile")
    for ext in ["png","svg"]:fig.savefig(figures/("Top20_encoding_agreement."+ext),dpi=180)
    plt.close(fig)
    counts={"main_genes":len(main_genes),"sensitivity_extra_genes":len(extra),"protein_pairs":len(catalog),
            "models":len(panel),"protein_model_predictions":len(raw),"gene_model_predictions":len(scores),
            "main_Kd_vs_IC50_consensus_Spearman":combined.Kd_consensus_percentile.corr(combined.IC50_consensus_percentile,method="spearman"),
            "main_CNN_truncated_genes":int(main_genes.any_CNN_truncation.sum()),
            "main_GENCODE_translation_only_genes":int(main_genes.any_GENCODE_translation_only.sum()),
            "main_nonstandard_aminoacid_genes":int(main_genes.any_nonstandard_amino_acid.sum()),
            "main_multiple_protein_genes":int(main_genes.n_input_protein_sequences.gt(1).sum()),
            "Kd_first10_display":all_tables["Kd"][["gene_symbol","canonical_id","consensus_rank","model_percentile_IQR"]].head(10).to_dict("records"),
            "IC50_first10_display":all_tables["IC50"][["gene_symbol","canonical_id","consensus_rank","model_percentile_IQR"]].head(10).to_dict("records"),
            "no_final_target_selected":True,"clinical_or_experimental_validation":False}
    (results/"Prediction_summary_metrics.json").write_text(json.dumps(counts,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(counts,ensure_ascii=True),flush=True)

if __name__=="__main__":main()
