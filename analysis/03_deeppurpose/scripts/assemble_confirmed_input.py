"""Freeze the human-confirmed 1,881-gene scope; retain all mapping limitations."""
import hashlib
import json
import pathlib
import pandas as pd

ROOT=pathlib.Path(__file__).resolve().parents[1]

def main():
    out=ROOT/"main_input"
    out.mkdir(parents=True,exist_ok=True)
    review=ROOT/"candidate_scope_review"
    genes=pd.read_csv(review/"celllevel_6biological_genes.tsv",sep="\t")
    evidence=pd.read_csv(review/"celllevel_6biological_gene_celltype_evidence.tsv",sep="\t",low_memory=False)
    extended=pd.read_csv(review/"celllevel_all8_including_diagnostic_genes.tsv",sep="\t")
    assert len(genes)==1881 and genes.canonical_id.is_unique
    assert not evidence.diagnostic_annotation.astype(str).str.lower().eq("true").any()
    annotations=[]
    for gene,group in evidence.groupby("canonical_id"):
        up=group.loc[group.SC_log2FC>0,"cell_type"]
        down=group.loc[group.SC_log2FC<0,"cell_type"]
        pb=group.loc[group.PB_DE_pass.astype(str).str.lower().eq("true"),"cell_type"]
        rare=group.loc[group.PB_status.eq("not_run_insufficient_donors"),"cell_type"]
        annotations.append({"canonical_id":gene,"SC_up_cell_types":";".join(sorted(up)),
                            "SC_down_cell_types":";".join(sorted(down)),"PB_DE_support_cell_types":";".join(sorted(pb)),
                            "donor_limited_SC_cell_types":";".join(sorted(rare)),
                            "has_only_donor_limited_SC_support":len(rare)==len(group)})
    genes=genes.merge(pd.DataFrame(annotations),on="canonical_id",validate="one_to_one")
    genes["scope_status"]="user_confirmed_main_input_2026-10-24_not_validated_drug_targets"
    genes.to_csv(out/"Disease_matching_candidates_1881.tsv",sep="\t",index=False,na_rep="NA")
    evidence.to_csv(out/"Disease_candidates_gene_celltype_evidence.tsv",sep="\t",index=False,na_rep="NA")
    additional=extended.loc[~extended.canonical_id.isin(genes.canonical_id)].copy()
    assert len(additional)==76
    additional["scope_status"]="diagnostic_annotation_sensitivity_only_not_main_input"
    additional.to_csv(out/"Additional_diagnostic_only_76genes.tsv",sep="\t",index=False,na_rep="NA")
    mapping=pd.read_csv(ROOT/"protein_mapping_preparation"/"Direct_mapping_with_symbol_checks.tsv",sep="\t",keep_default_na=False)
    merged=genes.merge(mapping.drop(columns="gene_symbol"),on="canonical_id",how="left",validate="one_to_many",suffixes=("","_mapping"))
    merged["mapping_available"]=merged.uniprot_accession.notna()
    counts=merged.loc[merged.mapping_available].groupby("canonical_id").size()
    merged["n_reviewed_entries_for_gene"]=merged.canonical_id.map(counts).fillna(0).astype(int)
    merged["inference_status"]="not_run_drug_identity_under_review"
    merged.to_csv(out/"Main_gene_protein_mapping_review.tsv",sep="\t",index=False,na_rep="NA")
    decisions={"date":"2026-10-24","human_selection":"细胞级六个生物学大类：1,881 基因；混合/未定另作敏感性（建议）",
               "main_genes":1881,"additional_diagnostic_only_genes":76,"original_Pool_A_813":"preserved_not_replaced",
               "single_cell_direction":"both_up_and_down_existing_within_type_thresholds",
               "bulk_rule":"existing_abs_log2FC_gt0.585_raw_P_lt0.05_both_directions_protein_coding",
               "PB_role":"support_not_main_input_hard_gate","novel_drug_predictions":False,
               "n_exact_ID_mapped_main_genes":int(merged.loc[merged.mapping_available,"canonical_id"].nunique()),
               "n_main_genes_without_exact_reviewed_mapping":int(genes.canonical_id.nunique()-merged.loc[merged.mapping_available,"canonical_id"].nunique()),
               "files_sha256":{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in out.glob("*.tsv")}}
    (out/"Confirmed_scope_and_provenance.json").write_text(json.dumps(decisions,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in decisions.items() if k!="files_sha256"},ensure_ascii=True),flush=True)

if __name__=="__main__":main()
