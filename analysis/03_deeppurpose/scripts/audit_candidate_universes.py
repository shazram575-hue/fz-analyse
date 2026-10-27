"""Construct reviewable alternative disease intersections without selecting a main pool."""
import hashlib
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
FOUNDATION = WORKSPACE / "outputs" / "单细胞基础补全_2026-10-23"
STAGE1 = WORKSPACE / "outputs" / "第一阶段_2026-10-23"

def flag(series):
    return series.astype(str).str.upper().eq("TRUE")

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    result = ROOT / "candidate_scope_review"
    result.mkdir(parents=True, exist_ok=True)
    evidence_path = FOUNDATION / "results" / "Within_type_combined_gene_evidence.tsv"
    bulk_path = STAGE1 / "results" / "Bulk_exploratory_DE_protein_coding.tsv"
    evidence = pd.read_csv(evidence_path, sep="\t", low_memory=False)
    bulk = pd.read_csv(bulk_path, sep="\t")
    assert bulk.ensembl_stable.is_unique
    assert len(bulk) == 7494
    assert (bulk.gene_type == "protein_coding").all()
    assert ((bulk.logFC.abs() > 0.585) & (bulk["P.Value"] < 0.05)).all()
    assert not evidence.duplicated(["cell_type", "canonical_id"]).any()
    biological = ~flag(evidence.diagnostic_annotation)
    protein = evidence.gene_type.eq("protein_coding")
    definitions = {
        "celllevel_6biological": flag(evidence.SC_DE_pass) & biological,
        "celllevel_all8_including_diagnostic": flag(evidence.SC_DE_pass),
        "pseudobulk_4biological_estimable": flag(evidence.PB_DE_pass) & biological,
        "pseudobulk_all6_estimable_including_diagnostic": flag(evidence.PB_DE_pass),
    }
    bulk_names = {"ensembl_stable": "canonical_id", "gene_name": "bulk_gene_symbol", "gene_type": "bulk_gene_type"}
    bulk_join = bulk.rename(columns={**bulk_names, **{c: "bulk_" + c for c in bulk if c not in bulk_names}})
    summaries = []
    unions = set()
    for name, selection in definitions.items():
        selected = evidence.loc[selection & protein].copy()
        matched = selected.merge(bulk_join, on="canonical_id", how="inner", validate="many_to_one")
        matched.to_csv(result / (name + "_gene_celltype_evidence.tsv"), sep="\t", index=False, na_rep="NA")
        pool = matched.groupby("canonical_id", sort=True).agg(
            gene_symbol=("gene_symbol", "first"),
            cell_types=("cell_type", lambda x: ";".join(sorted(set(x)))),
            n_cell_types=("cell_type", "nunique"),
            diagnostic_evidence_present=("diagnostic_annotation", lambda x: bool(flag(x).any())),
            bulk_log2FC=("bulk_logFC", "first"),
            bulk_P=("bulk_P.Value", "first"),
            bulk_BH=("bulk_adj.P.Val", "first"),
        ).reset_index()
        pool["scope_status"] = "alternative_for_user_review_not_selected_main_pool"
        pool.to_csv(result / (name + "_genes.tsv"), sep="\t", index=False, na_rep="NA")
        unions.update(pool.canonical_id)
        summaries.append({"definition": name, "sc_de_unique_all_gene_types": evidence.loc[selection, "canonical_id"].nunique(),
                          "sc_de_unique_protein_coding": selected.canonical_id.nunique(),
                          "intersection_unique_genes": len(pool), "gene_celltype_evidence_rows": len(matched),
                          "cell_types_in_intersection": ";".join(sorted(matched.cell_type.unique()))})
    summary = pd.DataFrame(summaries)
    summary.to_csv(result / "Alternative_scope_counts.tsv", sep="\t", index=False)
    universe = evidence.loc[evidence.canonical_id.isin(unions) & protein, ["canonical_id", "gene_symbol"]].drop_duplicates()
    assert universe.canonical_id.is_unique
    universe.sort_values("canonical_id").to_csv(result / "Union_for_preparation_only.tsv", sep="\t", index=False)
    paths = [evidence_path, bulk_path, STAGE1 / "results" / "Pool_A_candidate_genes.tsv"]
    (result / "Source_SHA256.json").write_text(json.dumps({str(p): sha256(p) for p in paths}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(summary.to_json(orient="records", force_ascii=True))

if __name__ == "__main__":
    main()
