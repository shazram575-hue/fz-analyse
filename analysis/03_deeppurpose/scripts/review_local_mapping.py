"""Annotate exact-ID mappings and optional symbol-only rescue evidence locally."""
import collections
import gzip
import json
import pathlib
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]

def main():
    out = ROOT / "protein_mapping_preparation"
    genes = pd.read_csv(ROOT / "candidate_scope_review" / "Union_for_preparation_only.tsv", sep="\t")
    direct = pd.read_csv(out / "Reviewed_human_protein_mappings.tsv", sep="\t", keep_default_na=False)
    direct = direct.merge(genes, on="canonical_id", validate="many_to_one")
    direct["primary_symbol_exact_match"] = [symbol in names.split(";") for symbol,names in zip(direct.gene_symbol,direct.uniprot_gene_names)]
    direct["mapping_tier"] = "exact_Ensembl_GeneId"
    direct.to_csv(out / "Direct_mapping_with_symbol_checks.tsv", sep="\t", index=False)
    entries = json.load(gzip.open(ROOT / "public_human_proteome" / "Reviewed_human_UniProt.json.gz", "rt", encoding="utf-8"))["results"]
    symbols = collections.defaultdict(dict)
    for protein in entries:
        for gene in protein.get("genes", []):
            symbol = gene.get("geneName", {}).get("value", "")
            if symbol:
                symbols[symbol][protein["primaryAccession"]] = protein
    missing = genes.loc[~genes.canonical_id.isin(direct.canonical_id)]
    rescue = []
    for row in missing.itertuples(index=False):
        matches = symbols[row.gene_symbol]
        rescue.append({"canonical_id":row.canonical_id,"gene_symbol":row.gene_symbol,
                       "primary_symbol_candidate_accessions":";".join(sorted(matches)),
                       "n_primary_symbol_matches":len(matches),
                       "mapping_tier":"unique_primary_symbol_only_requires_ID_reconciliation" if len(matches)==1 else "unresolved",
                       "automatic_inclusion":False})
    pd.DataFrame(rescue).to_csv(out / "Symbol_only_rescue_review_not_automatically_included.tsv",sep="\t",index=False)
    summaries=[]
    for path in sorted((ROOT/"candidate_scope_review").glob("*_genes.tsv")):
        pool=pd.read_csv(path,sep="\t")
        matched=direct[direct.canonical_id.isin(pool.canonical_id)]
        n_per_gene=matched.groupby("canonical_id").size()
        summaries.append({"definition":path.stem.removesuffix("_genes"),"input_genes":len(pool),
                          "exact_id_mapped_genes":matched.canonical_id.nunique(),"gene_protein_rows":len(matched),
                          "unmapped_genes":len(pool)-matched.canonical_id.nunique(),
                          "genes_with_multiple_reviewed_entries":int((n_per_gene>1).sum()),
                          "CNN_truncated_gene_protein_rows":int(matched.CNN_truncates_after1000.astype(str).str.lower().eq("true").sum()),
                          "symbol_conflict_gene_protein_rows":int((~matched.primary_symbol_exact_match).sum())})
    pd.DataFrame(summaries).to_csv(out/"Coverage_by_alternative_scope.tsv",sep="\t",index=False)
    print(json.dumps(summaries),flush=True)

if __name__=="__main__":
    main()
