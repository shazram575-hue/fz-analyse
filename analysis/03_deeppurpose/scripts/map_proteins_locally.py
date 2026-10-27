"""Strictly local candidate mapping; this module performs no network operations."""
import csv
import gzip
import hashlib
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "protein_mapping_preparation"

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with (ROOT / "candidate_scope_review" / "Union_for_preparation_only.tsv").open(encoding="utf-8") as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))
    ids = {row["canonical_id"] for row in genes}
    with gzip.open(ROOT / "public_human_proteome" / "Reviewed_human_UniProt.json.gz", "rt", encoding="utf-8") as handle:
        entries = json.load(handle)["results"]
    rows = []
    for protein in entries:
        if protein.get("organism", {}).get("taxonId") != 9606:
            continue
        matched = set()
        for crossref in protein.get("uniProtKBCrossReferences", []):
            if crossref.get("database") == "Ensembl":
                for prop in crossref.get("properties", []):
                    if prop.get("key") == "GeneId":
                        gene = prop["value"].split(".")[0]
                        if gene in ids:
                            matched.add(gene)
        sequence = protein.get("sequence", {}).get("value", "")
        if matched and not sequence:
            raise ValueError("Matched UniProt entry lacks sequence")
        accession = protein["primaryAccession"]
        names = [g.get("geneName", {}).get("value", "") for g in protein.get("genes", [])]
        for gene in sorted(matched):
            rows.append({"canonical_id": gene, "uniprot_accession": accession, "uniprot_gene_names": ";".join(names),
                         "organism_id": 9606, "entry_type": protein.get("entryType", ""), "sequence_length": len(sequence),
                         "protein_sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                         "CNN_truncates_after1000": len(sequence) > 1000,
                         "nonstandard_amino_acids": ";".join(sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"))),
                         "mapping_method": "local_exact_Ensembl_GeneId_from_public_reviewed_human_dataset",
                         "scope_status": "preparation_union_not_main_pool_selection"})
    rows.sort(key=lambda row: (row["canonical_id"], row["uniprot_accession"]))
    if rows:
        with (OUT / "Reviewed_human_protein_mappings.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    mapped = {row["canonical_id"] for row in rows}
    with (OUT / "Unmapped_in_reviewed_human.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["canonical_id", "gene_symbol"], delimiter="\t")
        writer.writeheader()
        writer.writerows(row for row in genes if row["canonical_id"] not in mapped)
    counts = {gene: sum(row["canonical_id"] == gene for row in rows) for gene in mapped}
    summary = {"input_genes_preparation_union": len(ids), "mapped_genes": len(mapped), "gene_protein_rows": len(rows),
               "unmapped_reviewed_human": len(ids - mapped), "genes_with_multiple_reviewed_entries": sum(n > 1 for n in counts.values()),
               "CNN_truncated_rows": sum(row["CNN_truncates_after1000"] for row in rows),
               "candidate_ids_transmitted": False, "ambiguous_entries": "all_retained_no_forced_choice"}
    (OUT / "Mapping_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary), flush=True)

if __name__ == "__main__":
    main()
