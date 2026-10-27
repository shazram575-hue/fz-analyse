"""Resumable exact-Ensembl mapping to reviewed UniProt records for scope alternatives."""
import csv
import hashlib
import json
import pathlib
import time
import urllib.parse
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "protein_mapping_preparation"
BASE = "https://rest.uniprot.org"

def request_json(url, data=None):
    request = urllib.request.Request(url, data=data, headers={"User-Agent": "Workflow3-DeepPurpose/1.0"})
    with urllib.request.urlopen(request, timeout=40) as response:
        body = response.read(150_000_001)
        if len(body) > 150_000_000:
            raise ValueError("UniProt response exceeds 150 MB limit")
        headers = {k: v for k, v in response.headers.items() if k.lower() in ["x-uniprot-release", "x-uniprot-release-date", "x-total-results"]}
    return json.loads(body), headers

def main():
    raise RuntimeError("DISABLED: external candidate-ID submission was rejected. Use download_public_human_proteome.py followed by map_proteins_locally.py; candidate lists must remain local.")
    OUT.mkdir(parents=True, exist_ok=True)
    input_path = ROOT / "candidate_scope_review" / "Union_for_preparation_only.tsv"
    with input_path.open(encoding="utf-8") as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))
    ids = sorted({row["canonical_id"] for row in genes})
    manifest_path = OUT / "Mapping_job.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert manifest["ids_sha256"] == hashlib.sha256("\n".join(ids).encode()).hexdigest()
    else:
        payload = urllib.parse.urlencode({"from": "Ensembl", "to": "UniProtKB-Swiss-Prot", "ids": ",".join(ids)}).encode()
        job, headers = request_json(BASE + "/idmapping/run", payload)
        manifest = {"jobId": job["jobId"], "n_input_ids": len(ids), "source": "Ensembl", "target": "UniProtKB-Swiss-Prot",
                    "scope": "union_of_alternative_disease_inputs_for_preparation_only",
                    "ids_sha256": hashlib.sha256("\n".join(ids).encode()).hexdigest()}
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    deadline = time.monotonic() + 45
    while time.monotonic() < deadline:
        status, _ = request_json(BASE + "/idmapping/status/" + manifest["jobId"])
        (OUT / "Mapping_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
        if status.get("jobStatus") in ["NEW", "RUNNING"]:
            time.sleep(3)
            continue
        if status.get("jobStatus"):
            raise ValueError(status)
        break
    else:
        print("MAPPING_PENDING", manifest["jobId"], flush=True)
        return
    fields = "accession,id,gene_primary,organism_id,length,sequence,xref_ensembl"
    query = urllib.parse.urlencode({"format": "json", "fields": fields})
    results_path = OUT / "UniProt_mapping_results.json"
    if results_path.exists():
        result = json.loads(results_path.read_text(encoding="utf-8"))
    else:
        result, headers = request_json(BASE + "/idmapping/uniprotkb/results/stream/" + manifest["jobId"] + "?" + query)
        results_path.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
        (OUT / "UniProt_response_headers.json").write_text(json.dumps(headers, indent=2), encoding="utf-8")
    rows = []
    for pair in result.get("results", []):
        protein = pair["to"]
        if not isinstance(protein, dict):
            raise ValueError("UniProt result did not include protein records")
        organism = protein.get("organism", {}).get("taxonId")
        sequence = protein.get("sequence", {}).get("value", "")
        if organism != 9606:
            continue
        accession = protein["primaryAccession"]
        names = [g.get("geneName", {}).get("value", "") for g in protein.get("genes", [])]
        rows.append({"canonical_id": pair["from"], "uniprot_accession": accession, "uniprot_gene_names": ";".join(names),
                     "organism_id": organism, "entry_type": protein.get("entryType", ""), "sequence_length": len(sequence),
                     "protein_sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                     "CNN_truncates_after1000": len(sequence) > 1000,
                     "nonstandard_amino_acids": ";".join(sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"))),
                     "mapping_method": "exact_Ensembl_UniProt_ID_mapping",
                     "scope_status": "preparation_union_not_main_pool_selection"})
    if rows:
        keys = list(rows[0])
        with (OUT / "Reviewed_human_protein_mappings.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    mapped = {row["canonical_id"] for row in rows}
    missing = [row for row in genes if row["canonical_id"] not in mapped]
    with (OUT / "Unmapped_in_reviewed_human.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["canonical_id", "gene_symbol"], delimiter="\t")
        writer.writeheader()
        writer.writerows(missing)
    summary = {"n_input_genes": len(ids), "n_mapped_genes": len(mapped), "n_gene_protein_rows": len(rows),
               "n_unmapped_reviewed_human": len(missing), "n_truncated_gene_protein_rows": sum(row["CNN_truncates_after1000"] for row in rows),
               "n_mapping_multiple_genes": sum(sum(r["canonical_id"] == gene for r in rows) > 1 for gene in mapped),
               "sequence_selection": "all_exact_reviewed_human_matches_preserved_no_ambiguous_choice"}
    (OUT / "Mapping_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary), flush=True)

if __name__ == "__main__":
    main()
