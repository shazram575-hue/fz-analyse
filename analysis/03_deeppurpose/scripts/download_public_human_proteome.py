"""Download the entire public reviewed human dataset. No candidate IDs are read or transmitted."""
import gzip
import hashlib
import json
import pathlib
import time
import urllib.parse
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "public_human_proteome"

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    query = urllib.parse.urlencode({"query": "(organism_id:9606) AND (reviewed:true)", "format": "json", "compressed": "true",
                                   "fields": "accession,id,gene_primary,organism_id,length,sequence,xref_ensembl"})
    url = "https://rest.uniprot.org/uniprotkb/stream?" + query
    destination = OUT / "Reviewed_human_UniProt.json.gz"
    if not destination.exists():
        partial = destination.with_suffix(".gz.part")
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "Workflow3-public-human-proteome-download/1.0"})
                with urllib.request.urlopen(request, timeout=60) as response, partial.open("wb") as out:
                    headers = {k: v for k, v in response.headers.items() if k.lower() in ["x-uniprot-release", "x-uniprot-release-date", "x-total-results", "content-type"]}
                    total = 0
                    while chunk := response.read(1024 * 1024):
                        total += len(chunk)
                        if total > 100_000_000:
                            raise ValueError("Download exceeds 100 MB limit")
                        out.write(chunk)
                with gzip.open(partial, "rt", encoding="utf-8") as handle:
                    data = json.load(handle)
                if not data.get("results"):
                    raise ValueError("No public human entries in downloaded response")
                partial.replace(destination)
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(3)
        audit = {"source_url": url, "query_scope": "all_reviewed_human_UniProt_records", "candidate_ids_transmitted": False,
                 "headers": headers, "n_entries": len(data["results"]), "bytes": destination.stat().st_size,
                 "sha256": hashlib.sha256(destination.read_bytes()).hexdigest()}
        (OUT / "Download_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    else:
        audit = json.loads((OUT / "Download_audit.json").read_text(encoding="utf-8"))
    print(json.dumps(audit), flush=True)

if __name__ == "__main__":
    main()
