"""Retrieve public primary-source structure evidence, preserving failures."""
import concurrent.futures
import hashlib
import json
import pathlib
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCES = [
    ("Fang2022_Fig1.jpg", "https://ars.els-cdn.com/content/image/1-s2.0-S094471132200229X-gr1.jpg"),
    ("Fang2022_SI_mmc1.pdf", "https://ars.els-cdn.com/content/image/1-s2.0-S094471132200229X-mmc1.pdf"),
    ("Kikuchi1975_primary.pdf", "https://www.jstage.jst.go.jp/article/cpb1958/23/2/23_2_416/_pdf"),
    ("Shang2025_crystal.html", "https://www.degruyterbrill.com/document/doi/10.1515/ncrs-2025-0148/html"),
    ("Shang2025_crystal.pdf", "https://www.degruyterbrill.com/document/doi/10.1515/ncrs-2025-0148/pdf"),
    ("Fang2022_publisher.html", "https://www.sciencedirect.com/science/article/pii/S094471132200229X"),
]

def fetch(source):
    name, url = source
    record = {"name": name, "url": url}
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Workflow3 scientific source verification)"})
        with urllib.request.urlopen(request, timeout=25) as response:
            content = response.read(20_000_001)
            record.update(status=response.status, content_type=response.headers.get("Content-Type"))
        if len(content) > 20_000_000:
            raise ValueError("Source exceeds 20 MB limit")
        if name.endswith(".pdf") and not content.startswith(b"%PDF"):
            raise ValueError("Response is not a PDF")
        if name.endswith(".jpg") and not content.startswith(b"\xff\xd8"):
            raise ValueError("Response is not a JPEG")
        path = ROOT / "structure_sources" / name
        path.write_bytes(content)
        record.update(saved=str(path), bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
    except Exception as exc:
        record["error"] = str(exc)
    return record

def main():
    (ROOT / "structure_sources").mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(fetch, SOURCES))
    (ROOT / "structure_sources" / "Retrieval_audit.json").write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(records, ensure_ascii=True))

if __name__ == "__main__":
    main()
