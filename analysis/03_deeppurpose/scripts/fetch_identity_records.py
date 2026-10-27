"""Public chemical and primary-reference retrieval; no research gene list is read."""
import concurrent.futures
import hashlib
import json
import pathlib
import urllib.request

ROOT=pathlib.Path(__file__).resolve().parents[1]
SOURCES=[
    ("PubChem_209254_properties.json","https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/209254/property/CanonicalSMILES,IsomericSMILES,InChI,InChIKey,MolecularFormula,IUPACName/JSON"),
    ("PubChem_500127_properties.json","https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/500127/property/CanonicalSMILES,IsomericSMILES,InChI,InChIKey,MolecularFormula,IUPACName/JSON"),
    ("Chang1998_primary.pdf","https://pubs.acs.org/doi/pdf/10.1021/np980162x"),
    ("Chang1998_primary.html","https://pubs.acs.org/doi/10.1021/np980162x"),
    ("Chiu1990_primary.html","https://www.sciencedirect.com/science/article/pii/003194229085362J"),
]

def fetch(item):
    name,url=item
    record={"name":name,"url":url}
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (Workflow3 structure verification)"})
        with urllib.request.urlopen(req,timeout=25) as response:
            content=response.read(20_000_001)
            record.update(status=response.status,content_type=response.headers.get("Content-Type"))
        if len(content)>20_000_000:raise ValueError("Response exceeds 20 MB")
        if name.endswith(".pdf") and not content.startswith(b"%PDF"):raise ValueError("Response is not a PDF")
        if name.endswith(".json"):json.loads(content)
        if not content:raise ValueError("Empty response")
        path=ROOT/"structure_sources"/name
        path.write_bytes(content)
        record.update(saved=str(path),bytes=len(content),sha256=hashlib.sha256(content).hexdigest())
    except Exception as exc:record["error"]=str(exc)
    return record

def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,SOURCES))
    (ROOT/"structure_sources"/"Identity_record_retrieval.json").write_text(json.dumps(records,indent=2),encoding="utf-8")
    print(json.dumps(records),flush=True)

if __name__=="__main__":main()
