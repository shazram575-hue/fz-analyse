"""Retrieve public chemical identity evidence only; no study gene lists transmitted."""
import concurrent.futures
import hashlib
import json
import pathlib
import re
import urllib.request

ROOT=pathlib.Path(__file__).resolve().parents[1]
SOURCES=[
    ("ChEMBL456512_molecule.json","https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL456512.json"),
    ("ChEMBL_Fang2022_document.json","https://www.ebi.ac.uk/chembl/api/data/document.json?doi=10.1016%2Fj.phymed.2022.154151"),
    ("Shang2025_crossref.json","https://api.crossref.org/works/10.1515/ncrs-2025-0148"),
    ("Shang2025_RG.html","https://www.researchgate.net/publication/392401026_Crystal_structure_of_spiropachysine_C31H46N2O"),
    ("Shang2025_COD.json","https://www.crystallography.net/cod/result?doi=10.1515%2Fncrs-2025-0148&format=json"),
    ("CCDC2433196.html","https://www.ccdc.cam.ac.uk/structures/StructureSummary/form?pid=2433196"),
    ("PubChem209254_full_record.json","https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/209254/JSON"),
    ("NP_MRD_NP0066471.html","https://np-mrd.org/natural_products/NP0066471"),
]

def fetch(item):
    name,url=item
    record={"name":name,"url":url}
    try:
        request=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (Workflow3 chemical identity research)","Accept":"application/json,text/html,application/pdf,*/*"})
        with urllib.request.urlopen(request,timeout=30) as response:
            body=response.read(15_000_001)
            record.update(status=response.status,content_type=response.headers.get("Content-Type"))
        if len(body)>15_000_000:raise ValueError("Source exceeds 15 MB limit")
        if not body:raise ValueError("Empty response")
        (ROOT/"sources"/name).write_bytes(body)
        record.update(bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
        if name.endswith(".json"):
            try:
                parsed=json.loads(body)
                if name=="ChEMBL456512_molecule.json":record["key_content"]={k:parsed.get(k) for k in ["molecule_chembl_id","pref_name","molecule_structures","molecule_synonyms"]}
                elif name=="ChEMBL_Fang2022_document.json":record["key_content"]=parsed
                elif name=="Shang2025_crossref.json":record["key_content"]={k:parsed.get("message",{}).get(k) for k in ["link","resource","license","URL"]}
                elif name=="Shang2025_COD.json":record["key_content"]=parsed
            except ValueError:record["parse_error"]="Response is not JSON"
        else:
            text=body.decode("utf-8",errors="replace")
            links=re.findall(r'''(?:href|src)\s*=\s*["']([^"']+)["']''',text,re.I)
            record["relevant_links"]=[link for link in links if any(word in link.lower() for word in ["download","fulltext",".pdf",".cif","2433196"] )][:25]
    except Exception as exc:record["error"]=str(exc)
    print(json.dumps(record,ensure_ascii=True),flush=True)
    return record

def main():
    (ROOT/"sources").mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,SOURCES))
    (ROOT/"Source_retrieval_audit.json").write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding="utf-8")

if __name__=="__main__":main()
