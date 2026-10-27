"""Retrieve public primary nomenclature sources and expose source-linked PDF URLs."""
import concurrent.futures
import hashlib
import json
import pathlib
import re
import urllib.parse
import urllib.request

ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/"nomenclature_review"
SOURCES=[
    ("Qiu1990_primary.html","https://sioc-journal.cn/Jwk_yjhx/EN/abstract/abstract330661.shtml"),
    ("Qiu1992_primary.html","https://www.integrativebiology.ac.cn/EN/abstract/abstract31331.shtml"),
    ("Qiu1992_primary_cn.html","https://www.integrativebiology.ac.cn/CN/abstract/abstract31331.shtml"),
    ("Qiu1992_plantdiversity.html","https://journal.kib.ac.cn/EN/abstract/abstract31331.shtml"),
    ("Kikuchi1968_pubmed.html","https://pubmed.ncbi.nlm.nih.gov/5651213/"),
    ("Chang1998_primary.html","https://pubs.acs.org/doi/10.1021/np980162x"),
]

def fetch(item):
    name,url=item
    record={"name":name,"url":url}
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (Workflow3 scientific nomenclature check)"})
        with urllib.request.urlopen(req,timeout=25) as response:
            raw=response.read(10_000_000)
            record.update(status=response.status,content_type=response.headers.get("Content-Type"))
        (OUT/name).write_bytes(raw)
        text=raw.decode("utf-8",errors="replace")
        links=re.findall(r'''(?:href|src)\s*=\s*["']([^"']+)["']''',text,re.I)
        record.update(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),
                      relevant_links=[urllib.parse.urljoin(url,link) for link in links if any(x in link.lower() for x in ["pdf","download","filedown","31331","330661"])])
        scripts=[line[:1600] for line in text.splitlines() if any(x in line.lower() for x in ["pdfdown","downloadpdf","downloadfile","pdffile","article_pdf"])][:15]
        record["download_script_lines"]=scripts
    except Exception as exc:record["error"]=str(exc)
    return record

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,SOURCES))
    (OUT/"Retrieval_audit.json").write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(records,ensure_ascii=True),flush=True)

if __name__=="__main__":main()
