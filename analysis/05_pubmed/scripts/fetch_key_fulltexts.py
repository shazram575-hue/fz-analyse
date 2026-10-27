import requests,json,time,xml.etree.ElementTree as ET
from pathlib import Path
B=Path(__file__).resolve().parents[1]
s=requests.Session()
for pmc in ['PMC12049414','PMC4747093','PMC5929433','PMC6481349','PMC9385710']:
    p=B/'sources'/f'{pmc}.xml'
    if p.exists():continue
    for url in [f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML',f'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmc}/unicode']:
        try:
            r=s.get(url,timeout=60);r.raise_for_status();ET.fromstring(r.content)
            p.write_bytes(r.content);print(pmc,len(r.content),flush=True);break
        except Exception as e:print(type(e).__name__,pmc,flush=True)
for url,name in [('https://www.nature.com/articles/s41594-024-01237-2','SAMD8_structure_paper.html'),('https://cgp.iiarjournals.org/content/20/1/88.full','SAMD8_SH003_paper.html'),('https://www.nature.com/articles/s41598-018-26429-x','Anacardic_breast_transcriptome.html'),('https://www.spandidos-publications.com/2187/downloadSupplementary','FEN1_supplement.pdf')]:
    p=B/'sources'/name
    if not p.exists():
        try:
            r=s.get(url,timeout=60);r.raise_for_status();p.write_bytes(r.content);print(name,len(r.content),flush=True)
        except Exception as e:print(name,type(e).__name__,flush=True)
