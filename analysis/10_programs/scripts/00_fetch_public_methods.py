from pathlib import Path
import sys,json,urllib.request,xml.etree.ElementTree as ET
out=Path(sys.argv[1])/'outputs/MFSD10_程序发现与独立验证_2026-10-26'
urls={
 'Wu2021_fulltext.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9044823/fullTextXML',
 'GSE176078_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/soft/GSE176078_family.soft.gz',
 'BrCa_cell_atlas_README.md':'https://raw.githubusercontent.com/Swarbricklab-code/BrCa_cell_atlas/main/README.md'}
records=[]
for name,url in urls.items():
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (public research methods)'})
        with urllib.request.urlopen(req,timeout=60) as r:b=r.read()
        (out/'sources'/name).write_bytes(b);records.append(dict(name=name,url=url,bytes=len(b),success=True));print(name,len(b),flush=True)
        if name.endswith('.xml'):
            t=ET.fromstring(b);(out/'sources/Wu2021_fulltext.txt').write_text('\n'.join(' '.join(n.itertext()) for n in t.findall('.//p')),encoding='utf-8')
    except Exception as e:records.append(dict(name=name,url=url,success=False,error=str(e)));print(name,str(e),flush=True)
(out/'sources/Download_records.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
