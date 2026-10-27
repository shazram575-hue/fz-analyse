import json,hashlib,tarfile,zipfile
from pathlib import Path
import requests
BASE=Path(__file__).resolve().parents[1]
software=BASE/'software';software.mkdir(exist_ok=True)
jobs=[]
for repo,asset in [('AutoDock-Vina','vina_1.2.7_win.exe'),('p2rank','p2rank_2.5.1.tar.gz')]:
    d=json.loads((BASE/'sources'/f'{repo}_release.json').read_text())
    a=next(x for x in d['assets'] if x['name']==asset)
    jobs.append((a['browser_download_url'],asset,a['size'],a.get('digest')))
j=json.loads((BASE/'sources/java_release.json').read_text())[0]['binary']['package']
jobs.append((j['link'],j['name'],j['size'],'sha256:'+j['checksum']))
manifest=[]
for url,name,size,digest in jobs:
    path=software/name
    if not path.exists() or path.stat().st_size!=size:
        with requests.get(url,stream=True,timeout=120) as r:
            r.raise_for_status()
            with path.open('wb') as f:
                for b in r.iter_content(1024*1024):f.write(b)
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    assert path.stat().st_size==size
    if digest:assert digest=='sha256:'+sha
    manifest.append({'url':url,'file':name,'bytes':size,'sha256':sha,'official_digest_checked':bool(digest)})
    print('DOWNLOADED',name,size,flush=True)
    if name.endswith('.zip'):
        with zipfile.ZipFile(path) as z:z.extractall(software)
    elif name.endswith('.tar.gz'):
        with tarfile.open(path) as t:t.extractall(software,filter='data')
(BASE/'audit/software_download_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
print('TOOLS_DONE',flush=True)
