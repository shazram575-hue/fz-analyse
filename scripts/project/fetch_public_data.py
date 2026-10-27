"""Plan-only by default. Explicitly stream verified public bytes; never overwrite originals."""
from pathlib import Path
import argparse,csv,json,hashlib,shutil,urllib.request,urllib.parse,os
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--dataset',required=True);p.add_argument('--asset');p.add_argument('--file-id');p.add_argument('--all',action='store_true');p.add_argument('--execute',action='store_true');p.add_argument('--destination',type=Path);a=p.parse_args();root=a.root.resolve()
config=json.loads((root/'config/public_datasets.json').read_text(encoding='utf-8'))
if a.dataset=='TCGA_BRCA_STAR':
 with (root/config['TCGA_BRCA_STAR']['manifest']).open(encoding='utf-8') as f:items=list(csv.DictReader(f,delimiter='\t'))
 if a.file_id:items=[r for r in items if r['file_id']==a.file_id]
elif a.dataset in config:
 items=config[a.dataset].get('assets',[])
 if a.asset:items=[r for r in items if r['id']==a.asset]
else:raise SystemExit('Unknown dataset; inspect config/public_datasets.json')
if not items:raise SystemExit('No matching assets')
dest=(a.destination or Path(os.environ.get('BIOINFO_DATA_ROOT',str(root/'data/external'))))/a.dataset
total=sum(int(r['bytes']) for r in items);print(json.dumps(dict(dataset=a.dataset,files=len(items),bytes=total,destination=str(dest),mode='execute' if a.execute else 'plan',example=items[:2]),ensure_ascii=False,indent=2))
if not a.execute:raise SystemExit(0)
if len(items)>1 and not a.all:raise SystemExit('Multiple files selected: explicit --all required')
dest.mkdir(parents=True,exist_ok=True)
if shutil.disk_usage(dest).free<total+1024**3:raise SystemExit('Insufficient disk: require selected bytes plus 1 GiB reserve')
allowed={'ftp.ncbi.nlm.nih.gov','api.gdc.cancer.gov','static-content.springer.com'}
def file_sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
for r in items:
 url=r['url'];assert urllib.parse.urlparse(url).scheme=='https' and urllib.parse.urlparse(url).hostname in allowed
 filename=r['filename'];assert Path(filename).name==filename
 target=dest/(r['file_id'] if a.dataset=='TCGA_BRCA_STAR' else '')/filename;target.parent.mkdir(parents=True,exist_ok=True)
 if target.exists():
  assert file_sha(target)==r['sha256'],'Existing file differs; refusing overwrite'
  print('Verified existing',target);continue
 part=target.with_name(target.name+'.part');assert not part.exists(),'Partial file exists; review it before restarting'
 h=hashlib.sha256();m=hashlib.md5(usedforsecurity=False);size=0
 request=urllib.request.Request(url,headers={'User-Agent':'Workflow3-public-data-fetch/1.0'})
 with urllib.request.urlopen(request,timeout=120) as response,part.open('xb') as f:
  for block in iter(lambda:response.read(4194304),b''):f.write(block);h.update(block);m.update(block);size+=len(block)
 assert size==int(r['bytes']) and h.hexdigest()==r['sha256'],'Size/SHA mismatch: .part retained, not promoted'
 if r.get('md5'):assert m.hexdigest()==r['md5'],'Official/local MD5 mismatch'
 part.rename(target);print('Downloaded and verified',target)
