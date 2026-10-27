"""Stream raw-input hashes; inspect small-file duplicates; never edit scientific files."""
from pathlib import Path
import csv,hashlib,json,argparse,collections,time
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--inventory',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root.resolve();out=a.output;out.mkdir(parents=True,exist_ok=True)
files=list(csv.DictReader(a.inventory.open(encoding='utf-8'),delimiter='\t'));hashes={}
def digest(path,md5=False):
 h=hashlib.sha256();m=hashlib.md5(usedforsecurity=False) if md5 else None
 with path.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b);m.update(b) if m else None
 return h.hexdigest(),m.hexdigest() if m else ''
def write(name,rows,fields):
 with (out/name).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
gdc={}
with (root/'outputs/第一阶段_2026-10-23/metadata/TCGA_verified_file_sample_metadata.tsv').open(encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):gdc[r['file_id']]=r
raw=[];downloads=[]
for n,r in enumerate(f for f in files if f['role']=='original_data'):
 path=root/r['path'];s,m=digest(path,True);hashes[r['path']]=s
 fid=path.parent.name;g=gdc.get(fid);ok=(m.lower()==g['gdc_md5'].lower() and path.stat().st_size==int(g['gdc_bytes'])) if g else None
 raw.append(dict(path=r['path'],bytes=path.stat().st_size,sha256=s,md5=m,gdc_file_id=fid if g else '',official_gdc_md5_matches=ok if g else 'not_applicable'))
 if g:downloads.append(dict(dataset='TCGA_BRCA_STAR',file_id=fid,filename=path.name,url='https://api.gdc.cancer.gov/data/'+fid,bytes=path.stat().st_size,md5=g['gdc_md5'],sha256=s,legacy_relative_path=r['path']))
 if n%200==0:print('RAW_HASH_PROGRESS',n+1,flush=True)
write('RAW_INPUT_CHECKSUMS.tsv',raw,['path','bytes','sha256','md5','gdc_file_id','official_gdc_md5_matches'])
write('TCGA_DOWNLOAD_MANIFEST.tsv',downloads,['dataset','file_id','filename','url','bytes','md5','sha256','legacy_relative_path'])
same_size=collections.defaultdict(list)
for r in files:
 if r['role']=='environment_or_dependency' or int(r['bytes'])==0 or int(r['bytes'])>50*1024**2:continue
 same_size[int(r['bytes'])].append(r)
dups=collections.defaultdict(list)
for size,group in same_size.items():
 if len(group)<2:continue
 for r in group:
  if r['path'] not in hashes:hashes[r['path']]=digest(root/r['path'])[0]
  dups[(size,hashes[r['path']])].append(r['path'])
dup_rows=[]
for (size,h),paths in dups.items():
 if len(paths)>1:
  for path in paths:dup_rows.append(dict(group_sha256=h,bytes=size,copies=len(paths),path=path,action='Retained; no deduplication or deletion'))
write('DUPLICATE_SMALL_FILES.tsv',dup_rows,['group_sha256','bytes','copies','path','action'])
summary=dict(raw_files_hashed=len(raw),raw_bytes=sum(r['bytes'] for r in raw),TCGA_files_rechecked=len(downloads),TCGA_MD5_failures=sum(r['official_gdc_md5_matches'] is False for r in raw),small_duplicate_groups=sum(len(v)>1 for v in dups.values()),small_duplicate_members=len(dup_rows),duplicate_scope='Non-environment files <=50 MiB; large-file duplication not verified',original_files_deleted=0)
(out/'INTEGRITY_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
