"""Restore only hash-verified small review files; never overwrite different content."""
from pathlib import Path
import argparse,csv,hashlib,shutil,json
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);mode=p.add_mutually_exclusive_group();mode.add_argument('--apply',action='store_true');mode.add_argument('--dry-run',action='store_true');a=p.parse_args();root=a.root.resolve();records=list(csv.DictReader((root/'docs/PATH_MAP.tsv').open(encoding='utf-8'),delimiter='\t'));counts=dict(already_present=0,would_restore=0,restored=0);plans=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for r in records:
 source=(root/r['repository_path']).resolve();target=(root/r['source']).resolve()
 if not(source.is_relative_to(root) and target.is_relative_to(root) and r['source'].startswith('outputs/')):raise ValueError('Path outside expected project outputs')
 if not(source.is_file() and sha(source)==r['sha256']):raise ValueError('Review copy changed: '+r['repository_path'])
 if target.exists():
  if sha(target)!=r['sha256']:raise ValueError('Refuse overwrite: '+r['source'])
  counts['already_present']+=1;continue
 counts['would_restore']+=1
 plans.append((source,target))
if a.apply:
 for source,target in plans:
  target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);counts['restored']+=1
print(json.dumps(dict(mode='apply' if a.apply else 'dry-run',**counts),indent=2))
