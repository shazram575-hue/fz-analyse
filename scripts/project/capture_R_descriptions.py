"""Fallback for sandbox R enumeration failures; metadata is NOT a load test."""
from pathlib import Path
import argparse,csv,json
p=argparse.ArgumentParser();p.add_argument('--library',action='append',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();rows=[]
for library in a.library:
 for fn in sorted(Path(library).glob('*/DESCRIPTION')):
  fields={};key=None
  for line in fn.read_text(encoding='utf-8',errors='replace').splitlines():
   if line.startswith((' ','\t')) and key:fields[key]+=' '+line.strip()
   elif ':' in line:key,value=line.split(':',1);fields[key]=value.strip()
  rows.append({k:fields.get(k,'') for k in ['Package','Version','Built','NeedsCompilation','License','Repository']}|{'LibPath':library,'verification':'DESCRIPTION metadata only; not package load validation'})
assert rows,'No package descriptions found; do not create a fake lockfile'
with a.output.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print(json.dumps(dict(description_records=len(rows),libraries=a.library,package_load_verified=False),indent=2))
