from pathlib import Path
import csv,json,collections,sys,re
root=Path(__file__).resolve().parents[2];out=root/'docs/audit/2026-10-26';inv=list(csv.DictReader((out/'before/FILE_INVENTORY.tsv').open(encoding='utf-8'),delimiter='\t'))
def write(name,rows,fields):
 with (out/name).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
markers=[]
for r in inv:
 p=root/r['path'];name=p.name.lower()
 if r['role']=='environment_or_dependency' or p.suffix.lower() in {'.py','.log','.html'}:continue
 if not (name.endswith('.complete') or 'completed' in name or name=='final_cleanup_validation.json'):continue
 if int(r['bytes'])>2*1024**2:continue
 status='marker_exists_not_independent_scientific_proof';details=''
 if p.suffix=='.json':
  try:
   d=json.loads(p.read_text(encoding='utf-8-sig'));details=json.dumps({k:v for k,v in d.items() if isinstance(v,(str,int,float,bool,type(None)))},ensure_ascii=False)[:2200]
  except Exception as e:status='json_parse_failed';details=str(e)
 else:details=p.read_text(encoding='utf-8',errors='replace')[:500]
 markers.append(dict(path=r['path'],status=status,details=details))
write('COMPLETION_MARKERS.tsv',markers,['path','status','details'])
links=list(csv.DictReader((out/'before/DOCUMENT_LINKS.tsv').open(encoding='utf-8'),delimiter='\t'));missing=[]
for r in links:
 if r['exists']=='True':continue
 path=r['document'];target=r['target']
 if target.startswith('mailto:'):cat='URI_false_positive_not_local_file'
 elif '[C@H]' in target or '@@' in target:cat='regex_false_positive_in_SMILES_not_a_link'
 elif 'snapshot' in path.lower() or '快照' in path:cat='historical_snapshot_relative_base_or_removed_history'
 elif '/software/' in path:cat='third_party_document_relative_base'
 elif '/sources/DeepPurpose/' in path:cat='third_party_partial_source_snapshot_tutorial_missing'
 else:cat='legacy_document_link_requires_review'
 missing.append(dict(document=path,target=target,classification=cat,action='Original retained; canonical new docs use verified review links'))
write('LEGACY_LINK_TRIAGE.tsv',missing,['document','target','classification','action'])
sizes=collections.Counter()
for r in inv:
 if r['role']=='environment_or_dependency':continue
 suffix=Path(r['path']).suffix.lower()
 if suffix in {'.r','.ps1'} and r['role']=='source_code':sizes[suffix]+=1
R=[str(root/r['path']) for r in inv if r['role']=='source_code' and Path(r['path']).suffix.lower()=='.r']
(out/'R_SCRIPT_PATHS.txt').write_text('\n'.join(R),encoding='utf-8')
write('SCIENCE_FILE_INVENTORY.tsv',[r for r in inv if r['role']!='environment_or_dependency' and not r['path'].startswith(('docs/audit/','scripts/project/'))],list(inv[0]))
print(json.dumps(dict(completion_markers=len(markers),legacy_link_triage=dict(collections.Counter(r['classification'] for r in missing)),script_parse_targets=dict(sizes)),ensure_ascii=False,indent=2))
