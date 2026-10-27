"""Finish documentation bookkeeping after engineering acceptance; no scientific execution."""
from pathlib import Path
import csv,json,hashlib,subprocess,sys,re,collections
root=Path(__file__).resolve().parents[2];out=root/'docs/audit/2026-10-26'
summary=json.loads((out/'AUDIT_COMPLETED.json').read_text());assert summary['status']=='complete' and summary['scientific_analysis_rerun'] is False
before={r['path']:r for r in csv.DictReader((out/'before/FILE_INVENTORY.tsv').open(encoding='utf-8'),delimiter='\t')}
after=list(csv.DictReader((out/'after/inventory/FILE_INVENTORY.tsv').open(encoding='utf-8'),delimiter='\t'));changes=[]
for r in after:
 old=before.get(r['path'])
 if old is None:action='created'
 elif int(r['bytes'])!=int(old['bytes']) or int(r['mtime_ns'])!=int(old['mtime_ns']):action='modified'
 else:continue
 if r['role']=='environment_or_dependency':continue
 category='ignored_audit_fixture_or_inventory' if '/after/' in r['path'] or '/before/' in r['path'] else 'project_preparation_or_documentation'
 changes.append(dict(path=r['path'],action=action,bytes=r['bytes'],category=category))
with (out/'CREATED_MODIFIED_FILES.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['path','action','bytes','category'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(changes)
links=[]
for p in [root/'README.md',root/'AGENTS.md',*list((root/'docs').glob('*.md'))]:
 for s in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
  s=s.strip('<>').split('#')[0]
  if not s or re.match(r'^[a-z]+:',s):continue
  links.append(dict(document=p.relative_to(root).as_posix(),target=s,exists=(p.parent/s).exists()))
assert all(r['exists'] for r in links),[r for r in links if not r['exists']]
with (out/'FINAL_DOCUMENT_LINKS.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['document','target','exists'],delimiter='\t');w.writeheader();w.writerows(links)
manifest=[]
for name in ['analysis','results','config','scripts','environments','docs','data','tests']:
 for p in (root/name).rglob('*'):
  if not p.is_file() or any(x in p.parts for x in ['__pycache__','local_snapshots','before','after']):continue
  if p.name in ['REPOSITORY_FILE_SHA256_LOCAL.tsv','GITHUB_CANDIDATES.tsv','GITHUB_READINESS.json','AUDIT_CLOSEOUT_COMPLETED.json']:continue
  manifest.append(dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
with (out/'REPOSITORY_FILE_SHA256_LOCAL.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(manifest)
r=subprocess.run([sys.executable,str(root/'scripts/project/audit_git_payload.py')],cwd=root,capture_output=True,text=True,encoding='utf-8');assert r.returncode==0,r.stderr
ready=json.loads((out/'GITHUB_READINESS.json').read_text());assert not ready['project_git_repository'] and ready['candidate_files_over_100MiB']==0 and ready['credential_pattern_flags']==0
result=dict(status='complete',handoff_section=39,new_scientific_analysis=False,final_canonical_document_links=len(links),review_copy_count=summary['review_copies'],review_copy_bytes=summary['review_copy_bytes'],candidate_files=ready['candidate_files'],candidate_bytes_at_listing=ready['candidate_bytes'],git_push=False,linux_full_reproduction_verified=False,stop_route_unchanged=True,created_modified_list_scope='Before/after inventory snapshots; final bookkeeping itself not recursively self-enumerated',local_checksum_scope='Project preparation view, excludes ignored caches/local snapshots; immutable review payload uses PATH_MAP cross-platform hashes')
(out/'AUDIT_CLOSEOUT_COMPLETED.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
