"""Engineering acceptance only: no scientific pipelines or statistical refits."""
from pathlib import Path
import csv,json,hashlib,subprocess,sys,shutil,re,ast
root=Path(__file__).resolve().parents[2];out=root/'docs/audit/2026-10-26';checks=[]
def read(path):
 with Path(path).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
def check(name,condition,details=''):
 checks.append(dict(check=name,passed=bool(condition),details=str(details)))
 if not condition:raise RuntimeError(name+': '+str(details))
before=read(out/'before/FILE_INVENTORY.tsv');allowed={'AGENTS.md','工作流3_交接.md','新对话接手_当前状态.md','outputs/project_audit_inventory_2026-10-26.log'};changed=[];missing=[];checked=0
for r in before:
 if r['role']=='environment_or_dependency' or r['path'] in allowed:continue
 p=root/r['path'];checked+=1
 if not p.is_file():missing.append(r['path']);continue
 s=p.stat()
 if s.st_size!=int(r['bytes']) or s.st_mtime_ns!=int(r['mtime_ns']):changed.append(r['path'])
check('Existing_scientific_files_still_present',not missing,str(missing[:10]))
check('Existing_scientific_bytes_mtime_unchanged',not changed,str(changed[:10]))
raw=json.loads((out/'integrity/INTEGRITY_SUMMARY.json').read_text());check('All1214_raw_inputs_hashed',raw['raw_files_hashed']==1214);check('All1208_official_GDC_MD5_match',raw['TCGA_files_rechecked']==1208 and raw['TCGA_MD5_failures']==0)
man=json.loads((out/'STAGE_MANIFEST_SUMMARY.json').read_text());check('No_missing_stage_manifest_files',man['missing']==0);check('No_nonroot_stage_hash_mismatch',man['nonroot_mismatches']==0,'Root-document historical drift retained, not rewritten')
mapping=read(root/'docs/PATH_MAP.tsv')
check('All_review_copies_identical',all(sha(root/r['repository_path'])==r['sha256']==sha(root/r['source']) for r in mapping))
edges=json.loads((root/'config/stages.json').read_text());ids={s['id'] for s in edges};visited=set()
def visit(s,active):
 if s['id'] in visited:return
 if s['id'] in active:raise ValueError('Cyclic stage dependency')
 for dep in s['depends_on']:visit(next(d for d in edges if d['id']==dep),active|{s['id']})
 visited.add(s['id'])
for s in edges:visit(s,set())
check('Stage_graph_11_nodes_acyclic',len(visited)==11)
check('Declared_completion_evidence_exists',all((root/p).exists() for s in edges for p in s['completion_markers']))
links=[]
active_docs=[root/'README.md',root/'AGENTS.md',root/'analysis/README.md',root/'results/README.md',root/'data/README.md',root/'environments/README.md',*list((root/'docs').glob('*.md'))]
for p in active_docs:
 for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
  target=target.strip('<>').split('#')[0]
  if not target or re.match(r'^[a-z]+:',target):continue
  links.append(dict(document=p.relative_to(root).as_posix(),target=target,exists=(p.parent/target).exists()))
check('All_new_canonical_document_links',all(r['exists'] for r in links),str([r for r in links if not r['exists']]))
with (out/'CANONICAL_DOCUMENT_LINKS.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['document','target','exists'],delimiter='\t');w.writeheader();w.writerows(links)
check('Original_R_scripts_parse',all(r['parse_ok']=='TRUE' for r in read(out/'R_PARSE.tsv')))
check('Original_PS_scripts_parse',all(r['parse_ok']=='True' for r in read(out/'POWERSHELL_PARSE.tsv')))
for p in [*list((root/'scripts/project').glob('*.py')),*list((root/'tests').glob('*.py'))]:ast.parse(p.read_text(encoding='utf-8'))
check('New_Python_syntax',True)
runtime=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=root,capture_output=True,text=True,encoding='utf-8');(out/'NEW_TOOL_TESTS.txt').write_text(runtime.stdout+runtime.stderr,encoding='utf-8');check('Five_restore_safety_tests',runtime.returncode==0 and 'Ran 5 tests' in runtime.stderr)
git=Path(shutil.which('git') or '');bash=git.parent.parent/'bin/bash.exe';bash_status='not_available'
record=out/'BASH_SYNTAX.json'
if record.exists():
 b=json.loads(record.read_text(encoding='utf-8-sig'));check('Cloud_setup_shell_syntax',b['returncode']==0 and b['sha256']==sha(root/'scripts/cloud/setup.sh'));bash_status=b['check']
elif bash.is_file():
 r=subprocess.run([str(bash),'--noprofile','--norc','-n','scripts/cloud/setup.sh'],cwd=root,capture_output=True,text=True);check('Cloud_setup_shell_syntax',r.returncode==0,r.stderr);bash_status='Git Bash syntax only; Linux runtime NOT validated'
# Data-free review clone in a new ignored audit directory; originals are never moved.
r=subprocess.run([sys.executable,str(root/'scripts/project/audit_git_payload.py')],cwd=root,capture_output=True,text=True,encoding='utf-8');check('Git_payload_policy_tests',r.returncode==0,r.stderr)
ready=json.loads((out/'GITHUB_READINESS.json').read_text());check('No_oversized_Git_candidates',ready['candidate_files_over_100MiB']==0);check('No_credential_pattern_flags',ready['credential_pattern_flags']==0);check('Project_Git_not_initialized',not ready['project_git_repository'])
clone=out/'after/cloud review 克隆';check('Clone_fixture_confined',clone.resolve().is_relative_to((out/'after').resolve()))
clone.mkdir(exist_ok=True)
for item in read(out/'GITHUB_CANDIDATES.tsv'):
 p=root/item['path'];dest=clone/item['path'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
clone_results=[]
base_python=Path(sys.base_prefix)/'python.exe' if sys.platform=='win32' else Path(sys.base_prefix)/'bin/python3'
clone_python=str(base_python) if base_python.is_file() else sys.executable
for command in [['scripts/project/doctor.py','--evidence-only'],['scripts/project/verify_project.py'],['scripts/project/restore_legacy_layout.py','--dry-run']]:
 r=subprocess.run([clone_python,str(clone/command[0]),'--root',str(clone),*command[1:]],cwd=clone,capture_output=True,text=True,encoding='utf-8');check('Data_free_clone_'+Path(command[0]).stem,r.returncode==0,r.stdout+r.stderr);clone_results.append(dict(command=command,python=clone_python,stdout=r.stdout,stderr=r.stderr))
(out/'DATA_FREE_CLONE_TEST.json').write_text(json.dumps(clone_results,ensure_ascii=False,indent=2),encoding='utf-8')
with (out/'FINAL_VERIFICATION.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['check','passed','details'],delimiter='\t');w.writeheader();w.writerows(checks)
summary=dict(status='complete',scope='engineering_audit_and_repository_preparation',scientific_analysis_rerun=False,originals_deleted=0,originals_moved=0,baseline_scientific_files_metadata_checked=checked,unexpected_modified_files=changed,missing_original_files=missing,raw_files_hashed=1214,TCGA_MD5_failures=0,stage_manifest_records=man['records'],nonroot_stage_hash_mismatches=0,historical_mutable_root_document_drift=man['root_document_drift'],review_copies=len(mapping),review_copy_bytes=sum(int(r['bytes']) for r in mapping),canonical_links_checked=len(links),checks_passed=len(checks),data_free_review_clone_passed=True,linux_full_reproduction_verified=False,cloud_setup_syntax=bash_status,project_Git_initialized=False,git_commit=False,git_push=False,scientific_stop_route_preserved=True,remaining_blockers=['Project repository and Cloud environment not configured','Large derived data external storage not designated','Linux analysis environments and legacy paths not fully ported','Third-party literature/weight redistribution licenses not fully verified'])
(out/'AUDIT_COMPLETED.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
