"""Use an ignored scratch Git directory to test rules; never init/stage/commit the project."""
from pathlib import Path
import subprocess,csv,json,re,os
root=Path(__file__).resolve().parents[2];out=root/'docs/audit/2026-10-26';scratch=out/'after/ignore-test.git';scratch.parent.mkdir(parents=True,exist_ok=True)
if not scratch.exists():subprocess.run(['git','init','--bare','--quiet',str(scratch)],check=True)
git=['git','--git-dir='+str(scratch),'--work-tree='+str(root),'-c','core.bare=false']
def ignored(paths):
 r=subprocess.run(git+['check-ignore','--no-index','--stdin','-z'],input=('\0'.join(paths)+'\0').encode('utf-8'),capture_output=True)
 if r.returncode not in [0,1]:raise RuntimeError(r.stderr.decode(errors='replace'))
 return set(r.stdout.decode('utf-8').strip('\0').split('\0')) if r.stdout else set()
samples={'outputs/test.rds':True,'原始数据/test.tsv':True,'data/external/test.tsv':True,'environments/local_snapshots/programs.freeze.txt':True,'docs/audit/2026-10-26/before/FILE_INVENTORY.tsv':True,'.env':True,'private.pem':True,'analysis/10_programs/scripts/program_analysis.py':False,'results/10_programs/frozen/Programs_v1.json':False,'docs/RESEARCH_PROGRESS.md':False,'README.md':False,'requirements.txt':False}
seen=ignored(list(samples));assert all((p in seen)==expected for p,expected in samples.items())
def write(name,rows,fields):
 with (out/name).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
write('GITIGNORE_TESTS.tsv',[dict(path=p,expected_ignored=v,actual_ignored=p in seen) for p,v in samples.items()],['path','expected_ignored','actual_ignored'])
paths=[]
for base,dirs,files in os.walk(root):
 # Prune the large namespaces that are explicitly ignored; do not inspect dependency contents.
 dirs[:]=[d for d in dirs if d not in ['outputs','原始数据','env','.venv','__pycache__','.git','node_modules','local_snapshots','before','after']]
 for name in files:paths.append((Path(base)/name).relative_to(root).as_posix())
skip=ignored(paths);candidates=[];secrets=[]
patterns={'private_key':r'-----BEGIN [A-Z ]*PRIVATE KEY-----','github_token':r'\bgh[pousr]_[A-Za-z0-9]{30,}\b','openai_like_token':r'\bsk-[A-Za-z0-9_-]{24,}\b','credential_literal':r'''(?i)(?:api_key|access_token|auth_token|client_secret|secret_access_key|password)\s*[=:]\s*["'][^"'\n]{20,}["']'''}
for rel in paths:
 if rel in skip:continue
 p=root/rel;size=p.stat().st_size;candidates.append(dict(path=rel,bytes=size,over_50MiB=size>50*1024**2,over_100MiB=size>100*1024**2,publication_note='Prepared only; no upload or Git staging'))
 if size<=5*1024**2 and p.suffix.lower() in ['.py','.r','.ps1','.sh','.md','.json','.txt','.yml','.yaml','.toml','.ini']:
  try:text=p.read_text(encoding='utf-8-sig')
  except UnicodeError:continue
  for i,line in enumerate(text.splitlines(),1):
   for kind,pat in patterns.items():
    if re.search(pat,line) and not any(v in line.lower() for v in ['your_api_key','example_token','placeholder','credential_literal']):secrets.append(dict(path=rel,line=i,kind=kind,value='REDACTED_NOT_CAPTURED',status='Manual review required; no raw credential retained'))
write('GITHUB_CANDIDATES.tsv',sorted(candidates,key=lambda r:r['path']),['path','bytes','over_50MiB','over_100MiB','publication_note'])
write('SECRET_PATTERN_REVIEW.tsv',secrets,['path','line','kind','value','status'])
r=subprocess.run(['git','rev-parse','--show-toplevel'],cwd=root,capture_output=True,text=True)
summary=dict(project_git_repository=r.returncode==0,project_git_root=r.stdout.strip() if r.returncode==0 else None,project_git_initialized=False,scratch_git_metadata_created_for_tests=True,gitignore_cases_passed=len(samples),candidate_files=len(candidates),candidate_bytes=sum(r['bytes'] for r in candidates),candidate_files_over_100MiB=sum(r['over_100MiB'] for r in candidates),credential_pattern_flags=len(secrets),secret_scan_scope='Eligible small text files only; pattern scan is not a full security or privacy assessment',git_add=False,git_commit=False,git_push=False)
(out/'GITHUB_READINESS.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
