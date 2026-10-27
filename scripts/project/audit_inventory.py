"""Read-only recursive inventory; outputs only to the requested audit directory."""
from pathlib import Path
import argparse,csv,json,os,re,collections,hashlib,ast,shutil,subprocess,sys

PARSER=argparse.ArgumentParser();PARSER.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);PARSER.add_argument('--output',type=Path,required=True);args=PARSER.parse_args()
root=args.root.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
ENV_PARTS={'env','.venv','venv','site-packages','node_modules','__pycache__','.git'}
CODE_EXT={'.py','.r','.ps1','.sh','.bash','.bat','.cmd','.ipynb'}
def role(rel):
 parts=rel.parts;low={p.lower() for p in parts}
 if low&ENV_PARTS:return 'environment_or_dependency'
 if parts[0]=='原始数据':return 'original_data'
 if 'scripts' in low or rel.suffix.lower() in CODE_EXT:return 'source_code'
 if low&{'cache','r_tmp','tmp','prepared','prepared_by_library','canonical_by_library','epithelial_by_library','cnmf_tmp','model_cache','package_cache','pip','numba','predictions','inference_cache'}:return 'cache_or_intermediate'
 if low&{'models','model_audit','model_cache'} or rel.suffix.lower() in {'.pt','.pth','.ckpt'}:return 'model_weights_or_artifact'
 if 'objects' in low or rel.suffix.lower() in {'.rds','.h5ad','.h5','.hdf5','.loom'}:return 'processed_object'
 if 'figures' in low or rel.suffix.lower() in {'.png','.svg','.pdf','.jpg'}:return 'figure_or_document'
 if 'logs' in low or rel.suffix.lower()=='.log':return 'log'
 if 'audit' in low:return 'audit_evidence'
 if low&{'results','main_input','main_input1881'}:return 'result_or_input_table'
 if 'metadata' in low:return 'metadata'
 if rel.suffix.lower() in {'.md','.html'}:return 'documentation'
 return 'other'
files=[];dirs=[];errors=[];linked=[]
for base,subdirs,names in os.walk(root,followlinks=False):
 base=Path(base)
 for d in list(subdirs):
  p=base/d
  if p.is_symlink() or (hasattr(p,'is_junction') and p.is_junction()):linked.append(str(p.relative_to(root)));subdirs.remove(d)
 dirs.append(base.relative_to(root).as_posix())
 for name in names:
  p=base/name;rel=p.relative_to(root)
  try:
   s=p.stat();files.append(dict(path=rel.as_posix(),bytes=s.st_size,mtime_ns=s.st_mtime_ns,role=role(rel),stage=rel.parts[1] if rel.parts[0]=='outputs' and len(rel.parts)>1 else rel.parts[0],extension=p.suffix.lower()))
  except OSError as e:errors.append(dict(path=rel.as_posix(),error=str(e)))
def tsv(name,rows,fields):
 with (out/name).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
tsv('FILE_INVENTORY.tsv',files,['path','bytes','mtime_ns','role','stage','extension'])
by=collections.defaultdict(lambda:[0,0])
for f in files:by[(f['stage'],f['role'])][0]+=1;by[(f['stage'],f['role'])][1]+=f['bytes']
tsv('STAGE_STORAGE.tsv',[dict(stage=k[0],role=k[1],files=v[0],bytes=v[1]) for k,v in sorted(by.items())],['stage','role','files','bytes'])
tsv('LARGE_FILES.tsv',sorted([f for f in files if f['bytes']>=5*1024**2 and f['role']!='environment_or_dependency'],key=lambda f:-f['bytes']),['path','bytes','mtime_ns','role','stage','extension'])
science=[f for f in files if f['role']!='environment_or_dependency' and not f['path'].startswith(out.relative_to(root).as_posix()+'/')]
scripts=[f for f in science if Path(f['path']).suffix.lower() in CODE_EXT and f['role'] not in ['cache_or_intermediate','original_data']]
findings=[];imports=[];links=[];references=[]
patterns={'windows_drive':r'(?i)(?<![a-z])[a-z]:[\\/][^\s\"\'\)\];,]+','windows_R_library':r'\.libPaths\s*\(', 'powershell_or_windows_executable':r'(?i)(Rscript\.exe|python\.exe|\.dll\b|Get-CimInstance|Start-Process|cmd /c)', 'explicit_random_seed':r'(?i)(set\.seed|random_state|random_seed|--seed|seed\s*=)', 'integrated_expression_mention':r'(?i)(assay\s*=\s*[\"\']integrated|X_umap|IntegrateData|Harmony)', 'shell_copy_in_python':r'subprocess\.(?:call|run)\(.{0,100}shell\s*=\s*True'}
for rec in scripts:
 p=root/rec['path']
 if rec['bytes']>10*1024**2:continue
 try:text=p.read_text(encoding='utf-8-sig')
 except UnicodeError:
  try:text=p.read_text(encoding='cp1252')
  except UnicodeError:errors.append(dict(path=rec['path'],error='Could not decode script; not inspected semantically'));continue
 for i,line in enumerate(text.splitlines(),1):
  for kind,pat in patterns.items():
   if re.search(pat,line):findings.append(dict(path=rec['path'],line=i,kind=kind,excerpt=line[:360]))
  # Literal stage-relative paths provide hints, not proof of dynamic dependencies.
  for match in re.finditer(r'''["']([^"'\n]{1,260}\.(?:tsv|csv|rds|h5ad|mtx|json|py|R|ps1|txt|gz))["']''',line):
   s=match.group(1)
   if '/' in s or '\\' in s:references.append(dict(script=rec['path'],line=i,literal=s,note='Static literal; resolution may depend on root/out variables'))
 if p.suffix.lower()=='.py':
  try:
   tree=ast.parse(text,filename=rec['path'])
   for n in ast.walk(tree):
    if isinstance(n,ast.Import):imports.extend(dict(path=rec['path'],module=a.name.split('.')[0]) for a in n.names)
    elif isinstance(n,ast.ImportFrom) and n.module:imports.append(dict(path=rec['path'],module=n.module.split('.')[0]))
  except (SyntaxError,ValueError) as e:findings.append(dict(path=rec['path'],line=getattr(e,'lineno',0),kind='python_parse_error',excerpt=str(e)))
 elif p.suffix.lower()=='.r':
  imports.extend(dict(path=rec['path'],module=m.group(1)) for m in re.finditer(r'(?:library|require)\s*\(\s*[\"\']?([\w.]+)',text))
for rec in science:
 p=root/rec['path']
 if p.suffix.lower()!='.md' or rec['bytes']>2*1024**2:continue
 try:text=p.read_text(encoding='utf-8-sig')
 except UnicodeError:continue
 for match in re.finditer(r'\]\(([^)]+)\)',text):
  s=match.group(1).strip('<>');s=s.split('#')[0]
  if not s or re.match(r'^[a-z]+://',s):continue
  target=Path(s) if re.match(r'^[A-Za-z]:',s) else p.parent/s
  links.append(dict(document=rec['path'],target=s,exists=target.exists(),note='Historical missing links require stage review; not auto-repaired'))
tsv('SCRIPT_INVENTORY.tsv',scripts,['path','bytes','mtime_ns','role','stage','extension'])
tsv('STATIC_FINDINGS.tsv',findings,['path','line','kind','excerpt'])
tsv('SCRIPT_IMPORTS.tsv',imports,['path','module'])
tsv('SCRIPT_PATH_LITERALS.tsv',references,['script','line','literal','note'])
tsv('DOCUMENT_LINKS.tsv',links,['document','target','exists','note'])
summary=dict(root=str(root),file_count=len(files),directory_count=len(dirs),bytes_total=sum(f['bytes'] for f in files),science_files=len(science),bytes_excluding_environments=sum(f['bytes'] for f in science),original_data_files=sum(f['role']=='original_data' for f in files),original_data_bytes=sum(f['bytes'] for f in files if f['role']=='original_data'),source_script_count=len(scripts),large_files_over_100MiB=sum(f['bytes']>100*1024**2 and f['role']!='environment_or_dependency' for f in files),hardcoded_path_scripts=len(set(r['path'] for r in findings if r['kind']=='windows_drive')),python_parse_errors=sum(r['kind']=='python_parse_error' for r in findings),broken_markdown_links=sum(not r['exists'] for r in links),disk_free_bytes=shutil.disk_usage(root).free,skipped_linked_directories=linked,errors=errors)
(out/'INVENTORY_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
