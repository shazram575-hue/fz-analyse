"""Copy a size-bounded, hash-traceable review view; legacy scientific files stay in place."""
from pathlib import Path
import csv,json,hashlib,shutil,sys,re
root=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_tsv(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
specs=[
 ('01_disease','第一阶段_2026-10-23',[],['audit/PHASE1_COMPLETED.txt'],['第一阶段_结果报告.md','第一阶段_分析方案与验收.md']),
 ('02_cell_atlas','单细胞基础补全_2026-10-23',['01_disease'],['audit/BASICS_COMPLETED.txt'],['单细胞基础补全_结果报告.md']),
 ('03_deeppurpose','DeepPurpose_预测_2026-10-24',['01_disease','02_cell_atlas'],['audit/PREDICTION_STAGE_COMPLETED.txt'],[]),
 ('04_proteome','身份深核与独立交集_2026-10-24',['03_deeppurpose'],['audit/INDEPENDENT_INTERSECTION_STAGE_COMPLETED.txt'],[]),
 ('05_pubmed','docking_audit_2026-10-24',['04_proteome'],['pubmed_comprehensive_recheck/COMPREHENSIVE_SEARCH_COMPLETED.json'],[]),
 ('06_docking47','reverse_docking47_2026-10-24',['04_proteome'],['audit/ALL47_DOCKING_DELIVERY_COMPLETED.json','audit/ALL47_TECHNICAL_CHECKS_COMPLETED.json'],[]),
 ('07_candidate','候选推进讨论_2026-10-25',['05_pubmed','06_docking47'],[],[]),
 ('08_expression','MFSD10_表达分布摸底_2026-10-26',['01_disease','02_cell_atlas'],['audit/EXPRESSION_SURVEY_COMPLETED.json'],[]),
 ('09_feasibility','恶性上皮_表达程序可行性_2026-10-26',['08_expression'],['audit/FEASIBILITY_COMPLETED.json'],[]),
 ('10_programs','MFSD10_程序发现与独立验证_2026-10-26',['09_feasibility','01_disease'],['audit/PIPELINE_COMPLETED.json'],[]),
 ('11_aggregation','MFSD10_稀疏与聚合探索_2026-10-26',['10_programs'],['audit/COMPLETED.json','audit/CLOSEOUT_COMPLETED.json','audit/STOP_THIS_ROUTE.json'],[])]
config=[];copies=[];skipped=[]
for sid,folder,deps,markers,reports in specs:
 base=root/'outputs'/folder
 stage=dict(id=sid,legacy_directory='outputs/'+folder,depends_on=deps,completion_markers=[('outputs/'+folder+'/'+p) for p in markers],computational_state='user_selection_record' if sid=='07_candidate' else 'completed_evidence_available',scientific_state='Not a validated MFSD10 causal mechanism',legacy_entrypoints_ported_to_linux=False)
 config.append(stage)
 selected=[]
 for p in base.rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(base);parts={s.lower() for s in rel.parts}
  if parts&{'env','cache','__pycache__','cnmf','sensitivity','models','package_cache','public_reference_code','public_reference','model_cache','inference_cache','prepared','poses','docking','controls','sources','author_reference'}:continue
  size=p.stat().st_size;suffix=p.suffix.lower();kind=None
  if rel.parts[0]=='scripts' and suffix in {'.py','.r','.ps1','.sh','.txt'}:kind='code'
  elif len(rel.parts)==1 and suffix=='.md':kind='report'
  elif rel.parts[0] in {'results','frozen','main_input','discovery','validation'} and suffix in {'.tsv','.csv','.json','.py','.md'}:kind='evidence'
  elif rel.parts[0]=='metadata' and suffix in {'.tsv','.json'} and any(s in p.name.lower() for s in ['patient','library','gene_mapping','scope','provenance']):kind='evidence'
  elif sid=='01_disease' and rel.as_posix() in {'metadata/TCGA_local_file_ids.tsv','metadata/TCGA_verified_file_sample_metadata.tsv','metadata/TCGA_model_sample_order.tsv','metadata/GSE161529_author_QC_and_analysis_membership.tsv','metadata/GSE161529_features.tsv'}:kind='input_provenance'
  elif rel.parts[0] in {'structure_verification_CID209254','structure_comparison'} and suffix in {'.json','.tsv','.sdf','.mol','.smi','.txt','.png','.svg','.md'}:kind='input_provenance'
  elif sid=='06_docking47' and rel.parts[0]=='inputs' and len(rel.parts)==2 and suffix in {'.json','.tsv','.sdf','.mol','.smi','.pdbqt'}:kind='input_provenance'
  elif rel.parts[0]=='figures' and suffix in {'.png','.svg'} and len(rel.parts)==2:kind='figure'
  elif rel.parts[0]=='audit' and (rel.as_posix() in markers or any(x in p.name for x in ['validation','Validation','COMPLETED','COMPLETION','STOP_','Summary'])) and suffix in {'.tsv','.json','.txt'}:kind='audit'
  elif rel.as_posix() in {'audit/方法实现与核读说明.md','audit/图形与收尾说明.md','audit/Figure_visual_QA.md'}:kind='method_limitations'
  elif rel.parts[0] in {'pubmed_comprehensive_recheck','pubmed_zero_recheck'} and suffix in {'.tsv','.json','.md'} and len(rel.parts)==2:kind='evidence'
  elif sid=='03_deeppurpose' and rel.as_posix()=='predictions/Frozen_prediction_design.json':kind='evidence'
  elif sid=='04_proteome' and rel.as_posix()=='background/Frozen_independent_drug_side_protocol.json':kind='evidence'
  elif sid=='06_docking47' and rel.as_posix()=='inputs/47_actual_search_regions.tsv':kind='evidence'
  if not kind:continue
  if size>2*1024**2:skipped.append(dict(path=p.relative_to(root).as_posix(),bytes=size,reason='Review export >2 MiB; kept in legacy directory'));continue
  selected.append((p,rel,kind))
 for p,rel,kind in selected:
  target=root/('analysis' if kind=='code' else 'results')/sid/rel;target.parent.mkdir(parents=True,exist_ok=True)
  h=sha(p)
  if target.exists():assert sha(target)==h,'Refuse to overwrite different review copy: '+str(target)
  else:shutil.copy2(p,target)
  copies.append(dict(stage=sid,kind=kind,source=p.relative_to(root).as_posix(),repository_path=target.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=h,copy_identical=sha(target)==h))
 # Root-level original reports remain in place and are referenced, not moved.
 stage['root_reports']=reports
 stage['review_directory']='results/'+sid
 stage['code_directory']='analysis/'+sid+'/scripts'
(root/'config').mkdir(exist_ok=True)
(root/'config/stages.json').write_text(json.dumps(config,ensure_ascii=False,indent=2),encoding='utf-8')
(root/'config/project_paths.json').write_text(json.dumps(dict(project_root_environment='BIOINFO_PROJECT_ROOT',external_data_environment='BIOINFO_DATA_ROOT',default_external_data='data/external',legacy_raw_directory='原始数据',legacy_output_directory='outputs',review_code_directory='analysis',review_evidence_directory='results',read_only_originals=True),ensure_ascii=False,indent=2),encoding='utf-8')
write_tsv(root/'docs/PATH_MAP.tsv',copies,['stage','kind','source','repository_path','bytes','sha256','copy_identical'])
write_tsv(root/'docs/audit/2026-10-26/EXPORT_EXCLUSIONS.tsv',skipped,['path','bytes','reason'])
# Publish only names/versions; keep absolute local runtime paths and complete raw snapshots local.
env=root/'environments';snap=env/'local_snapshots'
for name in ['programs','deeppurpose','docking']:
 packages=json.loads((snap/(name+'.pip_list.json')).read_text(encoding='utf-8'))
 write_tsv(env/(name+'.installed.tsv'),packages,['name','version'])
 runtime=json.loads((snap/(name+'.runtime.json')).read_text(encoding='utf-8'));runtime.pop('executable',None)
 (env/(name+'.runtime.json')).write_text(json.dumps(runtime,indent=2),encoding='utf-8')
ip=list(csv.DictReader((snap/'R_installed_packages_local.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
if ip:write_tsv(env/'R-installed-observed.tsv',[{k:r.get(k,'') for k in ['Package','Version','Built','Priority','NeedsCompilation']} for r in ip],['Package','Version','Built','Priority','NeedsCompilation'])
else:assert (env/'R-installed-observed.tsv').exists(),'R enumeration empty: first run capture_R_descriptions.py; do not invent a dependency lock'
shutil.copyfile(snap/'R_sessionInfo.txt',env/'R-sessionInfo-observed.txt')
pk={r['name'].lower():r['version'] for r in json.loads((snap/'programs.pip_list.json').read_text())}
core=['numpy','scipy','pandas','scikit-learn','anndata','statsmodels','matplotlib','h5py','cnmf','scanpy','pyyaml']
req=['# Actual observed Python-program environment versions on Windows, 2026-10-26.','# Not a complete lockfile. Linux resolution and binary compatibility are NOT yet validated.','# DeepPurpose and docking use separate environments; do not merge them into this one.']+[f'{p}=={pk[p]}' for p in core]
(root/'requirements.txt').write_text('\n'.join(req)+'\n',encoding='utf-8')
print(json.dumps(dict(copied_files=len(copies),copied_bytes=sum(r['bytes'] for r in copies),excluded_large_review_files=len(skipped),stages=len(config),originals_moved=0),indent=2))
