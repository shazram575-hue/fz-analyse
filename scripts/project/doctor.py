"""Portable, read-only resource and evidence doctor; never starts an analysis."""
from pathlib import Path
import argparse,json,os,platform,shutil,sys,importlib.util,csv
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--evidence-only',action='store_true');a=p.parse_args();root=a.root.resolve()
required=['AGENTS.md','README.md','docs/RESEARCH_PROGRESS.md','docs/PROJECT_AUDIT.md','docs/DATA_MANIFEST.md','docs/CLOUD_MIGRATION.md','docs/PATH_MAP.tsv','config/stages.json','config/public_datasets.json']
missing=[v for v in required if not(root/v).is_file()]
stages=json.loads((root/'config/stages.json').read_text()) if(root/'config/stages.json').exists() else []
raw=root/'原始数据';large_objects=[root/'outputs/MFSD10_程序发现与独立验证_2026-10-26/objects'/n for n in ['GSE161529_confirmed_malignant_singlets_counts.h5ad','GSE176078_malignant_candidates_author_counts_clinical_v2.h5ad']]
ram=None
if Path('/proc/meminfo').exists():
 for line in Path('/proc/meminfo').read_text().splitlines():
  if line.startswith('MemAvailable:'):ram=int(line.split()[1])*1024
report=dict(project_root=str(root),python=sys.version.split()[0],platform=platform.system(),available_RAM_bytes=ram,disk_free_bytes=shutil.disk_usage(root).free,missing_review_files=missing,stage_count=len(stages),raw_directory_present=raw.is_dir(),large_program_objects_present=[p.is_file() for p in large_objects],Rscript_on_PATH=shutil.which('Rscript'),modules={m:bool(importlib.util.find_spec(m)) for m in ['numpy','scipy','pandas','sklearn','anndata','statsmodels','cnmf','torch','rdkit','DeepPurpose']},mode='evidence-only' if a.evidence_only else 'inspection',scientific_run_attempted=False)
cg=Path('/sys/fs/cgroup/memory.max')
report['cgroup_memory_limit_bytes']=int(cg.read_text().strip()) if cg.exists() and cg.read_text().strip().isdigit() else None
report['review_ready']=not missing and len(stages)==11
report['linux_full_reproduction_verified']=False
report['note']='Review can pass without raw data or R. Review readiness is not scientific runtime readiness. Read CLOUD_MIGRATION before any large run.'
print(json.dumps(report,ensure_ascii=False,indent=2));sys.exit(0 if report['review_ready'] else 1)
