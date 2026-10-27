"""Evidence-only checks that also work in a data-free repository clone."""
from pathlib import Path
import argparse,json,csv,hashlib,ast,sys
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);a=p.parse_args();root=a.root.resolve()
with (root/'docs/PATH_MAP.tsv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
 p=root/r['repository_path'];assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['repository_path']
stages=json.loads((root/'config/stages.json').read_text());assert len(stages)==11
bank=json.loads((root/'results/10_programs/frozen/Programs_v1.json').read_text())
assert len(bank['meta_programs'])==5 and all(len(m['genes'])==50 and 'ENSG00000109736' not in m['genes'] for m in bank['meta_programs'])
stop=json.loads((root/'results/11_aggregation/audit/STOP_THIS_ROUTE.json').read_text());assert stop['stop_automatic_expression_program_association_optimization']
with (root/'results/11_aggregation/results/Aggregation_main_all10.tsv').open(encoding='utf-8') as f:agg=list(csv.DictReader(f,delimiter='\t'))
assert len(agg)==10 and all(float(r['FDR_main10'])>.05 for r in agg if r['cohort']=='discovery')
for p in (root/'scripts/project').glob('*.py'):ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
assert (root/'docs/RESEARCH_PROGRESS.md').is_file()
print(json.dumps(dict(status='pass',review_files_verified=len(rows),programs=5,aggregation_route_stop_preserved=True,large_data_loaded=False,scientific_analysis_rerun=False),indent=2))
