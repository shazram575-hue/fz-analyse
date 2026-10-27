"""Index historical error/limit messages without copying log content or running code."""
from pathlib import Path
import csv,re,json,collections
root=Path(__file__).resolve().parents[2];out=root/'docs/audit/2026-10-26'
with (out/'before/FILE_INVENTORY.tsv').open(encoding='utf-8') as f:files=list(csv.DictReader(f,delimiter='\t'))
patterns={'traceback':r'^Traceback','R_execution_halted':r'Execution halted','R_error':r'^Error(?: in|:)', 'permission_error':r'(?:PermissionError|Access is denied|拒绝访问)', 'missing_module':r'ModuleNotFoundError','iteration_limit':r'ConvergenceWarning|Maximum number of iterations','fatal_error':r'fatal error'}
rows=[];skipped=[]
for r in files:
 if r['role']!='log' or Path(r['path']).suffix.lower() not in ['.log','.txt']:continue
 if int(r['bytes'])>5*1024**2:skipped.append(r['path']);continue
 p=root/r['path']
 for line,text in enumerate(p.read_text(encoding='utf-8',errors='replace').splitlines(),1):
  for kind,pat in patterns.items():
   if re.search(pat,text):rows.append(dict(path=r['path'],line=line,diagnostic=kind,interpretation='Historical occurrence only; final stage status requires completion evidence, not this message count'))
with (out/'LOG_DIAGNOSTIC_INDEX.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['path','line','diagnostic','interpretation'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
(out/'LOG_DIAGNOSTIC_SUMMARY.json').write_text(json.dumps(dict(indexed_occurrences=len(rows),files_with_diagnostics=len(set(r['path'] for r in rows)),by_kind=dict(collections.Counter(r['diagnostic'] for r in rows)),large_logs_not_scanned=skipped,warning_counts_not_complete_convergence_counts=True),ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(files_with_diagnostics=len(set(r['path'] for r in rows)),occurrences=len(rows),large_logs_skipped=len(skipped)),indent=2))
