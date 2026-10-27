"""Local supervisor for this authorized computation; not a recurring automation."""
import json,subprocess,sys,time,os
from pathlib import Path
os.environ['PYTHONUTF8']='1'
B=Path(__file__).resolve().parents[1];ROOT=B.parents[1]
if (B/'audit/ALL47_DOCKING_DELIVERY_COMPLETED.json').exists():
    print('All47 delivery already complete. Read README and final report; no queue restart.')
    raise SystemExit(0)
def run(script):
    logfile=B/'logs'/f'{script}_supervised.log'
    with logfile.open('a',encoding='utf8') as f:
        p=subprocess.run([sys.executable,'-X','utf8','-u',str(B/'scripts'/f'{script}.py')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError(f'{script} exit={p.returncode}; see {logfile}')
last_status=None;last_checked=-12;last_controls=-1
while True:
    run('reuse_identical_context')
    try:p=json.loads((B/'audit/Run_progress.json').read_text(encoding='utf8'))
    except (FileNotFoundError,json.JSONDecodeError):time.sleep(5);continue
    signature=(p['main_completed'],p['control_completed'],p['failed'],p['ready_receptors'])
    if signature!=last_status:
        run('update_runtime_status');print(json.dumps({'main':p['main_completed'],'controls':p['control_completed'],'ready_receptors':p['ready_receptors'],'failed':p['failed']},ensure_ascii=False),flush=True)
        last_status=signature
    if p['failed']:
        raise RuntimeError('Docking failures reported; independent remaining jobs continue. Inspect Run_progress and docking logs.')
    if p['main_completed']>=last_checked+12 or p['control_completed']>last_controls or p['main_completed']==282:
        run('evaluate_completed');last_checked=p['main_completed'];last_controls=p['control_completed']
    if p['main_completed']==282 and p['control_completed']==9 and not p['running']:
        run('evaluate_completed');run('assemble_final_report');run('audit_final_delivery');run('update_runtime_status')
        print('ALL47_RAW_DATA_AND_REPORT_READY_FOR_VISUAL_REVIEW',flush=True);break
    time.sleep(30)
