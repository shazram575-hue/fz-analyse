"""Execute user-confirmed 282 main jobs and nine separately labelled redocking controls."""
import csv,json,time,subprocess,hashlib,sys,threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[3];B=Path(__file__).resolve().parents[1];OLD=ROOT/'outputs/docking_audit_2026-10-24'
V=OLD/'software/vina_1.2.7_win.exe'
SEEDS=[20261024,20261025,20261026]
for d in ['docking','controls']:(B/d).mkdir(exist_ok=True)
def tab(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
PLAN=tab(B/'inputs/47_structure_and_search_plan.tsv')
large_lock=threading.Semaphore(1)
combine_lock=threading.Lock()
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def ready(g):
    p=B/'audit'/f'{g}_receptor_preparation.json'
    if not p.exists():return False
    try:r=json.loads(p.read_text());return r.get('exit_code')==0 and r.get('connectivity_builder') in ['standard_residue_templates_v1','standard_residue_templates_v2']
    except Exception:return False
def combined(g):
    with combine_lock:
        base=B/'prepared/receptors'/f'{g}.pdbqt';co=B/'prepared/cofactors'/f'{g}_retained.pdbqt'
        out=B/'prepared/receptors'/f'{g}_docking.pdbqt'
        text=base.read_text(encoding='utf8')
        if co.exists():text+='\n'+co.read_text(encoding='utf8')
        if not out.exists() or out.read_text()!=text:out.write_text(text,encoding='ascii')
        return out
def build_jobs():
    jobs=[]
    for p in PLAN:
        for state in ['pH74','neutral']:
            for seed in SEEDS:jobs.append({'gene':p['gene'],'state':state,'seed':seed,'kind':'main',
              'ligand':str((B/'prepared/ligands'/f'spiropachysine_{state}.pdbqt').relative_to(ROOT))})
    controls=json.loads((B/'inputs/Redocking_controls.json').read_text())
    for c in controls:
        for seed in SEEDS:jobs.append({'gene':c['gene'],'state':c['CCD'],'seed':seed,'kind':'control','ligand':c['ligand']})
    return jobs
def run_one(job):
    g=job['gene'];folder=B/('docking' if job['kind']=='main' else 'controls')
    tag=f'{g}_{job["state"]}_ex64_s{job["seed"]}';out=folder/f'{tag}.pdbqt';meta=folder/f'{tag}.json'
    if not ready(g):return None
    box=json.loads((B/'prepared'/f'{g}_box.json').read_text())
    receptor=combined(g);ligand=ROOT/job['ligand']
    hashes={'receptor':sha(receptor),'ligand':sha(ligand),'box':sha(B/'prepared'/f'{g}_box.json'),'vina':sha(V)}
    if meta.exists() and out.exists():
        old=json.loads(meta.read_text())
        if old.get('exit_code')==0 and old.get('input_SHA256')==hashes and 'REMARK VINA RESULT:' in out.read_text():return old
    cmd=[str(V),'--receptor',str(receptor),'--ligand',str(ligand),'--cpu','2','--seed',str(job['seed']),
         '--exhaustiveness','64','--num_modes','20','--energy_range','5','--out',str(out)]
    cmd+=sum(([f'--center_{a}',str(round(v,5))] for a,v in zip('xyz',box['center'])),[])
    cmd+=sum(([f'--size_{a}',str(round(v,5))] for a,v in zip('xyz',box['size'])),[])
    cmd=[str(Path(x).relative_to(ROOT)) if Path(x).is_absolute() and str(x).startswith(str(ROOT)) else x for x in cmd]
    big=box['volume_A3']>1e6
    if big:large_lock.acquire()
    try:
        start=time.monotonic()
        print('START',tag,job['kind'],flush=True)
        with (folder/f'{tag}.log').open('w',encoding='utf8') as f:proc=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
        scores=[float(l.split()[3]) for l in out.read_text().splitlines() if l.startswith('REMARK VINA RESULT:')] if out.exists() else []
        data={**job,'tag':tag,'cmd':cmd,'exit_code':proc.returncode,'elapsed_seconds':time.monotonic()-start,
              'input_SHA256':hashes,'scores_kcal_mol':scores,'best_vina_kcal_mol':scores[0] if scores else None,
              'poses_returned':len(scores),'site_rule':box['site_rule'],'volume_A3':box['volume_A3']}
        meta.write_text(json.dumps(data,indent=2),encoding='utf8')
        print('DONE' if proc.returncode==0 and scores else 'FAILED',tag,round(data['elapsed_seconds'],1),data['best_vina_kcal_mol'],flush=True)
        return data
    finally:
        if big:large_lock.release()
def main():
    if (B/'audit/ALL47_DOCKING_DELIVERY_COMPLETED.json').exists():
        print('All47 delivery already completed; no new docking runs started.',flush=True)
        return
    assert json.loads((B/'audit/Protocol_authorized.json').read_text())['main_runs']==282
    jobs=build_jobs();assert sum(j['kind']=='main' for j in jobs)==282
    with (B/'inputs/All47_main_and_control_jobs.tsv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(jobs[0]),delimiter='\t');w.writeheader();w.writerows(jobs)
    pending=jobs[:];finished=[];failed=[];last_progress=time.monotonic()
    # Run controls early, then mix target jobs so each target is represented before consuming repeats.
    pending.sort(key=lambda j:(j['kind']!='control',j['seed'],j['state']!='pH74',j['gene']))
    with ThreadPoolExecutor(max_workers=6) as ex:
        active={}
        while pending or active:
            for j in pending[:]:
                if len(active)>=6:break
                if ready(j['gene']):active[ex.submit(run_one,j)]=j;pending.remove(j)
            completed=[f for f in active if f.done()]
            for fut in completed:
                job=active.pop(fut)
                try:
                    data=fut.result()
                    if data and data['exit_code']==0 and data['poses_returned']>0:finished.append(data)
                    else:failed.append(data or job)
                except Exception as e:failed.append({**job,'error':str(e)})
                last_progress=time.monotonic()
            status={'main_completed':sum(d['kind']=='main' for d in finished),'main_expected':282,
                    'control_completed':sum(d['kind']=='control' for d in finished),'control_expected':9,
                    'failed':len(failed),'pending':len(pending),'running':list(active.values()),
                    'ready_receptors':sum(ready(p['gene']) for p in PLAN),'updated_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
            (B/'audit/Run_progress.json').write_text(json.dumps(status,indent=2),encoding='utf8')
            if pending and not active and time.monotonic()-last_progress>1800:
                failed.extend({'gene':j['gene'],'kind':j['kind'],'error':'receptor not ready after30min'} for j in pending);break
            if pending or active:time.sleep(2)
    (B/'results/All47_raw_run_summary.json').write_text(json.dumps(finished,indent=2),encoding='utf8')
    (B/'audit/Run_errors.json').write_text(json.dumps(failed,indent=2),encoding='utf8')
    if failed:raise SystemExit(1)
    assert sum(d['kind']=='main' for d in finished)==282
    (B/'audit/RAW_282_MAIN_RUNS_COMPLETED.json').write_text(json.dumps({'main_runs':282,'control_runs':9,'next':'pose and evidence evaluation pending'},indent=2),encoding='utf8')
if __name__=='__main__':main()
