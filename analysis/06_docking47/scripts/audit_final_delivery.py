import csv,json,hashlib,time,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];B=Path(__file__).resolve().parents[1]
def tab(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
checks=[]
def check(name,ok,details):checks.append({'check':name,'pass':bool(ok),'details':str(details)})
pool=tab(ROOT/'outputs/身份深核与独立交集_2026-10-24/results/disease_intersections/top10percent_Kd_intersection_IC50_disease1881.tsv')
genes={r['gene_symbol'] for r in pool};summary=tab(B/'results/Target_state_summary.tsv');evidence=tab(B/'results/47_integrated_docking_disease_evidence.tsv')
main=[json.loads(p.read_text(encoding='utf8')) for p in (B/'docking').glob('*.json')]
controls=[json.loads(p.read_text(encoding='utf8')) for p in (B/'controls').glob('*.json')]
check('original47_target_set',len(genes)==47 and {r['gene'] for r in evidence}==genes,len(evidence))
keys={(r['gene'],r['state'],r['seed']) for r in main if r.get('exit_code')==0 and r.get('poses_returned',0)>0}
expected={(g,s,k) for g in genes for s in ['pH74','neutral'] for k in [20261024,20261025,20261026]}
check('exact282_unique_main_jobs',keys==expected and len(main)==282,f'{len(keys)}/282')
valid=True;bad=[]
for r in main:
    p=B/'docking'/f'{r["tag"]}.pdbqt'
    scores=[float(l.split()[3]) for l in p.read_text(encoding='ascii').splitlines() if l.startswith('REMARK VINA RESULT:')]
    if scores!=r['scores_kcal_mol']:valid=False;bad.append(r['tag'])
    receptor=B/'prepared/receptors'/f'{r["gene"]}_docking.pdbqt'
    ligand=ROOT/r['ligand'];box=B/'prepared'/f'{r["gene"]}_box.json'
    if any(sha(p)!=r['input_SHA256'][k] for k,p in [('receptor',receptor),('ligand',ligand),('box',box)]):valid=False;bad.append(r['tag']+'_input')
check('raw_scores_and_input_hashes',valid,bad)
check('9_separate_control_jobs',len(controls)==9 and all(r.get('exit_code')==0 for r in controls),len(controls))
check('94_complete_state_summaries',len(summary)==94 and all(int(r['completed_seeds'])==3 for r in summary),len(summary))
poses=tab(B/'results/All_completed_top_pose_checks.tsv')
check('282_pose_checks_present',len(poses)==282,len(poses))
check('ligand_identity_and_9_center_stereo_retained',all(r['stereochemistry_and_identity_match']=='True' for r in poses),sum(r['stereochemistry_and_identity_match']=='True' for r in poses))
check('all_raw_pose_SDFs_exported',all((ROOT/r['sdf_path']).exists() for r in poses),len(poses))
prep=[]
for g in genes:
    p=B/'audit'/f'{g}_receptor_preparation.json'
    prep.append(p.exists() and json.loads(p.read_text(encoding='utf8')).get('exit_code')==0)
check('47_receptors_prepared',all(prep),sum(prep))
regions=tab(B/'inputs/47_actual_search_regions.tsv')
check('47_actual_search_regions',{r['gene'] for r in regions}==genes,len(regions))
check('5_authorized_global_searches',{r['gene'] for r in regions if r['site_rule']=='global_search'}=={'ATP5MC1','FXYD3','MT-ND4L','NHSL2','SMIM22'},'explicit global group')
before=tab(B/'audit/Protected_inputs_before.tsv');changed=[]
for r in before:
    p=ROOT/r['path']
    if not p.exists() or sha(p)!=r['sha256']:changed.append(r['path'])
check('protected_upstream_public_inputs_unchanged',not changed,changed)
deleted=tab(B/'audit/Authorized_old_docking_deletion_manifest.tsv');remaining=[r['relative_path'] for r in deleted if (ROOT/r['relative_path']).exists()]
check('old_docking_outputs_absent',not remaining,remaining)
check('final_report_exists',(B/'结果报告_47靶点反向对接.md').exists(),'47 row report')
check('figures_PNG_SVG_present',all((B/'figures'/f'{stem}.{ext}').exists() for stem in ['01_all47_scores_and_consistency','02_redocking_controls'] for ext in ['png','svg']),'two figure sets')
with (B/'audit/Final_independent_validation.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(checks[0]),delimiter='\t');w.writeheader();w.writerows(checks)
result={'technical_checks_passed':sum(c['pass'] for c in checks),'technical_checks_total':len(checks),
        'technical_delivery_complete':all(c['pass'] for c in checks),'main_runs':282,'control_runs':9,
        'geometry_pass_top_poses':sum(r['PoseBusters_geometry_pass']=='True' for r in poses),
        'no_major_receptor_clash_top_poses':sum(r['no_major_receptor_clash']=='True' for r in poses),
        'scientific_pose_and_control_findings':'reported without replacing failed results',
        'protected_input_files':len(before),'deleted_old_output_files':len(deleted),'timestamp':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
(B/'audit/Final_validation_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(result,ensure_ascii=False,indent=2))
if not result['technical_delivery_complete']:raise SystemExit(1)
(B/'audit/ALL47_TECHNICAL_CHECKS_COMPLETED.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
subprocess.run([sys.executable,'-X','utf8',str(B/'scripts/update_runtime_status.py')],cwd=ROOT,check=True)
