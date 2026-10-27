"""Inspect local inputs; enumerate explicitly authorized obsolete docking outputs."""
import csv, hashlib, json, statistics
from pathlib import Path
import gemmi
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
B = Path(__file__).resolve().parents[1]
OLD = ROOT / 'outputs/docking_audit_2026-10-24'
if (B/'audit/Old_docking_deletion_completed.json').exists():
    raise SystemExit('Initial intake and authorized deletion already completed. Preserve the original135-file manifest; resume from README and Run_progress.json.')
for name in ['audit', 'inputs', 'results', 'sources', 'structures', 'prepared', 'logs', 'scripts']:
    (B / name).mkdir(parents=True, exist_ok=True)

def readtab(p):
    return list(csv.DictReader(p.open(encoding='utf-8-sig'), delimiter='\t'))

def write(p, rows):
    with p.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)

def sha(p):
    return hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()

poolpath = ROOT/'outputs/身份深核与独立交集_2026-10-24/results/disease_intersections/top10percent_Kd_intersection_IC50_disease1881.tsv'
pool = readtab(poolpath)
inventory = readtab(OLD/'results/47_structure_pocket_assessment.tsv')
background = {r['uniprot_accession']: r for r in readtab(ROOT/'outputs/身份深核与独立交集_2026-10-24/background/All_reviewed_human_protein_background.tsv')}
assert len(pool) == len(inventory) == 47
assert {r['gene_symbol'] for r in pool} == {r['gene'] for r in inventory}
genes = {r['gene_symbol']: r for r in pool}
rows = []
for r in inventory:
    g, acc = r['gene'], r['accession']
    af = OLD/'structures'/f'{acc}_AF.pdb'
    st = gemmi.read_structure(str(af))
    chain = st[0][0]
    residues = [x for x in chain if x.het_flag == 'A']
    seq = ''.join(gemmi.find_tabulated_residue(x.name).one_letter_code for x in residues)
    assert seq == background[acc]['protein_sequence'], (g, acc)
    xyz = [(a.pos.x,a.pos.y,a.pos.z) for x in residues for a in x if not a.element.is_hydrogen]
    span = [max(x[i] for x in xyz)-min(x[i] for x in xyz) for i in range(3)]
    folder = 'AF_MFSD10' if g == 'MFSD10' else 'AF47'
    pp = OLD/'pockets'/folder/f'{acc}_AF.pdb_predictions.csv'
    assert pp.exists(), pp
    pockets = [{k.strip(): v.strip() for k,v in x.items()} for x in csv.DictReader(pp.open())]
    assert len(pockets) == int(r['P2Rank_pocket_count'])
    rr = dict(r)
    rr.update({
        'AF_path':str(af.relative_to(ROOT)), 'AF_SHA256':sha(af),
        'AF_verified_residues':len(residues), 'AF_extent_x_A':round(span[0],2),
        'AF_extent_y_A':round(span[1],2), 'AF_extent_z_A':round(span[2],2),
        'pocket_predictions_path':str(pp.relative_to(ROOT)), 'pocket_predictions_SHA256':sha(pp),
        'bulk_log2FC':genes[g]['bulk_log2FC'], 'SC_up_cell_types':genes[g]['SC_up_cell_types'],
        'SC_down_cell_types':genes[g]['SC_down_cell_types'], 'PB_support_cell_types':genes[g]['PB_DE_support_cell_types'],
        'Kd_proteome_rank':genes[g]['Kd_best_protein_rank'], 'IC50_proteome_rank':genes[g]['IC50_best_protein_rank'],
        'new_docking_status':'awaiting_structure_and_region_plan',
    })
    rows.append(rr)
write(B/'inputs/47_verified_input_inventory.tsv', rows)

ligands=[]
for state in ['neutral','pH74']:
    p = OLD/'prepared'/f'spiropachysine_{state}.sdf'
    m = next(Chem.ForwardSDMolSupplier(p.open('rb'), removeHs=False))
    assert m and m.GetNumConformers()==1 and m.GetConformer().Is3D()
    centers=Chem.FindMolChiralCenters(Chem.RemoveHs(m), includeUnassigned=True)
    assert len(centers)==9 and all(x[1]!='?' for x in centers)
    ligands.append({'state':state, 'formal_charge':Chem.GetFormalCharge(m), 'specified_centers':len(centers),
                    'InChIKey':Chem.MolToInchiKey(m), 'sdf_path':str(p.relative_to(ROOT)), 'SHA256':sha(p)})
write(B/'inputs/Ligand_input_audit.tsv', ligands)
assert ligands[0]['InChIKey']=='SLGWGPQWJRVPAD-PZOAWPBASA-N'

protected = [poolpath]
for folder in ['outputs/第一阶段_2026-10-23/results', 'outputs/单细胞基础补全_2026-10-23/results',
               'outputs/DeepPurpose_预测_2026-10-24/main_input','outputs/身份深核与独立交集_2026-10-24/results',
               'outputs/docking_audit_2026-10-24/pubmed_comprehensive_recheck',
               'outputs/docking_audit_2026-10-24/pubmed_zero_recheck']:
    protected.extend(p for p in (ROOT/folder).rglob('*') if p.is_file())
protected.extend(OLD/'structures'/f'{r["accession"]}_AF.pdb' for r in inventory)
write(B/'audit/Protected_inputs_before.tsv', [{'path':str(p.relative_to(ROOT)), 'bytes':p.stat().st_size, 'sha256':sha(p)} for p in sorted(set(protected))])

obsolete={}
for folder in ['docking','figures','audit/报告修订前快照','audit/读图修正前快照','audit/本次接手文档同步']:
    for p in (OLD/folder).rglob('*'):
        if p.is_file(): obsolete[p]='旧对接输出或旧结论文档快照，按用户确认清除'
for pat in ['Pair_*','Pilot_*','TBZ_*','pilot_scores.json','validation_scores.json']:
    for p in (OLD/'results').glob(pat):
        if p.is_file(): obsolete[p]='旧两候选姿势、评分、接触或质量检查结果'
for name in ['两候选配对对接_结果报告.md','结果图逐张读法.md','试验规则_计算前固定.md']:
    obsolete[OLD/name]='旧对接报告或旧限定试算规则'
for name in ['PAIR_PILOT_DELIVERY_COMPLETED.txt','PILOT_RUNS_FINISHED.txt','VALIDATION_RUNS_FINISHED.txt',
             'Deliverable_SHA256_manifest.tsv','Figure_visual_QA.md','Independent_integrity_validation.tsv','Integrity_summary.json',
             '报告解释修订_原始计算保留.json','QC_limited_to_geometry_and_atom_distances.json',
             'QC_streamlined_to_official_dock_fast.json','serial_to_parallel_interrupted_processes.json']:
    obsolete[OLD/'audit'/name]='与已删除旧对接结果配套的完成标记、检查记录或旧交付清单'
for pat in ['pilot*','*pair*','figures*','TBZ*','validation_run.log']:
    for p in (OLD/'logs').glob(pat):
        if p.is_file(): obsolete[p]='旧对接、旧QA或旧绘图运行日志'
for name in ['run_docking.py','evaluate_docking.py','quality_check_pilot.py','make_figures.py',
             'summarize_pair_contacts.py','rewrite_decision_reports.py','finalize_pair_delivery.py',
             'close_delivery.py','sync_current_handoff.py']:
    obsolete[OLD/'scripts'/name]='只服务于旧两候选结果的运行或报告脚本，清除以免重建旧结论'
manifest=[]
for p,reason in sorted(obsolete.items()):
    if p.is_file():
        assert p.resolve().is_relative_to(OLD.resolve())
        manifest.append({'absolute_path':str(p.resolve()),'relative_path':str(p.relative_to(ROOT)),
                         'bytes':p.stat().st_size,'sha256':sha(p),'reason':reason})
write(B/'audit/Authorized_old_docking_deletion_manifest.tsv',manifest)
summary={'targets':len(rows),'AF_sequence_matches':47,'with_PDB':sum(int(r['experimental_PDB_count'])>0 for r in rows),
         'with_AF_pockets':sum(int(r['P2Rank_pocket_count'])>0 for r in rows),
         'without_AF_pockets':[r['gene'] for r in rows if int(r['P2Rank_pocket_count'])==0],
         'deletion_manifest_files':len(manifest),'deletion_manifest_bytes':sum(r['bytes'] for r in manifest),
         'new_docking_runs_completed':0,'scope_authorized':'all47_experimental_first_AF_fallback_include_global_search',
         'deletion_authorized':'old_docking_outputs_and_corresponding_conclusions_only'}
(B/'audit/Intake_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
