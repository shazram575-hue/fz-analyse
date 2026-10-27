import csv,json,statistics,re,hashlib,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import gemmi
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import rdMolAlign
ROOT=Path(__file__).resolve().parents[3];B=Path(__file__).resolve().parents[1];OLD=ROOT/'outputs/docking_audit_2026-10-24'
def table(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def truth(x):return str(x)=='True'
def num(x):return float(x) if x not in ['',None,'None'] else None
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def inspect_returned_target_poses(g,run_rows):
    """Diagnostic only: preserve primary rank1 scores, inspect existing returned poses for target attribution."""
    ext=json.loads((B/'audit'/f'{g}_extraction.json').read_text(encoding='utf8'))
    target_chains=set(ext['target_chains'])
    protein=gemmi.read_structure(str(B/'prepared/receptors'/f'{g}_prepared.pdb'))
    xyz=np.array([[a.pos.x,a.pos.y,a.pos.z] for c in protein[0] if c.name in target_chains for r in c for a in r if not a.element.is_hydrogen])
    tree=cKDTree(xyz)
    partner_xyz=np.array([[a.pos.x,a.pos.y,a.pos.z] for c in protein[0] if c.name not in target_chains for r in c for a in r if not a.element.is_hydrogen])
    partner_tree=cKDTree(partner_xyz) if len(partner_xyz) else None
    records=[];chosen=[]
    for row in sorted(run_rows,key=lambda x:int(x['seed'])):
        sdf=ROOT/row['sdf_path'];poses=[];molecules={}
        for rank,m in enumerate(Chem.ForwardSDMolSupplier(sdf.open('rb'),removeHs=True),1):
            if m is None:continue
            d=tree.query(m.GetConformer().GetPositions())[0]
            poses.append({'rank':rank,'target_contacts_within4A':int(sum(d<=4)),
                'partner_contacts_within4A':int(sum(partner_tree.query(m.GetConformer().GetPositions())[0]<=4)) if partner_tree else 0,
                'minimum_target_distance_A':float(min(d)),'score_kcal_mol':float(m.GetProp('vina_kcal_mol'))})
            molecules[rank]=m
        eligible=[x for x in poses if x['target_contacts_within4A']>=5]
        best=min(eligible,key=lambda x:x['score_kcal_mol']) if eligible else None
        if best:chosen.append(molecules[best['rank']])
        records.append({'gene':g,'state':'pH74','seed':int(row['seed']),'source_SDF_SHA256':sha(sdf),
            'primary_target_distance_A':poses[0]['minimum_target_distance_A'],
            'returned_poses_with_at_least1_target_contact':sum(x['target_contacts_within4A']>0 for x in poses),
            'returned_poses_with_at_least5_target_contacts':len(eligible),'best_target_associated_returned_pose':best,'all_poses':poses})
    rmsds=[float(rdMolAlign.CalcRMS(a,b,maxMatches=100000)) for i,a in enumerate(chosen) for b in chosen[i+1:]]
    result={'gene':g,'purpose':'returned-pose target-attribution diagnostic; primary scores and ranks unchanged',
        'runs':records,'best_target_associated_pose_pairwise_RMSD_A':rmsds,
        'all3_have_target_associated_returned_pose':len(chosen)==3}
    (B/'audit'/f'{g}_returned_pose_target_attribution_diagnostic.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    return result
states=table(B/'results/Target_state_summary.tsv');checks=table(B/'results/All_completed_top_pose_checks.tsv')
controls=table(B/'results/Redocking_control_results.tsv')
assert len(states)==94 and all(int(r['completed_seeds'])==3 for r in states)
assert len(checks)==282 and len(controls)==9
plan={r['gene']:r for r in table(B/'inputs/47_structure_and_search_plan.tsv')}
inv={r['gene']:r for r in table(B/'inputs/47_verified_input_inventory.tsv')}
original47={r['gene_symbol']:r for r in table(ROOT/'outputs/身份深核与独立交集_2026-10-24/results/disease_intersections/top10percent_Kd_intersection_IC50_disease1881.tsv')}
literature={r['gene']:r for r in table(OLD/'pubmed_comprehensive_recheck/47_comprehensive_PubMed_counts.tsv')}
disease=table(ROOT/'outputs/DeepPurpose_预测_2026-10-24/main_input/Disease_candidates_gene_celltype_evidence.tsv')
disease=[r for r in disease if r['gene_symbol'] in plan and r['diagnostic_annotation']=='False']
write(B/'results/47_target_celltype_disease_evidence.tsv',disease)
state={(r['gene'],r['state']):r for r in states}
evidence=[]
for g in sorted(plan):
    p=plan[g];r=inv[g];main=state[g,'pH74'];neutral=state[g,'neutral']
    good=all(truth(main[k]) for k in ['all_top_poses_geometry_pass','all_top_poses_identity_pass','all_top_poses_no_major_clash'])
    main_checks=[x for x in checks if x['gene']==g and x['state']=='pH74']
    target_counts=[int(x['ligand_heavy_atoms_within4A_target']) for x in main_checks]
    partner_counts=[int(x['ligand_heavy_atoms_within4A_partner']) for x in main_checks]
    repeat=num(main['pairwise_pose_RMSD_max_A'])<=2
    contacts=int(main['min_ligand_heavy_atoms_within4A_target'])
    conf=num(main['contact_residue_AF_pLDDT_mean']);global_search=p['site_rule']=='global_search'
    neutral_good=all(truth(neutral[k]) for k in ['all_top_poses_geometry_pass','all_top_poses_identity_pass','all_top_poses_no_major_clash'])
    neutral_conf=num(neutral['contact_residue_AF_pLDDT_mean'])
    neutral_stable=num(neutral['pairwise_pose_RMSD_max_A'])<=2
    neutral_contacts=int(neutral['min_ligand_heavy_atoms_within4A_target'])
    actual_box=json.loads((B/'prepared'/f'{g}_box.json').read_text(encoding='utf8'))
    pocket_probability=num(actual_box['basis'].get('pocket_probability'))
    known_site=p['site_rule'].startswith('experimental_ligand:')
    reasons=[]
    if not good:reasons.append('存在姿势几何或受体碰撞检查未通过项')
    if not repeat:reasons.append('最高排名姿势在种子间位置不同')
    if contacts<5:reasons.append('药物与目标链的近距离接触较少')
    if conf is not None and conf<70:reasons.append('实际接触区域模型可信度偏低')
    if global_search:reasons.append('全局探索，需要进一步定位可重复的局部位点')
    if good and contacts>=5 and (conf is None or conf>=70) and not global_search:
        support='有可重复的结构支持'
        action='可进入小规模结合与药效依赖性验证'
        if not known_site:
            support='有可重复姿势，位点为预测';action='可进入结合预实验，并核查预测位点'
        if not repeat:
            support='有结构可行性，最高姿势尚不一致';action='可做结合预实验，具体结合方式继续核验'
    elif global_search and good:
        support='获得全局探索线索';action='先核查局部位点和模型，再决定实验优先级'
    else:
        support='当前结构支持不足';action='先解决表内具体问题，再投入位点机制实验'
    if max(target_counts)==0:
        support='最高排名仅伙伴链有4Å接触';action='主评分不归属于单靶点，另查目标相关返回姿势'
        diagnostic=inspect_returned_target_poses(g,main_checks)
        dd=[x['primary_target_distance_A'] for x in diagnostic['runs']]
        reasons.append(f'主状态3/3最高排名姿势均无4Å内目标链接触，最近目标距离{min(dd):.2f}–{max(dd):.2f}Å')
        if diagnostic['all3_have_target_associated_returned_pose']:
            support='最高排名归于伙伴，另有目标相关候选';action='保留返回姿势线索，先明确具体结合亚基'
            ranks=[x['best_target_associated_returned_pose']['rank'] for x in diagnostic['runs']]
            reasons.append('三次返回中均另有至少5个药物重原子接近目标的姿势，排名分别'+('/'.join(map(str,ranks)))+'；独立附表保留，不替换主评分')
    elif min(target_counts)==0:
        support='部分最高排名姿势未接触目标链';action='先核对姿势归属与目标相关位点'
    elif any(partner_counts) and good and not global_search:
        support='有可重复的目标-伙伴界面姿势' if repeat else '存在多种目标-伙伴界面姿势'
        action='可检验复合物相关结合，并区分具体结合亚基'
    if any(partner_counts):reasons.append(f'主状态近距离接触的药物重原子数：目标{min(target_counts)}–{max(target_counts)}，伙伴{min(partner_counts)}–{max(partner_counts)}（两者可重叠）')
    if not global_search and (not good or (conf is not None and conf<70)) and neutral_good and neutral_stable and neutral_contacts>=5 and (neutral_conf is None or neutral_conf>=70):
        support='结果对药物状态敏感，中性姿势可用';action='保留中性状态结合假说，优先检验实际条件下的结合'
        reasons.append('预定+1主状态有距离或模型质量问题，中性状态三种子均通过相应检查；两状态分开保留')
    if not reasons:reasons.append('主状态三种子姿势收敛、几何可用且接触目标链')
    if pocket_probability is not None:reasons.append(f'实际采用的预测口袋P2Rank概率为{pocket_probability:.3f}')
    if original47[g]['has_only_donor_limited_SC_support']=='True':
        reasons.append('疾病侧目前仅有供者不足类型的细胞级支持，实验需匹配该细胞背景')
    if g=='MT-ND4':reasons.append('对接使用459aa P03905；原预测中的99aa C0HME5条目另保留')
    mutated_contacts=sorted({x['sequence_difference_residues_within5A'] for x in checks if x['gene']==g and x['state']=='pH74' and x.get('sequence_difference_residues_within5A')})
    if mutated_contacts:reasons.append('姿势5Å内有实验构建序列差异：'+' / '.join(mutated_contacts))
    ctr=[x for x in controls if x['gene']==g]
    if ctr and not any(truth(x['top_pose_within2A']) for x in ctr):
        support='有计算姿势，参照姿势验证未通过';action='可结合功能线索讨论预实验，精确位点结论暂不采用';reasons.append('本轮参照最高排名姿势0/3达到2Å')
    novelty='严格0' if int(literature[g]['expanded_domain_count'])==0 else ('机构缩写误命中，靶点相关0' if g=='SAMD8' else 'PubMed有命中')
    evidence.append({'gene':g,'accession':p['accession'],'structure':p['structure'],'source':p['source'],
        'site_rule':p['site_rule'],'pH74_median_kcal_mol':main['score_median_kcal_mol'],
        'actual_P2Rank_pocket_probability':pocket_probability,
        'pH74_min_kcal_mol':main['score_min_kcal_mol'],'pH74_max_kcal_mol':main['score_max_kcal_mol'],
        'neutral_median_kcal_mol':neutral['score_median_kcal_mol'],'pH74_max_seed_RMSD_A':main['pairwise_pose_RMSD_max_A'],
        'neutral_max_seed_RMSD_A':neutral['pairwise_pose_RMSD_max_A'],
        'neutral_geometry_identity_clash_all_pass':neutral_good,'neutral_min_target_contact_heavy_atoms':neutral_contacts,
        'neutral_contact_AF_pLDDT_mean':neutral_conf,
        'geometry_identity_clash_all_pass':good,'min_target_contact_heavy_atoms':contacts,'contact_AF_pLDDT_mean':conf,
        'max_target_contact_heavy_atoms':max(target_counts),'min_partner_contact_heavy_atoms':min(partner_counts),
        'max_partner_contact_heavy_atoms':max(partner_counts),
        'seeds_no_major_receptor_clash_pass':sum(truth(x['no_major_receptor_clash']) for x in main_checks),
        'sequence_difference_contacts_within5A':' / '.join(mutated_contacts),
        'bulk_log2FC':r['bulk_log2FC'],'SC_up_cell_types':r['SC_up_cell_types'],'SC_down_cell_types':r['SC_down_cell_types'],
        'PB_support_cell_types':r['PB_support_cell_types'],'Kd_proteome_rank':r['Kd_proteome_rank'],'IC50_proteome_rank':r['IC50_proteome_rank'],
        'has_only_donor_limited_SC_support':original47[g]['has_only_donor_limited_SC_support'],
        'donor_limited_SC_cell_types':original47[g]['donor_limited_SC_cell_types'],
        'Kd_prediction_proteins':original47[g]['Kd_reviewed_proteins'],
        'IC50_prediction_proteins':original47[g]['IC50_reviewed_proteins'],
        'Kd_CNN_sequence_truncation_flag':original47[g]['Kd_any_CNN_truncation'],
        'IC50_CNN_sequence_truncation_flag':original47[g]['IC50_any_CNN_truncation'],
        'PubMed_expanded_count':literature[g]['expanded_domain_count'],'novelty_annotation':novelty,
        'structure_support':support,'suggested_next_step':action,'specific_reason':'；'.join(reasons)})
write(B/'results/47_integrated_docking_disease_evidence.tsv',evidence)

for name in ['Microsoft YaHei','SimHei','Noto Sans CJK SC']:
    if any(f.name==name for f in font_manager.fontManager.ttflist):plt.rcParams['font.family']=name;break
plt.rcParams['axes.unicode_minus']=False
figdir=B/'figures';figdir.mkdir(exist_ok=True)
fig,axes=plt.subplots(1,2,figsize=(13,18),gridspec_kw={'width_ratios':[2.2,1]})
y=np.arange(47)
for i,r in enumerate(evidence):
    marker='s' if r['site_rule']=='global_search' else 'o'
    axes[0].plot([num(r['pH74_min_kcal_mol']),num(r['pH74_max_kcal_mol'])],[i,i],color='#246C9A',lw=2)
    axes[0].scatter(num(r['pH74_median_kcal_mol']),i,color='#246C9A',marker=marker,s=25)
    axes[0].scatter(num(r['neutral_median_kcal_mol']),i,color='#D39C31',marker=marker,s=22)
    axes[1].scatter(num(r['pH74_max_seed_RMSD_A']),i,color='#246C9A',marker=marker,s=25)
axes[0].set_yticks(y,[r['gene']+'  '+r['structure'] for r in evidence]);axes[1].set_yticks(y,[])
for ax in axes:ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
axes[0].set_title('+1（蓝）/中性（金）中位数；线为+1三种子范围')
axes[0].set_xlabel('Vina评分 (kcal/mol)');axes[1].set_title('+1最高排名姿势：种子间最大RMSD')
axes[1].set_xscale('symlog',linthresh=1)
right_limit=max(2.5,max(num(r['pH74_max_seed_RMSD_A']) for r in evidence)*1.25)
axes[1].set_xlim(0,right_limit)
right_ticks=[x for x in [0,1,2,5,10,20,40,80] if x<=right_limit]
axes[1].set_xticks(right_ticks,[str(x) for x in right_ticks])
axes[1].set_xlabel('RMSD (Å；1以下线性，以上压缩跨度)');axes[1].axvline(2,color='gray',ls='--',lw=1)
fig.suptitle('47靶点完整结果；方形为全局探索；按基因排列',y=.995)
fig.tight_layout();fig.savefig(figdir/'01_all47_scores_and_consistency.png',dpi=160,bbox_inches='tight');fig.savefig(figdir/'01_all47_scores_and_consistency.svg',bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4))
for j,g in enumerate(['PPARG','PNP','SLC18B1']):
    rr=[r for r in controls if r['gene']==g]
    ax.scatter([j-.1,j,j+.1],[num(r['top_pose_RMSD_A']) for r in rr],color='#A74343',label='最高排名姿势' if j==0 else None)
    ax.scatter([j-.1,j,j+.1],[num(r['best_returned_pose_RMSD_A']) for r in rr],color='#737373',marker='s',label='返回姿势中最接近实验者' if j==0 else None)
ax.axhline(2,ls='--',color='#246C9A');ax.set_xticks(range(3),['PPARG / BRL','PNP / JTC','SLC18B1 / YHL'])
ax.set_ylabel('对实验姿势的重原子RMSD (Å)');ax.legend();ax.set_title('本轮重新计算的9次参照复对接')
fig.tight_layout();fig.savefig(figdir/'02_redocking_controls.png',dpi=180);fig.savefig(figdir/'02_redocking_controls.svg');plt.close(fig)

lines=['# 47个共有靶点反向对接结果','',
f'完成日期：{time.strftime("%Y-%m-%d %H:%M:%S %z")}。按用户确认的方案完成47靶点×两状态×三种子＝282次主对接，另完成9次参照复对接。旧两候选对接结果已清除，本报告仅使用新一轮计算。','',
'## 这轮结果支持怎样继续','',
'47个对象均已有实际对接结果，下表逐项列出结构支持和建议的下一步。主状态为pH7.4规则生成的+1，中性状态单列；每个数值为三种子中位数。主攻靶点仍需结合这47项证据与你讨论确定。','',
'判断依据是姿势质量、三种子位置一致性、实际目标链接触、位点来源和局部结构可信度；表内不设统一Vina分数阳性线。全局探索单列解释。','',
'| 靶点 | +1中位数 | 中性中位数 | +1种子最大RMSD Å | 结构支持 | 后续建议 |',
'|---|---:|---:|---:|---|---|']
for r in evidence:lines.append(f'| {r["gene"]} | {num(r["pH74_median_kcal_mol"]):.3f} | {num(r["neutral_median_kcal_mol"]):.3f} | {num(r["pH74_max_seed_RMSD_A"]):.2f} | {r["structure_support"]} | {r["suggested_next_step"]} |')
lines+=['','282个最高排名姿势的身份/立体信息及基本分子几何检查均通过；268个通过本轮受体距离检查。ADRA2A与NR3C1在两状态均有距离问题，NR3C1当前口袋中的评分还为正值；S1PR1的+1姿势存在距离问题，而中性状态三次均通过，因此按状态敏感处理。GPRC5A、IFI27的主要问题是实际接触区域模型可信度偏低。各项原始指标均保留在检查表中。']
lines+=['','## 复合物中的姿势归属','',
'下表按+1状态三次最高排名姿势，分别统计4Å内接近目标链和伙伴链的药物重原子数。两类计数可以重叠；它们描述接触位置，不是能量分解。','',
'| 靶点 | 目标链接触数范围 | 伙伴链接触数范围 | 当前解释 |','|---|---:|---:|---|']
for r in evidence:
    if r['max_partner_contact_heavy_atoms']:
        lines.append(f'| {r["gene"]} | {r["min_target_contact_heavy_atoms"]}–{r["max_target_contact_heavy_atoms"]} | {r["min_partner_contact_heavy_atoms"]}–{r["max_partner_contact_heavy_atoms"]} | {r["structure_support"]} |')
for r in evidence:
    diagnostic_path=B/'audit'/f'{r["gene"]}_returned_pose_target_attribution_diagnostic.json'
    if not diagnostic_path.exists():continue
    diagnostic=json.loads(diagnostic_path.read_text(encoding='utf8'))
    distances=[x['primary_target_distance_A'] for x in diagnostic['runs']]
    lines+=['',f'**{r["gene"]}的具体情况：**最高排名姿势距目标链最近重原子仍为{min(distances):.2f}–{max(distances):.2f}Å，原主评分对应伙伴附近的位置。']
    if diagnostic['all3_have_target_associated_returned_pose']:
        candidates=[x['best_target_associated_returned_pose'] for x in diagnostic['runs']]
        ranks='/'.join(str(x['rank']) for x in candidates)
        candidate_scores='/'.join(f'{x["score_kcal_mol"]:.3f}' for x in candidates)
        rmsvals=diagnostic['best_target_associated_pose_pairwise_RMSD_A']
        lines.append(f'进一步核对每次已返回的20个姿势，三次均能找到至少5个药物重原子在4Å内接近目标的候选，分别排第{ranks}位，评分分别为{candidate_scores} kcal/mol；这三个候选姿势相互RMSD为{min(rmsvals):.2f}–{max(rmsvals):.2f}Å。它们作为目标相关返回姿势另列，保持原主评分与排名。')
    lines.append(f'[返回姿势逐项归属记录](audit/{diagnostic_path.name})。下一步判断应区分复合物结合与具体亚基归属。')
lines+=['','[完整47行证据与具体理由](results/47_integrated_docking_disease_evidence.tsv)包含Bulk方向、单细胞来源、患者PB、独立预测排名、PubMed结果及逐靶点理由。[逐细胞类型完整疾病证据](results/47_target_celltype_disease_evidence.tsv)保留统计量和PB方向。','',
'## 全47结果图','',
'![完整47靶点](figures/01_all47_scores_and_consistency.png)','',
'蓝点为+1状态、金点为中性状态，横线为+1三次分数范围；方形表示全局探索。右图比较种子间姿势，1Å以下为线性刻度，以上压缩跨度；2Å仅作位置一致性的描述参照。按基因排列，未把跨受体原始分数直接当成真实亲和力排序。','',
'## 本轮参照复对接','',
'![参照复对接](figures/02_redocking_controls.png)','',
'| 体系 | 最高排名姿势达到2Å | 返回姿势中至少一个达到2Å |','|---|---:|---:|']
for g in ['PPARG','PNP','SLC18B1']:
    rr=[r for r in controls if r['gene']==g];lines.append(f'| {g} | {sum(truth(r["top_pose_within2A"]) for r in rr)}/3 | {sum(truth(r["any_pose_within2A"]) for r in rr)}/3 |')
lines+=['','参照结果仅评价对应受体和设置；逐次RMSD、评分及计算方式见[参照表](results/Redocking_control_results.tsv)。','',
'## 五个全局探索与六个查新重点','']
for g in ['MFSD10','NHSL2','SLC18B1','TMEM254','TSPAN3','SAMD8']:
    r=next(r for r in evidence if r['gene']==g)
    lines.append(f'- **{g}：**{r["structure_support"]}。{r["specific_reason"]}。{r["suggested_next_step"]}。疾病背景：上调类型{r["SC_up_cell_types"] or "无"}；下调类型{r["SC_down_cell_types"] or "无"}；患者PB支持{r["PB_support_cell_types"] or "未列支持类型"}。')
lines+=['','五个严格PubMed0命中与SAMD8机构缩写消歧附列保留区分；查新信息用于讨论创新性，47个对象的计算范围不受它过滤。ATP5MC1、FXYD3、MT-ND4L、NHSL2、SMIM22的全局网格、姿势和质量数据均保留。','',
'## 方法与可复查文件','',
'使用Vina1.2.7、exhaustiveness64、seed20261024/25/26、最多20姿势。25个靶点采用实验结构，共23个不同PDB条目；其余22个采用AF模型。实际网格与结构选择在计算前固定。[完整固定方案](47靶点对接_计算方案.md)与[47实际网格](inputs/47_actual_search_regions.tsv)。','',
'受体逐残基映射，去融合段，按生物学组装保留目标和必要伙伴；模板补缺失重原子，已解析原子固定。金属与辅因子逐项留档。基本几何采用PoseBusters配置；受体/辅因子距离用全部保留重原子计算；种子RMSD和参照RMSD均在受体坐标系中作对称性校正，不先对齐配体。','',
'[282个最高排名姿势检查](results/All_completed_top_pose_checks.tsv)、[94组分状态汇总](results/Target_state_summary.tsv)、[种子RMSD](results/Seed_pose_pairwise_RMSD.tsv)。docking保留全部原始姿势、命令、输入SHA256和日志；poses保留全部返回姿势的三维SDF。','',
'方法依据：[Vina官方教程](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)、[Vina搜索空间与评分说明](https://autodock-vina.readthedocs.io/en/latest/faq.html)、[P2Rank作者文档](https://github.com/rdk/p2rank)。局部/全局搜索结果按其实际覆盖解释；当前研究选择的下一步是结合及药效依赖性实验。']
(B/'结果报告_47靶点反向对接.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
print('Created integrated47 evidence,full report and2 figure sets. Final checks still required.')
