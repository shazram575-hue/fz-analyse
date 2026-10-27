import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]; B=Path(__file__).resolve().parents[1]
def table(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
inv=table(B/'inputs/47_verified_input_inventory.tsv')
exp={r['gene']:r for r in table(B/'inputs/25_experimental_structure_review.tsv')}
LIG={'ADRA2A':'J59','CA12':'AZM','CYBRD1':'ASC','GPR132':'140','ITPR3':'I3P',
     'NR3C1':'CV7','PNP':'JTC','PPARG':'BRL','PTGER4':'7UR','S1PR1':'J8C',
     'SAMD8':'Z0P+OPE','SCD':'ST9','SLC18B1':'YHL','SSTR2':'GI9','VKORC1':'UA7'}
COMPLEX={'APH1A','ATP5MC1','LRRC26','MT-ND2','MT-ND4','MT-ND4L'}
HOMO={'ACKR1','ITPR3','KCNK1','PNP','SAMD8'}
rows=[]
for r in inv:
    g=r['gene'];e=exp.get(g)
    global_search=int(r['P2Rank_pocket_count'])==0
    site='global_search' if global_search else ('experimental_ligand:'+LIG[g] if g in LIG else 'P2Rank_first_target_associated_pocket')
    rows.append({'gene':g,'accession':r['accession'],'source':'PDB' if e else 'AlphaFold',
      'structure':e['pdb'] if e else r['accession']+'_AF','target_chains':e['canonical_chains'] if e else 'A',
      'resolved_residues_per_best_chain':e['resolved_best_chain_residues'] if e else r['frozen_length'],
      'canonical_span':e['canonical_min']+'-'+e['canonical_max'] if e else '1-'+r['frozen_length'],
      'assembly_policy':'target_and_native_human_complex_partners' if g in COMPLEX else ('target_homomer' if g in HOMO and e else 'target_chain'),
      'site_rule':site,'known_ligand_CCD':LIG.get(g,''),
      'sequence_differences':e['mutations'] if e else '',
      'AF_fraction_below50':r['AF_fraction_residues_pLDDT_lt50'],
      'new_run_status':'protocol_confirmed_282_main_runs'})
with (B/'inputs/47_structure_and_search_plan.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
lines=['# 47靶点重新对接方案（用户已确认，运行前固定）','',
'本轮固定原47集合；实验结构优先，AlphaFold补齐；5个未检出AF口袋者全局探索。下表尚非对接结果。','',
'## 配体与重复','',
'Spiropachysine A采用已核准CID209254，+1为主状态，中性状态作敏感性；九中心保持一致。Vina1.2.7、exhaustiveness64、3个随机种子20261024/25/26、20个返回姿势、energy_range5；两状态分别报告，不取两者更低者代替主状态。47×2×3＝282次主对接。','',
'## 结构与搜索区域','',
'优先保留实验坐标及其实际序列，去除融合段、表达标签和结晶辅助抗体。实验工程差异逐一列明，保持原坐标，不静默改回野生型。缺失重原子用模板补齐；缺失长环和未解析残基不凭空补建。涉及功能复合物者保留目标相关的天然人源伙伴，姿势统计分别记录目标链和伙伴接触。','',
'有共结构小分子/底物者，按其实验位置确定一个主搜索区域；这不预设药物一定竞争该位点。其余用P2Rank在实际采用结构上找目标相关第一口袋。ATP5MC1、FXYD3、MT-ND4L、NHSL2、SMIM22按用户决定全局探索。局部盒按位点重原子包围范围每侧外扩6Å、最小24Å；全局盒覆盖全部目标重原子并留配体空间。记录中心、边长及区域体积。','',
'全局探索与局部位点对接分组展示；尤其NHSL2全长模型约75%残基pLDDT<50，其全局分数同时标注实际接触区域可信度。若后续需要细分局部搜索，先列具体问题与工作量再讨论。','',
'Zn/血红素等结构性金属或辅因子按位点作用保留，结晶溶剂及待替换的实验配体移除；逐项记录。对金属/辅因子体系单列处理方式。','',
'## 验证与交付','',
'在本轮重新做已知小分子参照的复对接并报告全部重复；参照检验与本药282次计算分开。每个靶点报告：两状态三种子分数中位数及范围、姿势重复性、几何/碰撞、接触残基与局部可信度、疾病/细胞背景/预测/查新证据，以及是否值得进入结合和药效依赖性实验的明确理由。','',
'保留47行，不按分数先删对象。跨受体Vina分数用于展示计算结果，后续主攻选择结合位点依据、可重复性、疾病证据和实验可行性讨论；本轮不自动唯一锁靶，不做MD。','',
'具体结构审计见[25实验结构审阅](inputs/25_experimental_structure_review.tsv)，逐残基映射见[映射表](inputs/Experimental_residue_mapping.tsv)，可筛选的47行方案见[方案TSV](inputs/47_structure_and_search_plan.tsv)。','',
'| 基因 | 结构 | 已解析残基/链 | 搜索区域依据 |','|---|---|---:|---|']
for r in rows:lines.append(f'| {r["gene"]} | {r["structure"]} | {r["resolved_residues_per_best_chain"]} | {r["site_rule"]} |')
lines+=['','方法依据：[Vina官方教程](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)、[P2Rank作者文档](https://github.com/rdk/p2rank)。参数为本轮拟定设置，计算前固定。']
(B/'47靶点对接_计算方案.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
(B/'audit/Protocol_authorized.json').write_text(json.dumps({'date':'2026-10-24','user_reply':'按282次主方案执行','main_targets':47,'ligand_states':['pH74_plus1','neutral_sensitivity'],'seeds':[20261024,20261025,20261026],'exhaustiveness':64,'max_poses':20,'main_runs':282,'controls':'separate redocking runs','structure_preference':'experimental_first_AF_fallback','global_search_genes':['ATP5MC1','FXYD3','MT-ND4L','NHSL2','SMIM22']},ensure_ascii=False,indent=2),encoding='utf8')
print('Protocol draft:47 targets,25 experimental,22 AF,282 main docking jobs.')
