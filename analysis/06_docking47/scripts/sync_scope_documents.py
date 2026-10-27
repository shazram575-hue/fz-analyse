import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
B=Path(__file__).resolve().parents[1]
if (B/'audit/Run_progress.json').exists():
    raise SystemExit('Initial scope sync already applied. Use update_runtime_status.py; do not restore the initial zero-run snapshot.')
def save(p,s): p.write_text(s.strip()+'\n',encoding='utf8')
newhead='> **当前执行范围（2026-10-24，用户最新确认）：**47个前10%两终点共有∩疾病靶点全部与Spiropachysine A重新对接。实验结构优先、AlphaFold补齐；未检出口袋的5个也做全局探索并单列。用户确认仅清除旧对接结果及对应旧结论，保留原始数据、疾病分析、DeepPurpose、PubMed、公共结构和环境。旧两候选输出135份已删除，删除清单在outputs/reverse_docking47_2026-10-24/audit；新对接尚未运行。当前入口[新对话接手_当前状态](新对话接手_当前状态.md)及交接第28节。下文早期“不扩47”“首选SAMD8”等安排已被本次指令覆盖。'
p=ROOT/'AGENTS.md';s=p.read_text(encoding='utf8')
s=re.sub(r'\A> .*?\n\n',newhead+'\n\n',s,count=1,flags=re.S)
s=re.sub(r'## 2026-10-24查新复核与配对对接最新范围.*?(?=## 2026-10-24最新：PubMed)',
'''## 2026-10-24旧两候选对接清除与全47重做

用户最新明确要求全47对接，并确认只删除旧对接结果和对应结论。旧限两候选规则、旧对接分数/姿势/图、原优先建议及旧完成标记已清除。疾病分析、预测、PubMed记录、公共结构、输入准备及环境保留。正式结构/区域清单在新目录准备中，执行状态见交接第28节。

''',s,flags=re.S)
save(p,s)

current='''# 工作流3：全47靶点反向对接的接手状态

更新：2026-10-24。课题Spiropachysine A＋乳腺癌，目标Phytomedicine。

## 最新授权与实际状态

用户要求读取相关交接和具体信息后，对47个前10%两终点共有∩疾病靶点全部与药物对接，评估是否支持后续研究。用户已确认：实验结构优先，缺少合适实验结构时AlphaFold补齐；5个未检出口袋者仍做全局探索，结果单列。

用户确认仅删除旧对接结果及对应结论。已按逐文件SHA256清单删除135份旧姿势、分数、图、报告、配套完成标记、旧快照与旧两候选结果脚本；未把副本留作当前结果。原始数据、疾病分析、DeepPurpose预测、PubMed、公共结构、输入准备和环境保留。

新轮当前完成输入核对：47个AF模型逐位匹配冻结蛋白序列；25个有PDB条目，42个检出AF口袋。ATP5MC1、FXYD3、MT-ND4L、NHSL2、SMIM22的AF未检出口袋。新药物对接已完成0次；正在核对实验结构及搜索区域，不把准备完成写成对接完成。

入口：[新轮README](outputs/reverse_docking47_2026-10-24/README_当前状态.md)、[47行输入核验](outputs/reverse_docking47_2026-10-24/inputs/47_verified_input_inventory.tsv)、[删除清单](outputs/reverse_docking47_2026-10-24/audit/Authorized_old_docking_deletion_manifest.tsv)、[删除核对](outputs/reverse_docking47_2026-10-24/audit/Old_docking_deletion_completed.json)。

## 已完成且保留的上游证据

| 阶段 | 实际结果 | 入口 |
|---|---|---|
| 第一阶段 | Pool A813；CNV主分类采用17c coherent final calls | [结果报告](第一阶段_结果报告.md)、[冻结方案](第一阶段_分析方案与验收.md) |
| 单细胞补全 | 6生物学大类＋2诊断组，同类型Tumor-origin−Healthy双向DE；20项检查、22对图 | [补全报告](单细胞基础补全_结果报告.md) |
| 疾病主输入 | 六大类细胞级DE并集∩Bulk∩protein_coding＝1,881；诊断额外76另列 | [输入](outputs/DeepPurpose_预测_2026-10-24/main_input/Disease_matching_candidates_1881.tsv) |
| 早期池内预测 | 1,993配对×10模型＝19,930评分 | [报告](outputs/DeepPurpose_预测_2026-10-24/DeepPurpose_结果报告.md) |
| 独立全蛋白背景预测 | 20,431条×8非MPNN模型＝163,448评分；33/33验收 | [报告](outputs/身份深核与独立交集_2026-10-24/独立预测与疾病交集_结果报告.md) |

用户采用前10%范围，两终点共有∩疾病＝47；1%/5%及前10%并集296对照保留。Kd与IC50各4编码独立共识，共用BindingDB；MPNN因批次依赖单列。47为本轮完整对接对象，不再按查新先缩成两候选。[原47池](outputs/身份深核与独立交集_2026-10-24/results/disease_intersections/top10percent_Kd_intersection_IC50_disease1881.tsv)。

药物采用CID209254，完整InChIKey SLGWGPQWJRVPAD-PZOAWPBASA-N；原始合成研究Figure1的九中心重建与其一致，不再等待结构核准。[身份核对](outputs/身份深核与独立交集_2026-10-24/身份核对结论.md)。本轮重新核验既有中性、+1三维输入及九中心，结果见[配体核验](outputs/reverse_docking47_2026-10-24/inputs/Ligand_input_audit.tsv)。

## PubMed资料保留

47项完整领域检索47/47有效。严格0为MFSD10/SLC75A1、NHSL2、SLC18B1、TMEM254、TSPAN3；SAMD8补SMSr后1条为作者单位误命中，靶点相关0，单列。MT-ND4L有PMID24061460、SMIM22/CASIMO1有PMID29765154，均不属严格0。

[最新查新结论](outputs/docking_audit_2026-10-24/PubMed零命中标准_复核结论.md)、[47项完整查询与计数](outputs/docking_audit_2026-10-24/pubmed_comprehensive_recheck/47_comprehensive_PubMed_counts.tsv)。本轮对接不按文献热度跳过其余候选。

## 具体执行注意

- MFSD10按Q14728映射现名SLC75A1。MT-ND4结构用459aa P03905，额外99aa C0HME5的上游模型结果保留。
- NHSL2全长1225aa、约75%残基pLDDT<50；ITPR3全长2671aa。全局搜索区域、分域或分箱必须明示，不能静默截前1000aa。
- APH1A、ATP5MC1、LRRC26及线粒体ND亚基涉及复合物，结构核对需记录伙伴与目标链，结果逐一交代接触归属。
- 用户已确认可开展敲低/过表达、CETSA/DARTS及纯化蛋白结合检测，具体细胞系名称尚未给出。无需再次询问已确认能力。
- 本轮只重做对接，原813、1,881和47筛选依据不重算；未授权MD。

现有环境：outputs/docking_audit_2026-10-24/env/Scripts/python.exe；Vina1.2.7、Meeko0.7.1、RDKit2024.09.6、Molscrub0.3.2、P2Rank2.5.1、Java17。公共结构、序列映射和软件可复用；旧两候选结果已经清除，不再引用其数值或推进结论。新轮脚本和日志统一写outputs/reverse_docking47_2026-10-24。
'''
save(ROOT/'新对话接手_当前状态.md',current)

p=ROOT/'工作流3_交接.md';s=p.read_text(encoding='utf8')
s=re.sub(r'\A> .*?\n\n', '> **当前入口（2026-10-24）：**用户已授权全47重新对接并确认清除旧对接结果，优先读[新接手状态](新对话接手_当前状态.md)及第28节。第22–24节旧对接数值与结论按用户要求清除，第27节为前次接手记录。\n\n',s,count=1,flags=re.S)
s=re.sub(r'## 22\..*?(?=## 25\.)', '''## 22–24. 旧两候选对接及解释已清除（用户2026-10-24明确要求）

本段原两候选对接输出、分数、图和推进结论已随旧文件清除。原查新材料、公共结构和输入准备保留，最新PubMed事实见第25–26节。用户确认的敲低/过表达、CETSA/DARTS和纯化蛋白结合检测能力保留为实验条件信息。全47新对接见第28节。

''',s,flags=re.S)
s=re.sub(r'## 27\..*\Z', '''## 27. 前次接手记录已由全47重做安排更新

前次仅两候选且不扩47的范围已经被用户新指令覆盖。原两候选对接结果和对应旧结论按用户要求清除；上游数据、预测、PubMed检索保留。完整当前入口见根目录《新对话接手_当前状态.md》。

## 28. 2026-10-24：用户要求全47反向对接，旧对接清除、新轮准备中

用户要求完整接手相关交接和具体信息，对47个筛选出的共有靶点全部与Spiropachysine A进行对接，评估是否支撑后续研究，同步清理旧结果并修改MD。用户已明确答复两个问题：（1）只清除旧对接结果及对应旧结论，原始数据、疾病分析、DeepPurpose、PubMed、公共结构和软件环境保留；（2）实验结构优先，AlphaFold补齐，5个无预测口袋者也做全局探索。

已读取根接手文件、全交接相关阶段、冻结方案、独立预测/身份/疾病输入、反向对接调研、结构准备和旧执行脚本。已独立核对47原池与结构表集合相同、47个AF逐位匹配冻结蛋白序列、口袋CSV数量相同、中性/+1配体三维及九中心完整。现有25个有PDB条目、42个有AF口袋；5个无AF口袋为ATP5MC1、FXYD3、MT-ND4L、NHSL2、SMIM22。

已按绝对路径边界与逐文件SHA256检查删除135份旧两候选对接输出及相应报告、旧完成标记、图、旧快照和旧结果生成脚本。删除清单及核对在新目录audit/Authorized_old_docking_deletion_manifest.tsv、Old_docking_deletion_completed.json；原始数据未改，保留部分已建立Protected_inputs_before.tsv。旧对接结果不再作为新轮结论。

新轮目录outputs/reverse_docking47_2026-10-24，当前完成输入审计，尚无新对接分数（0次）。正在准备25个实验结构的具体覆盖/组装/配体核验及47个搜索区域清单。涉及研究取舍的问题集中向用户说明，常规技术直接解决；正式运行和完成状态将逐次追加，不能以脚本或结构下载代替计算完成。

恢复入口：新目录README_当前状态.md、inputs/47_verified_input_inventory.tsv、inputs/Ligand_input_audit.tsv、audit/Intake_summary.json。本轮不重跑疾病/DeepPurpose，不改变47原池，不自动做MD。PubMed严格0仍5个、SAMD8机构误命中单列；查新信息作为47行证据字段，不再作为只计算两者的过滤条件。
''',s,flags=re.S)
save(p,s)

save(B/'README_当前状态.md','''# 47靶点反向对接：新轮当前状态

用户2026-10-24确认实验结构优先、AlphaFold补齐，5个无预测口袋者仍做全局探索。47个共有靶点全部纳入，逐项给出结果及是否支持后续研究的依据。

实际完成：47个输入的序列、结构库存和口袋数量核对；药物中性/+1三维与九中心核验；旧对接135份输出及对应旧结论已清除。新对接完成0次，具体实验结构与搜索区域准备中。

- [47行输入核验与疾病/预测信息](inputs/47_verified_input_inventory.tsv)
- [配体输入核验](inputs/Ligand_input_audit.tsv)
- [输入核对汇总](audit/Intake_summary.json)
- [已执行删除清单](audit/Authorized_old_docking_deletion_manifest.tsv) · [删除核对](audit/Old_docking_deletion_completed.json)
- [保留输入SHA256基线](audit/Protected_inputs_before.tsv)

环境与公共结构复用旧目录；本轮脚本/新结果均写本目录。根目录[新接手文件](../../新对话接手_当前状态.md)和[总交接第28节](../../工作流3_交接.md)已同步。
''')
save(ROOT/'outputs/docking_audit_2026-10-24/README_当前入口.md','''# 保留资料入口：PubMed、公共结构及软件环境

用户已要求全47重新对接。旧两候选对接结果、图、报告和对应旧结论已清除；新轮转到[47靶点对接入口](../reverse_docking47_2026-10-24/README_当前状态.md)。

本目录继续保留：

- [最新PubMed完整领域检索](PubMed零命中标准_复核结论.md)
- [47项计数、查询与PMID](pubmed_comprehensive_recheck/47_comprehensive_PubMed_counts.tsv)
- [47结构库存及口袋准备](results/47_structure_pocket_assessment.tsv)
- structures、sources、pockets、prepared中的公共结构及输入准备，以及env、software环境。

不再以此前两候选结果指定主攻。删除清单见新轮audit，正式对接状态以新轮结果为准。
''')
save(ROOT/'outputs/docking_audit_2026-10-24/查新复核与候选结论.md','''# 查新复核：保留事实与当前用途

2026-10-24更新。当前标准为PubMed完整乳腺癌领域词＋正确靶点名称/别名联合检索0命中。47项扩充查询均有效，严格0为MFSD10/SLC75A1、NHSL2、SLC18B1、TMEM254、TSPAN3；SAMD8原始1条为SMSR作者单位误命中，靶点相关0，单列。MT-ND4L及SMIM22/CASIMO1有真实命中。

完整结论、检索式与原始响应见[最新PubMed复核](PubMed零命中标准_复核结论.md)。查新用于标注创新性，本轮47个靶点均进入对接，不再按此前两候选安排缩减。

原附件47查询计数可以复现，167篇PMID与题名真实，但部分归属有误：MALL不能混同MAL/MAL2；NHSL2不能混同KIAA2022/NEXMIF；SMIM22必须补CASIMO1。具体记录保留在[47项修订表](results/47_literature_evidence_revised.tsv)、[167篇归属核对表](results/167_representative_paper_attribution_audit.tsv)，PubMed当前计数以完整领域检索为准。

原基于两候选旧对接的优先建议和位点解释按用户要求清除。新对接与选择讨论见[新轮入口](../reverse_docking47_2026-10-24/README_当前状态.md)。
''')
p=ROOT/'outputs/README_结果索引.md';s=p.read_text(encoding='utf8');idx=s.find('---')
tail=s[idx:] if idx>=0 else ''
save(p,'''# 当前结果入口：全47反向对接重做（2026-10-24）

先读[当前接手状态](../新对话接手_当前状态.md)、[交接第28节](../工作流3_交接.md)和[47靶点新轮入口](reverse_docking47_2026-10-24/README_当前状态.md)。用户已确认实验结构优先、AlphaFold补齐、无口袋者也全局探索。旧两候选对接输出与结论已清除。

[保留的PubMed及公共结构资料](docking_audit_2026-10-24/README_当前入口.md)。原始数据、第一阶段、单细胞、DeepPurpose及47筛选依据保留。下方为上游阶段的历史索引，执行范围以本段及新轮入口为准。

'''+tail)
save(ROOT/'outputs/反向对接_调研_2026-10-24/README_新对话接手.md','''# 当前接手入口

用户最新已授权全47对接，采用实验结构优先、AlphaFold补齐，未检出口袋者也全局探索。旧两候选结果及对应结论已清除。

请读[根接手状态](../../新对话接手_当前状态.md)、[总交接第28节](../../工作流3_交接.md)、[新轮入口](../reverse_docking47_2026-10-24/README_当前状态.md)。本目录方法调研保留为技术参考，其早期三候选优先安排已经被本轮范围覆盖。
''')
print('Updated scope and handoff documents; obsolete docking conclusions removed.')
