import json,re,csv,time
from pathlib import Path
from runtime_io import read_json
ROOT=Path(__file__).resolve().parents[3];B=Path(__file__).resolve().parents[1]
if (B/'audit/ALL47_DOCKING_DELIVERY_COMPLETED.json').exists():
    print('Final delivery already recorded; preserve final documents.')
    raise SystemExit(0)
progress=read_json(B/'audit/Run_progress.json') if (B/'audit/Run_progress.json').exists() else {'main_completed':0,'control_completed':0,'failed':0,'running':[],'ready_receptors':0}
qa=read_json(B/'audit/Evaluation_progress.json') if (B/'audit/Evaluation_progress.json').exists() else {}
preps=[]
for p in (B/'audit').glob('*_receptor_preparation.json'):
    d=read_json(p)
    if d.get('exit_code')==0 and (B/'prepared/receptors'/f'{d["gene"]}.pdbqt').exists():preps.append(d['gene'])
stamp=time.strftime('%Y-%m-%d %H:%M:%S %z')
status=f'更新：{stamp}。受体已准备{len(preps)}/47；新药物主对接已完成{progress["main_completed"]}/282；参照复对接已完成{progress["control_completed"]}/9；当前运行{len(progress["running"])}项；运行失败记录{progress["failed"]}项。'
raw_done=progress['main_completed']==282 and progress['control_completed']==9 and not progress['running'] and not progress['failed']
phase=('原始计算已完成；完整报告和技术验收已生成，下一步只需核对两组实际PNG并完成最终文档。' if (B/'audit/ALL47_TECHNICAL_CHECKS_COMPLETED.json').exists()
       else '原始计算已完成；正在汇总姿势检查、报告和技术验收。' if raw_done else '当前主对接与检查持续运行。')
text=f'''# 全47靶点反向对接：本轮当前状态

{status}

{phase}

用户已确认：实验结构优先、AlphaFold补齐；每靶点一个主区域；ATP5MC1、FXYD3、MT-ND4L、NHSL2、SMIM22全局探索；药物+1/中性各3种子，共282次主对接，另做PPARG/PNP/SLC18B1的9次参照复对接。

旧两候选对接135份结果及对应旧结论已按用户确认清除，上游疾病、DeepPurpose、PubMed及公共结构保留。本轮结果从头计算，不混用旧分数。

## 方案与实际输入

- [固定计算方案](47靶点对接_计算方案.md)
- [47行结构与区域选择](inputs/47_structure_and_search_plan.tsv) · [实际网格坐标](inputs/47_actual_search_regions.tsv)
- [25个实验结构审计](inputs/25_experimental_structure_review.tsv) · [逐残基映射](inputs/Experimental_residue_mapping.tsv)
- [47行输入/疾病/预测证据](inputs/47_verified_input_inventory.tsv) · [药物结构核验](inputs/Ligand_input_audit.tsv)

## 执行与已有输出

- [实时进度JSON](audit/Run_progress.json) · [实际任务清单](inputs/All47_main_and_control_jobs.tsv)
- [计算日志](logs/all47_docking_6parallel.log) · [大部分受体准备日志](logs/prepare_all47_receptors_v2.log) · [AF准备日志](logs/prepare_AF_receptors.log)
- 原始本药姿势、命令和逐次日志：docking；独立参照复对接：controls；转换后的带键级三维SDF：poses。
- 已有质量检查及参照RMSD：results/All_completed_top_pose_checks.tsv、Redocking_control_results.tsv、Target_state_summary.tsv。它们是持续更新的中间表，未完成靶点不作最终评价。

## 输入准备的具体处理

SIFTS逐残基映射再与冻结序列核对；SLC18B1按mmCIF及确切序列限制在实际人源片段，去掉融合段。采用PDB生物学组装1，功能复合物保留天然人源伙伴；缺失长环不补建。PDBFixer1.12/OpenMM8.4补缺失重原子，已解析原子固定；不人为新增未解析链尾OXT。

Meeko0.7.1按标准氨基酸模板连接与参数化，避免把距离接近误认成化学键；肽键和二硫键用原共价半径/距离规则识别，大复合物使用等价直接查找减少两两遍历。该处理仅在本轮脚本进程内应用，没有改安装包或系统配置。源代码及逐靶点补原子记录在scripts、audit。

保留的Zn、血红素等使用实测坐标；CCD定义芳香性，Vina固定辅因子原子类型逐项记录于audit/Cofactor_preparation.json。使用Vina评分，不解释金属方向性或氧化还原反应。所有上下游保留检查以audit/Protected_inputs_before.tsv为基线。

## 清理与接手

[已删除文件清单](audit/Authorized_old_docking_deletion_manifest.tsv) · [删除完成核对](audit/Old_docking_deletion_completed.json)。根目录[新接手文件](../../新对话接手_当前状态.md)、[交接第28节及后续](../../工作流3_交接.md)同步本轮范围与实际状态。{phase}

用户最新要求：先把对接做完做好，额度可能用完，随时记录进度。本轮不再追问细胞系或展开实验条件。当前聊天中断时，已启动的本地队列与检查脚本仍可继续；新对话先读本页及audit/Run_progress.json，核对进程后再恢复，避免重复启动。已有ALL47_TECHNICAL_CHECKS_COMPLETED.json时，不重跑282次，先查看两组PNG，记录真实视觉检查后运行finalize_documents.py。
'''
(B/'README_当前状态.md').write_text(text,encoding='utf8')
p=ROOT/'AGENTS.md';s=p.read_text(encoding='utf8');s=s.replace('新对接尚未运行','新对接已启动，实际进度见新轮README与audit/Run_progress.json');p.write_text(s,encoding='utf8')
if raw_done:
    s=s.replace('新对接已启动，实际进度见新轮README与audit/Run_progress.json','282次主对接及9次参照已完成，报告和技术状态见新轮README；待实际图形复核及最终文档确认')
    p.write_text(s,encoding='utf8')
p=ROOT/'新对话接手_当前状态.md';s=p.read_text(encoding='utf8')
s=re.sub(r'新药物对接已完成0次；正在核对实验结构及搜索区域，不把准备完成写成对接完成。',
    '47个实际搜索区域已固定，用户已确认282次主方案；新药物对接与9次参照复对接正在运行。实时进度见新轮README及audit/Run_progress.json。',s)
s=re.sub(r'新一轮即时状态：.*?\n','',s)
if raw_done:
    s=s.replace('新药物对接与9次参照复对接正在运行。实时进度见新轮README及audit/Run_progress.json。','282次主对接及9次参照均已完成。'+phase)
s=s.replace('## 已完成且保留的上游证据',f'新一轮即时状态：{status}\n\n## 已完成且保留的上游证据')
p.write_text(s,encoding='utf8')
p=ROOT/'工作流3_交接.md';s=p.read_text(encoding='utf8')
s=s.replace('新轮目录outputs/reverse_docking47_2026-10-24，当前完成输入审计，尚无新对接分数（0次）。正在准备25个实验结构的具体覆盖/组装/配体核验及47个搜索区域清单。',
 '新轮目录outputs/reverse_docking47_2026-10-24，已完成25个实验结构的具体覆盖/组装/配体核验及47个实际搜索区域清单。用户进一步确认按282次主方案执行，正式计算与9次参照复对接已启动。')
marker='## 29. 2026-10-24：282次主方案已确认并执行（实时状态）'
if marker in s:s=s[:s.index(marker)].rstrip()+'\n\n'
s+=marker+f'''\n\n{status}

{phase} 用户最新明确先完成对接、减少额度消耗并随时记录；不再以细胞系名称作为前置，不展开新的实验条件讨论。

用户答复“按282次主方案执行”。结构25个实验、22个AF；每靶点一个主区域，5个无AF口袋者全局探索。两药物状态、3种子、exhaustiveness64，主对接282次；PPARG/PNP/SLC18B1各3种子参照另9次。实际任务和网格已固定，单次运行保存输入/引擎SHA256及完整命令，可复用完成结果。按16逻辑CPU使用6个任务×CPU2，大网格仅允许1个同时运行。一次并发调整中停止了命令行核对的本轮未完成子任务，已完成结果保留，记录见audit/Concurrency_reconfiguration_processes.json。

本轮输入准备新增PDBFixer1.12/OpenMM8.4于隔离环境，补缺失重原子而不补长环；已解析原子固定。Meeko在本轮进程内按标准模板确定残基内部连接，肽键/二硫键用原距离阈值直接查找，解决错误近距离连接及大复合物两两遍历耗时；安装包源码和系统环境未改。Zn/血红素等按实际坐标与CCD定义保留，处理清单在audit/Cofactor_preparation.json。所有技术修复在新scripts和逐靶点audit可追溯。

恢复先读新轮README、audit/Run_progress.json及logs/all47_docking_6parallel.log；受体准备见prepare_all47_receptors_v2.log、prepare_AF_receptors.log。evaluate_completed.py持续核对完成姿势的立体信息、基本几何、原子距离、目标/伙伴接触、种子RMSD及参照对实验RMSD，结果未齐前不生成最终候选结论或完成标记。旧两候选结果已清除，不用旧分数补齐。上游数据及47原池保留，MD未开展。
'''
p.write_text(s,encoding='utf8')
print(status)
