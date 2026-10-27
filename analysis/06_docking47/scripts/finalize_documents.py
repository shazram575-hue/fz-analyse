"""Run only after the actual figures have been visually reviewed."""
import csv,json,hashlib,time,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];B=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
validation=json.loads((B/'audit/ALL47_TECHNICAL_CHECKS_COMPLETED.json').read_text(encoding='utf8'))
assert validation['technical_delivery_complete']
assert (B/'audit/Figure_visual_QA.md').exists(),'Review both actual PNG figures before finalizing.'
stamp=time.strftime('%Y-%m-%d %H:%M:%S %z')
report=B/'结果报告_47靶点反向对接.md'
summary=f'47个靶点×两状态×三种子共282次主对接及9次参照复对接全部完成；{validation["technical_checks_passed"]}/{validation["technical_checks_total"]}项技术检查通过。'
top=f'''# 全47靶点反向对接已完成

更新：{stamp}。{summary}

优先阅读[47靶点完整结果报告](结果报告_47靶点反向对接.md)与[47行结构/疾病/预测整合证据](results/47_integrated_docking_disease_evidence.tsv)。报告逐项回答结构支持、具体问题和建议的下一步，唯一主攻仍待与用户讨论。

- [94组分状态汇总](results/Target_state_summary.tsv) · [282姿势质量检查](results/All_completed_top_pose_checks.tsv)
- [9次参照对实验RMSD](results/Redocking_control_results.tsv) · [种子间姿势RMSD](results/Seed_pose_pairwise_RMSD.tsv)
- [固定方案](47靶点对接_计算方案.md) · [47实际结构与网格](inputs/47_actual_search_regions.tsv)
- [最终技术验收](audit/Final_independent_validation.tsv) · [图形检查](audit/Figure_visual_QA.md)
- [旧135份对接文件删除记录](audit/Authorized_old_docking_deletion_manifest.tsv) · [上游保护检查](audit/Final_validation_summary.json)

原始本药姿势、命令与日志在docking；参照在controls；全部返回姿势的三维SDF在poses；输入、序列映射、修复与辅因子记录在inputs/prepared/audit。公共软件与结构沿用保留目录，旧两候选分数和对应旧结论已清除，没有混入新结果。

当前本轮主计算与检查结束，不重跑已完成282次，不自动增加口袋、分子动力学或锁定唯一靶点。恢复先读本页和根目录新接手状态；需要重建表图时使用本目录脚本，不使用旧两候选报告脚本。
'''
(B/'README_当前状态.md').write_text(top,encoding='utf8')
p=ROOT/'新对话接手_当前状态.md';s=p.read_text(encoding='utf8')
s=s.replace('更新：2026-10-24。课题',f'更新：{time.strftime("%Y-%m-%d")}。课题')
s=s.replace('原始计算已完成；完整报告和技术验收已生成，下一步只需核对两组实际PNG并完成最终文档。','原始计算、技术验收和实际图形复核均完成，最终文档已同步。')
s=re.sub(r'新一轮即时状态：.*?\n',f'新一轮即时状态：{stamp}。{summary}\n',s)
s=s.replace('47个实际搜索区域已固定，用户已确认282次主方案；新药物对接与9次参照复对接正在运行。实时进度见新轮README及audit/Run_progress.json。',
    f'47个实际搜索区域已固定，用户确认的282次主方案和9次参照均已完成。{validation["technical_checks_passed"]}/{validation["technical_checks_total"]}项技术检查通过；新轮报告已生成并完成两组结果图检查。')
if '本轮[最终47靶点报告]' not in s:
    s=s.replace('## 已完成且保留的上游证据','本轮[最终47靶点报告](outputs/reverse_docking47_2026-10-24/结果报告_47靶点反向对接.md)及[整合证据表](outputs/reverse_docking47_2026-10-24/results/47_integrated_docking_disease_evidence.tsv)为候选讨论的新依据。当前无本轮待运行主计算。\n\n## 已完成且保留的上游证据')
p.write_text(s,encoding='utf8')
p=ROOT/'AGENTS.md';s=p.read_text(encoding='utf8')
s=re.sub(r'\A> .*?\n\n',f'> **当前优先入口（{stamp}）：**用户已确认的全47反向对接已完成：282次主计算＋9次参照，{validation["technical_checks_passed"]}/{validation["technical_checks_total"]}项技术验收通过，两组结果图已检查。旧两候选135份对接输出及旧结论已清除，上游疾病、DeepPurpose、PubMed及公共结构保留。先读[新接手状态](新对话接手_当前状态.md)、交接第30节与[本轮47报告](outputs/reverse_docking47_2026-10-24/结果报告_47靶点反向对接.md)。唯一主攻待讨论，不重跑主体或自动加口袋/MD。下方历史范围以本段及用户新指令为准。\n\n',s,count=1,flags=re.S)
p.write_text(s,encoding='utf8')
p=ROOT/'工作流3_交接.md';s=p.read_text(encoding='utf8')
final_marker='## 30. 全47反向对接正式交付（本轮完成）'
if final_marker in s:s=s[:s.index(final_marker)].rstrip()+'\n'
s=s.replace('原始计算已完成；完整报告和技术验收已生成，下一步只需核对两组实际PNG并完成最终文档。','原始计算、技术验收和实际图形复核均完成，最终文档已同步。')
s=re.sub(r'\A> .*?\n\n', '> **当前总入口：**用户确认的全47反向对接已完成，优先读[新接手状态](新对话接手_当前状态.md)及第30节。旧两候选结果按用户要求清除，新结果只引用本轮282次主计算及9次参照。\n\n',s,count=1,flags=re.S)
s+='''\n\n## 30. 全47反向对接正式交付（本轮完成）

'''+f'更新：{stamp}。{summary}\n\n'+'''用户先确认仅删除旧对接及对应结论、保留上游，采用实验结构优先/AF补齐/5个无口袋者全局探索；随后确认282次主方案。实际25个实验结构、22个AF，所有47均已完成两状态三种子计算，未用旧分数补齐。PPARG/PNP/SLC18B1参照各3次重新执行，成败按实际RMSD逐项报告，不把一个体系结果推广到其他体系。

旧135份输出及对应MD结论、死链接、旧结果生成脚本和空目录均已清除。保留输入402份按SHA256核验不变。CA12作者残基有14个插入码，映射已补全完整SeqId并重新准备，其他实验结构无此类插入码；SLC18B1融合段按确切序列/原mmCIF边界去除。MT-ND2/ND4/ND4L保留同一天然复合物，只有原始蛋白原子记录完全一致时复用受体参数化，药物对接按各自区域独立运行，没有复用对接分数。

本轮PoseBusters基本几何、九中心/身份、全受体及辅因子原子距离、目标/伙伴接触、局部模型可信度、三种子RMSD均已记录。技术验收完成与科学支持程度分别报告，质量失败或参照失败保留在47行证据中。图形检查见新目录audit/Figure_visual_QA.md。未运行MD，未改变原813、1,881或47入池规则，未据原始Vina分数唯一锁靶。

最终入口：outputs/reverse_docking47_2026-10-24/结果报告_47靶点反向对接.md、results/47_integrated_docking_disease_evidence.tsv、README_当前状态.md。raw docking/controls、三维poses、prepared输入、逐残基映射、原始日志、软件版本与SHA256均保留。当前主计算与质量检查结束，下一步与用户围绕完整47证据确定是否推进及主攻/备选，不擅自扩大计算。
'''
s=s.replace('实际25个实验结构、22个AF','实际25个靶点采用23个不同PDB条目、22个靶点采用AF')
s=s.replace('CA12作者残基有14个插入码','CA12两链共有14条带插入码的残基映射')
p.write_text(s,encoding='utf8')
p=ROOT/'outputs/README_结果索引.md';s=p.read_text(encoding='utf8')
if s.startswith('# 最新完成：47靶点反向对接') and '以下为历史阶段索引。\n\n---\n\n' in s:
    s=s.split('以下为历史阶段索引。\n\n---\n\n',1)[1]
s='# 最新完成：47靶点反向对接\n\n'+summary+'\n\n[完整报告](reverse_docking47_2026-10-24/结果报告_47靶点反向对接.md) · [47行整合证据](reverse_docking47_2026-10-24/results/47_integrated_docking_disease_evidence.tsv) · [当前入口](reverse_docking47_2026-10-24/README_当前状态.md)。\n\n以下为历史阶段索引。\n\n---\n\n'+s
p.write_text(s,encoding='utf8')
manifest=[]
for p in sorted(B.glob('*.md')):
    manifest.append({'path':str(p.relative_to(B)),'bytes':p.stat().st_size,'sha256':sha(p)})
for p in [B/'audit/Figure_visual_QA.md',B/'audit/environment_pip_freeze.txt',ROOT/'新对话接手_当前状态.md',ROOT/'工作流3_交接.md',ROOT/'AGENTS.md',ROOT/'outputs/README_结果索引.md']:
    manifest.append({'path':str(p.relative_to(B)) if p.is_relative_to(B) else '../../'+str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
for folder in ['inputs','results','figures','scripts','docking','controls','poses','prepared','sources','structures']:
    for p in sorted((B/folder).rglob('*')):
        if p.is_file():manifest.append({'path':str(p.relative_to(B)),'bytes':p.stat().st_size,'sha256':sha(p)})
with (B/'audit/Deliverable_SHA256_manifest.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
validation['visual_QA_completed']=True;validation['deliverable_files_hashed']=len(manifest)
validation['final_delivery_timestamp']=time.strftime('%Y-%m-%dT%H:%M:%S%z')
(B/'audit/ALL47_DOCKING_DELIVERY_COMPLETED.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf8')
print(summary)
