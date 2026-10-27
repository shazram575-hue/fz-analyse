"""Write a factual report after inference, summaries and independent validation pass."""
import hashlib
import json
import pathlib
import pandas as pd

ROOT=pathlib.Path(__file__).resolve().parents[1]

def displayed_table(frame,endpoint):
    lines=["|基因|面板内排名|模型百分位中位数|模型百分位IQR|去掉一种编码的最大排名变化|",
           "|---|---:|---:|---:|---:|"]
    for row in frame.head(10).itertuples(index=False):
        lines.append(f"|{row.gene_symbol}|{row.consensus_rank}|{row.consensus_percentile:.3f}|{row.model_percentile_IQR:.3f}|{row.leave_one_model_out_max_abs_rank_shift:.0f}|")
    return "\n".join(lines)

def main():
    results=ROOT/"results"
    metrics=json.loads((results/"Prediction_summary_metrics.json").read_text(encoding="utf-8"))
    validation=json.loads((ROOT/"audit"/"Final_validation_summary.json").read_text(encoding="utf-8"))
    assert validation["n_checks"]==validation["n_pass"] and not validation["failed"]
    kd=pd.read_csv(results/"Main_Kd_consensus_ranked.tsv",sep="\t")
    ic50=pd.read_csv(results/"Main_IC50_consensus_ranked.tsv",sep="\t")
    stability=pd.read_csv(results/"Leave_one_encoding_out_stability.tsv",sep="\t")
    scope=json.loads((ROOT/"protein_mapping_preparation"/"Resolved_scope_summary.json").read_text(encoding="utf-8"))
    main_scope=next(item for item in scope if item["scope"]=="main")
    rho=metrics["main_Kd_vs_IC50_consensus_Spearman"]
    overlap=len(set(kd.head(20).canonical_id)&set(ic50.head(20).canonical_id))
    text=f'''# DeepPurpose探索性多模型预测：CID209254与1,881疾病候选

日期：2026-10-24；项目Spiropachysine A＋乳腺癌。状态：本轮原始评分、终点分开汇总和数据一致性验收已完成；技术验收不是预测准确性、直接结合或实验机制验证。未选定最终靶点。

## 1. 实际完成结果

- 用户确认的主输入：六个生物学大类细胞级双向DE∩既定Bulk探索DE∩protein_coding，**1,881基因**；诊断组额外76基因另作敏感性，合计1,957。
- 主输入1,881基因对应{main_scope['gene_protein_rows']:,}条去重参考蛋白序列；敏感性76条，总计**{metrics['protein_pairs']:,}个化合物—蛋白配对**。所有主输入基因均有参考序列，没有因未映射而偷偷删基因。
- 10个作者预训练模型，共**{metrics['protein_model_predictions']:,}条蛋白—模型原始评分**。按基因聚合为{metrics['gene_model_predictions']:,}条基因—模型评分；蛋白原分数及范围保留。
- Kd和IC50各用4个非MPNN编码形成独立排名，共8个主要模型；2个MPNN只作逐对运行敏感性。
- 独立一致性检查**{validation['n_pass']}/{validation['n_checks']}通过**；7份第一阶段冻结结果SHA256保持。

主要入口：[逐基因疾病与模型证据总表](results/Main1881_combined_disease_model_evidence.tsv)、[Kd面板排序](results/Main_Kd_consensus_ranked.tsv)、[IC50面板排序](results/Main_IC50_consensus_ranked.tsv)。两种终点的模型百分位共识Spearman为**{rho:.3f}**；各自前20展示名单交集**{overlap}/20**。这些是本轮计算的一致性描述，不是外部验证。

## 2. 具体分子身份与命名边界

本课题名称沿用Spiropachysine A；实际计算分子明确为**PubChem CID209254（Spiropachysine）**，InChIKey `SLGWGPQWJRVPAD-PZOAWPBASA-N`，C31H46N2O。9个可识别手性中心全指定，结构文件写入/读回保持同一InChIKey。

原始分离/结构资料支持早期spiropachysine、后来的B和ISO应分别区分。B为16β-羟基衍生物；ISO是涉及C-3构型的非对映异构体。尚未找到原文明确作“A”为无后缀名同义词的直接定义。用户已知此边界后明确接受CID209254进行计算，故本轮使用的是透明、固定的**文献推定＋用户采用的结构映射**，不宣称独立证明Fang2022实物样品身份。详见[命名证据](nomenclature_review/README_命名核对与采用决定.md)及[结构核验](structure_verification_CID209254/README_核验结论.md)。

用户/官方原字序SMILES长101字符，CNN药物编码上限100；采用经InChIKey验证等价的RDKit规范异构SMILES92字符，避免截断。未删除立体标记、改盐型或优化构型。实际字符串、权重哈希及预先固定的汇总规则见[predictions/Frozen_prediction_design.json](predictions/Frozen_prediction_design.json)。

## 3. 疾病候选与蛋白输入

本轮沿用已完成的47主文库、45供者键、236,309细胞的同类型Tumor-origin−Healthy补全结果。单细胞使用absFC>0.5、max(pct.tumor,pct.healthy)>0.15、全33,487特征Bonferroni<0.05；Bulk实际为TCGA STAR RNA-seq，absFC>0.585且原始P<0.05，上下调保留，并报告BH。未按药物、通路、PPI或预测结果再改变疾病集合。

该集合与原恶性上皮专属Pool A813不是同一比较，交集576基因；原813及其统计未修改。PB仅为支持列，不成为新硬门槛；供者不足的Mast/Perivascular_SMC证据限制保留。Mixed/Unresolved只进入诊断敏感性。细胞类型可能有组成变化，表达方向不自动等于药物作用方向。

蛋白先按公开UniProt2026_03（2026-09-02）全部人类已审校记录在本地作Ensembl GeneId映射，再以官方MD5核验过的GENCODEv36蛋白翻译参考对未直接映射者补充。精确全序列匹配可恢复已审校记录；其余基因保留全部唯一GENCODE翻译，标明“无精确已审校匹配”，不按分数挑最佳异构体。

主输入中**{metrics['main_GENCODE_translation_only_genes']}基因**含仅GENCODE翻译参考，**{metrics['main_multiple_protein_genes']}基因**有多条输入蛋白；**{metrics['main_CNN_truncated_genes']}基因**至少一条蛋白超过1,000氨基酸而在CNN中截断，**{metrics['main_nonstandard_aminoacid_genes']}基因**含U或X等非20标准字符。未将U/X替换成其他氨基酸，按作者原编码处理并保留标记。CNN对其词表内字符编码，AAC没有为这些字符单独建立20氨基酸组成通道；不能将这些结果视为完整保留全部原始序列信息。

输入采用的是基因产品参考序列，未声称已确认这些异构体在相应患者中实际翻译。蛋白来源、序列哈希、长序列、非标准字符和多条目均在目录和总表可追溯。蛋白数据采用公开整库下载后本地连接，候选名单未向外部服务发送。

## 4. 模型、评分和共识

实际环境：Python3.12.14、DeepPurpose0.1.5、PyTorch2.5.1+cpu、RDKit2024.09.6；其余版本见logs/pip_freeze.txt。无自训、无微调、无模型权重修改。作者权重原件、配置及SHA256保留。

|编码|Kd/IC50模型|本轮角色|
|---|---|---|
|CNN药物＋CNN蛋白|各1个|各终点主要共识|
|Morgan药物＋CNN蛋白|各1个|各终点主要共识|
|Morgan药物＋AAC蛋白|各1个|各终点主要共识|
|Daylight药物＋AAC蛋白|各1个|各终点主要共识|
|MPNN药物＋CNN蛋白|各1个|逐对敏感性，不并入主要共识|

Kd权重输出按作者对数尺度解释为预测pKd，IC50权重为预测pIC50，值较高表示模型预测效力较强。**它们不是P值、FDR、结合概率、实验测量Kd/IC50或已校准的临床效应。** individual模型原文件的内部列名predicted_p_value指p尺度评分，不是统计P；最终交付总表统一改名predicted_log_affinity。未把两终点原分数直接相加或混合换算成概率。

先按每个基因的全部保留唯一蛋白序列取模型分数中位数，保留最小/最大值与序列数；在每个终点内部，分别把4个主要模型的基因分数转为同一候选范围的百分位，再取4模型百分位中位数排序。百分位为rank/n，同分采用平均秩；最终同共识值保留相同最小名次，展示再以Ensembl ID稳定排序。没有观察后取最大有利分数、设Top N或阈值删基因。

每个主要面板做留一编码敏感性；逐基因报告百分位IQR和最大排名变化。多编码可揭示表示差异，**不能消除共享BindingDB的训练偏置**。这10个权重并不构成10个独立数据库验证；Kd/IC50终点不同，但仍共享数据库来源。训练版本、全部训练靶点和拆分未由权重名称完整重建，本轮不作cold-drug/cold-protein或无训练重叠的保证。

本轮Morgan/Daylight原生特征在所检查的同连接图立体变体间完全相同，因此6个指纹模型没有为这些变体提供立体区分。CNN特征随SMILES字符变化，也不等于已证实模型能正确识别物理手性选择性。采用完整立体输入解决的是计算分子的明确性，不自动解决各编码的表示能力。

## 5. MPNN问题与技术检验的实际范围

原版MPNN在相同公共配对上随批大小/位置出现明显变化，batch1与batch4最大差异约2.619（Kd）和7.217（IC50）个对数单位；详见public_validation/MPNN_batch_dependency_diagnostic.json。使用作者未改权重、逐对推理后重复一致，但这**不证明模型原准确性已恢复**，所以只单列敏感性。

10模型都通过严格state_dict加载、有限值和公开作者2配对重复4行的技术检查。实际药物输入又在各模型首/中/末配对上复核单项与批量分数差异≤1e-4。没有将这些检查当作预测准确性。当前候选没有配套的实验结合标签，无法给出这1,881基因上的真实AUROC、误差、阳性率或“命中率”。

## 6. 排序展示

以下每个面板仅展示首10以便查看输出，不是正式靶点入选名单。完整1,881条及全部模型分歧都保留。

### Kd面板

{displayed_table(kd,'Kd')}

### IC50面板

{displayed_table(ic50,'IC50')}

留一编码结果见[敏感性表](results/Leave_one_encoding_out_stability.tsv)。Kd面板前20与留一后前20交集范围{int(stability.loc[stability.endpoint.eq('Kd'),'top20_overlap'].min())}–{int(stability.loc[stability.endpoint.eq('Kd'),'top20_overlap'].max())}/20，IC50为{int(stability.loc[stability.endpoint.eq('IC50'),'top20_overlap'].min())}–{int(stability.loc[stability.endpoint.eq('IC50'),'top20_overlap'].max())}/20。它描述排名对编码选择的敏感性，不作为额外疾病候选门槛。

## 7. 文件、图和验收

- [蛋白级全部评分](results/All_protein_model_predictions.tsv)：19,930行；[基因×模型评分](results/All_gene_model_scores.tsv)：各蛋白中位数、范围和数目。
- [主证据总表](results/Main1881_combined_disease_model_evidence.tsv)：疾病来源、方向、PB支持、映射限制、Kd/IC50各自共识及分歧。
- [1,957范围Kd敏感性](results/Diagnostic_extended1957_Kd_ranked.tsv)、[IC50敏感性](results/Diagnostic_extended1957_IC50_ranked.tsv)：诊断组额外候选明确标记，不改主范围。
- [模型排序相关图](figures/Model_rank_agreement.png)、[终点共识散点](figures/Kd_vs_IC50_consensus.png)、[首20编码一致性展示](figures/Top20_encoding_agreement.png)，均另有SVG向量版。图形只作结果描述，不因图形观感调整统计或排名。
- [独立验收](audit/Final_prediction_validation.tsv)：{validation['n_checks']}项一致性检查；原冻结结果、输入、模型哈希、配对完整性、实际结构、单项复核及独立重算秩均通过。

scripts/run_confirmed_multimodel_prediction.py可从已完成单模型输出断点恢复；summarize_predictions.py重建汇总/图；audit_prediction_results.py独立验收；本报告由write_prediction_report.py按实际输出生成。AAC特征缓存保留，避免重复编码。

## 8. 当前可作出的结论

本轮得到的是**固定CID209254在预先确认的疾病候选蛋白中的多模型探索排序**。它支持后续讨论优先验证对象，尚不证明直接结合、调控方向、功能必要性、乳腺癌机制、治疗效果或创新性。表达、模型分数、患者复核与实验因果证据应分开。

本轮未启动GO/KEGG筛选、PPI定靶、对接、MD或空间计算。后续候选取舍需结合模型分歧、蛋白映射/截断、细胞来源、文献和可行实验，由用户与助手共同讨论，不把任何一份前列名单自动当最终机制。
'''
    (ROOT/"DeepPurpose_结果报告.md").write_text(text,encoding="utf-8")
    print(json.dumps({"report":str(ROOT/"DeepPurpose_结果报告.md"),"checks_pass":validation["n_pass"],"cross_endpoint_rho":rho,"top20_overlap":overlap}),flush=True)

if __name__=="__main__":main()
