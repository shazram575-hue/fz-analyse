# 科研进度档案

更新：2026-10-26。范围仅工作流3：Spiropachysine A＋乳腺癌，目标Phytomedicine。MFSD10（SLC75A1，Q14728）是候选方向，药物直接结合、功能必需性及因果机制均未建立。本次整理没有新增科学假设检验。

## 当前结论与停止状态

**计算筛选、对接、表达摸底、逐患者程序发现、外部程序评估和探索性聚合均已完成；MFSD10表达程序关联仍没有可靠跨队列证据。已按用户指令停止自动优化这一表达程序关联路线。** 此判断不是证明MFSD10无作用，也不是取消候选身份。不得通过换阈值、换signature、换粒度或选患者追逐显著性。

新接手者先看[三层关联对照](../results/11_aggregation/简洁对照表.md)、[聚合结果摘要](../results/11_aggregation/结果摘要_探索性聚合再分析.md)、[停止标记](../results/11_aggregation/audit/STOP_THIS_ROUTE.json)。复制件源路径和SHA256统一在[PATH_MAP](PATH_MAP.tsv)。完整旧目录仍在`outputs`；`analysis/results`只是Git审阅副本。

## 跨任务恢复与共享状态

上下文机制更新：2026-10-27。每个新任务从根[AGENTS](../AGENTS.md)和本页当前结论/导航恢复，再读相关方向记录、代码与结果；不假定知道其他聊天。实质进展写独立任务日志，并在本页补入口或状态变化，详细历史无需重复堆入AGENTS。

**2026-10-27用户已授权首次GitHub安全同步。** 仓库[shazram575-hue/fz-analyse](https://github.com/shazram575-hue/fz-analyse)归属已核实，实际提交、推送及同步状态以[本轮日志](TASK_LOGS/2026-10-27_github_sync.md)和当前Git为准；2026-10-26的`LOCAL_ONLY`/非仓库记录是历史快照。分别记录计算状态、科学结论和同步状态；只有经过审核、提交并合并到约定集成分支的变更才是`MERGED_REVIEWED`。上传本身不等于科研审核或外置数据恢复，未提交/未合并部分仍需标明。

具体流程与任务索引见[TASK_LOGS](TASK_LOGS/README.md)，新任务用[启动指令](CLOUD_TASK_TEMPLATE.md)和[交接模板](TASK_LOGS/TEMPLATE.md)。本次文档任务见[交接记录](TASK_LOGS/2026-10-26_context_handoff.md)，未新增科学计算、未推送或合并。以下“完成”均指已有本地计算交付，非假设成立或GitHub同步完成。

2026-10-27同步完成：用户明确允许公开上传公共分析结果，取消先改私有要求；代码和小型证据已同步到GitHub main，导入提交`9342813d8f45cee7228e5a1fa58d33ad76b6c046`，831文件树与本地审查版本完全一致，已独立fetch和读取远程必需文档。实际传输经已登录GitHub连接器完成；本机Git CLI push缺凭据的失败保留在[同步报告](GITHUB_SYNC_REPORT.md)。没有新增科学检验或改变下文结论，大型数据仍未上传/恢复到Cloud；未额外进行科研PR审核合并。

## 全项目模块导航

| 分析模块 | 定向接手记录 | 结果入口 | 代码 | 科学/工程状态 |
|---|---|---|---|---|
| 01 疾病侧Pool813 | [疾病/图谱](TASK_LOGS/01_DISEASE_AND_ATLAS.md) | [索引](../results/01_disease/README_结果索引.md) | [脚本](../analysis/01_disease/scripts) | 计算完成，关联筛选 |
| 02 单细胞图谱及1881输入 | [疾病/图谱](TASK_LOGS/01_DISEASE_AND_ATLAS.md) | [索引](../results/02_cell_atlas/README_结果索引.md) | [脚本](../analysis/02_cell_atlas/scripts) | 计算完成，少供者组仅描述 |
| 03 DeepPurpose池内预测 | [预测/对接](TASK_LOGS/02_PREDICTION_AND_DOCKING.md) | [报告](../results/03_deeppurpose/DeepPurpose_结果报告.md) | [脚本](../analysis/03_deeppurpose/scripts) | 计算完成，非实测结合 |
| 04 全蛋白预测及疾病交集 | [预测/对接](TASK_LOGS/02_PREDICTION_AND_DOCKING.md) | [报告](../results/04_proteome/独立预测与疾病交集_结果报告.md) | [脚本](../analysis/04_proteome/scripts) | 计算完成，偏置限制保留 |
| 05 查新 | [预测/对接](TASK_LOGS/02_PREDICTION_AND_DOCKING.md) | [入口](../results/05_pubmed/README_当前入口.md) | [脚本](../analysis/05_pubmed/scripts) | 历史检索完成，非永久创新性证明 |
| 06 全47反向对接 | [预测/对接](TASK_LOGS/02_PREDICTION_AND_DOCKING.md) | [报告](../results/06_docking47/结果报告_47靶点反向对接.md) | [脚本](../analysis/06_docking47/scripts) | 计算完成，无MD/实测亲和力 |
| 07 候选选择 | [预测/对接](TASK_LOGS/02_PREDICTION_AND_DOCKING.md) | [证据](../results/07_candidate/现有证据支持的候选优先级.md) | 无独立脚本，依据前述结果汇总 | MFSD10候选，非确证靶点 |
| 08 MFSD10表达摸底 | [表达/程序](TASK_LOGS/03_MFSD10_PROGRAMS.md) | [报告](../results/08_expression/MFSD10_表达分布摸底_结果报告.md) | [脚本](../analysis/08_expression/scripts) | 完成，稀疏/深度限制 |
| 09 程序可行性 | [表达/程序](TASK_LOGS/03_MFSD10_PROGRAMS.md) | [患者表](../results/09_feasibility/患者可用性表.md) | [脚本](../analysis/09_feasibility/scripts) | 历史前置检查完成，正式运行见10 |
| 10 程序发现与外部评估 | [表达/程序](TASK_LOGS/03_MFSD10_PROGRAMS.md) | [报告](../results/10_programs/结果报告_程序发现与独立验证.md) | [脚本](../analysis/10_programs/scripts) | 完成，程序协同性有支持，MFSD10关联不足 |
| 11 探索性聚合 | [表达/程序](TASK_LOGS/03_MFSD10_PROGRAMS.md) | [三层对照](../results/11_aggregation/简洁对照表.md) | [脚本](../analysis/11_aggregation/scripts) | 完成，停止自动优化关联路线 |
| 工程整理与云迁移 | [工程/云](TASK_LOGS/04_ENGINEERING_AND_CLOUD.md) | [审计](PROJECT_AUDIT.md)、[迁移](CLOUD_MIGRATION.md) | [工程工具](../scripts/project) | 审阅视图已有，Linux完整复现未验证 |

上述模块原接手记录的`LOCAL_ONLY`为当时快照；最新同步确认见本轮Git日志。大型数据负载仍未上传，代码合并后也须独立核实恢复位置及校验。详细方法、阴性结果与证据继续保留在下文和原报告中。

## 数据与来源

| 数据 | 实际范围 | 输入位置与来源证据 | 边界 |
|---|---|---|---|
| TCGA-BRCA STAR counts | 1208文件，1095肿瘤＋113癌旁；1095患者、113配对 | `原始数据/芯片/files`；[官方UUID/MD5下载清单](../data/TCGA_DOWNLOAD_MANIFEST.tsv) | 目录叫“芯片”，实际是RNA-seq，不得按microarray处理 |
| GSE161529 | 69文库、428024原始细胞；作者QC后340979；主47文库236309细胞 | `原始数据/单细胞/GSE161529_RAW.tar`及features/SOFT | 主34原发Total＋13健康Total；11正常Epi、4癌前、7LN分开 |
| GSE176078 | 100064细胞×29733作者基因，24489作者Cancer候选 | `原始数据/单细胞/GSE176078`；官方GEO及本地作者补表/方法 | 26临床Case键；20病例在最终矩阵有Cancer候选，主要16Naïve、4Treated另列 |
| 分子及蛋白 | CID209254，九中心完整结构；UniProt release2026_03、GENCODE v36 | `outputs/DeepPurpose_预测_2026-10-24`与`outputs/身份深核与独立交集_2026-10-24` | 计算分子与实物质量控制不是一回事；蛋白异构体及截断限制保留 |
| 结构/权重/文献 | PDB/AlphaFold、DeepPurpose作者预训练权重、PubMed原始查询 | `outputs/docking_audit_2026-10-24`及各阶段sources/models | 未取得或未核读的全文不能算已读；文献/权重重分发授权不默认成立 |

体积、格式、精确校验及恢复见[DATA_MANIFEST](DATA_MANIFEST.md)。本轮1214个原始文件已流式SHA256/MD5；1208个TCGA官方MD5全部一致。

## 已完成的流程、方法和证据

### 1. 疾病侧初筛：Pool A813

Bulk用edgeR TMM＋limma-voom，Tumor−Normal，`abs(log2FC)>0.585 & P.Value<0.05`，报告BH；双向保留。113真实配对校准block相关，另有配对敏感性。

单细胞主比较为确认恶性上皮vs健康来源正常上皮，`pct.1>0.15 & avg_log2FC>0.5 & Bonferroni<0.05`，校正覆盖33487规范特征。只保留protein_coding进入原Pool813；患者PB另作支持，不强制成为候选池门槛。≥20细胞是此阶段PB描述/可估计门槛，不是程序推断充分条件。

CNV唯一主分类为17c `Tumor_discovery_cells_CNV_coherent_final_calls.tsv`：基于每患者最多1000上皮的发现集，9486高置信恶性、19627已评分未定；未向其余细胞投射恶性标签。HMM=FALSE、denoise=TRUE、sdAmplifier=1.5，RMS>健康q99且同供者基因组模式一致。CNV阴性/未定不等于正常。

证据：[原Pool表](../results/01_disease/results/Pool_A_candidate_genes.tsv)、[冻结方案原件](../第一阶段_分析方案与验收.md)、[第一阶段报告](../第一阶段_结果报告.md)。代码：[Bulk](../analysis/01_disease/scripts/05_bulk_limma_voom.R)、[CNV17c](../analysis/01_disease/scripts/17c_coherent_CNV_final_calls.R)、[标记/PB](../analysis/01_disease/scripts/19_singlecell_markers_and_pseudobulk.R)。

### 2. 单细胞基础补全与疾病1881输入

全主图谱236309细胞；28558供者×类型平衡子集训练PCA，再投影全量Harmony/UMAP/聚类。tSNE仅子集，不伪造全量位置。6生物学大类＋2诊断组，不宣称精细免疫分型。

DE仅同类型Tumor-origin−Healthy，`abs(avg_log2FC)>0.5`、`max(pct.tumor,pct.healthy)>0.15`、全33487特征Bonferroni<0.05，双向保留；PB以供者聚合，主条件每组≥20细胞、每条件≥3供者，另年龄/性别敏感性。Mast与Perivascular_SMC一侧仅1供者，不能作群体推断。

后续药物疾病输入为六大类细胞级双向DE∩Bulk∩protein_coding，共1881；Mixed/Unresolved额外76单列。原813不改，不能把1881改称813。证据：[补全报告](../单细胞基础补全_结果报告.md)、[分类型摘要](../results/02_cell_atlas/results/Within_celltype_DE_summary.tsv)、[1881输入](../results/03_deeppurpose/main_input/Disease_matching_candidates_1881.tsv)。

### 3. DeepPurpose及独立蛋白背景

先完成1993化合物−蛋白配对×10模型＝19930原分数。Kd/IC50各4非MPNN编码分开共识；2个MPNN已复现批次依赖，逐对敏感性，不并入主共识。无自训/微调，不能称预测准确性已验证。

随后独立人类已审校背景20431蛋白×8模型＝163448分数；先排药物侧百分位，再与疾病1881相交。排名边界并列全保留；双终点共有∩疾病的1%/5%/10%结果为2/17/47基因，用户采用10%。这里“独立”指背景不由疾病候选预筛，不是独立训练数据库验证；模型共享BindingDB偏置。

实际分子为CID209254、InChIKey `SLGWGPQWJRVPAD-PZOAWPBASA-N`；后续原始Figure1九中心重建已支持A对应完整结构，覆盖早期仅命名推定的历史状态。输入采用等价92字符规范异构SMILES以避CNN100字符截断；蛋白CNN1000残基限制及异构体来源保留。预测pKd/pIC50不是统计P值、结合概率或实验测定值。

证据：[池内报告](../results/03_deeppurpose/DeepPurpose_结果报告.md)、[独立预测报告](../results/04_proteome/独立预测与疾病交集_结果报告.md)、[阈值总表](../results/04_proteome/results/Threshold_comparison_1_5_10percent.tsv)。

### 4. 查新、全47对接及MFSD10选择

2026-10-24扩充PubMed检索是带日期的操作性记录：严格返回0的5项含MFSD10/SLC75A1；不等于永久文献空白或机制创新已证明。SAMD8机构误命中消歧另列，MT-ND4L和CASIMO1相关真实文献保留。

全47靶点×两药物状态×3种子＝282主对接，9参照复对接；Vina1.2.7、exhaustiveness64，25靶点实验结构优先、22个AF。15项技术验收通过，但距离检查通过268/282，失败项没有清除；PPARG参照最高姿势2/3、PNP3/3、SLC18B1为0/3，不能推广参照有效性。MFSD10两状态评分约−11.683/−11.934是计算分数，不是实测亲和力。

用户选择MFSD10为研究方向；不是模型自动证明的直接靶点。没有本线正式MD结果可核验，不视作已完成。

证据：[查新](../results/05_pubmed/PubMed零命中标准_复核结论.md)、[47报告](../results/06_docking47/结果报告_47靶点反向对接.md)、[整合证据](../results/06_docking47/results/47_integrated_docking_disease_evidence.tsv)。

### 5. MFSD10表达摸底、程序发现与外部评估

MFSD10唯一映射ENSG00000109736。患者/组织/制备分层统计原始RNA counts，未用integrated表达做统计，未按MFSD10高低或阳阴分组。恶性9486细胞中60.3%零值，检出者55.6%仅1UMI；32患者中27名≥20恶性细胞，但这不是cNMF充分条件。健康背景、CNV未定19627、未评分肿瘤上皮93898分别保留。LN/癌前无粗注释，只作混合样本描述。

程序发现用已有8722确认恶性singlet中13个A/B候选患者（6071细胞）；完整未抽样图谱没有额外同标准已判恶性可直接补入。A/B只是工作分层；后续做了每患者K4/6/8、各100重复、2000HVG，共3900主分解及19场景380次敏感性。MFSD10物理移除出输入、归一化分母及末端全基因重拟合。5个固定50基因程序支持患者11/5/5/3/3；不按MFSD10相关性挑程序。

外部GSE176078未与发现集合并/整合：26病例中20在最终矩阵有作者Cancer细胞；主要16Naïve（21521细胞）、4Treated（2968）单列，另6无该标签不等于肿瘤无恶性。作者QC均满足，不追加剔除；doublet状态未知。CID3963 metadata为TNBC但补表IHC为ER+，保留冲突、按显式IHC调整，该例治疗史另列。未推PAM50。

程序/评分/分析公式/发现效应方向先冻结。固定基因覆盖94%–100%，5程序的外部基因残差协同性有患者级支持，不能仅据评分高称复现。发现13患者的MFSD10患者内/间10检验联合BH均未显著，最小q0.207；确认性关联清单为空。外部MP05患者内rho0.057、q0.00290仅探索性，不能倒改发现。

原模型：每患者partial rho控制目标剔除后UMI、基因数和线粒体，再患者等权Fisher-z、患者符号翻转P及患者bootstrap CI；患者间PB模型调整IHC/深度/线粒体，n13/16、df7/10。并非把细胞数作患者数。部分主重复到达迭代上限，局部延长仍有分量变化；0043 UMI减半的谱余弦/usage相关约0.729/0.662，不能称所有分量均稳定。

证据：[表达摸底](../results/08_expression/MFSD10_表达分布摸底_结果报告.md)、[程序报告](../results/10_programs/结果报告_程序发现与独立验证.md)、[冻结基因](../results/10_programs/frozen/Program_genes_v1.tsv)、[发现全表](../results/10_programs/discovery/MFSD10_patient_level_associations.tsv)、[外部全表](../results/10_programs/validation/Discovery_external_comparison_all.tsv)。

### 6. 探索性聚合再分析及停止

两个队列已经查看，故全部新增结果统一探索性。原20项模型直接复用；患者间PB标准化及程序均分与旧表逐项相符，未再跑NMF或筛signature。

每患者×样本约20细胞非重叠状态聚合，构邻域排除MFSD10，只用RNA counts/PCA；冻结评分保持，另有原细胞均分和随机同容量等6规格。发现全部8722细胞→425单元，验证24489→1215，UMI守恒。主要效应13/15患者；验证4465仅6单元，按既定规则不可估，同15患者旧效应对照保留。不是按目标/P值剔患者。

| 聚合患者内效应 | 发现13患者 rho [CI]；q | 验证15患者 rho [CI]；q |
|---|---|---|
| MP02 | 0.088 [−0.030,0.201]；0.379 | 0.285 [0.174,0.421]；0.00137 |
| MP05 | 0.136 [−0.037,0.285]；0.370 | 0.282 [0.164,0.405]；0.00137 |

两队列10项聚合主检验共同BH；旧FDR保留旧家族，不能只凭q数值变化认定改善。全部5发现聚合结果均未过校正；患者间完整和精简模型仍未支持。零值比例患者中位约53%/48%降至0%，但不等于每个细胞表达或恢复单细胞真值。MP01/MP03近零/不稳；MP04对评分/调整敏感；MP02/MP05外部增强而发现不精确。没有强贴“技术因素解释较多”标签。

**停止优化表达程序关联路线**，保留30项主对照、60规格、30患者间敏感性和所有阴性项。32技术核验、5473份前轮文件保护及文档收尾标记可查。研究展望需新授权，不能在云端继续悄悄调参。

## 未解决问题与允许的下一步

1. 云端工程：先通过证据核验、补外部数据位置与Linux环境；仅修复路径/资源问题，不重估已完成结果。见[CLOUD_MIGRATION](CLOUD_MIGRATION.md)。
2. 完整复现：Linux依赖、R混合包库、Windows对接二进制、模型权重来源/许可证及大型对象持久存储尚未全部验证。文献查新是日期快照，不能延用为绝对创新性保证。
3. 科研：MFSD10直接结合、功能依赖、作用方向及细胞背景仍待证实。当前结果不指定任何已验证的下游机制。
4. 湿实验仅展望：核实化合物身份/纯度；在实验可行的乳腺癌或相关细胞背景验证MFSD10表达、扰动与药效依赖及救援；独立评估结合或靶点占用。所需样本、细胞系和检测能力未核实，不是本轮执行授权。

后续可讨论的生信工作需先由用户选择新问题：基于现有47证据整理结构/序列/疾病证据缺口、核对公开功能资料及实验可行性；或者先做样本量与独立验证设计。它们尚未启动，不默认继续已停止的表达程序优化，也不预设具体通路机制。

历史完整记录见[工作流3交接](../工作流3_交接.md)。2026-10-25已授权清除的阶段不能恢复或引用；本轮之后新分析是独立新授权结果，不能与已清除内容混同。
