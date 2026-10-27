# 疾病侧筛选与单细胞图谱：接手记录

记录日期：2026-10-26；类型：既有结果索引迁入。本轮没有重跑分析。同步快照：`LOCAL_ONLY`，本目录尚无可核实Git/审核合并基线；后续确认按[同步规则](README.md)追加。

## 实际状态与证据

- 01疾病侧初筛、02单细胞基础补全已有完成交付；不是MFSD10功能或药物结合验证。[科研进度第1–2节](../RESEARCH_PROGRESS.md)保留完整方法与结果说明。
- 原Pool A813冻结：[结果索引](../../results/01_disease/README_结果索引.md)、[候选表](../../results/01_disease/results/Pool_A_candidate_genes.tsv)、[方案原件](../../第一阶段_分析方案与验收.md)。Bulk为TCGA STAR RNA-seq，以edgeR TMM＋limma-voom分析，`abs(log2FC)>0.585 & P.Value<0.05`，报告BH；上下调均保留。单细胞确认恶性vs健康上皮用`pct.1>0.15 & avg_log2FC>0.5 & 全特征Bonferroni<0.05`，protein_coding进入813。
- [基础图谱索引](../../results/02_cell_atlas/README_结果索引.md)、[方法参数](../../results/02_cell_atlas/方法与参数说明.md)；六大类同类型肿瘤来源vs健康来源DE形成后续疾病输入1881，Mixed/Unresolved额外76仅诊断。1881不是改写813，不做类型vs其余替代主比较。
- 代码入口：[01](../../analysis/01_disease/scripts)、[02](../../analysis/02_cell_atlas/scripts)；原件映射及SHA在[PATH_MAP](../PATH_MAP.tsv)，完成标记与依赖在[stages.json](../../config/stages.json)。副本不是已移植Linux入口。

## 必须沿用的解释边界

- GSE161529唯一CNV主分类17c：`Tumor_discovery_cells_CNV_coherent_final_calls.tsv`；未定/未评分不纳入确认恶性，也不归非恶性。不能把全部肿瘤来源上皮视为恶性。
- 患者映射优先作者临床/Case信息：0029多区域、0114不同分型病灶仍各为同患者；健康0064与癌症0064不同；0230.16/.17保留原供者区分。不得自行猜PAM50或未知治疗。
- 7个LN和4个癌前文库无粗注释，只作混合细胞表达描述；不按组织名称补细胞类型。
- 患者为独立重复；≥20细胞仅是既有PB描述/可估计门槛，不是NMF程序推断充分标准。Mast与Perivascular_SMC一侧仅1供者，不作患者群体推断。
- 原始数据、对象与所有统计结果保持原状；2026-10-25被授权清除的阶段不恢复/引用。后来明确新授权产物独立保留。

## 恢复与待办

本轮仅查阅索引/已有科研档案，没有启动本模块任务或验证全套环境。若用户后续授权工程复现，先按[数据清单](../DATA_MANIFEST.md)确认对象及原件，再核对对应脚本依赖；不得借移植改统计语义。未解决科研问题见统一进度，不能把表达差异称为因果作用。
