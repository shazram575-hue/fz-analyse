options(stringsAsFactors=FALSE)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","单细胞基础补全_2026-10-23")
stopifnot(all(file.exists(file.path(out,"cache",c("Global_embedding.complete","Within_celltype_DE.complete","Figures.complete","Independent_validation.complete")))))
ds<-fread(file.path(out,"results","Within_celltype_DE_summary.tsv"))
metrics<-fread(file.path(out,"audit","Completion_metrics.tsv"))
getn<-function(n)metrics[metric==n,value]
fmt<-function(x)ifelse(is.na(x),"不估计",format(x,big.mark=",",scientific=FALSE,trim=TRUE))
lab<-c(Endothelial="内皮",Epithelial="上皮（组织来源分组）",Fibroblast="成纤维",Immune="免疫大类",Mast="肥大细胞",Mixed_epithelial_immune="上皮/免疫混合〔诊断组〕",Perivascular_SMC="周细胞/平滑肌",Unresolved="未定〔诊断组〕")
tab<-vapply(seq_len(nrow(ds)),function(i){s<-ds[i];paste0("|",lab[s$cell_type],"|",s$n_healthy_donors_min20," / ",s$n_tumor_donors_min20,
 "|",fmt(s$PB_up)," / ",fmt(s$PB_down),"|",fmt(s$celllevel_up)," / ",fmt(s$celllevel_down),"|")},character(1))
report<-c("# 单细胞基础补全与同类型疾病差异结果", "", "项目：Spiropachysine A＋乳腺癌。主体计算日期：2026-10-23；收尾验收日期：2026-10-24。", "",
 "## 本轮范围与结果保护", "",
 "按用户确认，本轮只比较同一现有注释组内的肿瘤来源vs健康来源，没有进行细胞类型vs其余类型的DE。现有第一阶段813基因Pool A、恶性上皮专属比较和原始数据保留；新增内容全部进入outputs/单细胞基础补全_2026-10-23，原结果另有冻结摘要及Pool快照。", "",
 sprintf("主图谱沿用47份Total文库、45个供者键、%s个作者QC后细胞。当前注释为6个生物学大类及Mixed/Unresolved两个诊断组，没有把免疫大类自动命名为精细T/B/髓系亚群。正常Epi富集、癌前和淋巴结未混入主图谱。",fmt(getn("Main_cells"))), "",
 "## 补齐的基础分析", "",
 "- 逐文库QC分布、当前矩阵的QC保留数、疑似双细胞比例；作者QC和标记策略未更改。",
 sprintf("- 跨文库共识HVG、PCA载荷与碎石图。按供者×现有类型最多200细胞预先均匀抽样，共%s细胞训练PCA，再对全部主细胞显式投影。不是对全部细胞直接拟合PCA，也不将该抽样称为Seurat LeverageScore。",fmt(getn("PCA_training_landmarks"))),
 "- 未整合与Harmony整合后的全细胞UMAP对照，以及按现有注释、组织来源、全局cluster、文库、线粒体比例与细胞周期着色。",
 "- 全细胞邻居图与Louvain 0.2/0.4/0.8/1.2多分辨率聚类、cluster注释构成及分辨率转换计数。原逐细胞注释保持，cluster多数标签不自动替换身份。",
 "- 平衡子集tSNE；该图只含训练/描述子集，没有将tSNE位置外推到其余细胞。",
 "- 已知标志基因点图、逐供者细胞组成、细胞周期评分与比例，以及同类型疾病DE火山图/供者热图。已知标志点图是注释检查，不是类型vs其余的DE。", "",
 sprintf("全局cluster数量：resolution0.2=%s，0.4=%s，0.8=%s，1.2=%s。0.4为基线展示，不宣称是唯一最优分辨率。",fmt(getn("Global_clusters_resolution_0.2")),fmt(getn("Global_clusters_resolution_0.4")),fmt(getn("Global_clusters_resolution_0.8")),fmt(getn("Global_clusters_resolution_1.2"))), "",
 "## 同类型肿瘤来源vs健康来源DE", "",
 "患者pseudobulk使用全部主细胞的规范原始UMI；同供者重复文库合并，每供者×类型至少20细胞，每条件至少3位供者才估计患者疾病差异。主模型condition，年龄/性别调整另作敏感性；edgeR TMM＋limma-voom，BH在每个类型/模型内部控制。DE显示门槛为BH<0.05且abs(log2FC)>0.5，完整被检结果一并保存。", "",
 "细胞级结果使用预先选定、每供者×类型最多200细胞的集合；Seurat RNA LogNormalize10,000＋双向Wilcoxon，min.pct=0.15、abs(avg_log2FC)>0.5、max(pct.1,pct.2)>0.15、全33,487特征Bonferroni<0.05。细胞级P仅作描述，不把细胞数当患者数。未做蛋白编码或通路成员筛选。", "",
 "|现有注释组|PB合格供者：健康 / 肿瘤|PB上调 / 下调|细胞级上调 / 下调|", "|---|---:|---:|---:|",tab,"",
 "肥大细胞健康侧仅1位合格供者，周细胞/平滑肌肿瘤侧仅1位，因此没有生成可推广的pseudobulk P值；探索性细胞级表保留patient_evidence_limited标记。这是当前队列的生物重复限制，不能用增加同一供者细胞或改阈值解决。混合/未定组的差异是诊断信息，不作为已建立的细胞类型机制。", "",
 "本轮Epithelial按组织来源比较，包含肿瘤来源但CNV未定或未评分上皮；它与原9,486恶性vs12,264正常的比较对象不同。原1,838正向markers及813 Pool不替换、不重新取交集。免疫大类内的组成差异也可能影响表达，不能直接声称细胞内在改变。", "",
 "## 文件入口", "",
 "- [完整输出索引](outputs/单细胞基础补全_2026-10-23/README_结果索引.md)",
 "- [各注释组DE汇总](outputs/单细胞基础补全_2026-10-23/results/Within_celltype_DE_summary.tsv)",
 "- [逐基因并列证据表](outputs/单细胞基础补全_2026-10-23/results/Within_type_combined_gene_evidence.tsv)：8组×33,487特征，未检验项为NA",
 "- [患者模型调整与细胞级一致性](outputs/单细胞基础补全_2026-10-23/results/Within_type_DE_robustness_summary.tsv)",
 "- [患者pseudobulk全部统计](outputs/单细胞基础补全_2026-10-23/results/All_celltypes_pseudobulk_all_tested.tsv)",
 "- [年龄/性别调整全部统计](outputs/单细胞基础补全_2026-10-23/results/All_celltypes_pseudobulk_age_gender_all_tested.tsv)",
 "- [细胞级全部通过表达预筛的统计](outputs/单细胞基础补全_2026-10-23/results/All_celltypes_celllevel_all_prefiltered.tsv)",
 "- [全细胞元数据及cluster](outputs/单细胞基础补全_2026-10-23/metadata/All_main_metadata_with_global_clusters.tsv)",
 "- [独立验收](outputs/单细胞基础补全_2026-10-23/audit/Independent_validation.tsv)",
 "- [参数、参考与实现说明](outputs/单细胞基础补全_2026-10-23/方法与参数说明.md)","",
 "全局可恢复对象为cache/All_main_atlas_embedding_bundle.rds，含元数据、PCA/Harmony/UMAP/tSNE及规范counts分片路径。它是降维bundle，不冒充已合并全基因Seurat assay。原47份全基因Seurat对象仍在第一阶段缓存，新增规范全基因counts按文库分片，适应本机内存。", "",
 "## 收尾验收", "",
 "2026-10-24完成22组PDF/PNG。20项独立数据/模型检查通过；图形检查覆盖全部PNG概览及关键图完整尺寸，记录见输出目录audit/Figure_visual_QA.md。原结果保留校验、PDF完整性和交付文件SHA256分别见audit/Frozen_stage1_after.tsv、Published_stage1_output_preservation.tsv、PDF_completeness.tsv及Deliverable_SHA256_manifest.tsv。", "",
 "## 解释边界", "",
 "Harmony仅用于图谱几何与聚类，DE使用原始RNA/counts。供者与疾病来源存在嵌套，不能仅凭UMAP混合判定批次已完全消除，也不能把整合图当机制证明。细胞组成比例受取材、解离和QC影响；本轮未进行患者群体的细胞丰度差异检验。", "",
 "本轮没有按GO/KEGG筛基因，没有更新药物靶点或创新性结论。所有疾病DE均为同数据集内证据。")
writeLines(enc2utf8(report),file.path(root,"单细胞基础补全_结果报告.md"),useBytes=TRUE)
index<-c("# 单细胞基础补全输出索引", "", "状态：完成。用户口径：只做同类型肿瘤来源vs健康来源，不做类型vs其余DE。", "",
 "报告：[单细胞基础补全结果](../../单细胞基础补全_结果报告.md)。原Pool813位于第一阶段，未替换。", "",
 "## 结果表", "",
 "- results/Within_celltype_DE_summary.tsv：8现有组的样本量、可估计状态和上下调数。",
 "- results/within_celltype/：各组细胞级预筛全表/显著表、可估计组pseudobulk全表/显著表及年龄性别敏感性。",
 "- Within_type_combined_gene_evidence.tsv：267,896条基因×类型并列证据；Within_type_DE_robustness_summary.tsv：调整前后数量及方向一致性，无新筛选门槛。",
 "- All_celltypes_pseudobulk_all_tested.tsv和All_celltypes_pseudobulk_age_gender_all_tested.tsv：各类型合并统计。",
 "- All_celltypes_celllevel_all_prefiltered.tsv：双向Seurat描述；不是未经表达预筛的全基因双向检验表。",
 "- Library_QC_summary.tsv、Celltype_condition_counts.tsv、Donor_celltype_composition.tsv、Cellcycle_by_celltype_condition.tsv：基础统计。",
 "- PCA_variance_explained.tsv、PCA_top_absolute_loadings.tsv、Global_cluster_annotation_composition.tsv、Global_resolution_transition_counts.tsv：降维与聚类解释。", "",
 "## 元数据与缓存", "",
 "metadata/All_main_metadata_with_global_clusters.tsv含236,309细胞及全部分辨率；All_main_UMAP_pca.tsv和All_main_UMAP_harmony.tsv分别为未整合/整合全量坐标；Balanced_landmark_tSNE_Harmony.tsv仅28,558子集坐标。",
 "",
 "cache/All_main_atlas_embedding_bundle.rds为可恢复降维bundle；All_main_Harmony_SNN.rds为全局图；Balanced_landmark_PCA_model.rds包含均值、数值标准差、载荷和训练cell IDs；canonical_by_library为33,487全基因计数分片。原逐库完整Seurat不改写。", "",
 "## 图", "",
 "Foundation_01–02：逐文库QC/保留；03–04：PCA/Harmony/碎石；05–08：全局UMAP、整合比较、子集tSNE、按文库拆图；09–10：周期；11：已知marker点图；12：供者组成；13–14：聚类分辨率/注释；15–16：疾病DE汇总/火山；17_*：各可估计组供者热图。每图PDF及PNG。", "",
 "## 审计与恢复", "",
 "audit/Frozen_stage1_before.tsv及Frozen_stage1_after.tsv验证7份原主结果SHA256；Independent_validation.tsv为20项数据/模型口径检查；Completion_metrics.tsv为实际数量；Frozen_stage1文件校验需在本轮最后执行。",
 "",
 "scripts顺序01准备→03疾病DE/03b并列证据→02合并图谱→04出图→05验收→06报告。各模块有cache/*.complete及logs。原件不删；本轮只新增目录。稀有组/诊断组限制见报告。")
writeLines(enc2utf8(index),file.path(out,"README_结果索引.md"),useBytes=TRUE)
methods<-c("# 方法、参数与来源", "", "## 范围", "",
 "原47主文库及旧逐细胞大类注释保留。图谱不自动精细命名免疫亚群。新增DE只限同类型Tumor-origin−Healthy；原Pool和恶性专属对比冻结。", "",
 "## 降维与聚类", "",
 "按各文库原Seurat HVG的跨文库出现频率及中位排名选共识2,000 HVG；与疾病DE和候选名单独立。每供者×现有类型最多200细胞均匀随机，种子20261023，28,558 landmark用于PCA训练。所有细胞的基因标准化均除以全基因总UMI后乘10,000并log1p；不是按HVG总计数重新定库大小。",
 "",
 "显式中心/标准差＋irlba SVD，50PC；以保存的数值中心/尺度/载荷对全236,309细胞投影。前30PC用于下游。该代表子集训练/全量投影参考[Seurat可扩展分析思路](https://satijalab.org/seurat/articles/seurat5_sketch_analysis)，实际采用供者×类型分层均匀抽样，未运行LeverageScore。",
 "",
 "Harmony1.2.4：vars_use=gsm、theta=1、lambda=1、sigma=0.1、nclust=100、max_iter=10、ncores=2，保存未校正PCA。依据[Korsunsky等，Nature Methods2019](https://doi.org/10.1038/s41592-019-0619-0)及[Seurat整合教程](https://satijalab.org/seurat/articles/seurat5_integration)。仅低维坐标校正，原始表达不改。",
 "",
 "邻居图：Seurat FindNeighbors，k20、Annoy50棵树、Euclidean/SNN；Louvain algorithm1，res0.2/0.4/0.8/1.2，种子20261023。UMAP：原PCA与Harmony分别全量运行，cosine、neighbors30、min_dist0.3、spread1、200epochs、计算线程2、SGD线程1。tSNE仅28,558landmark，在Harmony30维上perplexity30、theta0.5、1000iterations、2线程，不向其余细胞投影。",
 "",
 "细胞周期：Seurat cc.genes.updated.2019映射Ensembl，CellCycleScoring逐库描述，无回归或剔除。基础QC沿用作者，不看图反复调阈值。", "",
 "## 疾病差异", "",
 "全量原始UMI按供者×类型聚合，合并重复文库；至少20细胞的组且每条件至少3供者，condition主模型及年龄/性别可识别模型，edgeR TMM/filterByExpr＋limma-voom。BH按每类型/模型全部被检基因；显著展示abs(log2FC)>0.5、BH<0.05。供者不足不捏造患者推断。",
 "",
 "Seurat5.4.0细胞级双向Wilcoxon用平衡抽样，RNA LogNormalize，min.pct0.15、logFC门槛0.5；最终absFC>0.5、max(pct.1,pct.2)>0.15、p_val_adj<0.05。校正分母33,487全特征，非按Bulk预选。[官方DE说明](https://satijalab.org/seurat/reference/findmarkers)。患者证据设计参考[Squair等2021](https://doi.org/10.1038/s41467-021-25960-2)。所有细胞级P为描述，未把细胞当独立患者。", "",
 "## 图形与复现", "",
 "2026-10-24收尾时，47文库分面仅调整点层渲染分辨率以加快出图；修复cluster0的factor显示编码，并改善火山图点的可读性。坐标、细胞分配、注释、表达和统计均未因图形修复更改。绘图支持复用完整已存在的图对；需要全部重画可传force参数。", "",
 "沿用Tableau系列固定类别色值：上皮#4E79A7、成纤维#F28E2B、免疫#59A14F、内皮#76B7B2、周细胞#B07AA1、肥大#EDC948、混合#E15759、未定#BAB0AC。并非复刻某篇Nature图的原配色。全细胞散点仅点层栅格化，其余PDF文字/轴保留向量；并存200dpi PNG便于检查。",
 "",
 "实现修复：历史逐库元数据列数不完全一致，按列名合并补NA且不改生物标签；data.table中位数显式转数值。当前prcomp_irlba封装在scale.=TRUE时保留逻辑标志而不是拟合标准差，独立投影校验及时检出；改用显式数值均值/标准差的irlba，保存并验算投影，未放宽错误模型。完整版本/日志在logs，原输入与第一阶段结果不改。")
writeLines(enc2utf8(methods),file.path(out,"方法与参数说明.md"),useBytes=TRUE)
cat("REPORTS WRITTEN",format(Sys.time()),"\n")
