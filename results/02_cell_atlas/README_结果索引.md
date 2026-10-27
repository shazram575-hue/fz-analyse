# 单细胞基础补全输出索引

状态：完成。用户口径：只做同类型肿瘤来源vs健康来源，不做类型vs其余DE。

报告：[单细胞基础补全结果](../../单细胞基础补全_结果报告.md)。原Pool813位于第一阶段，未替换。

## 结果表

- results/Within_celltype_DE_summary.tsv：8现有组的样本量、可估计状态和上下调数。
- results/within_celltype/：各组细胞级预筛全表/显著表、可估计组pseudobulk全表/显著表及年龄性别敏感性。
- Within_type_combined_gene_evidence.tsv：267,896条基因×类型并列证据；Within_type_DE_robustness_summary.tsv：调整前后数量及方向一致性，无新筛选门槛。
- All_celltypes_pseudobulk_all_tested.tsv和All_celltypes_pseudobulk_age_gender_all_tested.tsv：各类型合并统计。
- All_celltypes_celllevel_all_prefiltered.tsv：双向Seurat描述；不是未经表达预筛的全基因双向检验表。
- Library_QC_summary.tsv、Celltype_condition_counts.tsv、Donor_celltype_composition.tsv、Cellcycle_by_celltype_condition.tsv：基础统计。
- PCA_variance_explained.tsv、PCA_top_absolute_loadings.tsv、Global_cluster_annotation_composition.tsv、Global_resolution_transition_counts.tsv：降维与聚类解释。

## 元数据与缓存

metadata/All_main_metadata_with_global_clusters.tsv含236,309细胞及全部分辨率；All_main_UMAP_pca.tsv和All_main_UMAP_harmony.tsv分别为未整合/整合全量坐标；Balanced_landmark_tSNE_Harmony.tsv仅28,558子集坐标。

cache/All_main_atlas_embedding_bundle.rds为可恢复降维bundle；All_main_Harmony_SNN.rds为全局图；Balanced_landmark_PCA_model.rds包含均值、数值标准差、载荷和训练cell IDs；canonical_by_library为33,487全基因计数分片。原逐库完整Seurat不改写。

## 图

Foundation_01–02：逐文库QC/保留；03–04：PCA/Harmony/碎石；05–08：全局UMAP、整合比较、子集tSNE、按文库拆图；09–10：周期；11：已知marker点图；12：供者组成；13–14：聚类分辨率/注释；15–16：疾病DE汇总/火山；17_*：各可估计组供者热图。每图PDF及PNG。

## 审计与恢复

audit/Frozen_stage1_before.tsv及Frozen_stage1_after.tsv验证7份原主结果SHA256；Independent_validation.tsv为20项数据/模型口径检查；Completion_metrics.tsv为实际数量；Frozen_stage1文件校验需在本轮最后执行。

scripts顺序01准备→03疾病DE/03b并列证据→02合并图谱→04出图→05验收→06报告。各模块有cache/*.complete及logs。原件不删；本轮只新增目录。稀有组/诊断组限制见报告。
