# 第一阶段输出索引

日期：2026-10-23；状态：完成。课题Spiropachysine A＋乳腺癌。入口报告：[第一阶段结果报告](../../第一阶段_结果报告.md)；接手：[工作流3交接](../../工作流3_交接.md)。

## 主结果

| 文件 | 内容 |
|---|---|
| [Pool_A_candidate_genes.tsv](results/Pool_A_candidate_genes.tsv) | **813个唯一蛋白编码基因，正式Pool A**；Ensembl排序，不是靶点优先级 |
| [Pool_A_gene_ids_symbols.tsv](results/Pool_A_gene_ids_symbols.tsv) | 仅ID/符号，供后续输入 |
| [Pool_A_full_gene_evidence_universe.tsv](results/Pool_A_full_gene_evidence_universe.tsv) | 全33,487规范RNA特征的合并证据；没有通过正向预筛/未检验项为NA |
| [Pool_A_stage1_summary.tsv](results/Pool_A_stage1_summary.tsv) | 逐层数量 |
| [Pool_A_sensitivity_retention.tsv](results/Pool_A_sensitivity_retention.tsv) | 各替代Pool大小、主候选保留数、Jaccard |
| [Pool_A_expression_by_donor_celltype.tsv](results/Pool_A_expression_by_donor_celltype.tsv) | 813候选全主图谱逐供者/细胞类型表达；169,104行 |
| [Pool_A_expression_by_celltype_summary.tsv](results/Pool_A_expression_by_celltype_summary.tsv) | 按组织组和大类汇总，13,821行；供者组至少20细胞；细胞加权及供者均值检测率分列 |

主表中bulk_P为原始P，bulk_BH为BH；sc_Bonferroni为Seurat全特征Bonferroni；PB_*_BH为患者pseudobulk BH。PB支持需正向效应及BH<0.05，不是新增Pool硬门槛。marker_retained_*为各替代分析是否仍通过单细胞规则。bulk_sc_direction_conflict指Bulk下降、单细胞上升；保留不删除。

未通过可检测性/正向预筛的基因不等同不存在；统计表NA不能替换为P=1、logFC=0或非编码。novelty_status全为本阶段未评估，不代表已经文献查新。

## 差异统计

- [Bulk全统计](results/Bulk_all_tested_genes.tsv)，[探索DE](results/Bulk_exploratory_DE_P005_FC1p5.tsv)，[蛋白编码DE](results/Bulk_exploratory_DE_protein_coding.tsv)。
- [113配对全统计](results/Bulk_paired113_sensitivity_all_genes.tsv)，[主模型对配对](results/Bulk_main_vs_paired_sensitivity.tsv)，[模型摘要](results/Bulk_model_summary.tsv)。
- [单细胞主正向预筛统计](results/SingleCell_main_all_positive_prefiltered.tsv)，[通过阈值markers](results/SingleCell_main_eligible_markers.tsv)。该文件不是所有特征的未预筛双向检验表；完整特征和counts仍保留。
- [患者pseudobulk](results/Pseudobulk_main_all_tested_genes.tsv)，[年龄/性别调整](results/Pseudobulk_age_gender_all_tested_genes.tsv)，[扩正常调整](results/Pseudobulk_extended_normals_age_gender_all_tested_genes.tsv)。
- results/SingleCell_without_*、female_only、patient_half1/2、extended_normal_*为敏感性；half为同供者内分细胞，不是独立患者验证。

## CNV主结果：请勿混用旧诊断文件

**唯一主文件：[Tumor_discovery_cells_CNV_coherent_final_calls.tsv](results/Tumor_discovery_cells_CNV_coherent_final_calls.tsv)**；主汇总[CNV_coherent_final_summary_by_patient.tsv](results/CNV_coherent_final_summary_by_patient.tsv)，主规则[说明](audit/CNV_PRIMARY_CLASSIFICATION_METHOD.txt)，脚本scripts/17c_coherent_CNV_final_calls.R。

Tumor_discovery_cells_CNV_calls.tsv为旧整染色体严格门槛，Tumor_discovery_cells_CNV_final_calls.tsv为RMS-only诊断；均不是主RNA差异输入。含final字样不代表优先于coherent主版本。低CNV/不一致细胞为未定；未评分上皮不外推标签。

参考/校准/留出、冻结基因坐标和逐患者inferCNV缓存可在metadata、audit及cache中恢复。健康留出5供者/200个唯一细胞重复用于不同查询背景，不算32组独立验证。

## 供者、质控和注释

- [GDC样本映射](metadata/TCGA_verified_file_sample_metadata.tsv)：对应全部1,208文件。
- [单细胞统一供者映射](metadata/GSE161529_reconciled_patient_metadata.tsv)：69文库、54操作性供者键；主47文库45键。
- [主DE细胞表](metadata/SingleCell_main_DE_cell_metadata.tsv)：9,486恶性＋12,264正常。
- [患者聚合样本表](metadata/Patient_pseudobulk_sample_metadata.tsv)：至少20细胞纳入标记；实际27恶性＋13正常。
- [规范基因注释](metadata/SingleCell_canonical_gene_metadata.tsv)及[来源映射](metadata/SingleCell_canonical_gene_identity_map.tsv)：GENCODE v36；33,487规范特征。
- [主图谱大类](results/SingleCell_broad_annotation_summary.tsv)：全236,309细胞，混合/未定保留；全逐细胞表在cache/All_main_cells_broad_annotation.tsv.gz。
- author_reference下保存全文、补充文件、GitHub快照/SHA及临床提取物；genome_reference下保存GTF和官方校验依据。

## 图与验收

figures/Final_01_candidate_pool、Final_02_sensitivity_retention、Final_03_CNV_donor_classification、Final_04_cell_patient_evidence、Final_05_broad_cell_atlas均有PDF和PNG。总结图已检查；技术pilot及逐文库QC图另存，不能视为投稿最终图。

独立验收：[8项检查](audit/Final_independent_validation.tsv)；定义验收：[6项检查](audit/Pool_A_definition_checks.tsv)；[供者/留出等报告指标](audit/Final_report_metrics.tsv)；[年龄/性别](audit/Pseudobulk_demographic_summary.tsv)。

## 脚本与恢复

主要顺序：00–06输入/元数据/Bulk/QC；07/09/13逐库图谱；10/11作者补充及患者核对；12HPCA；14/18/18b基因注释；15技术CNV验证；16/16b发现集与独立健康控制；17逐患者CNV→17c主分类；19主markers/PB及敏感性；20扩正常；21组装Pool；22独立验收/总结图；23全图谱候选表达。

05使用bulk_block_gls.R；图谱注释使用annotate_library.R。完整会话、时间和原尝试日志在logs。旧文档历史快照在audit，不作为当前待办。常规不需重跑；若只更新显示，运行22而不重新估计生物学集合。

当前没有后台未完成分析。未开展药物预测、GO/KEGG筛选、空间、对接或MD；813候选未完成新颖性核验。
