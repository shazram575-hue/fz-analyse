# 探索性聚合再分析

[简洁对照表](简洁对照表.md)和[结果摘要](结果摘要_探索性聚合再分析.md)已生成。原模型已复用，固定20细胞左右非重叠聚合、全部5程序及技术敏感性已计算；没有重跑NMF或修改signature。

聚合降低了零值，并增强验证MP02/MP05，但发现仍未过校正、区间不精确，未建立可靠跨队列证据。按用户要求停止继续优化本表达程序关联路线，见audit/STOP_THIS_ROUTE.json。

本轮已完成，完成标记为audit/COMPLETED.json：32项核验通过，前轮5473份清单文件SHA256全部不变。当前无本轮后台任务。中途审批服务故障及额度限制造成文档收尾中断，现已补齐；未绕过审批，未重复主体计算。

- [完整30项效应、CI、原始P、FDR及患者数](results/Comparison_all30_effects_CI_P_FDR_N.tsv)
- [全部60个聚合规格](results/Aggregation_all60_fixed_specifications.tsv)
- [患者间精简敏感性模型](results/Patient_between_parsimonious_sensitivities.tsv)
- [逐患者效应及不可估原因](results/Aggregation_patient_effects_and_failures.tsv)
- [相同患者的原单细胞对照](results/Matched_original_singlecell_reused.tsv)
- [固定方案](探索性再分析_固定规则.md)
- [复现入口](scripts/run.ps1)

原始聚合counts、患者pseudobulk、成员映射及状态特征均保留。run.ps1遇完成标记会退出，不自动继续优化。本轮收尾校验另见audit/CLOSEOUT_COMPLETED.json。
