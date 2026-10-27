# 全47靶点反向对接已完成

更新：2026-10-25 04:37:35 +0800。47个靶点×两状态×三种子共282次主对接及9次参照复对接全部完成；15/15项技术检查通过。

优先阅读[47靶点完整结果报告](结果报告_47靶点反向对接.md)与[47行结构/疾病/预测整合证据](results/47_integrated_docking_disease_evidence.tsv)。报告逐项回答结构支持、具体问题和建议的下一步，唯一主攻仍待与用户讨论。

- [94组分状态汇总](results/Target_state_summary.tsv) · [282姿势质量检查](results/All_completed_top_pose_checks.tsv)
- [9次参照对实验RMSD](results/Redocking_control_results.tsv) · [种子间姿势RMSD](results/Seed_pose_pairwise_RMSD.tsv)
- [固定方案](47靶点对接_计算方案.md) · [47实际结构与网格](inputs/47_actual_search_regions.tsv)
- [最终技术验收](audit/Final_independent_validation.tsv) · [图形检查](audit/Figure_visual_QA.md)
- [旧135份对接文件删除记录](audit/Authorized_old_docking_deletion_manifest.tsv) · [上游保护检查](audit/Final_validation_summary.json)

原始本药姿势、命令与日志在docking；参照在controls；全部返回姿势的三维SDF在poses；输入、序列映射、修复与辅因子记录在inputs/prepared/audit。公共软件与结构沿用保留目录，旧两候选分数和对应旧结论已清除，没有混入新结果。

当前本轮主计算与检查结束，不重跑已完成282次，不自动增加口袋、分子动力学或锁定唯一靶点。恢复先读本页和根目录新接手状态；需要重建表图时使用本目录脚本，不使用旧两候选报告脚本。
