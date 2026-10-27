# 本轮入口：恶性上皮表达程序可行性

2026-10-26。已完成可行性检查；未运行NMF/cNMF或程序投射，未新增过滤或改动Seurat对象。

- [32患者可用性表](患者可用性表.md)
- [完整机器可读表](results/患者可用性表.tsv)，含GSM、全部原发病灶、确认恶性与singlet复杂度、双细胞数、用途及理由。
- [依据与未执行的最小试跑提案](方法依据与最小试运行方案.md)
- [11项计数/QC检查](audit/Counts_QC_validation.tsv)
- [按样本恶性细胞构成](results/Confirmed_malignant_sample_composition.tsv)
- [审计复现入口](scripts/run.ps1)

当前32患者9,486确认恶性，其中8,722 singlet/764 doublet标记；A发现候选10、B小K候选3、C仅投射候选5、D低复杂度描述优先4、E仅描述10。分类依据为本项目保守提议，不是公认阈值或已验证成功。当前无本轮后台任务；下一步只能按用户新授权试运行。
