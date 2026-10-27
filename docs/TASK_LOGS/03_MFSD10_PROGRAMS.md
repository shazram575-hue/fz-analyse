# MFSD10表达、程序与探索性聚合：接手记录

记录日期：2026-10-26；类型：既有交付摘要。本轮无新统计。同步快照：`LOCAL_ONLY`，无可核实审核合并基线。科学状态：程序重复性有支持；MFSD10关联不足，停止自动优化这一路线。

## 按实际阶段恢复

| 模块 | 证据入口 | 代码入口 |
|---|---|---|
| 08 表达分布摸底 | [报告](../../results/08_expression/MFSD10_表达分布摸底_结果报告.md) | [脚本](../../analysis/08_expression/scripts) |
| 09 恶性上皮程序可行性（历史前置阶段） | [患者可用性](../../results/09_feasibility/患者可用性表.md)、[原方案](../../results/09_feasibility/方法依据与最小试运行方案.md) | [脚本](../../analysis/09_feasibility/scripts) |
| 10 正式程序发现及独立外部评估 | [报告](../../results/10_programs/结果报告_程序发现与独立验证.md)、[固定基因](../../results/10_programs/frozen/Program_genes_v1.tsv)、[发现结果](../../results/10_programs/discovery/MFSD10_patient_level_associations.tsv)、[外部结果](../../results/10_programs/validation/Discovery_external_comparison_all.tsv) | [脚本](../../analysis/10_programs/scripts) |
| 11 两队列已查看后的探索聚合 | [三层对照](../../results/11_aggregation/简洁对照表.md)、[摘要](../../results/11_aggregation/结果摘要_探索性聚合再分析.md)、[固定规则](../../results/11_aggregation/探索性再分析_固定规则.md)、[停止标记](../../results/11_aggregation/audit/STOP_THIS_ROUTE.json) | [脚本](../../analysis/11_aggregation/scripts) |

原目录、完整参数/版本/患者映射和完成标记见[科研进度第5–6节](../RESEARCH_PROGRESS.md)、[PATH_MAP](../PATH_MAP.tsv)及[stages](../../config/stages.json)。09方案中的“尚未运行”是历史快照，不能覆盖10/11实际交付状态。

## 已有发现及不可越过的边界

- MFSD10在恶性上皮单细胞稀疏且受深度影响；不按高低或阳性/阴性分组。LN/癌前混合细胞与CNV未定/未评分边界沿用[疾病模块](01_DISEASE_AND_ATLAS.md)。
- cNMF发现用GSE161529确认恶性singlet的13患者6071细胞，2000HVG、K4/6/8各100重复起始方案，3900主分解和19场景380次敏感性已交付；A/B不是方法充分标准。5个固定程序支持患者11/5/5/3/3，不按MFSD10关联筛程序。
- MFSD10从拟合输入、归一化分母和末端全基因重拟合排除。5程序基因与评分冻结，两个队列未联合整合发现。统计不用integrated或UMAP代替counts，不用MAGIC/ALRA补值作为主要证据。
- GSE176078主对象clinical_v2；26病例映射均保留，其中20有作者Cancer候选，16Naïve主要验证、4Treated单列。doublet状态未知；CID3963的metadata/IHC分型冲突保留，不推PAM50。作者Cancer标签与已确认所有上皮恶性不是同一回事。
- 外部基因协同性有支持，不把“评分高”单独称程序复现。发现MFSD10关联未过BH（最小q0.207），没有确认性成功；外部MP05弱正相关只探索性。迭代上限及0043深度敏感性限制保留，不能称全部程序分量稳定。
- 后续固定约20细胞非重叠聚合全部探索性；每患者内估效应再患者层汇总，不把聚合单元当独立患者。主要发现13、外部15患者；4465仅6单元按既定规则不可估，旧外部16患者结果不抹去。
- 聚合减少零值但发现全部5程序均未过校正，患者间仍未支持；MP02/MP05外部增强、发现不精确。所有效应、CI、FDR及技术敏感性继续以全5程序对照表为准，不只看有利项。

## 停止与恢复

本轮只读档案/证据入口，没有重跑NMF、筛signature、调整模型或启动后台分析。既有计算完成与科学证据不足同时成立；停止不等于证明MFSD10无功能。未经新的明确科学授权，不继续优化表达程序关联、不自动追加GSEA/通讯/拟时序。后续工程任务先核实数据和运行条件，不删除完成标记强制重跑。
