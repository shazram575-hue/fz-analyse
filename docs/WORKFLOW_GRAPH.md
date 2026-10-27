# 输入输出依赖概览

这是主要数据依赖，不代表每个动态文件引用都已静态解析。精确旧目录、代码、报告和完成记录见[stages.json](../config/stages.json)，字面路径线索见本次审计。小文件恢复不等于大数据依赖已恢复。

```mermaid
flowchart TD
 TCGA[TCGA STAR counts] --> BULK[TMM + voom]
 SC[GSE161529 RNA counts] --> QC[作者QC与粗注释]
 QC --> CNV[抽样发现集inferCNV17c]
 CNV --> MAL[确认恶性 vs 健康上皮]
 BULK --> POOL[原Pool813冻结]
 MAL --> POOL
 QC --> TYPES[六大类同类型Tumor来源 vs Healthy]
 TYPES --> D1881[疾病1881 + 诊断76另列]
 BULK --> D1881
 MOL[CID209254完整结构] --> DP[独立20431蛋白 x 8模型]
 DP --> TOP[1% / 5% / 10% 双终点交并]
 TOP --> G47[用户采用10%共有与疾病相交47]
 D1881 --> G47
 G47 --> DOCK[282主对接 + 9参照]
 G47 --> PM[PubMed日期快照]
 DOCK --> MF[用户采用MFSD10候选]
 PM --> MF
 CNV --> DISC[13患者逐患者cNMF; 目标从拟合排除]
 MF --> EXPR[患者分层表达摸底]
 DISC --> FROZEN[5程序与分析规则冻结]
 VAL[GSE176078作者Cancer候选] --> SCORE[独立固定程序评分与基因协同性]
 FROZEN --> SCORE
 FROZEN --> AGG[已查看两队列后的探索性聚合]
 SCORE --> AGG
 AGG --> STOP[跨队列目标关联仍不足; 停止自动优化]
```

DeepPurpose1881池内10模型计算也已保留，但独立全蛋白筛选不是把池内分数直接当全背景百分位。Pool813与1881不是同一集合；程序共表达复现不等于MFSD10作用机制验证。
