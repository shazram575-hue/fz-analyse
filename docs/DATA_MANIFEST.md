# 数据清单与获取约定

审计日期2026-10-26。原始数据、对象、缓存和完整结果没有移动。全量文件位置/体积见[科学文件索引](audit/2026-10-26/SCIENCE_FILE_INVENTORY.tsv)，精确字节数以TSV为准。

GitHub同步补记（2026-10-27）：本页列出的原始归档、Seurat/AnnData、大型counts、NMF原始重复、模型权重及完整环境仍为**LOCAL_ONLY_DATA，不在GitHub，也未恢复到Cloud**。仅代码、文档、小型派生表图和校验/下载清单纳入受控上传。公开去标识化作者Case/患者映射随小表保留，不包含受控原始reads；上传不能推定数据负载已可用。实际同步证据见[任务日志](TASK_LOGS/2026-10-27_github_sync.md)。

## 原始数据

原始目录1214文件，共8052105862字节（约8.05 GB）。[每文件SHA256/MD5](audit/2026-10-26/integrity/RAW_INPUT_CHECKSUMS.tsv)覆盖全部原件；本次重新核对1208个TCGA文件，官方MD5失败0。[核验摘要](audit/2026-10-26/integrity/INTEGRITY_SUMMARY.json)

| 数据/格式 | 本地路径 | 体积/数量 | 来源及获取 |
|---|---|---|---|
| GSE161529归档 TAR | `原始数据/单细胞/GSE161529_RAW.tar` | 2359685120字节 | [官方GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529)，计数矩阵/条形码按既有脚本解包 |
| GSE161529基因与SOFT gzip | 同目录features.tsv.gz、family.soft.gz | 270380/8229字节 | 官方GEO补充与SOFT，和计数归档一起保留 |
| GSE176078计数归档 TAR.GZ | `原始数据/单细胞/GSE176078/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz` | 558829202字节 | [官方GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)；基因×细胞稀疏整数矩阵，另metadata |
| GSE176078作者方法PDF/补表XLSX | 同目录41588_2021_911_MOESM1_ESM.pdf、MOESM4_ESM.xlsx | 4865527/2467630字节 | [Wu等原文](https://doi.org/10.1038/s41588-021-00911-1)官方补充；临床与CNV依据。留本地，不默认整篇重分发 |
| TCGA-BRCA STAR counts TSV | `原始数据/芯片/files/<UUID>/<file_name>.tsv` | 1208文件 | [GDC](https://portal.gdc.cancer.gov/projects/TCGA-BRCA)；[精确UUID/URL/官方MD5/SHA清单](../data/TCGA_DOWNLOAD_MANIFEST.tsv)。不是芯片数据 |

公共下载配置在[public_datasets.json](../config/public_datasets.json)。URL来自官方GEO布局/既有GDC记录，本轮没有再次下载大型数据或探测每个下载端点。若公共内容版本已变化，必须报告校验不符，不能悄悄换数据或更新预期哈希。

```bash
# 只查看计划，默认不下载
python scripts/project/fetch_public_data.py --dataset GSE176078
python scripts/project/fetch_public_data.py --dataset TCGA_BRCA_STAR --file-id <清单中的UUID>
# 明确获准恢复公共数据且检查资源后才加 --execute
```

下载默认写入被Git忽略的`data/external/<dataset>`，不覆盖`原始数据`。流式写.part、校验字节数/哈希后再改名；不同内容的已有文件不会被覆盖。脚本不自动解包，需另估解压空间。GSE176078的公开UMI计数不是EGA受控原始测序reads；本项目没有自动获取受控reads。

## 处理对象、缓存和完整结果

| 类型 | 现有位置 | 云端方案 |
|---|---|---|
| TCGA全计数/模型RDS | `outputs/第一阶段_2026-10-23/cache` | 不入Git；外部持久存储或明确授权的重建。保持GENCODE/样本顺序 |
| 逐文库Seurat、上皮counts、CNV评分 | 同阶段cache/seurat_pilot、epithelial_by_library、CNV_discovery_by_patient | 大数据只读；17c最终判定不能由旧RMS-only替代 |
| 图谱规范counts分片、embedding bundle | `outputs/单细胞基础补全_2026-10-23/cache` | 总缓存约6.35 GB；bundle不是合并全基因Seurat assay |
| GSE161529程序发现对象 | `outputs/MFSD10_程序发现与独立验证_2026-10-26/objects/GSE161529_confirmed_malignant_singlets_counts.h5ad` | 原始RNA counts在X；只读，不能加未定/未评分 |
| GSE176078全/恶性候选对象 | 同目录objects，两个clinical_v2.h5ad | 主对象为clinical_v2；v1仅保留历史。对象和clinical map相配 |
| 聚合与患者PB对象 | `outputs/MFSD10_稀疏与聚合探索_2026-10-26/objects` | 原始聚合counts、成员映射和检查点留本地；不为复现自动重做已停止分析 |
| NMF原始随机重复与扰动 | 程序阶段discovery/cnmf、discovery/sensitivity | 不普通Git；固定signature/诊断/完整效应表已在results审阅视图 |
| 全蛋白分数/序列/大型全量DE | 预测、独立蛋白、疾病阶段results/cache | 小摘要已复制；大表按[导出排除表](audit/2026-10-26/EXPORT_EXCLUSIONS.tsv)与科学索引定位 |

上述大型派生文件没有统一公开下载URL，也没有指定外部对象存储地址，**这部分尚未解决**。不能把“可从原始数据重算”当作现已在云端可用，尤其不能无授权重跑大型分析。迁移前应为确需继续的对象建立用户指定的持久存储位置、权限与SHA256。

## 药物、蛋白、模型与结构来源

- 药物计算身份CID209254，InChIKey在科研进度中固定；小型MOL/SDF/核验JSON副本在`results/03_deeppurpose/structure_verification_CID209254`及`results/04_proteome/structure_comparison`。这是计算输入明确性，不是实物纯度/药效验证。
- UniProt2026_03人类已审校canonical背景与GENCODE v36参考保存在预测/独立蛋白阶段。原下载记录和版本冻结协议保留；不要用“最新UniProt”替换冻结版本。完整序列负载与部分下载URL当前只在旧outputs中，迁移清单需按实际继续阶段补齐，不能猜造URL。
- DeepPurpose作者预训练权重/配置/ZIP与原软件环境保留于`outputs/DeepPurpose_预测_2026-10-24`，大权重未入Git。来源核验与加载记录见该阶段scripts/logs/audit。重分发许可、未来可下载性和Linux加载兼容性未在本轮全部核实。
- PDB/AlphaFold结构、P2Rank与Vina本地准备位于docking_audit/reverse_docking47。固定网格与结构来源见`results/06_docking47/inputs/47_actual_search_regions.tsv`；原受体/伙伴/辅因子坐标需和SHA一致，不能只按蛋白名重建后称完全复现。
- 根目录用户提供的Fang2022 PDF及其他正文留本地，来源/重分发许可未全面确认，.gitignore排除PDF/XLSX和归档。查新条目与方法摘要可用于审阅，不能冒称云端已读取未提供全文。

## 上传与校验范围

新仓库视图仅包含2 MiB以内受控副本、脚本、方法/进度、图形及下载清单；不包含原始/大型对象/环境。副本SHA是[PATH_MAP](PATH_MAP.tsv)；原始数据SHA是本页第一节。现有阶段自己的完整交付SHA仍保留在旧outputs，不能因为整理而重写科学文件哈希。
