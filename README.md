# Spiropachysine A＋乳腺癌：工作流3

候选靶点MFSD10/SLC75A1的计算筛选与表达分析项目，目标期刊Phytomedicine。**候选身份不等于直接结合或因果机制已被验证。** 最新表达程序关联经探索性聚合再分析仍缺少可靠跨队列支持，已停止自动优化该路线。

新任务先读[AI工作规范](AGENTS.md)、[科研进度的结论与模块导航](docs/RESEARCH_PROGRESS.md)、[任务记录规则](docs/TASK_LOGS/README.md)，再定向读相关代码、日志和结果；无需每次重审整个项目。新接手智能体不得假定知道其他聊天，或只凭完成标记推断科学假设成立。工程审计按需查[项目审计](docs/PROJECT_AUDIT.md)。

新开云端聊天可复制[启动指令](docs/CLOUD_TASK_TEMPLATE.md)。任务独立留日志，实质进展更新统一进度，收尾说明实际输出、问题和恢复入口。2026-10-27经用户明确授权，代码、文档及小型结果已同步到[GitHub仓库main](https://github.com/shazram575-hue/fz-analyse)，**远程文件树校验通过，详见[同步报告](docs/GITHUB_SYNC_REPORT.md)**。大型数据仍只在本地。只有已审核、提交并合并到约定集成分支的修改才是共享最新科研状态；上传成功不自动代表新的科研结论通过审核。

## 目录与保留策略

```text
AGENTS.md                    当前项目级规范
docs/                        科研进度、数据清单、迁移及审计
docs/TASK_LOGS/              分方向接手摘要、独立任务记录与交接模板
docs/CLOUD_TASK_TEMPLATE.md  可复制的云端任务启动指令
docs/archive/                整理前根文档快照（历史，不作为执行规范）
analysis/<stage>/scripts/    逐字节代码副本，保留源路径和SHA256
results/<stage>/            报告、冻结程序、主要结果表和图形副本
config/                     阶段依赖、路径约定、公共下载配置
environments/               实测依赖清单；不是虚构的跨平台锁文件
scripts/project/            轻量审计、校验、恢复小文件与下载计划工具
scripts/cloud/              仅环境准备，不启动生信分析
data/                       下载清单；external/为未来外置数据位置
原始数据/                    本地原件，未移动，Git忽略
outputs/                     本地阶段对象/缓存/日志/完整结果，未移动，Git忽略
```

原件和旧代码没有迁移或统计改写。新`analysis`和`results`是受控审阅副本，映射与哈希见[PATH_MAP](docs/PATH_MAP.tsv)。大数据、二进制环境、原始NMF重复等仍在本地，见[数据清单](docs/DATA_MANIFEST.md)。复制的历史报告可能含旧路径或早期状态，以科研进度为当前入口。

## 安全的首次运行

只需Python 3.10+；下面命令不会运行NMF、R差异分析、对接或关联检验：

```bash
python scripts/project/doctor.py --evidence-only
python scripts/project/verify_project.py
python scripts/project/restore_legacy_layout.py --dry-run
python scripts/project/fetch_public_data.py --dataset GSE176078
```

最后一条只打印下载计划，不下载。数据恢复和Linux移植见[CLOUD_MIGRATION](docs/CLOUD_MIGRATION.md)。2026-10-26整理任务未初始化或推送；2026-10-27首次Git同步的新授权和实际结果见[同步日志](docs/TASK_LOGS/2026-10-27_github_sync.md)。尚未创建或验证Cloud环境。Windows默认非UTF-8时使用`python -X utf8`执行上述工程命令。

## 当前科研结果

- TCGA STAR RNA-seq＋GSE161529疾病侧筛选完成：原Pool A813冻结；后续六大类疾病输入1881另存，二者不是同一候选池。
- DeepPurpose池内10模型及独立20431蛋白×8模型预测完成；用户选择前10%双终点共有∩疾病47池。预测共识不等于实验结合证据。
- 全47靶点282次主对接＋9次参照完成；MFSD10为用户采用的候选方向，未运行MD。
- GSE161529逐患者cNMF及GSE176078独立固定程序评估完成，未联合两个队列发现；5程序共表达重复性有支持，MFSD10关联未获得发现校正后再独立确认。
- 后续固定约20细胞聚合统一标探索性：验证MP02/MP05增强，但发现不精确、未过校正；患者间结果未支持。**到此停止继续调粒度、signature或模型，不自动追加其他机制分析。**

## 上传与复现边界

本目录原先约35.48 GB，不可整目录上传。仅考虑本README、AGENTS、config、scripts、environments公开版本清单、analysis及results受控副本和docs；先检查[上传清单](docs/GITHUB_FILE_POLICY.md)。原始数据和Windows环境不适合普通Git；GitHub普通Git会阻止超过100 MiB文件。[GitHub官方说明](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)

当前可供云端审阅和轻量核验，**不等于所有旧Windows流水线已经Linux化或一键复现**。没有伪造renv.lock/environment.yml；requirements.txt只列实测Python程序环境核心版本，Linux重建未验证。
