# GitHub＋Codex Cloud迁移

## 当前可做与未完成项

本项目准备的是**可审阅、可轻量校验的源码与证据仓库视图**，不声称Windows生信流水线已在Linux重现。2026-10-26整理阶段未初始化Git；2026-10-27用户新授权首次同步到[指定公开仓库](https://github.com/shazram575-hue/fz-analyse)，最新Git/上传状态见[同步日志](TASK_LOGS/2026-10-27_github_sync.md)。尚未创建或验证Cloud环境。

2026-10-27最新：用户明确允许公开上传，代码和小型结果已通过GitHub连接器同步至main并独立fetch核验，原本地草稿保存在archive分支。此前拦截/浏览器失败是历史，不再等待私有化。可开始创建Cloud环境；实际传输及校验证据见[同步报告](GITHUB_SYNC_REPORT.md)。本机Git CLI仍未配置登录凭据，不把其失败的push命令写成成功；Cloud自己的仓库访问权限另行实测。

原目录约35.48 GB，其中原始数据8.05 GB；Windows环境与缓存不能原封不动搬入GitHub或小型云VM。官方资料说明Cloud使用可复用环境，需选择仓库、准备依赖/访问并完成测试后发布；本轮没有实际执行这些账户操作。[OpenAI官方Cloud环境说明](https://learn.chatgpt.com/docs/environments/cloud-environments)

## 建议迁移顺序

每次Cloud任务使用[启动指令](CLOUD_TASK_TEMPLATE.md)，按[任务记录与同步规则](TASK_LOGS/README.md)恢复上下文。历史`LOCAL_ONLY`不冒充已同步，代码上传或PR存在仍不足以称科研状态已审核合并。代码同步也不代表Git忽略的数据/对象已恢复。全项目状态统一见[科研进度](RESEARCH_PROGRESS.md)。

1. 本地先运行`python scripts/project/doctor.py --evidence-only`和`python scripts/project/verify_project.py`。检查[上传政策与候选表](GITHUB_FILE_POLICY.md)，不要直接将整个本地目录拖到GitHub。
2. 用户已指定仓库；首次初始化、审阅/提交/推送结果按同步日志核实。后续仍须检查`git status`、分支/上游、远程新增、大文件和密钥，不能按历史状态盲推。
3. 在Cloud选择该仓库并创建环境；只先做审阅与轻量校验。设置/安装脚本可参考`scripts/cloud/setup.sh`，默认不安装大型分析栈，不下载数据，不启动分析。
4. 只有需要某分析阶段且有明确新授权时，重建对应Linux环境、恢复其数据负载，再做最小兼容性测试。不要因“移植”而重估已完成科学结果。

当前Cloud官方文档区分新环境与Legacy环境；网络访问、账户权限和资源应在所用环境实测，不能沿用旧教程假设。旧环境脚本的默认离线/缓存行为不应套用于所有新环境。[当前文档](https://learn.chatgpt.com/docs/environments/cloud-environments)、[Legacy文档](https://learn.chatgpt.com/docs/environments/cloud-environment)

## 代码与证据路径

`analysis/<stage>/scripts`、`results/<stage>`是原件的字节相同副本。旧代码没有被强行改路径；[PATH_MAP](PATH_MAP.tsv)列出旧outputs位置。

```bash
python scripts/project/restore_legacy_layout.py --dry-run
# 新克隆确需旧小文件布局时才执行，全部先验hash/冲突检查
python scripts/project/restore_legacy_layout.py --apply
```

此工具只恢复映射中的小文件到项目内outputs；已有内容不一致就拒绝，不覆盖。它不会恢复大型对象、模型或文献，不能自动解决硬编码解释器、R库和系统依赖。恢复后outputs仍被Git忽略，避免把本地工作负载误上传。

新工具从自身位置推断项目根，可检查`BIOINFO_DATA_ROOT`数据位置；旧脚本不一定读此变量。`config/project_paths.json`是新约定，不是对旧硬编码的神奇修复。各阶段移植需在保留旧代码基础上建立新入口、明确输入/输出及版本，改前解释依据。

## 环境与Linux兼容问题

| 事项 | 已核实状态 | 云端待办 |
|---|---|---|
| Python | 三个独立Windows Python3.12.14环境已采集包清单 | 分开重建，不能把programs/DeepPurpose/docking版本混装 |
| 程序环境 | requirements.txt为实际核心版本，非完整锁 | Linux包解析/二进制wheel和核心导入仍需实测；不称lock已验证 |
| DeepPurpose | Torch2.5.1+cpu/DeepPurpose0.1.5等见observed清单 | PyTorch CPU来源、RDKit、DGL/MPNN及旧编码接口分别核验；不能把安装成功当准确性验证 |
| R | 实测R4.5.3；两个R库DESCRIPTION共1028条，含重复包及不同Built版本 | Linux从源/对应仓库重建；不要复制Windows DLL；不伪造renv.lock |
| R包枚举 | 本次installed.packages返回空，直接DESCRIPTION回退可读 | 属本机会话可见性/枚举限制，不是“无依赖”；加载验证与历史成功日志分开 |
| 路径/启动 | 基线57个脚本含Windows盘符；10个PS、57个R脚本语法检查通过 | `.libPaths(D:...)`、`.exe`、PowerShell、English_United States locale等需适配，语法通过不代表可运行 |
| 结构软件 | 旧Vina/P2Rank/Java及准备脚本有Windows工件 | 获取对应Linux可执行文件，固定版本及输入SHA，先小规模技术核验 |
| 数据/内存 | 大对象外置，本轮未创建云数据持久层 | 实测VM/cgroup内存、磁盘和下载/解压峰值，逐库稀疏处理，禁止完整巨矩阵密集化 |
| 权限/网络 | 没有验证用户Cloud账户、仓库或持久存储权限 | 依所用Cloud环境配置官方GEO/GDC/包源及必要私有访问，不把凭证写仓库 |

官方当前默认VM资源随计划而异，本地目录已超过若干默认磁盘规格；应以环境实测为准，不能把整目录35.48 GB作为启动负载。[官方资源说明](https://learn.chatgpt.com/docs/environments/cloud-environments)

## 数据恢复

公共数据按[DATA_MANIFEST](DATA_MANIFEST.md)及`config/public_datasets.json`规划。默认下载命令只列计划，明确`--execute`后流式下载并校验，绝不覆盖原件；大对象外部存储地址目前未知。新鲜公共下载若与冻结哈希不同，报告差异并停下，不能换成“最新版”继续旧比较。

Windows包环境、模型权重、原始NMF重复、Seurat/AnnData对象不普通Git。GitHub会阻止超过100 MiB的普通文件；LFS/外部存储是另行设计，不是本轮已配置完成。[GitHub官方大文件说明](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)

## 下一位Cloud智能体的首条任务

本次只在Windows实测Python3.12.14工程入口；默认GBK首次失败，使用`python -X utf8`复核通过。Cloud setup已设置`PYTHONUTF8=1`；Linux环境本身、R依赖和生信入口仍未实测。工程核验只需Python3.10+标准库；无需为了读取证据安装整套cNMF/DeepPurpose环境。

数据下载/依赖安装时按实际阶段申请对应网络访问：官方GEO/NCBI、GDC、PyPI及相应R/Bioconductor包源；预测/对接另核实冻结权重/结构来源。不是所有域名都要无条件开放；不得为继续任务上传凭证或自动下载受控EGA reads。

读AGENTS和RESEARCH_PROGRESS；运行证据级doctor/verify；确认5程序和停止标记；列明当前VM缺哪些数据、依赖和系统库。可继续做工程移植/恢复规划、检查文档与结果一致性。**不要重跑NMF、重选signature，或继续优化已停止的MFSD10表达程序关联。** 新的生物学目标或湿实验路线需要用户新指令。
