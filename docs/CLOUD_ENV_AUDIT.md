# Codex Cloud首次启动：科研上下文与环境审计

检查日期：2026-10-27（Asia/Shanghai）。项目：Spiropachysine A＋乳腺癌，工作流3，目标期刊Phytomedicine。

**结论：已准确恢复仓库中的研究上下文，云端可安全开展文档审阅、结果表核对及轻量工程工作；尚不具备完整生信分析或原流程复现条件。** 728份受控副本哈希通过，冻结程序和停止标记完整。原始/大型数据未恢复，R及主要分析依赖未就绪，GEO/GDC/R包源访问受限，持久存储未证实。本次没有执行新科学假设验证，没有重跑NMF、预测、对接或关联检验，没有提交、推送、建PR或合并。

本任务交接及命令记录见[TASK_LOGS](TASK_LOGS/2026-10-27_cloud_first_start_audit.md)。本报告和本任务文档增量均为`LOCAL_ONLY`，独立任务不能假设这些未提交文件自动可见。

后续保存说明：上述同步状态是首次审计收尾快照。用户随后另行授权独立分支提交、推送与PR，实际保存结果以[同步交接](TASK_LOGS/2026-10-27_cloud_audit_github_sync.md)及远程Git记录为准，不改变本报告的检查事实或科学结论。

## 1. 仓库、基线与检查范围

| 项目 | 本次观测 |
|---|---|
| 启动目录 | `/workspace`，此目录本身不是Git工作树 |
| 实际仓库 | `/workspace/fz-analyse` |
| origin | `https://github.com/shazram575-hue/fz-analyse.git`，与README及GitHub同步报告指定地址一致 |
| 当前分支 | `work`，没有配置上游；本地有`origin/main`引用 |
| HEAD | `655cebaf3857ae0bff107e1448dc95649d4c9264` |
| 远程实际状态 | 只读`git ls-remote --symref origin HEAD refs/heads/main`成功：默认分支main，HEAD及main均为上述SHA |
| 起始工作区 | 干净；832个已跟踪路径 |
| 相关历史 | `58f54b6`初始提交 → `9342813`代码/证据导入 → `655ceba`同步收尾文档 |
| 共享与审核 | 可确认当前checkout与远程main同SHA；没有取得额外科研PR审核/合并证据，不标为`MERGED_REVIEWED`。按项目规则保持`BRANCH_UNMERGED`，含“直接同步main但未取得科研审核证据”的含义 |
| 仓库写权限 | 未测试；只读成功不证明可push，本次也无写远端授权 |

依次读取了AGENTS、README、RESEARCH_PROGRESS、PROJECT_AUDIT、DATA_MANIFEST、CLOUD_MIGRATION、CLOUD_TASK_TEMPLATE、TASK_LOGS规则/模板、全部既有模块与独立任务日志；并进一步核对实际表格、程序冻结文件、停止标记、方法实现说明、依赖清单、下载配置和工程脚本。历史审计与同步记录只作为各自日期的证据，不将其中“尚未建仓/尚未创建Cloud”沿用为本次实况。

未读取本地未上传的论文全文、对象、权重或原始大数据；未逐一核查每份历史科学日志，未重新估计统计量。本次表格计数、读取既有效应/FDR和哈希核对属于证据审计。

## 2. 真实研究进度与证据边界

核心问题是Spiropachysine A可能通过何种靶点和细胞背景影响乳腺癌；当前MFSD10/SLC75A1（ENSG00000109736、Q14728）是用户采用的候选，直接结合、功能必需性、作用方向与因果机制均未建立。

| 已有工作 | 仓库证据与本次核对 | 能支持与不能支持的结论 |
|---|---|---|
| 疾病侧筛选 | [原Pool表](../results/01_disease/results/Pool_A_candidate_genes.tsv)实际813行；TCGA STAR counts的TMM＋limma-voom，与GSE161529确认恶性/健康上皮比较 | 疾病表达筛选；不是药效或因果机制证据。TCGA目录虽叫“芯片”，数据实际为RNA-seq |
| 图谱及后续疾病输入 | [1881表](../results/03_deeppurpose/main_input/Disease_matching_candidates_1881.tsv)实际1881行；六大类同类型Tumor-origin−Healthy双向差异与Bulk交集 | 1881和原813分别冻结，不互换。少供者细胞组保留描述性限制 |
| 药物预测 | 池内1993配对×10模型；独立20431蛋白×8非MPNN模型；[阈值表](../results/04_proteome/results/Threshold_comparison_1_5_10percent.tsv)核对双终点共有∩疾病1%/5%/10%为2/17/47 | 用户采用10%得到47；不是模型准确率验证。MPNN批次依赖、共享BindingDB偏置、输入截断限制保留 |
| 全47反向对接 | [整合表](../results/06_docking47/results/47_integrated_docking_disease_evidence.tsv)实际47行；报告记录282主对接＋9参照，Vina1.2.7、exhaustiveness64 | MFSD10采用6S4M实验结构，但口袋为预测；两状态分数−11.683/−11.934为计算分数，不是实测亲和力。没有本线MD交付 |
| 表达摸底 | [报告](../results/08_expression/MFSD10_表达分布摸底_结果报告.md)：9486确认恶性细胞中60.3%零值，检出细胞55.6%仅1UMI | 测量稀疏且受深度影响，不按MFSD10高低/阳阴划组推机制 |
| 程序发现及外部评估 | 13患者6071细胞逐患者cNMF，3900主分解及380次局部敏感性；[冻结程序](../results/10_programs/frozen/Programs_v1.json)为5×50基因，目标排除；外部主要16 Naïve患者 | 程序共表达重复性有支持；MFSD10关联未获得发现后独立确认。两队列未联合整合发现 |
| 探索性聚合 | [主10项表](../results/11_aggregation/results/Aggregation_main_all10.tsv)、[完整30项表](../results/11_aggregation/results/Comparison_all30_effects_CI_P_FDR_N.tsv)及[停止标记](../results/11_aggregation/audit/STOP_THIS_ROUTE.json)已读取 | 两队列已查看后的探索性工作。外部部分增强，发现仍不精确；停止自动优化，不能反称确认成功 |

关键阴性/不稳定结果完整保留：

- [发现关联10项](../results/10_programs/discovery/MFSD10_patient_level_associations.tsv)最小BH q为0.207087；[外部对照10项](../results/10_programs/validation/Discovery_external_comparison_all.tsv)的confirmatory与confirmatory_supported均为False。外部MP05单细胞患者内rho约0.057、q约0.00290只属探索性。
- 聚合主家族为两队列共10项BH。MP02：发现rho 0.088，95%CI [−0.030, 0.201]，q0.379，n13；外部rho0.285，[0.174, 0.421]，q0.00137，n15。MP05：发现0.136，[−0.037, 0.285]，q0.370；外部0.282，[0.164, 0.405]，q0.00137。全部五个发现程序未过校正，患者间模型也未获得支持。旧q与聚合q家族不同，不只比较q大小宣称改善。
- 主要原外部n16；聚合n15是4465仅6单元、未达固定8单元条件，不是按结果挑患者。患者是生物学重复，细胞/聚合单元不增加患者n。
- cNMF部分重复触及迭代上限；0043深度扰动稳定性有限。P值符号翻转与患者bootstrap CI构造不同，不能只挑CI不跨零的项。[实现说明](../results/10_programs/audit/方法实现与核读说明.md)保留K选择描述与实际代码的差别。
- 对接报告保留距离检查268/282通过及未通过项；参照PPARG最高姿势2/3、PNP3/3、SLC18B1为0/3，不将技术验收概括为所有位点均已验证。
- PubMed零命中是2026-10-24检索快照，不是永久创新性证明。冻结UniProt2026_03、GENCODE v36不能换成“最新版本”后仍声称原比较复现。

细胞背景需要单独澄清：47整合表中MFSD10的SC_up_cell_types及PB_support_cell_types均为Fibroblast；[候选证据](../results/07_candidate/现有证据支持的候选优先级.md)记载成纤维PB log2FC约0.948、BH约7.40×10⁻⁵、调整后约1.01×10⁻⁴。这是既有疾病背景线索，不表示MFSD10仅在成纤维表达，也不等于药物通过成纤维作用。后来的恶性上皮程序阴性结果不能覆盖此背景差异。

GSE161529唯一主CNV分类仍为17c；未定/未评分不当正常或恶性。GSE176078的clinical_v2、16 Naïve/4 Treated分层、doublet未知及CID3963注释冲突保留。停止标记明确不授权任何额外分析；本次未调整粒度、患者、signature、模型，未追加GSEA/通讯/拟时序，未恢复被清除的历史阶段。

## 3. 云端可访问与缺失的资料

| 类别 | 本次确认 |
|---|---|
| 文档与代码 | 所需入口全部可读；analysis目录139文件、772785字节；各阶段R/Python/PowerShell脚本可审阅 |
| 小型结果 | results目录591文件、67080373字节；候选、阈值、对接、患者映射、冻结程序、关联/敏感性表及图形可访问 |
| 原件/副本 | verify逐项验证PATH_MAP的728份受控副本；与记录SHA256一致。不是与云端不存在的旧outputs现场再次对照 |
| 下载与校验清单 | data目录含README及1208项TCGA UUID/URL/MD5/SHA清单；GEO URL/体积/哈希配置可读；1214原始文件校验表及旧INTEGRITY_SUMMARY可读 |
| 原始数据 | `原始数据/`不存在，未恢复；本次未重验原始文件哈希。旧审计“TCGA MD5失败0”是历史证据 |
| 大型派生负载 | `outputs/`不存在；两个doctor指定H5AD均缺失。Seurat/RDS、counts、NMF原重复、完整全蛋白分数/序列、权重/结构工件不在当前可用负载 |
| 外置数据 | BIOINFO_DATA_ROOT未配置；未发现项目专用持久数据挂载证据或已给出的对象存储URI/权限/校验映射；不推断未检查的其他位置有数据 |

只读恢复dry-run显示728份小文件可恢复、已恢复0。没有执行`--apply`，因为本次审阅不需要重建旧目录；该工具不会恢复大对象或修复Windows解释器。原始8.05 GB、原项目约35.48 GB为历史负载，不是当前97 MiB工作树的大小。

## 4. Python、R与Shell环境

实测平台Debian GNU/Linux 13.6、Linux x86_64、glibc2.41。Python3.12.14（Clang构建），可执行文件位于`/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python`；与本地Windows同Python版本号不表示同包/ABI环境。Bash5.2.37、Git2.52.0、Java OpenJDK21.0.12.1可用。

| Python包 | requirements记录的Windows版本 | 当前Linux元数据 |
|---|---|---|
| numpy | 2.5.4 | 2.3.5 |
| scipy | 1.18.1 | 1.17.0 |
| pandas | 3.1.0 | 2.2.3 |
| scikit-learn | 1.9.1 | 1.8.0 |
| matplotlib | 3.11.2 | 3.10.8 |
| pyyaml | 6.0.3 | 6.0.3 |
| anndata / statsmodels | 0.13.4 / 0.15.0 | 均无已安装distribution元数据，doctor亦未找到模块 |
| h5py / cnmf / scanpy | 3.16.0 / 1.7.1 / 1.12.4 | 均无已安装distribution元数据 |

numpy/scipy/pandas/sklearn/matplotlib/yaml在每命令单线程条件下成功导入，2×2稀疏矩阵构造检查通过；没有分析科研数据。Torch/RDKit/DeepPurpose/DGL/Biopython无已安装distribution元数据。没有把find_spec或元数据存在当成全部依赖可加载的证明。

R和Rscript不在PATH，未取得可用R运行时；不能执行包加载、R语法或统计模型检查。仓库内R4.5.3与1028条DESCRIPTION是Windows历史证据，不是此云端R环境。后续按阶段重建data.table/Matrix/edgeR/limma及所需Seurat/harmony/infercnv，分别核验版本、系统库和ABI；无需为轻量表格工作安装全栈。

Conda/mamba/micromamba、PowerShell、Vina/P2Rank未在PATH找到。Java21存在，但旧对接脚本调用Java17 Windows EXE，未核验替代兼容性。Docker命令存在，未检测守护进程或创建容器；nvidia-smi不存在，GPU能力未证实，后续规划按CPU考虑。

Matplotlib初次导入尝试写只读`/home/agent/.config/matplotlib`，回退临时缓存且导入成功；按命令设置`MPLCONFIGDIR=/tmp/cloud-first-audit-mpl`后配置路径可用。后续可选项目可写缓存目录；临时目录不当持久结果目录。未修改系统环境、安装软件或重写requirements；Linux固定版本可获取性、完整依赖解析/二进制兼容仍待专项验证。

## 5. Windows/Linux兼容与资源

### 路径和启动问题

本次对analysis中`.py/.R/.ps1/.sh`使用正则`[A-Za-z]:[\\/]`静态查找，78个文件命中盘符字面量。它包含脚本内字符串/文案，不是78个经动态验证的失败；与历史全目录“57个源脚本”统计范围和规则不同，不解释为新增错误。

实际例子：

- `analysis/01_disease/scripts/05_bulk_limma_voom.R:3`调用D盘`.libPaths`。
- `analysis/10_programs/scripts/run_pipeline.ps1:4-14`调用`env/Scripts/python.exe`、C盘Python和D盘Rscript.exe。
- `analysis/06_docking47/scripts/run_all47.py:6`调用`vina_1.2.7_win.exe`；`define_search_regions.py:9`调用Java17的`java.exe`。

另有旧outputs相对位置和locale等依赖。`config/project_paths.json`/BIOINFO_DATA_ROOT是新约定，不能自动修复旧脚本。后续仅为选定新任务建立Linux入口，保留源脚本、输入SHA和统计语义，先用微型技术输入测试；不批量重构。`bash -n scripts/cloud/setup.sh`通过，但只证明Shell语法；本次未执行依赖安装分支。

### 实际资源与可承受任务

| 资源 | 实测及判断 |
|---|---|
| 内存 | cgroup memory.max=8589934592字节，即8 GiB；检查时memory.current约120 MB。主机/proc显示可用约9–10 GB，不能用其代替8 GiB任务上限 |
| CPU | affinity可见3 CPU；cpu.max=`200000 100000`，约2核时间配额，规划以约2核为准 |
| swap | /proc显示0；不能依赖swap避免OOM |
| 磁盘 | /workspace总33770192896字节，剩余约32034521088字节（29.8 GiB，随运行变化）；工作树du约97 MiB |
| 临时空间 | /tmp约4.9 GiB tmpfs，消耗内存，不用来存大型下载/解压对象 |
| 持久性 | /workspace是可写overlay，/workspace/shared目录存在，但未见项目专用持久卷证据。目录存在不证明跨独立任务持久，未做重启/跨任务试验 |
| 线程 | OMP/OPENBLAS/MKL线程变量初始未设置；本次导入测试限制为1。后续命令显式限制并记录 |

可支撑文档/小表、指定基因的列级提取与小型稀疏技术测试。完整GSE161529/全图谱Seurat、全矩阵标准化、多患者并行NMF、全蛋白批量预测不获资源可行性保证。以GSE176078的100064×29733矩阵为例，仅一个float64稠密副本约23.8 GB，已经超过8 GiB，尚不计中间对象；禁止完整稠密化。原项目35.48 GB也大于当前剩余空间。未来按患者/文库分块或磁盘支持的稀疏方案测峰值，达到预算再换更大环境；本次未跑压力测试。

## 6. 网络、公开恢复与存储

已按cloud-environment-runtime技能读取environment_status和`/etc/codex/network-policy.json`：当前观测有效、running/connected，restricted策略已enforced，仅package_managers预设；没有自定义科学数据域名，没有配置的secret/runtime变量/outbound身份，VPN未配置、TCP目标列表为空。没有输出凭据值，也不由空secret列表推断Git公共读取一定不可用。

默认命令沙箱的首次网络探测出现proxy连接失败/Operation not permitted；启用工具支持的命令网络权限后，保留代理和TLS再做同样只读检查，结果如下。初次沙箱错误与下表目的域策略拒绝分别保留，未关闭代理、绕过网络策略或下载正文数据。

| 有界只读请求 | 本次结果 |
|---|---|
| GitHub origin ls-remote | 成功，远端main与本地HEAD一致 |
| `https://pypi.org/simple/anndata/` HEAD | HTTP200；只证明此索引可达，不证明所有固定wheel可安装 |
| GEO GSE161529 features URL HEAD | `Tunnel connection failed: 403 Forbidden` |
| GEO GSE176078归档URL HEAD | 同上；未下载归档 |
| GDC清单首UUID下载URL HEAD | 同上；未下载counts |
| cloud.r-project.org包索引HEAD | 同上 |
| bioconductor.org包索引HEAD | 同上 |

403发生于代理隧道，不可写作官方源404/数据消失；这也不是工具自动审批拒绝。策略中确未包含上述科学域名，本次实际无法恢复它们。解决方案是在环境配置流程增加所需官方域名并复测；如官方重定向需CDN，再依据真实Location补对应访问。不要无条件开放全部域名；本次未修改网络配置。

默认plan-only下载方案可正常运行：

| 数据计划 | 数量 / 字节 | 本次状态 |
|---|---|---|
| GSE161529 | 3文件 / 2359963729 | 配置、哈希和计划可读；负载未恢复 |
| GSE176078 | 1文件 / 558829202 | 同上 |
| TCGA_BRCA_STAR | 1208文件 / 5125979774 | UUID/官方MD5/SHA清单可读；负载未恢复 |

以上公开计数负载合计约8.04 GB，另有作者补充资料、解压副本和中间对象。download脚本只提供流式`.part`、体积/SHA/MD5核对与不覆盖保护，不自动解包；plan成功不证明官方URL内容和冻结哈希今天仍一致。实际恢复须在后续获授权任务中分批执行，校验不符保留差异并停止，不替换成新版本。公开处理counts与受控原始reads分开，本次不访问EGA。

大型派生对象、模型权重/完整结构没有已配置的统一远程恢复地址。应先建立项目专用持久存储：明确URI/目录、权限、版本、bytes/SHA256、只读原件与新结果位置、保留期/预算；再测试所需小样例跨任务访问。无需把整个35 GB项目搬入Git或当前VM。不因对象缺失而重跑已停止的科学分析。BIOINFO_DATA_ROOT只是选择器，旧脚本还须验证接线。

## 7. 技术问题优先级与完成标准

P0表示相应正式分析前必须解决，P1表示按所选阶段解决；不要求为所有方向一次性安装/恢复全部资源。

| 优先级 | 问题及影响 | 建议与可核验完成标准 |
|---|---|---|
| P0 | 科学数据域名被代理403；原始数据不可恢复、R包源不可用 | 按选定任务配置官方域名，保留代理/TLS；HEAD及获授权的小样例校验成功后再恢复所需数据 |
| P0 | 大对象缺失、持久层/地址未知；独立任务无法取得负载 | 用户指定或配置持久位置，建立对象清单与权限，逐件SHA校验；不声称存在即已恢复 |
| P0 | 未提交文档增量不会自动共享；当前无科研审核合并证据 | 先审阅本次文档；用户后续明确授权后走普通提交/审核同步流程。同步后用具体commit恢复；本次不代为推送 |
| P1 | R不可用、Python版本不符/缺包；旧流水线不能运行 | 选定模块后分环境重建，保留实际解析/导入日志，生成真实Linux环境记录；不伪造lock，不混装DeepPurpose/程序/对接栈 |
| P1 | Windows启动/路径和旧outputs依赖 | 保存原件，建立新版本入口，逐项核对映射、哈希、命令与微型兼容性测试；路径修复不改统计语义 |
| P1（大型任务前为P0） | 8 GiB/2核/29.8 GiB无法保证全量峰值 | 先按输入维度和中间对象估预算，分块稀疏、限制线程，任务级测峰值；必要时升级资源或外置大对象 |
| P1（预测/对接专属） | 权重、结构负载/许可及Linux二进制未验证 | 冻结来源/版本/SHA和使用条件；对应Linux工具仅做技术加载检查，不重新计算科研结果 |
| P2 | Matplotlib默认配置不可写、线程变量未设 | 每命令选可写缓存并设线程；本次MPLCONFIGDIR替代路径已验证 |

## 8. 跨独立任务连续性检查

AGENTS → RESEARCH_PROGRESS模块导航 → TASK_LOGS模块/独立记录 → 实际代码/结果的恢复链条可用；模板具备授权、方法/版本、失败与阴性、后台/断点、Git同步和外置数据字段。verify同时检查冻结程序和停止标记，适合无大数据的独立checkout核验。

限制是Git只传已同步文件，聊天记忆、未提交文档与外置数据不会自动传递。旧模块LOCAL_ONLY是历史快照，不应追改成当时已同步；当前远程main文件可读亦不等于科研审核通过。独立任务应记录实际分支/HEAD/差异，不盲目pull、切分支或复用未合并结果。

本次新增唯一任务日志，并给工程模块、任务索引、科研进度和迁移/README入口增加当前审计链接；只改变工程实况说明，不改既有科学正文或原件。交接中明确本任务没有后台生信进程；不替其他任务宣布无运行。现有机制足以恢复已提交研究上下文，但本次没有创建第二任务或重启VM来宣称端到端连续性验证。

## 9. 下一阶段科研计划（仅建议，均未执行）

现有证据较充分的范围是候选池/流程的计算交付、固定程序的部分重复性及输入身份记录；证据不足的是跨队列MFSD10表达程序关联、药物靶向特异性和作用细胞背景；未验证的是直接结合、功能依赖、药效因果与任何下游机制。停止关联优化不等于MFSD10无功能。

| 建议方向 | 科学问题与价值 | 数据/方法与资源 | 限制与准入条件 |
|---|---|---|---|
| **优先：独立细胞背景与遗传依赖评估** | MFSD10是否在可实验建模的乳腺癌细胞或相关基质背景有可检测蛋白和功能依赖？使后续实验对象有证据依据 | 后续选择固定release的DepMap/CCLE表达、CRISPR gene-effect/QC及HPA蛋白/定位资料；现有成纤维疾病证据作背景。Python表格/列级提取、1–2核，目标≤4 GiB RAM，下载量先查manifest | 覆盖和版本目前未核实，需数据域名/许可。肿瘤细胞系不代表CAF；RNA不等于蛋白，CRISPR依赖不等于药物直接靶点。无覆盖/无依赖同样交付 |
| 独立患者层的细胞来源验证 | 现有Bulk升高和成纤维差异是否在真正独立样本中重复，细胞组成是否足以解释Bulk趋势？ | 新队列须排除重复患者/既有两队列；先冻结细胞类型、患者映射、分组、目标及检验家族，再患者级PB/必要协变量。优先取得作者处理矩阵；按患者稀疏处理，资源先评估，完整图谱可能需16–32 GiB以上 | 新问题需明确授权；不重启恶性上皮程序优化，不新选signature。分型/治疗/平台混杂、CAF供者不足可导致不可估，不补造注释 |
| MFSD10结构与药理可检验性审阅 | 当前预测口袋是否便于设计结合/占用验证，序列、亚型和结构环境有哪些缺口？ | 核对冻结6S4M/序列、原网格与坐标SHA、已知配体/功能资料；小结构/表格，CPU和少量内存即可。获取原坐标前只做证据清单 | 不是新增对接排名或MD；未校准分数不证明亲和力，膜环境和实验结构状态限制需保留 |
| 扰动机制数据可行性与验证设计 | 是否有匹配细胞背景的MFSD10遗传扰动数据，能提供比共表达更直接的功能证据？ | 先查有无公开扰动及化合物覆盖；若存在再另行预设对照、时间/剂量、统计单位、全部终点和多重校正，使用独立验证。资源随数据预估 | 当前未确认存在匹配数据；没有就交付证据空缺，不替换近似基因/化合物后冒称靶点验证。不为显著性增加通路检验 |

所有建议先明确问题与估计量，报告效应/不确定性、样本数、全部预定比较、不可估原因和负面结果。新分析不因结果不显著改变患者、阈值或模型。湿实验展望仍是化合物身份/纯度、靶点占用或独立结合测定、遗传扰动及救援/药效依赖，尚无执行或能力核实。

### 最值得启动的第一个正式生信任务

建议题目：**MFSD10在乳腺癌细胞与基质背景中的独立表达/蛋白证据及遗传依赖评估。** 这是候选可实验验证性的评估，不预设某条通路或重新优化现有程序关联。

1. 先冻结问题：MFSD10是否具备支持选择实验细胞背景的表达/蛋白与遗传依赖证据；以ENSG00000109736/Q14728核对跨库映射，固定数据release、访问日期、SHA、覆盖范围和QC规则。
2. 乳腺来源细胞系为遗传依赖主范围，保留全部按来源/QC可用模型及缺失原因；同一细胞系的重复/筛选批次不当独立模型。成纤维背景若无相应功能筛选，明确“未覆盖”，用独立蛋白/表达证据分栏描述，不从癌细胞系外推。
3. 初始以数据覆盖、逐模型gene-effect、表达/蛋白及发布者QC为主要交付；任何正式组间比较、对照集和校正家族先预设，不选择最有利癌种/阈值。基因依赖阈值采用所选发布版本有依据的规则，现不虚构数值。
4. 交付全部模型表、证据分层、技术限制、可考虑的实验背景/不可判定清单，以及下一步验证设计。无依赖、蛋白证据缺失、跨来源冲突均属于有效结论；不增加NMF、GSEA、通讯、拟时序或MD来“补显著”。
5. 开始前完成该任务实际所需数据的网络/持久存储与最小Python环境；不必先恢复全部8 GB原始数据或整个旧R栈。当前8 GiB/2核预期可用于列级表格任务，但须先核实实际文件大小。

本次仅提出方案，未检索或下载这些新增数据，未核实DepMap/HPA具体覆盖，也未执行上述科学判断。需在后续任务明确授权新问题和数据范围。

## 10. 本次轻量验证记录与交接

在`/workspace/fz-analyse`运行，未加下载`--execute`或恢复`--apply`：

```bash
python scripts/project/doctor.py --evidence-only
python scripts/project/verify_project.py
python scripts/project/restore_legacy_layout.py --dry-run
python scripts/project/fetch_public_data.py --dataset GSE161529
python scripts/project/fetch_public_data.py --dataset GSE176078
python scripts/project/fetch_public_data.py --dataset TCGA_BRCA_STAR
bash -n scripts/cloud/setup.sh
git ls-remote --symref origin HEAD refs/heads/main
```

doctor：review_ready=true、missing_review_files=[]、11阶段、原始目录不存在、两个大对象不存在、linux_full_reproduction_verified=false。verify：pass、728副本、5程序、停止边界保留、large_data_loaded=false、scientific_analysis_rerun=false。restore：would_restore=728/restored=0。三个下载计划正常退出。Shell语法通过；网络失败与按支持方式复测结果见第6节。

附加探测使用Python标准库读取版本元数据、/proc及cgroup、shutil.disk_usage、只列所需项目变量、按固定正则静态查路径、CSV读取已有结果；包导入测试设OMP/OPENBLAS/MKL各1线程，25秒超时；网络HEAD逐URL6秒超时，ls-remote25秒超时。没有卸载/安装包或运行setup的安装分支。

恢复入口：本报告 → [任务日志](TASK_LOGS/2026-10-27_cloud_first_start_audit.md) → [RESEARCH_PROGRESS](RESEARCH_PROGRESS.md)及对应模块。没有本任务后台科学计算或待续分析断点；外置数据未恢复。报告变更须经后续审核并获明确Git操作授权后同步，不能将本地完成直接标为已共享。

收尾检查：`git diff --check`通过，7份本任务文档的137个本地链接可解析；科研进度既有科学正文SHA256仍为`855515075b90e91ed556c738ef8a83db11aa13fe1d5acc10e7a9a8a586875b4f`，与HEAD一致。PATH_MAP、阶段配置、冻结程序及停止标记亦逐字节未变。修改仅涉及7份文档；未创建原始数据、outputs或data/external目录。
