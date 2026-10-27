# 2026-10-26：长期科研上下文与任务交接机制

任务标识沿用起始日期；文档收尾日期：2026-10-27（Asia/Shanghai）。未把旧科学结果日期改为本轮日期。

## 任务与基线

- 用户授权：在已整理项目基础上完善AGENTS、统一进度入口、方向日志及云端启动模板；只修改文档，不重审全项目、不改数据/统计/核心参数、不开始分析、不push。
- 已读：根AGENTS/README、科研进度当前状态与相关章节、迁移入口、现有模块报告文件名及阶段配置。科学摘要依据已有档案/结果入口，不声称本轮重新审阅全部原始证据或重算结果。
- 起始/收尾Git：项目根`git rev-parse --show-toplevel`返回不是Git仓库；分支/HEAD/工作区Git diff不适用。未核实远端仓库、审核或合并证据，已共享基线未知。
- 本任务新增修改同步状态：`LOCAL_ONLY`。这是本地记录，不是已发布云端状态。
- 计算状态：文档工作完成；科学状态：不适用，无新科学计算。原有阴性结果和停止边界保留。

## 实际修改

- [AGENTS](../../AGENTS.md)：精简恢复顺序、定向阅读、独立任务日志、审核合并共享基线及交接要求；保留科学与数据规则，细化内容指向模块日志。原版保存于[archive](../archive/AGENTS.pre_context_handoff_2026-10-26.md)。
- [RESEARCH_PROGRESS](../RESEARCH_PROGRESS.md)：增加统一模块导航和本地/共享状态说明，既有科学正文不改。
- [README](../../README.md)及[CLOUD_MIGRATION](../CLOUD_MIGRATION.md)：增加任务规则、日志与启动模板入口。
- 新建[TASK_LOGS规则/索引](README.md)、[疾病图谱](01_DISEASE_AND_ATLAS.md)、[预测对接](02_PREDICTION_AND_DOCKING.md)、[表达程序](03_MFSD10_PROGRAMS.md)、[工程云迁移](04_ENGINEERING_AND_CLOUD.md)、[空白交接模板](TEMPLATE.md)、本记录及[云端启动指令](../CLOUD_TASK_TEMPLATE.md)。方向摘要是历史证据接手索引，不是本轮运行成果。

## 验证、限制与恢复

- 已做定向文档验证；本轮未执行全项目审计、依赖安装或生信程序。
- `python --version`遇WindowsApps别名无法启动，使用PowerShell/.NET完成定向文档验证，未改变系统环境。
- 无本任务启动的后台分析或待恢复进程；其他会话任务未核实。恢复入口为AGENTS→科研进度导航→本目录相关模块/任务。
- 未进行Git初始化/暂存/提交/推送/建PR/合并。未经审核合并前，不升级共享状态。外置数据位置和Linux完整复现仍沿用前轮未核实事项。
- 下一步建议：用户在新任务中粘贴[CLOUD_TASK_TEMPLATE](../CLOUD_TASK_TEMPLATE.md)并填写具体目标；仓库建立后由获授权任务核实审核/合并状态并追加同步确认。此建议不是自动执行Git操作或新分析授权。

## 收尾验证（2026-10-27）

- PowerShell/.NET读取13份本轮相关文档（含原AGENTS快照），178个本地Markdown链接均存在，未发现冲突标记。科研进度导航与`config/stages.json`的11个阶段逐项对应。
- 首次链接检查发现候选选择阶段不存在独立`analysis/07_candidate/scripts`目录，已将两处新导航改为“无独立脚本、依据已有结果汇总”，没有创建虚构代码或运行记录。
- 人工核对恢复顺序、停止边界、患者独立性、阴性记录和共享状态定义一致。历史审计/报告的“本轮”仍是其原任务快照，不冒充本次新增验证。
- 科研进度从“数据与来源”到文末的原科学正文SHA256前后相同：`855515075B90E91ED556C738EF8A83DB11AA13FE1D5ACC10E7A9A8A586875B4F`。
- 定向保护检查前后相同：`config/stages.json`、`docs/PATH_MAP.tsv`、`results/10_programs/frozen/Program_genes_v1.tsv`、`results/11_aggregation/audit/STOP_THIS_ROUTE.json`。这是有限范围校验，不宣称重新校验全部数据和结果。
- `.gitignore`中无排除新任务Markdown文档的规则；未修改它。Git根检查仍为非仓库，无法提供项目Git diff/HEAD，没有提交或推送。

## 后续同步确认

尚无。本记录及本轮全部修改仍为`LOCAL_ONLY`，无commit/PR/审核/合并SHA可填写。
