# 工程整理与云迁移：接手记录

记录日期：2026-10-26；类型：前轮工程成果索引＋本轮机制入口。同步快照：`LOCAL_ONLY`；本轮实查项目根不是Git仓库，未核实远端、commit、审核或合并。不能据此宣称其他位置不存在仓库。

## 既有成果

- [README](../../README.md)、[审计报告](../PROJECT_AUDIT.md)、[数据清单](../DATA_MANIFEST.md)、[迁移说明](../CLOUD_MIGRATION.md)、[上传政策](../GITHUB_FILE_POLICY.md)已存在，本轮保留。
- `原始数据/`和`outputs/`保留原位置；`analysis/`与`results/`为受控审阅副本，源位置/校验在[PATH_MAP](../PATH_MAP.tsv)，不是全流程已完成Linux移植的声明。
- 阶段关系见[工作流图](../WORKFLOW_GRAPH.md)、[配置](../../config/stages.json)。数据校验及工程测试数是前轮审计记录，本轮没有重新执行全项目审计或科学复算。
- 真实环境清单在[environments](../../environments)，不是已验证Linux锁文件。未伪造renv.lock，未确认外置对象云存储、权限与全套Linux环境已恢复。

## 本轮新增机制

[任务规则](README.md)、[全项目进度](../RESEARCH_PROGRESS.md)、[云端启动模板](../CLOUD_TASK_TEMPLATE.md)和[本轮日志](2026-10-26_context_handoff.md)提供定向恢复、任务独立记录、同步标记及交接。AGENTS保留简明强制规则，细节在模块记录和原科学报告中。旧版AGENTS另存[修订前快照](../archive/AGENTS.pre_context_handoff_2026-10-26.md)。

## 尚未解决与下一步

- 未执行git初始化、add、commit、push、建PR、合并或Cloud部署；下一位应实查其checkout和审核合并证据，不沿用本地时间戳当共享状态。
- 先按实际任务做轻量核验；大型处理先查RAM/cgroup、磁盘、稀疏格式、线程与依赖，禁止完整大矩阵转稠密。下载须核对官方来源/校验且有任务授权，外置数据恢复独立记录。
- 本轮Shell中的`python`命令解析到不可启动的WindowsApps别名，文档验证改用PowerShell/.NET，无安装或环境修改；不是宣称项目全部Python环境失效。
- 后续云端工程可在授权范围内补数据持久层、对应Linux环境及最小兼容性测试；不能借移植重估冻结科学结果。未经授权不推送或合并。
