# GitHub同步报告

更新：2026-10-27。**代码、研究文档和已审查小型结果已成功同步至指定公开仓库main，远程读取及完整文件树哈希核验通过。可以开始创建Codex Cloud环境；大型数据及Linux分析运行环境仍需另行准备。**

| 项目 | 实际状态 |
|---|---|
| 用户指定仓库 | [shazram575-hue/fz-analyse](https://github.com/shazram575-hue/fz-analyse)，public，owner及push权限已核实 |
| 本地分支/上游 | main → origin/main |
| 已核实的远程代码导入提交 | `9342813d8f45cee7228e5a1fa58d33ad76b6c046`，831个已检查文件；后续收尾文档提交通过git log核实 |
| 起始远程提交 | `58f54b6a1f3859086893f5f524580fc799de8a70`，仅占位README；保留父历史 |
| 已同步分支 | GitHub main；本地main跟踪origin/main，内容一致 |
| 同步方式 | 本机Git CLI push因无登录凭据失败；随后通过已登录GitHub连接器上传Git对象、创建提交并非强制更新main，最后用Git fetch独立验证成功。不将失败的push命令写成成功 |
| 具体清单 | [PATIENT_FIELD_REVIEW.tsv](audit/2026-10-27/github-sync/PATIENT_FIELD_REVIEW.tsv)，列出全部133路径、字段、行数及来源说明，不包含临床值 |
| 当前下一步 | 在Codex Cloud选择该仓库main，先做工程证据核验和资源/数据缺口检查，不重跑冻结分析 |

已完成：逐文件路径/体积/哈希和凭据模式检查；829暂存路径与批准清单一致；728份科学副本在Git blob中逐项SHA256一致；无原始/大对象/密钥暂存，11项ignore用例和6项预检安全测试通过。图像30项人工标记已检查为绘图软件/分子元数据及内部嵌图。以上不是绝对隐私或许可证保证。

可供审阅的主要内容是README、AGENTS、科研进度、数据清单、迁移说明、任务日志、R/Python/Shell脚本、环境版本清单、冻结程序和小表图。原始数据约8.05 GB、完整outputs、RDS/H5AD/reads、权重、Windows环境、作者PDF/XLSX均未纳入提交，仍只在本地。

Python3.12.14的UTF-8模式下doctor和verify通过；Windows默认GBK首次失败已记录。Linux整套依赖、R/结构软件、绝对路径、网络和大型对象外部存储未验证。不能将本地工程通过写为Cloud运行成功。

已通过本地独立Git clone的无大型数据测试：clone不含原始目录及大型对象，doctor和728份副本校验仍通过；这是Windows本地Git传输测试，不是Linux或GitHub上传测试。

此前尝试私有化时电脑操作因URL识别不足停止，没有修改仓库可见性。随后用户明确指示“不需要设置为隐私的，里面都是公共数据的公共分析，你就直接上传就好了”，更新为公开上传授权。重新检查owner、public可见性、push权限及远程历史后，继续原先已审查文件的普通推送；不通过其他工具绕过审批，不上传受控reads或原始大型数据。

代码导入的远程树为`f0ce64f1be1cb1777c58eb48f28be51fe10ff264`，与本地已审查提交`bad8bb1bbbd8453c4861bf6de43d8bca26c4a08e`的完整Git树完全一致。独立fetch确认远程main，远程接口按提交读取README、AGENTS和科研进度成功；另确认PROJECT_AUDIT、DATA_MANIFEST和CLOUD_MIGRATION均在树中。另存[机器可读同步证据](audit/2026-10-27/github-sync/SYNC_VERIFIED.json)。

连接器分批大文件请求曾超时，已从保存断点降低并发重试；130个唯一较大blob及10批树对象核验后才更新main。原远程初始提交保留为父；原本地三份草稿提交保存在`archive/local-drafts-before-github-api-2026-10-27`，没有删除历史或强推。API导入提交因创建者/时间元数据不同而SHA不同，文件树保持相同。

本次文件范围及公开上传由用户明确批准；没有新建科研PR、合并科研结论或自动升级`MERGED_REVIEWED`。后续文档收尾提交通过git log定位，不要求本报告填写包含自身最终文本的commit SHA。Cloud尚未实际创建、未运行Linux完整分析；本机以后使用Git CLI写入仍需正常登录认证。
