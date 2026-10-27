# GitHub文件政策

2026-10-26整理阶段仅准备清单，未初始化或推送。2026-10-27用户授权安全同步并指定公开仓库[shazram575-hue/fz-analyse](https://github.com/shazram575-hue/fz-analyse)；本轮以真实Git ignore和逐文件预检为准，实际状态见[同步日志](TASK_LOGS/2026-10-27_github_sync.md)。

## 适合受控上传的内容

- README、当前AGENTS、docs科研进度/方法/审计摘要、config及新工程工具。
- analysis中逐字节代码副本、results中2 MiB以内的主要表/冻结程序/报告/PNG/SVG，均有PATH_MAP及哈希。
- environments中的实际包名/版本、R DESCRIPTION清单、核心requirements。它们不是Linux已验证锁文件。
- data/TCGA_DOWNLOAD_MANIFEST.tsv与原始文件校验/获取说明；仅包含公共文件UUID及校验元数据，不包含大型RNA计数。

旧候选清单留作2026-10-26快照；当前逐文件大小/SHA在[本轮候选清单](audit/2026-10-27/github-sync/UPLOAD_CANDIDATES.tsv)。公共去标识化Case ID不是原始身份索引；本轮核对作者公开来源、直接身份字段和原始reads限制。未确认可公开的新增临床数据不得自动沿用本轮许可。

## 不适合普通Git的内容

- 原始数据8.05 GB、完整outputs、大型RDS/H5AD/MTX、全量全蛋白/DE表、NMF随机重复及大缓存。
- Windows虚拟环境、site-packages、DLL/EXE/JAR、本机R包库和本地工具状态。
- `.env`、私钥/证书、服务账号凭证；真实凭证用云平台secret机制，不提交仓库。
- PDF全文、原作者XLSX、模型权重/配置pickle及归档：默认本地保留或外置，许可/平台兼容未全面核实，不自动重分发。
- 审计时完整83k文件盘点before/after及带本机绝对路径的local_snapshots默认忽略；可上传精简科学索引和摘要。

`.gitignore`默认忽略原始数据和outputs，故不能只写ignore而把所有代码/结果漏掉；analysis/results正是可审核的独立副本。`.gitattributes`对这些副本关闭EOL转换以保持科学证据SHA；新工程脚本用LF。

GitHub对普通Git文件的硬限制为100 MiB，超过50 MiB也会提示；不要用Git充当整个35.48 GB目录的备份。LFS/Release/对象存储各有不同限制，未在本轮自动配置。[GitHub官方说明](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)

## 上传前仍需的动作

首次同步按本轮日志记录；用户已经授权本次安全提交/推送，不需重复索取常规确认。大型派生对象位置及Linux环境/访问测试仍待解决；每次新增上传仍须复查实际索引、敏感字段和来源。本轮没有改写或清理既有远程历史。
