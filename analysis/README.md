# 可审阅代码副本

各阶段scripts目录是旧outputs代码的逐字节副本，未改变参数或统计行为。源路径、目标路径及SHA256见[路径映射](../docs/PATH_MAP.tsv)。原件仍在outputs。

这些副本不是已移植的Linux流水线。许多旧脚本依赖outputs下的相对位置、Windows R库或本机解释器路径。**不要直接批量运行本目录。**先读[科研进度](../docs/RESEARCH_PROGRESS.md)及[迁移说明](../docs/CLOUD_MIGRATION.md)。`scripts/project/restore_legacy_layout.py`可在新克隆中按原路径恢复小文件副本，但不能替代大型数据、模型权重或环境重建。
