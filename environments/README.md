# 实测环境，不是跨平台锁

本轮从三个实际Python环境采集包元数据：programs、deeppurpose、docking的installed.tsv/runtime.json。完整pip freeze和含本机绝对路径的原始输出放local_snapshots（Git忽略）。requirements.txt仅列programs环境核心实测版本，既不是完整transitive lock，也未在Linux重建验证；不将三个阶段依赖混装。

R运行时为4.5.3。受限会话installed.packages返回空，因此R-installed-observed.tsv采用两个实际R库的DESCRIPTION回退，1028条记录含重复包及不同Built版本。优先库为D:/R4.5.2/R-4.5.2/library；其中Matrix1.7-5与另一库Matrix1.7-4并存。infercnv1.20.0的Built为R4.4.0/Windows。该清单证明文件/版本存在，不证明全部包在当前会话可加载；既往分析加载版本见原日志。

没有生成猜测的renv.lock或environment.yml，也没有修改系统R/Python配置。Linux需对应系统库、包源和ABI复核。元数据回退可用scripts/project/capture_R_descriptions.py复现；包枚举空表不得覆盖有效回退表。
