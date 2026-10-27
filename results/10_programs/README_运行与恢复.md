# MFSD10程序发现与独立验证：本轮运行入口

用户已授权连续完成GSE176078导入/独立验证对象、GSE161529正式逐患者cNMF、共享程序与规则冻结、独立外部验证。此前第35节“仅可行性/不自动启动”已被本轮新授权更新。两个队列不得联合整合或共同发现程序。

本轮2026-10-26已完成计算、程序与规则冻结、外部统计及图形。最终交付完成必须以audit/PIPELINE_COMPLETED.json为准；audit/PROGRESS.json保留阶段恢复状态，不能只据脚本存在或启动日志判断完成。

优先读[结果报告](结果报告_程序发现与独立验证.md)、[冻结程序](frozen/Program_genes_v1.tsv)、[冻结规则](frozen/Analysis_freeze_v1.json)和[发现/外部完整比较](validation/Discovery_external_comparison_all.tsv)。已完成13患者3900次主分解及19场景380次局部敏感性；5程序外部基因协同性得到支持，但发现MFSD10关联均未过BH，无确认性关联成功。外部MP05弱正相关仅探索性，未倒改发现。

主验证对象为objects/GSE176078_malignant_candidates_author_counts_clinical_v2.h5ad；全细胞参考为objects/GSE176078_all_cells_author_counts_clinical_v2.h5ad。26病例均保留对应表，20病例有作者Cancer细胞，其中16Naïve主要验证、4Treated单列；其余6病例的该最终矩阵没有作者Cancer细胞，并非否定其肿瘤恶性。原v1导入对象保留，主统计使用clinical_v2。

原始下载包实际位置：原始数据/单细胞/GSE176078，含作者tar.gz、补表xlsx和方法PDF；用户给出的原始数据/GSE176078不存在。仅向本目录cache解压，原文件保留。

已核查此前只有可行性/表达摸底，无成功cNMF拟合文件，不复用其他工作流的程序、脚本、阈值或交接包。完整未抽样图谱仍只有原发现集17c的9486个确认恶性，主拟合使用其中8722既有singlet；未定/未评分不补入。

可恢复脚本与输出均在本目录。所有患者、K与随机种子的完整输出逐项保存，可断点跳过完成项。MFSD10从所有拟合输入、归一化分母及末端重拟合排除；其原始表达单列。程序筛选先于MFSD10关联，验证集只在冻结后评分。
