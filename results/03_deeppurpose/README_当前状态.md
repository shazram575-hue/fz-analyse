# DeepPurpose：当前状态与恢复入口

更新：2026-10-24（北京时间）。**用户已接受CID209254作为明确计算结构并保留命名限制；本轮10模型探索性预测已完成。** 1,881主基因＋76诊断额外基因共1,993条参考蛋白序列，19,930条原始评分；27项数据/复算检查通过，未作预测准确性或实验验证，也未选定最终靶点。当前无本轮遗留任务。

最新入口：[结果报告](DeepPurpose_结果报告.md)、[主证据总表](results/Main1881_combined_disease_model_evidence.tsv)、[Kd排序](results/Main_Kd_consensus_ranked.tsv)、[IC50排序](results/Main_IC50_consensus_ranked.tsv)。两个终点共识Spearman0.269，前20交集3/20；保留分歧，不混合原始终点分数。以下准备记录用于追溯，完成状态以本段和结果报告为准。

## 1. 已同步的疾病输入及用户决定

单细胞基础补全已于2026-10-24验收；当前是6个生物学大类＋Mixed/Unresolved两个诊断组，非精细免疫分型。最新同类型DE为Tumor-origin−Healthy，与原恶性专属比较不同。目录称“芯片”的Bulk实际仍为TCGA STAR RNA-seq。

用户已明确选择：六个生物学大类的细胞级双向DE取并集，与既定Bulk双向探索DE相交并限定protein_coding，得到**1,881个主输入基因**。原单细胞门槛沿用本次补全结果：abs(avg_log2FC)>0.5、max(pct.1,pct.2)>0.15、全33,487特征Bonferroni<0.05；Bulk为abs(log2FC)>0.585且原始P<0.05，并保留BH。

Mixed/Unresolved纳入后的全集为1,957，其中76个是主输入之外的诊断组额外候选，仅作敏感性。患者pseudobulk保留为支持信息，不加为主输入硬门槛。Mast/Perivascular_SMC供者不足标记、上下调方向、来源类型均保留。原Pool A813及恶性专属比较未改，新1,881集合不冒称原Pool A。

正式入口：

- main_input/Disease_matching_candidates_1881.tsv：已确认范围的主基因表；不是已验证药物靶点。
- main_input/Disease_candidates_gene_celltype_evidence.tsv：3,599条基因×类型证据，保留各统计和方向。
- main_input/Additional_diagnostic_only_76genes.tsv：仅诊断组带来的额外候选。
- main_input/Main_gene_protein_mapping_review.tsv：主基因全部保留的蛋白映射审阅表，未映射项NA。
- main_input/Confirmed_scope_and_provenance.json：用户选择、定义、数量和文件SHA256。

candidate_scope_review的4种备选表仅为最初讨论材料；当前主范围以main_input确认文件为准。患者PB4类交集4,172、含诊断PB6组交集4,255均未被选为主输入。

## 2. 本地环境与10模型实际检查

独立环境：env/Scripts/python.exe；Python3.12.14、DeepPurpose0.1.5、PyTorch2.5.1+cpu、RDKit2024.09.6、numpy1.26.4、pandas2.2.3、scikit-learn1.5.2。完整版本见logs/pip_freeze.txt。未修改系统Python或R；此处CPU构建不等于确认本机没有GPU。与本次编码无关的DGL/调参分支依赖未安装，不宣称整个工具库全部功能均已验证。

已下载/审计10个作者预训练模型：CNN_CNN、Morgan_CNN、Morgan_AAC、Daylight_AAC、MPNN_CNN各自的BindingDB Kd和BindingDB IC50版本。model_archives保留ZIP及本地校验，models保留权重与安全读取的纯配置。一个IC50 ZIP含Mac资源副本，脚本按精确文件名排除__MACOSX，未改模型权重。

10个模型均严格加载state_dict，在作者公开toy_data的2个配对重复成4行后，通过有限值、重复运行、单项/批量和重复行一致性检查。**这是技术运行检查，不是预测准确性、泛化或外部验证。** public_validation/Public_example_scores.tsv仅为公开样例分数，绝不能误用为Spiropachysine A结果。

发现的真实问题：未改动的DeepPurpose0.1.5 MPNN在batch1与batch4时公共样例最大差异分别达约2.619个Kd对数单位和7.217个IC50对数单位。原始诊断、分数与失败日志保留；逐对batch_size1后重复一致。没有改权重或放宽断言。逐对计算不证明已恢复原模型准确性，因此这两个MPNN模型暂不直接并入正式共识，后续仍须核验；其余8个模型技术检查通过。初次QA预期4条作者样例，实际仅2条，已改成2条重复为4行并明示，不伪称4个独立配对。

后续组合原则：多编码可检查表示差异，不能保证消除数据库偏置。上述模型都来自BindingDB；Kd与IC50是不同终点，应分别记录和排名，不直接平均成结合概率。Davis/KIBA的激酶背景也不能自动推广到所有人类蛋白。最终共识规则/准确性评估尚未冻结，当前没有Top N或亲和力硬阈值。

## 3. 本地蛋白映射与数据保护（直接映射准备快照，后续已补齐）

首次外部UniProt ID提交被自动审批拒绝：候选是研究衍生名单，尚无明确授权发送给该外部服务。**被拒脚本未执行，没有候选名单外传；脚本已禁用。** 已采用安全替代：下载全部公开人类已审校UniProt数据，随后仅在本地映射。

公开版本UniProt2026_03（2026-09-02），20,431条记录，原始压缩JSON与SHA256保留。为兼容当时尚未选择的备选范围，先映射4,335基因并集：4,308个按Ensembl GeneId直接匹配，4,318个基因—蛋白条目，27个未匹配；10个基因对应多个已审校条目，全保留。

已确认的1,881主范围：1,869个基因有直接映射，共1,874个蛋白条目，12个基因未有直接已审校映射；5个基因多条目。223条蛋白超过CNN的1,000氨基酸上限；22条映射的当前主要符号与历史输入符号不完全一致，待结合稳定ID核对，不静默删改。未映射/多条目/符号不同均不是阴性药物结论。primary_symbol救援候选另表保留，没有自动替代正式ID映射。

还已下载GENCODE v36蛋白翻译参考，MD5与第一阶段已保存的官方MD5清单一致：20e8a4eb0f54be8e7e85b293021da08a。尚未执行该参考的额外序列救援，不把“已下载”写成“映射问题已解决”。

## 4. 结构证据与早期等待安排（已由后续用户接受CID决定覆盖）

已读取目录中的14页Fang2022完整正文，定位来源/提取鉴定、试剂、Figure1A与参考文献；不是只读摘要。源PDF原件未改，全文提取物在structure_sources/Fang2022_user_fulltext.txt。正文明确将NMR鉴定参照指向Chang等1998研究（DOI10.1021/np980162x），另引Chiu等1990的Isospiropachysine论文（DOI10.1016/0031-9422(90)85362-J）。不能仅凭近似名称将不同立体记录当同一分子。

作者公开Figure1A已取得并查看：连接关系可辨，但未画出完整立体标记。公开mmc1.pdf已实际读取，内容为语言润色证明，不是NMR或结构鉴定补充；其余试探的mmc2/3/4请求为404。2025晶体文献的下载响应为空/非PDF，不算取得全文。1975原始结构研究全文已取得，已视觉核读首页与结构图相关页；未把它自动等同于2022实验样品。

PubChem209254及500127原始属性JSON已取得。二者连接图相同；500127记录缺少立体指定，不能把其名称当作已提供异构体的完整立体结构。对同图的立体/去立体/翻转测试，DeepPurpose Morgan及Daylight特征完全相同，CNN字符特征不同。此处只做特征核查，**未计算研究蛋白的药物分数**。

用户在得知可先用6个指纹模型探索的选项后，明确选择：**先等其给出全部结构核准，再进行。** 因此不以无立体图、PubChem近似名称或任一推测构型启动正式或探索性匹配。

用户需要准备的优先材料和验收口径见[结构核准材料清单](结构核准_需要准备的材料.md)。2022正文和1975全文已在本地，不要重复索取；下一份最有针对性的文献是正文指定的Chang1998鉴定依据，或直接与目标样品一致的带立体结构文件/原作者鉴定材料。

## 5. 验收、恢复与未完成项

audit/Preparation_validation.tsv中7项均通过：4种集合从各类型显著表独立重建、7份第一阶段冻结文件不变、10模型技术运行检查及公开下载无候选名单。最终药物结果没有完成标记。

scripts入口：audit_candidate_universes→assemble_confirmed_input；download_public_human_proteome→map_proteins_locally→review_local_mapping；prepare_model_panel→validate_model_runtime；fetch_structure_sources/fetch_identity_records/check_stereochemistry_sensitivity均为结构资料与特征准备。map_protein_sequences.py已禁用，不得重试候选名单外发。

运行时曾由库默认参数创建根目录空result文件夹，经确认空目录后已移入本线outputs；加载器现显式把内部输出定位到public_validation/model_internal_results，不再向根目录写库结果。

以上为准备期恢复记录。后续用户已接受具体CID并允许继续：参考序列已补齐，10模型推理和结果汇总/一致性验收已完成；实验准确性尚未验证。原813不改，不干预其他会话。需要恢复时优先阅读最新结果报告，不重复执行已完成主体。

## 2026-10-24后续：用户提供结构已完成CID一致性核验

用户提供的SMILES仅需去掉括号前的Markdown反斜杠转义，即与官方CID209254异构SMILES逐字相同；9个手性中心全指定，InChIKey比对及MOL立体往返通过。文件与核验见structure_verification_CID209254/README_核验结论.md。该CID记录名为Spiropachysine，当前同义名未列Spiropachysine A；研究样品与数据库结构的对应尚未自动宣布完成。本轮未开始药物预测，第4节的等待安排保留；不应再把用户提供的CID结构说成“没有立体信息”。
