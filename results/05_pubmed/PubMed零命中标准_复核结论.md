# PubMed乳腺癌领域扩充检索：最终复核结果

执行日期：2026-10-24。用户要求立即用完整、正确的疾病领域词组重查。**本轮已对47个候选全部执行扩充检索，并核查旧零命中组的新增记录。**

## 实际结论

**严格返回0条的候选更新为5个：MFSD10/SLC75A1、NHSL2、SLC18B1、TMEM254、TSPAN3。**

SAMD8继续单列：基因符号直接查询为0；补SMSr后原始返回1条，为已核实的机构缩写误命中，靶点相关0。SMIM22/CASIMO1有真实论文，不纳入零命中组。

| 靶点 | 上版扩名计数 | 本轮完整领域词组计数 | 结论 |
|---|---:|---:|---|
| MFSD10 | 0 | 0 | 当前扩充式仍返回0 |
| MT-ND4L | 0 | 1 | 新增真实乳腺癌记录，剔除 |
| NHSL2 | 0 | 0 | 当前扩充式仍返回0 |
| SAMD8 | 1 | 1 | SMSR机构缩写误命中，相关0 |
| SLC18B1 | 0 | 0 | 当前扩充式仍返回0 |
| SMIM22 | 1 | 1 | CASIMO1真实乳腺癌记录，剔除 |
| TMEM254 | 0 | 0 | 当前扩充式仍返回0 |
| TSPAN3 | 0 | 0 | 当前扩充式仍返回0 |

**与上一版相比，MT-ND4L必须剔除。** 新增文献：Chintha等，*Mitochondria and tumors: a new perspective*，Indian J Cancer，2013；[PMID24061460](https://pubmed.ncbi.nlm.nih.gov/24061460/)，DOI10.4103/0019-509X.118732。正式摘要明确说明样本包含乳腺良恶性组织，并报告ND4L、ND6、COX-II变异只见于恶性乳腺组织。这是实质性漏检，已经纠正；不能再称MT-ND4L为乳腺癌PubMed零命中靶点。

这篇摘要使用“breast”和“malignant tissues”等分散表达，没有依赖“breast cancer”固定词组；MeSH也为较笼统的Neoplasms。新加入的“乳腺部位词 AND 肿瘤词干”找回了这条记录。

SAMD8的1条为[PMID36381921](https://pubmed.ncbi.nlm.nih.gov/36381921/)，SMSR出现在作者单位School of Medical Sciences and Research，未研究SAMD8；扩充疾病词后没有新增其他记录。SMIM22的真实命中为[PMID29765154](https://pubmed.ncbi.nlm.nih.gov/29765154/)，研究对象别名CASIMO1。

## 本轮检索具体补了什么

靶点端沿用47项各自已核实的符号、现名、历史别名、蛋白全名与短名。疾病端在上一版全部词项之外，增加以下模块，模块之间用OR合并，再与靶点端AND：

1. **MeSH及下位分类：**Breast Neoplasms、Mammary Neoplasms Animal、Mammary Neoplasms Experimental。未关闭自动下位词扩展。
2. **部位与肿瘤词干：**`(breast* OR mammary) AND (cancer* OR neoplas* OR carcinom* OR tumor* OR tumour* OR malignan* OR adenocarcinom*)`，检索题名/摘要。该模块不要求“breast cancer”紧邻出现。
3. **病理表述：**导管/小叶原位癌、浸润性导管/小叶癌、炎性乳腺癌、乳头Paget病，以及带乳腺背景限制的DCIS/LCIS/IDC/ILC。
4. **分型：**TNBC、triple-negative/triple negative breast，以及带乳腺背景限制的basal-like、luminal A/B、HER2阳性/阴性/富集、ER阳性/阴性、激素受体阳性。
5. **常用乳腺癌模型及名称变体：**MCF-7/MCF7、MDA-MB-231/468/453/436/157、T47D、BT474、BT549、BT20、SKBR3、ZR75、多个HCC及SUM系列、Hs578T；4T1、EMT6、E0771/EO771、MMTV-PyMT/neu/Wnt1等动物模型。
6. **队列及模型背景：**TCGA-BRCA、METABRIC、SCAN-B，以及带乳腺背景的PDX/类器官。

HER2、ER、IDC、ILC等跨领域缩写均有乳腺背景约束；不把BRCA或BC单独当作乳腺癌，不把正常乳腺细胞MCF10A单独当作癌症模型。无时间、语言或物种限制。

上述模块依据[NLM PubMed检索规则](https://pubmed.ncbi.nlm.nih.gov/help/#mesh-terms)、[乳腺肿瘤MeSH分类](https://www.ncbi.nlm.nih.gov/mesh/D001943)及[实验性乳腺肿瘤MeSH定义](https://www.ncbi.nlm.nih.gov/mesh/68008325)整理。每个候选实际执行的完整字符串、NCBI翻译和警告均已保存，不能只凭模块标题代替执行记录。

## 完整可复核文件

- [47项扩充前后实际计数、PMID、完整检索式和PubMed链接](pubmed_comprehensive_recheck/47_comprehensive_PubMed_counts.tsv)
- [疾病端全部检索模块](pubmed_comprehensive_recheck/Disease_search_blocks.json)
- [旧零命中及近零组命中文献的正式元数据与摘要](pubmed_comprehensive_recheck/focused_paper_metadata.tsv)
- `pubmed_comprehensive_recheck/*_expanded_domain.json`：47份NCBI原始ESearch响应。

47次查询全部有效；每项均保留旧式作为OR分支，计数不应下降，本轮全部满足。旧零命中组新增记录已逐条核对：仅MT-ND4L增加PMID24061460，为真实乳腺组织研究；SAMD8原机构误命中及SMIM22原真实命中性质不变。这里报告PubMed索引检索结果；旧版“6个严格0”的结论由本版“5个严格0＋SAMD8误命中单列”覆盖。

原预测、疾病分析、对接分数和姿势均未修改。没有启动新的对接或MD。
