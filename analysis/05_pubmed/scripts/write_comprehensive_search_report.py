import csv,json,hashlib
from pathlib import Path
B=Path(__file__).resolve().parents[1];ROOT=B.parents[1];O=B/'pubmed_comprehensive_recheck'
rows=list(csv.DictReader((O/'47_comprehensive_PubMed_counts.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
assert len(rows)==47 and all(r['status']=='ok' for r in rows)
zero=[r['gene'] for r in rows if r['expanded_domain_count']=='0']
assert zero==['MFSD10','NHSL2','SLC18B1','TMEM254','TSPAN3'],zero
papers=list(csv.DictReader((O/'focused_paper_metadata.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
nd=next(r for r in papers if r['PMID']=='24061460')
assert 'malignant breast tissues' in nd['abstract'] and 'ND4L' in nd['abstract']
focus=['MFSD10','MT-ND4L','NHSL2','SAMD8','SLC18B1','SMIM22','TMEM254','TSPAN3']
table='\n'.join('| '+r['gene']+' | '+r['previous_count']+' | '+r['expanded_domain_count']+' | '+('新增真实乳腺癌记录，剔除' if r['gene']=='MT-ND4L' else 'CASIMO1真实乳腺癌记录，剔除' if r['gene']=='SMIM22' else 'SMSR机构缩写误命中，相关0' if r['gene']=='SAMD8' else '当前扩充式仍返回0')+' |' for r in rows if r['gene'] in focus)
report='''# PubMed乳腺癌领域扩充检索：最终复核结果

执行日期：2026-10-24。用户要求立即用完整、正确的疾病领域词组重查。**本轮已对47个候选全部执行扩充检索，并核查旧零命中组的新增记录。**

## 实际结论

**严格返回0条的候选更新为5个：MFSD10/SLC75A1、NHSL2、SLC18B1、TMEM254、TSPAN3。**

SAMD8继续单列：基因符号直接查询为0；补SMSr后原始返回1条，为已核实的机构缩写误命中，靶点相关0。SMIM22/CASIMO1有真实论文，不纳入零命中组。

| 靶点 | 上版扩名计数 | 本轮完整领域词组计数 | 结论 |
|---|---:|---:|---|
'''+table+'''

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
'''
(B/'PubMed零命中标准_复核结论.md').write_text(report,encoding='utf8')
(O/'完整检索复核报告.md').write_text(report.replace('](pubmed_comprehensive_recheck/',']('),encoding='utf8')
full=[]
for r in rows:
    full.append('## '+r['gene']+'\n\n```text\n'+r['query']+'\n```\n\n返回数：'+r['expanded_domain_count']+'；上一版：'+r['previous_count']+'\n\n新增PMID：'+(r['new_PMIDs'] or '无')+'\n\n[在PubMed复核]('+r['PubMed_link']+')\n')
(O/'47项完整检索式.md').write_text('# 47项已执行的PubMed完整检索式\n\n2026-10-24；来源为实际ESearch请求。\n\n'+'\n'.join(full),encoding='utf8')
header='## 26. 2026-10-24：完整乳腺癌领域词组已重查，零命中更新为5个（覆盖25节旧计数）'
hp=ROOT/'工作流3_交接.md'
if header not in hp.read_text(encoding='utf8'):
    with hp.open('a',encoding='utf8') as f:
        f.write('\n\n'+header+'\n\n用户质疑乳腺癌领域词覆盖并明确要求直接重查。本轮对47个候选全部执行扩充PubMed式：保留旧式，增加人/动物/实验乳腺肿瘤MeSH，breast*/mammary与cancer*/neoplas*/carcinom*/tumor*/tumour*/malignan*等词干组合，病理/亚型，常用乳腺癌细胞及动物模型变体，METABRIC/TCGA-BRCA/SCAN-B等。跨癌种亚型词有乳腺背景限制，不用BRCA/BC或正常MCF10A独立代表癌症。47/47查询有效，记录querytranslation/警告，全部计数不低于旧式。\n\n新增真实漏检：MT-ND4L由0→1，PMID24061460，Chintha2013 Mitochondria and tumors: a new perspective（10.4103/0019-509X.118732），正式摘要明示恶性乳腺组织ND4L突变。旧式只含固定乳腺癌词组且该文MeSH为泛Neoplasms；部位词与肿瘤词干组合找回，故必须剔除MT-ND4L。\n\n当前严格返回0为5：MFSD10/SLC75A1、NHSL2、SLC18B1、TMEM254、TSPAN3。SAMD8扩充前后均原始1，为PMID36381921作者单位SMSR缩写误命中，相关靶点0，继续单列；SMIM22仍有CASIMO1真实论文，剔除。旧零/近零组5篇实际返回文献的元数据/摘要已保存并核读相应条目。\n\n现行入口outputs/docking_audit_2026-10-24/PubMed零命中标准_复核结论.md已重写；pubmed_comprehensive_recheck/47_comprehensive_PubMed_counts.tsv、Disease_search_blocks.json、47项完整检索式.md、47份原始JSON及focused_paper_metadata.tsv可复核。第25节6个0是旧检索式快照，被本节覆盖。未改疾病/预测/对接结果，无本轮后台任务。\n')
ag=ROOT/'AGENTS.md';tag='## 2026-10-24完整领域词复核最新计数'
if tag not in ag.read_text(encoding='utf8'):
    with ag.open('a',encoding='utf8') as f:
        f.write('\n\n'+tag+'\n\n已按用户要求完成47候选全面PubMed扩充词重查（MeSH＋乳腺/肿瘤词干＋病理/亚型＋癌细胞/动物模型＋队列词），47/47有效。新增MT-ND4L真实乳腺组织研究PMID24061460，故当前严格0为5：MFSD10/SLC75A1、NHSL2、SLC18B1、TMEM254、TSPAN3。SAMD8原始1仍为SMSR机构缩写误命中，相关0单列；SMIM22有CASIMO1论文排除。此计数覆盖此前6个0；入口PubMed零命中标准_复核结论.md及pubmed_comprehensive_recheck，交接26节。原疾病/预测/对接不改，无需重跑。\n')
checks={'date':'2026-10-24 Asia/Shanghai','n_valid_queries':47,'n_failures':0,'old_query_preserved_as_OR_branch':True,'all_counts_monotonic':all(int(r['expanded_domain_count'])>=int(r['previous_count']) for r in rows),'strict_zero_genes':zero,'new_true_breast_hit':{'gene':'MT-ND4L','PMID':'24061460'},'SAMD8_false_hit_only':'36381921','new_docking_or_MD':False}
(O/'COMPREHENSIVE_SEARCH_COMPLETED.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf8')
manifest=[]
for p in sorted(B.rglob('*')):
    if not p.is_file() or any(t in p.relative_to(B).parts for t in ['env','software','__pycache__']):continue
    if p.name=='Deliverable_SHA256_manifest.tsv':continue
    manifest.append({'relative_path':p.relative_to(B).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
with (B/'audit/Deliverable_SHA256_manifest.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
print('Complete: 47 valid queries;5 literal zeros;MT-ND4L true new hit excluded;reports and handover updated')
