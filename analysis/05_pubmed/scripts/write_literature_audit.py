import csv,json,hashlib
from pathlib import Path
B=Path(__file__).resolve().parents[1];ROOT=B.parents[1]
def read(p,sep=','):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter=sep))
original=read(ROOT/'附件1_靶点检索结果汇总.csv')
papers=read(B/'results/Representative_paper_metadata_audit.tsv','\t')
special={
 'AGPAT2':('已有表达/药物反应/关联证据','20363836研究AGPAT11/LPCAT2，不应归为AGPAT2的癌组织上调；34756920和33380715仍支持非空白','34756920;33380715'),
 'APH1A':('已有患者组织检测与预后证据','不只是生信模型：24223915有55份乳腺癌组织qPCR；41463075为405患者Notch比较','24223915;41463075'),
 'ATP5MC1':('已有组学记录，专门机制未确认','转录组记录不等于基因本身功能验证','28701793;39667070'),
 'GPR160':('已有临床药物遗传及功能后续线索','34048027不只是预后模型：SNP×AI比较并描述EMT/生长通路后续；未完整核读其功能实验，不作无功能断言','34048027'),
 'MALL':('已有组学/预后记录，原归属错误','MALL/BENE与MAL、MAL2是不同基因；原6代表文献均不能支撑所写MALL专门机制。2018转录组有MALL，2023七基因模型也含MALL','10.1038/s41598-018-26429-x;10.1155/2023/6242355'),
 'MFSD10':('已有乳腺癌组学记录，专门机制未检出','现名SLC75A1；原0命中不等于无研究。TNBC miRNA研究表中有MFSD10，当前未见独立乳腺癌扰动/救援机制','41877801;6S4M'),
 'MFSD12':('已有明确乳腺癌机制，必须纠正','2023 Molecular Cell：MFSD12/mTORC1-T254、CTNS、氧化应激耐受、化疗反应与突变体功能；不能写无专门功能/临床研究','37751742'),
 'MT-ND2':('已有遗传/表达及功能关联','不能把所有ND2命中视为线粒体ND2，保留原文归属；原表中有乳腺癌遗传研究','31577800;24414975'),
 'MT-ND4':('已有乳腺癌机制','let-7a靶向线粒体ND4及代谢重编程已研究','34838007'),
 'MT-ND4L':('已有乳腺癌表达/变异记录','2018原发乳腺癌线粒体变异和2022幸存者疲劳RNA研究都含ND4L；不是完全空白','10.3390/cancers10120500;PMC9420746'),
 'MYADM':('已有组学/家族综述线索，专门机制未确认','不得从家族其他成员功能移用于MYADM；需要区分家族综述与本基因乳腺癌扰动','37345137;41454416'),
 'NHSL2':('已有乳腺癌预后/组学记录，专门机制未检出','原26870154为KIAA2022/NEXMIF，与NHSL2无关。2023 TME研究Table2列NHSL2；另有学位论文siRNA筛查线索，仅检索段落未取得有效PDF，不写为完整核读','36936167;10.3389/fendo.2023.1131525'),
 'ORMDL2':('已有乳腺癌患者预后证据','2018研究已分析ORMDL1/2表达与患者生存，2025综述不是唯一证据；不把ORMDL3机制移给ORMDL2','PMC5929433;40757047'),
 'PNP':('已有酶前药研究，需区分人源与细菌PNP','38039708是症状缩写PNP，不是基因。酶前药/自杀基因中的PNP多为外源微生物酶，不证明人内源PNP是同一乳腺癌机制靶点','24098491;33708463'),
 'PIGC':('已有乳腺癌分子轴研究','两篇corrigendum不是两项新的独立研究；保留原NEAT1/MALAT1轴研究但不提前宣称PIGC直接药物靶点','35378140;33171372'),
 'SAMD8':('已有乳腺癌组学记录，专门机制未检出','别名SMSr/FLJ25082必须补检。FEN1乳腺癌补充表实际有SAMD8；SH003表网页也列出该基因，但完整正文获取失败，不声称全文核读','本目录sources/FEN1_supplement.pdf;https://cgp.iiarjournals.org/content/20/1/88'),
 'SLC18B1':('已有乳腺癌表达记录，专门机制未检出','别名VPAT/C6orf192必须补检；GEO Profiles 50593035为乳腺癌细胞MTDH敲低表达数据。可优先探索机制，不写无任何乳腺癌数据','GEO50593035;10.1038/s41467-025-59549-w'),
 'SMIM22':('已有明确乳腺癌机制，撤销空白','CASIMO1为其别名：2018 Oncogene已有乳腺癌敲低、过表达、SQLE相互作用及救援','29765154'),
 'TMC4':('已有肿瘤表达/暴露关联，专门机制未确认','原代表研究是环境暴露及乳腺癌死亡关联，不能概括为仅线粒体组学','31593922;32341500'),
 'TMEM254':('已有乳腺癌补充组学记录，专门机制未检出','C10orf57/FLJ13263/bA369J21.6为别名；FEN1补充材料PDF第53页含TMEM254。TMEM254-AS1是不同转录本，不能当本基因机制','本目录sources/FEN1_supplement.pdf'),
 'TSPAN13':('已有明确乳腺癌机制，必须纠正','2005过表达与增殖侵袭、2019 miRNA-TSPAN13、2025在线/2026卷期免疫逃逸研究均已存在','15958618;30701690;41218148'),
 'TSPAN3':('已有泛癌/BRCA分析，专门机制未检出','2022泛癌TSPAN研究覆盖BRCA。2018 AnAc检索摘录的TSPAN33不是TSPAN3，已排除该误读','PMC9385710'),
 'VKORC1L1':('已有乳腺癌细胞检测与患者关联','2019 Oncotarget有多TNBC细胞系qPCR及患者数据；不能写仅泛癌间接提及。未证实VKORC1L1单基因因果机制','31040920'),
}
mechanism=set('ACKR1 CA12 CLDN7 CYBRD1 EMP1 F3 FXYD3 GPR132 GPRC5A IFI27 ITPR3 KCNK1 LRRC26 MAL2 NR3C1 PMP22 PPARG PTGER4 S1PR1 SCD SSTR2'.split())
out=[]
for row in original:
    gene=row['基因符号'];refs=';'.join(p['PMID'] for p in papers if p['gene']==gene)
    level,note,refs=special.get(gene,('已有乳腺癌功能/机制研究' if gene in mechanism else '已有乳腺癌临床/表达关联','非空白；代表文献题名摘要已核，未把所有列出正文标为全文核读',refs))
    out.append({'gene':gene,'submitted_level':row['证据等级'],'submitted_raw_count':row['乳腺癌相关文献数'],'submitted_query_reproduced':True,'revised_assessment':level,'fully_blank_confirmed':False,'evidence_and_correction':note,'source_identifiers':refs,'scope':'literature/structure prioritization; not confirmed drug target'})
with (B/'results/47_literature_evidence_revised.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
flags={('AGPAT2','20363836'):'AGPAT11/LPCAT2 differs from AGPAT2',('NHSL2','26870154'):'KIAA2022/NEXMIF differs from NHSL2',
 ('PNP','38039708'):'PNP = psychoneurophysiological symptoms, not the gene',('SCD','22245140'):'SCD = BRCA1 serine-containing domain, not stearoyl-CoA desaturase',
 ('F3','35146636'):'F3 = plant sub-fraction, not tissue factor',('F3','7648165'):'pancreatic cancer study; no breast-cancer-specific mechanism conclusion',
 ('PMP22','35178836'):'PERP is distinct from PMP22',('PMP22','41456178'):'PERP is distinct from PMP22',('PMP22','28925397'):'PERP is distinct from PMP22',
 ('CLDN7','21051248'):'canine mammary study; label species',('SSTR2','25985817'):'canine mammary study; label species',
 ('PIGC','37783583'):'corrigendum, not independent original study',('PIGC','37821336'):'corrigendum, not independent original study'}
for p in papers:
    p['attribution_flag']=('MAL/MAL2 or compound-name literature, not MALL' if p['gene']=='MALL' else flags.get((p['gene'],p['PMID']),'no clear name/species mismatch detected at title/abstract screen; not a full-text causal validation'))
with (B/'results/167_representative_paper_attribution_audit.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(papers[0]),delimiter='\t');w.writeheader();w.writerows(papers)
hashes=[]
for name in ['附件1_靶点检索结果汇总.csv','附件2_代表性文献清单.csv','乳腺癌47靶点_查新检索记录.html']:
    p=ROOT/name;hashes.append({'file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(B/'audit/user_supplied_file_hashes.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2),encoding='utf8')
report='''# 查新复核与候选采用结论

2026-10-24。已读用户HTML及两份CSV、工作流交接21节及相关最新入口。47原查询命中数逐一复现；167条代表记录的PMID与题名全部吻合。**这只能证明检索可复现和文献确实存在，不能证明基因归属、证据级别或“无研究”结论正确。原附件原样保留。**

## 本轮确定的采用意见

1. **优先计算SAMD8、SLC18B1**，使用“在本轮检索范围内，尚未检出以该基因为核心的乳腺癌专门机制研究”表述。二者已有乳腺癌表达/组学记录，禁止写“完全空白”“从未研究”。结构对接只提供假设。
2. **第三个优先补核候选建议MFSD10，现名SLC75A1**。其同一UniProt Q14728、原疾病/预测入池证据保留；有6S4M人源2.4Å实验结构，AF主口袋概率0.991且局部可信度高。已有TNBC miRNA研究中的组学列举，但本轮未检出乳腺癌单基因机制研究。此建议基于创新性与结构可评估性，不把它指定为已确认结合靶点，也不因建议而删除其余候选。
3. **NHSL2暂缓本批结构对接**。原指给它的PMID26870154实际上是KIAA2022/NEXMIF；2023原始乳腺癌TME研究又确实在Table2列出NHSL2（HR1.48、95%CI1.06–2.06、原始P0.0214），属于探索性组学预后关联。学位论文siRNA筛查有检索段落，下载收到HTML挑战页，未当作有效PDF/已读全文。其AF全长平均pLDDT47.89，74.78%残基低于50，P2Rank未检出位点；不适合把全长盲对接硬排在另两者旁边。这不是阴性结合证据。
4. **未确认任何一个候选具有“全领域/所有乳腺癌资料完全无记录”的绝对空白。** 应按专门机制、功能实验、临床/组学、数据库表达分别描述。原来用0–1条归为“空白”的标准撤销，不能唯一决定研究创新性。

## 影响选择的实质纠错

| 原记录 | 复核结论 | 原始证据 |
|---|---|---|
| SMIM22完全空白 | 错误；CASIMO1已有乳腺癌敲低、过表达和救援研究 | [Polycarpou-Schwarz2018](https://doi.org/10.1038/s41388-018-0281-5)，PMID29765154；原文摘要核读，非完整付费正文 |
| MFSD12无专门功能/临床研究 | 错误；mTORC1-T254调控、突变体抑瘤、患者化疗反应已研究 | [He2023](https://doi.org/10.1016/j.molcel.2023.08.032)，PMID37751742；官方摘要及公开全文段落 |
| TSPAN13机制缺乏 | 错误；2005与2019已有功能，2025在线/2026卷期有MHC-I降解免疫机制 | [2005](https://pubmed.ncbi.nlm.nih.gov/15958618/)、[2019](https://pubmed.ncbi.nlm.nih.gov/30701690/)、[2026](https://doi.org/10.1158/0008-5472.CAN-25-1223)；摘要核读、最新开放正文段落核读 |
| MALL=MAL | 错误；MALL/BENE、MAL、MAL2为不同基因，Pt-Mal也是化合物命名 | 原6条代表项归属均需剔出其专门机制依据；[MALL真实转录组记录](https://doi.org/10.1038/s41598-018-26429-x) |
| NHSL2=KIAA2022 | 错误；KIAA2022是NEXMIF，2015突变文献不支持NHSL2 | 原2015正文核读；[NEXMIF官方信息](https://www.ncbi.nlm.nih.gov/gene/340533)；[NHSL2真实乳腺癌预后](https://doi.org/10.3389/fendo.2023.1131525)全文Table2核读 |
| ORMDL2仅一综述 | 错误；2018已分析乳腺癌患者生存 | [原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC5929433/)全文核读相关段落 |
| VKORC1L1仅泛癌间接提及 | 错误；有TNBC细胞qPCR和乳腺癌患者数据，尚不等于单基因因果 | [2019原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC6481349/)全文核读相关方法结果 |
| MT-ND4L空白 | 错误；已有患者mtRNA/变异记录 | [2018原始研究](https://doi.org/10.3390/cancers10120500)、[2022乳腺癌疲劳研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC9420746/)公开正文/表格检索段落 |
| SAMD8空白 | 太强；乳腺癌FEN1补充表含SAMD8，SH003药物处理表也列出 | 本地FEN1补充PDF第46/65页已提取并渲染核验；[SH003](https://cgp.iiarjournals.org/content/20/1/88)仅网页检索表，403未取得完整正文 |
| SLC18B1空白 | 太强；已有乳腺癌细胞基因表达数据 | [NCBI原始GEO profile50593035](https://www.ncbi.nlm.nih.gov/geoprofiles/50593035)已保存；不能据此推断SLC18B1本身功能 |
| MFSD10空白 | 太强；现名补检必需，TNBC miRNA论文Table含MFSD10 | [原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13008373/)BioC全文表已核；[SLC75A1官方命名](https://www.ncbi.nlm.nih.gov/gene/10227) |
| TMEM254空白 | 太强；FEN1补充表第53页有记录 | 本地补充PDF已提取并渲染核验；AS1不是本基因 |
| TSPAN3空白 | 太强；泛癌研究已纳入BRCA | [原始泛癌分析](https://pmc.ncbi.nlm.nih.gov/articles/PMC9385710/)；2018文中的TSPAN33不能归给TSPAN3 |

另外，F3被混入植物分级F3，PNP被混入症状缩写PNP，PMP22被混入另一基因PERP，SCD被混入BRCA1蛋白片段缩写SCD；AGPAT2被借用AGPAT11/LPCAT2的癌组织上调。犬乳腺肿瘤研究单列物种。PIGC的两篇勘误不计作两项独立原始研究。详细标记见167条归属审计表。

“2020年至检索日”在2026年不能称严格近5年，应改“2020年以来”；原文按靶点累加1892条也不是所有靶点合并去重后的独立论文总数。PubMed AllFields仍不等于全文检索，0命中尤其不等于无补充表记录；短别名、作者同名、字段不存在与被忽略术语均可产生假阳性/假阴性。查询JSON保留了NCBI的errorlist/warnings与实际querytranslation，未经人工判读的扩展计数不作为纳入文献数。

## 可审阅输出与边界

- results/47_literature_evidence_revised.tsv：47项逐行修订，旧标签保留对照。
- results/PubMed_search_reproduction.tsv：原式与关键扩展式、PMID、翻译和警告；不当作人工相关性筛后的文献计数。
- results/Representative_paper_metadata_audit.tsv：167题名/PMID核对；元数据核对不代表全文审读。
- results/167_representative_paper_attribution_audit.tsv：名称/物种/勘误归属风险。
- audit/user_supplied_file_hashes.json：新增用户文件哈希，原件未修改。

本次为PubMed复现＋开放原始资料及公开网页补检，不是跨所有商业数据库、专利库和未发表资料的法定查新。未声称“绝对没有”某项研究。结构/表达/预测均不证明本药直接结合或因果作用。已获用户回复：基因敲低/过表达、CETSA/DARTS和纯化蛋白结合检测均可开展；具体细胞系未指定，本轮计算不以其为前置条件。
'''
# Do not claim an unverified DOI for SH003 (retained only by known title/URL).
report=report.replace(';10.21873/cgp.20359','')
report+='\n## 用户本轮最后收窄（优先）\n用户明确配对对接只需解释已核实候选与药物的相互作用可行性，不要过度扩展。因此本轮收尾范围为SAMD8、SLC18B1；不要求47统一正式对接，不为凑三个候选启动MFSD10。前述MFSD10是备用建议，不是新主靶点；NHSL2保留而暂缓结构对接。\n'
(B/'audit/初次查新复核_技术记录.md').write_text(report,encoding='utf8')
print('47 revised rows; 167 metadata/attribution rows; report written')
