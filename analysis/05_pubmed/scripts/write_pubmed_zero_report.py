import csv,json,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
B=Path(__file__).resolve().parents[1];ROOT=B.parents[1]
P=B/'pubmed_zero_recheck'
rows=list(csv.DictReader((P/'47_zero_hit_summary.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
raw=list(csv.DictReader((P/'47_uniform_PubMed_counts.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
assert len(rows)==47 and len(raw)==94
assert all(r['status']=='ok' for r in raw)
simple=[r['gene'] for r in rows if r['simple_symbol_count']=='0']
expanded=[r['gene'] for r in rows if r['expanded_names_count']=='0']
assert len(simple)==8 and len(expanded)==6
root=ET.parse(B/'sources/key_paper_metadata.xml')
paper=next(a for a in root.findall('.//PubmedArticle') if a.findtext('./MedlineCitation/PMID')=='36381921')
affiliations=[' '.join(a.itertext()) for a in paper.findall('.//Affiliation') if 'SMSR' in ' '.join(a.itertext())]
assert any('School of Medical Sciences and Research' in a for a in affiliations)
(P/'SAMD8_SMSr_false_hit_review.json').write_text(json.dumps({'PMID':'36381921','title':''.join(paper.find('.//ArticleTitle').itertext()),'SMSR_affiliations':affiliations,'finding':'SMSR is an institution abbreviation, not SAMD8/SMSr target research'},ensure_ascii=False,indent=2),encoding='utf8')
table='\n'.join('| '+r['gene']+' | '+r['simple_symbol_count']+' | '+r['expanded_names_count']+' | '+('机构缩写误命中，靶点相关0' if r['gene']=='SAMD8' else 'CASIMO1真实乳腺癌研究，剔除' if r['gene']=='SMIM22' else '两种检索均0')+' |' for r in rows if r['gene'] in simple)
links='\n'.join('- ['+r['gene']+'：符号检索]('+r['simple_PubMed_link']+')；[补全名称检索]('+r['expanded_PubMed_link']+')' for r in rows if r['gene'] in simple)
report='''# PubMed零命中标准：47个候选重新复核

日期：2026-10-24。按用户明确要求，研究空白的操作标准为**PubMed中乳腺癌与该靶点联合检索结果为0条**。不再把“有文献但未见专门机制”列入零命中组。

## 实际结果

47个候选均完成两种统一查询，共94次有效查询，无查询失败被记作0。

- 基因符号直接联合：8个0命中。
- 补充已核实别名、蛋白名称及乳腺癌同义词：6个字面0命中。
- SAMD8补名后有1条机构缩写误命中，靶点相关文献仍为0，单列说明；不混同为页面字面0。
- SMIM22补CASIMO1后有1篇真实乳腺癌研究，剔出零命中候选。

| 候选 | 符号＋breast cancer | 补全名称＋统一乳腺癌词组 | 核对结果 |
|---|---:|---:|---|
'''+table+'''

**补名后字面0的6个：MFSD10/SLC75A1、MT-ND4L、NHSL2、SLC18B1、TMEM254、TSPAN3。**

SAMD8需单列：[PMID36381921](https://pubmed.ncbi.nlm.nih.gov/36381921/)研究乳腺癌组织形态与ER/PR/HER2，作者单位SMSR是School of Medical Sciences and Research的简称，未研究SAMD8。其原符号查询为0；加入SMSr后字面计数为1，相关靶点文献为0。

SMIM22的命中是[PMID29765154](https://pubmed.ncbi.nlm.nih.gov/29765154/)，以CASIMO1为名，包含乳腺癌功能研究，是真实漏检，不能保留为零命中靶点。

此前提及MALL、ATP5MC1、MYADM、TMC4等“机制未明确”的候选，**不符合本次零命中标准**，不能因为未找到专门机制而纳入本组。

## 检索方法

第一种为每个符号相同的直接检索，例如：

```text
"SAMD8" AND "breast cancer"
```

第二种靶点端使用本次冻结UniProt条目的官方符号、基因别名、推荐及替代蛋白名称/短名，并补入此前已核实的历史别名。NHSL2不含错误的KIAA2022；SMIM22包含CASIMO1；SAMD8包含SMSr/FLJ25082；SLC18B1包含VPAT/C6orf192；MFSD10包含现名SLC75A1。

统一疾病端为：

```text
("breast neoplasms"[MeSH Terms]
 OR "breast cancer"[All Fields]
 OR "breast carcinoma"[All Fields]
 OR "breast tumor"[All Fields]
 OR "breast tumour"[All Fields]
 OR "mammary carcinoma"[All Fields]
 OR "mammary cancer"[All Fields]
 OR "mammary tumor"[All Fields]
 OR "mammary tumour"[All Fields]
 OR "TNBC"[All Fields])
```

不使用mammary或breast单独作为疾病词，以免把正常乳房、哺乳资料或鸡胸肌研究混入。NCBI实际querytranslation和短语未收录警告逐条保留，未收录别名没有被伪称已经有效检索。

完整记录：[47项计数汇总](pubmed_zero_recheck/47_zero_hit_summary.tsv)、[94条查询式、翻译与PMID](pubmed_zero_recheck/47_uniform_PubMed_counts.tsv)。这里统计PubMed索引查询，网页、GEO和论文补充表中的基因记录不直接计入PubMed命中数。

## 可点击复核

'''+links+'''

本次只核对零命中名单，不改既有疾病输入、预测和对接结果，也不据零命中自动指定最终作用靶点。
'''
(B/'PubMed零命中标准_复核结论.md').write_text(report,encoding='utf8')
notice='> 最新选题口径（2026-10-24）：用户明确采用PubMed“乳腺癌＋靶点”联合检索0命中标准。此前正文的“未见专门机制研究”分类不再作为零命中名单。当前计数、SAMD8/SMSr机构误命中及SMIM22/CASIMO1漏检，统一见[PubMed零命中复核](PubMed零命中标准_复核结论.md)。原始对接结果不变。\n\n'
for name in ['查新复核与候选结论.md','两候选配对对接_结果报告.md','README_当前入口.md']:
    p=B/name;t=p.read_text(encoding='utf8')
    if not t.startswith('> 最新选题口径'):p.write_text(notice+t,encoding='utf8')
header='## 25. 2026-10-24：用户明确PubMed零命中标准，重新统一检索（最新选题口径）'
hp=ROOT/'工作流3_交接.md'
if header not in hp.read_text(encoding='utf8'):
    with hp.open('a',encoding='utf8') as f:
        f.write('\n\n'+header+'\n\n用户质疑按研究空白缩窄范围，并进一步明确操作标准：PubMed综合乳腺癌与靶点联合检索结果须为0。必须按此口径，不得继续把“有文献但未见专门机制”的候选混入零命中组。此前第23/24节机制创新性表述保留历史，但不能当作当前零命中名单。\n\n本轮47候选×符号直接查询/补全UniProt已核实名称和历史别名的统一疾病词组，共94条有效PubMed ESearch查询。符号直接查询有8个0：MFSD10、MT-ND4L、NHSL2、SAMD8、SLC18B1、SMIM22、TMEM254、TSPAN3。补名后严格字面0为6：MFSD10/SLC75A1、MT-ND4L、NHSL2、SLC18B1、TMEM254、TSPAN3。SAMD8补SMSr后字面1，为PMID36381921作者单位SMSR（School of Medical Sciences and Research）误命中，靶点相关仍0，单列而不混同字面0。SMIM22补CASIMO1后命中PMID29765154真实乳腺癌研究，剔除。MALL/ATP5MC1/MYADM/TMC4等虽可能机制未充分研究，但PubMed非0，不纳入本组。\n\n入口outputs/docking_audit_2026-10-24/PubMed零命中标准_复核结论.md，pubmed_zero_recheck/47_zero_hit_summary.tsv及47_uniform_PubMed_counts.tsv含所有实际计数、PMID、可点击PubMed检索、查询翻译和警告。旧首轮仅基因名扩展表已标SUPERSEDED，正式扩展包括蛋白短名SMSr等。未使用mammary/breast单词独立疾病端，避免正常乳房/哺乳/鸡胸误命中。无查询失败当作0，无新疾病或对接计算，现有结果保留。\n')
ag=ROOT/'AGENTS.md';tag='## 2026-10-24最新：PubMed零命中操作标准'
if tag not in ag.read_text(encoding='utf8'):
    with ag.open('a',encoding='utf8') as f:
        f.write('\n\n'+tag+'\n\n用户明确“乳腺癌＋该靶点在PubMed综合检索0条”为当前研究空白标准；不能用未见专门机制代替0命中。94条统一复查后，补名字面0为MFSD10/SLC75A1、MT-ND4L、NHSL2、SLC18B1、TMEM254、TSPAN3。SAMD8符号0、补SMSr后1条为作者单位缩写误命中，靶点相关0，需单列；SMIM22/CASIMO1有真实乳腺癌论文应剔。新入口PubMed零命中标准_复核结论.md及pubmed_zero_recheck目录；详见交接25节。原预测/疾病/对接不改，不自动扩大计算。\n')
manifest=[]
for p in sorted(B.rglob('*')):
    if not p.is_file() or any(t in p.relative_to(B).parts for t in ['env','software','__pycache__']):continue
    if p.name=='Deliverable_SHA256_manifest.tsv':continue
    manifest.append({'relative_path':p.relative_to(B).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
with (B/'audit/Deliverable_SHA256_manifest.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
print('94/94 queries valid; symbol zeros8, expanded literal zeros6; SAMD8 false hit separately reviewed')
