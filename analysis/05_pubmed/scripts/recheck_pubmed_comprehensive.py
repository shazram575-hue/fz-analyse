"""Expanded, layered PubMed audit; preserve former queries and identify new hits."""
from pathlib import Path
import csv,json,time,urllib.parse,xml.etree.ElementTree as ET,re
import requests
B=Path(__file__).resolve().parents[1]
O=B/'pubmed_comprehensive_recheck';O.mkdir(exist_ok=True)
old=list(csv.DictReader((B/'pubmed_zero_recheck/47_uniform_PubMed_counts.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
old={r['gene']:r for r in old if r['variant']=='expanded_all_confirmed_names'}
assert len(old)==47
mesh='("Breast Neoplasms"[Mesh] OR "Mammary Neoplasms, Animal"[Mesh] OR "Mammary Neoplasms, Experimental"[Mesh])'
site='(breast*[tiab] OR mammary[tiab])'
tumour='(cancer*[tiab] OR neoplas*[tiab] OR carcinom*[tiab] OR tumor*[tiab] OR tumour*[tiab] OR malignan*[tiab] OR adenocarcinom*[tiab])'
histology='("ductal carcinoma in situ"[tiab] OR "lobular carcinoma in situ"[tiab] OR "invasive ductal carcinoma"[tiab] OR "invasive lobular carcinoma"[tiab] OR "inflammatory breast cancer"[tiab] OR "Paget disease of the nipple"[tiab] OR "Paget disease of the breast"[tiab] OR ((DCIS[tiab] OR LCIS[tiab] OR IDC[tiab] OR ILC[tiab]) AND '+site+'))'
subtypes='(TNBC[tiab] OR "triple negative breast"[tiab] OR "triple-negative breast"[tiab] OR ('+site+' AND ("basal-like"[tiab] OR "luminal A"[tiab] OR "luminal B"[tiab] OR "HER2-positive"[tiab] OR "HER2-negative"[tiab] OR "HER2-enriched"[tiab] OR "ER-positive"[tiab] OR "ER-negative"[tiab] OR "hormone receptor-positive"[tiab])))'
models=['MCF-7','MCF7','MDA-MB-231','MDAMB231','MDA-MB-468','MDAMB468','MDA-MB-453','MDAMB453','MDA-MB-436','MDA-MB-157','T-47D','T47D','BT-474','BT474','BT-549','BT549','BT-20','BT20','SK-BR-3','SKBR3','ZR-75-1','ZR75-1','HCC1937','HCC1806','HCC70','HCC38','HCC1954','HCC1569','HCC1428','HCC1500','Hs578T','Hs-578T','SUM149','SUM149PT','SUM159','SUM159PT','SUM1315','SUM1315MO2','SUM190','SUM225','4T1','EMT6','E0771','EO771','MMTV-PyMT','MMTV-neu','MMTV-Wnt1']
model_query='('+' OR '.join('"'+t+'"[tiab]' for t in models)+')'
cohorts='("TCGA-BRCA"[tiab] OR METABRIC[tiab] OR "SCAN-B"[tiab] OR ('+site+' AND ("patient-derived xenograft"[tiab] OR "patient-derived organoid"[tiab])))'
extra='('+mesh+' OR ('+site+' AND '+tumour+') OR '+histology+' OR '+subtypes+' OR '+model_query+' OR '+cohorts+')'
strategy={'MeSH':mesh,'site_plus_neoplasm_stems':site+' AND '+tumour,'histology':histology,'subtypes':subtypes,'cell_and_animal_models':model_query,'cohorts':cohorts,'notes':'Previous disease query retained as OR branch. No species/date/language restrictions. Cancer cell models included; MCF10A not used alone as cancer indicator. HER2/ER/IDC/ILC subtype terms constrained to breast context. No unqualified BRCA/BC.'}
(O/'Disease_search_blocks.json').write_text(json.dumps(strategy,ensure_ascii=False,indent=2),encoding='utf8')
s=requests.Session();s.headers['User-Agent']='PubMedBreastDomainAudit/2.0'
def retrieve(query,path,db='pubmed',xml=False):
    if path.exists():return path.read_bytes() if xml else json.loads(path.read_text(encoding='utf8'))
    params={'db':db,'retmode':'xml' if xml else 'json'}
    params.update({'id':query} if xml else {'term':query,'retmax':1000})
    endpoint='efetch.fcgi' if xml else 'esearch.fcgi'
    for attempt in range(4):
        try:
            url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+endpoint
            r=s.get(url,params=params,timeout=60)
            if r.status_code in [414,431]:r=s.post(url,data=params,timeout=60)
            r.raise_for_status()
            if xml:
                ET.fromstring(r.content);path.write_bytes(r.content);return r.content
            obj=r.json()
            if 'esearchresult' not in obj:raise ValueError(str(obj)[:200])
            path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf8');return obj
        except Exception:
            if attempt==3:raise
            time.sleep(2*(attempt+1))
rows=[];new_ids=set()
for gene,prior in sorted(old.items()):
    target,old_disease=prior['query'].split(' AND ',1)
    query=target+' AND ('+old_disease+' OR '+extra+')'
    try:
        obj=retrieve(query,O/f'{gene}_expanded_domain.json');d=obj['esearchresult']
        if d.get('ERROR'):raise ValueError(d['ERROR'])
        count=int(d['count']);ids=set(d.get('idlist',[]));prev=set(filter(None,prior['PMIDs'].split(';')))
        assert count>=int(prior['count']),(gene,prior['count'],count)
        added=ids-prev
        if int(prior['count'])<=1:new_ids.update(added|ids)
        row={'gene':gene,'previous_count':prior['count'],'expanded_domain_count':count,'status':'ok','new_PMIDs':';'.join(sorted(added,key=int)),'PMIDs':';'.join(sorted(ids,key=int)),'query':query,'query_translation':d.get('querytranslation',''),'warnings':json.dumps({k:d[k] for k in ['warninglist','errorlist'] if k in d},ensure_ascii=False),'PubMed_link':'https://pubmed.ncbi.nlm.nih.gov/?term='+urllib.parse.quote(query)}
        print(gene,prior['count'],'->',count,('NEW '+row['new_PMIDs']) if int(prior['count'])<=1 else '',flush=True)
    except Exception as e:
        row={'gene':gene,'previous_count':prior['count'],'expanded_domain_count':'','status':'FAILED_NOT_ZERO','new_PMIDs':'','PMIDs':'','query':query,'query_translation':'','warnings':str(e),'PubMed_link':'https://pubmed.ncbi.nlm.nih.gov/?term='+urllib.parse.quote(query)}
        print('FAILED',gene,str(e),flush=True)
    rows.append(row)
    with (O/'47_comprehensive_PubMed_counts.tsv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    time.sleep(.4)
ids=sorted(new_ids,key=int)
for i in range(0,len(ids),30):
    retrieve(','.join(ids[i:i+30]),O/f'focused_papers_{i:03}.xml',xml=True);time.sleep(.4)
papers=[]
for p in sorted(O.glob('focused_papers_*.xml')):
    for a in ET.parse(p).findall('.//PubmedArticle'):
        txt=lambda q:' '.join(''.join(x.itertext()) for x in a.findall(q))
        pid=txt('./MedlineCitation/PMID');aids={x.attrib['IdType']:x.text for x in a.findall('./PubmedData/ArticleIdList/ArticleId')}
        papers.append({'PMID':pid,'title':txt('.//ArticleTitle'),'abstract':txt('.//AbstractText'),'authors_affiliations':txt('.//Author'),'MeSH':txt('.//MeshHeading'),'PMCID':aids.get('pmc',''),'doi':aids.get('doi','')})
if papers:
    with (O/'focused_paper_metadata.tsv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(papers[0]),delimiter='\t');w.writeheader();w.writerows(papers)
print('ZERO',[r['gene'] for r in rows if r['expanded_domain_count']==0 and r['status']=='ok'],flush=True)
print('FOCUS',[(r['gene'],r['previous_count'],r['expanded_domain_count'],r['new_PMIDs']) for r in rows if int(r['previous_count'])<=1],flush=True)
print('COMPLETE',len(rows),len(papers),flush=True)
