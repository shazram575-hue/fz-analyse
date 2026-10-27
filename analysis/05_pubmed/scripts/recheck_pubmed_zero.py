"""Uniform PubMed zero-hit audit. Read public metadata only; no file uploads."""
import csv,json,time,sys,re,urllib.parse,xml.etree.ElementTree as ET
from pathlib import Path
import requests
B=Path(__file__).resolve().parents[1];ROOT=B.parents[1]
OUT=B/'pubmed_zero_recheck';OUT.mkdir(exist_ok=True)
def table(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
pool=table(ROOT/'outputs/身份深核与独立交集_2026-10-24/results/disease_intersections/top10percent_Kd_intersection_IC50_disease1881.tsv')
back=table(ROOT/'outputs/身份深核与独立交集_2026-10-24/background/All_reviewed_human_protein_background.tsv')
disease='("breast neoplasms"[MeSH Terms] OR "breast cancer"[All Fields] OR "breast carcinoma"[All Fields] OR "breast tumor"[All Fields] OR "breast tumour"[All Fields] OR "mammary carcinoma"[All Fields] OR "mammary cancer"[All Fields] OR "mammary tumor"[All Fields] OR "mammary tumour"[All Fields] OR "TNBC"[All Fields])'
s=requests.Session();s.headers['User-Agent']='PubMedZeroHitResearchAudit/1.0'
rows=[]
for item in sorted(pool,key=lambda r:r['gene_symbol']):
    gene=item['gene_symbol']
    accs=item['Kd_reviewed_proteins'].split(';')
    entry=max((r for r in back if r['uniprot_accession'] in accs),key=lambda r:int(r['sequence_length']))
    u=json.loads((B/'sources'/f"uniprot_{entry['uniprot_accession']}.json").read_text(encoding='utf8'))
    aliases=[gene]
    for g in u.get('genes',[]):
        aliases += [g['geneName']['value']] if 'geneName' in g else []
        aliases += [x['value'] for x in g.get('synonyms',[])]
    aliases=sorted(set(aliases),key=lambda x:(x!=gene,x))
    additional={'SAMD8':['SMSr','FLJ25082'],'SLC18B1':['VPAT','C6orf192'],
      'MFSD10':['SLC75A1','TETRAN','TETTRAN'],'TMEM254':['C10orf57','FLJ13263','bA369J21.6'],
      'TSPAN3':['TM4SF8','TM4-A','TSP-3'],'NHSL2':['NHS like 2','NHS-like 2'],'SMIM22':['CASIMO1']}
    aliases=list(dict.fromkeys(aliases+additional.get(gene,[])))
    protein_name=u.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value','')
    terms=list(aliases)
    if protein_name:terms.append(protein_name)
    description=u.get('proteinDescription',{})
    for named in [description.get('recommendedName',{})]+description.get('alternativeNames',[]):
        if named.get('fullName'):terms.append(named['fullName']['value'])
        terms += [x['value'] for x in named.get('shortNames',[])]
    terms=list(dict.fromkeys(terms))
    variants={
      'simple_symbol':f'"{gene}" AND "breast cancer"',
      'expanded_all_confirmed_names':'('+' OR '.join('"'+t+'"[All Fields]' for t in terms)+') AND '+disease}
    for variant,query in variants.items():
        path=OUT/f'{gene}_{variant}.json'
        try:
            if path.exists():obj=json.loads(path.read_text(encoding='utf8'))
            else:
                for attempt in range(4):
                    try:
                        response=s.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',params={'db':'pubmed','retmode':'json','term':query,'retmax':1000},timeout=50)
                        response.raise_for_status();obj=response.json()
                        if 'esearchresult' not in obj:raise ValueError(str(obj)[:200])
                        path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf8');break
                    except Exception:
                        if attempt==3:raise
                        time.sleep(2*(attempt+1))
            d=obj['esearchresult'];count=int(d['count'])
            status='ok'
            if 'ERROR' in d or 'error' in obj:raise ValueError(str(d))
            translation=d.get('querytranslation','')
            if not any(t.lower() in translation.lower() for t in terms):
                status='check_target_terms_in_translation'
            row={'gene':gene,'variant':variant,'count':count,'status':status,'PMIDs':';'.join(d.get('idlist',[])),
             'confirmed_gene_aliases':';'.join(aliases),'protein_name':protein_name,'query':query,'query_translation':translation,
             'warnings_errors':json.dumps({k:d[k] for k in ['warninglist','errorlist'] if k in d},ensure_ascii=False),
             'PubMed_link':'https://pubmed.ncbi.nlm.nih.gov/?term='+urllib.parse.quote(query),'date':'2026-10-24 Asia/Shanghai'}
            print(gene,variant,count,status,flush=True)
        except Exception as exc:
            row={'gene':gene,'variant':variant,'count':'','status':'query_failed_not_zero','PMIDs':'','confirmed_gene_aliases':';'.join(aliases),'protein_name':protein_name,'query':query,'query_translation':'','warnings_errors':str(exc),'PubMed_link':'https://pubmed.ncbi.nlm.nih.gov/?term='+urllib.parse.quote(query),'date':'2026-10-24 Asia/Shanghai'}
            print('FAILED',gene,variant,exc,flush=True)
        rows.append(row)
        with (OUT/'47_uniform_PubMed_counts.tsv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        time.sleep(.38)
summary=[]
for gene in sorted({r['gene'] for r in rows}):
    rr={r['variant']:r for r in rows if r['gene']==gene};a=rr['simple_symbol'];b=rr['expanded_all_confirmed_names']
    summary.append({'gene':gene,'simple_symbol_count':a['count'],'expanded_names_count':b['count'],'strict_expanded_zero':b['count']==0 and b['status']=='ok',
      'simple_status':a['status'],'expanded_status':b['status'],'expanded_PMIDs':b['PMIDs'],'confirmed_aliases':b['confirmed_gene_aliases'],
      'simple_PubMed_link':a['PubMed_link'],'expanded_PubMed_link':b['PubMed_link']})
with (OUT/'47_zero_hit_summary.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
print('ZERO_SIMPLE',[r['gene'] for r in summary if r['simple_symbol_count']==0],flush=True)
print('ZERO_EXPANDED',[r['gene'] for r in summary if r['strict_expanded_zero']],flush=True)
print('COMPLETE',len(summary),flush=True)
