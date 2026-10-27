"""Read-only public-source audit; never uploads local files or sequences."""
import csv, json, re, time, hashlib, sys, xml.etree.ElementTree as ET
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1]
for name in ['sources','results','structures','logs','audit','software']:
    (OUT/name).mkdir(exist_ok=True, parents=True)
S = requests.Session()
S.headers['User-Agent'] = 'PublicResearchEvidenceAudit/1.0'
def get(url, path, params=None, binary=False):
    path = OUT/path
    if path.exists():
        return path.read_bytes() if binary else json.loads(path.read_text(encoding='utf-8'))
    for attempt in range(4):
        try:
            r=S.get(url, params=params, timeout=90); r.raise_for_status()
            if binary:
                path.write_bytes(r.content); return r.content
            obj=r.json(); path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8'); return obj
        except Exception:
            if attempt==3: raise
            time.sleep(2*(attempt+1))

def table(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig'), delimiter='\t' if path.suffix=='.tsv' else ','))
def write_table(name, rows):
    if not rows:return
    with (OUT/'results'/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def literature():
    rows=table(ROOT/'附件1_靶点检索结果汇总.csv')
    breast='("breast neoplasms"[MeSH Terms] OR "breast cancer"[Title/Abstract] OR "breast carcinoma"[Title/Abstract] OR "breast tumor"[Title/Abstract] OR "breast tumour"[Title/Abstract] OR mammary[Title/Abstract] OR TNBC[Title/Abstract])'
    corrected={'NHSL2':['NHSL2','NHS like 2','NHS-like 2'], 'MALL':['MALL','BENE','MAL-like'],
      'SAMD8':['SAMD8','SMSr','FLJ25082','sphingomyelin synthase related','ceramide phosphoethanolamine synthase'],
      'SLC18B1':['SLC18B1','VPAT','vesicular polyamine transporter','C6orf192'],
      'TSPAN3':['TSPAN3','TM4SF8','TM4-A','TSP-3','tetraspanin 3'],
      'MFSD10':['MFSD10','TETRAN'],'SMIM22':['SMIM22','CASIMO1'],
      'TMEM254':['TMEM254','C10orf57'], 'MT-ND4L':['MT-ND4L','MTND4L','ND4L'],
      'ORMDL2':['ORMDL2'], 'VKORC1L1':['VKORC1L1','VKORC1-like 1']}
    summary=[]
    for row in rows:
        gene=row['基因符号']
        variants={'submitted':row['检索式']}
        if gene in corrected:
            terms=corrected[gene]
            prefix='('+' OR '.join('"'+t+'"[Title/Abstract]' for t in terms)+')'
            variants['corrected_tiab']=prefix+' AND '+breast
            variants['corrected_allfields']='('+' OR '.join('"'+t+'"[All Fields]' for t in terms)+') AND '+breast.replace('[Title/Abstract]','[All Fields]')
        for variant,query in variants.items():
            data=get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',Path('sources')/f'pubmed_{gene}_{variant}.json',{'db':'pubmed','term':query,'retmode':'json','retmax':1000})
            d=data.get('esearchresult',{})
            summary.append({'gene':gene,'variant':variant,'user_count':row['乳腺癌相关文献数'],'current_count':d.get('count','ERROR'),'PMIDs':';'.join(d.get('idlist',[])),'query':query,'translation':d.get('querytranslation',''),'warnings':json.dumps(d.get('warninglist',{}),ensure_ascii=False)})
            write_table('PubMed_search_reproduction.tsv',summary)
            print(gene,variant,d.get('count'),flush=True);time.sleep(.4)
    papers=table(ROOT/'附件2_代表性文献清单.csv')
    ids=sorted(set(p['PMID'] for p in papers))
    metadata={}
    for n in range(0,len(ids),25):
        chunk=ids[n:n+25]
        data=get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',Path('sources')/f'pubmed_papers_{n:03}.xml',{'db':'pubmed','id':','.join(chunk),'retmode':'xml'},binary=True)
        for p in ET.fromstring(data).findall('.//PubmedArticle'):
            text=lambda q: ' '.join(''.join(x.itertext()) for x in p.findall(q))
            pid=text('./MedlineCitation/PMID')
            ids2={x.attrib.get('IdType'):x.text for x in p.findall('./PubmedData/ArticleIdList/ArticleId')}
            metadata[pid]={'title':text('.//ArticleTitle'),'pubYear':text('.//JournalIssue/PubDate/Year'),'doi':ids2.get('doi',''),'pmcid':ids2.get('pmc',''),'abstractText':text('.//AbstractText')}
        time.sleep(.4)
    audited=[]
    norm=lambda s:re.sub(r'\W+','',s.lower())
    for p in papers:
        m=metadata.get(p['PMID'],{})
        audited.append({'gene':p['基因符号'],'PMID':p['PMID'],'found':bool(m),'title_matches':norm(p['标题'])==norm(m.get('title','')),'submitted_title':p['标题'],'official_title':m.get('title',''),'pubyear':m.get('pubYear',''),'doi':m.get('doi',''),'PMCID':m.get('pmcid',''),'abstract':re.sub('<[^>]+>',' ',m.get('abstractText',''))})
    write_table('Representative_paper_metadata_audit.tsv',audited)

def structures():
    pool=table(ROOT/'outputs/身份深核与独立交集_2026-10-24/results/disease_intersections/top10percent_Kd_intersection_IC50_disease1881.tsv')
    print('POOL FIELDS',list(pool[0]),flush=True)
    background=table(ROOT/'outputs/身份深核与独立交集_2026-10-24/background/All_reviewed_human_protein_background.tsv')
    genes={r['gene_symbol'] for r in pool}
    selected=[r for r in background if any(g in genes for g in r['primary_gene_names'].split(';'))]
    results=[]
    for r in selected:
        acc=r['uniprot_accession']; gene=r['primary_gene_names']
        try:
            u=get(f'https://rest.uniprot.org/uniprotkb/{acc}.json',Path('sources')/f'uniprot_{acc}.json')
            seq=u['sequence']['value']
            pdb=[x for x in u.get('uniProtKBCrossReferences',[]) if x['database']=='PDB']
            af=get(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}',Path('sources')/f'alphafold_{acc}.json')
            if af and isinstance(af,list):
                for record in af[:1]:
                    get(record['pdbUrl'],Path('structures')/f'{acc}_AF.pdb',binary=True)
                    if gene in ['SAMD8','SLC18B1','NHSL2']:
                        get(record['paeDocUrl'],Path('structures')/f'{acc}_PAE.json',binary=True)
            results.append({'gene':gene,'accession':acc,'length':len(seq),'sequence_matches_frozen':seq==r['protein_sequence'],'PDB_entries':';'.join(x['id'] for x in pdb),'PDB_properties':json.dumps(pdb),'AF_available':bool(af),'status':'metadata_downloaded_not_pocket_validated'})
            print(gene,acc,'PDB',len(pdb),'AF',bool(af),flush=True)
        except Exception as e:
            results.append({'gene':gene,'accession':acc,'length':r['sequence_length'],'sequence_matches_frozen':'','PDB_entries':'','PDB_properties':'','AF_available':'','status':str(e)})
            print('ERROR',gene,str(e),flush=True)
        write_table('47_structure_inventory.tsv',results)
    for pdb in ['9D7W','9D7X','9D7V','9D7U']:
        get(f'https://files.rcsb.org/download/{pdb}.cif',Path('structures')/f'{pdb}.cif',binary=True)
        get(f'https://files.rcsb.org/download/{pdb}.pdb',Path('structures')/f'{pdb}.pdb',binary=True)
    print('STRUCTURES_DONE',flush=True)

if __name__=='__main__':
    {'literature':literature,'structures':structures}[sys.argv[1]]()
