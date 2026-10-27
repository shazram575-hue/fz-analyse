import csv,json,statistics,re
from pathlib import Path
import gemmi
B=Path(__file__).resolve().parents[1];ROOT=B.parents[1]
def tab(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
pool=tab(ROOT/'outputs/身份深核与独立交集_2026-10-24/results/disease_intersections/top10percent_Kd_intersection_IC50_disease1881.tsv')
back=tab(ROOT/'outputs/身份深核与独立交集_2026-10-24/background/All_reviewed_human_protein_background.tsv')
rows=[]
for prec in sorted(pool,key=lambda x:x['gene_symbol']):
    g=prec['gene_symbol']
    # Map by frozen accession, not a current primary symbol (MFSD10 -> SLC75A1).
    accessions=prec['Kd_reviewed_proteins'].split(';')
    seqrec=max((x for x in back if x['uniprot_accession'] in accessions),key=lambda x:int(x['sequence_length']))
    acc=seqrec['uniprot_accession'];u=json.loads((B/'sources'/f'uniprot_{acc}.json').read_text(encoding='utf8'))
    pdb=[x for x in u.get('uniProtKBCrossReferences',[]) if x['database']=='PDB']
    feats=u.get('features',[]);af=B/'structures'/f'{acc}_AF.pdb'
    pd=[];match='';pockets=[]
    confidence={}
    if af.exists():
        chain=gemmi.read_structure(str(af))[0][0]
        residues=[r for r in chain if r.het_flag=='A']
        sequence=''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in residues)
        match=sequence==seqrec['protein_sequence']
        confidence={str(r.seqid.num):r['CA'][0].b_iso for r in residues}
        pd=list(confidence.values())
        file=B/('pockets/AF_MFSD10' if g=='MFSD10' else 'pockets/AF47')/f'{acc}_AF.pdb_predictions.csv'
        if file.exists():pockets=[{k.strip():v.strip() for k,v in x.items()} for x in csv.DictReader(file.open())]
    best=pockets[0] if pockets else {}
    pocket_pd=[confidence[r.split('_')[-1]] for r in best.get('residue_ids','').split() if r.split('_')[-1] in confidence]
    note=[]
    if g=='NHSL2':note+=['long protein; whole-protein docking not justified; local confidence and interdomain PAE required']
    if g=='SAMD8':note+=['experimental hexamer; resolved 144-407 only; inspect membrane reaction chamber']
    if g=='SLC18B1':note+=['experimental GFP/nanobody fusion; canonical-only extraction verified']
    if g=='MT-ND4':note+=['459-aa P03905 assessed; additional frozen 99-aa C0HME5 protein has no AF entry; experimental complex I entries exist']
    if g=='MFSD10':note+=['current primary symbol SLC75A1; accession-based mapping required; 6S4M has citrate/lipid, not validated drug ligand']
    if sum(x['type']=='Transmembrane' for x in feats):note+=['membrane protein: membrane-facing grooves are not automatically ligand sites']
    rows.append({'gene':g,'accession':acc,'frozen_length':len(seqrec['protein_sequence']),
      'current_UniProt_sequence_matches_frozen':u['sequence']['value']==seqrec['protein_sequence'],
      'experimental_PDB_count':len(pdb),'experimental_PDB_entries':';'.join(x['id'] for x in pdb),
      'PDB_ranges_method_resolution':json.dumps(pdb),'transmembrane_segments':sum(x['type']=='Transmembrane' for x in feats),
      'AF_downloaded':af.exists(),'AF_sequence_matches_frozen':match,'AF_mean_pLDDT':statistics.mean(pd) if pd else '',
      'AF_fraction_residues_pLDDT_ge70':sum(p>=70 for p in pd)/len(pd) if pd else '',
      'AF_fraction_residues_pLDDT_lt50':sum(p<50 for p in pd)/len(pd) if pd else '',
      'P2Rank_pocket_count':len(pockets),'AF_top_pocket_probability':best.get('probability',''),
      'AF_top_pocket_residue_mean_pLDDT':statistics.mean(pocket_pd) if pocket_pd else '',
      'AF_top_pocket_fraction_pLDDT_ge70':sum(p>=70 for p in pocket_pd)/len(pocket_pd) if pocket_pd else '',
      'AF_top_pocket_residues':best.get('residue_ids',''),
      'assessment':'experimental_structure_requires_site_assembly_review' if pdb else ('predicted_pocket_requires_local_review' if pockets else 'no_detected_AF_pocket_not_negative_binding_evidence'),
      'notes':'; '.join(note)})
with (B/'results/47_structure_pocket_assessment.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('ROWS',len(rows),'PDB proteins',sum(x['experimental_PDB_count']>0 for x in rows),'AF',sum(x['AF_downloaded'] for x in rows),'pocket proteins',sum(x['P2Rank_pocket_count']>0 for x in rows))
print(json.dumps([x for x in rows if x['gene'] in ['NHSL2','SAMD8','SLC18B1','TMEM254','MFSD10','SMIM22','MT-ND4']],ensure_ascii=False,indent=2))
