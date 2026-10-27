"""Download public coordinate records and inspect canonical coverage before receptor preparation."""
import csv, json, hashlib, shutil, time, gzip, xml.etree.ElementTree as ET
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests, gemmi
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[3]; B=Path(__file__).resolve().parents[1]
OLD=ROOT/'outputs/docking_audit_2026-10-24'
SELECT={
 'ACKR1':('8YNF','full-length receptor complex rather than N-terminal peptide'),
 'ADRA2A':('7EJ8','ligand-bound receptor complex without T4-lysozyme insertion'),
 'APH1A':('8KCS','high-resolution gamma-secretase complex'),
 'ATP5MC1':('8H9S','high-resolution ATP synthase with c-ring copies'),
 'CA12':('1JD0','well-characterized inhibitor-bound catalytic domain'),
 'CYBRD1':('5ZLG','human membrane redox protein experimental structure'),
 'F3':('1BOY','human extracellular domain without assigning factor-VII ligand to F3'),
 'GPR132':('8HQM','small-molecule-bound receptor complex'),
 'ITPR3':('8TKG','high-resolution full-channel structure'),
 'KCNK1':('3UKM','available human channel structure'),
 'LRRC26':('7YO3','auxiliary-subunit complex with broader LRRC26 coverage than 8VAZ'),
 'MFSD10':('6S4M','available human SLC75A1 experimental structure'),
 'MT-ND2':('9I4I','high-resolution human respiratory complex I'),
 'MT-ND4':('9I4I','high-resolution human respiratory complex I'),
 'MT-ND4L':('9I4I','high-resolution human respiratory complex I'),
 'NR3C1':('4UDD','high-resolution small-molecule-bound ligand-binding domain'),
 'PNP':('7ZSO','human trimeric enzyme; avoids engineered four-Trp 4EAR mutant'),
 'PPARG':('1FM6','well-characterized ligand-binding domain complex'),
 'PTGER4':('5YWY','small-molecule-bound human receptor'),
 'S1PR1':('7TD4','ligand-bound receptor complex without T4-lysozyme insertion'),
 'SAMD8':('8IJR','human hexamer with resolved reaction-chamber ligands'),
 'SCD':('4ZYO','human desaturase with catalytic cofactors'),
 'SLC18B1':('9D7W','human transporter with small-molecule ligand'),
 'SSTR2':('7XN9','small-molecule-bound human receptor'),
 'VKORC1':('6WVH','high-resolution inhibitor-bound human enzyme'),
}

def get(url,p):
    if p.exists() and p.stat().st_size>100: return
    for n in range(3):
        try:
            r=requests.get(url,timeout=55);r.raise_for_status()
            tmp=p.with_suffix(p.suffix+'.part');tmp.write_bytes(r.content);tmp.replace(p);return
        except Exception:
            if n==2:raise
            time.sleep(n+1)

def fetch(code):
    cif=B/'structures'/f'{code}.cif'
    local=next((p for p in [OLD/'structures'/f'{code}.cif',OLD/'sources'/f'{code}.cif'] if p.exists()),None)
    if not cif.exists() and local:shutil.copy2(local,cif)
    get(f'https://files.rcsb.org/download/{code}.cif',cif)
    get(f'https://data.rcsb.org/rest/v1/core/entry/{code}',B/'sources'/f'{code}_entry.json')
    get(f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{code.lower()}',B/'sources'/f'{code}_SIFTS.json')
    get(f'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{code.lower()}.xml.gz',B/'sources'/f'{code}_SIFTS.xml.gz')
    st=gemmi.read_structure(str(cif))
    print(code,'models',len(st),'chains',len(st[0]),flush=True)
    return code

def table(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def review():
    inv={r['gene']:r for r in table(B/'inputs/47_verified_input_inventory.tsv')}
    rows=[]; ligrows=[]; refs=[]
    for gene,(code,reason) in SELECT.items():
        acc=inv[gene]['accession']
        seq=json.loads((OLD/'sources'/f'uniprot_{acc}.json').read_text(encoding='utf8'))['sequence']['value']
        cif=B/'structures'/f'{code}.cif';block=gemmi.cif.read_file(str(cif)).sole_block()
        st=gemmi.read_structure(str(cif))
        # struct_ref id -> database accession; sequence mappings use author numbering.
        refmap={x[0]:x[1] for x in block.find(['_struct_ref.id','_struct_ref.pdbx_db_accession'])}
        segments=[]
        for x in block.find(['_struct_ref_seq.ref_id','_struct_ref_seq.pdbx_strand_id',
                '_struct_ref_seq.pdbx_auth_seq_align_beg','_struct_ref_seq.pdbx_auth_seq_align_end',
                '_struct_ref_seq.db_align_beg','_struct_ref_seq.db_align_end']):
            if refmap.get(x[0],'').split('-')[0]!=acc: continue
            for ch in x[1].split(','):
                segments.append((ch.strip(),int(x[2]),int(x[3]),int(x[4]),int(x[5])))
        sifts=json.loads((B/'sources'/f'{code}_SIFTS.json').read_text())[code.lower()]['UniProt']
        target=sifts.get(acc)
        if target is None:
            target=next((v for k,v in sifts.items() if k.split('-')[0]==acc),None)
        assert target,(gene,code,'accession missing in SIFTS')
        mapped={}; mismatches=[]
        xml=ET.fromstring(gzip.decompress((B/'sources'/f'{code}_SIFTS.xml.gz').read_bytes()))
        residue_map={}
        for item in xml.iter():
            if item.tag.split('}')[-1]!='residue': continue
            crosses=[c.attrib for c in item if c.tag.split('}')[-1]=='crossRefDb']
            ur=next((x for x in crosses if x.get('dbSource')=='UniProt' and x.get('dbAccessionId','').split('-')[0]==acc),None)
            pr=next((x for x in crosses if x.get('dbSource')=='PDB'),None)
            if ur is not None and pr is not None:
                residue_map[(pr['dbChainId'],pr['dbResNum'])]=int(ur['dbResNum'])
        for chain in st[0]:
            ch=chain.name
            for res in chain:
                key=(ch,str(res.seqid))
                if res.het_flag!='A' or key not in residue_map: continue
                number=residue_map[key]
                # SIFTS extends the Q6NT16 mapping into the attached nanobody;
                # mmCIF database alignment and exact native sequence define this boundary.
                if gene=='SLC18B1':
                    if ch!='A' or not 232<=res.seqid.num<=638:continue
                    number=res.seqid.num-207
                if not 1<=number<=len(seq): continue
                aa=gemmi.find_tabulated_residue(res.name).one_letter_code
                mapped.setdefault(ch,[]).append((res,number,aa))
                if aa!=seq[number-1]:mismatches.append(f'{ch}:{res.seqid.num}:{aa}>{seq[number-1]}')
        assert mapped,(gene,code,'no canonical mapping')
        best=max(mapped,key=lambda c:len(mapped[c]))
        residues=mapped[best];nums=sorted({n for _,n,_ in residues})
        xyz=np.array([[a.pos.x,a.pos.y,a.pos.z] for r,_,_ in residues for a in r if not a.element.is_hydrogen])
        tree=cKDTree(xyz)
        ligs=[]
        for chain in st[0]:
            for res in chain:
                if res.het_flag=='A' or res.name in ['HOH','WAT','DOD']:continue
                atoms=[a for a in res if not a.element.is_hydrogen]
                if not atoms:continue
                pos=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in atoms])
                distances=tree.query(pos)[0]
                if min(distances)>6:continue
                ligand={'gene':gene,'pdb':code,'target_chain':best,'ligand_chain':chain.name,
                  'ligand_resnum':res.seqid.num,'CCD':res.name,'heavy_atoms':len(atoms),
                  'min_target_distance_A':round(float(min(distances)),3),
                  'heavy_atoms_within4A':int(sum(distances<=4)),
                  'center_x':float(pos.mean(0)[0]),'center_y':float(pos.mean(0)[1]),'center_z':float(pos.mean(0)[2])}
                ligs.append(ligand);ligrows.append(ligand)
        info=json.loads((B/'sources'/f'{code}_entry.json').read_text())
        rows.append({'gene':gene,'accession':acc,'pdb':code,'selection_reason':reason,
          'title':info['struct']['title'],'method':';'.join(x['method'] for x in info.get('exptl',[])),
          'resolution_A':';'.join(map(str,info.get('rcsb_entry_info',{}).get('resolution_combined',[]))),
          'canonical_chains':';'.join(mapped),'best_chain':best,'resolved_best_chain_residues':len(nums),
          'canonical_min':min(nums),'canonical_max':max(nums),'canonical_length':len(seq),
          'coverage_fraction':round(len(nums)/len(seq),4),'mutations':';'.join(mismatches),
          'mapping_segments':json.dumps(target['mappings']),'mapping_source':'PDBe SIFTS XML per-residue canonical mapping',
          'nearby_nonwater_CCDs':';'.join(sorted({x['CCD'] for x in ligs})),
          'total_structure_chains':len(st[0]),'file_SHA256':hashlib.file_digest(cif.open('rb'),'sha256').hexdigest(),
          'preparation_status':'needs_site_and_partner_review'})
        refs += [{'gene':gene,'chain':ch,'author_resnum':res.seqid.num,'author_resid':str(res.seqid),'canonical_resnum':n,
                 'resname':res.name,'canonical_aa':seq[n-1],'matches':aa==seq[n-1]}
                 for ch,rs in mapped.items() for res,n,aa in rs]
    write(B/'inputs/25_experimental_structure_review.tsv',rows)
    write(B/'inputs/Experimental_nearby_heteroatoms.tsv',ligrows)
    write(B/'inputs/Experimental_residue_mapping.tsv',refs)
    print('reviewed',len(rows),'experimental targets',flush=True)
    print(json.dumps([{'gene':r['gene'],'pdb':r['pdb'],'resolved':r['resolved_best_chain_residues'],
                      'mutations':r['mutations'],'ligands':r['nearby_nonwater_CCDs']} for r in rows],ensure_ascii=False,indent=2))

if __name__=='__main__':
    errors=[]
    with ThreadPoolExecutor(max_workers=4) as ex:
        futures={ex.submit(fetch,c):c for c in sorted({v[0] for v in SELECT.values()})}
        for fut in as_completed(futures):
            try:fut.result()
            except Exception as e:errors.append({'pdb':futures[fut],'error':str(e)});print('ERROR',errors[-1],flush=True)
    (B/'audit/Experimental_download_errors.json').write_text(json.dumps(errors,indent=2),encoding='utf8')
    if errors:raise SystemExit(1)
    review()
