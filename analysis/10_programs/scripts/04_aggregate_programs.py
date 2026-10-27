from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform

root=Path(sys.argv[1]);out=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26'
fit=json.loads((out/'metadata/Discovery_patient_fit_list.json').read_text())['patients']
geneinfo=pd.read_csv(out/'metadata/GSE161529_canonical_genes.tsv',sep='\t').set_index('canonical_id')
nodes=[];diagnostics=[];choices=[]
for pid in fit:
    p=out/'discovery/cnmf'/pid;assert (p/'FACTORIZATION_COMPLETE.json').exists(),f'{pid} is not complete'
    d=pd.read_csv(p/'K_diagnostics.tsv',sep='\t');diagnostics.append(d)
    stable=d[(d.silhouette>=.5)&(d.fraction_spectra_kept>=.5)&d.consensus_success.astype(bool)].sort_values('k')
    selected=None
    for _,r in stable.iterrows():
        if selected is None:selected=r
        elif (selected.prediction_error-r.prediction_error)/selected.prediction_error>=.02:selected=r
    choices.append(dict(patient_id=pid,selected_K=None if selected is None else int(selected.k),n_stable_K=len(stable),reason='Largest stable K with >=2% additional reconstruction improvement; no MFSD10 information'))
    for _,r in stable.iterrows():
        k=int(r.k);sc=pd.read_csv(p/f'{pid}.gene_spectra_score.k_{k}.dt_0_5.txt',sep='\t',index_col=0)
        assert 'ENSG00000109736' not in sc.columns
        for program,row in sc.iterrows():
            z=row[np.isfinite(row)&(row>0)].sort_values(ascending=False).head(50)
            if len(z)<50:continue
            nodes.append(dict(node=f'{pid}_K{k}_P{program}',patient_id=pid,k=k,program=str(program),genes=z.index.tolist(),weights=z.to_dict(),silhouette=float(r.silhouette),selected_K=selected is not None and k==int(selected.k)))
pd.concat(diagnostics).to_csv(out/'discovery/All_patient_K_diagnostics.tsv',sep='\t',index=False)
pd.DataFrame(choices).to_csv(out/'discovery/Selected_K_by_patient.tsv',sep='\t',index=False)
# Within-patient robustness across genuinely different K values, then deduplicate same patient's programs.
robust=[]
for n in nodes:
    n['within_rank_overlap']=max([len(set(n['genes'])&set(m['genes'])) for m in nodes if m['patient_id']==n['patient_id'] and m['k']!=n['k']] or [0])
    if n['within_rank_overlap']>=35:robust.append(n)
dedup=[]
for pid in fit:
    candidates=sorted([n for n in robust if n['patient_id']==pid],key=lambda n:(not n['selected_K'],-n['silhouette'],n['node']))
    kept=[]
    for n in candidates:
        if all(len(set(n['genes'])&set(m['genes']))<35 for m in kept):kept.append(n)
    dedup+=kept
assert len(dedup)>=3,'Too few robust programs for shared-program estimation'
G=[set(n['genes']) for n in dedup];D=np.zeros((len(G),len(G)))
for i in range(len(G)):
    for j in range(i):D[i,j]=D[j,i]=1-len(G[i]&G[j])/len(G[i]|G[j])
lab=fcluster(linkage(squareform(D),method='complete'),t=.89,criterion='distance')
meta=[]
for cl in sorted(set(lab)):
    mem=[dedup[i] for i in np.flatnonzero(lab==cl)];patients=sorted(set(n['patient_id'] for n in mem))
    if len(patients)<3:continue
    # One vote per patient per gene, even if multiple nonredundant components from a patient occur in the cluster.
    freq={};weights={}
    for p in patients:
        pn=[n for n in mem if n['patient_id']==p];pg=set().union(*(set(n['genes']) for n in pn))
        for g in pg:
            freq[g]=freq.get(g,0)+1
            w=np.mean([n['weights'].get(g,0)/max(n['weights'].values()) for n in pn]);weights[g]=weights.get(g,0)+w
    top=sorted(freq,key=lambda g:(-freq[g],-weights[g],g))[:50]
    if len(top)<50:continue
    symbols=geneinfo.loc[top,'gene_symbol'].fillna('').astype(str)
    meta.append(dict(genes=top,symbols=symbols.tolist(),n_support_patients=len(patients),support_patients=patients,support_nodes=[n['node'] for n in mem],gene_patient_frequency={g:freq[g] for g in top},MT_fraction=float(symbols.str.startswith('MT-').mean()),ribosomal_fraction=float(symbols.str.match(r'^RP[SL]\d').mean())))
meta=sorted(meta,key=lambda m:(-m['n_support_patients'],','.join(m['genes'])))
for i,m in enumerate(meta):m['program_id']=f'MP{i+1:02d}';m['technical_flag']=bool(m['MT_fraction']>.2 or m['ribosomal_fraction']>.4)
assert meta,'No shared meta-program passed prespecified criteria; do not force a program'
bank={'version':'discovery_programs_v1','cohort_discovery':'GSE161529','cohort_validation':'GSE176078','target_excluded':'ENSG00000109736/MFSD10/SLC75A1','score_method':'Mean log1p(CP10K) of fixed 50 genes; denominator all RNA UMI minus MFSD10; require >=80% mapped genes','selection_without_target':True,'meta_programs':meta}
fn=out/'frozen/Programs_v1.json';assert not fn.exists(),'Frozen programs already exist; do not overwrite'
fn.write_text(json.dumps(bank,indent=2),encoding='utf-8')
rows=[]
for m in meta:
    for g,s in zip(m['genes'],m['symbols']):rows.append(dict(program_id=m['program_id'],canonical_id=g,gene_symbol=s,n_support_patients=m['n_support_patients'],gene_patient_support=m['gene_patient_frequency'][g],support_patients=';'.join(m['support_patients']),technical_flag=m['technical_flag']))
pd.DataFrame(rows).to_csv(out/'frozen/Program_genes_v1.tsv',sep='\t',index=False)
pd.DataFrame([{k:v for k,v in n.items() if k not in ['genes','weights']} for n in nodes]).to_csv(out/'discovery/All_components_robustness.tsv',sep='\t',index=False)
(out/'discovery/Robust_program_membership.json').write_text(json.dumps(dedup,indent=2),encoding='utf-8')
summary=pd.DataFrame([dict(program_id=m['program_id'],n_support_patients=m['n_support_patients'],support_patients=';'.join(m['support_patients']),MT_fraction=m['MT_fraction'],ribosomal_fraction=m['ribosomal_fraction'],technical_flag=m['technical_flag'],genes=';'.join(m['symbols'])) for m in meta])
summary.to_csv(out/'frozen/Program_summary_v1.tsv',sep='\t',index=False)
status=json.loads((out/'audit/PROGRESS.json').read_text());status.update(programs_frozen=True,stage='discovery_target_association',n_meta_programs=len(meta));(out/'audit/PROGRESS.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
print(summary[['program_id','n_support_patients','technical_flag']].to_string(index=False));print('PROGRAMS_FROZEN_SHA256',hashlib.sha256(fn.read_bytes()).hexdigest())
