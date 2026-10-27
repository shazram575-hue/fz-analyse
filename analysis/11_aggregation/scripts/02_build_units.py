import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,time,hashlib,gc
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from reused_original_models import prepare

TARGET='ENSG00000109736'
def log(*s):print(time.strftime('%H:%M:%S'),*s,flush=True)
def center_groups(Z,lab,k):
    ans=np.zeros((k,Z.shape[1]));np.add.at(ans,lab,Z);return ans/np.bincount(lab,minlength=k)[:,None]
def partition(X,obs,var,cohort,seed):
    n=X.shape[0];k=max(1,n//20);caps=np.full(k,n//k,dtype=int);caps[:n%k]+=1
    ids=var.index.astype(str).to_numpy() if cohort=='discovery' else var.unified_gene_id.astype(str).to_numpy()
    symbols=var.gene_symbol.fillna('').astype(str).to_numpy() if cohort=='discovery' else var.author_gene.astype(str).to_numpy()
    target=(ids==TARGET)|(symbols=='MFSD10')|(symbols=='SLC75A1');assert target.sum()==1
    other=X[:,~target].astype(np.float64).tocsr();den=np.asarray(other.sum(axis=1)).ravel();assert (den>0).all()
    sym=symbols[~target];mt=np.char.startswith(sym.astype(str),'MT-');mtp=np.asarray(other[:,mt].sum(axis=1)).ravel()/den*100
    rng=np.random.default_rng(seed)
    if n<40:return np.zeros(n,int),np.zeros(n,int),np.zeros(n),np.zeros(n),[],0
    other=other.multiply((1e4/den)[:,None]).tocsr();other.data=np.log1p(other.data)
    mu=np.asarray(other.mean(axis=0)).ravel();va=np.asarray(other.power(2).mean(axis=0)).ravel()-mu*mu
    nz=np.asarray((other>0).sum(axis=0)).ravel();eligible=(nz>=max(3,int(np.ceil(n*.01))))&(va>1e-8)
    order=np.flatnonzero(eligible);order=order[np.argsort(-va[order],kind='stable')[:2000]]
    assert len(order)>=10,'Insufficient non-target state features'
    A=other[:,order].toarray();A/=np.sqrt(va[order])[None,:]
    pcs=min(30,n-1,len(order)-1);Z=PCA(n_components=pcs,svd_solver='randomized',random_state=seed).fit_transform(A)
    C=np.column_stack([np.ones(n),np.log1p(den),mtp]);Z-=C@np.linalg.lstsq(C,Z,rcond=None)[0]
    km=KMeans(n_clusters=k,n_init=10,max_iter=300,random_state=seed).fit(Z);centers=km.cluster_centers_
    slots=np.repeat(np.arange(k),caps)
    for _ in range(2):
        cost=np.repeat(cdist(Z,centers,metric='sqeuclidean'),caps,axis=1)
        rr,cc=linear_sum_assignment(cost);lab=np.empty(n,int);lab[rr]=slots[cc];centers=center_groups(Z,lab,k)
    randomlab=rng.permutation(slots);rc=center_groups(Z,randomlab,k)
    dist=((Z-centers[lab])**2).sum(axis=1);rdist=((Z-rc[randomlab])**2).sum(axis=1)
    hvg=ids[~target][order].tolist();assert TARGET not in hvg
    return lab,randomlab,dist,rdist,hvg,pcs

def aggregate(obj,obs,scores,labels,compact,cohort,kind):
    keys=pd.unique(labels);idx=pd.Index(keys).get_indexer(labels);G=sp.csr_matrix((np.ones(len(labels),np.int32),(idx,np.arange(len(labels)))),shape=(len(keys),len(labels)))
    A=(G@obj.X).tocsr();A.eliminate_zeros();assert np.array_equal(np.asarray(A.sum(axis=0)),np.asarray(obj.X.sum(axis=0)))
    rows=[]
    for i,key in enumerate(keys):
        j=np.flatnonzero(idx==i);d=obs.iloc[j];total=d.total_without_target.to_numpy()+d.MFSD10_counts.to_numpy()
        mt_raw=np.sum(d.mt_pct.to_numpy()*d.total_without_target.to_numpy())/total.sum()
        a=dict(unit_id=key,patient_id=d.patient_id.iloc[0],sample_id=d.unit_sample.iloc[0],clinical_group=str(d.clinical_group.iloc[0]),treatment_status=str(d.unit_treatment.iloc[0]),cohort=cohort,partition=kind,n_cells=len(j),nCount_RNA=int(total.sum()),nFeature_RNA=int(A[i].nnz),**{'percent.mt':mt_raw,'mt_pct_recomputed':mt_raw},median_cell_UMI=float(np.median(total)),mean_state_squared_distance=float(np.mean(compact[j])))
        for mp in scores:a['cellmean_'+mp]=float(scores.iloc[j][mp].mean())
        rows.append(a)
    md=pd.DataFrame(rows).set_index('unit_id');result=ad.AnnData(X=A,obs=md,var=obj.var.copy())
    u,sc,cov,_,_=prepare(result,bank,cohort);u['n_features_aggregate']=np.diff(A.indptr)
    for mp in sc:u[mp]=sc[mp]
    # Mapped gene coverage and observed detection coverage have different meanings.
    ids=obj.var.index.astype(str).to_numpy() if cohort=='discovery' else obj.var.unified_gene_id.astype(str).to_numpy()
    for m in bank['meta_programs']:
        jj=np.flatnonzero(np.isin(ids,m['genes']));u[m['program_id']+'_mapped_fraction']=len(jj)/50
        u[m['program_id']+'_detected_fraction']=np.asarray((A[:,jj]>0).sum(axis=1)).ravel()/50
    result.obs=u;result.uns['analysis_type']='Exploratory reanalysis after both cohorts viewed';result.uns['fixed_unit_target_size']=20
    result.write_h5ad(out/'objects'/f'{cohort}_{kind}_raw_count_units.h5ad',compression='gzip')
    u.to_csv(out/'results'/f'{cohort}_{kind}_unit_summary.tsv.gz',sep='\t');return u,cov

def main(root):
    global out,bank
    old=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26';out=root/'outputs/MFSD10_稀疏与聚合探索_2026-10-26';bank=json.loads((old/'frozen/Programs_v1.json').read_text())
    cache=out/'objects/partition_checkpoints';cache.mkdir(exist_ok=True)
    fit=json.loads((old/'metadata/Discovery_patient_fit_list.json').read_text())['patients'];external=json.loads((old/'audit/External_validation_summary.json').read_text())['primary_naive_patient_ids']
    audits=[];roles=[]
    for cohort,fn,primary in [('discovery','GSE161529_confirmed_malignant_singlets_counts.h5ad',fit),('validation','GSE176078_malignant_candidates_author_counts_clinical_v2.h5ad',external)]:
        done=out/'audit'/f'{cohort}_AGGREGATION_COMPLETE.json'
        if done.exists():log('REUSE_AGGREGATION',cohort);continue
        obj=ad.read_h5ad(old/'objects'/fn);obs,scores,coverage,lookup,denom=prepare(obj,bank,cohort)
        obs['unit_sample']=obs['gsm'].astype(str) if cohort=='discovery' else obs.sample_id.astype(str)
        obs['unit_treatment']=obs.treatment_status.astype(str)
        labels=np.full(len(obs),'',dtype=object);randomlabels=labels.copy();compact=np.zeros(len(obs));randomcompact=np.zeros(len(obs));features=[]
        for (patient,sample),positions in obs.groupby(['patient_id','unit_sample'],observed=True).indices.items():
            positions=np.asarray(positions);key=cohort+'__'+patient+'__'+sample;path=cache/(key+'.npz')
            seed=20261026+int(hashlib.sha256(key.encode()).hexdigest()[:7],16)
            if path.exists():
                z=np.load(path,allow_pickle=False);assert np.array_equal(z['cells'],obs.index[positions].to_numpy(dtype=str));lab,rlab,dist,rdist,hvg,pcs=z['lab'],z['randomlab'],z['dist'],z['randomdist'],z['hvg'].tolist(),int(z['pcs'])
            else:
                log('PARTITION',cohort,patient,sample,len(positions))
                lab,rlab,dist,rdist,hvg,pcs=partition(obj.X[positions],obs.iloc[positions],obj.var,cohort,seed)
                np.savez_compressed(path,cells=obs.index[positions].to_numpy(dtype=str),lab=lab,randomlab=rlab,dist=dist,randomdist=rdist,hvg=np.asarray(hvg,dtype=str),pcs=np.array(pcs))
            labels[positions]=[key+f'__state_{i:04d}' for i in lab];randomlabels[positions]=[key+f'__random_{i:04d}' for i in rlab];compact[positions]=dist;randomcompact[positions]=rdist
            features.extend([dict(cohort=cohort,patient_id=patient,sample_id=sample,feature_id=g) for g in hvg])
        assert (labels!='').all() and len(np.unique(obs.index))==len(labels)
        obs[['patient_id','unit_sample']].assign(state_unit=labels,random_unit=randomlabels).to_csv(out/'metadata'/f'{cohort}_cell_unit_membership.tsv.gz',sep='\t')
        pd.DataFrame(features,columns=['cohort','patient_id','sample_id','feature_id']).to_csv(out/'metadata'/f'{cohort}_neighborhood_features.tsv.gz',sep='\t',index=False)
        st,cov=aggregate(obj,obs,scores,labels,compact,cohort,'state');ra,_=aggregate(obj,obs,scores,randomlabels,randomcompact,cohort,'random')
        # Full-gene patient pseudobulk is an audit of the already performed patient-level analysis.
        pats=pd.unique(obs.patient_id);P=sp.csr_matrix((np.ones(len(obs),np.int32),(pd.Index(pats).get_indexer(obs.patient_id),np.arange(len(obs)))),shape=(len(pats),len(obs)))
        PB=(P@obj.X).tocsr();pbmd=[]
        for p in pats:
            d=obs[obs.patient_id==p];us=st[st.patient_id==p]
            a=dict(patient_id=p,clinical_group=str(d.clinical_group.iloc[0]),treatment_status=str(d.unit_treatment.iloc[0]),n_cells=len(d),n_units=len(us),primary_patient=p in primary,MFSD10_total_UMI=float(d.MFSD10_counts.sum()),RNA_UMI_excluding_target=float(d.total_without_target.sum()),MFSD10_log2_CPM=float(np.log2(1+1e6*d.MFSD10_counts.sum()/d.total_without_target.sum())),median_depth=float(d.total_without_target.median()),median_mt=float(d.mt_pct.median()),singlecell_zero_fraction=float((d.MFSD10_counts==0).mean()),aggregate_zero_fraction=float((us.MFSD10_counts==0).mean()))
            for mp in scores:a[mp]=float(scores.loc[d.index,mp].mean())
            pbmd.append(a)
        pbmd=pd.DataFrame(pbmd).set_index('patient_id');ad.AnnData(X=PB,obs=pbmd,var=obj.var.copy()).write_h5ad(out/'objects'/f'{cohort}_patient_pseudobulk_raw_counts.h5ad',compression='gzip');pbmd.to_csv(out/'results'/f'{cohort}_patient_pseudobulk_audit.tsv',sep='\t')
        oldfn=old/('discovery/MFSD10_between_patient_inputs.tsv' if cohort=='discovery' else 'validation/MFSD10_between_patient_inputs_Naive.tsv');oldpb=pd.read_csv(oldfn,sep='\t').set_index('patient_id')
        for column in ['MFSD10_log2_CPM',*scores.columns]:assert np.allclose(pbmd.loc[oldpb.index,column],oldpb[column],atol=1e-10,rtol=1e-10),column
        for p,row in pbmd.iterrows():roles.append(dict(cohort=cohort,patient_id=p,n_cells=int(row.n_cells),n_units=int(row.n_units),primary_patient=bool(row.primary_patient),treatment_status=row.treatment_status,reason='Primary patient unchanged' if row.primary_patient else 'Previously outside main discovery / treated secondary; not selected by MFSD10'))
        audit={'cohort':cohort,'n_cells':obj.n_obs,'state_units':len(st),'random_units':len(ra),'nonoverlap':True,'UMI_conserved':True,'patient_pseudobulk_matches_existing':True,'MFSD10_excluded_from_neighborhood':True,'state_mean_squared_distance':float(compact.mean()),'random_mean_squared_distance':float(randomcompact.mean())}
        done.write_text(json.dumps(audit,indent=2),encoding='utf-8');audits.append(audit);log('AGGREGATION_COMPLETE',audit)
        del obj,obs,st,ra,PB;gc.collect()
    roleparts=[]
    for cohort in ['discovery','validation']:
        a=pd.read_csv(out/'results'/f'{cohort}_patient_pseudobulk_audit.tsv',sep='\t');a['cohort']=cohort;roleparts.append(a)
    pd.concat(roleparts).to_csv(out/'results/All_patient_unit_availability.tsv',sep='\t',index=False)
    p=out/'audit/PROGRESS.json';s=json.loads(p.read_text());s.update(aggregation_complete=True,stage='patient_effect_estimation');p.write_text(json.dumps(s,indent=2),encoding='utf-8')

if __name__=='__main__':main(Path(sys.argv[1]))
