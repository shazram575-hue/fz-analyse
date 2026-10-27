"""Fixed cohort-separate scores, patient-level associations, and gene-gene coherence."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import json,itertools
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.stats import rankdata,t as student_t,norm,wilcoxon
from statsmodels.stats.multitest import multipletests

TARGET='ENSG00000109736'
def bh(p):
    p=np.asarray(p,float);ans=np.full(p.shape,np.nan);ok=np.isfinite(p)
    if ok.any():ans[ok]=multipletests(p[ok],method='fdr_bh')[1]
    return ans
def resid(y,c):return y-c@np.linalg.lstsq(c,y,rcond=None)[0]
def ranks_covariates(obs):
    a=np.c_[np.log1p(obs['total_without_target'].to_numpy(float)),np.log1p(obs['features_without_target'].to_numpy(float)),obs['mt_pct'].to_numpy(float)]
    cols=[np.ones(len(a))]
    for j in range(a.shape[1]):
        r=rankdata(a[:,j]);r-=r.mean();s=r.std()
        if s>0:cols.append(r/s)
    return np.column_stack(cols)
def partial_rho(x,y,c):
    rx=resid(rankdata(x),c);ry=resid(rankdata(y),c)
    if np.std(rx)<1e-10 or np.std(ry)<1e-10:return np.nan
    return float(np.clip(np.corrcoef(rx,ry)[0,1],-.999999,.999999))
def sign_meta(rhos,seed=20261026,direction=0):
    rhos=np.asarray(rhos,float);rhos=rhos[np.isfinite(rhos)];n=len(rhos)
    if n==0:return dict(n_patients=0,effect_rho=np.nan,p=np.nan,ci_low=np.nan,ci_high=np.nan)
    z=np.arctanh(np.clip(rhos,-.999,.999));mu=z.mean();effect=np.tanh(mu)
    if n<6:return dict(n_patients=n,effect_rho=effect,p=np.nan,ci_low=np.nan,ci_high=np.nan)
    rng=np.random.default_rng(seed)
    if n<=16:signs=np.array(list(itertools.product([-1,1],repeat=n)),dtype=np.int8)
    else:signs=rng.choice([-1,1],size=(50000,n))
    null=signs@z/n
    hits=np.sum(np.abs(null)>=abs(mu)-1e-12) if direction==0 else np.sum(direction*null>=direction*mu-1e-12)
    p=(hits+1)/(len(null)+1)
    boot=np.tanh(z[rng.integers(0,n,size=(5000,n))].mean(axis=1));lo,hi=np.quantile(boot,[.025,.975])
    return dict(n_patients=n,effect_rho=float(effect),p=float(p),ci_low=float(lo),ci_high=float(hi))
def column_sum(X,indices):return np.asarray(X[:,indices].sum(axis=1)).ravel()
def extract_normalized(X,indices,denom):
    # Only requested gene subsets (never the complete cohort matrix) become dense.
    z=X[:,indices].toarray().astype(np.float64);return np.log1p(z/denom[:,None]*1e4)
def prepare(obj,bank,cohort):
    obs=obj.obs.copy();X=obj.X.tocsr()
    if cohort=='discovery':
        ids=obj.var_names.to_numpy();symbols=obj.var['gene_symbol'].fillna('').astype(str).to_numpy();obs['mt_pct']=obs['percent.mt'].to_numpy(float)
    else:
        ids=obj.var['unified_gene_id'].astype(str).to_numpy();symbols=obj.var['author_gene'].astype(str).to_numpy();obs['mt_pct']=obs['mt_pct_recomputed'].to_numpy(float)
    target=np.flatnonzero((ids==TARGET)|(symbols=='MFSD10')|(symbols=='SLC75A1'));assert len(target)==1
    mf=column_sum(X,target);total=np.asarray(X.sum(axis=1)).ravel();denom=total-mf;assert (denom>0).all()
    obs['MFSD10_counts']=mf;obs['total_without_target']=denom;obs['features_without_target']=np.diff(X.indptr)-(mf>0)
    obs['MFSD10_log1p_CP10K']=np.log1p(mf/denom*1e4);obs['mt_pct']=obs['mt_pct']*total/denom
    columns={}
    for i,g in enumerate(ids):columns.setdefault(g,[]).append(i)
    # Exclude ambiguous multi-row mappings from program scoring rather than silently counting them twice.
    lookup={g:ii[0] for g,ii in columns.items() if len(ii)==1 and g!=TARGET}
    scores=pd.DataFrame(index=obs.index);coverage=[]
    for m in bank['meta_programs']:
        genes=[g for g in m['genes'] if g in lookup];cov=len(genes)/len(m['genes'])
        coverage.append(dict(cohort=cohort,program_id=m['program_id'],n_genes_expected=50,n_genes_mapped=len(genes),fraction_mapped=cov,missing_genes=';'.join(g for g in m['genes'] if g not in lookup)))
        scores[m['program_id']]=extract_normalized(X,[lookup[g] for g in genes],denom).mean(axis=1) if cov>=.8 else np.nan
    return obs,scores,pd.DataFrame(coverage),lookup,denom

def associations(obs,scores,patients,seed=20261026,directions=None):
    obs=obs[obs.patient_id.isin(patients)].copy();scores=scores.loc[obs.index]
    rows=[];patient_rows=[];between_rows=[]
    for p,d in obs.groupby('patient_id',observed=True):
        idx=d.index
        if len(d)<50:continue
        c=ranks_covariates(d)
        for mp in scores:
            y=scores.loc[idx,mp].to_numpy(float);ok=np.isfinite(y)
            rho=partial_rho(d.MFSD10_log1p_CP10K.to_numpy()[ok],y[ok],c[ok]) if ok.sum()>=50 else np.nan
            patient_rows.append(dict(patient_id=p,program_id=mp,n_cells=len(d),rho=rho,clinical_group=str(d.clinical_group.iloc[0]),n_MFSD10_detected=int((d.MFSD10_counts>0).sum()),median_UMI=float(d.total_without_target.median())))
    per=pd.DataFrame(patient_rows)
    for mp in scores:
        rr=per[per.program_id==mp];direction=0 if directions is None else directions.get((mp,'within'),0)
        s=sign_meta(rr.rho,seed=seed,direction=direction);s.update(program_id=mp,association='within',n_cells=int(rr.n_cells.sum()),note='Patient equal-weight Fisher-z aggregation; patient sign-flip P and patient bootstrap CI')
        rows.append(s)
    # Between-patient summaries; no cell-level independent replicate assumptions.
    means=[]
    for p,d in obs.groupby('patient_id',observed=True):
        if len(d)<50:continue
        a=dict(patient_id=p,clinical_group=str(d.clinical_group.iloc[0]),MFSD10_log2_CPM=np.log2(1+1e6*d.MFSD10_counts.sum()/d.total_without_target.sum()),median_depth=float(d.total_without_target.median()),median_mt=float(d.mt_pct.median()))
        a.update({mp:float(scores.loc[d.index,mp].mean()) for mp in scores});means.append(a)
    means=pd.DataFrame(means)
    for mp in scores:
        b=means[np.isfinite(means[mp])].copy();n=len(b)
        cats=pd.get_dummies(b.clinical_group,drop_first=True,dtype=float).to_numpy()
        cov=np.column_stack([np.ones(n),cats,rankdata(np.log1p(b.median_depth)),rankdata(b.median_mt)])
        rank=np.linalg.matrix_rank(cov);rho=partial_rho(b.MFSD10_log2_CPM.to_numpy(),b[mp].to_numpy(),cov) if n>rank+2 else np.nan;df=n-rank-1
        estimable=n>=12 and df>=5 and np.isfinite(rho)
        stat=rho*np.sqrt(df/max(1-rho*rho,1e-12)) if estimable else np.nan
        direction=0 if directions is None else directions.get((mp,'between'),0)
        pv=float(2*student_t.sf(abs(stat),df)) if estimable else np.nan
        if estimable and direction:pv=float(student_t.sf(direction*stat,df))
        se=1/np.sqrt(max(n-rank-3,1))
        lo,hi=np.tanh(np.arctanh(np.clip(rho,-.999,.999))+np.array([-1,1])*1.96*se) if estimable else [np.nan,np.nan]
        rows.append(dict(program_id=mp,association='between',n_patients=n,effect_rho=rho,p=pv,ci_low=lo,ci_high=hi,residual_df=df,n_cells=int(obs.shape[0]),note='Patient-level partial rank rho; explicit clinical IHC, depth and mitochondrial covariates; descriptive if n<12 or df<5'))
    result=pd.DataFrame(rows);result['BH_q']=bh(result.p)
    return result,per,means

def mean_offdiag(residual):
    residual=residual-residual.mean(axis=0);length=np.sqrt((residual**2).sum(axis=0));keep=length>1e-8
    if keep.sum()<3:return np.nan
    Z=residual[:,keep]/length[keep];g=Z.shape[1]
    return float(((Z.sum(axis=1)**2).sum()-g)/(g*(g-1)))

def coherence(obj,obs,bank,lookup,denom,patients,seed=20261026):
    X=obj.X.tocsr();rng=np.random.default_rng(seed)
    detected=np.asarray((X>0).sum(axis=0)).ravel();mean_linear=np.asarray(X.T@(1/denom)).ravel()/len(obs)
    union=set(g for m in bank['meta_programs'] for g in m['genes'])
    universe=[g for g,i in lookup.items() if detected[i]>=max(5,int(len(obs)*.03))]
    vals=np.array([np.log1p(mean_linear[lookup[g]]*1e4) for g in universe]);bins=np.unique(np.quantile(vals,np.linspace(0,1,11)))
    gene_bin={g:int(np.searchsorted(bins,v,side='right')) for g,v in zip(universe,vals)}
    pools={}
    for g in universe:
        if g not in union:pools.setdefault(gene_bin[g],[]).append(g)
    nulls={}
    for m in bank['meta_programs']:
        genes=[g for g in m['genes'] if g in lookup]
        choices=[]
        for _ in range(25):
            selected=[]
            for g in genes:
                bn=gene_bin.get(g,5);pool=pools.get(bn,universe);selected.append(pool[rng.integers(len(pool))])
            choices.append(list(dict.fromkeys(selected)))
        nulls[m['program_id']]=choices
    rows=[]
    for p in patients:
        ix=np.flatnonzero(obs.patient_id.to_numpy()==p)
        if len(ix)<50:continue
        dd=obs.iloc[ix];C=ranks_covariates(dd)
        need=sorted(set(g for m in bank['meta_programs'] for g in m['genes'] if g in lookup)|set(g for sets in nulls.values() for ss in sets for g in ss))
        E=extract_normalized(X[ix],[lookup[g] for g in need],denom[ix]);R=resid(E,C);pos={g:i for i,g in enumerate(need)}
        for m in bank['meta_programs']:
            genes=[g for g in m['genes'] if g in pos]
            if len(genes)<40:continue
            val=mean_offdiag(R[:,[pos[g] for g in genes]])
            ns=[mean_offdiag(R[:,[pos[g] for g in ss]]) for ss in nulls[m['program_id']]]
            null=float(np.nanmean(ns));rows.append(dict(patient_id=p,program_id=m['program_id'],n_cells=len(ix),n_genes=len(genes),residual_gene_coherence=val,matched_random_coherence=null,coherence_excess=val-null))
    per=pd.DataFrame(rows);summary=[]
    for mp,d in per.groupby('program_id',observed=True):
        v=d.coherence_excess.to_numpy();v=v[np.isfinite(v)]
        pv=float(wilcoxon(v,alternative='greater',zero_method='wilcox').pvalue) if len(v)>=6 and np.any(v!=0) else np.nan
        summary.append(dict(program_id=mp,n_patients=len(v),median_coherence=float(d.residual_gene_coherence.median()),median_matched_null=float(d.matched_random_coherence.median()),median_excess=float(np.median(v)),patients_positive=int((v>0).sum()),p=pv))
    summary=pd.DataFrame(summary);summary['BH_q']=bh(summary.p)
    return per,summary
