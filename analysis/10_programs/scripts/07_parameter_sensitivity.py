import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,json,warnings,contextlib,concurrent.futures
import numpy as np,pandas as pd,anndata as ad
from cnmf import cNMF
from sklearn.decomposition import non_negative_factorization
from scipy.optimize import linear_sum_assignment

class AuditedCNMF(cNMF):
    def _nmf(self,X,nmf_kwargs):
        params=nmf_kwargs
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always');U,H,n=non_negative_factorization(X,**params)
        rec={'seed':str(params.get('random_state')),'K':int(params['n_components']),'n_iter':int(n),'max_iter':int(params.get('max_iter',200)),'limit_reached':bool(n>=params.get('max_iter',200)),'warnings':[str(w.message) for w in ws]}
        p=Path(self.output_dir)/self.name/'iteration_audit.jsonl'
        with p.open('a',encoding='utf-8') as f:f.write(json.dumps(rec)+'\n')
        return H,U

def task(args):
    outstr,pid,scenario,k=args;out=Path(outstr);name=pid+'__'+scenario;base=out/'discovery/sensitivity';marker=base/name/'SENSITIVITY_COMPLETE.json'
    if marker.exists():return json.loads(marker.read_text())
    c=AuditedCNMF(output_dir=str(base),name=name)
    with (out/'logs'/f'sensitivity_{name}.log').open('a',encoding='utf-8',buffering=1) as f,contextlib.redirect_stdout(f),contextlib.redirect_stderr(f):
        indir='discovery_all_confirmed_inputs' if scenario=='with_flagged_doublets' else 'discovery_patient_inputs'
        obj=ad.read_h5ad(out/'cache'/indir/f'{pid}.h5ad');rng=np.random.default_rng(20261026+sum(map(ord,scenario)))
        if scenario.startswith('sub80'):obj=obj[np.sort(rng.choice(obj.n_obs,int(obj.n_obs*.8),replace=False))].copy()
        if scenario=='UMI50':obj.X.data=rng.binomial(obj.X.data,.5).astype(np.int32);obj.X.eliminate_zeros()
        keep=np.asarray(obj.X.sum(axis=0)).ravel()>0;obj=obj[:,keep].copy();assert 'ENSG00000109736' not in obj.var_names
        hvg=1000 if scenario=='HVG1000' else 3000 if scenario=='HVG3000' else 2000
        fn=base/name/'perturbed_input.h5ad';obj.write_h5ad(fn,compression='gzip')
        if not Path(c.paths['nmf_replicate_parameters']).exists():c.prepare(str(fn),components=[int(k)],n_iter=20,seed=20261026,num_highvar_genes=hvg,max_NMF_iter=2000)
        c.update_nmf_iter_params();c.factorize(skip_completed_runs=True);c.combine();c.consensus(int(k),density_threshold=.5,show_clustering=False)
        orig=out/'discovery/cnmf'/pid
        a=pd.read_csv(orig/f'{pid}.gene_spectra_score.k_{k}.dt_0_5.txt',sep='\t',index_col=0)
        b=pd.read_csv(base/name/f'{name}.gene_spectra_score.k_{k}.dt_0_5.txt',sep='\t',index_col=0)
        common=a.columns.intersection(b.columns);A=np.maximum(a[common].to_numpy(),0);B=np.maximum(b[common].to_numpy(),0)
        A/=np.maximum(np.linalg.norm(A,axis=1,keepdims=True),1e-10);B/=np.maximum(np.linalg.norm(B,axis=1,keepdims=True),1e-10)
        cos=A@B.T;ii,jj=linear_sum_assignment(-cos)
        ua=pd.read_csv(orig/f'{pid}.usages.k_{k}.dt_0_5.consensus.txt',sep='\t',index_col=0);ub=pd.read_csv(base/name/f'{name}.usages.k_{k}.dt_0_5.consensus.txt',sep='\t',index_col=0)
        cells=ua.index.intersection(ub.index);ua=ua.loc[cells].div(ua.loc[cells].sum(axis=1),axis=0);ub=ub.loc[cells].div(ub.loc[cells].sum(axis=1),axis=0)
        rows=[]
        for i,j in zip(ii,jj):
            ga=set(a.iloc[i].nlargest(50).index);gb=set(b.iloc[j].nlargest(50).index)
            rows.append(dict(patient_id=pid,scenario=scenario,K=int(k),base_program=str(a.index[i]),matched_program=str(b.index[j]),positive_loading_cosine=float(cos[i,j]),top50_overlap=len(ga&gb),top50_Jaccard=len(ga&gb)/len(ga|gb),usage_correlation=float(np.corrcoef(ua.iloc[:,i],ub.iloc[:,j])[0,1]),n_cells=obj.n_obs,n_HVG=hvg,replicates=20))
        pd.DataFrame(rows).to_csv(base/name/'Component_comparison.tsv',sep='\t',index=False)
        r={'patient_id':pid,'scenario':scenario,'median_cosine':float(np.median([r['positive_loading_cosine'] for r in rows])),'median_top50_overlap':float(np.median([r['top50_overlap'] for r in rows])),'median_usage_correlation':float(np.median([r['usage_correlation'] for r in rows])),'status':'complete'}
        marker.write_text(json.dumps(r,indent=2),encoding='utf-8');return r

if __name__=='__main__':
    out=Path(sys.argv[1])/'outputs/MFSD10_程序发现与独立验证_2026-10-26';chosen=pd.read_csv(out/'discovery/Selected_K_by_patient.tsv',sep='\t').set_index('patient_id')
    if (out/'cache/GSE161529_three_patients_all_confirmed.mtx').exists():
        from scipy.io import mmread
        idir=out/'cache/discovery_all_confirmed_inputs';idir.mkdir(exist_ok=True)
        if len(list(idir.glob('*.h5ad')))<3:
            md=pd.read_csv(out/'cache/GSE161529_three_patients_all_confirmed_metadata.tsv',sep='\t').set_index('cell_id');genes=pd.read_csv(out/'metadata/GSE161529_canonical_genes.tsv',sep='\t').set_index('canonical_id')
            with (out/'cache/GSE161529_three_patients_all_confirmed.mtx').open('rb') as f:X=mmread(f).T.tocsr()
            obj=ad.AnnData(X=X,obs=md,var=genes)
            for pid in ['Cancer__0126','Cancer__0043','Cancer__0031']:
                a=obj[obj.obs.patient_id==pid,obj.var_names!='ENSG00000109736'].copy();nz=np.asarray((a.X>0).sum(axis=0)).ravel();a=a[:,nz>=max(5,int(np.ceil(a.n_obs*.03)))].copy();a.write_h5ad(idir/f'{pid}.h5ad',compression='gzip')
    jobs=[]
    for pid in ['Cancer__0126','Cancer__0043','Cancer__0031']:
        k=chosen.loc[pid,'selected_K'];assert np.isfinite(k)
        for scenario in ['HVG1000','HVG3000','sub80_A','sub80_B','UMI50']:jobs.append((str(out),pid,scenario,int(k)))
        if (out/'cache/discovery_all_confirmed_inputs'/f'{pid}.h5ad').exists():jobs.append((str(out),pid,'with_flagged_doublets',int(k)))
    # Targeted local convergence check, motivated by an observed primary-run iteration-limit warning.
    pid='Cancer__0173';k=chosen.loc[pid,'selected_K']
    if np.isfinite(k):jobs.append((str(out),pid,'maxiter2000',int(k)))
    summary=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        for r in pool.map(task,jobs):summary.append(r);print(json.dumps(r),flush=True)
    pd.DataFrame(summary).to_csv(out/'discovery/Parameter_sensitivity_summary.tsv',sep='\t',index=False)
    parts=list((out/'discovery/sensitivity').glob('*/Component_comparison.tsv'));pd.concat([pd.read_csv(p,sep='\t') for p in parts]).to_csv(out/'discovery/Parameter_sensitivity_components.tsv',sep='\t',index=False)
    print('PARAMETER_SENSITIVITY_COMPLETE')
