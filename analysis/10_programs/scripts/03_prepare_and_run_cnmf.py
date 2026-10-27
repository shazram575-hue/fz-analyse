import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,time,contextlib,traceback,concurrent.futures
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.io import mmread
import anndata as ad
from cnmf import cNMF
from cnmf.cnmf import load_df_from_npz

def log(*x):print(time.strftime('%H:%M:%S'),*x,flush=True)
def one_patient(args):
    path,pid=args;out=Path(path);pdir=out/'discovery/cnmf';done=pdir/pid/'FACTORIZATION_COMPLETE.json'
    if done.exists():return json.loads(done.read_text())
    with (out/'logs'/f'cnmf_{pid}.log').open('a',encoding='utf-8',buffering=1) as f,contextlib.redirect_stdout(f),contextlib.redirect_stderr(f):
        c=cNMF(output_dir=str(pdir),name=pid)
        try:
            start=time.time()
            if not Path(c.paths['nmf_replicate_parameters']).exists():
                c.prepare(str(out/'cache/discovery_patient_inputs'/f'{pid}.h5ad'),components=[4,6,8],n_iter=100,seed=20261026+int(pid.split('__')[-1]),num_highvar_genes=2000,max_NMF_iter=1000)
            c.update_nmf_iter_params();c.factorize(worker_i=0,total_workers=1,skip_completed_runs=True)
            c.combine(skip_missing_files=False)
            rows=[]
            for k in [4,6,8]:
                stats=c.consensus(k,skip_density_and_return_after_stats=True,show_clustering=False)
                d={str(idx):float(v) for idx,v in stats.iloc[:,0].items()};d.update(patient_id=pid)
                try:
                    c.consensus(k,density_threshold=.5,show_clustering=True,close_clustergram_fig=True)
                    density=load_df_from_npz(c.paths['local_density_cache']%k).iloc[:,0]
                    d['fraction_spectra_kept']=float((density<.5).mean());d['consensus_success']=True
                except Exception as ex:
                    d['consensus_success']=False;d['consensus_error']=str(ex);d['fraction_spectra_kept']=0.
                rows.append(d)
            pd.DataFrame(rows).to_csv(pdir/pid/'K_diagnostics.tsv',sep='\t',index=False)
            r={'patient_id':pid,'status':'complete','replicates':300,'elapsed_seconds':time.time()-start,'K_diagnostics':rows}
            done.write_text(json.dumps(r,indent=2),encoding='utf-8');return r
        except Exception:
            traceback.print_exc();return {'patient_id':pid,'status':'failed','traceback':traceback.format_exc()}

def main(root):
    out=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26';inp=out/'cache/discovery_patient_inputs';inp.mkdir(exist_ok=True)
    objpath=out/'objects/GSE161529_confirmed_malignant_singlets_counts.h5ad'
    if not objpath.exists():
        md=pd.read_csv(out/'metadata/GSE161529_singlet_metadata.tsv',sep='\t').set_index('cell_id')
        genes=pd.read_csv(out/'metadata/GSE161529_canonical_genes.tsv',sep='\t').set_index('canonical_id')
        with (out/'cache/GSE161529_singlet_counts_gene_by_cell.mtx').open('rb') as f:x=mmread(f).T.tocsr()
        assert x.shape==(len(md),len(genes));x.data=x.data.astype(np.int32)
        obj=ad.AnnData(X=x,obs=md,var=genes);obj.write_h5ad(objpath,compression='gzip')
    else:obj=ad.read_h5ad(objpath)
    availability=pd.read_csv(root/'outputs/恶性上皮_表达程序可行性_2026-10-26/results/患者可用性表.tsv',sep='\t')
    candidates=availability[availability['allocation'].str.startswith(('A_','B_'))].patient_id.tolist()
    targetmask=(obj.var_names=='ENSG00000109736')|obj.var.gene_symbol.isin(['MFSD10','SLC75A1']).to_numpy()
    assert targetmask.sum()==1
    for pid in candidates:
        fn=inp/f'{pid}.h5ad'
        if fn.exists():continue
        pp=obj[obj.obs.patient_id==pid,~targetmask].copy()
        nz=np.asarray((pp.X>0).sum(axis=0)).ravel();keep=nz>=max(5,int(np.ceil(pp.n_obs*.03)))
        pp=pp[:,keep].copy();assert pp.n_vars>=2000 and 'ENSG00000109736' not in pp.var_names
        assert (np.asarray(pp.X.sum(axis=1)).ravel()>0).all()
        pp.write_h5ad(fn,compression='gzip');log('PATIENT_INPUT',pid,pp.shape)
    (out/'metadata/Discovery_patient_fit_list.json').write_text(json.dumps({'patients':candidates,'n':len(candidates),'selection':'Existing A/B cell-count+complexity candidates; stability may fail; no MFSD10 selection','target_excluded':'ENSG00000109736 and MFSD10/SLC75A1 aliases physically removed before normalization and all refits','components':[4,6,8],'n_iter':100,'num_highvar_genes':2000,'max_NMF_iter':1000},indent=2),encoding='utf-8')
    progress=json.loads((out/'audit/PROGRESS.json').read_text());progress.update(discovery_counts_exported=True,stage='per_patient_cnmf_running',discovery_candidates=candidates);(out/'audit/PROGRESS.json').write_text(json.dumps(progress,indent=2),encoding='utf-8')
    results=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(one_patient,(str(out),pid)):pid for pid in candidates}
        for f in concurrent.futures.as_completed(futures):
            r=f.result();results.append(r);log('PATIENT_FACTOR_RESULT',r['patient_id'],r['status'],r.get('elapsed_seconds'))
            (out/'audit/CNMF_patient_run_progress.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    assert all(r['status']=='complete' for r in results), 'One or more patients failed; see individual logs and progress for restart.'
    progress=json.loads((out/'audit/PROGRESS.json').read_text());progress.update(discovery_factorization_complete=True,stage='meta_program_aggregation');(out/'audit/PROGRESS.json').write_text(json.dumps(progress,indent=2),encoding='utf-8')
    log('ALL_PATIENT_FACTORIZATIONS_COMPLETE')

if __name__=='__main__':main(Path(sys.argv[1]))
