import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,re,warnings,yaml
import numpy as np,pandas as pd,anndata as ad
from cnmf import cNMF
from cnmf.cnmf import load_df_from_npz
from sklearn.decomposition import non_negative_factorization
from scipy.optimize import linear_sum_assignment
out=Path(sys.argv[1])/'outputs/MFSD10_程序发现与独立验证_2026-10-26';rows=[]
for log in sorted((out/'logs').glob('cnmf_Cancer__*.log')):
    pid=log.stem.removeprefix('cnmf_');task=None;flag=[]
    for line in log.read_text(encoding='utf-8',errors='replace').splitlines():
        m=re.search(r'Starting task (\d+)',line)
        if m:task=int(m.group(1))
        if 'ConvergenceWarning' in line and task is not None:flag.append(task)
    if not flag:continue
    c=cNMF(output_dir=str(out/'discovery/cnmf'),name=pid);params=load_df_from_npz(c.paths['nmf_replicate_parameters'])
    with open(c.paths['nmf_run_parameters']) as f:base=yaml.safe_load(f)
    obj=ad.read_h5ad(c.paths['normalized_counts'])
    for ix in sorted(set(flag)):
        p=params.iloc[ix];kw=dict(base);kw.update(n_components=int(p.n_components),random_state=int(p.nmf_seed),max_iter=3000)
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always');U,H,it=non_negative_factorization(obj.X,**kw)
        before=load_df_from_npz(c.paths['iter_spectra']%(p.n_components,p['iter'])).to_numpy()
        a=before/np.maximum(np.linalg.norm(before,axis=1,keepdims=True),1e-12);b=H/np.maximum(np.linalg.norm(H,axis=1,keepdims=True),1e-12);cos=a@b.T;ii,jj=linear_sum_assignment(-cos)
        rows.append(dict(patient_id=pid,task=ix,K=int(p.n_components),seed=int(p.nmf_seed),original_max_iter=1000,followup_max_iter=3000,followup_n_iter=int(it),followup_limit_reached=it>=3000,median_component_cosine=float(np.median(cos[ii,jj])),min_component_cosine=float(np.min(cos[ii,jj])),note='Only explicitly logged primary-run limits revisited; warnings may be suppressed across workers, so this is not an exhaustive convergence count. Original/frozen files unchanged.'))
pd.DataFrame(rows).to_csv(out/'discovery/Convergence_limit_followup.tsv',sep='\t',index=False);print(pd.DataFrame(rows).to_string(index=False))
