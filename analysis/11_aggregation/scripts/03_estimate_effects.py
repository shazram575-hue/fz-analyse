import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from scipy.stats import rankdata,t as student_t
from reused_original_models import resid,partial_rho,sign_meta,bh

def covariates(d,cols):
    ans=[np.ones(len(d))]
    for col in cols:
        v=d[col].to_numpy(float);v=np.log1p(v) if col!='mt_pct' else v
        v=rankdata(v);v-=v.mean()
        if v.std()>0:ans.append(v/v.std())
    return np.column_stack(ans)
def between(d,mp,adj):
    n=len(d);cols=[np.ones(n)]
    if adj!='unadjusted':cols.extend(pd.get_dummies(d.clinical_group,drop_first=True,dtype=float).to_numpy().T)
    if adj=='IHC_depth':cols.append(rankdata(np.log1p(d.median_depth)))
    C=np.column_stack(cols);rank=np.linalg.matrix_rank(C);df=n-rank-1
    rho=partial_rho(d.MFSD10_log2_CPM.to_numpy(),d[mp].to_numpy(),C);ok=n>=12 and df>=5 and np.isfinite(rho)
    p=2*student_t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-12)),df) if ok else np.nan
    se=1/np.sqrt(max(n-rank-2,1));ci=np.tanh(np.arctanh(np.clip(rho,-.999,.999))+np.array([-1,1])*1.96*se) if ok else [np.nan,np.nan]
    return dict(program_id=mp,specification=adj,n_patients=n,effect_rho=rho,ci_low=ci[0],ci_high=ci[1],p=p,residual_df=df,analysis_unit='patient',note='Exploratory fixed parsimonious sensitivity; no cell-level sample size')

root=Path(sys.argv[1]);old=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26';out=root/'outputs/MFSD10_稀疏与聚合探索_2026-10-26'
programs=[m['program_id'] for m in json.loads((old/'frozen/Programs_v1.json').read_text())['meta_programs']]
specs={
 'state_full':('state',['total_without_target','features_without_target','mt_pct'],'countscore'),
 'state_unadjusted':('state',[],'countscore'),
 'state_depth':('state',['total_without_target'],'countscore'),
 'state_depth_mt':('state',['total_without_target','mt_pct'],'countscore'),
 'state_cellmean_full':('state',['total_without_target','features_without_target','mt_pct'],'cellmean'),
 'random_full':('random',['total_without_target','features_without_target','mt_pct'],'countscore')}
effect_rows=[];patient_rows=[];paired_rows=[];matched_rows=[];between_rows=[]
rng=np.random.default_rng(20261026)
for cohort in ['discovery','validation']:
    availability=pd.read_csv(out/'results'/f'{cohort}_patient_pseudobulk_audit.tsv',sep='\t');primary=availability.loc[availability.primary_patient.eq(True),'patient_id'].tolist()
    frames={k:pd.read_csv(out/'results'/f'{cohort}_{k}_unit_summary.tsv.gz',sep='\t',index_col=0) for k in ['state','random']}
    old_per=pd.read_csv(old/('discovery/MFSD10_within_patient_effects.tsv' if cohort=='discovery' else 'validation/MFSD10_within_patient_effects_Naive.tsv'),sep='\t')
    for name,(partition,cols,scoring) in specs.items():
        u=frames[partition];local=[]
        for pid,d in u.groupby('patient_id',observed=True):
            C=covariates(d,cols);rank=np.linalg.matrix_rank(C);df=len(d)-rank-1
            reason='Fewer_than_8_units' if len(d)<8 else 'Residual_df_less_than_4' if df<4 else ''
            for mp in programs:
                scorecol=mp if scoring=='countscore' else 'cellmean_'+mp
                rho=np.nan
                if not reason:
                    rho=partial_rho(d.MFSD10_log1p_CP10K.to_numpy(float),d[scorecol].to_numpy(float),C)
                why=reason or ('Constant_or_degenerate_residual_variable' if not np.isfinite(rho) else '')
                row=dict(cohort=cohort,specification=name,patient_id=pid,program_id=mp,primary_patient=pid in primary,n_units=len(d),n_cells=int(d.n_cells.sum()),nuisance_rank=rank,residual_df=df,rho=rho,reason_not_estimated=why)
                patient_rows.append(row);local.append(row)
        pp=pd.DataFrame(local)
        for mp in programs:
            rr=pp[(pp.program_id==mp)&pp.primary_patient&pp.rho.notna()]
            result=sign_meta(rr.rho,seed=20261026);result.update(cohort=cohort,program_id=mp,specification=name,n_units_total=int(rr.n_units.sum()),n_cells_total=int(rr.n_cells.sum()),analysis_type='Exploratory_reanalysis',analysis_unit='patient_effect_not_unit',adjustment=';'.join(cols) or 'none')
            effect_rows.append(result)
            if name=='state_full':
                ids=rr.patient_id.tolist();prev=old_per[(old_per.program_id==mp)&old_per.patient_id.isin(ids)].set_index('patient_id')
                aa=rr.set_index('patient_id').loc[prev.index,'rho'].to_numpy();bb=prev.rho.to_numpy();ok=np.isfinite(aa)&np.isfinite(bb);aa=aa[ok];bb=bb[ok]
                match=sign_meta(bb,seed=20261026);match.update(cohort=cohort,program_id=mp,note='Reuses saved original patient rhos, restricts to aggregate-estimable patients; no single-cell model refit');matched_rows.append(match)
                if len(aa)>=6:
                    take=rng.integers(0,len(aa),size=(5000,len(aa)));za=np.arctanh(np.clip(aa,-.999,.999));zb=np.arctanh(np.clip(bb,-.999,.999));delta=np.tanh(za[take].mean(axis=1))-np.tanh(zb[take].mean(axis=1));lo,hi=np.quantile(delta,[.025,.975]);point=np.tanh(za.mean())-np.tanh(zb.mean())
                else:point=lo=hi=np.nan
                paired_rows.append(dict(cohort=cohort,program_id=mp,n_matched_patients=len(aa),delta_meta_rho=point,ci_low=lo,ci_high=hi,note='Paired patient bootstrap; aggregation changes resolution/estimand, not a causal effect'))
    b=availability[availability.primary_patient.eq(True)].copy()
    for adj in ['unadjusted','IHC_only','IHC_depth']:
        for mp in programs:r=between(b,mp,adj);r['cohort']=cohort;between_rows.append(r)

effects=pd.DataFrame(effect_rows)
effects['q_specification_10']=np.nan
for spec,idx in effects.groupby('specification').groups.items():effects.loc[idx,'q_specification_10']=bh(effects.loc[idx,'p'])
effects['q_all60_exploratory']=bh(effects.p)
main=effects[effects.specification=='state_full'].copy();main['FDR_main10']=bh(main.p)
main['FDR_within_cohort5']=np.nan
for coh,idx in main.groupby('cohort').groups.items():main.loc[idx,'FDR_within_cohort5']=bh(main.loc[idx,'p'])
effects.to_csv(out/'results/Aggregation_all60_fixed_specifications.tsv',sep='\t',index=False)
main.to_csv(out/'results/Aggregation_main_all10.tsv',sep='\t',index=False)
pd.DataFrame(patient_rows).to_csv(out/'results/Aggregation_patient_effects_and_failures.tsv',sep='\t',index=False)
pd.DataFrame(paired_rows).to_csv(out/'results/Matched_patient_delta_vs_singlecell.tsv',sep='\t',index=False)
matched=pd.DataFrame(matched_rows);matched['FDR_matched10']=bh(matched.p);matched.to_csv(out/'results/Matched_original_singlecell_reused.tsv',sep='\t',index=False)
b=pd.DataFrame(between_rows);b['q_all30_sensitivities']=bh(b.p);b['q_specification10']=np.nan
for spec,idx in b.groupby('specification').groups.items():b.loc[idx,'q_specification10']=bh(b.loc[idx,'p'])
b.to_csv(out/'results/Patient_between_parsimonious_sensitivities.tsv',sep='\t',index=False)
p=out/'audit/PROGRESS.json';s=json.loads(p.read_text());s.update(patient_effects_complete=True,stage='comparison_and_report');p.write_text(json.dumps(s,indent=2),encoding='utf-8')
print(main[['cohort','program_id','n_patients','effect_rho','ci_low','ci_high','p','FDR_main10','q_all60_exploratory']].to_string(index=False))
