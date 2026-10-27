from pathlib import Path
import sys,json
import anndata as ad,pandas as pd
out=Path(sys.argv[1])/'outputs/MFSD10_程序发现与独立验证_2026-10-26'
files=['GSE176078_all_cells_author_counts','GSE176078_malignant_candidates_author_counts']
for name in files:
    target=out/'objects'/f'{name}_clinical_v2.h5ad'
    if target.exists():continue
    obj=ad.read_h5ad(out/'objects'/f'{name}.h5ad',backed='r')
    obj.obs['metadata_subtype_author']=obj.obs['subtype'].astype(str)
    obj.obs['clinical_group']=obj.obs['clinical_subtype_author_IHC'].astype(str).map({'ER+':'ER+','HER2+':'HER2+','HER2+/ER+':'HER2+','TNBC':'TNBC'}).fillna('Unknown')
    obj.obs['clinical_label_conflict']=obj.obs['clinical_group']!=obj.obs['metadata_subtype_author']
    obj.uns['clinical_covariate_source']='Supplementary Table 1 explicit contemporary IHC; metadata subtype separately retained; CID3963 ER+ versus metadata TNBC discrepancy retained'
    obj.write_h5ad(target,compression='gzip');obj.file.close();print('CLINICAL_VERSION_SAVED',name,flush=True)
obj=ad.read_h5ad(out/'objects'/f'{files[1]}_clinical_v2.h5ad',backed='r')
conf=obj.obs.groupby('sample_id',observed=True)[['clinical_group','metadata_subtype_author','clinical_subtype_author_IHC','treatment_status','clinical_label_conflict']].first()
conf.to_csv(out/'audit/GSE176078_clinical_label_reconciliation.tsv',sep='\t')
summary=pd.read_csv(out/'validation/GSE176078_malignant_patient_summary.tsv',sep='\t')
summary['clinical_group']=summary.sample_ids.map(conf['clinical_group'])
summary['metadata_subtype_author']=summary.sample_ids.map(conf['metadata_subtype_author'])
summary['clinical_label_conflict']=summary.sample_ids.map(conf['clinical_label_conflict'])
summary.to_csv(out/'validation/GSE176078_malignant_patient_summary_clinical_v2.tsv',sep='\t',index=False)
obj.file.close()
p=out/'audit/PROGRESS.json';d=json.loads(p.read_text());d.update(validation_clinical_reconciled=True,validation_primary_object='objects/GSE176078_malignant_candidates_author_counts_clinical_v2.h5ad');p.write_text(json.dumps(d,indent=2),encoding='utf-8')
