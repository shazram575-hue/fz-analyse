import os
os.environ.setdefault('OMP_NUM_THREADS','2'); os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
from pathlib import Path
import sys,json,gc,hashlib,time
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.io import mmread
import anndata as ad

root=Path(sys.argv[1]);out=root/'outputs/MFSD10_程序发现与独立验证_2026-10-26';data=out/'cache/GSE176078_unpacked/Wu_etal_2021_BRCA_scRNASeq'
def status(**kw):
    p=out/'audit/PROGRESS.json';d=json.loads(p.read_text(encoding='utf-8-sig'));d.update(kw);p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
def log(*s):print(time.strftime('%H:%M:%S'),*s,flush=True)
bars=pd.read_csv(data/'count_matrix_barcodes.tsv',header=None,sep='\t',dtype=str)[0].to_numpy()
genes=pd.read_csv(data/'count_matrix_genes.tsv',header=None,sep='\t',dtype=str)[0].to_numpy()
md=pd.read_csv(data/'metadata.csv',index_col=0);assert md.index.is_unique and pd.Index(bars).is_unique
assert set(md.index)==set(bars);md=md.loc[bars].copy()
assert 'orig.ident' in md and 'celltype_major' in md
log('Reading sparse matrix; no full dense conversion')
with (data/'count_matrix_sparse.mtx').open('rb') as stream:
    mm=mmread(stream) # Binary file handle avoids Windows Unicode-path limitations in the C++ reader.
if mm.shape==(len(genes),len(bars)): orientation='genes_by_cells';mm=mm.T
elif mm.shape==(len(bars),len(genes)):orientation='cells_by_genes'
else:raise ValueError('Matrix dimensions do not match actual gene/barcode tables')
assert np.isfinite(mm.data).all() and (mm.data>=0).all() and (mm.data==np.rint(mm.data)).all()
assert mm.data.max()<np.iinfo(np.int32).max
mm.data=mm.data.astype(np.int32);X=mm.tocsr();X.sum_duplicates();X.eliminate_zeros();del mm;gc.collect()
log('Sparse matrix ready',X.shape,X.nnz)
ncount=np.asarray(X.sum(axis=1)).ravel();nfeature=np.diff(X.indptr)
assert np.array_equal(ncount,md['nCount_RNA'].to_numpy())
assert np.array_equal(nfeature,md['nFeature_RNA'].to_numpy())
gm=pd.read_csv(out/'metadata/GSE161529_canonical_genes.tsv',sep='\t')
valid=gm.dropna(subset=['gene_symbol']).groupby('gene_symbol')['canonical_id'].agg(list).to_dict()
unified=[];mapstatus=[]
for g in genes:
    hits=valid.get(g,[])
    if len(hits)==1:unified.append(hits[0]);mapstatus.append('exact_unique_symbol_to_GENCODEv36')
    else:unified.append('AUTHOR_SYMBOL::'+g);mapstatus.append('ambiguous' if len(hits)>1 else 'unmapped_preserved')
var=pd.DataFrame({'author_gene':genes,'unified_gene_id':unified,'mapping_status':mapstatus},index=pd.Index([f'author_feature_{i:05d}' for i in range(len(genes))],name='feature_id'))
target=np.flatnonzero((genes=='MFSD10')|(genes=='SLC75A1')|(np.array(unified)=='ENSG00000109736'))
assert len(target)==1, 'MFSD10 mapping ambiguous; manual review required'
var.iloc[target].to_csv(out/'metadata/GSE176078_MFSD10_mapping.tsv',sep='\t')
var.to_csv(out/'metadata/GSE176078_author_to_unified_gene_map.tsv',sep='\t')
# The actual clinical table header is on row four, identified during inventory.
clinical=pd.read_csv(out/'metadata/supplement_Supplementary Table 1.tsv',sep='\t',skiprows=3,dtype=str).dropna(subset=['Case ID'])
assert clinical['Case ID'].is_unique
alias={'CID4290A':'4290','CID44041':'4404-1','CID44971':'4497-1','CID44991':'4499-1','CID45171':'4517-1','CID4530N':'4530'}
records=[]
for sid in sorted(md['orig.ident'].unique()):
    case=alias.get(sid,sid.removeprefix('CID'));hit=clinical[clinical['Case ID']==case];assert len(hit)==1,(sid,case)
    r=hit.iloc[0].to_dict();r.update(sample_id=sid,patient_id='GSE176078_Case__'+case,mapping_basis='Unique author clinical Case ID; explicit suffix/hyphen reconciliation' if sid in alias else 'Exact author clinical Case ID')
    records.append(r)
pt=pd.DataFrame(records);assert len(pt)==26 and pt.patient_id.nunique()==26
pt['treatment_status']=pt['Treatment status'].map({'Naïve':'Naive','Treated':'Treated'}).fillna('Unknown')
pt['clinical_subtype_author_IHC']=pt['Subtype by IHC'];pt['tissue_origin']='Primary_breast_tumor';pt['tissue_source_evidence']='GEO series and sample tissue characteristic; treated/recurrent details retained separately'
pt.to_csv(out/'metadata/GSE176078_patient_sample_clinical_map.tsv',sep='\t',index=False)
lookup=pt.set_index('sample_id')
for c in ['patient_id','treatment_status','clinical_subtype_author_IHC','tissue_origin','mapping_basis']:
    md[c]=md['orig.ident'].map(lookup[c]);assert md[c].notna().all()
md['sample_id']=md['orig.ident'];md['metadata_subtype_author']=md['subtype']
md['clinical_group']=md['clinical_subtype_author_IHC'].map({'ER+':'ER+','HER2+':'HER2+','HER2+/ER+':'HER2+','TNBC':'TNBC'}).fillna('Unknown')
md['clinical_label_conflict']=md['clinical_group']!=md['metadata_subtype_author']
md['author_celltype_major']=md['celltype_major'];md['author_celltype_minor']=md['celltype_minor'];md['author_celltype_subset']=md['celltype_subset']
md['malignancy_status']=np.where(md['celltype_major']=='Cancer Epithelial','Author_neoplastic_candidate',np.where(md['celltype_major']=='Normal Epithelial','Author_normal_epithelial','Not_epithelial'))
md['doublet_status']='Unknown_not_in_published_metadata'
md['RNA_counts_total']=ncount;md['RNA_counts_features']=nfeature
mt=np.flatnonzero(np.char.startswith(genes.astype(str),'MT-'));mtp=100*np.asarray(X[:,mt].sum(axis=1)).ravel()/ncount
md['mt_pct_recomputed']=mtp
md['existing_author_QC_consistent']=(ncount>250)&(nfeature>200)&(mtp<=20+1e-7)
mf=np.asarray(X[:,target[0]].toarray()).ravel();md['MFSD10_counts']=mf
md['RNA_UMI_excluding_MFSD10']=ncount-mf
md['MFSD10_log1p_CP10K']=np.log1p(mf/np.maximum(ncount-mf,1)*1e4)
marker_names=['EPCAM','KRT8','KRT18','KRT19','KRT5','KRT14','PTPRC','CD3D','LYZ','LST1','COL1A1','DCN','VWF','PECAM1']
marker_counts={g:np.asarray(X[:,np.flatnonzero(genes==g)].sum(axis=1)).ravel() for g in marker_names}
md['lineage_mixture_flag']=(marker_counts['PTPRC']>0)&((marker_counts['CD3D']>0)|(marker_counts['LST1']>0))&(md['celltype_major']=='Cancer Epithelial')
# Flags are descriptive; co-detection is not sufficient to call a doublet or revoke author CNV status.
marker_rows=[]
for ct,idx in md.groupby('celltype_major',observed=True).indices.items():
    for g,cnt in marker_counts.items():marker_rows.append(dict(author_celltype=ct,gene=g,n_cells=len(idx),pct_detected=100*np.mean(cnt[idx]>0),mean_log1p_CP10K=np.log1p(cnt[idx]/ncount[idx]*1e4).mean()))
pd.DataFrame(marker_rows).to_csv(out/'audit/GSE176078_lineage_marker_check.tsv',sep='\t',index=False)
q=[]
for (patient,ct),d in md.groupby(['patient_id','malignancy_status'],observed=True):
    q.append(dict(patient_id=patient,status=ct,n_cells=len(d),n_authorQC_inconsistent=int((~d.existing_author_QC_consistent).sum()),n_lineage_flags=int(d.lineage_mixture_flag.sum()),UMI_median=float(d.nCount_RNA.median()),genes_median=float(d.nFeature_RNA.median()),mt_pct_median=float(d.mt_pct_recomputed.median())))
pd.DataFrame(q).to_csv(out/'audit/GSE176078_existing_quality_distribution.tsv',sep='\t',index=False)
allobj=ad.AnnData(X=X,obs=md,var=var)
allobj.uns['provenance']={'matrix_orientation_detected':orientation,'author_raw_gene_names_preserved':True,'integrated_expression_used':False,'author_malignancy_basis':'inferCNV with immune/endothelial references, genomic instability and CNV profile correlation; see local methods PDF','doublet_handling':'Not specified in supplied metadata/methods; no new algorithm applied','original_processed_input':'GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz'}
log('Writing all-cell reference object')
allobj.write_h5ad(out/'objects/GSE176078_all_cells_author_counts.h5ad',compression='gzip')
candidate=md['malignancy_status'].eq('Author_neoplastic_candidate').to_numpy()
# Only reject demonstrable author-rule inconsistencies, preserving them in the reference object and QC table.
eligible=candidate&md.existing_author_QC_consistent.to_numpy()
malobj=allobj[eligible].copy();malobj.uns['candidate_count_before_QC']=int(candidate.sum());malobj.uns['candidate_count_after_QC']=int(eligible.sum())
malobj.write_h5ad(out/'objects/GSE176078_malignant_candidates_author_counts.h5ad',compression='gzip')
sm=[]
for pat,d in md[eligible].groupby('patient_id',observed=True):
    sm.append(dict(patient_id=pat,sample_ids=';'.join(sorted(d.sample_id.unique())),treatment_status=d.treatment_status.iloc[0],clinical_group=d.clinical_group.iloc[0],clinical_subtype_author_IHC=d.clinical_subtype_author_IHC.iloc[0],n_malignant=len(d),UMI_median=float(d.nCount_RNA.median()),genes_median=float(d.nFeature_RNA.median()),MFSD10_n_detected=int((d.MFSD10_counts>0).sum()),MFSD10_pct_detected=float(100*(d.MFSD10_counts>0).mean()),MFSD10_CPM=float(1e6*d.MFSD10_counts.sum()/d.RNA_UMI_excluding_MFSD10.sum()),lineage_flags=int(d.lineage_mixture_flag.sum()),doublet_status='Unknown_author_metadata_missing'))
pd.DataFrame(sm).to_csv(out/'validation/GSE176078_malignant_patient_summary.tsv',sep='\t',index=False)
md.to_csv(out/'metadata/GSE176078_all_cells_unified_metadata.tsv.gz',sep='\t')
audit={'shape_cells_genes':list(X.shape),'nnz':int(X.nnz),'orientation_detected':orientation,'matrix_nonnegative_integer':True,'barcodes_unique':True,'metadata_one_to_one':True,'UMI_matches_author_metadata':True,'nFeature_matches_author_metadata':True,'MFSD10_author_gene':str(genes[target[0]]),'MFSD10_n_feature_rows':1,'author_gene_duplicates':int(pd.Index(genes).duplicated().sum()),'Cancer_Epithelial_before_QC':int(candidate.sum()),'Cancer_Epithelial_after_QC':int(eligible.sum()),'QC_inconsistent_all_cells':int((~md.existing_author_QC_consistent).sum()),'all_cells':len(md),'patients_clinical_table':len(pt),'naive_patients':int((pt.treatment_status=='Naive').sum()),'treated_patients':int((pt.treatment_status=='Treated').sum())}
(out/'audit/GSE176078_import_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
status(validation_object_saved=True,validation_import_audit=audit,stage='discovery_fitting_preparation');log('VALIDATION_OBJECTS_COMPLETE',audit)
