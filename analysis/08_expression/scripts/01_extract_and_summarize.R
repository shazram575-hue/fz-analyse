# Read-only reuse of frozen annotations/counts; no Seurat reprocessing or CNV projection.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c('D:/R4.5.2/R-4.5.2/library',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(Seurat);library(edgeR)})
setDTthreads(2L);set.seed(20261026)
root<-normalizePath(commandArgs(TRUE)[1],winslash='/');old<-file.path(root,'outputs/第一阶段_2026-10-23')
base<-file.path(root,'outputs/单细胞基础补全_2026-10-23');out<-file.path(root,'outputs/MFSD10_表达分布摸底_2026-10-26')
for(d in c('metadata','results','figures','logs','audit','cache/r_tmp','cache/by_library'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
sink(file.path(out,'logs/01_extract.log'),split=TRUE)
map<-fread(file.path(old,'metadata/SingleCell_canonical_gene_identity_map.tsv'))
genes<-fread(file.path(old,'metadata/SingleCell_canonical_gene_metadata.tsv'))
target<-genes[gene_symbol=='MFSD10'];stopifnot(nrow(target)==1L,target$canonical_id=='ENSG00000109736')
fwrite(map[canonical_id==target$canonical_id],file.path(out,'metadata/MFSD10_gene_mapping.tsv'),sep='\t')
collapse<-sparseMatrix(i=match(map$canonical_id,genes$canonical_id),j=seq_len(nrow(map)),x=1,dims=c(nrow(genes),nrow(map)))
m<-fread(file.path(out,'metadata/Library_inventory.tsv'),colClasses=c(patient='character'))
m[,tissue_context:=fcase(analysis_role=='lymph_node_separate','Lymph_node',analysis_role=='preneoplastic_separate','Preneoplastic',cancer.type=='Normal','Healthy_breast',default='Primary_tumor')]
m[,preparation:=fcase(cell.population=='Epithelial','Epithelial_enriched',cell.population=='Involved lymph node','LN_unspecified',default='Total')]
m[,clinical_label:=fcase(cancer.type=='Normal','Normal',cancer.type=='ER+ tumour','ER+',cancer.type=='PR+ tumour','PR+',cancer.type=='HER2+ tumour','HER2+',cancer.type=='Triple negative tumour','TNBC',cancer.type=='Triple negative BRCA1 tumour','TNBC_BRCA1',default='BRCA1_preneoplastic')]
# Exact specimen joins to archived clinical tables; no receptor/subtype imputation.
ref<-file.path(old,'metadata/author_reference');clean<-function(z)gsub('[[:space:]]+','',z)
er<-fread(file.path(ref,'EMBJ-40-e107333-s001.docx_table1.tsv'),skip=1,header=FALSE,fill=TRUE,sep='\t')
tn<-fread(file.path(ref,'EMBJ-40-e107333-s005.docx_table1.tsv'),skip=1,header=FALSE,fill=TRUE,sep='\t')
clin<-rbindlist(list(data.table(clinical_specimen=clean(er[[1]]),published_tumor_type=er[[3]],histology=er[[4]],ER_text=er[[5]],grade_text=er[[6]],size_mm_text=er[[7]],BRCA1_variant=NA_character_,HER2_text=NA_character_,clinical_source='Pal2021_EV3'),data.table(clinical_specimen=clean(tn[[1]]),published_tumor_type=tn[[3]],histology=NA_character_,ER_text=NA_character_,grade_text=tn[[5]],size_mm_text=tn[[4]],BRCA1_variant=tn[[6]],HER2_text=tn[[7]],clinical_source='Pal2021_EV2')),fill=TRUE)
stopifnot(!anyDuplicated(clin$clinical_specimen));m<-merge(m,clin,by='clinical_specimen',all.x=TRUE,sort=FALSE);setorder(m,gsm)
fwrite(m,file.path(out,'metadata/Library_clinical_context.tsv'),sep='\t')
cnv<-fread(file.path(old,'results/Tumor_discovery_cells_CNV_coherent_final_calls.tsv'))
stopifnot(!anyDuplicated(cnv$cell_id),nrow(cnv)==29113,sum(cnv$cnv_high_q99)==9486)
kept<-readRDS(file.path(old,'cache/GSE161529_author_QC_kept_barcodes_by_gsm.rds'))
pairs<-fread(file.path(old,'audit/GSE161529_tar_file_pairs.tsv'))
qc<-as.data.table(read.delim(gzfile(file.path(old,'cache/GSE161529_all_cell_QC_membership.tsv.gz')),check.names=FALSE))
sources<-c(file.path(old,c('metadata/SingleCell_canonical_gene_identity_map.tsv','metadata/SingleCell_canonical_gene_metadata.tsv','metadata/GSE161529_reconciled_patient_metadata.tsv','results/Tumor_discovery_cells_CNV_coherent_final_calls.tsv','cache/GSE161529_author_QC_kept_barcodes_by_gsm.rds','audit/GSE161529_tar_file_pairs.tsv','cache/GSE161529_all_cell_QC_membership.tsv.gz')),file.path(ref,c('EMBJ-40-e107333-s001.docx_table1.tsv','EMBJ-40-e107333-s005.docx_table1.tsv')))
allcells<-list();pbparts<-list();libaudit<-list()
for(i in seq_len(nrow(m))){
 g<-m$gsm[i];cat('LIBRARY',i,'/69',g,format(Sys.time()),'\n');flush.console()
 af<-file.path(old,'audit/cell_annotation',paste0(g,'_cell_identity.tsv'));cf<-file.path(base,'cache/canonical_by_library',paste0(g,'.rds'));sf<-file.path(old,'cache/seurat_pilot',paste0(g,'_pilot_seurat.rds'))
 if(file.exists(af)){
  md<-fread(af);sources<-c(sources,af)
  if(file.exists(cf)){
   x<-readRDS(cf);sources<-c(sources,cf);raw_source<-'RNA_counts_canonical_cache';object_fields<-'cached canonical RNA counts; annotations in cell_identity.tsv'
  }else{
   so<-readRDS(sf);x<-GetAssayData(so,assay='RNA',layer='counts');stopifnot(identical(rownames(x),map$ensembl_id))
   object_fields<-paste(names(so@meta.data),collapse=';');x<-collapse%*%x;rownames(x)<-genes$canonical_id;rm(so);gc(FALSE)
   sources<-c(sources,sf);raw_source<-'Seurat_RNA_counts'
  }
  md<-md[match(colnames(x),cell_id)];stopifnot(identical(md$cell_id,colnames(x)))
 }else{
  pp<-pairs[gsm==g];mf<-file.path(old,'cache/GSE161529_unpacked',pp$matrix_member);bf<-file.path(old,'cache/GSE161529_unpacked',pp$barcode_member)
  con<-gzfile(mf,'rt');raw<-as(readMM(con),'CsparseMatrix');close(con);bc<-readLines(gzfile(bf),warn=FALSE)
  stopifnot(nrow(raw)==nrow(map),length(bc)==ncol(raw));colnames(raw)<-bc;raw<-raw[,kept[[g]],drop=FALSE];colnames(raw)<-paste(g,kept[[g]],sep='__')
  x<-collapse%*%raw;rownames(x)<-genes$canonical_id;rm(raw);gc(FALSE)
  md<-qc[match(colnames(x),cell_id)];stopifnot(identical(md$cell_id,colnames(x)),all(md$pass_author_qc))
  md[,`:=`(cell_type='Not_annotated',scDblFinder_class=NA_character_)];sources<-c(sources,mf,bf);raw_source<-'raw_matrix_author_QC';object_fields<-'No existing Seurat object or cell annotation'
 }
 stopifnot(identical(rownames(x),genes$canonical_id),all(abs(Matrix::colSums(x)-md$nCount_RNA)<1e-6),all(x@x>=0),all(x@x==round(x@x)))
 count<-as.numeric(x[target$canonical_id,]);raw_umi<-Matrix::colSums(x)
 cc<-data.table(cell_id=md$cell_id,gsm=g,patient_id=m$patient_id[i],clinical_specimen=m$clinical_specimen[i],clinical_label=m$clinical_label[i],tissue_context=m$tissue_context[i],preparation=m$preparation[i],cell_type=md$cell_type,nCount_RNA=md$nCount_RNA,nFeature_RNA=md$nFeature_RNA,scDblFinder_class=md$scDblFinder_class,MFSD10_counts=count,MFSD10_detected=count>0,MFSD10_log1p_CP10K=log1p(count/raw_umi*10000),RNA_total_UMI=raw_umi)
 cc[,epithelial_status:='Not_epithelial_or_unannotated']
 cc[cell_type=='Epithelial' & tissue_context=='Healthy_breast',epithelial_status:='Healthy_context_no_CNV_call']
 cc[cell_type=='Epithelial' & tissue_context=='Primary_tumor',epithelial_status:='Tumor_epithelial_unscored']
 ix<-match(cc$cell_id,cnv$cell_id);has<-which(!is.na(ix));stopifnot(all(cc$cell_type[has]=='Epithelial'),all(cc$patient_id[has]==cnv$patient_id[ix[has]]))
 cc[has,epithelial_status:=ifelse(cnv$cnv_high_q99[ix[has]],'Malignant_CNV_coherent','Uncertain_CNV_scored')]
 cc[,`:=`(nCount_without_MFSD10=nCount_RNA-MFSD10_counts,nFeature_without_MFSD10=nFeature_RNA-as.integer(MFSD10_detected))]
 # Two disjoint partitions separately: broad lineage, then epithelial CNV/context state.
 memb<-rbindlist(list(cc[,.(cell_index=seq_len(.N),partition='broad',cell_type,epithelial_status='All')],cc[cell_type=='Epithelial',.(cell_index=which(cc$cell_type=='Epithelial'),partition='epithelial',cell_type,epithelial_status)]))
 memb[,local_key:=paste(partition,cell_type,epithelial_status,sep='::')];keys<-unique(memb$local_key)
 D<-sparseMatrix(i=memb$cell_index,j=match(memb$local_key,keys),x=1,dims=c(nrow(cc),length(keys)))
 pb<-as.matrix(x%*%D);colnames(pb)<-paste(g,keys,sep='::');rownames(pb)<-genes$canonical_id
 sm<-memb[,.(n_cells=.N,n_expressing=sum(cc$MFSD10_detected[cell_index]),sum_log1p_CP10K=sum(cc$MFSD10_log1p_CP10K[cell_index]),MFSD10_UMI=sum(cc$MFSD10_counts[cell_index]),RNA_UMI=sum(cc$RNA_total_UMI[cell_index]),n_doublets=sum(cc$scDblFinder_class[cell_index]=='doublet',na.rm=TRUE)),by=.(local_key,partition,cell_type,epithelial_status)]
 sm[,`:=`(column_id=paste(g,local_key,sep='::'),gsm=g,patient_id=m$patient_id[i],clinical_specimen=m$clinical_specimen[i],clinical_label=m$clinical_label[i],tissue_context=m$tissue_context[i],preparation=m$preparation[i])]
 stopifnot(all(sm$MFSD10_UMI==pb[target$canonical_id,sm$column_id]),all(sm$RNA_UMI==colSums(pb[,sm$column_id,drop=FALSE])))
 saveRDS(list(counts=pb,samples=sm),file.path(out,'cache/by_library',paste0(g,'.rds')),compress=FALSE)
 allcells[[g]]<-cc;pbparts[[g]]<-list(counts=pb,samples=sm)
 libaudit[[g]]<-data.table(gsm=g,n_cells=nrow(cc),n_features=nrow(x),MFSD10_UMI=sum(count),RNA_UMI=sum(raw_umi),source=raw_source,object_metadata_fields=object_fields,n_CNV_scored=length(has))
 rm(x,md,pb,cc,D);gc(FALSE)
}
cells<-rbindlist(allcells);rm(allcells);gc(FALSE)
stopifnot(nrow(cells)==340979L,!anyDuplicated(cells$cell_id),cells[epithelial_status=='Malignant_CNV_coherent',.N]==9486,cells[epithelial_status=='Uncertain_CNV_scored',.N]==19627)
saveRDS(cells,file.path(out,'cache/MFSD10_cell_level.rds'),compress=FALSE)
fwrite(cells,file.path(out,'results/MFSD10_cell_level.tsv.gz'),sep='\t')
fwrite(rbindlist(libaudit),file.path(out,'audit/Library_count_checks.tsv'),sep='\t')
sm<-rbindlist(lapply(pbparts,`[[`,'samples'));bigpb<-do.call(cbind,lapply(pbparts,`[[`,'counts'));rm(pbparts);gc(FALSE)
sm[,`:=`(pct_detected=100*n_expressing/n_cells,CPM=1e6*MFSD10_UMI/RNA_UMI,mean_log1p_CP10K=sum_log1p_CP10K/n_cells,descriptive_min20=n_cells>=20)]
fwrite(sm,file.path(out,'results/Library_celltype_and_epithelial_summary.tsv'),sep='\t')
aggregate_pb<-function(part){
 s<-copy(sm[partition==part]);byvars<-c('patient_id','tissue_context','preparation','cell_type','epithelial_status');if(part=='epithelial')byvars<-c(byvars,'clinical_label')
 s[,patient_group:=do.call(paste,c(.SD,sep='::')),.SDcols=byvars];keys<-unique(s$patient_group)
 D<-sparseMatrix(i=seq_len(nrow(s)),j=match(s$patient_group,keys),x=1,dims=c(nrow(s),length(keys)))
 pb<-bigpb[,s$column_id,drop=FALSE]%*%D;pb<-as.matrix(pb);colnames(pb)<-keys
 a<-s[,.(n_cells=sum(n_cells),n_expressing=sum(n_expressing),MFSD10_UMI=sum(MFSD10_UMI),RNA_UMI=sum(RNA_UMI),sum_log1p_CP10K=sum(sum_log1p_CP10K),n_doublets=sum(n_doublets),n_libraries=uniqueN(gsm),GSMs=paste(sort(unique(gsm)),collapse=';'),specimens=paste(sort(unique(clinical_specimen)),collapse=';'),clinical_labels=paste(sort(unique(clinical_label)),collapse=';')),by=c('patient_group',byvars)]
 a<-a[match(keys,patient_group)];stopifnot(all(a$MFSD10_UMI==pb[target$canonical_id,]),all(a$RNA_UMI==colSums(pb)))
 a[,`:=`(pct_detected=100*n_expressing/n_cells,CPM=1e6*MFSD10_UMI/RNA_UMI,mean_log1p_CP10K=sum_log1p_CP10K/n_cells,descriptive_min20=n_cells>=20,log2_CPM_plus1=log2(1+1e6*MFSD10_UMI/RNA_UMI))]
 a[,`:=`(TMM_factor=NA_real_,TMM_CPM=NA_real_,TMM_eligible=FALSE)]
 # Within lineage x anatomical site x preparation; small groups retained descriptively.
 a[,normalization_stratum:=paste(cell_type,ifelse(tissue_context=='Lymph_node','LN','Breast'),preparation,sep='::')]
 norms<-list()
 for(st in unique(a$normalization_stratum)){
  j<-which(a$normalization_stratum==st & a$n_cells>=20);if(length(j)<2)next
  y<-calcNormFactors(DGEList(pb[,j,drop=FALSE]),method='TMM');a[j,`:=`(TMM_factor=y$samples$norm.factors,TMM_CPM=as.numeric(cpm(y,log=FALSE)[target$canonical_id,]),TMM_eligible=TRUE)]
  norms[[st]]<-data.table(stratum=st,n_groups=length(j),n_unique_patients=uniqueN(a$patient_id[j]),zero_fraction=mean(pb[,j,drop=FALSE]==0))
 }
 saveRDS(list(counts=pb,samples=a,genes=genes),file.path(out,'cache',paste0('Patient_',part,'_pseudobulk_counts.rds')),compress=FALSE)
 fwrite(a,file.path(out,'results',paste0('Patient_',part,'_MFSD10_summary.tsv')),sep='\t')
 fwrite(rbindlist(norms),file.path(out,'audit',paste0('TMM_',part,'_strata.tsv')),sep='\t');a
}
broad<-aggregate_pb('broad');epi<-aggregate_pb('epithelial')
fwrite(cells[,.(n_cells=.N,n_patients=uniqueN(patient_id),n_libraries=uniqueN(gsm),MFSD10_UMI=sum(MFSD10_counts),n_detected=sum(MFSD10_detected)),by=.(tissue_context,preparation,cell_type,epithelial_status)],file.path(out,'results/Coverage_by_context_celltype_state.tsv'),sep='\t')
fwrite(broad[n_cells>=20,.(n_patient_groups=.N,n_patients=uniqueN(patient_id),n_cells=sum(n_cells),n_detected=sum(n_expressing),cell_weighted_detection=100*sum(n_expressing)/sum(n_cells),patient_median_detection=median(pct_detected),patient_q25_detection=quantile(pct_detected,.25),patient_q75_detection=quantile(pct_detected,.75),patient_median_CPM=median(CPM),patient_median_TMM_CPM=median(TMM_CPM,na.rm=TRUE),patients_with_detection=sum(n_expressing>0),largest_patient_UMI_share=max(MFSD10_UMI)/max(1,sum(MFSD10_UMI)),top3_patient_UMI_share=sum(head(sort(MFSD10_UMI,decreasing=TRUE),3))/max(1,sum(MFSD10_UMI)),LOO_median_CPM_min=if(.N>1)min(vapply(seq_len(.N),function(j)median(CPM[-j]),numeric(1))) else NA_real_,LOO_median_CPM_max=if(.N>1)max(vapply(seq_len(.N),function(j)median(CPM[-j]),numeric(1))) else NA_real_),by=.(tissue_context,preparation,cell_type)],file.path(out,'results/Patient_balanced_celltype_overview.tsv'),sep='\t')
# A leave-one-patient-out statistic is undefined for a one-patient stratum.
ovfile<-file.path(out,'results/Patient_balanced_celltype_overview.tsv');ov<-fread(ovfile)
for(v in c('LOO_median_CPM_min','LOO_median_CPM_max'))ov[!is.finite(get(v)),(v):=NA_real_]
fwrite(ov,ovfile,sep='\t')
sources<-unique(normalizePath(sources,winslash='/',mustWork=TRUE));fi<-file.info(sources)
fwrite(data.table(path=sources,bytes=fi$size,mtime=as.character(fi$mtime)),file.path(out,'audit/Read_only_input_manifest.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'logs/01_sessionInfo.txt'))
cat('EXTRACTION_COMPLETE',format(Sys.time()),'CELLS',nrow(cells),'\n');sink()
