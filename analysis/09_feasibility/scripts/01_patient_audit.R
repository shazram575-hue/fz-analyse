options(stringsAsFactors=FALSE,warn=1)
.libPaths(c('D:/R4.5.2/R-4.5.2/library',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix)})
setDTthreads(2L)
root<-normalizePath(commandArgs(TRUE)[1],winslash='/');old<-file.path(root,'outputs/第一阶段_2026-10-23');survey<-file.path(root,'outputs/MFSD10_表达分布摸底_2026-10-26');out<-file.path(root,'outputs/恶性上皮_表达程序可行性_2026-10-26')
for(d in c('results','audit','logs','sources'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
sink(file.path(out,'logs/01_patient_audit.log'),split=TRUE)
cnv<-fread(file.path(old,'results/Tumor_discovery_cells_CNV_coherent_final_calls.tsv'))
md<-cnv[cnv_high_q99==TRUE];stopifnot(nrow(md)==9486,!anyDuplicated(md$cell_id),all(md$cell_type=='Epithelial'))
m<-fread(file.path(survey,'metadata/Library_clinical_context.tsv'))
md<-merge(md,m[,.(gsm,clinical_specimen,clinical_label,preparation,tissue_context)],by='gsm',all.x=TRUE,sort=FALSE)
obj<-readRDS(file.path(old,'cache/Main_malignant_normal_raw_counts.rds'));x<-obj$counts
stopifnot(!anyDuplicated(colnames(x)),!anyDuplicated(rownames(x)),all(md$cell_id%in%colnames(x)))
x<-x[,md$cell_id,drop=FALSE];stopifnot(identical(rownames(x),obj$genes$canonical_id))
tot<-Matrix::colSums(x);features<-Matrix::colSums(x>0);all_integer<-all(x@x>=0 & is.finite(x@x) & x@x==round(x@x))
stopifnot(all_integer,all(tot==md$nCount_RNA),all(tot>0))
md[,`:=`(raw_counts_UMI=as.numeric(tot),canonical_nFeature=as.numeric(features),is_singlet=scDblFinder_class=='singlet')]
qc<-as.data.table(read.delim(gzfile(file.path(old,'cache/GSE161529_all_cell_QC_membership.tsv.gz')),check.names=FALSE))
qq<-qc[match(md$cell_id,cell_id)];stopifnot(identical(qq$cell_id,md$cell_id),all(qq$pass_author_qc))
fwrite(md[,.(cell_id,gsm,patient_id,clinical_specimen,clinical_label,nCount_RNA,nFeature_RNA,canonical_nFeature,percent.mt,scDblFinder_class,is_singlet)],file.path(out,'audit/Confirmed_malignant_cell_inventory.tsv'),sep='\t')
fnum<-function(z)as.numeric(median(z))
parts<-md[,.(confirmed_malignant_cells=.N,singlet_flag_cells=sum(is_singlet),doublet_flag_cells=sum(scDblFinder_class=='doublet'),UMI_median=fnum(nCount_RNA),nFeature_RNA_median=fnum(nFeature_RNA),canonical_nFeature_median=fnum(canonical_nFeature),mitochondrial_pct_median=fnum(percent.mt)),by=.(patient_id,gsm,clinical_specimen,clinical_label)]
fwrite(parts,file.path(out,'results/Confirmed_malignant_sample_composition.tsv'),sep='\t')
pat<-md[,.(n_confirmed=.N,n_singlet=sum(is_singlet),n_doublet=sum(scDblFinder_class=='doublet'),doublet_pct=100*mean(scDblFinder_class=='doublet'),UMI_median=fnum(nCount_RNA),UMI_q25=as.numeric(quantile(nCount_RNA,.25)),UMI_q75=as.numeric(quantile(nCount_RNA,.75)),genes_median=fnum(nFeature_RNA),canonical_genes_median=fnum(canonical_nFeature),genes_q25=as.numeric(quantile(nFeature_RNA,.25)),genes_q75=as.numeric(quantile(nFeature_RNA,.75)),mt_pct_median=fnum(percent.mt),singlet_UMI_median=if(any(is_singlet))fnum(nCount_RNA[is_singlet]) else NA_real_,singlet_genes_median=if(any(is_singlet))fnum(nFeature_RNA[is_singlet]) else NA_real_,n_libraries_with_malignant=uniqueN(gsm),malignant_clinical_labels=paste(sort(unique(clinical_label)),collapse=';'),malignant_sample_composition=paste(paste0(gsm,'[',clinical_specimen,']:',.N),collapse=';')),by=patient_id]
# Explicit per-library counts rather than repeating a patient's total for every library.
comp<-parts[,.(malignant_sample_composition=paste(paste0(gsm,'[',clinical_specimen,']:',confirmed_malignant_cells,'(S=',singlet_flag_cells,',D=',doublet_flag_cells,')'),collapse=';')),by=patient_id]
pat[,malignant_sample_composition:=comp$malignant_sample_composition[match(patient_id,comp$patient_id)]]
allprimary<-m[tissue_context=='Primary_tumor',.(all_primary_clinical_labels=paste(sort(unique(clinical_label)),collapse=';'),all_primary_GSMs=paste(gsm,collapse=';'),all_primary_specimens=paste(clinical_specimen,collapse=';'),n_primary_libraries=.N),by=patient_id]
pat<-merge(pat,allprimary,by='patient_id',all.x=TRUE,sort=FALSE)
# Complexity of the existing singlet-flag subset: audit only; no cells are removed or saved as a new analysis object.
pat[,`:=`(singlet_genes_detected_ge1=0L,singlet_genes_detected_ge5=0L,singlet_genes_detected_ge3pct=0L,MFSD10_present_in_raw_counts='ENSG00000109736'%in%rownames(x))]
for(p in pat$patient_id){
 j<-which(md$patient_id==p & md$is_singlet);if(!length(j))next
 ng<-Matrix::rowSums(x[,j,drop=FALSE]>0)
 # MFSD10 is excluded from these prospective fitting-universe diagnostics.
 ng<-ng[rownames(x)!='ENSG00000109736']
 pat[patient_id==p,`:=`(singlet_genes_detected_ge1=sum(ng>=1),singlet_genes_detected_ge5=sum(ng>=5),singlet_genes_detected_ge3pct=sum(ng>=max(5L,ceiling(length(j)*.03))))]
}
setorder(pat,-n_singlet,patient_id)
fwrite(pat,file.path(out,'results/Patient_counts_complexity_audit.tsv'),sep='\t')
checks<-data.table(check=c('Confirmed_9486','All_32_patient_keys','All_cell_IDs_unique','Every_confirmed_cell_in_raw_counts','Canonical_features_33487','Nonnegative_finite_integer_counts','UMI_equal_metadata','All_cells_pass_previous_author_QC','No_CNV_uncertain_or_unscored','Confirmed_all_primary_Total','Singlet_plus_doublet_equals_total'),pass=c(nrow(md)==9486,uniqueN(md$patient_id)==32,!anyDuplicated(md$cell_id),all(md$cell_id%in%colnames(x)),nrow(x)==33487,all_integer,all(tot==md$nCount_RNA),all(qq$pass_author_qc),all(md$cnv_high_q99==TRUE),all(md$tissue_context=='Primary_tumor' & md$preparation=='Total'),all(pat$n_singlet+pat$n_doublet==pat$n_confirmed)))
fwrite(checks,file.path(out,'audit/Counts_QC_validation.tsv'),sep='\t');stopifnot(all(checks$pass))
capture.output(sessionInfo(),file=file.path(out,'logs/R_sessionInfo.txt'))
print(pat[,.(patient_id,n_confirmed,n_singlet,n_doublet,UMI_median,genes_median,singlet_UMI_median,singlet_genes_median,singlet_genes_detected_ge3pct,malignant_clinical_labels,malignant_sample_composition)])
cat('QC_AUDIT_COMPLETE; NMF_NOT_RUN\n');sink()
