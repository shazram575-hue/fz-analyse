options(stringsAsFactors=FALSE,warn=1)
.libPaths(c('D:/R4.5.2/R-4.5.2/library',.libPaths()));suppressPackageStartupMessages({library(data.table);library(Matrix)})
setDTthreads(2L)
root<-normalizePath(commandArgs(TRUE)[1],winslash='/');out<-file.path(root,'outputs/MFSD10_程序发现与独立验证_2026-10-26');old<-file.path(root,'outputs/第一阶段_2026-10-23')
sink(file.path(out,'logs/01_discovery_input.log'),split=TRUE)
z<-readRDS(file.path(old,'cache/Main_malignant_normal_raw_counts.rds'))
cn<-fread(file.path(old,'results/Tumor_discovery_cells_CNV_coherent_final_calls.tsv'))
full<-fread(file.path(root,'outputs/单细胞基础补全_2026-10-23/metadata/All_main_cell_metadata.tsv'))
m<-fread(file.path(old,'metadata/GSE161529_reconciled_patient_metadata.tsv'))
mal<-cn[cnv_high_q99==TRUE];stopifnot(nrow(mal)==9486,!anyDuplicated(mal$cell_id))
mal<-merge(mal,m[,.(gsm,clinical_specimen,cancer.type,treatment)],by='gsm',all.x=TRUE,sort=FALSE)
mal[,clinical_group:=fcase(cancer.type%in%c('Triple negative tumour','Triple negative BRCA1 tumour'),'TNBC',cancer.type=='HER2+ tumour','HER2+',cancer.type=='PR+ tumour','PR+',default='ER+')]
mal[,treatment_status:=ifelse(treatment=='none','none_reported','Unknown')]
fwrite(mal,file.path(out,'metadata/GSE161529_confirmed_malignant_cells.tsv'),sep='\t')
md<-mal[scDblFinder_class=='singlet'];stopifnot(nrow(md)==8722,all(md$cell_id%in%colnames(z$counts)))
x<-z$counts[,md$cell_id,drop=FALSE];stopifnot(all(Matrix::colSums(x)==md$nCount_RNA),all(x@x>=0 & x@x==round(x@x)))
Matrix::writeMM(x,file.path(out,'cache/GSE161529_singlet_counts_gene_by_cell.mtx'))
fwrite(md,file.path(out,'metadata/GSE161529_singlet_metadata.tsv'),sep='\t');fwrite(z$genes,file.path(out,'metadata/GSE161529_canonical_genes.tsv'),sep='\t')
inv<-full[cell_type=='Epithelial',.(full_atlas_epithelial_cells=.N),by=patient_id]
disc<-cn[,.(CNV_scored=.N,confirmed_malignant=sum(cnv_high_q99),uncertain=sum(!cnv_high_q99)),by=patient_id]
inv<-merge(inv,disc,by='patient_id',all.x=TRUE);inv[,unscored:=full_atlas_epithelial_cells-CNV_scored]
fwrite(inv,file.path(out,'audit/Full_unsampled_atlas_CNV_coverage.tsv'),sep='\t')
fwrite(md[,.(n_singlets=.N,median_UMI=as.numeric(median(nCount_RNA)),median_genes=as.numeric(median(nFeature_RNA)),clinical_labels=paste(sort(unique(cancer.type)),collapse=';'),GSMs=paste(unique(gsm),collapse=';')),by=patient_id],file.path(out,'metadata/GSE161529_patient_input_summary.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'logs/R_export_versions.txt'));cat('EXPORTED',dim(x),'\n');sink()
