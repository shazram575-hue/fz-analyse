options(stringsAsFactors=FALSE);.libPaths(c('D:/R4.5.2/R-4.5.2/library',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix)})
root<-normalizePath(commandArgs(TRUE)[1],winslash='/');out<-file.path(root,'outputs/MFSD10_程序发现与独立验证_2026-10-26')
z<-readRDS(file.path(root,'outputs/第一阶段_2026-10-23/cache/Main_malignant_normal_raw_counts.rds'))
md<-fread(file.path(out,'metadata/GSE161529_confirmed_malignant_cells.tsv'))[patient_id%in%c('Cancer__0126','Cancer__0043','Cancer__0031')]
x<-z$counts[,md$cell_id,drop=FALSE];stopifnot(all(Matrix::colSums(x)==md$nCount_RNA))
Matrix::writeMM(x,file.path(out,'cache/GSE161529_three_patients_all_confirmed.mtx'));fwrite(md,file.path(out,'cache/GSE161529_three_patients_all_confirmed_metadata.tsv'),sep='\t')
