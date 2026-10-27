# Descriptive lineage expression for all 813 frozen candidates; no new gate or DE test.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
setDTthreads(2L)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
sink(file.path(out,"logs","23_candidate_celltype_expression.log"),split=TRUE)
m<-fread(file.path(out,"metadata","GSE161529_reconciled_patient_metadata.tsv"))[include_main==TRUE]
pool<-fread(file.path(out,"results","Pool_A_candidate_genes.tsv"))
map<-fread(file.path(out,"metadata","SingleCell_canonical_gene_identity_map.tsv"))
cnv<-fread(file.path(out,"results","Tumor_discovery_cells_CNV_coherent_final_calls.tsv"))
malig<-cnv[cnv_high_q99==TRUE,cell_id]
fp<-file.path(out,"results","Pool_A_expression_by_donor_celltype.tsv")
parts<-list()
for(g in m$gsm) {
 cat("LINEAGE EXPRESSION",g,format(Sys.time()),"\n");flush.console()
 so<-readRDS(file.path(out,"cache","seurat_pilot",paste0(g,"_pilot_seurat.rds")))
 raw<-GetAssayData(so,assay="RNA",layer="counts")
 stopifnot(identical(rownames(raw),map$ensembl_id))
 ix<-which(map$canonical_id%in%pool$canonical_id)
 collapse<-sparseMatrix(i=match(map$canonical_id[ix],pool$canonical_id),j=seq_along(ix),x=1,
 dims=c(nrow(pool),length(ix)))
 x<-collapse%*%raw[ix,,drop=FALSE]
 rownames(x)<-pool$canonical_id
 md<-fread(file.path(out,"audit","cell_annotation",paste0(g,"_cell_identity.tsv")))
 md<-md[match(colnames(raw),cell_id)]
 stopifnot(identical(md$cell_id,colnames(raw)),all(md$nCount_RNA>0))
 md[,description:=cell_type]
 md[cell_type=="Epithelial",description:=if(m[gsm==g,analysis_role]=="main_healthy_total")"Healthy_Epithelial" else "Tumor_epithelial_uncertain_or_unscored"]
 md[cell_id%in%malig,description:="Malignant_CNV_high"]
 linear<-x%*%Diagonal(x=10000/md$nCount_RNA)
 logn<-linear;logn@x<-log1p(logn@x)
 for(ct in unique(md$description)) {
  j<-which(md$description==ct)
  parts[[length(parts)+1L]]<-data.table(canonical_id=pool$canonical_id,gene_symbol=pool$gene_symbol,
   gsm=g,patient_id=m[gsm==g,patient_id],analysis_role=m[gsm==g,analysis_role],cell_type=ct,
   n_cells=length(j),n_expressing=as.numeric(Matrix::rowSums(x[,j,drop=FALSE]>0)),
   sum_log1p_normalized=as.numeric(Matrix::rowSums(logn[,j,drop=FALSE])),
   sum_UMI=as.numeric(Matrix::rowSums(x[,j,drop=FALSE])))
 }
 rm(so,raw,x,linear,logn,collapse);gc(verbose=FALSE)
}
all<-rbindlist(parts);rm(parts);gc(verbose=FALSE)
don<-all[,.(n_cells=sum(n_cells),n_expressing=sum(n_expressing),
 sum_log1p_normalized=sum(sum_log1p_normalized),sum_UMI=sum(sum_UMI)),
 by=.(canonical_id,gene_symbol,patient_id,analysis_role,cell_type)]
don[,`:=`(pct_expressing=n_expressing/n_cells,mean_log1p_normalized=sum_log1p_normalized/n_cells,
 descriptive_min20_donor=n_cells>=20)]
fwrite(don,fp,sep="\t")
sm<-don[descriptive_min20_donor==TRUE,.(n_donors=.N,n_cells=sum(n_cells),
 cell_weighted_pct=sum(n_expressing)/sum(n_cells),donor_mean_pct=mean(pct_expressing),
 donor_mean_log1p_normalized=mean(mean_log1p_normalized)),by=.(canonical_id,gene_symbol,analysis_role,cell_type)]
fwrite(sm,file.path(out,"results","Pool_A_expression_by_celltype_summary.tsv"),sep="\t")
stopifnot(setequal(don$canonical_id,pool$canonical_id))
totals<-unique(don[,.(patient_id,analysis_role,cell_type,n_cells)])
stopifnot(sum(totals$n_cells)==236309L,totals[cell_type=="Malignant_CNV_high",sum(n_cells)]==9486L)
fwrite(totals,file.path(out,"audit","Candidate_celltype_expression_cell_totals.tsv"),sep="\t")
writeLines(c("Descriptive only: no candidate exclusions or new P values.",
 "All main QC cells included; per-donor summaries combine repeat libraries.",
 "Summary uses donor/cell-type groups with >=20 cells, with cell-weighted and donor-mean detection separately reported.",
 "Healthy epithelium here uses the full atlas, not the 12,264-cell discovery cap.",
 "Tumor epithelial uncertain or unscored cells are never called normal or malignant by projection."),
 file.path(out,"audit","Candidate_celltype_expression_method.txt"))
capture.output(sessionInfo(),file=file.path(out,"logs","23_R_sessionInfo.txt"))
cat("COMPLETE",format(Sys.time()),"DONOR_GENE_ROWS",nrow(don),"SUMMARY_ROWS",nrow(sm),"\n")
sink()
