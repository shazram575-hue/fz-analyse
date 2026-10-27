options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table);library(Matrix)
setDTthreads(2L);set.seed(20261024)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
start<-Sys.time()
while(!file.exists(file.path(out,"cache","CNV_fixed_reference_design.complete"))) {
  if(as.numeric(difftime(Sys.time(),start,units="hours"))>5)stop("Discovery selection did not complete")
  Sys.sleep(15)
}
sel<-fread(file.path(out,"metadata","Discovery_patient_balanced_selected_cells.tsv"))
m<-fread(file.path(out,"metadata","GSE161529_reconciled_patient_metadata.tsv"))
norm<-m[include_main==TRUE & analysis_role=="main_healthy_total"]
ref_gsms<-c("GSM4909254","GSM4909253","GSM4909276")
remaining<-norm[!gsm%in%ref_gsms][order(clinical_age,gsm),gsm]
calib_gsms<-remaining[seq(1,length(remaining),by=2)]
validation_gsms<-setdiff(remaining,calib_gsms)
go<-fread(file.path(out,"metadata","genome_reference","infercnv_GENCODEv36_autosomal_gene_order.tsv"),
  header=FALSE,col.names=c("gene_id","chr","start","end"))
parts<-list();annotations<-list()
for(g in c(ref_gsms,calib_gsms,validation_gsms)) {
  md<-sel[gsm==g]
  role<-if(g%in%ref_gsms)"reference" else if(g%in%calib_gsms)"healthy_calibration" else "healthy_validation"
  cap<-if(role=="reference")100L else 40L
  ids<-sample(md$cell_id,min(nrow(md),cap))
  x<-readRDS(file.path(out,"cache","epithelial_by_library",paste0(g,"_epithelial_counts.rds")))$counts
  parts[[g]]<-x[go$gene_id,ids,drop=FALSE]
  annotations[[g]]<-data.table(cell_id=ids,gsm=g,source_patient_id=md$patient_id[match(ids,md$cell_id)],
        group=paste0(if(role=="reference")"ref__" else paste0(role,"__"),g),role=role)
}
x<-do.call(cbind,parts);ann<-rbindlist(annotations)
stopifnot(identical(colnames(x),ann$cell_id),length(intersect(calib_gsms,validation_gsms))==0L)
train_ix<-which(ann$role%in%c("reference","healthy_calibration"))
mean_count<-Matrix::rowMeans(x[,train_ix,drop=FALSE])
nz<-Matrix::rowSums(x[,train_ix,drop=FALSE]>0)
keep<-mean_count>=.1 & nz>=3L
stopifnot(sum(keep)>5000L)
genes<-go$gene_id[keep]
fwrite(data.table(gene_id=go$gene_id,reference_calibration_mean_count=mean_count,
       reference_calibration_nz_cells=nz,cnv_feature_keep=keep),
  file.path(out,"audit","CNV_final_reference_expression_panel.tsv"),sep="\t")
fwrite(go[keep],file.path(out,"metadata","genome_reference","infercnv_discovery_frozen_gene_order_v2.tsv"),sep="\t",col.names=FALSE)
fwrite(ann,file.path(out,"metadata","CNV_final_reference_calibration_validation_cells.tsv"),sep="\t")
saveRDS(list(genes=genes,counts=x[genes,,drop=FALSE],annotation=ann,
   reference_groups=unique(ann[role=="reference",group]),
   feature_selection="Autosomal genes with mean raw count>=0.1 and detection>=3 in reference+calibration healthy cells; validation donors excluded from panel selection",
   seed=20261024L),file.path(out,"cache","CNV_fixed_reference_design_v2.rds"))
writeLines("complete",file.path(out,"cache","CNV_fixed_reference_v2.complete"))
print(ann[,.(cells=.N,donors=uniqueN(source_patient_id)),by=role])
cat("CNV FROZEN GENES",length(genes),format(Sys.time()),"\n")
