options(stringsAsFactors=FALSE,warn=1,scipen=100)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table);library(Matrix);suppressPackageStartupMessages(library(infercnv))
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
old<-readRDS(file.path(out,"cache","infercnv_technical_pilot","run.final.infercnv_obj"))
ann<-fread(file.path(out,"audit","CNV_technical_pilot_sampled_cells.tsv"))
for(g in setdiff(unique(ann$group),"healthy_reference")) {
  ids<-which(ann$group==g)
  ann$group[ids]<-paste0(g,"__technical_",ceiling(seq_along(ids)/50))
}
pd<-file.path(out,"cache","infercnv_technical_group_validation")
dir.create(pd,recursive=TRUE,showWarnings=FALSE)
fwrite(ann,file.path(pd,"annotations.tsv"),sep="\t",col.names=FALSE)
raw<-old@count.data
stopifnot(all(is.finite(raw)),all(raw>=0),all(raw==floor(raw)))
futile.logger::flog.appender(futile.logger::appender.file(file.path(pd,"infercnv_run.log")))
obj<-CreateInfercnvObject(raw_counts_matrix=raw,
 gene_order_file=file.path(out,"metadata","genome_reference","infercnv_GENCODEv36_autosomal_gene_order.tsv"),
 annotations_file=file.path(pd,"annotations.tsv"),ref_group_names="healthy_reference")
t0<-proc.time()
new<-infercnv::run(obj,cutoff=.1,out_dir=pd,HMM=FALSE,denoise=TRUE,
 analysis_mode="samples",cluster_by_groups=TRUE,cluster_references=FALSE,num_threads=2,
 plot_steps=FALSE,inspect_subclusters=FALSE,no_plot=TRUE,no_prelim_plot=TRUE,
 save_rds=FALSE,save_final_rds=FALSE,resume_mode=FALSE,write_expr_matrix=FALSE)
stopifnot(setequal(rownames(old@expr.data),rownames(new@expr.data)),setequal(colnames(old@expr.data),colnames(new@expr.data)))
a<-old@expr.data;b<-new@expr.data[rownames(a),colnames(a),drop=FALSE]
err<-max(abs(a-b));stopifnot(err<1e-10)
fwrite(data.table(test="Technical observation groups with HMM disabled do not change inferred expression",
 max_absolute_expression_difference=err,validation_runtime_seconds=unname((proc.time()-t0)["elapsed"])),
 file.path(out,"audit","CNV_technical_group_equivalence_check.tsv"),sep="\t")
cat("PASS TECHNICAL GROUPS",err,format(Sys.time()),"\n")
