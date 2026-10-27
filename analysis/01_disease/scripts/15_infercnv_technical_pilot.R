options(stringsAsFactors=FALSE,warn=1,scipen=100)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(infercnv))
set.seed(20261023)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
cd<-file.path(out,"cache","epithelial_by_library")
pd<-file.path(out,"cache","infercnv_technical_pilot")
dir.create(pd,recursive=TRUE,showWarnings=FALSE)
sink(file.path(out,"logs","15_infercnv_technical_pilot.log"),split=TRUE)
cat("START",format(Sys.time()),"\n")
getepi<-function(gsm)readRDS(file.path(cd,paste0(gsm,"_epithelial_counts.rds")))
ref<-getepi("GSM4909276");ctl<-getepi("GSM4909253");qry<-getepi("GSM4909281")
balanced<-function(x,n) {
  groups<-split(x$metadata$cell_id,as.character(x$metadata$RNA_snn_res.0.4))
  each<-ceiling(n/length(groups))
  ids<-unlist(lapply(groups,function(z)sample(z,min(length(z),each))),use.names=FALSE)
  if(length(ids)>n)ids<-sample(ids,n)
  ids
}
ri<-balanced(ref,150);ci<-balanced(ctl,150);qi<-balanced(qry,500)
stopifnot(identical(rownames(ref$counts),rownames(ctl$counts)),identical(rownames(ref$counts),rownames(qry$counts)))
x<-cbind(ref$counts[,ri,drop=FALSE],ctl$counts[,ci,drop=FALSE],qry$counts[,qi,drop=FALSE])
ann<-data.table(cell_id=colnames(x),group=c(rep("healthy_reference",length(ri)),rep("heldout_healthy",length(ci)),rep("tumor_query",length(qi))))
fwrite(ann,file.path(pd,"cell_annotations.tsv"),sep="\t",col.names=FALSE)
fwrite(ann,file.path(out,"audit","CNV_technical_pilot_sampled_cells.tsv"),sep="\t")
rm(ref,ctl,qry);gc(verbose=FALSE)
obj<-CreateInfercnvObject(raw_counts_matrix=x,
 gene_order_file=file.path(out,"metadata","genome_reference","infercnv_GENCODEv36_autosomal_gene_order.tsv"),
 annotations_file=file.path(pd,"cell_annotations.tsv"),ref_group_names="healthy_reference")
rm(x);gc(verbose=FALSE)
cat("RUN",format(Sys.time()),"CELLS",nrow(ann),"\n");flush.console()
t0<-proc.time()
res<-infercnv::run(obj,cutoff=.1,out_dir=pd,HMM=FALSE,denoise=TRUE,
 analysis_mode="samples",cluster_by_groups=TRUE,cluster_references=FALSE,
 num_threads=2,plot_steps=FALSE,inspect_subclusters=FALSE,no_plot=TRUE,no_prelim_plot=TRUE,
 save_rds=FALSE,save_final_rds=TRUE,resume_mode=FALSE,write_expr_matrix=FALSE)
cat("RUNTIME_SECONDS",(proc.time()-t0)["elapsed"],"\n")
e<-res@expr.data
cat("FINAL DIM",dim(e),"RANGE",range(e),"MEDIAN",median(e),"\n")
stopifnot(all(is.finite(e)),all(e>0),median(e)>.5,median(e)<1.5)
score<-sqrt(colMeans((e-1)^2))
ann<-ann[match(colnames(e),cell_id)]
ann[,cnv_rms:=score]
thr<-quantile(ann[group=="heldout_healthy",cnv_rms],.99,names=FALSE)
ann[,above_heldout_healthy_q99:=cnv_rms>thr]
fwrite(ann,file.path(out,"audit","CNV_technical_pilot_scores.tsv"),sep="\t")
fwrite(ann[,.(n_cells=.N,median_cnv_rms=median(cnv_rms),above_q99=mean(above_heldout_healthy_q99)),by=group],
 file.path(out,"audit","CNV_technical_pilot_group_summary.tsv"),sep="\t")
pdf(file.path(out,"figures","CNV_technical_pilot_score_distributions.pdf"),width=8,height=6)
boxplot(cnv_rms~group,data=as.data.frame(ann),ylab="RNA-inferred CNV RMS around reference",las=2,
 col=c("#4E79A7","#59A14F","#E15759"),main="Technical pilot; not final malignant-cell calls")
abline(h=thr,lty=2)
dev.off()
chr<-as.character(res@gene_order[rownames(e),1])
chrlev<-paste0("chr",1:22)
chrmed<-do.call(rbind,lapply(chrlev,function(k)apply(e[chr==k,,drop=FALSE],2,median)))
rownames(chrmed)<-chrlev
saveRDS(list(chromosome_medians=chrmed,scores=ann,heldout_q99=thr,
             note="Balanced technical pilot sample; not Pool A classification"),file.path(pd,"pilot_summary.rds"))
png(file.path(out,"figures","CNV_technical_pilot_chromosome_heatmap.png"),width=1500,height=650,res=150,bg="white")
par(mar=c(5,5,3,5))
limits<-max(abs(chrmed-1));limits<-max(limits,.01)
image(seq_len(ncol(chrmed)),1:22,t(chrmed-1),zlim=c(-limits,limits),
 col=colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101),axes=FALSE,
 xlab="",ylab="",main="Pilot chromosome medians relative to reference; not final calls")
axis(2,at=1:22,labels=1:22,las=2)
axis(1,at=c(length(ri)/2,length(ri)+length(ci)/2,length(ri)+length(ci)+length(qi)/2),
 labels=c("Healthy reference","Independent healthy control","Tumor epithelial query"))
abline(v=c(length(ri)+.5,length(ri)+length(ci)+.5),lty=2)
box();dev.off()
capture.output(sessionInfo(),file=file.path(out,"logs","15_R_sessionInfo.txt"))
cat("COMPLETE",format(Sys.time()),"\n")
sink()
