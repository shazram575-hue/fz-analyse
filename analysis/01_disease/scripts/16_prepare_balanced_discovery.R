# Patient-balanced discovery cohort; full atlas remains preserved separately.
# Fixed cap is independent of candidate genes and differential-expression results.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
setDTthreads(2L)
set.seed(20261023)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
sink(file.path(out,"logs","16_prepare_discovery.log"),split=TRUE)
m<-fread(file.path(out,"metadata","GSE161529_reconciled_patient_metadata.tsv"))
m<-m[include_main==TRUE]
start<-Sys.time()
repeat {
  s<-tryCatch(fread(file.path(out,"audit","Annotation_run_status.tsv")),error=function(e)NULL)
  if(!is.null(s) && nrow(s)==47L && all(s$state=="complete"))break
  if(as.numeric(difftime(Sys.time(),start,units="hours"))>5)stop("Waiting for annotation exceeded five hours")
  Sys.sleep(15)
}
cat("START",format(Sys.time()),"\n")
epidir<-file.path(out,"cache","epithelial_by_library")
library_meta<-rbindlist(lapply(m$gsm,function(gsm){
  x<-readRDS(file.path(epidir,paste0(gsm,"_epithelial_counts.rds")))
  d<-as.data.table(x$metadata)
  d[,gsm:=gsm]
  d
}),fill=TRUE)
library_meta<-merge(library_meta,m[,.(gsm,patient_id,analysis_role,clinical_age,gender)],by="gsm",all.x=TRUE)
stopifnot(!anyDuplicated(library_meta$cell_id))
selected<-rbindlist(lapply(sort(unique(library_meta$patient_id)),function(p){
  md<-library_meta[patient_id==p]
  md[sample.int(nrow(md),min(1000L,nrow(md)))]
}))
fwrite(selected,file.path(out,"metadata","Discovery_patient_balanced_selected_cells.tsv"),sep="\t")
fwrite(library_meta[,.(total_epithelial_cells=.N,selected_cells=sum(cell_id%in%selected$cell_id)),by=.(patient_id,analysis_role)],
 file.path(out,"audit","Discovery_sampling_by_patient.tsv"),sep="\t")
fwrite(data.table(seed=20261023L,max_cells_per_patient=1000L,
                  sampling="Uniform random over epithelial-origin cells within public donor key; not equal-cluster sampling"),
 file.path(out,"audit","Discovery_sampling_parameters.tsv"),sep="\t")
f<-fread(file.path(out,"metadata","GSE161529_features.tsv"))
rs<-numeric(nrow(f));nz<-integer(nrow(f));n<-0L
for(loop_gsm in m$gsm) {
  ids<-selected[gsm==loop_gsm,cell_id]
  if(!length(ids))next
  x<-readRDS(file.path(epidir,paste0(loop_gsm,"_epithelial_counts.rds")))$counts
  z<-x[,ids,drop=FALSE]
  stopifnot(identical(rownames(z),f$ensembl_id))
  rs<-rs+Matrix::rowSums(z);nz<-nz+Matrix::rowSums(z>0);n<-n+ncol(z)
  rm(x,z);gc(verbose=FALSE)
}
go<-fread(file.path(out,"metadata","genome_reference","infercnv_GENCODEv36_autosomal_gene_order.tsv"),
          header=FALSE,col.names=c("gene_id","chr","start","end"))
ga<-copy(f)
ga[,`:=`(mean_raw_count_balanced=rs/n,n_cells_nonzero=nz,autosomal_coordinates=ensembl_id%in%go$gene_id)]
ga[,cnv_feature_panel:=autosomal_coordinates & mean_raw_count_balanced>=.1 & n_cells_nonzero>=3L]
fwrite(ga,file.path(out,"audit","CNV_frozen_expression_feature_panel.tsv"),sep="\t")
genes<-go[gene_id%in%ga[cnv_feature_panel==TRUE,ensembl_id],gene_id]
stopifnot(length(genes)>5000L)
go<-go[gene_id%in%genes]
fwrite(go,file.path(out,"metadata","genome_reference","infercnv_discovery_frozen_gene_order.tsv"),sep="\t",col.names=FALSE)

normal_gsms<-m[analysis_role=="main_healthy_total",gsm]
ref_gsms<-c("GSM4909254","GSM4909253","GSM4909276")
ctl_gsms<-setdiff(normal_gsms,ref_gsms)
fixed_ids<-list();fixed_ann<-list();fixed_counts<-list()
for(loop_gsm in c(ref_gsms,ctl_gsms)) {
  md<-selected[gsm==loop_gsm]
  nmax<-if(loop_gsm%in%ref_gsms)100L else 20L
  ids<-sample(md$cell_id,min(nrow(md),nmax))
  if(!length(ids))next
  x<-readRDS(file.path(epidir,paste0(loop_gsm,"_epithelial_counts.rds")))$counts
  fixed_counts[[loop_gsm]]<-x[genes,ids,drop=FALSE]
  fixed_ann[[loop_gsm]]<-data.table(cell_id=ids,gsm=loop_gsm,
              group=paste0(if(loop_gsm%in%ref_gsms)"ref__" else "healthy_control__",loop_gsm),
              role=if(loop_gsm%in%ref_gsms)"reference" else "healthy_control")
}
fixed_ann<-rbindlist(fixed_ann)
fixed_counts<-do.call(cbind,fixed_counts)
stopifnot(identical(colnames(fixed_counts),fixed_ann$cell_id),sum(fixed_ann$role=="healthy_control")>=150L)
fwrite(fixed_ann,file.path(out,"metadata","CNV_fixed_reference_and_healthy_controls.tsv"),sep="\t")
saveRDS(list(genes=genes,counts=fixed_counts,annotation=fixed_ann,
             reference_groups=unique(fixed_ann[role=="reference",group]),
             gene_feature_cutoff=.1,seed=20261023L),file.path(out,"cache","CNV_fixed_reference_design.rds"))
writeLines("complete",file.path(out,"cache","CNV_fixed_reference_design.complete"))
cat("COMPLETE",format(Sys.time()),"SELECTED CELLS",nrow(selected),"PATIENTS",uniqueN(selected$patient_id),
    "CNV PANEL",length(genes),"REFERENCE",sum(fixed_ann$role=="reference"),"HEALTHY CONTROLS",sum(fixed_ann$role=="healthy_control"),"\n")
sink()
