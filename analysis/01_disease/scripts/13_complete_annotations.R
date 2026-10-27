options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
setDTthreads(2L)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
source(file.path(out,"scripts","annotate_library.R"),encoding="UTF-8")
f<-fread(file.path(out,"metadata","GSE161529_features.tsv"))
m<-fread(file.path(out,"metadata","GSE161529_reconciled_patient_metadata.tsv"))
m<-m[include_main==TRUE]
ad<-file.path(out,"audit","cell_annotation")
completed<-character()
start<-Sys.time()
repeat {
  state<-tryCatch(fread(file.path(out,"audit","Seurat_main_library_run_status.tsv")),error=function(e)NULL)
  ready<-if(is.null(state))character() else state[state=="complete",gsm]
  for(gsm in setdiff(intersect(m$gsm,ready),completed)) {
    clusterfile<-file.path(ad,paste0(gsm,"_cluster_identity.tsv"))
    current<-FALSE
    if(file.exists(clusterfile)) {
      ct<-fread(clusterfile)
      current<-"annotation_version"%in%names(ct) && all(ct$annotation_version=="HPCA_marker_v1.1") &&
         file.exists(file.path(out,"cache","epithelial_by_library",paste0(gsm,"_epithelial_counts.rds")))
    }
    if(!current) {
      cat("ANNOTATE",gsm,format(Sys.time()),"\n");flush.console()
      so<-readRDS(file.path(out,"cache","seurat_pilot",paste0(gsm,"_pilot_seurat.rds")))
      annotate_library(so,f,root,gsm)
      rm(so);gc(verbose=FALSE)
    }
    completed<-union(completed,gsm)
    fwrite(data.table(gsm=completed,state="complete",version="HPCA_marker_v1.1"),
      file.path(out,"audit","Annotation_run_status.tsv"),sep="\t")
  }
  if(length(completed)==nrow(m))break
  if(as.numeric(difftime(Sys.time(),start,units="hours"))>5)stop("Annotation wait exceeded five hours; inspect upstream status")
  Sys.sleep(10)
}
sm<-rbindlist(lapply(m$gsm,function(gsm)fread(file.path(ad,paste0(gsm,"_celltype_counts.tsv")))))
sm<-merge(sm,m[,.(gsm,patient_id,analysis_role)],by="gsm",all.x=TRUE)
fwrite(sm,file.path(out,"results","SingleCell_broad_celltypes_by_library.tsv"),sep="\t")
cell_tables<-lapply(m$gsm,function(gsm) {
  fread(file.path(ad,paste0(gsm,"_cell_identity.tsv")))
})
allmd<-rbindlist(cell_tables,fill=TRUE)
allmd<-merge(allmd,m[,.(gsm,patient_id,analysis_role,clinical_age,gender)],by="gsm",all.x=TRUE)
stopifnot(!anyDuplicated(allmd$cell_id),nrow(allmd)==236309L)
con<-gzfile(file.path(out,"cache","All_main_cells_broad_annotation.tsv.gz"),"wt")
write.table(allmd,con,sep="\t",row.names=FALSE,quote=FALSE);close(con)
fwrite(allmd[,.(n_cells=.N,n_patients=uniqueN(patient_id)),by=.(analysis_role,cell_type)],
 file.path(out,"results","SingleCell_broad_annotation_summary.tsv"),sep="\t")
cat("ALL 47 ANNOTATIONS COMPLETE",format(Sys.time()),"CELLS",nrow(allmd),"\n")
