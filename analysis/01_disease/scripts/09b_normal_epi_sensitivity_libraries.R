options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
m<-fread(file.path(out,"metadata","GSE161529_author_QC_and_analysis_membership.tsv"))
m<-m[analysis_role=="normal_epi_sensitivity"]
stopifnot(nrow(m)==11L)
statuses<-list()
for(i in seq_len(nrow(m))) {
  gsm<-m$gsm[i]
  path<-file.path(out,"cache","seurat_pilot",paste0(gsm,"_pilot_seurat.rds"))
  if(!file.exists(path)) {
    cat("NORMAL EPI SENSITIVITY",i,"/11",gsm,format(Sys.time()),"\n")
    options(workflow3.library.args=c(root,gsm))
    env<-new.env(parent=globalenv())
    tryCatch(source(file.path(out,"scripts","07_seurat_library_pilot.R"),local=env,encoding="UTF-8"),
      error=function(e){while(sink.number())sink();stop(gsm,": ",conditionMessage(e))})
    rm(env);gc(verbose=FALSE)
  }
  statuses[[i]]<-data.table(gsm=gsm,state="complete")
  fwrite(rbindlist(statuses),file.path(out,"audit","Normal_Epi_sensitivity_library_status.tsv"),sep="\t")
}
writeLines("complete",file.path(out,"cache","Normal_Epi_sensitivity_libraries.complete"))
cat("ALL NORMAL EPI LIBRARIES COMPLETE",format(Sys.time()),"\n")
