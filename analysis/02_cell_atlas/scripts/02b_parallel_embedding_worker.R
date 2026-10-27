# UMAP/tSNE do not depend on clustering. Compute in parallel with Louvain using
# the identical frozen inputs/parameters; publish only fully serialized files.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table);setDTthreads(2L)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/")
out<-file.path(root,"outputs","单细胞基础补全_2026-10-23")
sink(file.path(out,"logs","02b_parallel_embedding_worker.log"),split=TRUE)
md<-fread(file.path(out,"metadata","All_main_cell_metadata.tsv"))
hm<-readRDS(file.path(out,"cache","All_main_Harmony_embeddings.rds"))
stopifnot(identical(rownames(hm),md$cell_id))
publish<-function(x,path) {
 temp<-paste0(path,".worker_ready")
 saveRDS(x,temp,compress=FALSE)
 if(!file.exists(path)) {
  ok<-file.rename(temp,path)
  cat("ATOMIC PUBLISH",basename(path),ok,format(Sys.time()),"\n")
 } else cat("PRIMARY FILE ALREADY EXISTS",basename(path),"\n")
}
for(reduction in c("harmony","pca")) {
 path<-file.path(out,"cache",paste0("All_main_UMAP_",reduction,".rds"))
 if(file.exists(path))next
 cat("PARALLEL UMAP",reduction,format(Sys.time()),"\n");flush.console()
 xx<-if(reduction=="harmony")hm else readRDS(file.path(out,"cache","All_main_PCA_embeddings.rds"))[,1:30,drop=FALSE]
 set.seed(20261023)
 um<-uwot::umap(xx,n_neighbors=30,n_components=2,metric="cosine",n_epochs=200,
   min_dist=.3,spread=1,scale=FALSE,init="spectral",n_threads=2,n_sgd_threads=1,ret_model=FALSE,verbose=TRUE)
 rownames(um)<-md$cell_id;publish(um,path)
 rm(um,xx);gc(verbose=FALSE)
}
path<-file.path(out,"cache","Balanced_landmark_tSNE_Harmony.rds")
if(!file.exists(path)) {
 ids<-readRDS(file.path(out,"cache","Balanced_landmark_PCA_model.rds"))$landmark_cell_ids
 cat("PARALLEL LANDMARK tSNE",format(Sys.time()),"\n");flush.console();set.seed(20261023)
 ti<-Rtsne::Rtsne(hm[match(ids,md$cell_id),,drop=FALSE],dims=2,perplexity=30,theta=.5,pca=FALSE,
 check_duplicates=FALSE,max_iter=1000,num_threads=2,verbose=TRUE)
 publish(data.table(cell_id=ids,tSNE_1=ti$Y[,1],tSNE_2=ti$Y[,2]),path)
}
cat("WORKER COMPLETE",format(Sys.time()),"\n");sink()
