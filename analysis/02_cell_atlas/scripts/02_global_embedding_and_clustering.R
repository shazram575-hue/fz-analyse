# Full main-cohort embedding and clustering. PCA is trained on pre-DE donor/type
# balanced landmarks and every main cell is explicitly projected. RNA DE is separate.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(harmony))
setDTthreads(2L);future::plan("sequential");set.seed(20261023)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","单细胞基础补全_2026-10-23")
stopifnot(file.exists(file.path(out,"cache","Preparation.complete")))
sink(file.path(out,"logs","02_global_embedding.log"),split=TRUE)
md<-fread(file.path(out,"metadata","All_main_cell_metadata.tsv"))
hv<-fread(file.path(out,"metadata","Consensus_HVG_ranking.tsv"))[selected_for_PCA==TRUE]
features<-hv$canonical_id
model_file<-file.path(out,"cache","Balanced_landmark_PCA_model.rds")
if(!file.exists(model_file)) {
 cat("BUILD LANDMARK PCA",format(Sys.time()),"\n");flush.console()
 parts<-list()
 for(g in unique(md$gsm)) {
  z<-readRDS(file.path(out,"cache","balanced_by_library",paste0(g,".rds")))
  x<-z$counts[features,,drop=FALSE]%*%Diagonal(x=10000/z$metadata$nCount_RNA)
  x@x<-log1p(x@x);colnames(x)<-z$metadata$cell_id
  parts[[g]]<-x;rm(z,x)
 }
 x<-do.call(cbind,parts);rm(parts);gc(verbose=FALSE)
 av<-Matrix::rowMeans(x);sq<-x;sq@x<-sq@x^2
 sdv<-sqrt(pmax((Matrix::rowSums(sq)-ncol(x)*av^2)/(ncol(x)-1L),0));rm(sq)
 keep<-sdv>1e-8;x<-x[keep,,drop=FALSE];features<-features[keep];av<-av[keep];sdv<-sdv[keep]
 set.seed(20261023)
 # Explicit numerical scaling: the local prcomp_irlba wrapper retains a logical
 # scale flag rather than the fitted scale vector when scale.=TRUE.
 sv<-irlba::irlba(A=t(x),nv=50,nu=50,center=av,scale=sdv,tol=1e-7,maxit=1000)
 pc<-list(rotation=sv$v,sdev=sv$d/sqrt(ncol(x)-1L),center=av,scale=sdv,
          x=sweep(sv$u,2,sv$d,"*"))
 rownames(pc$rotation)<-features;colnames(pc$rotation)<-paste0("PC",1:50)
 rownames(pc$x)<-colnames(x);colnames(pc$x)<-paste0("PC",1:50);class(pc)<-"prcomp"
 # Ensure sparse projection reproduces the fitted scores before projecting others.
 projected<-as.matrix(t(x)%*%(pc$rotation/pc$scale))
 projected<-sweep(projected,2,colSums(pc$rotation*(pc$center/pc$scale)),"-")
 err<-max(abs(projected-pc$x));relative_err<-sqrt(sum((projected-pc$x)^2)/sum(pc$x^2))
 cat("PCA PROJECTION ERROR",err,"RELATIVE",relative_err,"\n");flush.console()
 stopifnot(err<1e-4,relative_err<1e-6)
 saveRDS(list(model=pc,features=features,landmark_cell_ids=colnames(x),projection_max_error=err),model_file,compress=FALSE)
 fwrite(data.table(check=c("PCA_landmark_projection_max_error","PCA_landmark_projection_relative_error"),value=c(err,relative_err)),file.path(out,"audit","PCA_projection_validation.tsv"),sep="\t")
 rm(pc,x,projected);gc(verbose=FALSE)
}
pcmod<-readRDS(model_file);pc<-pcmod$model;features<-pcmod$features
fwrite(data.table(PC=seq_along(pc$sdev),sdev=pc$sdev,variance_fraction=pc$sdev^2/length(features)),
 file.path(out,"results","PCA_variance_explained.tsv"),sep="\t")
loads<-rbindlist(lapply(seq_len(ncol(pc$rotation)),function(i) {
 z<-data.table(canonical_id=rownames(pc$rotation),PC=i,loading=pc$rotation[,i])
 z[order(-abs(loading))][seq_len(20)]
}))
loads<-merge(loads,hv[,.(canonical_id,gene_symbol)],by="canonical_id",all.x=TRUE,sort=FALSE)
fwrite(loads,file.path(out,"results","PCA_top_absolute_loadings.tsv"),sep="\t")
project_file<-file.path(out,"cache","All_main_PCA_embeddings.rds")
if(!file.exists(project_file)) {
 parts<-list();L<-pc$rotation/pc$scale;offset<-colSums(pc$rotation*(pc$center/pc$scale))
 for(g in unique(md$gsm)) {
  cat("PROJECT ALL CELLS",g,format(Sys.time()),"\n");flush.console()
  x<-readRDS(file.path(out,"cache","canonical_by_library",paste0(g,".rds")))
  ids<-colnames(x);mm<-md[match(ids,cell_id)]
  x<-x[features,,drop=FALSE]%*%Diagonal(x=10000/mm$nCount_RNA);x@x<-log1p(x@x)
  scores<-as.matrix(t(x)%*%L);scores<-sweep(scores,2,offset,"-");rownames(scores)<-ids
  parts[[g]]<-scores;rm(x,scores);gc(verbose=FALSE)
 }
 pca<-do.call(rbind,parts);rm(parts)
 pca<-pca[match(md$cell_id,rownames(pca)),,drop=FALSE]
 stopifnot(all(is.finite(pca)),identical(rownames(pca),md$cell_id))
 saveRDS(pca,project_file,compress=FALSE)
} else pca<-readRDS(project_file)
rm(pc);gc(verbose=FALSE)
harm_file<-file.path(out,"cache","All_main_Harmony_embeddings.rds")
if(!file.exists(harm_file)) {
 cat("HARMONY FULL COHORT",format(Sys.time()),"\n");flush.console();set.seed(20261023)
 harm<-harmony::RunHarmony(data_mat=pca[,1:30,drop=FALSE],meta_data=as.data.frame(md),vars_use="gsm",
 theta=1,lambda=1,sigma=.1,nclust=100,max_iter=10,early_stop=TRUE,ncores=2,verbose=TRUE)
 if(nrow(harm)!=nrow(md))harm<-t(harm)
 rownames(harm)<-md$cell_id;stopifnot(ncol(harm)==30L,all(is.finite(harm)))
 saveRDS(harm,harm_file,compress=FALSE)
} else harm<-readRDS(harm_file)
graphfile<-file.path(out,"cache","All_main_Harmony_SNN.rds")
if(!file.exists(graphfile)) {
 cat("FULL CELL NEIGHBORS",format(Sys.time()),"\n");flush.console()
 nn<-FindNeighbors(object=harm,k.param=20,compute.SNN=TRUE,nn.method="annoy",n.trees=50,verbose=TRUE)
 snn<-nn$snn;saveRDS(snn,graphfile,compress=FALSE);rm(nn);gc(verbose=FALSE)
} else snn<-readRDS(graphfile)
clusterfile<-file.path(out,"metadata","All_main_global_cluster_assignments.tsv")
if(!file.exists(clusterfile)) {
 cat("GLOBAL CLUSTERING",format(Sys.time()),"\n");flush.console()
 cl<-FindClusters(object=snn,resolution=c(.2,.4,.8,1.2),algorithm=1,random.seed=20261023,verbose=FALSE)
 cl<-as.data.table(cl,keep.rownames="cell_id")
 stopifnot(nrow(cl)==nrow(md));setnames(cl,setdiff(names(cl),"cell_id"),paste0("global_cluster_res_",c("0.2","0.4","0.8","1.2")))
 if(identical(cl$cell_id,as.character(seq_len(nrow(cl)))))cl$cell_id<-md$cell_id
 stopifnot(setequal(cl$cell_id,md$cell_id))
 fwrite(cl,clusterfile,sep="\t")
} else cl<-fread(clusterfile)
md<-merge(md,cl,by="cell_id",all.x=TRUE,sort=FALSE);md<-md[match(rownames(pca),cell_id)]
fwrite(md,file.path(out,"metadata","All_main_metadata_with_global_clusters.tsv"),sep="\t")
cluster_summary<-rbindlist(lapply(c("0.2","0.4","0.8","1.2"),function(res) {
 col<-paste0("global_cluster_res_",res)
 d<-md[,.(n=.N),by=.(cluster=get(col),cell_type)]
 d<-d[,.(n_cells=sum(n),majority_annotation=cell_type[which.max(n)],purity=max(n)/sum(n)),by=cluster]
 d[,resolution:=as.numeric(res)];d
}))
fwrite(cluster_summary,file.path(out,"results","Global_cluster_annotation_composition.tsv"),sep="\t")
trans<-rbindlist(lapply(seq_len(3),function(i) {
 a<-paste0("global_cluster_res_",c("0.2","0.4","0.8")[i]);b<-paste0("global_cluster_res_",c("0.4","0.8","1.2")[i])
 d<-md[,.(n_cells=.N),by=.(from_cluster=get(a),to_cluster=get(b))]
 d[,`:=`(from_resolution=a,to_resolution=b)];d
}))
fwrite(trans,file.path(out,"results","Global_resolution_transition_counts.tsv"),sep="\t")
rm(snn);gc(verbose=FALSE)
for(reduction in c("pca","harmony")) {
 path<-file.path(out,"cache",paste0("All_main_UMAP_",reduction,".rds"))
 if(!file.exists(path)) {
  cat("FULL CELL UMAP",reduction,format(Sys.time()),"\n");flush.console();set.seed(20261023)
  xx<-if(reduction=="pca")pca[,1:30,drop=FALSE] else harm
  um<-uwot::umap(xx,n_neighbors=30,n_components=2,metric="cosine",n_epochs=200,
    min_dist=.3,spread=1,scale=FALSE,init="spectral",n_threads=2,n_sgd_threads=1,
    ret_model=FALSE,verbose=TRUE)
  rownames(um)<-md$cell_id;saveRDS(um,path,compress=FALSE)
 } else um<-readRDS(path)
 tab<-data.table(cell_id=md$cell_id,UMAP_1=um[,1],UMAP_2=um[,2])
 fwrite(tab,file.path(out,"metadata",paste0("All_main_UMAP_",reduction,".tsv")),sep="\t")
 rm(um);gc(verbose=FALSE)
}
tsfile<-file.path(out,"cache","Balanced_landmark_tSNE_Harmony.rds")
if(!file.exists(tsfile)) {
 cat("LANDMARK tSNE",format(Sys.time()),"\n");flush.console();set.seed(20261023)
 ids<-pcmod$landmark_cell_ids
 ti<-Rtsne::Rtsne(harm[match(ids,md$cell_id),,drop=FALSE],dims=2,perplexity=30,theta=.5,
   pca=FALSE,check_duplicates=FALSE,max_iter=1000,num_threads=2,verbose=TRUE)
 ts<-data.table(cell_id=ids,tSNE_1=ti$Y[,1],tSNE_2=ti$Y[,2])
 saveRDS(ts,tsfile,compress=FALSE)
} else ts<-readRDS(tsfile)
fwrite(ts,file.path(out,"metadata","Balanced_landmark_tSNE_Harmony.tsv"),sep="\t")
saveRDS(list(metadata=md,PCA=pca,Harmony=harm,
 UMAP_PCA=readRDS(file.path(out,"cache","All_main_UMAP_pca.rds")),
 UMAP_Harmony=readRDS(file.path(out,"cache","All_main_UMAP_harmony.rds")),tSNE_landmarks=ts,
 PCA_model_file=model_file,full_gene_counts_by_library=file.path(out,"cache","canonical_by_library"),
 note="Full-cell embedding bundle, not a merged all-gene Seurat assay; original per-library full Seurat objects preserved."),
 file.path(out,"cache","All_main_atlas_embedding_bundle.rds"),compress=FALSE)
capture.output(sessionInfo(),file=file.path(out,"logs","02_sessionInfo.txt"))
writeLines("complete",file.path(out,"cache","Global_embedding.complete"))
cat("GLOBAL EMBEDDING COMPLETE",format(Sys.time()),"\n");sink()
