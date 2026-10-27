# Broad annotation only. Epithelial origin is not a malignancy call.
annotate_library <- function(so, features, root, gsm) {
  out<-file.path(root,"outputs","第一阶段_2026-10-23")
  adir<-file.path(out,"audit","cell_annotation")
  edir<-file.path(out,"cache","epithelial_by_library")
  dir.create(adir,recursive=TRUE,showWarnings=FALSE)
  dir.create(edir,recursive=TRUE,showWarnings=FALSE)
  cl<-as.character(so$RNA_snn_res.0.4)
  lev<-sort(unique(cl))
  z<-Matrix::sparseMatrix(i=seq_along(cl),j=match(cl,lev),x=1,dims=c(length(cl),length(lev)))
  dat<-Seurat::GetAssayData(so,assay="RNA",layer="data")
  avg<-as.matrix(dat%*%z);avg<-sweep(avg,2,Matrix::colSums(z),"/")
  pct<-as.matrix((dat>0)%*%z);pct<-sweep(pct,2,Matrix::colSums(z),"/")
  rownames(avg)<-rownames(pct)<-rownames(dat);colnames(avg)<-colnames(pct)<-lev
  symbols<-features$gene_symbol[match(rownames(avg),features$ensembl_id)]
  ord<-order(symbols,-rowSums(avg),na.last=TRUE)
  use<-ord[!duplicated(symbols[ord]) & !is.na(symbols[ord]) & nzchar(symbols[ord])]
  test<-avg[use,,drop=FALSE];rownames(test)<-symbols[use]
  trained<-getOption("workflow3.hpca.trained")
  if(is.null(trained)) {
    ref<-readRDS(file.path(out,"cache","celldex_reference","HPCA_2024-02-26_reference.rds"))
    shared<-intersect(rownames(test),rownames(ref))
    stopifnot(length(shared)>10000L)
    trained<-SingleR::trainSingleR(ref=ref[shared,],labels=ref$label.main,
                 BPPARAM=BiocParallel::SerialParam())
    options(workflow3.hpca.trained=trained,workflow3.hpca.shared=shared)
    rm(ref);gc(verbose=FALSE)
  }
  shared<-getOption("workflow3.hpca.shared")
  pred<-SingleR::classifySingleR(test=test[shared,,drop=FALSE],trained=trained,
                                BPPARAM=BiocParallel::SerialParam())
  ptab<-data.table::data.table(cluster=rownames(pred),SingleR_label=as.character(pred$labels),
                               SingleR_pruned=as.character(pred$pruned.labels),SingleR_delta=pred$delta.next)
  panels<-list(Epithelial=c("EPCAM","KRT8","KRT18","KRT19","KRT7","KRT5","KRT14","KRT17","TP63"),
               Immune=c("PTPRC","LST1","TYROBP","FCER1G","CD3D","CD3E","TRAC","CD79A","MS4A1","NKG7"),
               Fibroblast=c("DCN","LUM","COL1A1","COL1A2","COL3A1","PDGFRA"),
               Endothelial=c("PECAM1","VWF","CLDN5","KDR","FLT1","EMCN"),
               Perivascular=c("RGS5","PDGFRB","CSPG4","ACTA2","TAGLN","MYH11"),
               Mast=c("TPSAB1","TPSB2","CPA3","KIT"),Erythroid=c("HBB","HBA1","HBA2","GYPA"))
  # Marker test is for cell identity, not for candidate-gene selection.
  lookup<-setNames(features$ensembl_id,features$gene_symbol)
  markerstats<-list();support<-matrix(0L,length(lev),length(panels),dimnames=list(lev,names(panels)))
  for(type in names(panels)) {
    sy<-panels[[type]];ids<-unname(lookup[sy]);ok<-!is.na(ids)&ids%in%rownames(avg)
    ids<-ids[ok];sy<-sy[ok]
    for(k in lev) {
      support[k,type]<-sum(pct[ids,k]>.20 & avg[ids,k]>.20)
      markerstats[[length(markerstats)+1L]]<-data.table::data.table(gsm=gsm,cluster=k,panel=type,
          gene_symbol=sy,ensembl_id=ids,mean_log_normalized=avg[ids,k],pct=pct[ids,k])
    }
  }
  ptab<-merge(ptab,data.table::as.data.table(support,keep.rownames="cluster"),by="cluster",sort=FALSE)
  ptab[,cell_type:="Unresolved"]
  ptab[Erythroid>=2L,cell_type:="Erythroid"]
  ptab[Immune>=2L,cell_type:="Immune"]
  ptab[Fibroblast>=2L & Immune<2L,cell_type:="Fibroblast"]
  ptab[Perivascular>=3L & Fibroblast<3L & Immune<2L & Epithelial<2L,cell_type:="Perivascular_SMC"]
  ptab[Endothelial>=2L & Immune<2L & (Epithelial<2L | SingleR_label=="Endothelial_cells"),cell_type:="Endothelial"]
  ptab[Mast>=2L & Epithelial<2L,cell_type:="Mast"]
  ptab[(Epithelial>=3L | (Epithelial>=2L & SingleR_label%in%c("Epithelial_cells","Keratinocytes"))) & Immune<2L,
       cell_type:="Epithelial"]
  ptab[Epithelial>=2L & Immune>=2L,cell_type:="Mixed_epithelial_immune"]
  ptab[,reference_marker_agreement:=fifelse(cell_type=="Epithelial",
    SingleR_label%in%c("Epithelial_cells","Keratinocytes"),NA)]
  ptab[,`:=`(gsm=gsm,annotation_cluster_resolution=.4,annotation_version="HPCA_marker_v1.1")]
  data.table::fwrite(ptab,file.path(adir,paste0(gsm,"_cluster_identity.tsv")),sep="\t")
  data.table::fwrite(data.table::rbindlist(markerstats),file.path(adir,paste0(gsm,"_lineage_markers.tsv")),sep="\t")
  md<-data.table::as.data.table(so[[]],keep.rownames="cell_id")
  md[,cell_type:=ptab$cell_type[match(as.character(RNA_snn_res.0.4),ptab$cluster)]]
  md[,SingleR_label:=ptab$SingleR_label[match(as.character(RNA_snn_res.0.4),ptab$cluster)]]
  md[,annotation_resolution:=.4]
  md[,annotation_version:="HPCA_marker_v1.1"]
  data.table::fwrite(md,file.path(adir,paste0(gsm,"_cell_identity.tsv")),sep="\t")
  epi<-md[cell_type=="Epithelial",cell_id]
  raw<-Seurat::GetAssayData(so,assay="RNA",layer="counts")
  saveRDS(list(counts=raw[,epi,drop=FALSE],metadata=md[cell_type=="Epithelial"],
               features=features,note="Epithelial-origin candidates; malignancy pending CNV"),
          file.path(edir,paste0(gsm,"_epithelial_counts.rds")),compress="gzip")
  summary<-md[,.(n_cells=.N,n_doublet_flag=sum(scDblFinder_class=="doublet")),by=cell_type]
  summary[,gsm:=gsm]
  data.table::fwrite(summary,file.path(adir,paste0(gsm,"_celltype_counts.tsv")),sep="\t")
  cat("BROAD ANNOTATION",gsm,"EPITHELIAL",length(epi),"OF",ncol(so),format(Sys.time()),"\n")
  invisible(summary)
}
