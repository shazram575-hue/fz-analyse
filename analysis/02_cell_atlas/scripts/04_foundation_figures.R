options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
force<-"force"%in%commandArgs(TRUE)
out<-file.path(root,"outputs","单细胞基础补全_2026-10-23")
old<-file.path(root,"outputs","第一阶段_2026-10-23")
stopifnot(file.exists(file.path(out,"cache","Global_embedding.complete")),file.exists(file.path(out,"cache","Within_celltype_DE.complete")))
sink(file.path(out,"logs","04_figures.log"),split=TRUE)
rd<-function(f,dir="results")fread(file.path(out,dir,f))
md<-rd("All_main_metadata_with_global_clusters.tsv","metadata")
for(n in grep("^global_cluster_res_",names(md),value=TRUE))md[,(n):=factor(as.character(get(n))) ]
colors<-c(Epithelial="#4E79A7",Fibroblast="#F28E2B",Immune="#59A14F",Endothelial="#76B7B2",Perivascular_SMC="#B07AA1",Mast="#EDC948",Mixed_epithelial_immune="#E15759",Unresolved="#BAB0AC")
cohort<-c(Healthy="#4E79A7",Tumor="#E15759")
lab<-c(Epithelial="Epithelial",Fibroblast="Fibroblast",Immune="Immune",Endothelial="Endothelial",Perivascular_SMC="Perivascular / SMC",Mast="Mast",Mixed_epithelial_immune="Mixed epithelial / immune",Unresolved="Unresolved")
theme_set(theme_classic(base_size=10,base_family="sans")+theme(plot.title=element_text(face="bold",size=12),
 plot.subtitle=element_text(size=9,color="#444444"),plot.caption=element_text(size=8,color="#555555",hjust=0),
 axis.text=element_text(color="#333333"),legend.position="bottom",legend.text=element_text(size=8),
 plot.margin=margin(10,12,10,10)))
savefig<-function(p,name,w=9,h=6) {
 if(!force && all(file.exists(file.path(out,"figures",paste0(name,c(".pdf",".png"))))))return(invisible(NULL))
 ggsave(file.path(out,"figures",paste0(name,".pdf")),p,width=w,height=h,device="pdf",useDingbats=FALSE,bg="white",limitsize=FALSE)
 ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=200,bg="white",limitsize=FALSE)
 cat("FIGURE",name,format(Sys.time()),"\n");flush.console()
}
# QC summaries show distributions without hiding outlier cells by deleting them.
lib<-md[,.(n_cells=.N,n_donors=uniqueN(patient_id),median_UMI=as.numeric(median(nCount_RNA)),median_features=as.numeric(median(nFeature_RNA)),
 median_mt=as.numeric(median(percent.mt)),doublet_fraction=mean(scDblFinder_class=="doublet")),by=.(gsm,author_sample,condition)]
raw<-fread(file.path(old,"audit","GSE161529_unfiltered_qc_summary.tsv"))
lib<-merge(lib,raw[,.(gsm,n_raw=n_cells)],by="gsm",all.x=TRUE);lib[,retained_fraction:=n_cells/n_raw]
fwrite(lib,file.path(out,"results","Library_QC_summary.tsv"),sep="\t")
lev<-lib[order(condition,author_sample),author_sample];md[,sample_plot:=factor(author_sample,levels=lev)]
qc<-md[,.(sample_plot,condition,`log10 UMI`=log10(nCount_RNA),`log10 detected genes`=log10(nFeature_RNA),`Mitochondrial UMIs (%)`=percent.mt)]
qc<-melt(qc,id.vars=c("sample_plot","condition"),variable.name="metric",value.name="value")
p<-ggplot(qc,aes(sample_plot,value,fill=condition))+geom_boxplot(outlier.shape=NA,width=.7,linewidth=.2)+
 facet_wrap(~metric,ncol=1,scales="free_y")+scale_fill_manual(values=cohort)+
 labs(title="Post-QC distributions across all 47 main libraries",x=NULL,y=NULL,fill="Source",
 caption="Author QC retained; no cells removed after inspecting these plots. Boxplots omit outlier symbols only.\nDoublet flags remain in the baseline. Counts refer to the full gene assay.")+
 theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6))
savefig(p,"Foundation_01_library_QC",13,9);rm(qc);gc(verbose=FALSE)
ret<-melt(lib[,.(gsm,author_sample,condition,Retained=n_cells,Excluded=n_raw-n_cells)],id.vars=c("gsm","author_sample","condition"),variable.name="QC_status",value.name="cells")
ret[,author_sample:=factor(author_sample,levels=lev)]
p<-ggplot(ret,aes(author_sample,cells,fill=QC_status))+geom_col(width=.8)+
 scale_fill_manual(values=c(Retained="#4E79A7",Excluded="#D0D0D0"))+scale_y_continuous(labels=scales::comma)+
 labs(title="Author-QC retention of current input matrices",x=NULL,y="Cells",fill=NULL,
 caption="Counts are from the provided matrices, including the two documented differences from the historical author table.")+
 theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6))
savefig(p,"Foundation_02_QC_retention",13,4.8)
# All-cell embeddings are rasterized only at the point layer in otherwise vector PDFs.
pca<-readRDS(file.path(out,"cache","All_main_PCA_embeddings.rds"))
hm<-readRDS(file.path(out,"cache","All_main_Harmony_embeddings.rds"))
up<-rd("All_main_UMAP_pca.tsv","metadata");uh<-rd("All_main_UMAP_harmony.tsv","metadata")
plotmd<-copy(md);plotmd[,`:=`(PC1=pca[,1],PC2=pca[,2],Harmony1=hm[,1],Harmony2=hm[,2])]
plotmd<-merge(plotmd,uh,by="cell_id",all.x=TRUE,sort=FALSE)
set.seed(20261023);plotmd<-plotmd[sample(.N)]
scatter<-function(d,x,y,group,title,subtitle=NULL,continuous=FALSE) {
 z<-ggplot(d,aes(x=.data[[x]],y=.data[[y]],color=.data[[group]]))+
  scattermore::geom_scattermore(pointsize=1,pixels=c(1600,1200))+
  labs(title=title,subtitle=subtitle,x=x,y=y,color=NULL)+theme(axis.text=element_blank(),axis.ticks=element_blank())
 if(group=="cell_type")z<-z+scale_color_manual(values=colors,labels=lab)+guides(color=guide_legend(ncol=2,override.aes=list(size=3)))
 else if(group=="condition")z<-z+scale_color_manual(values=cohort)+guides(color=guide_legend(override.aes=list(size=3)))
 else if(group=="Phase")z<-z+scale_color_manual(values=c(G1="#BAB0AC",S="#4E79A7",G2M="#E15759"))
 else if(continuous)z<-z+scale_color_viridis_c()
 if(startsWith(group,"global_cluster_")) {
  centers<-d[,.(xx=median(get(x)),yy=median(get(y))),by=.(cluster=get(group))]
  z<-z+guides(color="none")+ggrepel::geom_text_repel(data=centers,aes(xx,yy,label=cluster),inherit.aes=FALSE,
    size=2.2,color="#222222",seed=20261023,max.overlaps=Inf,box.padding=.15,min.segment.length=0)
 }
 z
}
p<-scatter(plotmd,"PC1","PC2","cell_type","Uncorrected PCA","All236,309 cells explicitly projected")+
 scatter(plotmd,"Harmony1","Harmony2","cell_type","Harmony coordinates","Library correction for geometry only")+plot_layout(guides="collect")
savefig(p,"Foundation_03_PCA_Harmony",12,5.7)
vv<-rd("PCA_variance_explained.tsv")
p<-ggplot(vv,aes(PC,variance_fraction*100))+geom_line(color="#4E79A7")+geom_point(color="#4E79A7",size=1.2)+
 geom_vline(xintercept=30,linetype=2,color="#888888")+
 labs(title="PCA scree plot",subtitle="50 PCs fitted to28,558 donor/type-balanced landmarks; first30 used downstream",
 x="Principal component",y="Variance in standardized HVGs (%)",caption="Consensus HVGs selected by across-library recurrence, independently of disease DE. No expression was corrected for the DE tests.")
savefig(p,"Foundation_04_PCA_scree",7.2,4.7)
p<-scatter(plotmd,"UMAP_1","UMAP_2","cell_type","Global UMAP: existing annotation","236,309 main cells; original labels preserved")+
 scatter(plotmd,"UMAP_1","UMAP_2","condition","Global UMAP: tissue source")+
 scatter(plotmd,"UMAP_1","UMAP_2","global_cluster_res_0.4","Global UMAP: clusters","Louvain resolution0.4; no automatic relabeling")+
 scatter(plotmd,"UMAP_1","UMAP_2","percent.mt","Global UMAP: mitochondrial fraction",continuous=TRUE)+plot_layout(ncol=2)
savefig(p,"Foundation_05_global_UMAP",12,10)
rawmd<-merge(md,up,by="cell_id",all.x=TRUE,sort=FALSE);set.seed(20261023);rawmd<-rawmd[sample(.N)]
p<-scatter(rawmd,"UMAP_1","UMAP_2","cell_type","Uncorrected UMAP: annotation")+
 scatter(plotmd,"UMAP_1","UMAP_2","cell_type","Harmony UMAP: annotation")+
 scatter(rawmd,"UMAP_1","UMAP_2","condition","Uncorrected UMAP: source")+
 scatter(plotmd,"UMAP_1","UMAP_2","condition","Harmony UMAP: source")+plot_layout(ncol=2,guides="collect")
savefig(p,"Foundation_06_integration_comparison",12,10)
ts<-rd("Balanced_landmark_tSNE_Harmony.tsv","metadata");ts<-merge(ts,md,by="cell_id",all.x=TRUE,sort=FALSE)
p<-scatter(ts,"tSNE_1","tSNE_2","cell_type","Balanced landmark tSNE",sprintf("%s sampled cells; no tSNE projection to the remaining cells",format(nrow(ts),big.mark=",")))+
 scatter(ts,"tSNE_1","tSNE_2","condition","Balanced landmark tSNE: source")+plot_layout(ncol=2)
savefig(p,"Foundation_07_landmark_tSNE",12,5.8)
# Per-library facets make donor segregation visible without assigning45 similar legend colors.
p<-ggplot(plotmd,aes(UMAP_1,UMAP_2,color=condition))+scattermore::geom_scattermore(pointsize=1,pixels=c(400,350))+
 facet_wrap(~author_sample,ncol=6)+scale_color_manual(values=cohort)+
 labs(title="Global Harmony UMAP, split by main library",x="UMAP1",y="UMAP2",color="Source",
 caption="Same global coordinates in every facet. Differences can reflect lineage, biology and technical composition.")+
 theme(axis.text=element_blank(),axis.ticks=element_blank(),strip.text=element_text(size=7))
savefig(p,"Foundation_08_UMAP_by_library",13,15)
p<-scatter(plotmd,"UMAP_1","UMAP_2","Phase","Cell-cycle phase on global UMAP","Descriptive module scores; no cell-cycle regression")
savefig(p,"Foundation_09_cell_cycle_UMAP",8.5,6.5)
cy<-md[,.(n=.N),by=.(condition,cell_type,Phase)];cy[,fraction:=n/sum(n),by=.(condition,cell_type)]
fwrite(cy,file.path(out,"results","Cellcycle_by_celltype_condition.tsv"),sep="\t")
p<-ggplot(cy,aes(cell_type,fraction,fill=Phase))+geom_col()+facet_wrap(~condition,ncol=1)+
 scale_fill_manual(values=c(G1="#BAB0AC",S="#4E79A7",G2M="#E15759"))+scale_y_continuous(labels=scales::percent)+
 scale_x_discrete(labels=lab)+labs(title="Descriptive cell-cycle composition",x=NULL,y="Cell fraction",fill=NULL)+
 theme(axis.text.x=element_text(angle=30,hjust=1))
savefig(p,"Foundation_10_cell_cycle_summary",10,6.5)
# Known marker dotplot is an annotation diagnostic, not a type-vs-rest DE test.
pan<-rd("Canonical_annotation_panel_by_donor.tsv")
pp<-pan[n_cells>=20,.(n_donors=.N,donor_mean_pct=mean(pct),donor_mean_log1p=mean(mean_log1p_normalized)),by=.(gene_symbol,condition,cell_type)]
fwrite(pp,file.path(out,"results","Annotation_dotplot_source.tsv"),sep="\t")
panel_order<-c("EPCAM","KRT8","KRT18","KRT19","KRT5","KRT14","TP63","ESR1","PGR","ELF5","KIT","PTPRC","CD3D","CD3E","TRAC","NKG7","GNLY","MS4A1","CD79A","JCHAIN","MZB1","LYZ","LST1","FCER1G","CD68","DCN","COL1A1","COL1A2","LUM","PECAM1","VWF","CLDN5","KDR","RGS5","PDGFRB","CSPG4","ACTA2","TAGLN","MYH11","TPSAB1","TPSB2","CPA3","MKI67","TOP2A","PCNA")
pp[,gene_symbol:=factor(gene_symbol,levels=panel_order)]
p<-ggplot(pp,aes(gene_symbol,cell_type,size=donor_mean_pct*100,color=donor_mean_log1p))+geom_point()+
 facet_wrap(~condition,ncol=1)+scale_color_gradient(low="#D9E6F2",high="#B2182B")+
 scale_size_area(max_size=5,breaks=c(0,25,50,75,100))+scale_y_discrete(labels=lab)+
 labs(title="Canonical-marker check of existing broad annotation",x=NULL,y=NULL,size="Donor-mean detection (%)",color="Donor-mean log1p expression",
 caption="Donor/type groups with>=20 cells. Expression is normalized by full-assay UMI totals. This is a descriptive marker panel, not a type-vs-rest DE analysis.")+
 theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7))
savefig(p,"Foundation_11_annotation_marker_dotplot",15,8.5)
comp<-md[,.(n=.N),by=.(patient_id,condition,cell_type)];comp[,fraction:=n/sum(n),by=patient_id]
fwrite(comp,file.path(out,"results","Donor_celltype_composition.tsv"),sep="\t")
p<-ggplot(comp,aes(patient_id,fraction,fill=cell_type))+geom_col(width=.9)+facet_grid(~condition,scales="free_x",space="free_x")+
 scale_fill_manual(values=colors,labels=lab)+scale_y_continuous(labels=scales::percent)+
 labs(title="Observed cell-type composition per donor",x=NULL,y="QC-cell fraction",fill=NULL,
 caption="Repeated libraries of the same donor are combined. Observed proportions depend on sampling and processing; no population abundance test is implied.")+
 theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7))+guides(fill=guide_legend(nrow=2))
savefig(p,"Foundation_12_donor_composition",13,6)
cc<-rd("Global_cluster_annotation_composition.tsv");ns<-cc[,.(n_clusters=.N),by=resolution]
p<-ggplot(ns,aes(resolution,n_clusters))+geom_line(color="#4E79A7")+geom_point(color="#4E79A7",size=2)+
 scale_x_continuous(breaks=c(.2,.4,.8,1.2))+labs(title="Global clustering across resolutions",x="Louvain resolution",y="Number of clusters",
 caption="Resolution0.4 is the baseline display. Original cell labels are preserved; majority labels do not automatically replace per-cell identities.")
savefig(p,"Foundation_13_resolution_summary",7,4.8)
cross<-md[,.(n=.N),by=.(cluster=global_cluster_res_0.4,cell_type)];cross[,fraction:=n/sum(n),by=cluster]
cross[,cluster:=factor(as.character(cluster),levels=as.character(sort(unique(as.integer(as.character(cluster))))))]
p<-ggplot(cross,aes(cluster,cell_type,fill=fraction))+geom_tile()+scale_fill_gradient(low="white",high="#4E79A7",labels=scales::percent)+
 scale_y_discrete(labels=lab)+labs(title="Annotation composition of global clusters",subtitle="Resolution0.4 | fraction within each cluster",x="Global cluster",y=NULL,fill="Fraction")+
 theme(axis.text.x=element_text(angle=90,size=7))
savefig(p,"Foundation_14_cluster_annotation_check",13,4.6)
ds<-rd("Within_celltype_DE_summary.tsv")
dd<-melt(ds[,.(cell_type,Up=PB_up,Down=PB_down)],id.vars="cell_type",variable.name="direction",value.name="genes")
p<-ggplot(dd[!is.na(genes)],aes(cell_type,genes,fill=direction))+geom_col(position="dodge",width=.75)+
 scale_fill_manual(values=c(Up="#E15759",Down="#4E79A7"))+scale_x_discrete(labels=lab)+
 labs(title="Within-type disease DE: patient pseudobulk",subtitle="BH<0.05 and absolute log2FC>0.5; all eligible original cells",x=NULL,y="Genes",fill=NULL,
 caption="Mast and Perivascular/SMC have inadequate donors on one side and no patient-level DE claim.\nMixed and Unresolved are diagnostic groups. Pseudobulk FDR is controlled within each cell type.")+
 theme(axis.text.x=element_text(angle=25,hjust=1))
savefig(p,"Foundation_15_within_type_PB_DE",10,6)
allpb<-rd("All_celltypes_pseudobulk_all_tested.tsv")
allpb[,direction:=ifelse(!DE_pass,"Other",ifelse(logFC>0,"Up","Down"))]
allpb[,neg_log10_BH:=pmin(-log10(pmax(adj.P.Val,1e-300)),50)]
p<-ggplot(allpb,aes(logFC,neg_log10_BH,color=direction))+scattermore::geom_scattermore(pointsize=1.3,pixels=c(700,600))+
 facet_wrap(~cell_type,ncol=3,scales="free_x",labeller=as_labeller(lab))+scale_color_manual(values=c(Up="#E15759",Down="#4E79A7",Other="#D0D0D0"))+
 geom_vline(xintercept=c(-.5,.5),linetype=2,color="#999999",linewidth=.3)+geom_hline(yintercept=-log10(.05),linetype=2,color="#999999",linewidth=.3)+
 labs(title="Within-type Tumor-origin vs Healthy: donor-level effects",x="Pseudobulk log2 fold change",y="-log10(BH), capped at50",color=NULL,
 caption="Condition-only model. Epithelial includes tumor-origin uncertain/unscored cells; it is not the frozen malignant-only contrast.")
savefig(p,"Foundation_16_pseudobulk_volcanoes",13,8.2)
# Heatmaps are display-only and never redefine the tested gene universe.
pb<-readRDS(file.path(out,"cache","All_celltypes_patient_pseudobulk_counts.rds"))
for(ct in ds[PB_status=="complete",cell_type]) {
 st<-allpb[cell_type==ct & DE_pass==TRUE]
 top<-unique(c(st[logFC>0][order(adj.P.Val,-abs(logFC))]$canonical_id[seq_len(min(10,sum(st$logFC>0)))],
   st[logFC<0][order(adj.P.Val,-abs(logFC))]$canonical_id[seq_len(min(10,sum(st$logFC<0)))]))
 top<-top[!is.na(top)];if(length(top)<2)next
 ss<-pb$samples[cell_type==ct & eligible_min20==TRUE][order(condition,patient_id)]
 y<-edgeR::calcNormFactors(edgeR::DGEList(pb$counts[,ss$donor_celltype,drop=FALSE]),method="TMM")
 z<-edgeR::cpm(y,log=TRUE,prior.count=2)[top,,drop=FALSE]
 z<-t(scale(t(z)));z[!is.finite(z)]<-0;z<-pmax(pmin(z,2.5),-2.5)
 ann<-data.frame(Source=ss$condition,row.names=ss$donor_celltype)
 sy<-pb$genes$gene_symbol[match(top,pb$genes$canonical_id)];sy[is.na(sy)|sy==""]<-top[is.na(sy)|sy==""]
 labels<-make.unique(sy)
 stem<-paste0("Foundation_17_",ct,"_donor_DE_heatmap")
 if(!force && all(file.exists(file.path(out,"figures",paste0(stem,c(".pdf",".png"))))))next
 drawheat<-function()pheatmap::pheatmap(z,color=colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100),
   breaks=seq(-2.5,2.5,length.out=101),cluster_cols=FALSE,cluster_rows=TRUE,labels_row=labels,
   show_colnames=FALSE,annotation_col=ann,annotation_colors=list(Source=cohort),border_color=NA,
   fontsize=9,main=paste(ct,": donor pseudobulk DE (row z-scores)"))
 pdf(file.path(out,"figures",paste0(stem,".pdf")),width=9,height=6.2,useDingbats=FALSE);drawheat();dev.off()
 png(file.path(out,"figures",paste0(stem,".png")),width=1800,height=1240,res=200,bg="white");drawheat();dev.off()
}
writeLines(c("All Foundation figures saved as PDF plus PNG.",
 "Scatter point layers rasterized for manageable vector-PDF sizes; labels remain vector.",
 "HVG/PCA-landmark sampling and tSNE subsets are explicitly labelled.",
 "Dotplot is a known-marker diagnostic, not a type-vs-rest DE test."),file.path(out,"audit","Figure_generation_notes.txt"))
capture.output(sessionInfo(),file=file.path(out,"logs","04_sessionInfo.txt"))
writeLines("complete",file.path(out,"cache","Figures.complete"))
cat("FIGURES COMPLETE",format(Sys.time()),"\n");sink()
