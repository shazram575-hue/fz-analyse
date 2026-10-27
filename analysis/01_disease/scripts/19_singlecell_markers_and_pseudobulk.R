# User-defined marker criteria, with patient-level checks kept separate from the
# hard Pool A definition. No GO/KEGG or drug prior is used in this script.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
setDTthreads(2L);future::plan("sequential");set.seed(20261023)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
start<-Sys.time()
while(!file.exists(file.path(out,"cache","CNV_coherent_final_calls.complete"))) {
 if(as.numeric(difftime(Sys.time(),start,units="hours"))>6)stop("CNV classification not complete")
 Sys.sleep(15)
}
sink(file.path(out,"logs","19_markers_and_pseudobulk.log"),split=TRUE)
cat("START",format(Sys.time()),"\n")
sel<-fread(file.path(out,"metadata","Discovery_patient_balanced_selected_cells.tsv"))
cnv<-fread(file.path(out,"results","Tumor_discovery_cells_CNV_coherent_final_calls.tsv"))
normal<-sel[analysis_role=="main_healthy_total"]
normal[,condition:="Normal_Epithelial"]
malig<-cnv[cnv_high_q99==TRUE]
malig[,condition:="Malignant_Epithelial"]
md<-rbindlist(list(normal,malig),fill=TRUE)
stopifnot(!anyDuplicated(md$cell_id),nrow(malig)>100L,nrow(normal)>100L)
mp<-fread(file.path(out,"metadata","SingleCell_canonical_gene_identity_map.tsv"))
genes<-fread(file.path(out,"metadata","SingleCell_canonical_gene_metadata.tsv"))
geneids<-genes$canonical_id
collapse<-sparseMatrix(i=match(mp$canonical_id,geneids),j=seq_len(nrow(mp)),x=1,
                       dims=c(length(geneids),nrow(mp)))
parts<-list()
for(g in unique(md$gsm)) {
 ids<-md[gsm==g,cell_id]
 x<-readRDS(file.path(out,"cache","epithelial_by_library",paste0(g,"_epithelial_counts.rds")))$counts
 stopifnot(identical(rownames(x),mp$ensembl_id))
 parts[[g]]<-collapse%*%x[,ids,drop=FALSE]
 rm(x);gc(verbose=FALSE)
}
counts<-do.call(cbind,parts)
rownames(counts)<-geneids
rm(parts);gc(verbose=FALSE)
md<-md[match(colnames(counts),cell_id)]
stopifnot(all(abs(Matrix::colSums(counts)-md$nCount_RNA)<1e-6))
so<-CreateSeuratObject(counts=counts,project="Workflow3_patient_balanced_discovery",min.cells=0,min.features=0)
so$condition<-md$condition;so$patient_id<-md$patient_id;so$gsm<-md$gsm
so<-NormalizeData(so,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
Idents(so)<-"condition"
fwrite(md,file.path(out,"metadata","SingleCell_main_DE_cell_metadata.tsv"),sep="\t")
cat("MARKERS CELLS",ncol(so),"GENE UNIVERSE",nrow(so),"\n");flush.console()
run_markers<-function(object,name) {
 z<-FindMarkers(object,ident.1="Malignant_Epithelial",ident.2="Normal_Epithelial",assay="RNA",slot="data",
      test.use="wilcox",logfc.threshold=.5,min.pct=.15,only.pos=TRUE,base=2,
      random.seed=20261023,verbose=FALSE)
 z<-as.data.table(z,keep.rownames="canonical_id")
 stopifnot(all(c("p_val","p_val_adj","avg_log2FC","pct.1","pct.2")%in%names(z)))
 z[,marker_pass:=pct.1>.15 & avg_log2FC>.5 & p_val_adj<.05]
 z<-merge(genes,z,by="canonical_id",all.y=TRUE,sort=FALSE)
 z[,contrast:="Malignant_Epithelial_minus_Normal_Epithelial"]
 z[,p_adjustment:="Seurat Bonferroni over full canonical feature universe"]
 fwrite(z,file.path(out,"results",paste0(name,"_all_positive_prefiltered.tsv")),sep="\t")
 fwrite(z[marker_pass==TRUE],file.path(out,"results",paste0(name,"_eligible_markers.tsv")),sep="\t")
 cat("MARKER RESULT",name,nrow(z),"PASS",sum(z$marker_pass),format(Sys.time()),"\n")
 invisible(z)
}
main<-run_markers(so,"SingleCell_main")

# Aggregate raw UMIs by biological donor; cells are not independent patients.
groupkeys<-unique(paste(md$patient_id,md$condition,sep="::"))
design_map<-sparseMatrix(i=seq_len(nrow(md)),j=match(paste(md$patient_id,md$condition,sep="::"),groupkeys),
                      x=1,dims=c(nrow(md),length(groupkeys)))
pb<-as.matrix(counts%*%design_map)
colnames(pb)<-groupkeys
pm<-md[,.(n_cells=.N,patient_id=patient_id[1],condition=condition[1],
           clinical_age=clinical_age[1],gender=gender[1]),by=.(pb_id=paste(patient_id,condition,sep="::"))]
pm<-pm[match(colnames(pb),pb_id)]
pm[,include_pseudobulk:=n_cells>=20L]
fwrite(pm,file.path(out,"metadata","Patient_pseudobulk_sample_metadata.tsv"),sep="\t")
saveRDS(list(counts=pb,samples=pm,genes=genes),file.path(out,"cache","Patient_pseudobulk_raw_counts.rds"))
fit_pb<-function(pb,pm,name,adjust=FALSE) {
 ix<-which(pm$include_pseudobulk)
 s<-as.data.frame(pm[ix])
 s$condition<-factor(s$condition,levels=c("Normal_Epithelial","Malignant_Epithelial"))
 s$gender<-factor(s$gender)
 if(adjust && nlevels(s$gender)>1L) {
  design<-model.matrix(~condition+clinical_age+gender,data=s)
 } else if(adjust) design<-model.matrix(~condition+clinical_age,data=s)
 else design<-model.matrix(~condition,data=s)
 stopifnot(qr(design)$rank==ncol(design))
 y<-DGEList(pb[,ix,drop=FALSE]);keep<-filterByExpr(y,design=design)
 y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method="TMM")
 v<-voom(y,design,plot=FALSE)
 fit<-eBayes(lmFit(v,design))
 tt<-as.data.table(topTable(fit,coef="conditionMalignant_Epithelial",number=Inf,sort.by="none",confint=TRUE),keep.rownames="canonical_id")
 tt<-merge(genes,tt,by="canonical_id",all.y=TRUE,sort=FALSE)
 tt[,`:=`(model=name,n_normal_patients=sum(s$condition=="Normal_Epithelial"),
          n_malignant_patients=sum(s$condition=="Malignant_Epithelial"))]
 fwrite(tt,file.path(out,"results",paste0(name,"_all_tested_genes.tsv")),sep="\t")
 tt
}
pbmain<-fit_pb(pb,pm,"Pseudobulk_main")
pbadj<-fit_pb(pb,pm,"Pseudobulk_age_gender",TRUE)
cat("PSEUDOBULK COMPLETE",format(Sys.time()),"\n")

# Reversible robustness checks preserve the same gene thresholds and main result.
singlets<-md[scDblFinder_class!="doublet",cell_id]
if(length(singlets)>100L)run_markers(subset(so,cells=singlets),"SingleCell_without_doublet_flags")
without_changed<-md[!gsm%in%c("GSM4909296","GSM4909319"),cell_id]
run_markers(subset(so,cells=without_changed),"SingleCell_without_two_changed_libraries")
female<-md[gender=="Female",cell_id]
run_markers(subset(so,cells=female),"SingleCell_female_only")
set.seed(20261025)
half1<-unlist(lapply(split(md$cell_id,paste(md$patient_id,md$condition,sep="::")),function(z)sample(z,max(1L,floor(length(z)/2)))),use.names=FALSE)
half2<-setdiff(md$cell_id,half1)
run_markers(subset(so,cells=half1),"SingleCell_patient_half1")
run_markers(subset(so,cells=half2),"SingleCell_patient_half2")
saveRDS(so,file.path(out,"cache","Main_malignant_normal_DE_Seurat.rds"),compress="gzip")
saveRDS(list(counts=counts,metadata=md,genes=genes),file.path(out,"cache","Main_malignant_normal_raw_counts.rds"),compress="gzip")
capture.output(sessionInfo(),file=file.path(out,"logs","19_R_sessionInfo.txt"))
writeLines("complete",file.path(out,"cache","SingleCell_main_markers_and_pseudobulk.complete"))
cat("COMPLETE",format(Sys.time()),"\n")
sink()
