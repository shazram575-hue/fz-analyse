options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
setDTthreads(2L);future::plan("sequential");set.seed(20261026)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
start<-Sys.time()
while(!file.exists(file.path(out,"cache","SingleCell_main_markers_and_pseudobulk.complete")) ||
      !file.exists(file.path(out,"cache","Normal_Epi_sensitivity_libraries.complete"))) {
 if(as.numeric(difftime(Sys.time(),start,units="hours"))>6)stop("Sensitivity inputs did not finish")
 Sys.sleep(15)
}
sink(file.path(out,"logs","20_extended_normal_sensitivity.log"),split=TRUE)
m<-fread(file.path(out,"metadata","GSE161529_reconciled_patient_metadata.tsv"))
norm<-m[analysis_role%in%c("main_healthy_total","normal_epi_sensitivity")]
nmd<-rbindlist(lapply(norm$gsm,function(g){
 z<-readRDS(file.path(out,"cache","epithelial_by_library",paste0(g,"_epithelial_counts.rds")))
 d<-as.data.table(z$metadata);d[,gsm:=g];d
}),fill=TRUE)
nmd<-merge(nmd,norm[,.(gsm,patient_id,clinical_age,gender,analysis_role)],by="gsm",all.x=TRUE)
nmd<-rbindlist(lapply(sort(unique(nmd$patient_id)),function(p){d<-nmd[patient_id==p];d[sample.int(nrow(d),min(nrow(d),1000L))]}))
nmd[,condition:="Normal_Epithelial"]
mp<-fread(file.path(out,"metadata","SingleCell_canonical_gene_identity_map.tsv"))
genes<-fread(file.path(out,"metadata","SingleCell_canonical_gene_metadata.tsv"))
A<-sparseMatrix(i=match(mp$canonical_id,genes$canonical_id),j=seq_len(nrow(mp)),x=1,
 dims=c(nrow(genes),nrow(mp)))
parts<-list()
for(g in unique(nmd$gsm)) {
 ids<-nmd[gsm==g,cell_id]
 z<-readRDS(file.path(out,"cache","epithelial_by_library",paste0(g,"_epithelial_counts.rds")))$counts
 parts[[g]]<-A%*%z[,ids,drop=FALSE]
 rm(z);gc(verbose=FALSE)
}
normcounts<-do.call(cbind,parts);rownames(normcounts)<-genes$canonical_id
rm(parts);gc(verbose=FALSE)
base<-readRDS(file.path(out,"cache","Main_malignant_normal_raw_counts.rds"))
cmd<-as.data.table(base$metadata)[condition=="Malignant_Epithelial"]
counts<-cbind(base$counts[,cmd$cell_id,drop=FALSE],normcounts)
rm(base,normcounts);gc(verbose=FALSE)
md<-rbindlist(list(cmd,nmd),fill=TRUE);md<-md[match(colnames(counts),cell_id)]
stopifnot(!anyDuplicated(md$cell_id),all(abs(Matrix::colSums(counts)-md$nCount_RNA)<1e-6))
so<-CreateSeuratObject(counts,min.cells=0,min.features=0)
so$condition<-md$condition
so<-NormalizeData(so,scale.factor=10000,verbose=FALSE);Idents(so)<-"condition"
z<-FindMarkers(so,ident.1="Malignant_Epithelial",ident.2="Normal_Epithelial",assay="RNA",slot="data",
 test.use="wilcox",logfc.threshold=.5,min.pct=.15,only.pos=TRUE,base=2,verbose=FALSE,random.seed=20261023)
z<-as.data.table(z,keep.rownames="canonical_id")
z[,marker_pass:=pct.1>.15 & avg_log2FC>.5 & p_val_adj<.05]
z<-merge(genes,z,by="canonical_id",all.y=TRUE,sort=FALSE)
fwrite(z,file.path(out,"results","SingleCell_extended_normal_Total_plus_Epi_all_positive_prefiltered.tsv"),sep="\t")
fwrite(z[marker_pass==TRUE],file.path(out,"results","SingleCell_extended_normal_Total_plus_Epi_eligible_markers.tsv"),sep="\t")
fwrite(md,file.path(out,"metadata","Extended_normal_sensitivity_cell_metadata.tsv"),sep="\t")
main<-fread(file.path(out,"results","SingleCell_main_all_positive_prefiltered.tsv"))
compare<-merge(main[,.(canonical_id,main_pct=pct.1)],z[,.(canonical_id,extended_pct=pct.1)],by="canonical_id")
stopifnot(max(abs(compare$main_pct-compare$extended_pct))<1e-8)
fwrite(data.table(check="Malignant cells held fixed; detection fractions identical",max_abs_difference=max(abs(compare$main_pct-compare$extended_pct))),
 file.path(out,"audit","Extended_normal_sensitivity_fixed_malignant_check.tsv"),sep="\t")
keys<-unique(paste(md$patient_id,md$condition,sep="::"))
G<-sparseMatrix(i=seq_len(nrow(md)),j=match(paste(md$patient_id,md$condition,sep="::"),keys),x=1,
 dims=c(nrow(md),length(keys)))
pb<-as.matrix(counts%*%G);colnames(pb)<-keys
pm<-md[,.(n_cells=.N,condition=condition[1],clinical_age=clinical_age[1],gender=gender[1]),
 by=.(pb_id=paste(patient_id,condition,sep="::"))]
pm<-pm[match(colnames(pb),pb_id)]
ix<-which(pm$n_cells>=20L)
s<-as.data.frame(pm[ix]);s$condition<-factor(s$condition,levels=c("Normal_Epithelial","Malignant_Epithelial"));s$gender<-factor(s$gender)
design<-if(nlevels(s$gender)>1L)model.matrix(~condition+clinical_age+gender,s) else model.matrix(~condition+clinical_age,s)
stopifnot(qr(design)$rank==ncol(design))
y<-DGEList(pb[,ix,drop=FALSE]);keep<-filterByExpr(y,design=design)
y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v<-voom(y,design,plot=FALSE)
fit<-eBayes(lmFit(v,design))
tt<-as.data.table(topTable(fit,coef="conditionMalignant_Epithelial",number=Inf,sort.by="none",confint=TRUE),keep.rownames="canonical_id")
tt<-merge(genes,tt,by="canonical_id",all.y=TRUE,sort=FALSE)
fwrite(tt,file.path(out,"results","Pseudobulk_extended_normals_age_gender_all_tested_genes.tsv"),sep="\t")
fwrite(pm,file.path(out,"metadata","Extended_normal_pseudobulk_samples.tsv"),sep="\t")
cat("EXTENDED NORMALS COMPLETE",format(Sys.time()),"NORMAL DONORS",uniqueN(nmd$patient_id),"NORMAL CELLS",nrow(nmd),
 "MALIGNANT CELLS",nrow(cmd),"MARKERS",sum(z$marker_pass),"\n")
capture.output(sessionInfo(),file=file.path(out,"logs","20_R_sessionInfo.txt"))
writeLines("complete",file.path(out,"cache","Extended_normal_sensitivity.complete"));sink()
