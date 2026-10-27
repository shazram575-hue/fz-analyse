# User scope: Tumor-origin vs Healthy within each existing annotation only.
# No identity-vs-rest tests. No revision to the frozen 813-gene Pool A.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
setDTthreads(2L);future::plan("sequential");set.seed(20261023)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","单细胞基础补全_2026-10-23")
stopifnot(file.exists(file.path(out,"cache","Preparation.complete")))
sink(file.path(out,"logs","03_within_celltype_DE.log"),split=TRUE)
p<-readRDS(file.path(out,"cache","All_celltypes_patient_pseudobulk_counts.rds"))
md<-fread(file.path(out,"metadata","All_main_cell_metadata.tsv"))
genes<-p$genes;m<-fread(file.path(out,"metadata","Main_library_patient_metadata.tsv"))
dir.create(file.path(out,"results","within_celltype"),showWarnings=FALSE)
cts<-sort(unique(md$cell_type));status<-list();allsc<-list();allpb<-list();alladj<-list()
for(ct in cts) {
 cat("TYPE START",ct,format(Sys.time()),"\n");flush.console()
 diagnostic<-ct%in%c("Mixed_epithelial_immune","Unresolved")
 s<-p$samples[cell_type==ct & eligible_min20==TRUE]
 nH<-sum(s$condition=="Healthy");nT<-sum(s$condition=="Tumor")
 pbstatus<-"not_run_insufficient_donors";nup<-ndown<-NA_integer_
 fit_pb<-function(adjust=FALSE) {
  ss<-as.data.frame(s);ss$condition<-factor(ss$condition,levels=c("Healthy","Tumor"));ss$gender<-factor(ss$gender)
  design<-if(!adjust)model.matrix(~condition,data=ss) else if(nlevels(ss$gender)>1L)
   model.matrix(~condition+clinical_age+gender,data=ss) else model.matrix(~condition+clinical_age,data=ss)
  if(qr(design)$rank<ncol(design) || nrow(design)-ncol(design)<3L)return(NULL)
  y<-DGEList(p$counts[,ss$donor_celltype,drop=FALSE]);keep<-filterByExpr(y,design=design)
  y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method="TMM");v<-voom(y,design,plot=FALSE)
  fit<-eBayes(lmFit(v,design));tt<-as.data.table(topTable(fit,coef="conditionTumor",number=Inf,sort.by="none",confint=TRUE),keep.rownames="canonical_id")
  tt<-merge(genes,tt,by="canonical_id",all.y=TRUE,sort=FALSE)
  tt[,`:=`(cell_type=ct,contrast="Tumor_origin_minus_Healthy_within_same_annotation",n_healthy_donors=nH,n_tumor_donors=nT,
   n_healthy_cells=sum(s[condition=="Healthy",n_cells]),n_tumor_cells=sum(s[condition=="Tumor",n_cells]),
   model=if(adjust)"condition_age_gender_if_variable" else "condition_only",
   diagnostic_annotation=diagnostic,FDR_family="BH within each cell type and model",
   DE_pass=abs(logFC)>.5 & adj.P.Val<.05)]
  tt
 }
 if(nH>=3L && nT>=3L) {
  tt<-fit_pb(FALSE)
  if(!is.null(tt)) {
   pbstatus<-"complete";allpb[[ct]]<-tt;nup<-tt[DE_pass==TRUE & logFC>0,.N];ndown<-tt[DE_pass==TRUE & logFC<0,.N]
   fwrite(tt,file.path(out,"results","within_celltype",paste0(ct,"_pseudobulk_all_tested.tsv")),sep="\t")
   fwrite(tt[DE_pass==TRUE],file.path(out,"results","within_celltype",paste0(ct,"_pseudobulk_DE_BH005_absFC05.tsv")),sep="\t")
  } else pbstatus<-"not_estimable_design_rank_or_df"
  ta<-fit_pb(TRUE)
  if(!is.null(ta)) {
   alladj[[ct]]<-ta
   fwrite(ta,file.path(out,"results","within_celltype",paste0(ct,"_pseudobulk_age_gender_all_tested.tsv")),sep="\t")
  }
 }
 # Balanced descriptive cell-level evidence: at most200 cells per donor/type,
 # selected before DE. All-cell pseudobulk above retains every eligible cell.
 cmeta<-md[cell_type==ct & balanced_selected==TRUE]
 parts<-list()
 for(g in unique(cmeta$gsm)) {
  z<-readRDS(file.path(out,"cache","balanced_by_library",paste0(g,".rds")))
  ids<-cmeta[gsm==g,cell_id]
  parts[[g]]<-z$counts[,ids,drop=FALSE]
  rm(z)
 }
 x<-do.call(cbind,parts);rm(parts);gc(verbose=FALSE)
 cmeta<-cmeta[match(colnames(x),cell_id)]
 so<-CreateSeuratObject(x,project=paste0("Within_",ct),min.cells=0,min.features=0)
 so$condition<-cmeta$condition;Idents(so)<-"condition"
 so<-NormalizeData(so,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
 st<-FindMarkers(so,ident.1="Tumor",ident.2="Healthy",assay="RNA",slot="data",test.use="wilcox",
   only.pos=FALSE,logfc.threshold=.5,min.pct=.15,base=2,random.seed=20261023,verbose=FALSE)
 st<-as.data.table(st,keep.rownames="canonical_id")
 stopifnot(max(abs(pmin(st$p_val*nrow(genes),1)-st$p_val_adj))<1e-10)
 st<-merge(genes,st,by="canonical_id",all.y=TRUE,sort=FALSE)
 st[,`:=`(cell_type=ct,contrast="Tumor_origin_minus_Healthy_within_same_annotation",
  DE_pass=abs(avg_log2FC)>.5 & pmax(pct.1,pct.2)>.15 & p_val_adj<.05,
  n_healthy_cells=sum(cmeta$condition=="Healthy"),n_tumor_cells=sum(cmeta$condition=="Tumor"),
  n_healthy_donors=uniqueN(cmeta[condition=="Healthy",patient_id]),n_tumor_donors=uniqueN(cmeta[condition=="Tumor",patient_id]),
  patient_evidence_limited=nH<3L | nT<3L,diagnostic_annotation=diagnostic,
  p_adjustment="Seurat Bonferroni across all33487 RNA features within contrast",
  celllevel_role="Descriptive donor-balanced cells; cells are not independent biological replicates")]
 allsc[[ct]]<-st
 fwrite(st,file.path(out,"results","within_celltype",paste0(ct,"_celllevel_all_prefiltered.tsv")),sep="\t")
 fwrite(st[DE_pass==TRUE],file.path(out,"results","within_celltype",paste0(ct,"_celllevel_DE_Bonf005_absFC05.tsv")),sep="\t")
 status[[ct]]<-data.table(cell_type=ct,diagnostic_annotation=diagnostic,n_healthy_donors_min20=nH,n_tumor_donors_min20=nT,
  PB_status=pbstatus,PB_up=nup,PB_down=ndown,celllevel_prefiltered=nrow(st),
  celllevel_up=st[DE_pass==TRUE & avg_log2FC>0,.N],celllevel_down=st[DE_pass==TRUE & avg_log2FC<0,.N],
  celllevel_healthy_cells=sum(cmeta$condition=="Healthy"),celllevel_tumor_cells=sum(cmeta$condition=="Tumor"))
 fwrite(rbindlist(status),file.path(out,"results","Within_celltype_DE_summary.tsv"),sep="\t",na="NA")
 rm(x,so,st);gc(verbose=FALSE)
 cat("TYPE COMPLETE",ct,pbstatus,format(Sys.time()),"\n");flush.console()
}
fwrite(rbindlist(allsc),file.path(out,"results","All_celltypes_celllevel_all_prefiltered.tsv"),sep="\t")
fwrite(rbindlist(allpb),file.path(out,"results","All_celltypes_pseudobulk_all_tested.tsv"),sep="\t")
fwrite(rbindlist(alladj),file.path(out,"results","All_celltypes_pseudobulk_age_gender_all_tested.tsv"),sep="\t")
writeLines(c("Only within-cell-type Tumor-origin vs Healthy tests were performed.",
 "Cell-level RNA tests use<=200 randomly selected cells per donor/type; reference counts are canonical all-gene counts, not integrated data.",
 "Pseudobulk uses all original main-cohort cells with donor/type>=20 cells and>=3 donors per condition.",
 "Mast/Perivascular have inadequate independent donors: descriptive cell-level evidence must not be treated as population inference.",
 "Mixed and Unresolved are diagnostic annotation groups, not established biological cell types.",
 "Tumor-origin Epithelial includes CNV-uncertain/unscored epithelium; not equivalent to the frozen malignant-only contrast.",
 "No gene biotype or pathway filter used. Pool A remains frozen."),file.path(out,"audit","Disease_DE_definition.txt"))
capture.output(sessionInfo(),file=file.path(out,"logs","03_sessionInfo.txt"))
writeLines("complete",file.path(out,"cache","Within_celltype_DE.complete"))
cat("ALL TYPES COMPLETE",format(Sys.time()),"\n");sink()
