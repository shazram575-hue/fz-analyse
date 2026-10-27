# Assemble already fitted contrasts; no new tests or candidate gates.
options(stringsAsFactors=FALSE)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table);setDTthreads(2L)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/")
out<-file.path(root,"outputs","单细胞基础补全_2026-10-23")
genes<-fread(file.path(out,"metadata","Canonical_gene_annotation.tsv"))
ds<-fread(file.path(out,"results","Within_celltype_DE_summary.tsv"))
sc<-fread(file.path(out,"results","All_celltypes_celllevel_all_prefiltered.tsv"))
pb<-fread(file.path(out,"results","All_celltypes_pseudobulk_all_tested.tsv"))
pa<-fread(file.path(out,"results","All_celltypes_pseudobulk_age_gender_all_tested.tsv"))
universe<-rbindlist(lapply(ds$cell_type,function(ct)copy(genes)[,cell_type:=ct]))
sc<-sc[,.(cell_type,canonical_id,SC_log2FC=avg_log2FC,SC_pct_tumor=pct.1,SC_pct_healthy=pct.2,
 SC_P=p_val,SC_Bonferroni=p_val_adj,SC_DE_pass=DE_pass)]
pp<-pb[,.(cell_type,canonical_id,PB_log2FC=logFC,PB_P=P.Value,PB_BH=adj.P.Val,PB_DE_pass=DE_pass)]
aa<-pa[,.(cell_type,canonical_id,PB_adjusted_log2FC=logFC,PB_adjusted_P=P.Value,PB_adjusted_BH=adj.P.Val,PB_adjusted_DE_pass=DE_pass)]
full<-merge(universe,sc,by=c("cell_type","canonical_id"),all.x=TRUE,sort=FALSE)
full<-merge(full,pp,by=c("cell_type","canonical_id"),all.x=TRUE,sort=FALSE)
full<-merge(full,aa,by=c("cell_type","canonical_id"),all.x=TRUE,sort=FALSE)
full<-merge(full,ds[,.(cell_type,diagnostic_annotation,n_healthy_donors_min20,n_tumor_donors_min20,PB_status)],by="cell_type",all.x=TRUE,sort=FALSE)
full[,`:=`(SC_tested=!is.na(SC_P),PB_tested=!is.na(PB_P),PB_adjusted_tested=!is.na(PB_adjusted_P))]
full[,SC_PB_same_direction:=ifelse(SC_tested&PB_tested,sign(SC_log2FC)==sign(PB_log2FC),NA)]
stopifnot(nrow(full)==nrow(genes)*nrow(ds),!anyDuplicated(full[,.(cell_type,canonical_id)]))
fwrite(full,file.path(out,"results","Within_type_combined_gene_evidence.tsv"),sep="\t",na="NA")
sm<-full[,.(n_RNA_features=.N,SC_DE=sum(SC_DE_pass%in%TRUE),PB_DE=sum(PB_DE_pass%in%TRUE),
 PB_adjusted_DE=sum(PB_adjusted_DE_pass%in%TRUE),
 SC_PB_DE_same_direction=sum(SC_DE_pass%in%TRUE & PB_DE_pass%in%TRUE & SC_PB_same_direction%in%TRUE),
 PB_and_adjusted_DE_same_direction=sum(PB_DE_pass%in%TRUE & PB_adjusted_DE_pass%in%TRUE & sign(PB_log2FC)==sign(PB_adjusted_log2FC),na.rm=TRUE)),by=cell_type]
sm<-merge(ds,sm,by="cell_type",all.x=TRUE,sort=FALSE)
sm[PB_status!="complete",c("PB_DE","PB_adjusted_DE","SC_PB_DE_same_direction","PB_and_adjusted_DE_same_direction"):=NA_integer_]
fwrite(sm,file.path(out,"results","Within_type_DE_robustness_summary.tsv"),sep="\t",na="NA")
cat("Evidence rows",nrow(full),"\n");print(sm[,.(cell_type,PB_DE,PB_adjusted_DE,SC_PB_DE_same_direction)])
