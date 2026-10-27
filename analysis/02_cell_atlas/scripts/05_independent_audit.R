options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table))
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","单细胞基础补全_2026-10-23")
rd<-function(f,d="results")fread(file.path(out,d,f))
md<-rd("All_main_metadata_with_global_clusters.tsv","metadata")
p<-readRDS(file.path(out,"cache","All_celltypes_patient_pseudobulk_counts.rds"))
pc<-readRDS(file.path(out,"cache","All_main_PCA_embeddings.rds"))
hm<-readRDS(file.path(out,"cache","All_main_Harmony_embeddings.rds"))
u<-rd("All_main_UMAP_harmony.tsv","metadata");up<-rd("All_main_UMAP_pca.tsv","metadata")
ts<-rd("Balanced_landmark_tSNE_Harmony.tsv","metadata")
sc<-rd("All_celltypes_celllevel_all_prefiltered.tsv");pb<-rd("All_celltypes_pseudobulk_all_tested.tsv")
pa<-rd("All_celltypes_pseudobulk_age_gender_all_tested.tsv");ds<-rd("Within_celltype_DE_summary.tsv")
gm<-rd("Canonical_gene_annotation.tsv","metadata")
bk<-md[,.(total_UMI=sum(nCount_RNA),n_cells=.N),by=donor_celltype]
bk<-bk[match(colnames(p$counts),donor_celltype)]
sc_err<-max(abs(pmin(sc$p_val*nrow(gm),1)-sc$p_val_adj))
pb_err<-max(pb[,max(abs(p.adjust(P.Value,method="BH")-adj.P.Val)),by=cell_type]$V1)
pa_err<-max(pa[,max(abs(p.adjust(P.Value,method="BH")-adj.P.Val)),by=cell_type]$V1)
sampling<-rd("Balanced_sampling_counts.tsv","audit")
proj<-rd("PCA_projection_validation.tsv","audit")
all_types<-sort(unique(md$cell_type))
checks<-data.table(check=c("All_main_cells_unique_236309","47_libraries_45_donor_keys","Balanced_subset_prespecified_cap200",
 "All_gene_counts_33487_canonical","Pseudobulk_UMI_conservation","Pseudobulk_cell_count_conservation",
 "All8_annotations_have_within_type_celllevel_results","Celllevel_Bonferroni_full_universe",
 "PB_BH_recomputed_per_type","Adjusted_PB_BH_recomputed_per_type","Celllevel_thresholds_recomputed",
 "PB_thresholds_recomputed","All_contrasts_within_type_Tumor_vs_Healthy","Rare_groups_not_given_patient_inference",
 "PCA_all_cells_finite_50dimensions","Harmony_all_cells_finite_30dimensions","All_cell_UMAP_complete_both_reductions",
 "Landmark_tSNE_only_prespecified_cells","PCA_projection_independently_checked","All_resolutions_assigned"),
 pass=c(nrow(md)==236309L && !anyDuplicated(md$cell_id),uniqueN(md$gsm)==47L && uniqueN(md$patient_id)==45L,
 all(sampling$n_selected<=200L) && sum(sampling$n_selected)==sum(md$balanced_selected),
 nrow(p$counts)==33487L && identical(rownames(p$counts),gm$canonical_id),
 max(abs(colSums(p$counts)-bk$total_UMI))<1e-6,sum(p$samples$n_cells)==nrow(md),
 identical(sort(unique(sc$cell_type)),all_types),sc_err<1e-10,pb_err<1e-10,pa_err<1e-10,
 all(sc$DE_pass==(abs(sc$avg_log2FC)>.5 & pmax(sc$pct.1,sc$pct.2)>.15 & sc$p_val_adj<.05)),
 all(pb$DE_pass==(abs(pb$logFC)>.5 & pb$adj.P.Val<.05)),
 all(sc$contrast=="Tumor_origin_minus_Healthy_within_same_annotation") && all(pb$contrast=="Tumor_origin_minus_Healthy_within_same_annotation"),
 all(ds[cell_type%in%c("Mast","Perivascular_SMC"),PB_status]=="not_run_insufficient_donors") && !any(pb$cell_type%in%c("Mast","Perivascular_SMC")),
 nrow(pc)==nrow(md) && ncol(pc)==50L && all(is.finite(pc)) && identical(rownames(pc),md$cell_id),
 nrow(hm)==nrow(md) && ncol(hm)==30L && all(is.finite(hm)) && identical(rownames(hm),md$cell_id),
 setequal(u$cell_id,md$cell_id) && setequal(up$cell_id,md$cell_id) && !anyNA(u) && !anyNA(up),
 setequal(ts$cell_id,md[balanced_selected==TRUE,cell_id]) && !anyNA(ts),
 proj[check=="PCA_landmark_projection_max_error",value]<1e-4 && proj[check=="PCA_landmark_projection_relative_error",value]<1e-6,
 all(vapply(grep("^global_cluster_res_",names(md),value=TRUE),function(n)!anyNA(md[[n]]),logical(1)))))
fwrite(checks,file.path(out,"audit","Independent_validation.tsv"),sep="\t")
stopifnot(all(checks$pass))
res<-rd("Global_cluster_annotation_composition.tsv")
metrics<-data.table(metric=c("Main_cells","Libraries","Donor_keys","Biological_broad_types","Diagnostic_groups",
 "PCA_training_landmarks","Canonical_RNA_features","PCA_HVGs","Celllevel_within_type_contrasts",
 "Patient_PB_contrasts_completed","Patient_PB_biological_types","Patient_PB_diagnostic_groups","Patient_PB_insufficient_donor_types",
 paste0("Global_clusters_resolution_",c(.2,.4,.8,1.2)),"PCA_projection_max_error","Celllevel_Bonferroni_max_error","PB_BH_max_error"),
 value=c(nrow(md),uniqueN(md$gsm),uniqueN(md$patient_id),6L,2L,sum(md$balanced_selected),nrow(gm),
 nrow(readRDS(file.path(out,"cache","Balanced_landmark_PCA_model.rds"))$model$rotation),
 uniqueN(sc$cell_type),ds[PB_status=="complete",.N],ds[PB_status=="complete" & diagnostic_annotation==FALSE,.N],
 ds[PB_status=="complete" & diagnostic_annotation==TRUE,.N],ds[PB_status!="complete",.N],
 vapply(c(.2,.4,.8,1.2),function(r)res[resolution==r,.N],integer(1)),
 proj[check=="PCA_landmark_projection_max_error",value],sc_err,pb_err))
fwrite(metrics,file.path(out,"audit","Completion_metrics.tsv"),sep="\t")
capture.output(sessionInfo(),file=file.path(out,"logs","05_sessionInfo.txt"))
writeLines("complete",file.path(out,"cache","Independent_validation.complete"))
print(checks);print(metrics)
