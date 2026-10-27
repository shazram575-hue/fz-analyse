options(stringsAsFactors=FALSE)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
start<-Sys.time()
while(!file.exists(file.path(out,"cache","Extended_normal_sensitivity.complete"))) {
 if(as.numeric(difftime(Sys.time(),start,units="hours"))>6)stop("Marker/sensitivity generation not complete")
 Sys.sleep(15)
}
genes<-fread(file.path(out,"metadata","SingleCell_canonical_gene_metadata.tsv"))
b<-fread(file.path(out,"results","Bulk_all_tested_genes.tsv"))
b[,canonical_id:=sub("\\.[0-9]+(_PAR_Y)?$","\\1",gene_id)]
stopifnot(!anyDuplicated(b$canonical_id))
b<-b[,.(canonical_id,bulk_gene_id=gene_id,bulk_gene_symbol=gene_name,bulk_gene_type=gene_type,
 bulk_log2FC=logFC,bulk_P=P.Value,bulk_BH=adj.P.Val,bulk_AveExpr=AveExpr,
 bulk_DE_pass=abs(logFC)>.585 & P.Value<.05)]
sc<-fread(file.path(out,"results","SingleCell_main_all_positive_prefiltered.tsv"))
expected<-pmin(sc$p_val*nrow(genes),1)
stopifnot(max(abs(expected-sc$p_val_adj),na.rm=TRUE)<1e-10)
sc<-sc[,.(canonical_id,sc_avg_log2FC=avg_log2FC,sc_pct_malignant=pct.1,sc_pct_normal=pct.2,
 sc_raw_P=p_val,sc_Bonferroni=p_val_adj,sc_marker_pass=marker_pass)]
full<-merge(genes,b,by="canonical_id",all.x=TRUE,sort=FALSE)
full<-merge(full,sc,by="canonical_id",all.x=TRUE,sort=FALSE)
paired<-fread(file.path(out,"results","Bulk_paired113_sensitivity_all_genes.tsv"))
paired<-paired[,.(canonical_id=ensembl_stable,bulk_paired_log2FC=logFC,
 bulk_paired_P=P.Value,bulk_paired_BH=adj.P.Val)]
stopifnot(!anyDuplicated(paired$canonical_id))
full<-merge(full,paired,by="canonical_id",all.x=TRUE,sort=FALSE)
full[,bulk_paired_direction_agrees:=!is.na(bulk_paired_log2FC)&sign(bulk_log2FC)==sign(bulk_paired_log2FC)]
full[,pool_A:=bulk_DE_pass%in%TRUE & sc_marker_pass%in%TRUE & gene_type%in%"protein_coding" & bulk_gene_type%in%"protein_coding"]
full[,bulk_sc_direction_conflict:=!is.na(bulk_log2FC) & bulk_log2FC<0 & sc_marker_pass%in%TRUE]
add_pb<-function(path,prefix) {
 d<-fread(path)
 d<-d[,.(canonical_id,logFC,P.Value,adj.P.Val)]
 setnames(d,c("logFC","P.Value","adj.P.Val"),paste0(prefix,c("_log2FC","_P","_BH")))
 d
}
for(v in list(c("Pseudobulk_main_all_tested_genes.tsv","PB_main"),
              c("Pseudobulk_age_gender_all_tested_genes.tsv","PB_age_gender"),
              c("Pseudobulk_extended_normals_age_gender_all_tested_genes.tsv","PB_extended_normal"))) {
 full<-merge(full,add_pb(file.path(out,"results",v[1]),v[2]),by="canonical_id",all.x=TRUE,sort=FALSE)
}
full[,PB_main_support:=!is.na(PB_main_BH)&PB_main_BH<.05&PB_main_log2FC>0]
full[,PB_adjusted_support:=!is.na(PB_age_gender_BH)&PB_age_gender_BH<.05&PB_age_gender_log2FC>0]
sens<-list(without_doublet_flags="SingleCell_without_doublet_flags_eligible_markers.tsv",
 without_changed_libraries="SingleCell_without_two_changed_libraries_eligible_markers.tsv",
 female_only="SingleCell_female_only_eligible_markers.tsv",patient_half1="SingleCell_patient_half1_eligible_markers.tsv",
 patient_half2="SingleCell_patient_half2_eligible_markers.tsv",extended_normal="SingleCell_extended_normal_Total_plus_Epi_eligible_markers.tsv")
for(n in names(sens)) {
 d<-fread(file.path(out,"results",sens[[n]]))
 full[,(paste0("marker_retained_",n)):=canonical_id%in%d$canonical_id]
}
full[,novelty_status:="Not evaluated in stage1; no claim of a novel drug target"]
md<-fread(file.path(out,"metadata","SingleCell_main_DE_cell_metadata.tsv"))
pm<-fread(file.path(out,"metadata","Patient_pseudobulk_sample_metadata.tsv"))
full[,`:=`(gencode_version="v36",sc_n_malignant_cells=sum(md$condition=="Malignant_Epithelial"),
 sc_n_normal_cells=sum(md$condition=="Normal_Epithelial"),
 sc_n_malignant_donor_keys=uniqueN(md[condition=="Malignant_Epithelial",patient_id]),
 sc_n_normal_donor_keys=uniqueN(md[condition=="Normal_Epithelial",patient_id]),
 PB_n_malignant_patients=sum(pm$include_pseudobulk & pm$condition=="Malignant_Epithelial"),
 PB_n_normal_patients=sum(pm$include_pseudobulk & pm$condition=="Normal_Epithelial"),
 cnv_classification="17c coherent q99; inferred RNA CNV, not DNA validation")]
full[,evidence_level:="Bulk differential RNA expression plus CNV-supported epithelial marker; functional/drug binding evidence not provided by this stage"]
fwrite(full,file.path(out,"results","Pool_A_full_gene_evidence_universe.tsv"),sep="\t",na="NA")
pool<-full[pool_A==TRUE][order(canonical_id)]
stopifnot(!anyDuplicated(pool$canonical_id),all(pool$gene_type=="protein_coding"),
 all(abs(pool$bulk_log2FC)>.585),all(pool$bulk_P<.05),all(pool$sc_avg_log2FC>.5),
 all(pool$sc_pct_malignant>.15),all(pool$sc_Bonferroni<.05))
fwrite(pool,file.path(out,"results","Pool_A_candidate_genes.tsv"),sep="\t",na="NA")
fwrite(pool[,.(canonical_id,gene_symbol)],file.path(out,"results","Pool_A_gene_ids_symbols.tsv"),sep="\t")
summary<-data.table(metric=c("Canonical_RNA_features","Bulk_tested","Bulk_exploratory_DE","Bulk_DE_protein_coding",
 "SingleCell_eligible_positive_markers","Pool_A","Pool_A_PB_main_supported","Pool_A_PB_age_gender_supported",
 "Pool_A_bulk_down_SC_up"),
 value=c(nrow(genes),nrow(b),sum(b$bulk_DE_pass),sum(b$bulk_DE_pass&b$bulk_gene_type=="protein_coding"),
 sum(sc$sc_marker_pass),nrow(pool),sum(pool$PB_main_support),sum(pool$PB_adjusted_support),sum(pool$bulk_sc_direction_conflict)))
fwrite(summary,file.path(out,"results","Pool_A_stage1_summary.tsv"),sep="\t")
sr<-rbindlist(lapply(names(sens),function(n) {
 ids<-fread(file.path(out,"results",sens[[n]]))$canonical_id
 alt<-full[bulk_DE_pass%in%TRUE & gene_type%in%"protein_coding" & bulk_gene_type%in%"protein_coding" & canonical_id%in%ids,canonical_id]
 data.table(sensitivity=n,n_pool=nrow(pool),n_alternative_pool=length(alt),
 n_retained=sum(pool$canonical_id%in%alt),jaccard=length(intersect(pool$canonical_id,alt))/length(union(pool$canonical_id,alt)))
}))
fwrite(sr,file.path(out,"results","Pool_A_sensitivity_retention.tsv"),sep="\t")
independent_expected<-intersect(b[bulk_DE_pass==TRUE & bulk_gene_type=="protein_coding",canonical_id],
 intersect(sc[sc_marker_pass==TRUE,canonical_id],genes[gene_type=="protein_coding",canonical_id]))
checks<-data.table(check=c("Pool_A_unique_ID","protein_coding_only","Bulk_thresholds","SC_thresholds","Bonferroni_definition","No_PB_hard_gate"),
 pass=c(!anyDuplicated(pool$canonical_id),all(pool$gene_type=="protein_coding"),
 all(abs(pool$bulk_log2FC)>.585&pool$bulk_P<.05),
 all(pool$sc_avg_log2FC>.5&pool$sc_pct_malignant>.15&pool$sc_Bonferroni<.05),
 max(abs(expected-sc$sc_Bonferroni),na.rm=TRUE)<1e-10,setequal(pool$canonical_id,independent_expected)))
stopifnot(all(checks$pass))
fwrite(checks,file.path(out,"audit","Pool_A_definition_checks.tsv"),sep="\t")
writeLines("complete",file.path(out,"cache","Pool_A_assembled.complete"))
print(summary);cat("POOL A ASSEMBLED",format(Sys.time()),"\n")
