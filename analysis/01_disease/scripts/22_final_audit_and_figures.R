# Descriptive figures and independent reconstruction; does not alter Pool A.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
rd<-function(n,where="results")fread(file.path(out,where,n))
pool<-rd("Pool_A_candidate_genes.tsv")
b<-rd("Bulk_all_tested_genes.tsv")
sc<-rd("SingleCell_main_all_positive_prefiltered.tsv")
gm<-rd("SingleCell_canonical_gene_metadata.tsv","metadata")
pm<-rd("Patient_pseudobulk_sample_metadata.tsv","metadata")
md<-rd("SingleCell_main_DE_cell_metadata.tsv","metadata")
cs<-rd("CNV_coherent_final_summary_by_patient.tsv")
val<-rd("CNV_heldout_validation_coherent_profiles.tsv","audit")
sn<-rd("Pool_A_sensitivity_retention.tsv")
bulk_ids<-sub("\\.[0-9]+(_PAR_Y)?$","\\1",b[abs(logFC)>.585 & P.Value<.05 & gene_type=="protein_coding",gene_id])
sc_ids<-sc[pct.1>.15 & avg_log2FC>.5 & p_val_adj<.05,canonical_id]
expected<-intersect(intersect(bulk_ids,sc_ids),gm[gene_type=="protein_coding",canonical_id])
cnv<-rd("Tumor_discovery_cells_CNV_coherent_final_calls.tsv")
checks<-data.table(check=c("independent_Pool_A_set_reconstruction","unique_Pool_A_ids",
 "CNV_calls_equal_DE_malignant_cells","malignant_count_matches_summary","normal_DE_count",
 "PB_donors_not_cells","full_feature_Bonferroni","all_sensitivity_primary_counts_match"),
 pass=c(setequal(expected,pool$canonical_id),!anyDuplicated(pool$canonical_id),
 setequal(cnv[cnv_high_q99==TRUE,cell_id],md[condition=="Malignant_Epithelial",cell_id]),
 sum(cs$n_malignant)==sum(md$condition=="Malignant_Epithelial"),sum(md$condition=="Normal_Epithelial")==12264,
 !anyDuplicated(pm$patient_id),max(abs(pmin(sc$p_val*nrow(gm),1)-sc$p_val_adj))<1e-10,
 all(sn$n_pool==nrow(pool))))
fwrite(checks,file.path(out,"audit","Final_independent_validation.tsv"),sep="\t")
stopifnot(all(checks$pass))
stats<-rbindlist(list(
 data.table(metric=c("Pool_A","Bulk_paired_same_direction_in_Pool_A","SC_malignant_donor_keys","SC_normal_donor_keys",
 "PB_malignant_donors","PB_normal_donors","Healthy_validation_unique_cells","Healthy_validation_unique_donors",
 "Healthy_validation_context_FPR_min","Healthy_validation_context_FPR_max","Healthy_validation_context_FPR_median",
 "Tumor_selected_cells","Tumor_CNV_uncertain_cells","Pool_A_PB_not_tested","Pool_A_PB_adjusted_not_tested"),
 value=c(nrow(pool),sum(pool$bulk_paired_direction_agrees,na.rm=TRUE),
 uniqueN(md[condition=="Malignant_Epithelial",patient_id]),uniqueN(md[condition=="Normal_Epithelial",patient_id]),
 sum(pm$include_pseudobulk & pm$condition=="Malignant_Epithelial"),
 sum(pm$include_pseudobulk & pm$condition=="Normal_Epithelial"),
 uniqueN(val$cell_id),uniqueN(val$source_patient_id),
 min(cs$validation_coherent_high_fraction),max(cs$validation_coherent_high_fraction),
 median(cs$validation_coherent_high_fraction),sum(cs$n_selected),sum(cs$n_uncertain),
 sum(is.na(pool$PB_main_BH)),sum(is.na(pool$PB_age_gender_BH))))))
fwrite(stats,file.path(out,"audit","Final_report_metrics.tsv"),sep="\t")
age<-pm[include_pseudobulk==TRUE,.(n=.N,age_min=min(clinical_age),age_median=median(clinical_age),
 age_max=max(clinical_age),female=sum(gender=="Female"),male=sum(gender=="Male")),by=condition]
fwrite(age,file.path(out,"audit","Pseudobulk_demographic_summary.tsv"),sep="\t")
# Tableau-family categorical colors, explicitly a local design choice.
pal<-c(blue="#4E79A7",red="#E15759",green="#59A14F",grey="#BAB0AC",orange="#F28E2B",purple="#B07AA1")
theme_set(theme_classic(base_size=11,base_family="sans")+theme(
 plot.title=element_text(face="bold",size=13),plot.subtitle=element_text(size=10,color="#444444"),
 plot.caption=element_text(size=8,color="#555555",hjust=0),legend.position="bottom",
 axis.text=element_text(color="#333333"),plot.margin=margin(12,16,12,12)))
savefig<-function(p,name,w,h) {
 ggsave(file.path(out,"figures",paste0(name,".pdf")),p,width=w,height=h,device="pdf",useDingbats=FALSE,bg="white")
 ggsave(file.path(out,"figures",paste0(name,".png")),p,width=w,height=h,dpi=220,bg="white")
}
z<-data.table(stage=factor(c("Bulk protein-coding DE","Malignant epithelial markers","Pool A"),
 levels=c("Bulk protein-coding DE","Malignant epithelial markers","Pool A")),n=c(7494,1838,nrow(pool)))
p<-ggplot(z,aes(stage,n,fill=stage))+geom_col(width=.58)+geom_text(aes(label=format(n,big.mark=",")),vjust=-.5,size=4)+
 scale_fill_manual(values=unname(c(pal["blue"],pal["red"],pal["green"])))+
 scale_y_continuous(expand=expansion(mult=c(0,.13)),labels=scales::comma)+
 labs(title="Expression-supported candidate pool",subtitle="TCGA-BRCA and GSE161529 | discovery stage",
 x=NULL,y="Number of genes",caption="Pool A: exact intersection plus GENCODE v36 protein-coding status.\nBulk DE is bidirectional; single-cell markers are positive. No pathway or drug-target filter.")+theme(legend.position="none")
savefig(p,"Final_01_candidate_pool",8,4.7)
labels<-c(without_doublet_flags="Exclude doublet flags",without_changed_libraries="Exclude two changed libraries",
 female_only="Female cells only",patient_half1="Within-donor cell half 1",patient_half2="Within-donor cell half 2",extended_normal="Add normal Epi libraries")
sn[,`:=`(label=unname(labels[sensitivity]),retention=100*n_retained/n_pool)]
sn[,label:=factor(label,levels=rev(unname(labels)))]
p<-ggplot(sn,aes(retention,label))+geom_segment(aes(x=0,xend=retention,yend=label),color="#D9E1E8",linewidth=3)+
 geom_point(color=pal["blue"],size=3)+geom_text(aes(label=sprintf("%d / %d (%.1f%%)",n_retained,n_pool,retention)),hjust=0,nudge_x=1,size=3.4)+
 scale_x_continuous(limits=c(0,126),breaks=c(0,25,50,75,100))+
 labs(title="Retention of the primary Pool A",subtitle="All alternatives use the same gene thresholds",x="Primary candidates retained (%)",y=NULL,
 caption="Within-donor halves divide cells, not patients: they assess sampling stability, not external replication.\nAlternative pools can contain additional genes; see the TSV for Jaccard similarity.")
savefig(p,"Final_02_sensitivity_retention",9,4.8)
csl<-melt(cs[,.(patient_id,Malignant=n_malignant,Uncertain=n_uncertain)],id.vars="patient_id",variable.name="call",value.name="cells")
csl[,patient_id:=factor(patient_id,levels=cs[order(n_malignant/n_selected),patient_id])]
p<-ggplot(csl,aes(cells,patient_id,fill=call))+geom_col(width=.76)+
 scale_fill_manual(values=c(Malignant=unname(pal["red"]),Uncertain=unname(pal["grey"])))+
 labs(title="CNV-supported epithelial classification",subtitle="32 tumor donor keys | at most 1,000 sampled epithelial cells per donor",
 x="Scored epithelial cells",y=NULL,fill=NULL,
 caption=sprintf("Primary call: RNA CNV amplitude plus coherent genomic profile; uncertain cells remain uncertain.\nHealthy holdout: 5 donors / 200 unique cells, reused across contexts; context call rate %.1f%%-%.1f%%.",
 100*min(cs$validation_coherent_high_fraction),100*max(cs$validation_coherent_high_fraction)))
savefig(p,"Final_03_CNV_donor_classification",8,9.2)
pool[,pb_category:=ifelse(is.na(PB_main_BH),"Not tested",ifelse(PB_main_support,"BH < 0.05, positive","Other tested"))]
p<-ggplot(pool[!is.na(PB_main_log2FC)],aes(sc_avg_log2FC,PB_main_log2FC,color=pb_category,shape=bulk_sc_direction_conflict))+
 geom_hline(yintercept=0,color="#AAAAAA",linetype=2)+geom_point(alpha=.65,size=1.65)+
 scale_color_manual(values=c("BH < 0.05, positive"=unname(pal["blue"]),"Other tested"=unname(pal["grey"])))+
 scale_shape_manual(values=c(`FALSE`=16,`TRUE`=17),labels=c(`FALSE`="Bulk up / SC up",`TRUE`="Bulk down / SC up"))+
 labs(title="Cell-level and patient-level expression evidence",subtitle="Pool A | patient pseudobulk supports are descriptive, not an extra inclusion gate",
 x="Single-cell average log2 fold change",y="Patient pseudobulk log2 fold change",color="Patient pseudobulk",shape="Expression direction",
 caption=sprintf("Pseudobulk: %d malignant and %d normal donors with at least 20 cells. %d Pool A genes not tested.\nCNV and expression evidence share the same RNA modality and do not establish drug binding or causality.",
 sum(pm$include_pseudobulk & pm$condition=="Malignant_Epithelial"),sum(pm$include_pseudobulk & pm$condition=="Normal_Epithelial"),sum(is.na(pool$PB_main_BH))))+
 theme(legend.box="vertical",legend.text=element_text(size=9),legend.title=element_text(size=10))+
 guides(color=guide_legend(order=1),shape=guide_legend(order=2))
savefig(p,"Final_04_cell_patient_evidence",8.4,6)
ct<-rd("SingleCell_broad_annotation_summary.tsv")
ct[,cohort:=factor(ifelse(analysis_role=="main_healthy_total","Healthy Total","Primary tumor Total"),levels=c("Healthy Total","Primary tumor Total"))]
cols<-c(Epithelial="#4E79A7",Fibroblast="#F28E2B",Immune="#59A14F",Endothelial="#76B7B2",Perivascular_SMC="#B07AA1",Mast="#EDC948",Mixed_epithelial_immune="#E15759",Unresolved="#BAB0AC")
p<-ggplot(ct,aes(cohort,n_cells,fill=cell_type))+geom_col(width=.55)+scale_fill_manual(values=cols)+
 scale_y_continuous(labels=scales::comma)+labs(title="Main-library epithelial and microenvironment atlas",
 subtitle="47 Total libraries after author QC | 45 donor keys | 236,309 cells",x=NULL,y="Cells",fill="Broad identity",
 caption="Library-wise multi-marker annotation with SingleR/HPCA support; mixed and unresolved labels retained.\nAtlas counts are not estimates of in vivo cell composition, and epithelial identity alone is not malignancy.")
savefig(p,"Final_05_broad_cell_atlas",8.4,5.5)
capture.output(sessionInfo(),file=file.path(out,"logs","22_R_sessionInfo.txt"))
print(stats);print(age);print(checks)
cat("FINAL AUDIT AND FIGURES COMPLETE",format(Sys.time()),"\n")
