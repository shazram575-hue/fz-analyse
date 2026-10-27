# Primary RMS calibration does not require whole-chromosome shifts, which would
# exclude arm-level/focal CNAs. The whole-chromosome gate is retained as sensitivity.
# This clarification is made before any scRNA DE/Pool A analysis, not to tune its size.
options(stringsAsFactors=FALSE)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
start<-Sys.time()
while(!file.exists(file.path(out,"cache","CNV_all_discovery_patients.complete"))) {
  if(as.numeric(difftime(Sys.time(),start,units="hours"))>5)stop("CNV completion wait exceeded five hours")
  Sys.sleep(15)
}
sc<-fread(file.path(out,"results","Tumor_discovery_cells_CNV_calls.tsv"))
ct<-fread(file.path(out,"audit","CNV_heldout_control_calibration_all_batches.tsv"))
q95<-ct[role=="healthy_calibration",.(rms_q95=quantile(cnv_rms,.95)),by=query_context_patient]
sc[,cnv_high_rms_and_whole_chromosome:=cnv_high_q99]
sc[,cnv_high_q99:=!is.na(cnv_rms) & cnv_rms>rms_control_q99]
sc[,rms_control_q95:=q95$rms_q95[match(patient_id,q95$query_context_patient)]]
sc[,cnv_high_q95:=!is.na(cnv_rms) & cnv_rms>rms_control_q95]
sc[,malignancy_call:=fifelse(is.na(cnv_rms),"Uncertain_insufficient_CNV_counts",
 fifelse(cnv_high_q99,"Malignant_CNV_high","Uncertain_tumor_epithelial"))]
sc[,classification_basis:="Epithelial identity + tumor origin + RNA-inferred CNV RMS above healthy-calibration q99; no claim of DNA validation"]
fwrite(sc,file.path(out,"results","Tumor_discovery_cells_CNV_final_calls.tsv"),sep="\t")
ct[,cnv_rms_above_q99:=cnv_rms>rms_control_q99]
qc<-ct[,.(n_cells=.N,n_healthy_donors=uniqueN(source_patient_id),
           above_rms_threshold_fraction=mean(cnv_rms_above_q99)),by=.(query_context_patient,role)]
fwrite(qc,file.path(out,"audit","CNV_independent_healthy_validation_by_context.tsv"),sep="\t")
summary<-sc[,.(n_CNV_evaluated=sum(!is.na(cnv_rms)),n_selected=.N,n_malignant=sum(cnv_high_q99),
 n_uncertain=sum(!cnv_high_q99),n_looser_q95=sum(cnv_high_q95),n_strict_whole_chr=sum(cnv_high_rms_and_whole_chromosome)),by=patient_id]
fwrite(summary,file.path(out,"results","CNV_final_classification_by_patient.tsv"),sep="\t")
writeLines(c("Primary cell calls use calibrated CNV RMS; whole-chromosome median shifts are diagnostic/sensitivity only.",
"A whole-chromosome requirement would miss arm-level and focal CNAs and was removed before gene DE/Pool A computation.",
"Calibration and validation use different healthy donors. Reused control cells across patient contexts are not independent biological replicates.",
"Low-CNV tumor-origin epithelial cells remain uncertain. Unscored atlas cells are not assigned projected malignancy labels."),
 file.path(out,"audit","CNV_classification_definition.txt"))
writeLines("complete",file.path(out,"cache","CNV_final_calls.complete"))
cat("FINAL CALLS",nrow(sc),"MALIGNANT",sum(sc$cnv_high_q99),format(Sys.time()),"\n")
