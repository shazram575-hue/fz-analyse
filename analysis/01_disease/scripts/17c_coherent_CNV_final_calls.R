# Final conservative RNA-CNV annotation. A whole-chromosome shift is not required:
# focal/arm CNAs can have a near-zero chromosome median. Instead require a high
# calibrated CNV burden AND a coherent genomic profile, with leave-one-cell-out
# centroids and a separate healthy-validation donor set. No gene DE has run yet.
options(stringsAsFactors=FALSE,warn=1,scipen=100)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table);suppressPackageStartupMessages(library(infercnv))
setDTthreads(2L)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
t0<-Sys.time()
while(!file.exists(file.path(out,"cache","CNV_all_discovery_patients.complete"))) {
  if(as.numeric(difftime(Sys.time(),t0,units="hours"))>5)stop("CNV signal generation did not finish")
  Sys.sleep(15)
}
allsc<-fread(file.path(out,"results","Tumor_discovery_cells_CNV_calls.tsv"))
primary<-list();validations<-list();summaries<-list()
for(pid in unique(allsc$patient_id)) {
  cat("COHERENCE",pid,format(Sys.time()),"\n");flush.console()
  pd<-file.path(out,"cache","CNV_discovery_by_patient",pid)
  obj<-readRDS(file.path(pd,"run.final.infercnv_obj"))
  delta<-obj@expr.data-1
  ann<-fread(file.path(pd,"input_cell_annotations.tsv"))
  ann<-ann[match(colnames(delta),cell_id)]
  sc<-allsc[patient_id==pid]
  qix<-match(sc$cell_id,colnames(delta))
  aidx<-which(ann$role=="healthy_calibration")
  bidx<-which(ann$role=="healthy_validation")
  sc[,`:=`(cnv_profile_cosine=NA_real_,cosine_calibration_q99=NA_real_,cosine_calibration_q95=NA_real_,
            profile_support_group=NA_character_,n_cells_in_support_group=0L)]
  key<-paste(sc$gsm,sc$RNA_snn_res.0.4,sep="__")
  v_any99<-v_any95<-rep(FALSE,length(bidx))
  for(k in unique(key)) {
    members<-which(key==k & !is.na(qix))
    if(!length(members))next
    support<-members
    support_name<-k
    if(length(support)<20L) {
      support<-which(!is.na(qix));support_name<-paste0(pid,"__patient_fallback")
    }
    if(length(support)<20L)next
    x<-delta[,qix[support],drop=FALSE]
    ss<-rowSums(x)
    norm2<-colSums(x*x)
    dot<-colSums(x*ss)
    loo_cos<-(dot-norm2)/sqrt(pmax(norm2*(sum(ss*ss)-2*dot+norm2),.Machine$double.eps))
    mu<-ss/length(support)
    munorm<-sqrt(sum(mu*mu))
    ca<-delta[,aidx,drop=FALSE];va<-delta[,bidx,drop=FALSE]
    ccos<-colSums(ca*mu)/pmax(sqrt(colSums(ca*ca))*munorm,.Machine$double.eps)
    vcos<-colSums(va*mu)/pmax(sqrt(colSums(va*va))*munorm,.Machine$double.eps)
    cut99<-max(.20,quantile(ccos,.99,na.rm=TRUE,names=FALSE))
    cut95<-max(.15,quantile(ccos,.95,na.rm=TRUE,names=FALSE))
    sc$cnv_profile_cosine[members]<-loo_cos[match(members,support)]
    sc$cosine_calibration_q99[members]<-cut99
    sc$cosine_calibration_q95[members]<-cut95
    sc$profile_support_group[members]<-support_name
    sc$n_cells_in_support_group[members]<-length(support)
    v_any99<-v_any99 | vcos>cut99
    v_any95<-v_any95 | vcos>cut95
  }
  control<-fread(file.path(pd,"CNV_heldout_control_scores.tsv"))
  cal<-control[role=="healthy_calibration"]
  q99<-quantile(cal$cnv_rms,.99,names=FALSE)
  q95<-quantile(cal$cnv_rms,.95,names=FALSE)
  sc[,`:=`(cnv_strict_whole_chromosome=cnv_high_q99,
             cnv_rms_only_q99=!is.na(cnv_rms)&cnv_rms>q99)]
  sc[,cnv_high_q99:=!is.na(cnv_rms)&!is.na(cnv_profile_cosine)&cnv_rms>q99&cnv_profile_cosine>cosine_calibration_q99]
  sc[,cnv_high_q95:=!is.na(cnv_rms)&!is.na(cnv_profile_cosine)&cnv_rms>q95&cnv_profile_cosine>cosine_calibration_q95]
  sc[,malignancy_call:=fifelse(cnv_high_q99,"Malignant_CNV_high","Uncertain_tumor_epithelial")]
  sc[,classification_basis:="Tumor-origin epithelial + CNV RMS>healthy q99 + leave-one-out genomic-profile coherence>healthy q99 (minimum cosine0.20)"]
  val_ids<-colnames(delta)[bidx]
  vm<-control[match(val_ids,cell_id)]
  vm[,`:=`(coherent_q99=cnv_rms>q99&v_any99,coherent_q95=cnv_rms>q95&v_any95,
             rms_only_q99=cnv_rms>q99)]
  validations[[pid]]<-vm
  primary[[pid]]<-sc
  summaries[[pid]]<-data.table(patient_id=pid,n_selected=nrow(sc),n_scored=sum(!is.na(sc$cnv_rms)),
     n_malignant=sum(sc$cnv_high_q99),n_uncertain=sum(!sc$cnv_high_q99),n_q95=sum(sc$cnv_high_q95),
     n_rms_only=sum(sc$cnv_rms_only_q99),n_whole_chromosome_strict=sum(sc$cnv_strict_whole_chromosome),
     healthy_validation_cells=nrow(vm),healthy_validation_donors=uniqueN(vm$source_patient_id),
     validation_coherent_high_fraction=mean(vm$coherent_q99),validation_rms_only_high_fraction=mean(vm$rms_only_q99))
  rm(obj,delta,x,ca,va);gc(verbose=FALSE)
}
calls<-rbindlist(primary,fill=TRUE)
stopifnot(!anyDuplicated(calls$cell_id))
fwrite(calls,file.path(out,"results","Tumor_discovery_cells_CNV_coherent_final_calls.tsv"),sep="\t")
fwrite(rbindlist(summaries),file.path(out,"results","CNV_coherent_final_summary_by_patient.tsv"),sep="\t")
fwrite(rbindlist(validations),file.path(out,"audit","CNV_heldout_validation_coherent_profiles.tsv"),sep="\t")
writeLines(c("PRIMARY: epithelial identity, tumor origin, CNV RMS above calibration healthy q99, and coherent genomic profile.",
"Profile cosine uses leave-one-cell-out centroids of within-library RNA clusters; groups with <20 cells use the same patient's pooled profile.",
"Cosine cutoffs are calibration-healthy q99 with a pre-DE minimum0.20; q95/minimum0.15 is sensitivity only.",
"Whole-chromosome median gating is a strict sensitivity and does not define the primary calls, preserving arm/focal CNA sensitivity.",
"No differential genes or Pool A sizes were used to choose cell-classification thresholds.",
"Uncertain tumor epithelial cells are retained as uncertain; no malignancy labels are projected to unscored atlas cells.",
"Healthy validation comprises5 donors/200 unique cells, reused across query contexts; contexts are not independent validation cohorts."),
 file.path(out,"audit","CNV_PRIMARY_CLASSIFICATION_METHOD.txt"))
writeLines("complete",file.path(out,"cache","CNV_coherent_final_calls.complete"))
cat("COHERENT CNV COMPLETE",nrow(calls),"HIGH",sum(calls$cnv_high_q99),format(Sys.time()),"\n")
