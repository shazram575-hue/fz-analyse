options(stringsAsFactors=FALSE,warn=1,scipen=100)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(infercnv))
setDTthreads(2L)
set.seed(20261023)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
start<-Sys.time()
while(!file.exists(file.path(out,"cache","CNV_fixed_reference_v2.complete"))) {
  if(as.numeric(difftime(Sys.time(),start,units="hours"))>5)stop("Discovery design did not complete")
  Sys.sleep(15)
}
ref<-readRDS(file.path(out,"cache","CNV_fixed_reference_design_v2.rds"))
sel<-fread(file.path(out,"metadata","Discovery_patient_balanced_selected_cells.tsv"))
cases<-sort(unique(sel[analysis_role=="main_primary_total",patient_id]))
mainout<-file.path(out,"cache","CNV_discovery_by_patient")
dir.create(mainout,recursive=TRUE,showWarnings=FALSE)
statuses<-list();scores_all<-list();controls_all<-list();summaries<-list()
for(ip in seq_along(cases)) {
  pid<-cases[ip]
  pd<-file.path(mainout,pid)
  dir.create(pd,recursive=TRUE,showWarnings=FALSE)
  final<-file.path(pd,"CNV_cell_scores.tsv")
  if(file.exists(file.path(pd,"complete.flag"))) {
    cat("EXISTS CNV",pid,"\n")
    sc<-fread(final)
    ct<-fread(file.path(pd,"CNV_heldout_control_scores.tsv"))
    sm<-fread(file.path(pd,"CNV_patient_summary.tsv"))
  } else {
    cat("START CNV",ip,"/",length(cases),pid,format(Sys.time()),"\n");flush.console()
    md<-sel[patient_id==pid]
    parts<-list()
    for(g in unique(md$gsm)) {
      ids<-md[gsm==g,cell_id]
      z<-readRDS(file.path(out,"cache","epithelial_by_library",paste0(g,"_epithelial_counts.rds")))$counts
      parts[[g]]<-z[ref$genes,ids,drop=FALSE]
    }
    q<-do.call(cbind,parts)
    md<-md[match(colnames(q),cell_id)]
    # Technical groups only bound dendrogram cost with HMM=FALSE; no clone claims.
    qann<-data.table(cell_id=colnames(q),gsm=md$gsm,
       source_patient_id=pid,group=paste0("query__",pid,"__technical_",sprintf("%03d",ceiling(seq_len(ncol(q))/50))),role="query")
    ann<-rbindlist(list(ref$annotation,qann),fill=TRUE)
    x<-cbind(ref$counts,q)
    fwrite(ann,file.path(pd,"input_cell_annotations.tsv"),sep="\t")
    fwrite(ann[,.(cell_id,group)],file.path(pd,"infercnv_annotations.tsv"),sep="\t",col.names=FALSE)
    futile.logger::flog.appender(futile.logger::appender.file(file.path(pd,"infercnv_run.log")))
    obj<-CreateInfercnvObject(raw_counts_matrix=x,
       gene_order_file=file.path(out,"metadata","genome_reference","infercnv_discovery_frozen_gene_order_v2.tsv"),
       annotations_file=file.path(pd,"infercnv_annotations.tsv"),ref_group_names=ref$reference_groups)
    rm(x,q,z,parts);gc(verbose=FALSE)
    t0<-proc.time()
    # Healthy reference/calibration mean-count >=0.1 selection was frozen in16b.
    # cutoff=0 avoids selecting a different gene universe in every patient batch.
    res<-infercnv::run(obj,cutoff=0,min_cells_per_gene=1,out_dir=pd,HMM=FALSE,denoise=TRUE,
       analysis_mode="samples",cluster_by_groups=TRUE,cluster_references=FALSE,
       num_threads=2,plot_steps=FALSE,inspect_subclusters=FALSE,no_plot=TRUE,no_prelim_plot=TRUE,
       save_rds=FALSE,save_final_rds=TRUE,resume_mode=FALSE,write_expr_matrix=FALSE)
    elapsed<-(proc.time()-t0)["elapsed"]
    e<-res@expr.data
    stopifnot(all(is.finite(e)),all(e>0),setequal(rownames(e),ref$genes))
    ann<-ann[match(colnames(e),cell_id)]
    chr<-as.character(res@gene_order[rownames(e),1])
    chrlev<-paste0("chr",1:22)
    cm<-do.call(rbind,lapply(chrlev,function(k)apply(e[chr==k,,drop=FALSE],2,median)))
    rownames(cm)<-chrlev
    ann[,`:=`(query_context_patient=pid,cnv_rms=sqrt(colMeans((e-1)^2)),
              max_chromosome_shift=apply(abs(cm-1),2,max))]
    ctl<-ann[role=="healthy_calibration"]
    validation<-ann[role=="healthy_validation"]
    stopifnot(nrow(ctl)>=150L,nrow(validation)>=150L)
    q99r<-quantile(ctl$cnv_rms,.99,names=FALSE)
    q99c<-quantile(ctl$max_chromosome_shift,.99,names=FALSE)
    q95r<-quantile(ctl$cnv_rms,.95,names=FALSE)
    q95c<-quantile(ctl$max_chromosome_shift,.95,names=FALSE)
    ann[,`:=`(rms_control_q99=q99r,chromosome_control_q99=q99c,
               cnv_high_q99=cnv_rms>q99r & max_chromosome_shift>q99c,
               cnv_high_q95=cnv_rms>q95r & max_chromosome_shift>q95c)]
    sc<-merge(md,ann[role=="query",.(cell_id,cnv_rms,max_chromosome_shift,rms_control_q99,
                    chromosome_control_q99,cnv_high_q99,cnv_high_q95)],by="cell_id",all.x=TRUE,sort=FALSE)
    sc[is.na(cnv_high_q99),`:=`(cnv_high_q99=FALSE,cnv_high_q95=FALSE)]
    sc[,malignancy_call:=fifelse(is.na(cnv_rms),"Uncertain_insufficient_CNV_counts",
              fifelse(cnv_high_q99,"Malignant_CNV_high","Uncertain_tumor_epithelial"))]
    stopifnot(nrow(sc)==nrow(md),!anyDuplicated(sc$cell_id))
    ct<-ann[role%in%c("healthy_calibration","healthy_validation")]
    sm<-data.table(patient_id=pid,n_query=nrow(sc),n_malignant=sum(sc$cnv_high_q99),
         n_uncertain=sum(!sc$cnv_high_q99),n_q95=sum(sc$cnv_high_q95),
         n_reference=sum(ann$role=="reference"),n_calibration=sum(ct$role=="healthy_calibration"),
         n_validation=sum(ct$role=="healthy_validation"),
         calibration_above_threshold=mean(ct[role=="healthy_calibration",cnv_high_q99]),
         validation_above_threshold=mean(ct[role=="healthy_validation",cnv_high_q99]),
         rms_control_q99=q99r,chromosome_control_q99=q99c,run_seconds=unname(elapsed),n_cnv_genes=nrow(e))
    fwrite(sc,final,sep="\t")
    fwrite(ct,file.path(pd,"CNV_heldout_control_scores.tsv"),sep="\t")
    fwrite(sm,file.path(pd,"CNV_patient_summary.tsv"),sep="\t")
    saveRDS(list(chromosome_medians=cm,cell_annotation=ann,summary=sm),file.path(pd,"CNV_chromosome_summary.rds"))
    writeLines("complete",file.path(pd,"complete.flag"))
    cat("COMPLETE CNV",pid,"HIGH",sum(sc$cnv_high_q99),"/",nrow(sc),"SECONDS",elapsed,format(Sys.time()),"\n")
    rm(e,res,obj,cm,ann,md,ctl);gc(verbose=FALSE)
  }
  scores_all[[pid]]<-sc;controls_all[[pid]]<-ct;summaries[[pid]]<-sm
  statuses[[pid]]<-data.table(patient_id=pid,state="complete")
  fwrite(rbindlist(statuses),file.path(out,"audit","CNV_patient_run_status.tsv"),sep="\t")
  fwrite(rbindlist(summaries),file.path(out,"results","CNV_patient_summary.tsv"),sep="\t")
}
scores<-rbindlist(scores_all,fill=TRUE)
stopifnot(!anyDuplicated(scores$cell_id))
fwrite(scores,file.path(out,"results","Tumor_discovery_cells_CNV_calls.tsv"),sep="\t")
fwrite(rbindlist(controls_all),file.path(out,"audit","CNV_heldout_control_calibration_all_batches.tsv"),sep="\t")
writeLines("complete",file.path(out,"cache","CNV_all_discovery_patients.complete"))
cat("ALL CNV COMPLETE",format(Sys.time()),"QUERY",nrow(scores),"MALIGNANT",sum(scores$cnv_high_q99),"\n")
