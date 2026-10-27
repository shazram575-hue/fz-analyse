options(stringsAsFactors=FALSE,warn=1)
.libPaths(c('D:/R4.5.2/R-4.5.2/library',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(Matrix)})
setDTthreads(2L);set.seed(20261026)
root<-normalizePath(commandArgs(TRUE)[1],winslash='/');out<-file.path(root,'outputs/MFSD10_表达分布摸底_2026-10-26')
sink(file.path(out,'logs/02_figures.log'),split=TRUE)
cells<-readRDS(file.path(out,'cache/MFSD10_cell_level.rds'));b<-fread(file.path(out,'results/Patient_broad_MFSD10_summary.tsv'));e<-fread(file.path(out,'results/Patient_epithelial_MFSD10_summary.tsv'));m<-fread(file.path(out,'metadata/Library_clinical_context.tsv'))
ovfile<-file.path(out,'results/Patient_balanced_celltype_overview.tsv');ov<-fread(ovfile)
for(v in c('LOO_median_CPM_min','LOO_median_CPM_max'))ov[!is.finite(get(v)),(v):=NA_real_]
fwrite(ov,ovfile,sep='\t')
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),strip.text=element_text(size=9),legend.position='bottom',plot.caption=element_text(hjust=0,size=8)))
color_max<-max(b[cell_type!='Not_annotated',log2_CPM_plus1])
savefig<-function(p,name,w,h){ggsave(file.path(out,'figures',paste0(name,'.png')),p,width=w,height=h,dpi=160,bg='white',limitsize=FALSE);ggsave(file.path(out,'figures',paste0(name,'.pdf')),p,width=w,height=h,limitsize=FALSE)}
cts<-c('Epithelial','Fibroblast','Endothelial','Immune','Mast','Perivascular_SMC','Mixed_epithelial_immune','Unresolved')
pal<-c('Malignant_CNV_coherent'='#B34043','Uncertain_CNV_scored'='#A77F22','Tumor_epithelial_unscored'='#8B8B8B','Healthy_context_no_CNV_call'='#237DAD')
state_labels<-c('Malignant_CNV_coherent'='CNV-coherent malignant','Uncertain_CNV_scored'='CNV-scored uncertain','Tumor_epithelial_unscored'='Tumor epithelium: unscored','Healthy_context_no_CNV_call'='Healthy tissue: no CNV call')
for(prep in c('Total','Epithelial_enriched')){
 d<-copy(b[preparation==prep & tissue_context%in%c('Primary_tumor','Healthy_breast') & cell_type!='Not_annotated'])
 ctxlabels<-m[,.(patient_clinical_labels=paste(sort(unique(clinical_label)),collapse=';')),by=.(patient_id,tissue_context,preparation)]
 d<-merge(d,ctxlabels,by=c('patient_id','tissue_context','preparation'),all.x=TRUE,sort=FALSE)
 d[,row_label:=paste(sub('.*__','',patient_id),ifelse(tissue_context=='Healthy_breast','Normal',patient_clinical_labels),sep=' | ')]
 oo<-unique(d[order(tissue_context,patient_clinical_labels,patient_id),row_label]);stopifnot(length(oo)==uniqueN(d$patient_id));d[,row_label:=factor(row_label,levels=rev(oo))];d[,cell_type:=factor(cell_type,levels=cts)]
 p<-ggplot(d,aes(cell_type,row_label))+geom_point(aes(size=pct_detected,color=log2_CPM_plus1),alpha=.95)+geom_point(data=d[pct_detected==0],shape=1,size=1,color='grey50')+geom_point(data=d[n_cells<20],shape=4,size=2.8,color='black',stroke=.5)+scale_color_viridis_c(option='C',name='log2(CPM + 1)',limits=c(0,color_max))+scale_size_area(max_size=6,limits=c(0,100),breaks=c(5,20,50,100),name='Detected (%)')+labs(title=paste('MFSD10 by patient and broad cell type:',prep),subtitle='Each row is one patient within tissue / preparation; no pooled-cell inference',x=NULL,y='Patient | published clinical label',caption='Cross = fewer than 20 cells. Small grey ring = sampled but zero detected. Blank = no sampled cells.\nRNA counts summed per patient / cell type; CPM uses all RNA UMIs. Cross-cell-type differences are descriptive.')+theme(axis.text.x=element_text(angle=35,hjust=1),axis.text.y=element_text(size=8))
 savefig(p,paste0('01_',prep,'_patient_celltype'),11,ifelse(prep=='Total',12,6))
}
# Retain specimen-level clinical labels for the multi-lesion patient 0114.
e[,row_label:=paste(sub('.*__','',patient_id),clinical_label,ifelse(preparation=='Total','Total','Epi'),sep=' | ')]
oo<-unique(e[order(clinical_label,patient_id,preparation),row_label]);e[,row_label:=factor(row_label,levels=rev(oo))]
el<-rbindlist(list(e[,.(row_label,patient_id,epithelial_status,n_cells,metric='Detection (%)',value=pct_detected)],e[,.(row_label,patient_id,epithelial_status,n_cells,metric='log2(CPM + 1)',value=log2_CPM_plus1)]))
p<-ggplot(el,aes(value,row_label,color=epithelial_status))+geom_point(aes(shape=n_cells>=20),size=2,position=position_dodge(width=.65))+scale_shape_manual(values=c('TRUE'=16,'FALSE'=1),name='At least 20 cells')+scale_color_manual(values=pal,labels=state_labels,name=NULL)+facet_grid(.~metric,scales='free_x')+labs(title='Epithelial MFSD10: every patient, clinical label and CNV status',subtitle='Healthy context is separate from CNV status; no verified tumor nonmalignant class is available',x=NULL,y='Patient | clinical label | preparation',caption='Patient 0114 appears under ER+ and TNBC because these are distinct lesions of the same patient.\nCNV labels cover only the discovery subset; uncertain and unscored cells are retained. Open symbols: fewer than 20 cells.')+theme(axis.text.y=element_text(size=7),legend.text=element_text(size=8),legend.box='vertical')
savefig(p,'02_epithelial_clinical_patient_states',12,14)
mal<-copy(cells[epithelial_status=='Malignant_CNV_coherent'])
patientstats<-mal[,.(n_cells=.N,n_detected=sum(MFSD10_detected),pct_detected=100*mean(MFSD10_detected),MFSD10_UMI=sum(MFSD10_counts),CPM=1e6*sum(MFSD10_counts)/sum(nCount_RNA),mean_log1p_CP10K=mean(MFSD10_log1p_CP10K),median_log1p_CP10K=median(MFSD10_log1p_CP10K),q90_log1p_CP10K=quantile(MFSD10_log1p_CP10K,.9),q99_log1p_CP10K=quantile(MFSD10_log1p_CP10K,.99),maximum_counts=max(MFSD10_counts),n_libraries=uniqueN(gsm),clinical_labels=paste(sort(unique(clinical_label)),collapse=';')),by=patient_id]
patientstats[,panel:=sprintf('%s | %s\nn=%d; detected=%.1f%%',sub('.*__','',patient_id),clinical_labels,n_cells,pct_detected)]
fwrite(patientstats,file.path(out,'results/Malignant_patient_continuous_summary.tsv'),sep='\t')
lev<-patientstats[order(clinical_labels,patient_id),panel];mal[,panel:=factor(patientstats$panel[match(patient_id,patientstats$patient_id)],levels=lev)]
ed<-mal[, {xx<-sort(unique(c(0,MFSD10_log1p_CP10K,max(mal$MFSD10_log1p_CP10K))));.(expression=xx,cumulative=ecdf(MFSD10_log1p_CP10K)(xx))},by=panel]
p<-ggplot(ed,aes(expression,cumulative))+geom_step(linewidth=.5,color='#8F3052')+geom_point(data=ed[expression==0],size=1,color='#8F3052')+facet_wrap(~panel,ncol=4)+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+labs(title='MFSD10 continuous expression within malignant epithelium',subtitle='One empirical cumulative distribution per patient; all cells and all zeros retained',x='log1p(MFSD10 RNA UMI / total RNA UMI x 10,000)',y='Fraction of patient cells at or below expression',caption='The value at zero is the non-detection fraction. No MFSD10 high / low split or cell-level significance tests.\nPatient keys combine available malignant cells; small patient groups are displayed without population inference.')
savefig(p,'03_malignant_patient_ECDF',12,12)
rho<-function(x,y){if(length(x)<20 || length(unique(x))<2 || length(unique(y))<2)return(NA_real_);suppressWarnings(cor(x,y,method='spearman'))}
cp<-list();curves<-list()
for(scope in c('All_malignant','Singlets_only')){
 z<-if(scope=='All_malignant')copy(mal) else copy(mal[scDblFinder_class=='singlet'])
 for(v in c('nCount_RNA','nFeature_RNA','nCount_without_MFSD10','nFeature_without_MFSD10')){
  tt<-z[,.(n_cells=.N,n_detected=sum(MFSD10_detected),rho_raw=rho(MFSD10_counts,get(v)),rho_log_normalized=rho(MFSD10_log1p_CP10K,get(v)),rho_detection=rho(as.integer(MFSD10_detected),get(v)),rho_positive_only=rho(MFSD10_log1p_CP10K[MFSD10_detected],get(v)[MFSD10_detected]),adequate_for_description=.N>=20 & sum(MFSD10_detected)>=5),by=patient_id];tt[,`:=`(scope=scope,covariate=v)];cp[[paste(scope,v)]]<-tt
  if(scope=='All_malignant' && v%in%c('nCount_RNA','nFeature_RNA')){
   z[,depth_bin:=pmin(5L,ceiling(frank(get(v),ties.method='average')/.N*5)),by=patient_id]
   cr<-z[,.(n_cells=.N,median_depth=as.numeric(median(get(v))),detection=100*mean(MFSD10_detected),mean_log_normalized=mean(MFSD10_log1p_CP10K)),by=.(patient_id,depth_bin)];cr[,covariate:=v];curves[[v]]<-cr
  }
 }
}
corr<-rbindlist(cp);curve<-rbindlist(curves)
fwrite(corr,file.path(out,'results/Malignant_patient_QC_correlations.tsv'),sep='\t');fwrite(curve,file.path(out,'results/Malignant_patient_depth_quintiles.tsv'),sep='\t')
cl<-rbindlist(list(curve[,.(patient_id,covariate,depth_bin,median_depth,metric='Detection (%)',value=detection,n_cells)],curve[,.(patient_id,covariate,depth_bin,median_depth,metric='Mean log1p(CP10K)',value=mean_log_normalized,n_cells)]))
p1<-ggplot(cl,aes(log10(median_depth),value,group=patient_id))+geom_line(color='#5680A0',alpha=.35,linewidth=.4)+geom_point(aes(shape=n_cells>=20),color='#5680A0',alpha=.4,size=1)+facet_grid(metric~covariate,scales='free')+scale_shape_manual(values=c('TRUE'=16,'FALSE'=1),name='Bin has >=20 cells')+labs(title='Technical depth relationship inside each malignant patient',subtitle='Each line is one patient; five within-patient depth bins, without expression grouping',x='log10(median RNA depth / detected genes in bin)',y=NULL,caption='These are descriptive within-patient curves, not independent patient replicates or adjusted causal effects.')
savefig(p1,'04_malignant_depth_relationship',10,7)
# Patient-level rho display provides patient identity beside technical curves.
qq<-corr[scope=='All_malignant' & covariate%in%c('nCount_RNA','nFeature_RNA')]
qq<-melt(qq,id.vars=c('patient_id','covariate','adequate_for_description'),measure.vars=c('rho_raw','rho_log_normalized','rho_detection'),variable.name='measure',value.name='rho')
qq[,patient:=factor(sub('.*__','',patient_id),levels=rev(sub('.*__','',patientstats[order(clinical_labels,patient_id),patient_id])))]
p<-ggplot(qq,aes(rho,patient,color=covariate))+geom_vline(xintercept=0,color='grey75')+geom_point(aes(shape=adequate_for_description),position=position_dodge(width=.45),size=2,na.rm=TRUE)+facet_grid(.~measure)+scale_color_manual(values=c('nCount_RNA'='#257BA6','nFeature_RNA'='#CF7A23'))+scale_shape_manual(values=c('TRUE'=16,'FALSE'=1))+coord_cartesian(xlim=c(-1,1))+labs(title='Within-patient MFSD10 association with RNA QC metrics',x='Spearman rho (descriptive; no cell-level P values)',y='Patient',color=NULL,shape='>=20 cells and >=5 detected',caption='Raw, normalized and detection associations shown separately. Blank = insufficient cells or constant variable. Part-whole and singlet-only checks are in the table.')
savefig(p,'04b_malignant_patient_QC_rho',11,9)
# Unannotated samples remain sample-level descriptive mixtures.
un<-b[cell_type=='Not_annotated'];fwrite(un,file.path(out,'results/Unannotated_LN_preneoplastic_mixture_only.tsv'),sep='\t')
# Completeness and missingness are explicit, with '-' treated as missing in supplements.
fields<-c('patient_id','clinical_specimen','clinical_age','gender','clinical_label','tissue_context','preparation','published_tumor_type','histology','ER_text','grade_text','size_mm_text','BRCA1_variant','HER2_text')
miss<-rbindlist(lapply(fields,function(f){a<-copy(m);a[,available_field:=!is.na(get(f)) & trimws(as.character(get(f)))!='' & trimws(as.character(get(f)))!='-'];a[,.(field=f,n_libraries=.N,n_available=sum(available_field)),by=tissue_context]}))
fwrite(miss,file.path(out,'metadata/Clinical_field_availability.tsv'),sep='\t')
checks<-data.table(check=c('All_QC_cells_340979','Unique_cell_ids','All_69_libraries','Annotated_58_libraries','Main_Total_cells_236309','Malignant_CNV_9486','Uncertain_CNV_19627','CNV_scored_only_Epithelial','No_CNV_labels_in_LN_or_preneoplastic','No_inferred_tumor_nonmalignant','Broad_patient_UMI_conservation','Broad_patient_cell_conservation','Epithelial_patient_UMI_conservation','Epithelial_patient_cell_conservation','Repeated_0029_regions_collapsed','Mixed_0114_one_broad_patient','MFSD10_id_exact'),pass=c(nrow(cells)==340979,!anyDuplicated(cells$cell_id),uniqueN(cells$gsm)==69,uniqueN(cells[cell_type!='Not_annotated',gsm])==58,cells[preparation=='Total' & tissue_context%in%c('Primary_tumor','Healthy_breast'),.N]==236309,nrow(mal)==9486,cells[epithelial_status=='Uncertain_CNV_scored',.N]==19627,all(cells[epithelial_status%in%c('Malignant_CNV_coherent','Uncertain_CNV_scored'),cell_type]=='Epithelial'),cells[tissue_context%in%c('Lymph_node','Preneoplastic') & epithelial_status%in%names(pal),.N]==0,!any(grepl('nonmalignant',cells$epithelial_status)),sum(b$MFSD10_UMI)==sum(cells$MFSD10_counts),sum(b$n_cells)==nrow(cells),sum(e$MFSD10_UMI)==cells[cell_type=='Epithelial',sum(MFSD10_counts)],sum(e$n_cells)==cells[cell_type=='Epithelial',.N],b[patient_id=='Cancer__0029' & cell_type=='Epithelial',n_libraries]==2,b[patient_id=='Cancer__0114' & cell_type=='Epithelial',.N]==1,fread(file.path(out,'metadata/MFSD10_gene_mapping.tsv'))$canonical_id=='ENSG00000109736'))
fwrite(checks,file.path(out,'audit/Analysis_validation.tsv'),sep='\t');stopifnot(all(checks$pass))
capture.output(sessionInfo(),file=file.path(out,'logs/02_sessionInfo.txt'));cat('FIGURES_AND_CHECKS_COMPLETE',format(Sys.time()),'\n');sink()
