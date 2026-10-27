options(stringsAsFactors=FALSE)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
ref<-file.path(out,"metadata","author_reference")
m<-fread(file.path(out,"metadata","GSE161529_author_QC_and_analysis_membership.tsv"),colClasses="character")
clean<-function(x)gsub("[[:space:]]+","",x)
normal<-fread(file.path(ref,"EMBJ-40-e107333-s004.docx_table1.tsv"),skip=1,header=FALSE,fill=TRUE,sep="\t")
er<-fread(file.path(ref,"EMBJ-40-e107333-s001.docx_table1.tsv"),skip=1,header=FALSE,fill=TRUE,sep="\t")
tn<-fread(file.path(ref,"EMBJ-40-e107333-s005.docx_table1.tsv"),skip=1,header=FALSE,fill=TRUE,sep="\t")
ages<-rbindlist(list(data.table(specimen=clean(normal[[1]]),age=as.numeric(clean(normal[[2]])),source="Pal2021_EV1"),
                    data.table(specimen=clean(er[[1]]),age=as.numeric(clean(er[[2]])),source="Pal2021_EV3"),
                    data.table(specimen=clean(tn[[1]]),age=as.numeric(clean(tn[[2]])),source="Pal2021_EV2")))
ages<-ages[!is.na(age) & grepl("^(N-|B1-|ER-|TN-|HER2-|mER-)",specimen)]
stopifnot(!anyDuplicated(ages$specimen))
m[,clinical_specimen:=sub("-(total|epi)$","",SampleName)]
m[,clinical_age:=ages$age[match(clean(clinical_specimen),ages$specimen)]]
m[,clinical_table:=ages$source[match(clean(clinical_specimen),ages$specimen)]]
if(anyNA(m$clinical_age)) { print(m[is.na(clinical_age),.(gsm,SampleName,clinical_specimen)]); stop("Unmapped clinical specimens") }
m[,donor_namespace:=fifelse(cancer.type=="Normal","Healthy",fifelse(cancer.type=="BRCA1 pre-neoplastic","BRCA1_preneo","Cancer"))]
m[,patient_id:=paste(donor_namespace,patient,sep="__")]
m[,patient_mapping_basis:="Published SOFT/EV4 patient field within cohort; EV1-EV3 clinical age cross-check"]
m[patient=="0064",patient_mapping_basis:="Healthy N-0064 age44 versus cancer ER-0064 age65: separate donors; T/LN share cancer donor"]
m[patient=="0114",patient_mapping_basis:="Published patient0114, lesions T2/T3, both age84: aggregate by patient; split sensitivity"]
m[patient=="0029",patient_mapping_basis:="Published patient0029, regions7C/9C, age59: aggregate by patient; split sensitivity"]
stopifnot(all(m[,uniqueN(clinical_age),by=patient_id]$V1==1L))
fwrite(m,file.path(out,"metadata","GSE161529_reconciled_patient_metadata.tsv"),sep="\t")
summary<-m[,.(n_profiles=.N,age=unique(clinical_age),gender=paste(unique(gender),collapse=";"),specimens=paste(SampleName,collapse=";"),GSMs=paste(gsm,collapse=";")),by=patient_id]
fwrite(summary,file.path(out,"audit","GSE161529_patient_identity_reconciliation.tsv"),sep="\t")
print(m[,.(n_profiles=.N,n_patient_keys=uniqueN(patient_id)),by=donor_namespace])
cat("MAIN PROFILES",sum(m$include_main=="TRUE"),"MAIN PATIENT KEYS",uniqueN(m[include_main=="TRUE",patient_id]),"\n")
