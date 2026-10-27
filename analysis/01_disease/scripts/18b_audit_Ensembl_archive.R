options(stringsAsFactors=FALSE)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table);library(jsonlite)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
records<-fromJSON(file.path(out,"metadata","Ensembl_archive_unresolved_feature_ids.json"),simplifyVector=FALSE)
mp<-fread(file.path(out,"metadata","SingleCell_canonical_gene_identity_map.tsv"))
g<-fread(file.path(out,"metadata","genome_reference","GENCODE_v36_gene_annotation.tsv"))
audit<-list();nrescue<-0L
for(r in records) {
 reps<-r$possible_replacement
 ids<-if(length(reps))vapply(reps,function(x)if(is.null(x$stable_id))"" else as.character(x$stable_id),character(1)) else character()
 scores<-if(length(reps))vapply(reps,function(x)if(is.null(x$score))NA_real_ else as.numeric(x$score),numeric(1)) else numeric()
 canon<-sub("\\.[0-9]+(_PAR_Y)?$","\\1",ids)
 eligible<-which(canon%in%g$ensembl_stable & !is.na(scores)&scores==1)
 accepted<-NA_character_
 if(length(eligible)==1L && length(ids)==1L) {
   accepted<-canon[eligible]
   j<-which(mp$ensembl_id==r$id & mp$mapping_status=="unmapped_original_id")
   if(length(j)) {
     z<-match(accepted,g$ensembl_stable)
     mp$canonical_id[j]<-accepted;mp$canonical_symbol[j]<-g$gene_name[z];mp$gene_type[j]<-g$gene_type[z]
     mp$mapping_status[j]<-"Ensembl_archive_exact_replacement"
     nrescue<-nrescue+length(j)
   }
 }
 audit[[length(audit)+1L]]<-data.table(source_id=as.character(r$id),
     archive_latest=if(is.null(r$latest))NA_character_ else r$latest,
     archive_last_release=if(is.null(r$release))NA_character_ else as.character(r$release),
     is_current=if(is.null(r$is_current))NA_character_ else as.character(r$is_current),
     n_replacements=length(ids),replacement_ids=paste(ids,collapse=";"),
     replacement_scores=paste(scores,collapse=";"),accepted_GENCODE36_replacement=accepted)
}
fwrite(rbindlist(audit,fill=TRUE),file.path(out,"audit","Ensembl_archive_legacy_ID_audit.tsv"),sep="\t",na="NA")
fwrite(mp,file.path(out,"metadata","SingleCell_canonical_gene_identity_map.tsv"),sep="\t",na="NA")
canon<-mp[,.(gene_symbol=canonical_symbol[1],gene_type=gene_type[1],
    source_feature_ids=paste(ensembl_id,collapse=";"),n_source_features=.N,
    mapping_status=paste(sort(unique(mapping_status)),collapse=";")),by=canonical_id]
fwrite(canon,file.path(out,"metadata","SingleCell_canonical_gene_metadata.tsv"),sep="\t",na="NA")
fwrite(mp[,.(n_source_features=.N),by=mapping_status],file.path(out,"audit","Gene_identity_mapping_summary.tsv"),sep="\t")
writeLines(c("Unresolved legacy IDs were checked with the official Ensembl archive API.",
"Only unequivocal score1 single replacements present in the TCGA-matched GENCODE36 model were accepted.",
"An archive peptide:null field is not evidence that a gene is noncoding and was not used as such.",
"Remaining legacy features stay in the RNA test universe; they are not relabeled as lncRNA/pseudogene.",
"Pool A requires an annotated protein-coding gene with a measured Bulk contrast; no gene identity is invented for unmatched records."),
 file.path(out,"audit","Gene_identity_legacy_ID_handling.txt"))
cat("ARCHIVE QUERIED",length(records),"ADDITIONAL EXACT RESCUES",nrescue,"CANONICAL ROWS",nrow(canon),"\n")
