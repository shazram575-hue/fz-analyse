options(stringsAsFactors=FALSE,warn=2)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
out<-file.path(root,"outputs","第一阶段_2026-10-23")
rd<-file.path(out,"metadata","genome_reference")
path<-file.path(rd,"gencode.v36.annotation.gtf.gz")
md5lines<-readLines(file.path(rd,"GENCODE_release36_MD5SUMS.txt"),warn=FALSE)
expectedLine<-md5lines[grepl("gencode.v36.annotation.gtf.gz",md5lines,fixed=TRUE)]
stopifnot(length(expectedLine)==1L)
expected<-strsplit(trimws(expectedLine),"[[:space:]]+")[[1]][1]
actual<-unname(tools::md5sum(path))
stopifnot(identical(tolower(expected),tolower(actual)))
fwrite(data.table(source="GENCODE release36 official MD5SUMS",file=path,expected_md5=expected,actual_md5=actual,match=TRUE),
       file.path(rd,"GENCODE_reference_md5_check.tsv"),sep="\t")
con<-gzfile(path,"rt")
blocks<-list();i<-0L
repeat {
  x<-readLines(con,n=100000L,warn=FALSE)
  if(!length(x))break
  x<-x[grepl("\tgene\t",x,fixed=TRUE)]
  if(!length(x))next
  i<-i+1L
  z<-as.data.table(tstrsplit(x,"\t",fixed=TRUE))
  att<-z[[9]]
  attr<-function(key)sub(paste0('.*',key,' "([^"]+)".*'),'\\1',att)
  blocks[[i]]<-data.table(gene_id=attr("gene_id"),gene_name=attr("gene_name"),gene_type=attr("gene_type"),
          chromosome=z[[1]],start=as.integer(z[[4]]),end=as.integer(z[[5]]),strand=z[[7]])
}
close(con)
genes<-rbindlist(blocks)
genes[,ensembl_stable:=sub("\\.[0-9]+(_PAR_Y)?$","\\1",gene_id)]
stopifnot(nrow(genes)>60000L,!anyDuplicated(genes$gene_id))
fwrite(genes,file.path(rd,"GENCODE_v36_gene_annotation.tsv"),sep="\t")
f<-fread(file.path(out,"metadata","GSE161529_features.tsv"))
auto<-genes[chromosome%in%paste0("chr",1:22)]
auto[,chr_number:=as.integer(sub("chr","",chromosome))]
setorder(auto,chr_number,start,end)
stopifnot(!anyDuplicated(auto$ensembl_stable))
go<-auto[ensembl_stable%in%f$ensembl_id,.(gene_id=ensembl_stable,chromosome,start,end)]
stopifnot(nrow(go)>28000L)
fwrite(go,file.path(rd,"infercnv_GENCODEv36_autosomal_gene_order.tsv"),sep="\t",col.names=FALSE)
map<-merge(f,genes[,.(ensembl_stable,gencode_gene_name=gene_name,gene_type,chromosome,start,end)],
              by.x="ensembl_id",by.y="ensembl_stable",all.x=TRUE,sort=FALSE,allow.cartesian=TRUE)
fwrite(map,file.path(out,"metadata","SingleCell_features_GENCODEv36_biotype_coordinates.tsv"),sep="\t",na="NA")
cat("GENCODE GENES",nrow(genes),"SC AUTOSOMAL COORDINATES",nrow(go),"SC MAPPED",sum(!is.na(map$gene_type)),"\n")
