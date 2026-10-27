# Read-only extraction of original published clinical/QC supplementary tables.
options(stringsAsFactors=FALSE,warn=1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(data.table)
library(xml2)
library(readxl)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
ref<-file.path(root,"outputs","第一阶段_2026-10-23","metadata","author_reference")
zip<-file.path(ref,"Pal2021_supplementaryFiles.zip")
contents<-unzip(zip,list=TRUE)
selected<-contents$Name[grepl("\\.(docx|xlsx)$",contents$Name)]
stopifnot(!any(grepl("(^|/)\\.\\.(/|$)|^/|^[A-Za-z]:",selected)))
unzip(zip,files=selected,exdir=ref)
for(file in selected) {
  path<-file.path(ref,file)
  if(grepl("\\.docx$",file)) {
    con<-unz(path,"word/document.xml")
    doc<-read_xml(con)
    ns<-xml_ns(doc)
    paragraphs<-xml_find_all(doc,".//w:p",ns)
    text<-vapply(paragraphs,function(p)paste(xml_text(xml_find_all(p,".//w:t",ns)),collapse=""),character(1))
    writeLines(text,file.path(ref,paste0(file,"_paragraphs.txt")),useBytes=TRUE)
    tabs<-xml_find_all(doc,".//w:tbl",ns)
    for(i in seq_along(tabs)) {
      rows<-lapply(xml_find_all(tabs[[i]],"./w:tr",ns),function(tr)
        vapply(xml_find_all(tr,"./w:tc",ns),function(tc)
          paste(xml_text(xml_find_all(tc,".//w:t",ns)),collapse=" "),character(1)))
      nc<-max(lengths(rows));mat<-matrix("",length(rows),nc)
      for(j in seq_along(rows))mat[j,seq_along(rows[[j]])]<-rows[[j]]
      fwrite(as.data.table(mat),file.path(ref,paste0(file,"_table",i,".tsv")),sep="\t",col.names=FALSE)
    }
    cat(file,"TABLES",length(tabs),"\n")
  } else {
    for(sheet in excel_sheets(path)) {
      tab<-as.data.table(read_excel(path,sheet=sheet,col_names=FALSE,.name_repair="minimal"))
      fwrite(tab,file.path(ref,paste0(file,"_",gsub("[^[:alnum:]_]","_",sheet),".tsv")),sep="\t",col.names=FALSE)
      cat(file,sheet,"ROWS",nrow(tab),"COLS",ncol(tab),"\n")
    }
  }
}
