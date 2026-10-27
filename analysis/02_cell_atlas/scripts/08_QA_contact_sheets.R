options(stringsAsFactors=FALSE)
.libPaths(c("D:/R4.5.2/R-4.5.2/library",.libPaths()))
library(png);library(grid)
root<-normalizePath(commandArgs(TRUE)[1],winslash="/")
out<-file.path(root,"outputs","单细胞基础补全_2026-10-23")
files<-list.files(file.path(out,"figures"),pattern="^Foundation_.*\\.png$",full.names=TRUE)
stopifnot(length(files)>0)
for(page in seq_len(ceiling(length(files)/6))) {
 path<-file.path(out,"audit",paste0("Figure_contact_sheet_",page,".png"))
 png(path,width=2000,height=2100,res=150,bg="white")
 grid.newpage();pushViewport(viewport(layout=grid.layout(3,2)))
 ids<-seq.int((page-1)*6+1,min(page*6,length(files)))
 for(j in seq_along(ids)) {
  pushViewport(viewport(layout.pos.row=ceiling(j/2),layout.pos.col=(j-1)%%2+1))
  grid.text(basename(files[ids[j]]),y=.98,gp=gpar(fontsize=8))
  img<-png::readPNG(files[ids[j]])
  ratio<-dim(img)[2]/dim(img)[1];target_ratio<-1000/650
  ww<-if(ratio>target_ratio).96 else .96*ratio/target_ratio
  hh<-if(ratio>target_ratio).9*target_ratio/ratio else .9
  grid.raster(img,x=.5,y=.47,width=ww,height=hh,interpolate=TRUE)
  popViewport()
 }
 popViewport();dev.off()
}
cat("Contact sheets for",length(files),"PNGs\n")
