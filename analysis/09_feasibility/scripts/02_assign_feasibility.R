options(stringsAsFactors=FALSE)
.libPaths(c('D:/R4.5.2/R-4.5.2/library',.libPaths()));library(data.table)
root<-normalizePath(commandArgs(TRUE)[1],winslash='/');out<-file.path(root,'outputs/恶性上皮_表达程序可行性_2026-10-26')
p<-fread(file.path(out,'results/Patient_counts_complexity_audit.tsv'))
# Proposed study-specific allocation; not a published universal minimum or an executed cell filter.
p[,allocation:=fcase(n_singlet<50,'E_仅描述',n_singlet<200 & singlet_UMI_median>=1000 & singlet_genes_median>=1000,'C_仅投射候选',singlet_UMI_median<1000 | singlet_genes_median<1000,'D_低复杂度描述优先',n_singlet<300,'B_小K试跑候选',default='A_发现试跑候选')]
p[,allocation_reason:=fcase(allocation=='E_仅描述','独立细胞数太少，不逐患者发现程序',allocation=='C_仅投射候选','50–199 singlets；不加入当前主发现集，投射仍需覆盖/重构检查',allocation=='D_低复杂度描述优先','singlet UMI中位<1000或基因中位<1000；若投射需验证覆盖与深度偏差',allocation=='B_小K试跑候选','200–299 singlets且复杂度达本轮提议线；仅低秩试跑，稳定性通过后再评估',default='>=300 singlets且复杂度达本轮提议线；仅候选，未验证程序可辨识性')]
setorder(p,allocation,-n_singlet)
fwrite(p,file.path(out,'results/患者可用性表.tsv'),sep='\t')
q<-p[,.(patients=.N,confirmed_cells=sum(n_confirmed),singlet_flags=sum(n_singlet),doublet_flags=sum(n_doublet),patient_ids=paste(sub('Cancer__','',patient_id),collapse=';')),by=allocation]
fwrite(q,file.path(out,'results/Allocation_summary.tsv'),sep='\t');print(q)
cat('TOTAL CONFIRMED',sum(p$n_confirmed),'SINGLETS',sum(p$n_singlet),'DOUBLETS',sum(p$n_doublet),'\n')
fmt<-function(z)ifelse(z==round(z),format(z,scientific=FALSE,trim=TRUE),format(round(z,1),scientific=FALSE,trim=TRUE))
short<-c('A_发现试跑候选'='A发现候选','B_小K试跑候选'='B小K候选','C_仅投射候选'='C仅投射','D_低复杂度描述优先'='D低复杂度','E_仅描述'='E仅描述')
samples<-vapply(seq_len(nrow(p)),function(i){g<-sub('.*__','',p$patient_id[i]);if(g=='0029')return('ER-0029-7C:19；9C:126');if(g=='0114')return('TN-0114-T2:1；ER-T3:0');paste0(sub('.*\\[([^]]+)\\].*','\\1',p$malignant_sample_composition[i]),':',p$n_confirmed[i])},character(1))
lines<-c('# 恶性上皮跨患者表达程序：患者可用性表','',
'UMI/基因中位数按全部确认恶性细胞统计；分类依据singlet子集数量及复杂度，均未执行新过滤。完整TSV另含singlet中位数、四分位数、双细胞数、GSM及所有原发病灶。',
'',
'| 患者 | 确认恶性/其中singlet | UMI中位 | 基因中位 | 临床标签（恶性来源） | 样本:恶性数 | 建议用途 |',
'|---|---:|---:|---:|---|---|---|',
sprintf('| %s | %d/%d | %s | %s | %s | %s | %s |',sub('.*__','',p$patient_id),p$n_confirmed,p$n_singlet,fmt(p$UMI_median),fmt(p$genes_median),p$malignant_clinical_labels,samples,short[p$allocation]),'',
'同一患者多文库不增加患者数；0029两区域合并；0114另有ER+病灶但无确认恶性细胞。以上仅原发Total；LN/癌前因无粗注释不纳入。临床标签不推定PAM50。')
writeLines(lines,file.path(out,'患者可用性表.md'),useBytes=TRUE)
