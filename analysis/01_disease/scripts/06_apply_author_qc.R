options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
root <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
sink(file.path(out, "logs", "06_apply_author_qc.log"), split = TRUE)
cat("START", format(Sys.time()), "\n")
meta <- fread(file.path(out, "metadata", "GSE161529_sample_metadata_unfiltered.tsv"), colClasses = "character")
author <- fread(file.path(out, "metadata", "author_reference", "Tables__SampleStats.txt"))
stopifnot(!anyDuplicated(author$GEO_ID), setequal(meta$gsm, author$GEO_ID))
meta <- merge(meta, author, by.x = "gsm", by.y = "GEO_ID", all.x = TRUE, sort = FALSE)
meta[, analysis_role := fifelse(cancer.type == "Normal" & cell.population == "Total", "main_healthy_total",
                  fifelse(cancer.type == "Normal" & cell.population == "Epithelial", "normal_epi_sensitivity",
                  fifelse(cancer.type == "BRCA1 pre-neoplastic", "preneoplastic_separate",
                  fifelse(cell.population == "Involved lymph node", "lymph_node_separate", "main_primary_total"))))]
meta[, include_main := analysis_role %in% c("main_healthy_total", "main_primary_total")]
stopifnot(sum(meta$analysis_role == "main_healthy_total") == 13L, sum(meta$analysis_role == "main_primary_total") == 34L)
meta[, patient_id_SOFT_unresolved := patient]
meta[, patient_identity_status := "SOFT identifiers require reconciliation; do not infer clinical pairing"]
fwrite(meta, file.path(out, "metadata", "GSE161529_author_QC_and_analysis_membership.tsv"), sep = "\t")
summaries <- list()
allQc <- list()
qcd <- file.path(out, "audit", "sc_unfiltered_cell_qc")
for (i in seq_len(nrow(meta))) {
  gsm <- meta$gsm[i]
  qfiles <- list.files(qcd, pattern = paste0("^", gsm, "_.*_raw_qc\\.tsv\\.gz$"), full.names = TRUE)
  stopifnot(length(qfiles) == 1L)
  qc <- as.data.table(read.delim(gzfile(qfiles), check.names = FALSE))
  stopifnot(uniqueN(qc$barcode) == nrow(qc))
  if (nrow(qc) != meta$CellNum[i]) cat("RAW CELL NUMBER DISAGREEMENT", gsm, "current", nrow(qc), "author table", meta$CellNum[i], "\n")
  qc[, `:=`(pass_mito = percent_mt/100 < meta$Mito[i],
             pass_gene_lower = nFeature_RNA > meta$GeneLower[i],
             pass_gene_upper = nFeature_RNA < meta$GeneUpper[i],
             pass_umi_upper = nCount_RNA < meta$LibSize[i])]
  qc[, pass_author_qc := pass_mito & pass_gene_lower & pass_gene_upper & pass_umi_upper]
  qc[, `:=`(cell_id = paste(gsm, barcode, sep = "__"), author_sample = meta$SampleName[i],
             analysis_role = meta$analysis_role[i], include_main = meta$include_main[i])]
  kept <- sum(qc$pass_author_qc)
  summaries[[i]] <- data.table(gsm = gsm, author_sample = meta$SampleName[i], analysis_role = meta$analysis_role[i],
                     include_main = meta$include_main[i], raw_cells = nrow(qc), author_reported_raw = meta$CellNum[i],
                     raw_difference_from_author = nrow(qc) - meta$CellNum[i],
                     author_reported_kept = meta$CellNumAfter[i], reproduced_kept = kept,
                     kept_difference_from_author = kept - meta$CellNumAfter[i],
                     fail_mito = sum(!qc$pass_mito), fail_gene_lower = sum(!qc$pass_gene_lower),
                     fail_gene_upper = sum(!qc$pass_gene_upper), fail_umi_upper = sum(!qc$pass_umi_upper))
  allQc[[i]] <- qc
}
sm <- rbindlist(summaries)
all <- rbindlist(allQc)
stopifnot(!anyDuplicated(all$cell_id))
fwrite(sm, file.path(out, "audit", "GSE161529_author_QC_reproduction.tsv"), sep = "\t")
gz <- gzfile(file.path(out, "cache", "GSE161529_all_cell_QC_membership.tsv.gz"), "wt")
write.table(all, gz, sep = "\t", row.names = FALSE, quote = FALSE); close(gz)
saveRDS(split(all[pass_author_qc == TRUE, barcode], all[pass_author_qc == TRUE, gsm]),
        file.path(out, "cache", "GSE161529_author_QC_kept_barcodes_by_gsm.rds"))
role <- sm[, .(profiles = .N, raw_cells = sum(raw_cells), qc_cells = sum(reproduced_kept)), by = analysis_role]
fwrite(role, file.path(out, "results", "SingleCell_author_QC_by_analysis_role.tsv"), sep = "\t")
print(role)
cat("AUTHOR QC DIFFERENCES", sum(sm$kept_difference_from_author != 0), "MAX ABS DIFF", max(abs(sm$kept_difference_from_author)), "\n")
pdf(file.path(out, "figures", "SingleCell_author_QC_reproduction.pdf"), width = 15, height = 7)
par(mar = c(12, 5, 3, 1))
barplot(rbind(sm$reproduced_kept, sm$raw_cells - sm$reproduced_kept), names.arg = sm$author_sample,
        las = 2, cex.names = .55, col = c("#4E79A7", "#BAB0AC"), border = NA,
        main = "GSE161529: published sample-specific QC reproduced", ylab = "Cells")
legend("topright", legend = c("Retained", "Filtered"), fill = c("#4E79A7", "#BAB0AC"), bty = "n")
dev.off()
cat("COMPLETE", format(Sys.time()), "\n")
capture.output(sessionInfo(), file = file.path(out, "logs", "06_R_sessionInfo.txt"))
sink()
