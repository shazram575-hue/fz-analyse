options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
meta <- fread(file.path(out, "metadata", "GSE161529_author_QC_and_analysis_membership.tsv"))
meta <- meta[include_main == TRUE]
statuspath <- file.path(out, "audit", "Seurat_main_library_run_status.tsv")
statuses <- list()
for (i in seq_len(nrow(meta))) {
  gsm <- meta$gsm[i]
  final <- file.path(out, "cache", "seurat_pilot", paste0(gsm, "_pilot_seurat.rds"))
  if (file.exists(final)) {
    cat("EXISTS", gsm, "\n")
    statuses[[i]] <- data.table(gsm = gsm, state = "complete", completed = as.character(file.info(final)$mtime))
  } else {
    cat("START LIBRARY", i, "/", nrow(meta), gsm, format(Sys.time()), "\n")
    options(workflow3.library.args = c(root, gsm))
    stage_env <- new.env(parent = globalenv())
    tryCatch(source(file.path(out, "scripts", "07_seurat_library_pilot.R"), local = stage_env, encoding = "UTF-8"),
             error = function(e) { while (sink.number()) sink(); stop("Library ", gsm, ": ", conditionMessage(e)) })
    rm(stage_env); gc(verbose = FALSE)
    stopifnot(file.exists(final))
    statuses[[i]] <- data.table(gsm = gsm, state = "complete", completed = format(Sys.time()))
  }
  fwrite(rbindlist(statuses), statuspath, sep = "\t")
}
cat("ALL MAIN LIBRARIES COMPLETE", format(Sys.time()), "\n")
