# Compute per-cell QC on all libraries; deliberately applies no cell/sample exclusion.
options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
setDTthreads(2L)
root <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
sink(file.path(out, "logs", "03_sc_unfiltered_qc.log"), split = TRUE)
cat("START", format(Sys.time()), "\n")
features <- fread(file.path(out, "metadata", "GSE161529_features.tsv"))
pairs <- fread(file.path(out, "audit", "GSE161529_tar_file_pairs.tsv"))
unpacked <- file.path(out, "cache", "GSE161529_unpacked")
qcdir <- file.path(out, "audit", "sc_unfiltered_cell_qc")
dir.create(unpacked, recursive = TRUE, showWarnings = FALSE)
dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)
needed <- c(pairs$matrix_member, pairs$barcode_member)
missing <- needed[!file.exists(file.path(unpacked, needed))]
if (length(missing)) {
  cat("EXTRACT", length(missing), "GZIP MEMBERS\n")
  utils::untar(file.path(root, "原始数据", "单细胞", "GSE161529_RAW.tar"), files = missing, exdir = unpacked)
}
mt <- startsWith(features$gene_symbol, "MT-")
stopifnot(sum(mt) > 0L)
summaries <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  key <- pairs$profile_key[i]
  qcp <- file.path(qcdir, paste0(key, "_raw_qc.tsv.gz"))
  sp <- file.path(qcdir, paste0(key, "_summary.tsv"))
  if (file.exists(qcp) && file.exists(sp)) {
    summaries[[i]] <- fread(sp)
    cat("RESUME", i, key, "\n"); next
  }
  cat("READ", i, "/", nrow(pairs), key, format(Sys.time()), "\n")
  flush.console()
  con <- gzfile(file.path(unpacked, pairs$matrix_member[i]), "rt")
  x <- readMM(con); close(con)
  stopifnot(nrow(x) == nrow(features), all(is.finite(x@x)), all(x@x >= 0), all(x@x == floor(x@x)))
  x <- as(x, "CsparseMatrix")
  bc <- readLines(gzfile(file.path(unpacked, pairs$barcode_member[i])), warn = FALSE)
  stopifnot(length(bc) == ncol(x), !anyDuplicated(bc))
  # QC needs only column summaries; retain no Seurat object or dense matrix.
  count <- Matrix::colSums(x)
  nf <- diff(x@p)
  stopifnot(!any(x@x == 0))
  mtcount <- Matrix::colSums(x[mt, , drop = FALSE])
  qc <- data.table(gsm = pairs$gsm[i], profile_key = key, barcode = bc,
                   nCount_RNA = count, nFeature_RNA = nf,
                   percent_mt = ifelse(count > 0, 100 * mtcount/count, NA_real_))
  con <- gzfile(qcp, "wt")
  write.table(qc, con, sep = "\t", row.names = FALSE, quote = FALSE); close(con)
  sm <- data.table(gsm = pairs$gsm[i], profile_key = key, n_genes = nrow(x), n_cells = ncol(x),
                   nnz = length(x@x), zero_library_cells = sum(count == 0),
                   count_q01 = unname(quantile(count, .01)), count_median = median(count), count_q99 = unname(quantile(count, .99)),
                   feature_q01 = unname(quantile(nf, .01)), feature_median = median(nf), feature_q99 = unname(quantile(nf, .99)),
                   mt_q01 = unname(quantile(qc$percent_mt, .01, na.rm = TRUE)), mt_median = median(qc$percent_mt, na.rm = TRUE),
                   mt_q99 = unname(quantile(qc$percent_mt, .99, na.rm = TRUE)))
  fwrite(sm, sp, sep = "\t")
  summaries[[i]] <- sm
  rm(x, count, nf, mtcount, bc, qc); gc(verbose = FALSE)
  fwrite(rbindlist(summaries[seq_len(i)]), file.path(out, "audit", "GSE161529_unfiltered_qc_summary.tsv"), sep = "\t")
}
total <- rbindlist(summaries)
fwrite(total, file.path(out, "audit", "GSE161529_unfiltered_qc_summary.tsv"), sep = "\t")
cat("COMPLETE", format(Sys.time()), "PROFILES", nrow(total), "RAW CELLS", sum(total$n_cells), "\n")
capture.output(sessionInfo(), file = file.path(out, "logs", "03_R_sessionInfo.txt"))
sink()
