options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
root <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
sink(file.path(out, "logs", "04_build_bulk_counts.log"), split = TRUE)
cat("START", format(Sys.time()), "\n")
meta <- fread(file.path(out, "metadata", "TCGA_verified_file_sample_metadata.tsv"))
gene <- fread(file.path(out, "metadata", "TCGA_first_file_GENCODEv36_gene_annotation.tsv"))
stopifnot(all(meta$md5_matches_gdc), all(meta$size_matches_gdc))
x <- matrix(0L, nrow(gene), nrow(meta), dimnames = list(gene$gene_id, meta$file_id))
audit <- vector("list", nrow(meta))
for (i in seq_len(nrow(meta))) {
  z <- fread(meta$path[i], skip = "gene_id", select = c("gene_id", "gene_name", "gene_type", "unstranded"), na.strings = "")
  z <- z[startsWith(gene_id, "ENSG")]
  stopifnot(identical(z$gene_id, gene$gene_id), identical(z$gene_name, gene$gene_name),
            identical(z$gene_type, gene$gene_type), !anyNA(z$unstranded),
            all(is.finite(z$unstranded)), all(z$unstranded >= 0), all(z$unstranded == floor(z$unstranded)),
            max(z$unstranded) <= .Machine$integer.max)
  x[, i] <- as.integer(z$unstranded)
  audit[[i]] <- data.table(file_id = meta$file_id[i], n_genes = nrow(z),
                           total_gene_reads = sum(z$unstranded), detected_genes = sum(z$unstranded > 0),
                           gene_order_matches = TRUE, annotation_matches = TRUE,
                           nonnegative_integer = TRUE)
  if (i %% 100L == 0L) { cat("READ", i, "/", nrow(meta), format(Sys.time()), "\n"); flush.console() }
}
fwrite(rbindlist(audit), file.path(out, "audit", "TCGA_all_file_count_checks.tsv"), sep = "\t")
geneAudit <- copy(gene)
geneAudit[, n_samples_nonzero := rowSums(x > 0)]
geneAudit[, total_counts := rowSums(x)]
fwrite(geneAudit, file.path(out, "audit", "TCGA_gene_detection_unfiltered.tsv"), sep = "\t")
saveRDS(list(counts = x, genes = gene, samples = meta,
             source = "unstranded raw gene counts, GENCODE v36, unfiltered"),
        file.path(out, "cache", "TCGA_unstranded_unfiltered_counts.rds"), compress = "gzip")
cat("COMPLETE", format(Sys.time()), "GENES", nrow(x), "FILES", ncol(x), "\n")
capture.output(sessionInfo(), file = file.path(out, "logs", "04_R_sessionInfo.txt"))
sink()
