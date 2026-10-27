# Input audit only: no sample exclusions, differential expression or target prediction.
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "Asia/Shanghai")
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
for (d in c("audit", "metadata", "cache", "results", "figures", "logs", "scripts"))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
logfile <- file.path(out, "logs", "00_input_audit.log")
sink(logfile, split = TRUE)
cat("START", format(Sys.time()), "\nROOT", root, "\n")
allfiles <- list.files(file.path(root, "原始数据"), recursive = TRUE, full.names = TRUE)
fi <- file.info(allfiles)
manifest <- data.table(path = allfiles, bytes = fi$size,
                       modified_local = format(fi$mtime, "%Y-%m-%d %H:%M:%S"))
fwrite(manifest, file.path(out, "audit", "raw_input_inventory.tsv"), sep = "\t")
cat("INPUT FILES", nrow(manifest), "BYTES", sum(manifest$bytes), "\n")

softpath <- file.path(root, "原始数据", "单细胞", "GSE161529_family.soft.gz")
s <- readLines(gzfile(softpath), warn = FALSE)
ss <- which(grepl("^\\^SAMPLE = ", s))
samplelist <- vector("list", length(ss))
charlist <- vector("list", length(ss))
getfield <- function(b, field) {
  v <- b[startsWith(b, paste0(field, " = "))]
  paste(sub("^[^=]+ = ", "", v), collapse = " | ")
}
for (i in seq_along(ss)) {
  b <- s[ss[i]:(if (i < length(ss)) ss[i + 1L] - 1L else length(s))]
  gsm <- sub("^\\^SAMPLE = ", "", b[1])
  ch <- b[startsWith(b, "!Sample_characteristics_ch1 = ")]
  ch <- sub("^!Sample_characteristics_ch1 = ", "", ch)
  ck <- sub(":.*$", "", ch)
  cv <- sub("^[^:]+: ?", "", ch)
  charlist[[i]] <- data.table(gsm = gsm, characteristic = ck, value = cv)
  one <- data.table(gsm = gsm, title = getfield(b, "!Sample_title"),
                   source = getfield(b, "!Sample_source_name_ch1"),
                   treatment = getfield(b, "!Sample_treatment_protocol_ch1"),
                   processing = getfield(b, "!Sample_data_processing"))
  for (j in seq_along(ck)) one[[make.names(ck[j])]] <- cv[j]
  samplelist[[i]] <- one
}
samples <- rbindlist(samplelist, fill = TRUE)
fwrite(samples, file.path(out, "metadata", "GSE161529_sample_metadata_unfiltered.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(charlist), file.path(out, "metadata", "GSE161529_characteristics_long.tsv"), sep = "\t")
writeLines(s[seq_len(ss[1] - 1L)], file.path(out, "metadata", "GSE161529_series_platform_metadata.txt"), useBytes = TRUE)
cat("SOFT PROFILES", nrow(samples), "\n")
for (key in intersect(c("cancer.type", "cell.population", "source"), names(samples))) {
  cat("FIELD", key, "\n"); print(samples[, .N, by = key])
}
if ("patient" %in% names(samples)) {
  dup <- samples[, .(n_profiles = .N, gsm = paste(gsm, collapse = ";"),
                    titles = paste(title, collapse = " | ")), by = patient][n_profiles > 1L]
  fwrite(dup, file.path(out, "audit", "GSE161529_patients_with_multiple_profiles.tsv"), sep = "\t")
  cat("UNIQUE PATIENT FIELD", uniqueN(samples$patient), "MULTI-PROFILE", nrow(dup), "\n")
}
features <- as.data.table(read.delim(gzfile(file.path(root, "原始数据", "单细胞", "GSE161529_features.tsv.gz")),
                  header = FALSE, col.names = c("ensembl_id", "gene_symbol", "feature_type")))
stopifnot(ncol(features) == 3L, !anyDuplicated(features$ensembl_id))
fwrite(features, file.path(out, "metadata", "GSE161529_features.tsv"), sep = "\t")
tarpath <- file.path(root, "原始数据", "单细胞", "GSE161529_RAW.tar")
members <- utils::untar(tarpath, list = TRUE)
stopifnot(!any(grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)|^[/\\\\]|^[A-Za-z]:", members)))
writeLines(members, file.path(out, "audit", "GSE161529_tar_members.txt"))
mat <- members[grepl("-matrix\\.mtx\\.gz$", members)]
bar <- members[grepl("-barcodes\\.tsv\\.gz$", members)]
matkey <- sub("-matrix\\.mtx\\.gz$", "", basename(mat))
barkey <- sub("-barcodes\\.tsv\\.gz$", "", basename(bar))
stopifnot(!anyDuplicated(matkey), !anyDuplicated(barkey), setequal(matkey, barkey))
pair <- data.table(profile_key = matkey, gsm = sub("_.*$", "", matkey), matrix_member = mat,
                   barcode_member = bar[match(matkey, barkey)])
stopifnot(setequal(pair$gsm, samples$gsm), length(mat) == nrow(samples))
fwrite(pair, file.path(out, "audit", "GSE161529_tar_file_pairs.tsv"), sep = "\t")
cat("FEATURES", nrow(features), "TAR MEMBERS", length(members), "MATRIX/BARCODE PAIRS", nrow(pair), "\n")

bulkpaths <- manifest$path[grepl("\\.rna_seq\\.augmented_star_gene_counts\\.tsv$", manifest$path)]
bulk <- data.table(file_id = basename(dirname(bulkpaths)), file_name = basename(bulkpaths), path = bulkpaths)
stopifnot(!anyDuplicated(bulk$file_id), !anyDuplicated(bulk$file_name))
fwrite(bulk, file.path(out, "metadata", "TCGA_local_file_ids.tsv"), sep = "\t")
first <- fread(bulkpaths[1], skip = "gene_id", na.strings = "")
gene <- first[grepl("^ENSG", gene_id), .(gene_id, gene_name, gene_type)]
gene[, ensembl_stable := sub("\\.[0-9]+$", "", gene_id)]
stopifnot(!anyDuplicated(gene$gene_id), !anyDuplicated(gene$ensembl_stable))
fwrite(gene, file.path(out, "metadata", "TCGA_first_file_GENCODEv36_gene_annotation.tsv"), sep = "\t")
fwrite(gene[, .N, by = gene_type], file.path(out, "audit", "TCGA_first_file_biotypes.tsv"), sep = "\t")
cat("BULK FILES", nrow(bulk), "FIRST FILE GENES", nrow(gene), "\n")
cat("Computing MD5 on original inputs; does not imply comparison with repository MD5.\n")
manifest[, md5_local := unname(tools::md5sum(path))]
fwrite(manifest, file.path(out, "audit", "raw_input_local_md5.tsv"), sep = "\t")
pkgs <- c("data.table", "Seurat", "SeuratObject", "Matrix", "limma", "edgeR", "infercnv", "copykat", "jsonlite", "ggplot2")
pv <- rbindlist(lapply(pkgs, function(x) data.table(package = x,
   version = tryCatch(as.character(packageVersion(x)), error = function(e) NA_character_),
   path = tryCatch(find.package(x), error = function(e) NA_character_))))
fwrite(pv, file.path(out, "audit", "R_package_locations.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "logs", "00_R_sessionInfo.txt"))
cat("COMPLETE", format(Sys.time()), "\n")
sink()
