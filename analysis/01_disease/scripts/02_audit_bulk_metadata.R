options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
setDTthreads(2L)
root <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
sink(file.path(out, "logs", "02_bulk_metadata_audit.log"), split = TRUE)
cat("START", format(Sys.time()), "\n")
local <- fread(file.path(out, "metadata", "TCGA_local_file_ids.tsv"))
checks <- fread(file.path(out, "audit", "raw_input_local_md5.tsv"))
response <- fromJSON(file.path(out, "metadata", "GDC_files_metadata_all.json"), simplifyVector = FALSE)
rows <- list()
val <- function(x, k) if (is.null(x[[k]])) NA_character_ else as.character(x[[k]])
for (hit in response) {
  for (case in hit$cases) {
    for (sample in case$samples) {
      rows[[length(rows) + 1L]] <- data.table(
        file_id = val(hit, "file_id"), file_name = val(hit, "file_name"),
        gdc_md5 = val(hit, "md5sum"), gdc_bytes = as.numeric(val(hit, "file_size")),
        state = val(hit, "state"), workflow = val(hit$analysis, "workflow_type"),
        case_uuid = val(case, "case_id"), patient_id = val(case, "submitter_id"),
        project = val(case$project, "project_id"), sample_uuid = val(sample, "sample_id"),
        sample_barcode = val(sample, "submitter_id"), sample_type = val(sample, "sample_type"),
        sample_type_id = val(sample, "sample_type_id"), tissue_type = val(sample, "tissue_type"),
        gender = val(case$demographic, "gender"), race = val(case$demographic, "race"),
        age_at_diagnosis_days = if (length(case$diagnoses)) val(case$diagnoses[[1]], "age_at_diagnosis") else NA_character_,
        primary_diagnosis = if (length(case$diagnoses)) val(case$diagnoses[[1]], "primary_diagnosis") else NA_character_,
        ajcc_stage = if (length(case$diagnoses)) val(case$diagnoses[[1]], "ajcc_pathologic_stage") else NA_character_
      )
    }
  }
}
meta <- rbindlist(rows, fill = TRUE)
stopifnot(nrow(meta) == nrow(local), !anyDuplicated(meta$file_id), setequal(meta$file_id, local$file_id))
meta <- merge(local, meta, by = c("file_id", "file_name"), all.x = TRUE, sort = FALSE)
meta <- merge(meta, checks[, .(path, md5_local, bytes_local = bytes)], by = "path", all.x = TRUE, sort = FALSE)
meta[, md5_matches_gdc := md5_local == gdc_md5]
meta[, size_matches_gdc := bytes_local == gdc_bytes]
meta[, sample_barcode_type := substr(sample_barcode, 14L, 15L)]
stopifnot(!anyNA(meta$md5_matches_gdc), all(meta$md5_matches_gdc), all(meta$size_matches_gdc))
fwrite(meta, file.path(out, "metadata", "TCGA_verified_file_sample_metadata.tsv"), sep = "\t", na = "NA")
summary <- meta[, .(n_files = .N, n_patients = uniqueN(patient_id)), by = .(project, sample_type, gender)]
fwrite(summary, file.path(out, "audit", "TCGA_sample_summary.tsv"), sep = "\t")
print(summary)
cat("FILES", nrow(meta), "UNIQUE PATIENTS", uniqueN(meta$patient_id),
    "MD5 VERIFIED", sum(meta$md5_matches_gdc), "\n")
dup <- meta[, .(n_files = .N, files = paste(file_id, collapse = ";"),
                sample_barcodes = paste(sample_barcode, collapse = ";")), by = .(patient_id, sample_type)][n_files > 1L]
fwrite(dup, file.path(out, "audit", "TCGA_patient_condition_multiple_files.tsv"), sep = "\t")
print(table(meta$sample_type))
cat("MULTIPLE FILES SAME PATIENT CONDITION", nrow(dup), "\n")
pair <- meta[, .(types = paste(sort(unique(sample_type)), collapse = ";"), n_files = .N), by = patient_id]
fwrite(pair, file.path(out, "audit", "TCGA_patient_condition_coverage.tsv"), sep = "\t")
cat("TUMOUR NORMAL PAIRED PATIENTS", sum(grepl("Primary Tumor", pair$types) & grepl("Solid Tissue Normal", pair$types)), "\n")
cat("COMPLETE", format(Sys.time()), "\n")
sink()
