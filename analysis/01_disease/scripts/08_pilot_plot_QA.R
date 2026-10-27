options(stringsAsFactors = FALSE)
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(TRUE)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
gsm <- args[2]
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
so <- readRDS(file.path(out, "cache", "seurat_pilot", paste0(gsm, "_pilot_seurat.rds")))
Idents(so) <- "RNA_snn_res.0.1"
f <- fread(file.path(out, "metadata", "GSE161529_features.tsv"))
panel <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "TP63", "KIT", "ELF5", "ESR1", "PGR",
           "PTPRC", "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "MS4A1", "CD79A", "JCHAIN", "MZB1", "LYZ", "LST1", "FCER1G", "CD68",
           "DCN", "COL1A1", "COL1A2", "LUM", "PECAM1", "VWF", "CLDN5", "RGS5", "PDGFRB", "CSPG4", "MKI67", "TOP2A")
mapping <- f[gene_symbol %in% panel][!duplicated(gene_symbol)]
mapping <- mapping[match(panel[panel %in% mapping$gene_symbol], gene_symbol)]
panel_ids <- setNames(mapping$ensembl_id, mapping$gene_symbol)
p <- DotPlot(so, features = unname(panel_ids), cols = c("#D9E6F2", "#B2182B")) +
     ggplot2::scale_x_discrete(labels = setNames(names(panel_ids), unname(panel_ids))) + RotatedAxis() +
     ggplot2::ggtitle(paste(gsm, "descriptive lineage markers; identities not finalized"))
ggplot2::ggsave(file.path(out, "figures", paste0(gsm, "_pilot_marker_panel.png")), p, width = 14, height = 6, dpi = 160, bg="white")
p <- DimPlot(so, label = TRUE) + ggplot2::ggtitle(paste(gsm, "pilot Seurat clusters, resolution 0.1"))
ggplot2::ggsave(file.path(out, "figures", paste0(gsm, "_pilot_UMAP.png")), p, width = 9, height = 6, dpi = 160, bg="white")
mat <- GetAssayData(so, assay = "RNA", layer = "data")[unname(panel_ids), , drop = FALSE]
cl <- as.character(Idents(so))
stats <- rbindlist(lapply(sort(unique(cl)), function(k) {
  ix <- cl == k
  data.table(cluster = k, gene_symbol = names(panel_ids), ensembl_id = unname(panel_ids), n_cells = sum(ix),
             mean_log_normalized = Matrix::rowMeans(mat[, ix, drop = FALSE]),
             pct_expressing = Matrix::rowMeans(mat[, ix, drop = FALSE] > 0))
}))
fwrite(stats, file.path(out, "audit", paste0(gsm, "_pilot_lineage_marker_summary.tsv")), sep = "\t")
cat("PLOT QA READY", gsm, "CELLS", ncol(so), "CLUSTERS", length(unique(cl)), "\n")
