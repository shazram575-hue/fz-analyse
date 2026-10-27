# Technical pilot on a named library, no final cell type or malignancy assignment.
# Reference: Pal/Chen HumanBreast10X per-library Seurat preparation; current Seurat 5.
options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
setDTthreads(2L)
future::plan("sequential")
set.seed(20261023)
args <- getOption("workflow3.library.args", commandArgs(TRUE))
stopifnot(length(args) == 2L)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
gsm <- args[2]
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
pout <- file.path(out, "cache", "seurat_pilot")
dir.create(pout, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out, "logs", paste0("07_pilot_", gsm, ".log")), split = TRUE)
cat("START PILOT", gsm, format(Sys.time()), "\n")
meta <- fread(file.path(out, "metadata", "GSE161529_author_QC_and_analysis_membership.tsv"), colClasses = "character")
one <- meta[gsm == args[2]]
stopifnot(nrow(one) == 1L, one$include_main == "TRUE" || one$analysis_role == "normal_epi_sensitivity")
pairs <- fread(file.path(out, "audit", "GSE161529_tar_file_pairs.tsv"))
pair <- pairs[gsm == args[2]]
f <- fread(file.path(out, "metadata", "GSE161529_features.tsv"))
kept <- readRDS(file.path(out, "cache", "GSE161529_author_QC_kept_barcodes_by_gsm.rds"))[[gsm]]
unpacked <- file.path(out, "cache", "GSE161529_unpacked")
con <- gzfile(file.path(unpacked, pair$matrix_member), "rt")
x <- readMM(con); close(con)
x <- as(x, "CsparseMatrix")
barcodes <- readLines(gzfile(file.path(unpacked, pair$barcode_member)), warn = FALSE)
stopifnot(length(barcodes) == ncol(x), nrow(x) == nrow(f), all(kept %in% barcodes))
rownames(x) <- f$ensembl_id
colnames(x) <- barcodes
x <- x[, kept, drop = FALSE]
colnames(x) <- paste(gsm, kept, sep = "__")
cat("AUTHOR QC CELLS", ncol(x), "\n")

cat("DOUBLET MARKING START", format(Sys.time()), "\n"); flush.console()
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scDblFinder))
sce <- SingleCellExperiment(list(counts = x))
set.seed(20261023)
sce <- scDblFinder(sce, BPPARAM = BiocParallel::SerialParam(RNGseed = 20261023))
doublets <- data.table(cell_id = colnames(sce), scDblFinder_score = sce$scDblFinder.score,
                        scDblFinder_class = as.character(sce$scDblFinder.class))
fwrite(doublets, file.path(out, "audit", paste0("scDblFinder_", gsm, "_pilot.tsv")), sep = "\t")
rm(sce); gc(verbose = FALSE)

so <- CreateSeuratObject(x, project = gsm, min.cells = 0, min.features = 0)
so$gsm <- gsm
so$author_sample <- one$SampleName
so$scDblFinder_class <- doublets$scDblFinder_class[match(colnames(so), doublets$cell_id)]
so[["percent.mt"]] <- PercentageFeatureSet(so, features = f$ensembl_id[startsWith(f$gene_symbol, "MT-")])
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
so <- RunPCA(so, features = VariableFeatures(so), npcs = 30, seed.use = 20261023, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = c(0.1, 0.4, 0.8), algorithm = 1, random.seed = 20261023, verbose = FALSE)
Idents(so) <- "RNA_snn_res.0.1"
so <- RunUMAP(so, dims = 1:30, seed.use = 20261023, verbose = FALSE)
cat("SEURAT PILOT DONE", format(Sys.time()), "\n")
md <- as.data.table(so[[]], keep.rownames = "cell_id")
fwrite(md, file.path(out, "audit", paste0("Seurat_", gsm, "_pilot_cell_metadata.tsv")), sep = "\t")
panel <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "TP63", "KIT", "ELF5", "ESR1", "PGR",
           "PTPRC", "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "MS4A1", "CD79A", "JCHAIN", "MZB1",
           "LYZ", "LST1", "FCER1G", "CD68", "DCN", "COL1A1", "COL1A2", "LUM",
           "PECAM1", "VWF", "CLDN5", "RGS5", "PDGFRB", "CSPG4", "MKI67", "TOP2A")
mapping <- f[gene_symbol %in% panel]
# Avoid ambiguous symbols in descriptive marker plot; no gene is excluded from counts.
mapping <- mapping[!duplicated(gene_symbol)]
mapping <- mapping[match(panel[panel %in% mapping$gene_symbol], gene_symbol)]
panel_ids <- setNames(mapping$ensembl_id, mapping$gene_symbol)
pdf(file.path(out, "figures", paste0("Seurat_", gsm, "_technical_pilot.pdf")), width = 11, height = 7)
print(DimPlot(so, reduction = "umap", label = TRUE) + ggplot2::ggtitle(paste(gsm, "pilot clusters; no final cell identities")))
print(DimPlot(so, reduction = "umap", group.by = "scDblFinder_class") + ggplot2::ggtitle("Doublet calls marked, no additional exclusion"))
print(DotPlot(so, features = unname(panel_ids), cols = c("#D9E6F2", "#B2182B")) +
        ggplot2::scale_x_discrete(labels = setNames(names(panel_ids), unname(panel_ids))) + RotatedAxis())
print(ElbowPlot(so, ndims = 30))
dev.off()
saveRDS(so, file.path(pout, paste0(gsm, "_pilot_seurat.rds")), compress = "gzip")
if(file.exists(file.path(out,"cache","celldex_reference","HPCA_2024-02-26_reference.rds"))) {
  source(file.path(out,"scripts","annotate_library.R"),local=TRUE,encoding="UTF-8")
  annotate_library(so,f,root,gsm)
}
capture.output(sessionInfo(), file = file.path(out, "logs", paste0("07_pilot_", gsm, "_sessionInfo.txt")))
cat("COMPLETE PILOT", gsm, format(Sys.time()), "\n")
sink()
