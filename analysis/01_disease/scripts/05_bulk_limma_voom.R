# User-confirmed 2026-10-23: all samples, patient correlation; paired sensitivity.
options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("D:/R4.5.2/R-4.5.2/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
setDTthreads(2L)
set.seed(20261023)
root <- normalizePath(commandArgs(TRUE)[1], winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs", "第一阶段_2026-10-23")
source(file.path(out,"scripts","bulk_block_gls.R"),encoding="UTF-8")
fwrite(validate_fast_pairblock_gls(), file.path(out,"audit","Bulk_block_GLS_vs_limma_validation.tsv"),sep="\t")
sink(file.path(out, "logs", "05_bulk_limma_voom.log"), split = TRUE)
cat("START", format(Sys.time()), "\n")
inp <- readRDS(file.path(out, "cache", "TCGA_unstranded_unfiltered_counts.rds"))
meta <- as.data.table(inp$samples)
genes <- as.data.table(inp$genes)
stopifnot(identical(colnames(inp$counts), meta$file_id), identical(rownames(inp$counts), genes$gene_id))
meta[, condition := factor(ifelse(sample_type == "Primary Tumor", "Tumor", "Normal"), levels = c("Normal", "Tumor"))]
meta[, patient_block := factor(patient_id)]
design <- model.matrix(~ condition, data = meta)
stopifnot(qr(design)$rank == ncol(design))
dge <- DGEList(counts = inp$counts)
keep <- filterByExpr(dge, group = meta$condition)
filt <- copy(genes)
filt[, filterByExpr_keep := keep]
# CPM summaries expose genes at the filtering boundary; not an additional filter.
cpm_cut <- 10 / median(dge$samples$lib.size) * 1e6
filt[, n_samples_at_count_cpm_cutoff := rowSums(cpm(dge) >= cpm_cut)]
filt[, total_count := rowSums(dge$counts)]
fwrite(filt, file.path(out, "audit", "TCGA_filterByExpr_decisions.tsv"), sep = "\t")
cat("FILTER", sum(keep), "/", length(keep), "CPM CUTOFF", cpm_cut, "\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
rm(inp); gc(verbose = FALSE)
pdf(file.path(out, "figures", "Bulk_QC_MDS_voom.pdf"), width = 8, height = 6)
condition_col <- c(Normal = "#4E79A7", Tumor = "#E15759")
cat("MDS COMMON TOP500 START", format(Sys.time()), "\n"); flush.console()
plotMDS(dge, top = 500, gene.selection = "common", labels = rep("", nrow(meta)), col = condition_col[as.character(meta$condition)], pch = 16, cex = .5,
        main = "TCGA-BRCA: unfiltered sample membership")
legend("topright", legend = names(condition_col), col = condition_col, pch = 16, bty = "n")
boxplot(log10(dge$samples$lib.size) ~ meta$condition, col = condition_col,
        ylab = "log10 library size after gene filtering", xlab = "", main = "Library size")
v0 <- voom(dge, design, plot = TRUE)
cat("DUPLICATE CORRELATION START", format(Sys.time()), "\n"); flush.console()
pairIdsForCorrelation <- meta[, .(n_condition = uniqueN(condition)), by = patient_id][n_condition == 2L, patient_id]
pairIxForCorrelation <- which(meta$patient_id %in% pairIdsForCorrelation)
# The within-patient correlation is estimated on the 113 actual paired blocks.
# All 1,208 samples remain in the main expression/effect model.
corfit <- duplicateCorrelation(v0[,pairIxForCorrelation], design[pairIxForCorrelation,,drop=FALSE],
                               block = droplevels(meta$patient_block[pairIxForCorrelation]))
cat("CONSENSUS CORRELATION", corfit$consensus.correlation, format(Sys.time()), "\n")
rm(v0); gc(verbose = FALSE)
v <- voom(dge, design, plot = TRUE, block = meta$patient_block, correlation = corfit$consensus.correlation)
fit <- fast_pairblock_gls(v$E,v$weights,as.integer(meta$condition=="Tumor"),meta$patient_block,corfit$consensus.correlation)
checkGenes <- unique(round(seq(1,nrow(v$E),length.out=10)))
standardCheck <- lmFit(v[checkGenes,],design,block=meta$patient_block,correlation=corfit$consensus.correlation)
realChecks <- data.table(quantity=c("coefficient","stdev_unscaled","sigma"),
 max_absolute_error=c(max(abs(fit$coefficients[checkGenes,]-standardCheck$coefficients)),
 max(abs(fit$stdev.unscaled[checkGenes,]-standardCheck$stdev.unscaled)),max(abs(fit$sigma[checkGenes]-standardCheck$sigma))))
stopifnot(all(realChecks$max_absolute_error<1e-7))
fwrite(realChecks,file.path(out,"audit","Bulk_block_GLS_real_genes_vs_limma.tsv"),sep="\t")
fit <- eBayes(fit)
dev.off()
tt <- as.data.table(topTable(fit, coef = "conditionTumor", number = Inf, sort.by = "none"), keep.rownames = "gene_id")
tt <- merge(genes[keep], tt, by = "gene_id", all.x = TRUE, sort = FALSE)
tt[, `:=`(contrast = "Tumor_minus_Normal", model = "voom_patient_block_duplicateCorrelation",
           exploratory_de = abs(logFC) > 0.585 & P.Value < .05,
           bh_significant = adj.P.Val < .05,
           direction = fifelse(logFC > 0, "Up_in_bulk_tumor", "Down_in_bulk_tumor"))]
stopifnot(nrow(tt) == sum(keep), !anyDuplicated(tt$gene_id), !anyNA(tt$P.Value))
fwrite(tt, file.path(out, "results", "Bulk_all_tested_genes.tsv"), sep = "\t")
fwrite(tt[exploratory_de == TRUE], file.path(out, "results", "Bulk_exploratory_DE_P005_FC1p5.tsv"), sep = "\t")
fwrite(tt[exploratory_de == TRUE & gene_type == "protein_coding"],
       file.path(out, "results", "Bulk_exploratory_DE_protein_coding.tsv"), sep = "\t")
modelinfo <- data.table(model = "all_1208_patient_block", n_normal = sum(meta$condition == "Normal"),
                       n_tumor = sum(meta$condition == "Tumor"), n_patients = uniqueN(meta$patient_id),
                       correlation = corfit$consensus.correlation, correlation_calibration_patients=113L,
                       n_genes_tested = nrow(tt),
                       n_exploratory_de = sum(tt$exploratory_de),
                       n_de_protein_coding = sum(tt$exploratory_de & tt$gene_type == "protein_coding"))
fwrite(modelinfo, file.path(out, "results", "Bulk_model_summary.tsv"), sep = "\t")
print(modelinfo)

# Reuse the main tested gene universe; refit TMM on paired samples.
pairIds <- meta[, .(n_condition = uniqueN(condition)), by = patient_id][n_condition == 2L, patient_id]
ix <- which(meta$patient_id %in% pairIds)
pm <- meta[ix]
stopifnot(length(pairIds) == 113L, length(ix) == 226L)
pd <- dge[, ix, keep.lib.sizes = FALSE]
pd <- calcNormFactors(pd, method = "TMM")
paired_design <- model.matrix(~ patient_block + condition, data = droplevels(as.data.frame(pm)))
stopifnot(qr(paired_design)$rank == ncol(paired_design))
cat("PAIRED SENSITIVITY START", format(Sys.time()), "\n"); flush.console()
pv <- voom(pd, paired_design, plot = FALSE)
pn <- match(pairIds,pm[condition=="Normal",patient_id])
ptum <- match(pairIds,pm[condition=="Tumor",patient_id])
normalIndices <- which(pm$condition=="Normal")[pn]
tumorIndices <- which(pm$condition=="Tumor")[ptum]
stopifnot(identical(pm$patient_id[normalIndices],pm$patient_id[tumorIndices]))
diffE <- pv$E[,tumorIndices,drop=FALSE]-pv$E[,normalIndices,drop=FALSE]
diffW <- 1/(1/pv$weights[,tumorIndices,drop=FALSE]+1/pv$weights[,normalIndices,drop=FALSE])
pairedFast <- lmFit(diffE,design=matrix(1,ncol(diffE),1,dimnames=list(NULL,"conditionTumor")),weights=diffW)
# Check the exact paired-difference reduction against the full patient fixed-effect model.
pairStandard <- lmFit(pv[checkGenes,],paired_design)
pairCheck <- data.table(quantity=c("coefficient","stdev_unscaled","sigma"),
 max_absolute_error=c(max(abs(pairedFast$coefficients[checkGenes,1]-pairStandard$coefficients[,"conditionTumor"])),
 max(abs(pairedFast$stdev.unscaled[checkGenes,1]-pairStandard$stdev.unscaled[,"conditionTumor"])),
 max(abs(pairedFast$sigma[checkGenes]-pairStandard$sigma))))
stopifnot(all(pairCheck$max_absolute_error<1e-7))
fwrite(pairCheck,file.path(out,"audit","Bulk_paired_difference_vs_full_model_validation.tsv"),sep="\t")
pf <- eBayes(pairedFast)
pt <- as.data.table(topTable(pf, coef = "conditionTumor", number = Inf, sort.by = "none"), keep.rownames = "gene_id")
pt <- merge(genes[keep], pt, by = "gene_id", all.x = TRUE, sort = FALSE)
pt[, exploratory_de := abs(logFC) > .585 & P.Value < .05]
pt[, model := "113_pairs_patient_fixed_effect"]
fwrite(pt, file.path(out, "results", "Bulk_paired113_sensitivity_all_genes.tsv"), sep = "\t")
comparison <- merge(tt[, .(gene_id, bulk_logFC = logFC, bulk_P = P.Value, bulk_BH = adj.P.Val, bulk_DE = exploratory_de)],
                    pt[, .(gene_id, paired_logFC = logFC, paired_P = P.Value, paired_BH = adj.P.Val, paired_DE = exploratory_de)], by = "gene_id")
comparison[, direction_agrees := sign(bulk_logFC) == sign(paired_logFC)]
fwrite(comparison, file.path(out, "results", "Bulk_main_vs_paired_sensitivity.tsv"), sep = "\t")
cat("PAIRED LOGFC SPEARMAN", cor(comparison$bulk_logFC, comparison$paired_logFC, method = "spearman"), "\n")

pdf(file.path(out, "figures", "Bulk_DE_exploratory_volcano.pdf"), width = 8, height = 6)
with(tt, plot(logFC, -log10(pmax(P.Value, .Machine$double.xmin)), pch = 16, cex = .35,
              col = ifelse(exploratory_de, ifelse(logFC > 0, "#E15759", "#4E79A7"), "#B7B7B7"),
              xlab = "Tumor - Normal log2 fold change", ylab = "-log10 raw P", main = "Exploratory threshold: raw P < 0.05, |log2FC| > 0.585"))
abline(v = c(-.585, .585), h = -log10(.05), lty = 2, col = "#777777")
dev.off()
saveRDS(list(model = modelinfo, genes = genes[keep], samples = meta,
             design = design, correlation = corfit$consensus.correlation,
             main_fit = fit, paired_fit = pf, filter = filt,
             note = "Exploratory raw P threshold as user specified; no enrichment filtering."),
        file.path(out, "cache", "Bulk_limma_fits_and_design.rds"), compress = "gzip")
fwrite(meta, file.path(out, "metadata", "TCGA_model_sample_order.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "logs", "05_R_sessionInfo.txt"))
cat("COMPLETE", format(Sys.time()), "\n")
sink()
