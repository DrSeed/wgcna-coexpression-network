#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(pheatmap)
})

allowWGCNAThreads()
options(stringsAsFactors = FALSE)

cat("=== WGCNA Co-expression Network Analysis ===\n")

# ── Load expression data (VST or log2-normalised) ─────────────────────
expr_file <- "data/vst_counts.csv"
traits_file <- "data/clinical.csv"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

datExpr0 <- read.csv(expr_file, row.names = 1, check.names = FALSE)
datExpr <- as.data.frame(t(datExpr0))

# ── QC: check for outlier samples ─────────────────────────────────────
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- hclust(dist(datExpr), method = "average")
png(file.path(output_dir, "sample_dendrogram.png"), width = 1000, height = 400)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()

# ── Soft-thresholding power selection ──────────────────────────────────
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
power <- sft$powerEstimate
if (is.na(power)) power <- 6
cat(sprintf("Selected soft-thresholding power: %d\n", power))

png(file.path(output_dir, "soft_threshold.png"), width = 800, height = 400)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold", ylab = "Scale Free Topology Model Fit",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold", ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# ── Network construction and module detection ──────────────────────────
net <- blockwiseModules(
  datExpr, power = power,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
cat(sprintf("Detected %d modules\n", length(unique(moduleColors)) - 1))

# ── Module-trait correlations ──────────────────────────────────────────
traits <- read.csv(traits_file, row.names = 1)
traits <- traits[rownames(datExpr), , drop = FALSE]

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

png(file.path(output_dir, "module_trait_heatmap.png"), width = 800, height = 600)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits),
  yLabels = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  main = "Module-trait relationships"
)
dev.off()

# ── Hub gene identification ────────────────────────────────────────────
hub_genes <- data.frame(gene = colnames(datExpr), module = moduleColors)
write.csv(hub_genes, file.path(output_dir, "gene_module_assignments.csv"), row.names = FALSE)

# Export for Cytoscape
TOM <- TOMsimilarityFromExpr(datExpr, power = power)
for (mod in unique(moduleColors)) {
  if (mod == "grey") next
  inModule <- moduleColors == mod
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(colnames(datExpr)[inModule], colnames(datExpr)[inModule])
  cyt <- exportNetworkToCytoscape(modTOM, edgeFile = file.path(output_dir, paste0("edges_", mod, ".txt")),
    nodeFile = file.path(output_dir, paste0("nodes_", mod, ".txt")),
    weighted = TRUE, threshold = 0.02)
}

cat("WGCNA pipeline complete.\n")
