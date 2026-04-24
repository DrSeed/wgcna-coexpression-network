#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(WGCNA); library(tidyverse); library(pheatmap) })
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
cat('=== WGCNA Co-expression Network Analysis ===\n')
datExpr0 <- read.csv('data/vst_counts.csv', row.names = 1, check.names = FALSE)
datExpr <- as.data.frame(t(datExpr0))
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
power <- sft$powerEstimate
if (is.na(power)) power <- 6
net <- blockwiseModules(datExpr, power = power, TOMType = 'unsigned', minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE, verbose = 3)
moduleColors <- labels2colors(net$colors)
traits <- read.csv('data/clinical.csv', row.names = 1)
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
moduleTraitCor <- cor(MEs, traits[rownames(datExpr), , drop = FALSE], use = 'p')
hub_genes <- data.frame(gene = colnames(datExpr), module = moduleColors)
write.csv(hub_genes, 'results/gene_module_assignments.csv', row.names = FALSE)
cat('WGCNA complete.\n')
