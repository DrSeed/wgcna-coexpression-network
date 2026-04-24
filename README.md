# WGCNA Co-expression Network Pipeline

> Differential expression tells you WHAT changed. Co-expression networks tell you WHAT CHANGED TOGETHER, which is usually the more interesting question.

## Why Co-expression Networks Matter

You ran DESeq2. You got 2,000 significant genes. Now what?

A list of genes is not biology. Biology is about relationships: which genes are co-regulated, which modules respond together to your condition, and which hub genes sit at the centre driving the response. That hub gene? That's your next paper. Or your drug target.

## How WGCNA Actually Works (The Intuition)

1. **Building a correlation matrix** between all gene pairs
2. **Raising correlations to a soft-thresholding power** to emphasise strong connections and suppress noise
3. **Clustering genes into modules** based on co-expression patterns
4. **Correlating modules with clinical traits** to find which module is associated with your phenotype

## Usage
```bash
Rscript wgcna_pipeline.R --expression data/vst_counts.csv --traits data/clinical.csv
```

## Critical Things People Get Wrong

**Sample size matters here more than in DEA.** WGCNA needs at least 15 samples. With fewer, your correlation estimates are noisy and your modules won't be reproducible.

**Use normalised, variance-stabilised data, not raw counts.** Feed VST or rlog-transformed values. Raw counts will give you garbage networks.

**Outlier samples will destroy your network.** One outlier can create spurious correlations across thousands of gene pairs. The pipeline includes a sample dendrogram so you can catch this.

## From Modules to Biology

Once you identify a module correlated with your trait, the hub genes are your candidates for follow-up. Export to Cytoscape, overlay with your DEGs, and you've got a systems-level view no single-gene analysis could give you.
