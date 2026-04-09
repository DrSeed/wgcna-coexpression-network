# WGCNA Co-expression Network Pipeline

Identify co-expressed gene modules and hub genes using Weighted Gene Co-expression Network Analysis.

## Features
- Soft-thresholding power selection with scale-free topology fitting
- Automatic module detection via dynamic tree cutting
- Module-trait correlation analysis
- Hub gene identification (module membership + gene significance)
- Network export for Cytoscape visualisation

## Usage
```bash
Rscript wgcna_pipeline.R --expression data/vst_counts.csv --traits data/clinical.csv
```

## Output
- Module assignments and eigengenes
- Module-trait correlation heatmap
- Hub gene lists per module
- Cytoscape-compatible edge/node tables
