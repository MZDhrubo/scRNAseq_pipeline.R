# scRNAseq_pipeline.R
#Complete Seurat workflow for PBMC3k scRNA-seq analysis: #   QC, normalization, PCA, clustering, UMAP, marker detection, #   cell-type annotation support, publication-ready plots, #   and interactive 3D UMAP/PCA visualization.
# PBMC3k Single-cell RNA-seq Analysis using Seurat

This repository contains a complete single-cell RNA-seq analysis pipeline using the PBMC3k dataset and the Seurat R package.

## Workflow

The pipeline includes:

1. Data loading using SeuratData
2. Quality control
3. Cell filtering
4. Normalization
5. Highly variable gene detection
6. Scaling
7. PCA
8. Clustering
9. 2D UMAP visualization
10. Marker gene discovery
11. Canonical PBMC marker visualization
12. Manual annotation template
13. Optional SingleR automatic annotation
14. Interactive 3D UMAP visualization
15. Interactive 3D PCA visualization

## Dataset

The analysis uses the PBMC3k dataset, a commonly used peripheral blood mononuclear cell single-cell RNA-seq dataset.

## Main Outputs

Results are saved in the `results/` directory.

### Figures

- QC violin plots
- QC scatter plots
- Variable feature plot
- PCA plot
- Elbow plot
- UMAP cluster plot
- Canonical marker DotPlot
- Marker FeaturePlots
- Marker violin plots
- Marker heatmap

### Tables

- Cluster sizes
- All cluster markers
- Top 10 markers per cluster
- Top 5 markers per cluster

### Interactive HTML

- 3D UMAP colored by clusters
- 3D UMAP colored by manual annotation
- 3D PCA colored by clusters

## How to Run

Open R or RStudio and run:

```r
source("PBMC3k_scRNAseq_pipeline.R")
