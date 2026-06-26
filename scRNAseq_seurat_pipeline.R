############################################################
# Generic Single-cell RNA-seq Analysis Pipeline using Seurat
# Author: Your Name
#
# Description:
#   This script performs a standard scRNA-seq workflow:
#   data loading, QC, filtering, normalization, PCA,
#   clustering, UMAP, marker discovery, and 3D visualization.
#
# Input supported:
#   1. 10x Genomics folder
#   2. 10x Genomics .h5 file
#   3. CSV count matrix
#   4. Existing Seurat .rds object
############################################################


############################################################
# 0. User Settings
############################################################

# Change this path for your dataset
input_path <- "data/filtered_feature_bc_matrix"

# Choose one:
# "10x_folder", "10x_h5", "csv", "seurat_rds"
input_type <- "10x_folder"

# Project/sample name
project_name <- "scRNAseq_project"

# QC thresholds
min_features <- 200
max_features <- 2500
max_percent_mt <- 5

# Analysis parameters
num_variable_features <- 2000
dims_use <- 1:10
cluster_resolution <- 0.5

# Species mitochondrial gene pattern
# Human: "^MT-"
# Mouse: "^mt-"
mito_pattern <- "^MT-"


############################################################
# 1. Setup
############################################################

options(timeout = 1000)

required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "plotly",
  "htmlwidgets"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(plotly)
library(htmlwidgets)

# Create output folders
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/interactive", showWarnings = FALSE, recursive = TRUE)


############################################################
# 2. Load Data
############################################################

if (input_type == "10x_folder") {
  
  counts <- Read10X(data.dir = input_path)
  
  seu <- CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = 3,
    min.features = 200
  )
  
} else if (input_type == "10x_h5") {
  
  counts <- Read10X_h5(filename = input_path)
  
  seu <- CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = 3,
    min.features = 200
  )
  
} else if (input_type == "csv") {
  
  counts <- read.csv(
    input_path,
    row.names = 1,
    check.names = FALSE
  )
  
  seu <- CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = 3,
    min.features = 200
  )
  
} else if (input_type == "seurat_rds") {
  
  seu <- readRDS(input_path)
  
} else {
  
  stop("Invalid input_type. Use: 10x_folder, 10x_h5, csv, or seurat_rds.")
}

print(seu)


############################################################
# 3. Quality Control
############################################################

seu[["percent.mt"]] <- PercentageFeatureSet(
  seu,
  pattern = mito_pattern
)

qc_violin <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

ggsave(
  filename = "results/figures/01_QC_violin.png",
  plot = qc_violin,
  width = 10,
  height = 5,
  dpi = 300
)

qc_scatter_1 <- FeatureScatter(
  seu,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

qc_scatter_2 <- FeatureScatter(
  seu,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

qc_scatter <- qc_scatter_1 + qc_scatter_2

ggsave(
  filename = "results/figures/02_QC_scatter.png",
  plot = qc_scatter,
  width = 10,
  height = 5,
  dpi = 300
)

# Filter cells
seu <- subset(
  seu,
  subset = nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < max_percent_mt
)

print(seu)


############################################################
# 4. Normalization and Variable Features
############################################################

seu <- NormalizeData(
  seu,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = num_variable_features
)

top10_variable_genes <- head(VariableFeatures(seu), 10)

variable_plot <- VariableFeaturePlot(seu)

variable_plot_labeled <- LabelPoints(
  plot = variable_plot,
  points = top10_variable_genes,
  repel = TRUE
)

ggsave(
  filename = "results/figures/03_variable_features.png",
  plot = variable_plot_labeled,
  width = 8,
  height = 6,
  dpi = 300
)


############################################################
# 5. Scaling and PCA
############################################################

all_genes <- rownames(seu)

seu <- ScaleData(
  seu,
  features = all_genes
)

seu <- RunPCA(
  seu,
  features = VariableFeatures(object = seu)
)

pca_plot <- DimPlot(
  seu,
  reduction = "pca"
) +
  ggtitle("PCA plot")

ggsave(
  filename = "results/figures/04_PCA.png",
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

elbow_plot <- ElbowPlot(
  seu,
  ndims = 30
)

ggsave(
  filename = "results/figures/05_elbow_plot.png",
  plot = elbow_plot,
  width = 7,
  height = 5,
  dpi = 300
)


############################################################
# 6. Clustering
############################################################

seu <- FindNeighbors(
  seu,
  dims = dims_use
)

seu <- FindClusters(
  seu,
  resolution = cluster_resolution
)

cluster_sizes <- as.data.frame(table(Idents(seu)))
colnames(cluster_sizes) <- c("cluster", "n_cells")

write.csv(
  cluster_sizes,
  "results/tables/cluster_sizes.csv",
  row.names = FALSE
)


############################################################
# 7. 2D UMAP
############################################################

seu <- RunUMAP(
  seu,
  dims = dims_use
)

umap_cluster <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("UMAP by Seurat clusters") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = "results/figures/06_UMAP_clusters.png",
  plot = umap_cluster,
  width = 8,
  height = 6,
  dpi = 300
)


############################################################
# 8. Marker Gene Discovery
############################################################

markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  "results/tables/all_cluster_markers.csv",
  row.names = FALSE
)

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(
  top10_markers,
  "results/tables/top10_markers_per_cluster.csv",
  row.names = FALSE
)

top5_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

write.csv(
  top5_markers,
  "results/tables/top5_markers_per_cluster.csv",
  row.names = FALSE
)


############################################################
# 9. General Marker Visualization
############################################################

# Edit this list according to your tissue/dataset.
# Example below is for immune/PBMC-like datasets.

marker_genes <- c(
  "CD3D", "CD3E", "IL7R", "CCR7",
  "CD8A", "CD8B",
  "MS4A1", "CD79A",
  "LYZ", "S100A8", "S100A9",
  "FCGR3A", "MS4A7",
  "GNLY", "NKG7", "GZMB",
  "FCER1A", "CST3",
  "PPBP", "PF4"
)

marker_genes <- marker_genes[marker_genes %in% rownames(seu)]

if (length(marker_genes) > 0) {
  
  dot_plot <- DotPlot(
    seu,
    features = marker_genes
  ) +
    RotatedAxis() +
    ggtitle("Marker gene expression by cluster")
  
  ggsave(
    filename = "results/figures/07_marker_dotplot.png",
    plot = dot_plot,
    width = 12,
    height = 6,
    dpi = 300
  )
  
  feature_plot <- FeaturePlot(
    seu,
    features = head(marker_genes, 6),
    ncol = 3
  )
  
  ggsave(
    filename = "results/figures/08_marker_featureplot.png",
    plot = feature_plot,
    width = 12,
    height = 8,
    dpi = 300
  )
}


############################################################
# 10. Manual Cell Type Annotation Template
############################################################

# At first, each cell type is named by cluster.
# After marker validation, manually edit this section.

seu$cluster_id <- as.character(Idents(seu))
seu$celltype_manual <- paste0("Cluster_", seu$cluster_id)

# Example:
# seu$celltype_manual[seu$seurat_clusters == "0"] <- "Naive CD4 T"
# seu$celltype_manual[seu$seurat_clusters == "1"] <- "CD14+ Monocytes"
# seu$celltype_manual[seu$seurat_clusters == "2"] <- "B cells"

umap_celltype <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "celltype_manual",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("UMAP by manual annotation") +
  theme_classic()

ggsave(
  filename = "results/figures/09_UMAP_manual_annotation.png",
  plot = umap_celltype,
  width = 9,
  height = 6,
  dpi = 300
)


############################################################
# 11. Heatmap of Top Marker Genes
############################################################

heatmap_genes <- unique(top5_markers$gene)
heatmap_genes <- heatmap_genes[heatmap_genes %in% rownames(seu)]

if (length(heatmap_genes) > 0) {
  
  marker_heatmap <- DoHeatmap(
    seu,
    features = heatmap_genes
  ) +
    NoLegend() +
    ggtitle("Top marker genes per cluster")
  
  ggsave(
    filename = "results/figures/10_top_marker_heatmap.png",
    plot = marker_heatmap,
    width = 12,
    height = 10,
    dpi = 300
  )
}


############################################################
# 12. 3D UMAP Interactive Plot
############################################################

seu <- RunUMAP(
  seu,
  dims = dims_use,
  n.components = 3,
  reduction.name = "umap3d",
  reduction.key = "UMAP3D_"
)

umap3d <- Embeddings(seu, "umap3d")

plot_df <- data.frame(
  UMAP_1 = umap3d[, 1],
  UMAP_2 = umap3d[, 2],
  UMAP_3 = umap3d[, 3],
  cluster = as.character(seu$seurat_clusters),
  celltype = as.character(seu$celltype_manual),
  nFeature_RNA = seu$nFeature_RNA,
  nCount_RNA = seu$nCount_RNA,
  percent.mt = seu$percent.mt
)

p3d_umap <- plot_ly(
  plot_df,
  x = ~UMAP_1,
  y = ~UMAP_2,
  z = ~UMAP_3,
  color = ~cluster,
  colors = "Set3",
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2, opacity = 0.75),
  text = ~paste(
    "Cluster:", cluster,
    "<br>Cell type:", celltype,
    "<br>Genes:", nFeature_RNA,
    "<br>UMIs:", nCount_RNA,
    "<br>MT%:", round(percent.mt, 2)
  ),
  hoverinfo = "text"
) %>%
  layout(
    title = "3D UMAP colored by cluster",
    scene = list(
      xaxis = list(title = "UMAP 1"),
      yaxis = list(title = "UMAP 2"),
      zaxis = list(title = "UMAP 3")
    )
  )

htmlwidgets::saveWidget(
  p3d_umap,
  "results/interactive/3D_UMAP_clusters.html",
  selfcontained = TRUE
)


############################################################
# 13. 3D PCA Interactive Plot
############################################################

pca_coords <- Embeddings(seu, "pca")

pca_df <- data.frame(
  PC_1 = pca_coords[, 1],
  PC_2 = pca_coords[, 2],
  PC_3 = pca_coords[, 3],
  cluster = as.character(seu$seurat_clusters),
  celltype = as.character(seu$celltype_manual),
  nFeature_RNA = seu$nFeature_RNA,
  nCount_RNA = seu$nCount_RNA,
  percent.mt = seu$percent.mt
)

p3d_pca <- plot_ly(
  pca_df,
  x = ~PC_1,
  y = ~PC_2,
  z = ~PC_3,
  color = ~cluster,
  colors = "Set3",
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2, opacity = 0.75),
  text = ~paste(
    "Cluster:", cluster,
    "<br>Cell type:", celltype,
    "<br>Genes:", nFeature_RNA,
    "<br>UMIs:", nCount_RNA,
    "<br>MT%:", round(percent.mt, 2)
  ),
  hoverinfo = "text"
) %>%
  layout(
    title = "3D PCA colored by cluster",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

htmlwidgets::saveWidget(
  p3d_pca,
  "results/interactive/3D_PCA_clusters.html",
  selfcontained = TRUE
)


############################################################
# 14. Save Processed Object and Session Info
############################################################

saveRDS(
  seu,
  "results/processed_seurat_object.rds"
)

sink("results/sessionInfo.txt")
sessionInfo()
sink()

cat("\nPipeline completed successfully!\n")
cat("All results saved in the 'results/' folder.\n")