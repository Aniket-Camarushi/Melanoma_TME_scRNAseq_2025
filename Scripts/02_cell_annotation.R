
# Enhanced variable feature detection
filtered_seurat <- FindVariableFeatures(
  filtered_seurat,
  selection.method = "vst",
  nfeatures = 3000  # Increased for comprehensive analysis
)

# Scale data
filtered_seurat <- ScaleData(filtered_seurat, features = rownames(filtered_seurat))

# Save variable features for downstream analysis
write.csv(
  HVFInfo(filtered_seurat),
  file = "../Results/highly_variable_features.csv"
)

# Part 6: Dimensionality Reduction ---------------------------------------
# Run PCA with more components
filtered_seurat <- RunPCA(
  filtered_seurat,
  npcs = 50,
  features = VariableFeatures(filtered_seurat)
)

# Harmony integration with parameters from the paper
filtered_seurat <- RunHarmony(
  filtered_seurat,
  group.by.vars = c("sample_id", "batch"),
  theta = c(2, 2),  # Default parameters
  lambda = c(1, 1),  # Default parameters
  max.iter.harmony = 20
)

# Final save point
saveRDS(filtered_seurat, file = "../Data/processed_seurat.rds")

# Generate and save PCA plot
pca_plot <- DimPlot(
  object = filtered_seurat,
  reduction = "harmony",
  group.by = "sample_id"
) + ggtitle("Harmony Integrated PCA")

ggsave(
  filename = "../Plots/harmony_integrated_pca.pdf",
  plot = pca_plot,
  width = 10,
  height = 8
)

# Harmony batch correction (since samples may have batch effects)
library(harmony)
filtered_seurat <- RunHarmony(filtered_seurat, group.by.vars = "orig.ident")

# UMAP and clustering
filtered_seurat <- RunUMAP(filtered_seurat, reduction = "harmony", dims = 1:30)
filtered_seurat <- FindNeighbors(filtered_seurat, reduction = "harmony", dims = 1:30)
filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.8)

# Visualization
DimPlot(filtered_seurat, reduction = "umap", label = TRUE)

# Cell type annotation using known markers (as in the paper)
markers <- list(
  T_cells = c("CD3D", "CD3E"),
  B_cells = c("MS4A1", "CD79A"),
  NK_cells = c("FGFBP2", "KLRD1"),
  Monocytes_Macrophages = c("LYZ", "CD68", "CD14"),
  Melanoma_cells = c("MLANA", "PMEL", "MITF", "DCT"),
  Endothelial = c("VWF", "PECAM1"),
  Fibroblasts = c("COL1A1", "COL3A1")
)

# Dot plot of marker genes
DotPlot(filtered_seurat, features = unlist(markers)) + RotatedAxis()

# Automated annotation with SingleR
library(SingleR)
library(celldex)

# Download reference dataset
ref <- HumanPrimaryCellAtlasData()

# Get single-cell expression matrix
sce <- as.SingleCellExperiment(filtered_seurat)
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

# Add annotations to Seurat object
filtered_seurat$celltype <- pred$labels

