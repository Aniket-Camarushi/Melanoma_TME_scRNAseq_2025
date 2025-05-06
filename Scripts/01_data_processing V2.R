# 01_data_processing.R
# Part 1: Initial Setup ---------------------------------------------------
setwd("C:/Users/cmbla/Bioinformatics Scripts/Recreate BINF Papers/Comprehensive scRNA-seq Analysis of Acral Melanoma/")

# Create directories using safer path handling
dirs_to_create <- c("Data", "Scripts", "Results", "Plots")
sapply(dirs_to_create, function(x) {
  dir.create(file.path(getwd(), x), showWarnings = FALSE)
})

# Load packages
setwd("./Scripts/")
source("Setup.R")
setup()

# Part 2: Data Download ---------------------------------------------------
geo_accession <- "GSE215120"
data_dir <- file.path("Data", "Melanoma_Cells_10X")

if (!dir.exists(data_dir)) {
  # Create more robust download function
  download_geo_data <- function(geo_id, dest_dir) {
    tryCatch({
      archive_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", 
                            geo_id, "&format=file")
      dest_file <- file.path("../Data", "Melanoma_Cells.tar")
      
      # Platform-independent download
      if (Sys.info()['sysname'] == "Windows") {
        download.file(archive_url, destfile = dest_file, mode = "wb")
      } else {
        download.file(archive_url, destfile = dest_file, method = "curl")
      }
      
      untar(dest_file, exdir = dest_dir)
      unlink(dest_file) # Remove tar after extraction
    }, error = function(e) {
      message("Download failed: ", e$message)
    })
  }
  
  download_geo_data(geo_accession, data_dir)
}

# Logical Breakpoint 1: Save raw data object
saveRDS(list.files(data_dir, full.names = TRUE), 
        file = "../Data/raw_file_paths.rds")

# Part 3: Data Loading and Merging ----------------------------------------
h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)

# Add parallel processing for large datasets
plan("multisession", workers = 4)

seurat_list <- lapply(h5_files, function(h5_file) {
  tryCatch({
    # Validate file
    if (!file.exists(h5_file)) {
      stop("File not found: ", h5_file)
    }
    
    sample_name <- tools::file_path_sans_ext(basename(h5_file))
    data <- Read10X_h5(h5_file)
    
    # Add basic QC metrics during creation
    CreateSeuratObject(
      counts = data,
      project = sample_name,
      min.cells = 3,
      min.features = 200,
      meta.data = list(
        sample_id = sample_name,
        batch = ifelse(grepl("AC", sample_name), "Acral", "Cutaneous")
      )
    )
    }, error = function(e) {
        message("Error processing ", h5_file, ": ", e$message)
        NULL
    })
  }
)
  
# Remove any failed samples
seurat_list <- seurat_list[!sapply(seurat_list, is.null)]

# Merge with improved metadata handling
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sapply(seurat_list, function(x) unique(x$sample_id))
)

# Logical Breakpoint 2: Save merged object
saveRDS(merged_seurat, file = "../Data/merged_unfiltered_seurat.rds")

# Part 4: Quality Control -------------------------------------------------
# Enhanced mitochondrial gene detection
mt_genes <- grep("^MT-", rownames(merged_seurat), value = TRUE)
merged_seurat[["percent_mt"]] <- PercentageFeatureSet(
  merged_seurat, 
  features = mt_genes
)

# Visualize QC metrics with plot saving
qc_plots <- VlnPlot(
  merged_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
  pt.size = 0.1,
  group.by = "sample_id"
) + patchwork::plot_annotation(title = "Pre-filter QC Metrics")

ggsave(
  filename = "../Plots/pre_filter_qc.pdf",
  plot = qc_plots,
  width = 12,
  height = 8
)

# Filter with more conservative thresholds
filtered_seurat <- subset(
  merged_seurat,
  subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &  
    nCount_RNA > 1000 & nCount_RNA < 30000 &          
    percent_mt < 15
)

# Logical Breakpoint 3: Save filtered object
saveRDS(filtered_seurat, file = "../Data/filtered_seurat.rds")

# Part 5: Normalization and Feature Selection -----------------------------
# SCTransform for improved normalization
filtered_seurat <- SCTransform(
  filtered_seurat,
  vars.to.regress = c("percent_mt", "nCount_RNA"),
  verbose = FALSE
)

# Enhanced variable feature detection
filtered_seurat <- FindVariableFeatures(
  filtered_seurat,
  selection.method = "vst",
  nfeatures = 3000  # Increased for comprehensive analysis
)

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
saveRDS(filtered_seurat, file = "Data/processed_seurat.rds")

# Generate and save PCA plot
pca_plot <- DimPlot(
  object = filtered_seurat,
  reduction = "harmony",
  group.by = "sample_id"
) + ggtitle("Harmony Integrated PCA")

ggsave(
  filename = "Plots/harmony_integrated_pca.pdf",
  plot = pca_plot,
  width = 10,
  height = 8
)