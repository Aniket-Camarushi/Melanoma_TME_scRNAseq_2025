  
setup <- function()
{
  # Install required packages
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  bioc_pkgs <- c("Seurat", "SingleCellExperiment", "scater", "scran", 
                         "harmony", "monocle3", "celldex", "SingleR", 
                         "infercnv", "CellChat", "GSVA", "clusterProfiler",
                         "org.Hs.eg.db", "GSEABase", "AUCell", "SCENIC")
  
  bioc2install <- bioc_pkgs[! bioc_pkgs %in% rownames(installed.packages())]
  
  if (length(bioc2install) > 0) 
  {
    BiocManager::install(bioc2install)
  }
  
  pkgs_list <- c("tidyverse", "ggplot2", "patchwork", "viridis", 
                     "ggridges", "ggrepel", "RColorBrewer", "devtools",
                     "remotes", "Matrix", "data.table")
  
  install.packages(pkgs_list[! pkgs_list %in% rownames(installed.packages())])
  
  # Install from GitHub if needed
  # remotes::install_github("satijalab/seurat-wrappers")
  # remotes::install_github("immunogenomics/harmony")
  # remotes::install_github("cole-trapnell-lab/monocle3")
  devtools::install_github('immunogenomics/presto')
  
  library(GEOquery)
  library(Seurat)
  library(monocle3)
  library(CellChat)
  library(harmony)
  library(SingleCellExperiment)
  library(SingleR)
  library(celldex)
  library(scater)
  library(GSVA)
  library(ggplot2)
  library(devtools)
  library(RColorBrewer)
  library(Matrix)
  library(tidyverse)
  library(rafalib)
  library(glue)
  library(dplyr)
  library(future)
  library(future.apply)
}