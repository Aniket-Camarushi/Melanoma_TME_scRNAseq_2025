# Single-Cell RNA-Seq Analysis of Melanoma Tumor Microenvironment  

**Goal**: This repository contains a complete reanalysis of the GSE215121 dataset using R and Seurat, reproducing the findings of [Zhang et al. (2022)](https://doi.org/10.1038/s41467-022-34877-3). 

## Key Features  
- **Data**: [GSE215121](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215121)  
- **Tools**: Seurat, CellChat, Harmony.
- **Key Findings**:  
  - Identified exhausted CD8+ T-cell subsets in PD-1-resistant tumors.  
  - Predicted tumor-immune crosstalk via ligand-receptor networks.  

## Project Structure
- `scripts/01_data_processing.R`: Download and preprocess scRNA-seq data
- `scripts/02_cell_annotation.R`: Cell type identification and clustering
- `scripts/03_tumor_heterogeneity.R`: Analysis of melanoma subpopulations
- `scripts/04_cell_communication.R`: CellChat analysis of interactions
- `scripts/05_differential_analysis.R`: Comparison of acral vs cutaneous melanoma

## Requirements  
- **R Version**: 4.5.0
- Seurat >= 5.2.1
- **Session Info**: [sessionInfo.txt](sessionInfo.txt)  

## Quick Start  
```bash
git clone https://github.com/yourname/melanoma-tme-scrnaseq
Rscript scripts/01_data_processing.R
```

## Results
Key findings reproduced from the original paper:
1. Identification of functional melanoma subpopulations
