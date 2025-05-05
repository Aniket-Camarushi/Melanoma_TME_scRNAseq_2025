# Single-Cell RNA-Seq Analysis of Melanoma Tumor Microenvironment  

**Goal**: This repository contains a complete reanalysis of the GSE215121 dataset using R and Seurat, reproducing the findings of [Zhang et al. (2022)](https://doi.org/10.1038/s41467-022-34877-3).
Replicate the pipeline from (Zhang, 2022) to analyze immune cell interactions in melanoma.  

## Key Features  
- **Data**: [GSE123814](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123814) (10x Genomics, 31 tumors).  
- **Tools**: Seurat, CellChat, Harmony.  
- **Key Findings**:  
  - Identified exhausted CD8+ T-cell subsets in PD-1-resistant tumors.  
  - Predicted tumor-immune crosstalk via ligand-receptor networks.  

## Reproducibility  
- **R Version**: 4.3.0  
- **Session Info**: [sessionInfo.txt](sessionInfo.txt)  
- **Conda Environment**: [environment.yml](environment.yml)  

## Quick Start  
```bash
git clone https://github.com/yourname/melanoma-tme-scrnaseq  
conda env create -f environment.yml  
Rscript code/01_qc_normalization.R
```
