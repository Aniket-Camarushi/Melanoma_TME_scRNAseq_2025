# Single-Cell RNA-Seq Analysis of Melanoma Tumor Microenvironment  

**Goal**: Replicate the pipeline from [Paper 1](https://doi.org/xxx) (JECCR, 2021) to analyze immune cell interactions and PD-1 resistance mechanisms in melanoma.  

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
