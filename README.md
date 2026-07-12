# GEO-Transcriptomic-Analysis

## Overview

This repository presents reproducible transcriptomic analysis workflows implemented in **R** and **Bioconductor** for the analysis of publicly available **Gene Expression Omnibus (GEO)** datasets. The project demonstrates complete pipelines for gene expression analysis, including data retrieval, preprocessing, quality assessment, normalization, differential gene expression analysis, and publication-quality visualization.

The repository was developed as part of a bioinformatics portfolio to demonstrate practical experience in transcriptomic data analysis using widely adopted Bioconductor packages and reproducible computational workflows.

---

## Project Objectives

The primary objectives of this repository are to:

* Retrieve publicly available transcriptomic datasets from GEO.
* Perform preprocessing and quality assessment of gene expression data.
* Normalize expression values using appropriate statistical methods.
* Identify differentially expressed genes (DEGs).
* Generate publication-quality visualizations.
* Produce reproducible analysis pipelines suitable for downstream biological interpretation.

---

## Repository Structure

```
GEO-Transcriptomic-Analysis/
│
├── data/
│   ├── raw/
│   └── processed/
│
├── scripts/
│   ├── microarray_analysis.R
│   └── rnaseq_analysis.R
│
├── results/
│   ├── figures/
│   ├── tables/
│   └── DEG_results/
│
├── README.md
└── LICENSE
```

---

## Analysis Workflow

### Microarray Analysis

The microarray workflow includes:

* GEO dataset retrieval using **GEOquery**
* Expression matrix extraction
* Data preprocessing
* Log2 transformation (when required)
* Experimental design matrix construction
* Linear modeling using **limma**
* Empirical Bayes statistical analysis
* Differentially expressed gene (DEG) identification
* Probe annotation
* Export of analysis results

---

### RNA-seq Analysis

The RNA-seq workflow includes:

* Import of raw count matrices
* Low-expression filtering using **edgeR**
* Library size normalization using **TMM**
* Experimental design matrix construction
* Mean-variance modeling using **voom**
* Linear modeling using **limma**
* Empirical Bayes statistics
* Differential gene expression analysis
* Export of normalized expression values and DEG tables

---

## Data Visualization

The repository includes scripts for generating publication-quality figures, including:

* Volcano plots (**EnhancedVolcano**)
* Heatmaps (**pheatmap**)
* Boxplots for quality assessment
* Expression distribution plots
* Differential expression summaries

---

## Software and Packages

### Programming Language

* R

### Bioconductor Packages

* GEOquery
* limma
* edgeR
* affy
* Biobase
* EnhancedVolcano

### CRAN Packages

* ggplot2
* pheatmap

---

## Workflow Summary

1. Download GEO datasets
2. Import expression data
3. Perform quality assessment
4. Normalize expression values
5. Construct the experimental design matrix
6. Fit linear models
7. Identify differentially expressed genes
8. Generate visualizations
9. Export analysis results

---

## Reproducibility

The workflows are organized as modular R scripts following reproducible research principles. File paths are project-relative, enabling execution across different operating systems with minimal modification. All statistical analyses rely on established Bioconductor methodologies widely used in transcriptomic research.

---

## Applications

These workflows can be adapted for a variety of transcriptomic studies, including:

* Differential gene expression analysis
* Comparative transcriptomics
* Biomarker discovery
* Candidate gene prioritization
* Functional genomics studies
* Exploratory analysis of publicly available transcriptomic datasets

---

## Requirements

The analyses were developed using:

* R
* Bioconductor
* RStudio (recommended)

Install required packages before running the scripts:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "GEOquery",
    "limma",
    "edgeR",
    "affy",
    "Biobase",
    "EnhancedVolcano"
))

install.packages(c(
    "ggplot2",
    "pheatmap"
))
```

---

## Citation

If you use these workflows in your research, please cite:

* The original GEO datasets used in the analyses.
* The corresponding Bioconductor packages.
* R and Bioconductor.

---

## Author

**Zeinab Rohani**

Bioinformatics | Computational Biology | Transcriptomics | Functional Genomics

GitHub: https://github.com/ZeinabRohani

---

## License

This repository is released for research and educational purposes.

---

## Acknowledgements

This project makes use of publicly available transcriptomic datasets deposited in the NCBI Gene Expression Omnibus (GEO) and relies on the Bioconductor ecosystem for reproducible statistical analysis of high-throughput genomic data.
