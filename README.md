# GEO-Transcriptomic-Analysis

## Overview

This repository provides reproducible R workflows for transcriptomic analysis of publicly available Gene Expression Omnibus (GEO) datasets. 
The project includes pipelines for both microarray and RNA-seq data, covering data acquisition, preprocessing, normalization, differential gene expression analysis, and visualization using Bioconductor packages.
The workflows were developed to support reproducible transcriptomic studies and can be adapted for the analysis of other publicly available gene expression datasets.

## Objectives

* Download and organize transcriptomic datasets from GEO.
* Perform preprocessing and quality assessment of expression data.
* Conduct differential gene expression analysis for both microarray and RNA-seq datasets.
* Generate publication-quality visualizations.
* Export reproducible analysis results for downstream biological interpretation.

## Methods

### Microarray Analysis

* Data retrieval using **GEOquery**
* Expression matrix extraction
* Log2 transformation
* Differential expression analysis using **limma**
* Probe annotation
* Export of differentially expressed genes (DEGs)

### RNA-seq Analysis

* Import of raw count matrix
* Low-expression filtering using **edgeR**
* Library size normalization using **TMM**
* Mean-variance modeling using **voom**
* Linear modeling and empirical Bayes statistics using **limma**
* Identification of differentially expressed genes

## Data Visualization

The repository includes scripts for generating publication-quality figures, including:

* Volcano plots using **EnhancedVolcano**
* Heatmaps using **pheatmap**
* Boxplots for quality assessment

## Software and Packages

* R
* Bioconductor
* GEOquery
* limma
* edgeR
* affy
* Biobase
* EnhancedVolcano
* pheatmap
* ggplot2

## Workflow

1. Download GEO datasets
2. Import expression data
3. Perform quality assessment
4. Normalize expression values
5. Construct experimental design matrix
6. Identify differentially expressed genes
7. Visualize results
8. Export analysis outputs

## Reproducibility

The workflows are organized as modular R scripts to facilitate reproducible analyses. 
File paths are intended to be project-relative, allowing the repository to be executed on different operating systems with minimal modification.

## Applications

These workflows can be adapted for:

* Differential gene expression analysis
* Comparative transcriptomics
* Biomarker discovery
* Candidate gene prioritization
* Exploratory analysis of publicly available transcriptomic datasets

## Citation

If you use this workflow in your research, please cite the original GEO datasets and the corresponding Bioconductor packages used in the analysis.

## Author

**Zeinab Rohani**

Research interests:

* Bioinformatics
* Computational Biology
* Transcriptomics
* Functional Genomics

**This repository is intended for research and educational purposes and demonstrates reproducible transcriptomic analysis workflows implemented in R using Bioconductor.**
