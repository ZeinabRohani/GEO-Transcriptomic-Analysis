# GEO-Transcriptomic-Analysis

## Overview

This repository presents a reproducible bioinformatics workflow for the analysis of public microarray gene expression datasets obtained from the **NCBI Gene Expression Omnibus (GEO)** using **R** and **Bioconductor**.

The project demonstrates a complete transcriptomic analysis pipeline, including GEO data retrieval, microarray preprocessing, differential gene expression analysis using **limma**, functional enrichment analysis, and publication-quality visualization.

This repository was developed as part of my bioinformatics portfolio to demonstrate practical experience in transcriptomic data analysis using reproducible computational workflows and widely adopted Bioconductor packages.

---

# Related Publication

The workflow implemented in this repository is based on the transcriptomic analysis strategy used in the following peer-reviewed publication.

**Rohani Z., Sazegar H., Rahimi E.** (2024).

*Unlocking the potential of Escherichia coli K-12: A novel approach for malignancy reduction in colorectal cancer through gene expression modulation.*

**Gene**, 906, 148266.

**DOI:** https://doi.org/10.1016/j.gene.2024.148266

---

## Relationship to this Repository

This repository provides a generalized and reproducible implementation of the transcriptomic analysis workflow used in the study, including:

* Downloading public GEO microarray datasets
* Expression matrix extraction
* Microarray preprocessing
* Differential gene expression analysis using **limma**
* Functional enrichment analysis using **Gene Ontology (GO)** and **KEGG pathways**
* Publication-quality visualization
* Reproducible analysis workflows implemented in **R/Bioconductor**

This repository is intended as an educational and research resource demonstrating reproducible transcriptomic analysis workflows and is **not intended to reproduce every analysis performed in the publication**.

---

# Project Objectives

The primary objectives of this repository are to:

* Retrieve publicly available GEO microarray datasets.
* Perform preprocessing of gene expression data.
* Normalize expression values using standard statistical methods.
* Identify differentially expressed genes (DEGs).
* Perform functional enrichment analysis.
* Generate publication-quality visualizations.
* Provide a reproducible transcriptomic analysis workflow.

---

# Analysis Workflow

```text
GEO Dataset
      │
      ▼
Download GEO Data
      │
      ▼
Microarray Preprocessing
      │
      ▼
Differential Expression Analysis (limma)
      │
      ▼
Functional Enrichment Analysis
      │
      ▼
Visualization
```

---

# Workflow Description

## 1. Download GEO Data

**Script**

`01_download_GEO_data.R`

This script downloads a GEO microarray dataset using **GEOquery** and extracts:

* Expression matrix
* Phenotype data
* Feature annotation

The extracted files are saved for downstream analysis.

---

## 2. Microarray Preprocessing

**Script**

`02_microarray_preprocessing.R`

This step prepares the expression matrix for downstream analysis by:

* Assessing expression value distribution
* Applying log2 transformation (when required)
* Removing missing values
* Performing quantile normalization
* Generating quality control boxplots

---

## 3. Differential Expression Analysis

**Script**

`03_limma_analysis.R`

Differential gene expression analysis is performed using the **limma** package.

The workflow includes:

* Experimental design matrix construction
* Linear model fitting
* Contrast analysis
* Empirical Bayes moderation
* Benjamini–Hochberg multiple testing correction
* Identification of significantly differentially expressed genes

---

## 4. Functional Enrichment Analysis

**Script**

`04_functional_enrichment.R`

Biological interpretation of differentially expressed genes is performed using **clusterProfiler**.

Analyses include:

* Gene Ontology (GO) Biological Process enrichment
* KEGG pathway enrichment

---

## 5. Visualization

**Script**

`05_visualization.R`

Publication-quality figures are generated for interpretation of transcriptomic results.

The visualization workflow includes:

* Volcano Plot
* Heatmap
* GO Dotplot
* KEGG Dotplot

---

# R Packages

The analysis pipeline uses the following R packages:

## Bioconductor

* GEOquery
* Biobase
* limma
* clusterProfiler
* org.Hs.eg.db

## CRAN

* ggplot2
* pheatmap

---


# Applications

This workflow can be adapted for a wide range of transcriptomic studies, including:

* Differential gene expression analysis
* Comparative transcriptomics
* Biomarker discovery
* Functional genomics
* Gene Ontology enrichment analysis
* KEGG pathway analysis
* Exploratory analysis of public GEO microarray datasets

---

# Requirements

The analyses were developed using:

* R
* Bioconductor
* RStudio (recommended)

Install the required packages before running the scripts:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "GEOquery",
    "Biobase",
    "limma",
    "clusterProfiler",
    "org.Hs.eg.db"
))

install.packages(c(
    "ggplot2",
    "pheatmap"
))
```

---

# Reproducibility

The workflow is organized as modular R scripts following reproducible research principles. Raw datasets are downloaded directly from the NCBI Gene Expression Omnibus (GEO), and all downstream analyses can be reproduced by executing the scripts in sequential order.

---

# Citation

If you use this workflow in your research, please cite:

* The original GEO datasets.
* The corresponding Bioconductor packages.
* R and Bioconductor.

---

# Author

**Zeinab Rohani**

Bioinformatics | Computational Biology | Transcriptomics | Functional Genomics

GitHub: https://github.com/ZeinabRohani

---

# License

This project is distributed under the **MIT License**.

---

# Acknowledgements

This project makes use of publicly available transcriptomic datasets deposited in the **NCBI Gene Expression Omnibus (GEO)** and relies on the **Bioconductor** ecosystem for reproducible statistical analysis of high-throughput gene expression data.
