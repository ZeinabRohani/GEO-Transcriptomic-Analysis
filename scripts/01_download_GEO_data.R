############################################################
# Project: GEO-Transcriptomic-Analysis
# Script : 01_download_GEO_data.R
#
# Description:
# Download a GEO dataset and extract:
# - Expression matrix
# - Phenotype data
# - Feature annotation
############################################################

## Load packages
library(GEOquery)
library(Biobase)

## GEO accession number
GEO_ID <- "GSE50040"

## Create directories
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

## Download GEO dataset
gset <- getGEO(
  GEO = GEO_ID,
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

## Select the first ExpressionSet
gset <- gset[[1]]

## Extract data
expression_matrix <- exprs(gset)
phenotype_data <- pData(gset)
feature_annotation <- fData(gset)

## Save expression matrix
write.table(
  expression_matrix,
  file = "data/expression_matrix.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

## Save phenotype data
write.table(
  phenotype_data,
  file = "data/phenotype_data.txt",
  sep = "\t",
  quote = FALSE
)

## Save feature annotation
write.table(
  feature_annotation,
  file = "data/feature_annotation.txt",
  sep = "\t",
  quote = FALSE
)

## Dataset summary
cat("\n")
cat("=====================================\n")
cat("GEO DATA DOWNLOAD COMPLETED\n")
cat("=====================================\n")
cat("Dataset :", GEO_ID, "\n")
cat("Genes   :", nrow(expression_matrix), "\n")
cat("Samples :", ncol(expression_matrix), "\n")
cat("=====================================\n")
