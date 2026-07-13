############################################################
# GEO-Transcriptomic-Analysis
# Script: 01_download_GEO_data.R
#
# Download GEO dataset and extract expression,
# phenotype and annotation data.
############################################################

library(GEOquery)
library(Biobase)

## Dataset
GEO_ID <- "GSE50040"

## Download dataset
gset <- getGEO(
  GEO = GEO_ID,
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

gset <- gset[[1]]

## Extract data
expression_matrix <- exprs(gset)
phenotype_data <- pData(gset)
feature_annotation <- fData(gset)

## Create folders
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

## Save files
write.table(
  expression_matrix,
  "data/expression_matrix.txt",
  sep = "\t",
  quote = FALSE
)

write.table(
  phenotype_data,
  "data/phenotype_data.txt",
  sep = "\t",
  quote = FALSE
)

write.table(
  feature_annotation,
  "data/feature_annotation.txt",
  sep = "\t",
  quote = FALSE
)

## Summary
cat("\n")
cat("Dataset:", GEO_ID, "\n")
cat("Genes:", nrow(expression_matrix), "\n")
cat("Samples:", ncol(expression_matrix), "\n")
cat("Download completed.\n")
