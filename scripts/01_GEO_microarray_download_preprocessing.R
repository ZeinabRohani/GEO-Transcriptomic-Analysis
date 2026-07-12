###############################################################
# Project: GEO-Transcriptomic-Analysis
#
# Script:
# 01_GEO_microarray_download_preprocessing.R
#
# Author:
# Zeinab Rohani
#
# Description:
# Download a GEO microarray dataset, extract expression matrix,
# phenotype information and annotation, perform preprocessing,
# and export processed data.
#
###############################################################

############################
# Load Required Packages
############################

library(GEOquery)
library(limma)
library(affy)
library(Biobase)
library(BiocGenerics)

############################
# Create Output Directories
############################

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

############################
# Download GEO Dataset
############################

geo_id <- "GSE50040"

gset <- getGEO(
  GEO = geo_id,
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

gset <- gset[[1]]

############################
# Extract Data
############################

expression_matrix <- exprs(gset)

annotation <- fData(gset)

phenotype <- pData(gset)

############################
# Save Raw Data
############################

write.table(
  expression_matrix,
  file = "data/expression_matrix.txt",
  sep = "\t",
  quote = FALSE
)

write.table(
  annotation,
  file = "data/annotation.txt",
  sep = "\t",
  quote = FALSE
)

write.table(
  phenotype,
  file = "data/phenotype.txt",
  sep = "\t",
  quote = FALSE
)

############################
# Log2 Transformation
############################

expression_matrix <- log2(expression_matrix + 1)

############################
# Quality Assessment
############################

pdf(
  "results/figures/boxplot_normalized_expression.pdf",
  width = 8,
  height = 6
)

boxplot(
  expression_matrix,
  las = 2,
  col = "lightblue",
  main = "Normalized Expression Values"
)

dev.off()

############################
# Save Normalized Matrix
############################

write.table(
  expression_matrix,
  file = "data/normalized_expression_matrix.txt",
  sep = "\t",
  quote = FALSE
)

############################
# Session Information
############################

writeLines(
  capture.output(sessionInfo()),
  "results/sessionInfo.txt"
)

###############################################################
# End of Script
###############################################################
