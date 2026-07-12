###############################################################
# Project : GEO-Transcriptomic-Analysis
# Script  : 01_GEO_microarray_download_preprocessing.R
# Author  : Zeinab Rohani
#
# Description:
# Download GEO microarray dataset
# Extract expression matrix
# Extract phenotype information
# Extract annotation
# Perform quality assessment
# Apply log2 transformation if required
# Save processed files
#
###############################################################

##############################
# Load Required Packages
##############################

required_packages <- c(
  "GEOquery",
  "Biobase",
  "BiocGenerics",
  "limma",
  "affy"
)

for(pkg in required_packages){

  if(!requireNamespace(pkg, quietly = TRUE)){

    stop(
      paste(
        "Package",
        pkg,
        "is not installed."
      )
    )

  }

  library(pkg, character.only = TRUE)

}

##############################
# Create Project Directories
##############################

dir.create(
  "data",
  showWarnings = FALSE
)

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

##############################
# GEO Accession
##############################

geo_accession <- "GSE50040"

message("-----------------------------------")
message("Downloading GEO dataset...")
message("-----------------------------------")

##############################
# Download Dataset
##############################

gset <- getGEO(

  GEO = geo_accession,

  GSEMatrix = TRUE,

  AnnotGPL = TRUE

)

if(length(gset) > 1){

  gset <- gset[[1]]

}else{

  gset <- gset[[1]]

}

message("Dataset successfully downloaded.")

##############################
# Extract Expression Matrix
##############################

expression_matrix <- exprs(gset)

##############################
# Extract Annotation
##############################

annotation <- fData(gset)

##############################
# Extract Phenotype Data
##############################

phenotype <- pData(gset)

##############################
# Save Raw Files
##############################

write.table(

  expression_matrix,

  file = "data/expression_matrix_raw.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

write.table(

  annotation,

  file = "data/annotation.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

write.table(

  phenotype,

  file = "data/phenotype.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

message("Raw files saved.")

##############################
# Expression Matrix Summary
##############################

cat("\n")
cat("-----------------------------------\n")
cat("Expression Matrix Summary\n")
cat("-----------------------------------\n")

cat(
  "Genes :",
  nrow(expression_matrix),
  "\n"
)

cat(
  "Samples :",
  ncol(expression_matrix),
  "\n"
)

cat("-----------------------------------\n")

##############################
# Check Whether
# Log2 Transformation
# Is Required
##############################

qx <- quantile(

  expression_matrix,

  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.99,
    1.00
  ),

  na.rm = TRUE

)

LogC <- (

  qx[5] > 100 ||

  (qx[6]-qx[1] > 50 &
     qx[2] > 0)

)

if(LogC){

  message("Applying log2 transformation...")

  expression_matrix[
    expression_matrix <= 0
  ] <- NA

  expression_matrix <-
    log2(expression_matrix)

}else{

  message(
    "Data appear to be already log-transformed."
  )

}

##############################
# Remove Missing Values
##############################

na_count <- sum(is.na(expression_matrix))

message(
  paste(
    "Missing values:",
    na_count
  )
)

##############################
# Quality Assessment
##############################

message("-----------------------------------")
message("Generating quality assessment plots...")
message("-----------------------------------")

pdf(
  file = "results/figures/Boxplot_Normalized_Expression.pdf",
  width = 10,
  height = 6
)

boxplot(
  expression_matrix,
  outline = FALSE,
  las = 2,
  col = "lightblue",
  main = paste(
    geo_accession,
    "- Normalized Expression"
  ),
  xlab = "Samples",
  ylab = "Log2 Expression"
)

dev.off()

##############################
# Save Normalized Matrix
##############################

write.table(

  expression_matrix,

  file = "data/expression_matrix_normalized.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Sample Information
##############################

sample_information <- data.frame(

  Sample = colnames(expression_matrix),

  stringsAsFactors = FALSE

)

write.table(

  sample_information,

  file = "results/tables/sample_information.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Dataset Summary
##############################

dataset_summary <- data.frame(

  GEO_Accession = geo_accession,

  Number_of_Genes = nrow(expression_matrix),

  Number_of_Samples = ncol(expression_matrix),

  Missing_Values = sum(is.na(expression_matrix)),

  stringsAsFactors = FALSE

)

write.table(

  dataset_summary,

  file = "results/tables/dataset_summary.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Save Session Information
##############################

writeLines(

  capture.output(sessionInfo()),

  con = "results/sessionInfo.txt"

)

##############################
# Final Messages
##############################

cat("\n")
cat("=============================================\n")
cat("Microarray preprocessing completed successfully.\n")
cat("=============================================\n")
cat("Dataset :", geo_accession, "\n")
cat("Genes   :", nrow(expression_matrix), "\n")
cat("Samples :", ncol(expression_matrix), "\n")
cat("=============================================\n")
cat("Output files\n")
cat("---------------------------------------------\n")
cat("data/expression_matrix_raw.txt\n")
cat("data/expression_matrix_normalized.txt\n")
cat("data/annotation.txt\n")
cat("data/phenotype.txt\n")
cat("results/figures/Boxplot_Normalized_Expression.pdf\n")
cat("results/tables/sample_information.txt\n")
cat("results/tables/dataset_summary.txt\n")
cat("results/sessionInfo.txt\n")
cat("=============================================\n")

###############################################################
# End of Script
###############################################################
