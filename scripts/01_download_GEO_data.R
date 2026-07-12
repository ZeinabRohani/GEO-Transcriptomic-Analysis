###############################################################
# GEO-Transcriptomic-Analysis
#
# Script : 01_download_GEO_data.R
# Author : Zeinab Rohani
#
# Description:
# Download a GEO dataset and extract:
#   - Expression matrix
#   - Phenotype data
#   - Feature annotation
#
# Compatible with:
#   - GEO Microarray datasets
#
###############################################################

##############################
# Load Required Packages
##############################

required_packages <- c(
  "GEOquery",
  "Biobase"
)

for(pkg in required_packages){

  if(!requireNamespace(pkg, quietly = TRUE)){

    stop(
      paste0(
        "Package '",
        pkg,
        "' is not installed.\n",
        "Please install it before running this script."
      )
    )

  }

  library(pkg, character.only = TRUE)

}

##############################
# User Settings
##############################

GEO_ID <- "GSE50040"

##############################
# Create Project Directories
##############################

directories <- c(

  "data",

  "data/raw",

  "data/processed",

  "results",

  "results/tables",

  "results/figures"

)

for(dir_name in directories){

  if(!dir.exists(dir_name)){

    dir.create(
      dir_name,
      recursive = TRUE
    )

  }

}

##############################
# Download GEO Dataset
##############################

message("---------------------------------------")
message("Downloading GEO dataset...")
message("---------------------------------------")

gset <- getGEO(

  GEO = GEO_ID,

  GSEMatrix = TRUE,

  AnnotGPL = TRUE

)

##############################
# Select ExpressionSet
##############################

if(length(gset) == 0){

  stop("No ExpressionSet was returned.")

}

if(length(gset) > 1){

  message(
    "Multiple platforms detected."
  )

  message(
    "Using the first ExpressionSet."
  )

}

gset <- gset[[1]]

message("Download completed.")

##############################
# Extract Data
##############################

expression_matrix <- exprs(gset)

phenotype_data <- pData(gset)

feature_annotation <- fData(gset)

##############################
# Basic Dataset Summary
##############################

cat("\n")

cat("=====================================\n")

cat("Dataset Summary\n")

cat("=====================================\n")

cat(
  "GEO accession : ",
  GEO_ID,
  "\n",
  sep = ""
)

cat(
  "Genes         : ",
  nrow(expression_matrix),
  "\n",
  sep = ""
)

cat(
  "Samples       : ",
  ncol(expression_matrix),
  "\n",
  sep = ""
)

cat("=====================================\n")

##############################
# Save Raw Expression Matrix
##############################

write.table(

  expression_matrix,

  file = "data/raw/expression_matrix.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Save Phenotype Data
##############################

write.table(

  phenotype_data,

  file = "data/raw/phenotype_data.txt",

  sep = "\t",

  quote = FALSE,

  row.names = TRUE

)

##############################
# Save Feature Annotation
##############################

write.table(

  feature_annotation,

  file = "data/raw/feature_annotation.txt",

  sep = "\t",

  quote = FALSE,

  row.names = TRUE

)

message("---------------------------------------")
message("Raw files successfully saved.")
message("---------------------------------------")

###############################################################
# Data Quality Assessment and Preprocessing
###############################################################

##############################
# Expression Value Summary
##############################

message("---------------------------------------")
message("Summarizing expression values...")
message("---------------------------------------")

summary(expression_matrix)

##############################
# Expression Quantiles
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

print(qx)

##############################
# Determine Whether Log2
# Transformation is Required
##############################

need_log2 <-

  (qx[5] > 100) ||

  ((qx[6] - qx[1]) > 50 && qx[2] > 0)

##############################
# Boxplot Before Processing
##############################

pdf(

  "results/figures/boxplot_raw_expression.pdf",

  width = 12,

  height = 6

)

boxplot(

  expression_matrix,

  outline = FALSE,

  las = 2,

  main = paste(
    GEO_ID,
    "Raw Expression"
  ),

  ylab = "Expression"

)

dev.off()

##############################
# Apply Log2 Transformation
##############################

if(need_log2){

  message("---------------------------------------")
  message("Applying log2 transformation...")
  message("---------------------------------------")

  expression_matrix[
    expression_matrix <= 0
  ] <- NA

  expression_matrix <-

    log2(expression_matrix)

}else{

  message("---------------------------------------")
  message("Dataset already appears to be log2 transformed.")
  message("---------------------------------------")

}

##############################
# Remove Genes with Missing Values
##############################

na_before <- sum(is.na(expression_matrix))

message(

  paste(

    "Missing values before filtering:",

    na_before

  )

)

expression_matrix <-

  expression_matrix[
    complete.cases(expression_matrix),
  ]

na_after <- sum(is.na(expression_matrix))

message(

  paste(

    "Missing values after filtering:",

    na_after

  )

)

##############################
# Boxplot After Processing
##############################

pdf(

  "results/figures/boxplot_processed_expression.pdf",

  width = 12,

  height = 6

)

boxplot(

  expression_matrix,

  outline = FALSE,

  las = 2,

  main = paste(
    GEO_ID,
    "Processed Expression"
  ),

  ylab = "Log2 Expression"

)

dev.off()

##############################
# Save Processed Matrix
##############################

write.table(

  expression_matrix,

  file = "data/processed/expression_matrix_processed.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Save Dataset Information
##############################

dataset_information <- data.frame(

  GEO_ID = GEO_ID,

  Genes = nrow(expression_matrix),

  Samples = ncol(expression_matrix),

  Log2_Transformation = need_log2,

  stringsAsFactors = FALSE

)

write.table(

  dataset_information,

  file = "results/tables/dataset_information.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

message("---------------------------------------")
message("Preprocessing completed successfully.")
message("---------------------------------------")

##############################
# Create Analysis Report
##############################

analysis_report <- c(

  "==========================================",

  "GEO Transcriptomic Analysis",

  "==========================================",

  paste("GEO Accession :", GEO_ID),

  paste("Number of Genes :", nrow(expression_matrix)),

  paste("Number of Samples :", ncol(expression_matrix)),

  paste("Log2 Transformation Applied :", need_log2),

  paste("Analysis Date :", Sys.Date()),

  "",

  "Files Generated:",

  "-------------------------",

  "data/raw/expression_matrix.txt",

  "data/raw/phenotype_data.txt",

  "data/raw/feature_annotation.txt",

  "data/processed/expression_matrix_processed.txt",

  "results/figures/boxplot_raw_expression.pdf",

  "results/figures/boxplot_processed_expression.pdf",

  "results/tables/dataset_information.txt"

)

writeLines(

  analysis_report,

  con = "results/analysis_report.txt"

)

##############################
# Save Session Information
##############################

writeLines(

  capture.output(sessionInfo()),

  con = "results/sessionInfo.txt"

)

##############################
# Final Console Summary
##############################

cat("\n")

cat("==========================================\n")

cat(" GEO DOWNLOAD COMPLETED SUCCESSFULLY\n")

cat("==========================================\n")

cat("Dataset :", GEO_ID, "\n")

cat("Genes   :", nrow(expression_matrix), "\n")

cat("Samples :", ncol(expression_matrix), "\n")

cat("==========================================\n")

cat("Output Directories\n")

cat("------------------------------------------\n")

cat("data/raw/\n")

cat("data/processed/\n")

cat("results/figures/\n")

cat("results/tables/\n")

cat("==========================================\n")

cat("Generated Files\n")

cat("------------------------------------------\n")

cat("✓ expression_matrix.txt\n")

cat("✓ phenotype_data.txt\n")

cat("✓ feature_annotation.txt\n")

cat("✓ expression_matrix_processed.txt\n")

cat("✓ dataset_information.txt\n")

cat("✓ analysis_report.txt\n")

cat("✓ sessionInfo.txt\n")

cat("✓ boxplot_raw_expression.pdf\n")

cat("✓ boxplot_processed_expression.pdf\n")

cat("==========================================\n")

message("Script completed successfully.")

###############################################################
# End of Script
###############################################################
