###############################################################
# GEO-Transcriptomic-Analysis
#
# Script : 02_microarray_preprocessing.R
# Author : Zeinab Rohani
#
# Description:
# Microarray preprocessing workflow
#   - Load expression matrix
#   - Quality assessment
#   - Detect log transformation
#   - Normalize if required
#   - Save processed expression matrix
#
###############################################################

##############################
# Load Required Packages
##############################

required_packages <- c(
  "limma"
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
# Create Output Directories
##############################

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

##############################
# Input Files
##############################

expression_file <-
  "data/raw/expression_matrix.txt"

##############################
# Check Input File
##############################

if(!file.exists(expression_file)){

  stop(
    paste(
      "File not found:",
      expression_file
    )
  )

}

##############################
# Read Expression Matrix
##############################

message("--------------------------------")
message("Loading expression matrix...")
message("--------------------------------")

expression_matrix <- read.delim(

  expression_file,

  row.names = 1,

  check.names = FALSE,

  stringsAsFactors = FALSE

)

expression_matrix <-
  as.matrix(expression_matrix)

mode(expression_matrix) <- "numeric"

##############################
# Dataset Summary
##############################

message("--------------------------------")
message("Expression Matrix Summary")
message("--------------------------------")

cat(
  "Genes   :",
  nrow(expression_matrix),
  "\n"
)

cat(
  "Samples :",
  ncol(expression_matrix),
  "\n"
)

##############################
# Remove Empty Rows
##############################

expression_matrix <-

  expression_matrix[
    rowSums(
      is.na(expression_matrix)
    ) != ncol(expression_matrix),
  ]

##############################
# Quality Assessment
##############################

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
    1
  ),

  na.rm = TRUE

)

print(qx)

##############################
# Determine Whether
# Log2 Transformation
# Is Required
##############################

need_log2 <-

(
  qx[5] > 100
)

||

(
  (qx[6]-qx[1]) > 50 &&
  qx[2] > 0
)

##############################
# Raw Boxplot
##############################

pdf(

  "results/figures/boxplot_before_preprocessing.pdf",

  width = 11,

  height = 6

)

boxplot(

  expression_matrix,

  outline = FALSE,

  las = 2,

  col = "lightgray",

  main = "Raw Expression",

  ylab = "Expression"

)

dev.off()

##############################
# Log2 Transformation
##############################

if(need_log2){

  message("--------------------------------")
  message("Applying log2 transformation...")
  message("--------------------------------")

  expression_matrix[
    expression_matrix <= 0
  ] <- NA

  expression_matrix <-
    log2(expression_matrix)

}else{

  message("--------------------------------")
  message("Expression data already appear to be log2 transformed.")
  message("--------------------------------")

}

##############################
# Remove Features with Missing Values
##############################

message("--------------------------------")
message("Removing missing values...")
message("--------------------------------")

expression_matrix <-

  expression_matrix[
    complete.cases(expression_matrix),
  ]

##############################
# Between-array Normalization
##############################

message("--------------------------------")
message("Performing quantile normalization...")
message("--------------------------------")

expression_matrix <-

  normalizeBetweenArrays(

    expression_matrix,

    method = "quantile"

  )

##############################
# Summary After Normalization
##############################

summary(expression_matrix)

##############################
# Boxplot After Normalization
##############################

pdf(

  "results/figures/boxplot_after_normalization.pdf",

  width = 11,

  height = 6

)

boxplot(

  expression_matrix,

  outline = FALSE,

  las = 2,

  col = "lightblue",

  main = "Normalized Expression",

  ylab = "Log2 Expression"

)

dev.off()

##############################
# Density Plot
##############################

pdf(

  "results/figures/density_plot.pdf",

  width = 10,

  height = 6

)

plotDensities(

  expression_matrix,

  main = "Expression Density"

)

dev.off()

##############################
# Mean-SD Plot
##############################

gene_mean <- rowMeans(expression_matrix)

gene_sd <- apply(
  expression_matrix,
  1,
  sd
)

pdf(

  "results/figures/mean_sd_plot.pdf",

  width = 7,

  height = 6

)

plot(

  gene_mean,

  gene_sd,

  pch = 16,

  cex = 0.5,

  xlab = "Mean Expression",

  ylab = "Standard Deviation",

  main = "Mean-SD Plot"

)

abline(
  h = median(gene_sd),
  lty = 2
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

message("--------------------------------")
message("Processed expression matrix saved.")
message("--------------------------------")

###############################################################
# Quality Metrics and Filtering
###############################################################

##############################
# Number of Genes
##############################

number_of_genes <- nrow(expression_matrix)

##############################
# Number of Samples
##############################

number_of_samples <- ncol(expression_matrix)

##############################
# Missing Values
##############################

missing_values <- sum(is.na(expression_matrix))

##############################
# Sample Means
##############################

sample_means <- colMeans(expression_matrix)

##############################
# Sample Standard Deviations
##############################

sample_sd <- apply(
  expression_matrix,
  2,
  sd
)

##############################
# Sample Medians
##############################

sample_medians <- apply(
  expression_matrix,
  2,
  median
)

##############################
# Sample Statistics
##############################

sample_statistics <- data.frame(

  Sample = colnames(expression_matrix),

  Mean = sample_means,

  Median = sample_medians,

  SD = sample_sd,

  stringsAsFactors = FALSE

)

write.table(

  sample_statistics,

  file = "results/tables/sample_statistics.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Gene Variance
##############################

gene_variance <- apply(

  expression_matrix,

  1,

  var

)

##############################
# Remove Zero Variance Genes
##############################

expression_matrix <-

  expression_matrix[
    gene_variance > 0,
  ]

##############################
# Save Filtered Matrix
##############################

write.table(

  expression_matrix,

  file = "data/processed/expression_matrix_filtered.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Histogram
##############################

pdf(

  "results/figures/expression_histogram.pdf",

  width = 8,

  height = 6

)

hist(

  expression_matrix,

  breaks = 100,

  main = "Expression Distribution",

  xlab = "Expression",

  col = "gray"

)

dev.off()

##############################
# Correlation Matrix
##############################

sample_cor <- cor(

  expression_matrix,

  method = "pearson"

)

write.table(

  sample_cor,

  file = "results/tables/sample_correlation_matrix.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Correlation Heatmap
##############################

pdf(

  "results/figures/sample_correlation_heatmap.pdf",

  width = 8,

  height = 8

)

heatmap(

  sample_cor,

  symm = TRUE,

  scale = "none"

)

dev.off()

##############################
# Principal Component Analysis
##############################

message("--------------------------------")
message("Running Principal Component Analysis...")
message("--------------------------------")

pca_result <- prcomp(

  t(expression_matrix),

  center = TRUE,

  scale. = TRUE

)

##############################
# PCA Plot
##############################

pdf(

  "results/figures/PCA_plot.pdf",

  width = 8,

  height = 6

)

plot(

  pca_result$x[,1],

  pca_result$x[,2],

  pch = 19,

  xlab = paste0(
    "PC1 (",
    round(
      100 * summary(pca_result)$importance[2,1],
      2
    ),
    "%)"
  ),

  ylab = paste0(
    "PC2 (",
    round(
      100 * summary(pca_result)$importance[2,2],
      2
    ),
    "%)"
  ),

  main = "Principal Component Analysis"

)

text(

  pca_result$x[,1],

  pca_result$x[,2],

  labels = colnames(expression_matrix),

  pos = 3,

  cex = 0.7

)

dev.off()

##############################
# Save PCA Coordinates
##############################

pca_coordinates <- data.frame(

  Sample = rownames(pca_result$x),

  PC1 = pca_result$x[,1],

  PC2 = pca_result$x[,2],

  PC3 = pca_result$x[,3],

  stringsAsFactors = FALSE

)

write.table(

  pca_coordinates,

  file = "results/tables/PCA_coordinates.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Save Processed Matrix
##############################

write.table(

  expression_matrix,

  file = "data/processed/expression_matrix_final.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Preprocessing Summary
##############################

summary_report <- data.frame(

  Total_Genes = nrow(expression_matrix),

  Total_Samples = ncol(expression_matrix),

  Missing_Values = missing_values,

  Log2_Transformation = need_log2,

  Normalization = "Quantile",

  Zero_Variance_Genes_Removed = sum(gene_variance == 0),

  stringsAsFactors = FALSE

)

write.table(

  summary_report,

  file = "results/tables/preprocessing_summary.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Save Session Information
##############################

writeLines(

  capture.output(sessionInfo()),

  con = "results/sessionInfo_preprocessing.txt"

)

##############################
# Final Console Output
##############################

cat("\n")

cat("=========================================\n")

cat(" MICROARRAY PREPROCESSING COMPLETED\n")

cat("=========================================\n")

cat("Genes   :", nrow(expression_matrix), "\n")

cat("Samples :", ncol(expression_matrix), "\n")

cat("=========================================\n")

cat("Output Files\n")

cat("-----------------------------------------\n")

cat("✓ expression_matrix_processed.txt\n")

cat("✓ expression_matrix_filtered.txt\n")

cat("✓ expression_matrix_final.txt\n")

cat("✓ preprocessing_summary.txt\n")

cat("✓ sample_statistics.txt\n")

cat("✓ sample_correlation_matrix.txt\n")

cat("✓ PCA_coordinates.txt\n")

cat("✓ boxplot_before_preprocessing.pdf\n")

cat("✓ boxplot_after_normalization.pdf\n")

cat("✓ density_plot.pdf\n")

cat("✓ expression_histogram.pdf\n")

cat("✓ sample_correlation_heatmap.pdf\n")

cat("✓ PCA_plot.pdf\n")

cat("=========================================\n")

message("Microarray preprocessing completed successfully.")

###############################################################
# End of Script
###############################################################
