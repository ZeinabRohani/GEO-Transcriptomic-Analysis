############################################################
# Project: GEO-Transcriptomic-Analysis
# Script : 02_microarray_preprocessing.R
#
# Description:
# Preprocess GEO microarray expression data:
# - Check expression distribution
# - Log2 transformation (if required)
# - Quantile normalization
# - Save normalized expression matrix
############################################################

## Load package
library(limma)

## Load expression matrix
expr <- read.delim(
  "data/expression_matrix.txt",
  row.names = 1,
  check.names = FALSE
)

############################################################
# Check expression distribution
############################################################

qx <- quantile(
  expr,
  probs = c(0, 0.25, 0.50, 0.75, 0.99, 1),
  na.rm = TRUE
)

need_log2 <-
  (qx[5] > 100) ||
  ((qx[6] - qx[1]) > 50 && qx[2] > 0)

############################################################
# Log2 transformation (if necessary)
############################################################

if (need_log2) {

  expr[expr <= 0] <- NA

  expr <- log2(expr)

}

############################################################
# Remove probes with missing values
############################################################

expr <- expr[
  complete.cases(expr),
]

############################################################
# Quantile normalization
############################################################

expr_norm <- normalizeBetweenArrays(
  expr,
  method = "quantile"
)

############################################################
# Save normalized matrix
############################################################

write.table(
  expr_norm,
  file = "data/expression_normalized.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

############################################################
# Boxplot before normalization
############################################################

pdf(
  "results/boxplot_before_normalization.pdf",
  width = 10,
  height = 6
)

boxplot(
  expr,
  outline = FALSE,
  las = 2,
  main = "Before Normalization",
  ylab = "Expression"
)

dev.off()

############################################################
# Boxplot after normalization
############################################################

pdf(
  "results/boxplot_after_normalization.pdf",
  width = 10,
  height = 6
)

boxplot(
  expr_norm,
  outline = FALSE,
  las = 2,
  main = "After Quantile Normalization",
  ylab = "Normalized Expression"
)

dev.off()

############################################################
# Summary
############################################################

cat("\n")
cat("=====================================\n")
cat("MICROARRAY PREPROCESSING COMPLETED\n")
cat("=====================================\n")
cat("Genes   :", nrow(expr_norm), "\n")
cat("Samples :", ncol(expr_norm), "\n")
cat("=====================================\n")
