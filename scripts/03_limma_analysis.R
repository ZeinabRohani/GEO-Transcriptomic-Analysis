############################################################
# Project: GEO-Transcriptomic-Analysis
# Script : 03_limma_analysis.R
#
# Description:
# Differential gene expression analysis using limma
############################################################

## Load package
library(limma)

############################################################
# Load expression matrix
############################################################

expr <- read.delim(
  "data/expression_normalized.txt",
  row.names = 1,
  check.names = FALSE
)

############################################################
# Define sample groups
############################################################
# Modify according to your GEO dataset

group <- factor(c(
  rep("Control", 37),
  rep("Tumor", 43)
))

############################################################
# Design matrix
############################################################

design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

############################################################
# Linear model fitting
############################################################

fit <- lmFit(expr, design)

############################################################
# Define contrast
############################################################

contrast.matrix <- makeContrasts(

  Tumor - Control,

  levels = design

)

############################################################
# Empirical Bayes statistics
############################################################

fit2 <- contrasts.fit(
  fit,
  contrast.matrix
)

fit2 <- eBayes(fit2)

############################################################
# Differentially expressed genes
############################################################

deg <- topTable(

  fit2,

  adjust.method = "BH",

  number = Inf,

  sort.by = "P"

)

############################################################
# Significant genes
############################################################

deg_sig <- subset(

  deg,

  adj.P.Val < 0.05 &
  abs(logFC) > 1

)

############################################################
# Save results
############################################################

write.table(

  deg,

  file = "results/limma_all_genes.txt",

  sep = "\t",

  quote = FALSE

)

write.table(

  deg_sig,

  file = "results/limma_significant_genes.txt",

  sep = "\t",

  quote = FALSE

)

############################################################
# Summary
############################################################

cat("\n")
cat("=====================================\n")
cat("LIMMA ANALYSIS COMPLETED\n")
cat("=====================================\n")
cat("Total genes       :", nrow(deg), "\n")
cat("Significant genes :", nrow(deg_sig), "\n")
cat("=====================================\n")
