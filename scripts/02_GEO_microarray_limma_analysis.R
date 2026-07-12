###############################################################
# Project: GEO-Transcriptomic-Analysis
#
# Script:
# 02_GEO_microarray_limma_analysis.R
#
# Author:
# Zeinab Rohani
#
# Description:
# Differential gene expression analysis of GEO microarray
# data using the limma package.
#
###############################################################

############################
# Load Required Packages
############################

library(limma)

############################
# Create Output Directories
############################

dir.create("results", showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

############################
# Load Normalized Expression Matrix
############################

expression_matrix <- read.delim(
  "data/normalized_expression_matrix.txt",
  row.names = 1,
  check.names = FALSE
)

############################
# Load Annotation File
############################

annotation <- read.delim(
  "data/annotation.txt",
  check.names = FALSE
)

############################
# Define Experimental Groups
############################
# Modify this vector according to your dataset.

group <- factor(c(
  rep("Cancer", 12),
  rep("Normal", 6)
))

############################
# Construct Design Matrix
############################

design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

############################
# Fit Linear Model
############################

fit <- lmFit(expression_matrix, design)

############################
# Define Contrast
############################

contrast_matrix <- makeContrasts(
  Cancer - Normal,
  levels = design
)

############################
# Apply Contrast
############################

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

############################
# Empirical Bayes Moderation
############################

fit2 <- eBayes(fit2)

############################
# Extract Differentially
# Expressed Genes
############################

deg_results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

############################
# Add Gene Symbols
############################

if ("probe" %in% colnames(annotation)) {

    rownames(annotation) <- annotation$probe

} else {

    rownames(annotation) <- annotation[,1]

}

gene_symbol_column <- grep(
    "Gene",
    colnames(annotation),
    value = TRUE
)[1]

if (!is.na(gene_symbol_column)) {

    deg_results$GeneSymbol <-
        annotation[
            rownames(deg_results),
            gene_symbol_column
        ]

}

############################
# Order by Adjusted P-value
############################

deg_results <-
    deg_results[
        order(deg_results$adj.P.Val),
    ]

############################
# Significant DEGs
############################

significant_deg <-
    subset(
        deg_results,
        adj.P.Val < 0.05 &
        abs(logFC) > 1
    )

############################
# Save Results
############################

write.table(
    deg_results,
    file = "results/tables/all_DEGs.txt",
    sep = "\t",
    quote = FALSE
)

write.table(
    significant_deg,
    file = "results/tables/significant_DEGs.txt",
    sep = "\t",
    quote = FALSE
)

############################
# Summary
############################

cat("\n")
cat("---------------------------------\n")
cat("Microarray Differential Expression\n")
cat("---------------------------------\n")
cat("Total genes:", nrow(deg_results), "\n")
cat("Significant DEGs:", nrow(significant_deg), "\n")
cat("---------------------------------\n")

###############################################################
# End of Script
###############################################################
