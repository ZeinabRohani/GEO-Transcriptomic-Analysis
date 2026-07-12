###############################################################
# GEO-Transcriptomic-Analysis
#
# Script : 05_visualization.R
# Author : Zeinab Rohani
#
# Description:
# Generate publication-quality figures for
# differential expression analysis.
#
###############################################################

##############################
# Load Required Packages
##############################

required_packages <- c(
  "EnhancedVolcano",
  "pheatmap",
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
# Create Output Directory
##############################

dir.create(

  "results/figures",

  recursive = TRUE,

  showWarnings = FALSE

)

##############################
# Input Files
##############################

deg_file <-
  "results/tables/significant_DEGs.txt"

expression_file <-
  "data/processed/expression_matrix_final.txt"

##############################
# Check Files
##############################

if(!file.exists(deg_file))
  stop("Differential expression file not found.")

if(!file.exists(expression_file))
  stop("Expression matrix not found.")

##############################
# Read Files
##############################

deg <- read.delim(

  deg_file,

  check.names = FALSE,

  stringsAsFactors = FALSE

)

expression_matrix <- read.delim(

  expression_file,

  row.names = 1,

  check.names = FALSE,

  stringsAsFactors = FALSE

)

expression_matrix <- as.matrix(expression_matrix)

mode(expression_matrix) <- "numeric"

##############################
# Volcano Plot
##############################

pdf(

  "results/figures/Volcano_plot.pdf",

  width = 8,

  height = 7

)

EnhancedVolcano(

  deg,

  lab = deg$ProbeID,

  x = "logFC",

  y = "adj.P.Val",

  pCutoff = 0.05,

  FCcutoff = 1,

  pointSize = 2.5,

  labSize = 3,

  title = "Differential Gene Expression",

  subtitle = NULL,

  legendPosition = "right",

  gridlines.major = FALSE,

  gridlines.minor = FALSE

)

dev.off()

##############################
# MA Plot
##############################

pdf(

  "results/figures/MA_plot.pdf",

  width = 7,

  height = 6

)

plotMA(

  deg,

  main = "MA Plot"

)

dev.off()

##############################
# Select Top DEGs
##############################

deg <- deg[
  order(deg$adj.P.Val),
]

top_n <- min(50, nrow(deg))

top_genes <- deg$ProbeID[1:top_n]

##############################
# Heatmap Matrix
##############################

heatmap_matrix <-

  expression_matrix[
    rownames(expression_matrix) %in% top_genes,
  ]

##############################
# Scale by Gene
##############################

heatmap_matrix <- t(

  scale(

    t(heatmap_matrix)

  )

)

##############################
# Remove Missing Values
##############################

heatmap_matrix <-

  heatmap_matrix[
    complete.cases(heatmap_matrix),
  ]

##############################
# Heatmap
##############################

pdf(

  "results/figures/Heatmap_Top50_DEGs.pdf",

  width = 8,

  height = 10

)

pheatmap(

  heatmap_matrix,

  scale = "none",

  cluster_rows = TRUE,

  cluster_cols = TRUE,

  show_rownames = TRUE,

  show_colnames = TRUE,

  border_color = NA,

  fontsize = 8,

  main = "Top Differentially Expressed Genes"

)

dev.off()

##############################
# Sample Correlation
##############################

sample_cor <- cor(

  expression_matrix,

  method = "pearson"

)

##############################
# Correlation Heatmap
##############################

pdf(

  "results/figures/Sample_Correlation.pdf",

  width = 8,

  height = 8

)

pheatmap(

  sample_cor,

  cluster_rows = TRUE,

  cluster_cols = TRUE,

  display_numbers = FALSE,

  border_color = NA,

  main = "Sample Correlation"

)

dev.off()

##############################
# Principal Component Analysis
##############################

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

  cex = 1.3,

  xlab = paste0(

    "PC1 (",

    round(

      100 *
      summary(pca_result)$importance[2,1],

      2

    ),

    "%)"

  ),

  ylab = paste0(

    "PC2 (",

    round(

      100 *
      summary(pca_result)$importance[2,2],

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

  cex = 0.8

)

dev.off()

##############################
# Expression Density Plot
##############################

pdf(

  "results/figures/Expression_Density.pdf",

  width = 8,

  height = 6

)

plotDensities(

  expression_matrix,

  main = "Expression Density"

)

dev.off()

##############################
# Expression Boxplot
##############################

pdf(

  "results/figures/Expression_Boxplot.pdf",

  width = 10,

  height = 6

)

boxplot(

  expression_matrix,

  outline = FALSE,

  las = 2,

  main = "Normalized Expression",

  ylab = "Expression"

)

dev.off()

##############################
# Expression Histogram
##############################

pdf(

  "results/figures/Expression_Histogram.pdf",

  width = 8,

  height = 6

)

hist(

  as.numeric(expression_matrix),

  breaks = 100,

  main = "Expression Distribution",

  xlab = "Expression"

)

dev.off()

##############################
# Figure Summary
##############################

figure_summary <- data.frame(

  Figure = c(

    "Volcano Plot",
    "MA Plot",
    "Heatmap (Top DEGs)",
    "Sample Correlation",
    "PCA Plot",
    "Expression Density",
    "Expression Boxplot",
    "Expression Histogram"

  ),

  File = c(

    "Volcano_plot.pdf",
    "MA_plot.pdf",
    "Heatmap_Top50_DEGs.pdf",
    "Sample_Correlation.pdf",
    "PCA_plot.pdf",
    "Expression_Density.pdf",
    "Expression_Boxplot.pdf",
    "Expression_Histogram.pdf"

  ),

  stringsAsFactors = FALSE

)

write.table(

  figure_summary,

  file = "results/tables/figure_summary.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Save Session Information
##############################

writeLines(

  capture.output(sessionInfo()),

  con = "results/sessionInfo_visualization.txt"

)

##############################
# Final Console Summary
##############################

cat("\n")

cat("=========================================\n")

cat(" VISUALIZATION COMPLETED\n")

cat("=========================================\n")

cat("Generated Figures\n")

cat("-----------------------------------------\n")

cat("✓ Volcano_plot.pdf\n")
cat("✓ MA_plot.pdf\n")
cat("✓ Heatmap_Top50_DEGs.pdf\n")
cat("✓ Sample_Correlation.pdf\n")
cat("✓ PCA_plot.pdf\n")
cat("✓ Expression_Density.pdf\n")
cat("✓ Expression_Boxplot.pdf\n")
cat("✓ Expression_Histogram.pdf\n")
cat("-----------------------------------------\n")
cat("Additional Files\n")
cat("-----------------------------------------\n")
cat("✓ figure_summary.txt\n")
cat("✓ sessionInfo_visualization.txt\n")
cat("=========================================\n")

message("Visualization completed successfully.")

###############################################################
# End of Script
###############################################################
