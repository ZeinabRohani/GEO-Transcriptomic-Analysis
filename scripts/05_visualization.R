############################################################
# Project: GEO-Transcriptomic-Analysis
# Script : 05_visualization.R
#
# Description:
# Generate publication-quality figures
# for differential expression and
# functional enrichment analysis.
############################################################

## Load packages
library(ggplot2)
library(pheatmap)
library(clusterProfiler)

############################################################
# Load expression matrix
############################################################

expr <- read.delim(
  "data/expression_normalized.txt",
  row.names = 1,
  check.names = FALSE
)

############################################################
# Load differential expression results
############################################################

deg <- read.delim(
  "results/limma_all_genes.txt",
  row.names = 1,
  check.names = FALSE
)

############################################################
# Volcano Plot
############################################################

deg$Significant <- "No"

deg$Significant[
  deg$adj.P.Val < 0.05 &
  abs(deg$logFC) > 1
] <- "Yes"

volcano <- ggplot(
  deg,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Significant
  )
) +
  geom_point(size = 1.5) +
  theme_bw() +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

ggsave(
  "results/VolcanoPlot.png",
  volcano,
  width = 7,
  height = 6,
  dpi = 300
)

############################################################
# Heatmap of Top 50 Differentially Expressed Genes
############################################################

top50 <- rownames(deg)[1:50]

heatmap_matrix <- expr[top50, ]

pheatmap(
  heatmap_matrix,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  filename = "results/Heatmap.png",
  width = 8,
  height = 10
)

############################################################
# Load Enrichment Results
############################################################

go <- read.delim(
  "results/GO_BP_results.txt"
)

kegg <- read.delim(
  "results/KEGG_results.txt"
)

############################################################
# GO Dotplot
############################################################

go_result <- readRDS("results/go_result.rds")

png(
  "results/GO_Dotplot.png",
  width = 2000,
  height = 1600,
  res = 300
)

dotplot(
  go_result,
  showCategory = 15
)

dev.off()

############################################################
# KEGG Dotplot
############################################################

kegg_result <- readRDS("results/kegg_result.rds")

png(
  "results/KEGG_Dotplot.png",
  width = 2000,
  height = 1600,
  res = 300
)

dotplot(
  kegg_result,
  showCategory = 15
)

dev.off()

############################################################
# Summary
############################################################

cat("\n")
cat("=====================================\n")
cat("VISUALIZATION COMPLETED\n")
cat("=====================================\n")
cat("Generated figures:\n")
cat("- VolcanoPlot.png\n")
cat("- Heatmap.png\n")
cat("- GO_Dotplot.png\n")
cat("- KEGG_Dotplot.png\n")
cat("=====================================\n")
