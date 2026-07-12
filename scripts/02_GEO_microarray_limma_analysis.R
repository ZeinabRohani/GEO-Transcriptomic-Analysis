###############################################################
# Project : GEO-Transcriptomic-Analysis
#
# Script  : 02_GEO_microarray_limma_analysis.R
#
# Author  : Zeinab Rohani
#
# Description:
# Differential gene expression analysis using limma
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
  "results",
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

##############################
# Load Expression Matrix
##############################

message("-----------------------------------")
message("Loading normalized expression matrix")
message("-----------------------------------")

expression_matrix <- read.delim(

  "data/expression_matrix_normalized.txt",

  row.names = 1,

  check.names = FALSE

)

##############################
# Load Annotation
##############################

annotation <- read.delim(

  "data/annotation.txt",

  check.names = FALSE,

  stringsAsFactors = FALSE

)

##############################
# Define Experimental Groups
##############################
#
# IMPORTANT
#
# Modify ONLY this section if your
# dataset contains different groups.
#
##############################

group <- factor(

  c(

    rep("Cancer",12),

    rep("Normal",6)

  )

)

##############################
# Check Sample Number
##############################

if(length(group) != ncol(expression_matrix)){

  stop(

    "Number of samples does not match group labels."

  )

}

##############################
# Design Matrix
##############################

design <- model.matrix(

  ~0 + group

)

colnames(design) <- levels(group)

message("-----------------------------------")
message("Design matrix")
message("-----------------------------------")

print(design)

##############################
# Linear Model
##############################

fit <- lmFit(

  expression_matrix,

  design

)

##############################
# Contrast Matrix
##############################

contrast_matrix <- makeContrasts(

  Cancer - Normal,

  levels = design

)

##############################
# Apply Contrast
##############################

fit2 <- contrasts.fit(

  fit,

  contrast_matrix

)

##############################
# Empirical Bayes
##############################

fit2 <- eBayes(

  fit2

)

##############################
# Extract All Genes
##############################

deg <- topTable(

  fit2,

  number = Inf,

  adjust.method = "BH",

  sort.by = "P"

)

message("Differential expression completed")

##############################
# Find Gene Symbol Column
##############################

possible_columns <- c(

  "Gene.symbol",
  "Gene Symbol",
  "GENE_SYMBOL",
  "Symbol",
  "SYMBOL",
  "GeneSymbol",
  "gene_symbol"

)

gene_column <- intersect(
  possible_columns,
  colnames(annotation)
)

if(length(gene_column) == 0){

  message("-----------------------------------")
  message("Gene symbol column was not found.")
  message("Probe IDs will be retained.")
  message("-----------------------------------")

  deg$GeneSymbol <- rownames(deg)

}else{

  rownames(annotation) <- annotation[,1]

  deg$GeneSymbol <- annotation[
    rownames(deg),
    gene_column[1]
  ]

}

##############################
# Remove Missing Symbols
##############################

deg$GeneSymbol[
  deg$GeneSymbol == ""
] <- NA

##############################
# Reorder Columns
##############################

deg <- deg[ ,

  c(

    "GeneSymbol",

    "logFC",

    "AveExpr",

    "t",

    "P.Value",

    "adj.P.Val",

    "B"

  )

]

##############################
# Order by Adjusted P-value
##############################

deg <- deg[
  order(deg$adj.P.Val),
]

##############################
# Significant Genes
##############################

significant_deg <- subset(

  deg,

  adj.P.Val < 0.05 &
    abs(logFC) >= 1

)

##############################
# Upregulated Genes
##############################

upregulated_deg <- subset(

  significant_deg,

  logFC >= 1

)

##############################
# Downregulated Genes
##############################

downregulated_deg <- subset(

  significant_deg,

  logFC <= -1

)

##############################
# Save All Results
##############################

write.table(

  deg,

  file = "results/tables/all_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

write.table(

  significant_deg,

  file = "results/tables/significant_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

write.table(

  upregulated_deg,

  file = "results/tables/upregulated_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

write.table(

  downregulated_deg,

  file = "results/tables/downregulated_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Summary Statistics
##############################

summary_table <- data.frame(

  Total_Genes = nrow(deg),

  Significant_Genes = nrow(significant_deg),

  Upregulated = nrow(upregulated_deg),

  Downregulated = nrow(downregulated_deg),

  LogFC_Cutoff = 1,

  Adjusted_Pvalue_Cutoff = 0.05

)

write.table(

  summary_table,

  file = "results/tables/analysis_summary.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Save Volcano Input
##############################

write.table(

  deg,

  file = "results/tables/volcano_input.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Save Heatmap Input
##############################

heatmap_genes <- rownames(significant_deg)

expression_heatmap <- expression_matrix[
  heatmap_genes,
]

write.table(

  expression_heatmap,

  file = "results/tables/heatmap_expression.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Console Report
##############################

cat("\n")

cat(" Differential Expression Analysis Summary\n")

cat(
  "Total genes          :",
  nrow(deg),
  "\n"
)

cat(
  "Significant genes    :",
  nrow(significant_deg),
  "\n"
)

cat(
  "Upregulated genes    :",
  nrow(upregulated_deg),
  "\n"
)

cat(
  "Downregulated genes  :",
  nrow(downregulated_deg),
  "\n"
)

cat("=========================================\n")

cat("Output files\n")

cat("-----------------------------------------\n")

cat("results/tables/all_DEGs.txt\n")

cat("results/tables/significant_DEGs.txt\n")

cat("results/tables/upregulated_DEGs.txt\n")

cat("results/tables/downregulated_DEGs.txt\n")

cat("results/tables/analysis_summary.txt\n")

cat("results/tables/volcano_input.txt\n")

cat("results/tables/heatmap_expression.txt\n")

###############################################################
# End of Script
###############################################################
