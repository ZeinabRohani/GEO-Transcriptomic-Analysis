###############################################################
# GEO-Transcriptomic-Analysis
#
# Script : 03_limma_analysis.R
# Author : Zeinab Rohani
#
# Description:
# Differential gene expression analysis
# using limma for GEO microarray datasets
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
# Input Files
##############################

expression_file <-
  "data/processed/expression_matrix_final.txt"

phenotype_file <-
  "data/raw/phenotype_data.txt"

annotation_file <-
  "data/raw/feature_annotation.txt"

##############################
# Check Input Files
##############################

input_files <- c(

  expression_file,

  phenotype_file,

  annotation_file

)

for(f in input_files){

  if(!file.exists(f)){

    stop(
      paste(
        "Missing input file:",
        f
      )
    )

  }

}

##############################
# Read Expression Matrix
##############################

message("--------------------------------")
message("Loading processed expression matrix...")
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
# Read Phenotype Data
##############################

phenotype <- read.delim(

  phenotype_file,

  row.names = 1,

  check.names = FALSE,

  stringsAsFactors = FALSE

)

##############################
# Read Annotation
##############################

annotation <- read.delim(

  annotation_file,

  row.names = 1,

  check.names = FALSE,

  stringsAsFactors = FALSE

)

##############################
# Check Sample Names
##############################

if(!all(colnames(expression_matrix) %in%
        rownames(phenotype))){

  stop(
    "Sample names in the expression matrix do not match phenotype data."
  )

}

##############################
# Reorder Phenotype
##############################

phenotype <-

  phenotype[
    colnames(expression_matrix),
  ]

##############################
# USER SETTINGS
##############################

###############################################################
# Replace the vector below with the
# experimental groups of your dataset.
###############################################################

group <- factor(

  c(

    rep("Cancer",12),

    rep("Normal",6)

  )

)

##############################
# Check Group Length
##############################

if(length(group) != ncol(expression_matrix)){

  stop(
    "Number of samples and group labels do not match."
  )

}

##############################
# Display Group Information
##############################

message("--------------------------------")
message("Experimental Groups")
message("--------------------------------")

print(table(group))

##############################
# Design Matrix
##############################

design <- model.matrix(

  ~0 + group

)

colnames(design) <- levels(group)

message("--------------------------------")
message("Design Matrix")
message("--------------------------------")

print(design)

##############################
# Save Design Matrix
##############################

write.table(

  design,

  file = "results/tables/design_matrix.txt",

  sep = "\t",

  quote = FALSE,

  row.names = TRUE

)

##############################
# User-defined Contrast
##############################

###############################################################
# Modify the comparison below according to your dataset.
#
# Examples:
# Cancer - Normal
# Tumor - Control
# Treated - Untreated
###############################################################

contrast_formula <- "Cancer-Normal"

##############################
# Create Contrast Matrix
##############################

contrast_matrix <- makeContrasts(

  contrasts = contrast_formula,

  levels = design

)

message("--------------------------------")
message("Contrast Matrix")
message("--------------------------------")

print(contrast_matrix)

##############################
# Save Contrast Matrix
##############################

write.table(

  contrast_matrix,

  file = "results/tables/contrast_matrix.txt",

  sep = "\t",

  quote = FALSE,

  row.names = TRUE

)

##############################
# Fit Linear Model
##############################

message("--------------------------------")
message("Fitting linear model...")
message("--------------------------------")

fit <- lmFit(

  expression_matrix,

  design

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
# Model Summary
##############################

message("--------------------------------")
message("Model Summary")
message("--------------------------------")

print(summary(decideTests(fit2)))

##############################
# Save Model Object
##############################

save(

  fit2,

  file = "results/tables/limma_fit.RData"

)

###############################################################
# Part 3
# Differential Expression Analysis
###############################################################

##############################
# Extract Differentially
# Expressed Genes
##############################

deg <- topTable(

  fit2,

  number = Inf,

  adjust.method = "BH",

  sort.by = "P"

)

##############################
# Add Gene IDs
##############################

deg$ProbeID <- rownames(deg)

##############################
# Add Gene Annotation
##############################

annotation$ProbeID <- rownames(annotation)

deg <- merge(

  deg,

  annotation,

  by = "ProbeID",

  all.x = TRUE,

  sort = FALSE

)

##############################
# Reorder Columns
##############################

first_columns <- c(

  "ProbeID",

  "logFC",

  "AveExpr",

  "t",

  "P.Value",

  "adj.P.Val",

  "B"

)

remaining_columns <-

  setdiff(

    colnames(deg),

    first_columns

  )

deg <- deg[ ,

  c(

    first_columns,

    remaining_columns

  )

]

##############################
# Order by Adjusted P-value
##############################

deg <-

  deg[
    order(
      deg$adj.P.Val
    ),
  ]

##############################
# User-defined Thresholds
##############################

logFC_cutoff <- 1

adjP_cutoff <- 0.05

##############################
# Significant DEGs
##############################

deg_sig <-

  subset(

    deg,

    abs(logFC) >= logFC_cutoff &
      adj.P.Val < adjP_cutoff

  )

##############################
# Upregulated Genes
##############################

deg_up <-

  subset(

    deg_sig,

    logFC >= logFC_cutoff

  )

##############################
# Downregulated Genes
##############################

deg_down <-

  subset(

    deg_sig,

    logFC <= -logFC_cutoff

  )

##############################
# Summary
##############################

summary_table <- data.frame(

  Total_Genes = nrow(deg),

  Significant = nrow(deg_sig),

  Upregulated = nrow(deg_up),

  Downregulated = nrow(deg_down),

  logFC_Cutoff = logFC_cutoff,

  Adjusted_Pvalue = adjP_cutoff,

  stringsAsFactors = FALSE

)

write.table(

  summary_table,

  file = "results/tables/DEG_summary.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export All Differentially
# Expressed Genes
##############################

write.table(

  deg,

  file = "results/tables/all_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Significant DEGs
##############################

write.table(

  deg_sig,

  file = "results/tables/significant_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Upregulated Genes
##############################

write.table(

  deg_up,

  file = "results/tables/upregulated_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Downregulated Genes
##############################

write.table(

  deg_down,

  file = "results/tables/downregulated_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Volcano Plot Input
##############################

write.table(

  deg,

  file = "results/tables/volcano_input.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Heatmap Input
##############################

heatmap_expression <- expression_matrix[
  rownames(expression_matrix) %in% deg_sig$ProbeID,
]

write.table(

  heatmap_expression,

  file = "results/tables/heatmap_expression.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Save Session Information
##############################

writeLines(

  capture.output(sessionInfo()),

  con = "results/sessionInfo_limma.txt"

)

##############################
# Console Summary
##############################

cat("\n")

cat("=========================================\n")

cat(" LIMMA DIFFERENTIAL EXPRESSION ANALYSIS\n")

cat("=========================================\n")

cat("Total genes          :", nrow(deg), "\n")

cat("Significant genes    :", nrow(deg_sig), "\n")

cat("Upregulated genes    :", nrow(deg_up), "\n")

cat("Downregulated genes  :", nrow(deg_down), "\n")

cat("=========================================\n")

cat("Output Files\n")

cat("-----------------------------------------\n")

cat("✓ all_DEGs.txt\n")

cat("✓ significant_DEGs.txt\n")

cat("✓ upregulated_DEGs.txt\n")

cat("✓ downregulated_DEGs.txt\n")

cat("✓ volcano_input.txt\n")

cat("✓ heatmap_expression.txt\n")

cat("✓ DEG_summary.txt\n")

cat("✓ design_matrix.txt\n")

cat("✓ contrast_matrix.txt\n")

cat("✓ limma_fit.RData\n")

cat("✓ sessionInfo_limma.txt\n")

cat("=========================================\n")

message("Differential expression analysis completed successfully.")

###############################################################
# End of Script
###############################################################
