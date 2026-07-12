###############################################################
# GEO-Transcriptomic-Analysis
#
# Script : 04_RNAseq_edgeR_analysis.R
# Author : Zeinab Rohani
#
# Description:
# Differential expression analysis of RNA-seq count data
# using edgeR + limma-voom
#
###############################################################

##############################
# Load Required Packages
##############################

required_packages <- c(
  "edgeR",
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
# User Settings
##############################

count_file <- "data/raw/count_matrix.txt"

###############################################################
# IMPORTANT
# Modify the group labels below according to your dataset.
###############################################################

group <- factor(

  c(

    rep("Control",3),

    rep("Treatment",3)

  )

)

##############################
# Check Input File
##############################

if(!file.exists(count_file)){

  stop(
    paste(
      "File not found:",
      count_file
    )
  )

}

##############################
# Read Count Matrix
##############################

message("--------------------------------")
message("Loading RNA-seq count matrix...")
message("--------------------------------")

counts <- read.delim(

  count_file,

  row.names = 1,

  check.names = FALSE,

  stringsAsFactors = FALSE

)

counts <- as.matrix(counts)

mode(counts) <- "numeric"

##############################
# Dataset Summary
##############################

cat(
  "Genes   :",
  nrow(counts),
  "\n"
)

cat(
  "Samples :",
  ncol(counts),
  "\n"
)

##############################
# Check Sample Number
##############################

if(length(group) != ncol(counts)){

  stop(
    "Number of samples and group labels do not match."
  )

}

##############################
# Create DGEList Object
##############################

dge <- DGEList(

  counts = counts,

  group = group

)

##############################
# Library Sizes
##############################

message("--------------------------------")
message("Library Sizes")
message("--------------------------------")

print(dge$samples)

##############################
# Save Library Information
##############################

write.table(

  dge$samples,

  file = "results/tables/library_information.txt",

  sep = "\t",

  quote = FALSE,

  row.names = TRUE

)

##############################
# Filter Lowly Expressed Genes
##############################

keep <- filterByExpr(

  dge,

  group = group

)

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

cat(
  "Genes retained :",
  nrow(dge),
  "\n"
)

##############################
# TMM Normalization
##############################

message("--------------------------------")
message("Performing TMM normalization...")
message("--------------------------------")

dge <- calcNormFactors(

  dge,

  method = "TMM"

)

##############################
# Save Normalization Factors
##############################

write.table(

  dge$samples,

  file = "results/tables/TMM_normalization_factors.txt",

  sep = "\t",

  quote = FALSE,

  row.names = TRUE

)

##############################
# Design Matrix
##############################

design <- model.matrix(

  ~0 + group

)

colnames(design) <- levels(group)

##############################
# Save Design Matrix
##############################

write.table(

  design,

  file = "results/tables/design_matrix_RNAseq.txt",

  sep = "\t",

  quote = FALSE,

  row.names = TRUE

)

##############################
# Voom Transformation
##############################

message("--------------------------------")
message("Running voom transformation...")
message("--------------------------------")

voom_object <- voom(

  dge,

  design,

  plot = TRUE

)

##############################
# Save Voom Plot
##############################

dev.copy(

  pdf,

  file = "results/figures/voom_mean_variance.pdf",

  width = 7,

  height = 6

)

dev.off()

##############################
# Extract Normalized Expression
##############################

normalized_expression <- voom_object$E

##############################
# Save Normalized Matrix
##############################

write.table(

  normalized_expression,

  file = "data/processed/RNAseq_voom_expression.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Boxplot
##############################

pdf(

  "results/figures/RNAseq_boxplot.pdf",

  width = 10,

  height = 6

)

boxplot(

  normalized_expression,

  outline = FALSE,

  las = 2,

  main = "Normalized RNA-seq Expression",

  ylab = "logCPM"

)

dev.off()

##############################
# Density Plot
##############################

pdf(

  "results/figures/RNAseq_density_plot.pdf",

  width = 8,

  height = 6

)

plotDensities(

  normalized_expression,

  main = "RNA-seq Density Plot"

)

dev.off()

##############################
# User-defined Contrast
##############################

###############################################################
# Modify the comparison below according to your dataset.
#
# Examples:
# Treatment-Control
# Tumor-Normal
# Disease-Healthy
###############################################################

contrast_formula <- "Treatment-Control"

##############################
# Create Contrast Matrix
##############################

contrast_matrix <- makeContrasts(

  contrasts = contrast_formula,

  levels = design

)

write.table(

  contrast_matrix,

  file = "results/tables/contrast_matrix_RNAseq.txt",

  sep = "\t",

  quote = FALSE,

  row.names = TRUE

)

##############################
# Linear Model
##############################

message("--------------------------------")
message("Fitting linear model...")
message("--------------------------------")

fit <- lmFit(

  voom_object,

  design

)

##############################
# Apply Contrast
##############################

fit <- contrasts.fit(

  fit,

  contrast_matrix

)

##############################
# Empirical Bayes
##############################

fit <- eBayes(

  fit

)

##############################
# Differential Expression Table
##############################

deg <- topTable(

  fit,

  number = Inf,

  adjust.method = "BH",

  sort.by = "P"

)

deg$GeneID <- rownames(deg)

##############################
# User-defined Thresholds
##############################

logFC_cutoff <- 1

adjP_cutoff <- 0.05

##############################
# Significant Genes
##############################

deg_sig <- subset(

  deg,

  abs(logFC) >= logFC_cutoff &
    adj.P.Val < adjP_cutoff

)

##############################
# Upregulated Genes
##############################

deg_up <- subset(

  deg_sig,

  logFC >= logFC_cutoff

)

##############################
# Downregulated Genes
##############################

deg_down <- subset(

  deg_sig,

  logFC <= -logFC_cutoff

)

##############################
# DEG Summary
##############################

deg_summary <- data.frame(

  Total_Genes = nrow(deg),

  Significant = nrow(deg_sig),

  Upregulated = nrow(deg_up),

  Downregulated = nrow(deg_down),

  stringsAsFactors = FALSE

)

write.table(

  deg_summary,

  file = "results/tables/RNAseq_DEG_summary.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export All DEGs
##############################

write.table(

  deg,

  file = "results/tables/RNAseq_all_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Significant DEGs
##############################

write.table(

  deg_sig,

  file = "results/tables/RNAseq_significant_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Upregulated Genes
##############################

write.table(

  deg_up,

  file = "results/tables/RNAseq_upregulated_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Downregulated Genes
##############################

write.table(

  deg_down,

  file = "results/tables/RNAseq_downregulated_DEGs.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Volcano Plot Input
##############################

write.table(

  deg,

  file = "results/tables/RNAseq_volcano_input.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

##############################
# Export Heatmap Input
##############################

heatmap_expression <- normalized_expression[
  rownames(normalized_expression) %in% deg_sig$GeneID,
]

write.table(

  heatmap_expression,

  file = "results/tables/RNAseq_heatmap_expression.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Save Normalized Expression
##############################

write.table(

  normalized_expression,

  file = "data/processed/RNAseq_normalized_expression.txt",

  sep = "\t",

  quote = FALSE,

  col.names = NA

)

##############################
# Save Fitted Model
##############################

save(

  fit,

  file = "results/tables/RNAseq_limma_fit.RData"

)

##############################
# Save Session Information
##############################

writeLines(

  capture.output(sessionInfo()),

  con = "results/sessionInfo_RNAseq.txt"

)

##############################
# Final Console Summary
##############################

cat("\n")

cat("=========================================\n")

cat(" RNA-seq DIFFERENTIAL EXPRESSION ANALYSIS\n")

cat("=========================================\n")

cat("Total genes          :", nrow(deg), "\n")

cat("Significant genes    :", nrow(deg_sig), "\n")

cat("Upregulated genes    :", nrow(deg_up), "\n")

cat("Downregulated genes  :", nrow(deg_down), "\n")

cat("=========================================\n")

cat("Output Files\n")

cat("-----------------------------------------\n")

cat("✓ RNAseq_all_DEGs.txt\n")
cat("✓ RNAseq_significant_DEGs.txt\n")
cat("✓ RNAseq_upregulated_DEGs.txt\n")
cat("✓ RNAseq_downregulated_DEGs.txt\n")
cat("✓ RNAseq_volcano_input.txt\n")
cat("✓ RNAseq_heatmap_expression.txt\n")
cat("✓ RNAseq_DEG_summary.txt\n")
cat("✓ RNAseq_normalized_expression.txt\n")
cat("✓ RNAseq_limma_fit.RData\n")
cat("✓ design_matrix_RNAseq.txt\n")
cat("✓ contrast_matrix_RNAseq.txt\n")
cat("✓ library_information.txt\n")
cat("✓ TMM_normalization_factors.txt\n")
cat("✓ voom_mean_variance.pdf\n")
cat("✓ RNAseq_boxplot.pdf\n")
cat("✓ RNAseq_density_plot.pdf\n")
cat("✓ sessionInfo_RNAseq.txt\n")

cat("=========================================\n")

message("RNA-seq differential expression analysis completed successfully.")

###############################################################
# End of Script
###############################################################
