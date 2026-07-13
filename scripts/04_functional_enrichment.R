############################################################
# Project: GEO-Transcriptomic-Analysis
# Script : 04_functional_enrichment.R
#
# Description:
# Functional enrichment analysis of
# differentially expressed genes using
# Gene Ontology (GO) and KEGG pathways.
############################################################

## Load packages
library(clusterProfiler)
library(org.Hs.eg.db)

############################################################
# Load DEG results
############################################################

deg <- read.delim(
  "results/limma_significant_genes.txt",
  row.names = 1,
  check.names = FALSE
)

############################################################
# Convert Gene Symbols to Entrez IDs
############################################################

gene.df <- bitr(

  rownames(deg),

  fromType = "SYMBOL",

  toType = "ENTREZID",

  OrgDb = org.Hs.eg.db

)

############################################################
# GO Biological Process
############################################################

go_bp <- enrichGO(

  gene          = gene.df$ENTREZID,

  OrgDb         = org.Hs.eg.db,

  keyType       = "ENTREZID",

  ont           = "BP",

  pAdjustMethod = "BH",

  pvalueCutoff  = 0.05,

  qvalueCutoff  = 0.05

)

############################################################
# KEGG Pathway Analysis
############################################################

kegg <- enrichKEGG(

  gene = gene.df$ENTREZID,

  organism = "hsa",

  pvalueCutoff = 0.05

)

############################################################
# Save Results
############################################################

write.table(

  as.data.frame(go_bp),

  file = "results/GO_BP_results.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

write.table(

  as.data.frame(kegg),

  file = "results/KEGG_results.txt",

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

############################################################
# Summary
############################################################

cat("\n")
cat("=====================================\n")
cat("FUNCTIONAL ENRICHMENT COMPLETED\n")
cat("=====================================\n")
cat("GO terms     :", nrow(as.data.frame(go_bp)), "\n")
cat("KEGG pathways:", nrow(as.data.frame(kegg)), "\n")
cat("=====================================\n")
