####Fernando Zamuner

# Load required packages
if (!requireNamespace("differential.coverage", quietly = TRUE)) {
  devtools::install_github("favorov/differential.coverage")
}
library(differential.coverage)
library(rtracklayer)
library(dplyr)
library(plyranges)
library(tidyr)
library(EnhancedVolcano)

setwd("~/data/")

# Get the genes
the_genes <- differential.coverage::get.Known.Gene.List(genome.annotation.id = "gencode19")
nonunique_names <- table(the_genes$ensembl) > 1
the_genes <- the_genes %>% filter((!ensembl %in% nonunique_names) | seqnames != "chrY")

# Function to attach gene names
attach_gene_names <- function(rds) {
  merge(as.data.frame(rds), as.data.frame(the_genes), by.x = "row.names", by.y = "ensembl")
}

# Read RDS files
mrna_O47 <- readRDS("JQ1_vs_DMSO_DE_results_O47_JQ1.rds")
mrna_O90 <- readRDS("JQ1_vs_DMSO_DE_results_O90_JQ1.rds")
mrna_pts <- readRDS("T_to_N_DE_results_hg19_ensembl_77_samples.rds")

# Build a named list for gene names
ensembl_to_genes <- setNames(the_genes$gene_name, the_genes$ensembl)


###Fig3_ABC
# Volcano plot function
volcano_plot <- function(rds) {
  this_rds <- deparse(substitute(rds))
  folder_path <- "volcano"
  
  # Check if the folder exists; if not, create it
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
  }
  
  pdfpath <- file.path(folder_path, paste0("volcano-", this_rds, ".pdf"))
  
  plot <- EnhancedVolcano(
    rds,
    maxoverlapsConnectors = Inf,
    lab = ensembl_to_genes[rownames(rds)],
    x = "log2FoldChange",
    y = "pvalue",
    title = this_rds,
    subtitle = NULL,
    pointSize = 1,
    pCutoff = 0.05,
    pCutoffCol = "padj",
    FCcutoff = 1.5,
    col = c("black", "black", "black", "red3"),
    legendPosition = "none",
    xlab = bquote(~Log[2]~ "fold change")
  )
  
  print(plot)
  ggsave(file = pdfpath, plot = plot)
}

# Call volcano_plot for each dataset
volcano_plot(mrna_O47)
volcano_plot(mrna_O90)
volcano_plot(mrna_pts)


###fig3_DE
########UP in Patients
# Volcano plot function
volcano_plot <- function(rds) {
  this_rds <- deparse(substitute(rds))
  folder_path <- "volcano"
  
  # Check if the folder exists; if not, create it
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
  }
  
  pdfpath <- file.path(folder_path, paste0("volcano_mod-", this_rds, ".pdf"))
  
  keyvals <- ifelse(
    rds$log2FoldChange < 0, 'lightgrey',
    ifelse(rds$log2FoldChange > 0 & rds$padj < 0.05, 'red', 'lightgrey')
  )
  
  keyvals[is.na(keyvals)] <- 'lightgrey'
  names(keyvals)[keyvals == 'red'] <- 'high'
  names(keyvals)[keyvals == 'lightgrey' & rds$log2FoldChange < 0] <- 'low'
  names(keyvals)[keyvals == 'lightgrey' & rds$log2FoldChange >= 0] <- 'mid'
  
  
  plot <- EnhancedVolcano(
    rds,
    maxoverlapsConnectors = Inf,
    lab = ensembl_to_genes[rownames(rds)],
    x = "log2FoldChange",
    y = "pvalue",
    title = this_rds,
    subtitle = NULL,
    pointSize = 1,
    colCustom = keyvals,
    labSize = 0,
    pCutoff = 0.05,
    pCutoffCol = "padj",
    FCcutoff = 0,
    legendPosition = "none",
    xlab = bquote(~Log[2]~ "Fold Change")
  )
  
  print(plot)
  ggsave(file = pdfpath, plot = plot)
}

# Call volcano_plot for each dataset
volcano_plot(mrna_pts)

########DOWN in Cell Lines
# Volcano plot function
volcano_plot <- function(rds) {
  this_rds <- deparse(substitute(rds))
  folder_path <- "volcano"
  
  # Check if the folder exists; if not, create it
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
  }
  
  pdfpath <- file.path(folder_path, paste0("volcano_mod-", this_rds, ".pdf"))
  
  keyvals <- ifelse(
    rds$log2FoldChange > 0, 'lightgrey',
    ifelse(rds$log2FoldChange < 0 & rds$padj < 0.05, 'green', 'lightgrey')
  )
  
  keyvals[is.na(keyvals)] <- 'lightgrey'
  names(keyvals)[keyvals == 'green'] <- 'high'
  names(keyvals)[keyvals == 'lightgrey' & rds$log2FoldChange < 0] <- 'low'
  names(keyvals)[keyvals == 'lightgrey' & rds$log2FoldChange >= 0] <- 'mid'
  
  plot <- EnhancedVolcano(
    rds,
    maxoverlapsConnectors = Inf,
    lab = ensembl_to_genes[rownames(rds)],
    x = "log2FoldChange",
    y = "pvalue",
    title = this_rds,
    subtitle = NULL,
    pointSize = 1,
    colCustom = keyvals,
    labSize = 0,
    pCutoff = 0.05,
    pCutoffCol = "padj",
    FCcutoff = 0,
    legendPosition = "none",
    xlab = bquote(~Log[2]~ "Fold Change")
  )
  
  print(plot)
  ggsave(file = pdfpath, plot = plot)
}

# Call volcano_plot for each dataset
volcano_plot(mrna_O47)
volcano_plot(mrna_O90)

