# Load required libraries
library(dplyr)
library(org.Hs.eg.db)
library(stringr)
library(RCurl)
library(fgsea)
library(igraph)
library(paletteer)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(qgraph)
library(gplots)
library(ComplexHeatmap)
library(circlize)

# Set working directory
setwd("~/data/")

# Load HNSCC data
df.tn <- as.data.frame(readRDS("fold_and_distance_77_samples_SED.rds"))
df.tn <- fold_and_distance_77_samples_SED

# Load JQ1 data for 047 and 090
df.047 <- readRDS("../RNAseq-gencode19/JQ1_vs_DMSO_DE_results_O47_JQ1.rds")
df.090 <- readRDS("../RNAseq-gencode19/JQ1_vs_DMSO_DE_results_O90_JQ1.rds")

# Preprocessing
names(df.047) <- str_c(names(df.047), "_047")
names(df.090) <- str_c(names(df.090), "_090")

mtch_047 <- match(rownames(df.tn), rownames(df.047))
mtch_090 <- match(rownames(df.tn), rownames(df.090))

df <- as.data.frame(cbind(df.tn, df.047[mtch_047,], df.090[mtch_090,]))

# Make custom gene lists representing TF target genes
TF.df <- readRDS("TF.enrichment_analysis.rds")
TFs <- unique(subset(TF.df, padj <= 0.05 & delta > 0)$TF)

trustt_URL <- "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
myfile <- getURL(trustt_URL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
TF_paths.df <- read.table(textConnection(myfile), header=F)

TF_gene_list <- list()
for (TF in TFs) {
  TF_genes <- unique(as.character(subset(TF_paths.df, V1 == TF)$V2))
  TF_gene_list[[TF]] <- data.frame(TF_genes = TF_genes)
}

expanded_TF_gene_list <- TF_gene_list
TF_in_trustt <- unique(as.character(TF_paths.df$V1))
for (TF in names(TF_gene_list)) {
  TF_df <- TF_gene_list[[TF]]
  genes_are_TFs <- TF_df$TF_genes %in% TF_in_trustt & !TF_df$TF_genes %in% names(TF_gene_list)
  if (sum(genes_are_TFs) > 0) {
    for (TF2 in TF_df$TF_genes[genes_are_TFs]) {
      TF2_genes <- as.character(subset(TF_paths.df, V1 == TF2)$V2)
      expanded_TF_gene_list[[TF2]] <- data.frame(TF_genes = TF2_genes)
    }
  }
}

only_TFs_list <- lapply(TF_gene_list, function(x) {
  hit <- x$TF_genes %in% names(TF_gene_list)
  return(data.frame(TF_genes = as.character(x$TF_genes[hit])))
})

# Prepare nodes and edges
edges <- dplyr::bind_rows(only_TFs_list, .id = "column_label")
names(edges) <- c("label", "pair")
edges <- edges[!edges$label == edges$pair,]

nodes <- as.data.frame(unique(c(names(only_TFs_list), as.character(edges$label), as.character(edges$pair))))
names(nodes) <- "label"

nodes$status <- "Gene"
nodes$status[nodes$label %in% names(TF_gene_list)] <- "TF"

mtch <- match(nodes$label, df$gene_name)
nodes$log2FC_patients <- df$log2FoldChange[mtch]
nodes$log2FC_JQ1047 <- df$log2FoldChange_047[mtch]
nodes$log2FC_JQ1090 <- df$log2FoldChange_090[mtch]
nodes$padj_patients <- df$padj[mtch]
nodes$padj_JQ1_047 <- df$padj_047[mtch]
nodes$padj_JQ1_090 <- df$padj_090[mtch]

# Rename columns for clarity
names(nodes)[3:8] <- c("log2FC_patients", "log2FC_JQ1047", "log2FC_JQ1090", "padj_patients", "padj_JQ1047", "padj_JQ1090")

# Rename row names
rownames(nodes) <- nodes$label

# Create heatmap matrix for log2FC, log2FC_JQ1047, and log2FC_JQ1090
heatmap_matrix <- as.matrix(nodes[, c("log2FC_patients", "log2FC_JQ1047", "log2FC_JQ1090")])

# Create annotations for significant padj values
significant_padj <- as.matrix(nodes[, c("padj_patients", "padj_JQ1047", "padj_JQ1090")])

# Ensure the rows in heatmap_matrix and significant_padj are ordered according to TF_genes_enr_order
TF_genes_enr_order <- c("TP63", "JUNB", "GRHL2", "JUND", "EPAS1", "KLF3", "FOS", "FOSL1", "SNAI2", "JUN", "ARNT", "TP53", "FOSL2", "KLF5", "TP73", "PPARG", "E2F3", "NKX2-1", "TFAP2C", "SMAD4", "SMAD3", "KLF4", "NRF1", "TEAD4", "ESR1", "HIF1A", "TEAD1", "TFAP2A", "ZBTB17", "E2F1", "SP1")

# Reorder by enriched t-TFs
heatmap_matrix <- heatmap_matrix[TF_genes_enr_order, , drop = FALSE]
significant_padj <- significant_padj[TF_genes_enr_order, , drop = FALSE]

# Convert significant padj to logical values for annotations
significant_padj_logical <- as.data.frame(significant_padj < 0.05)

# Replace NA with "NA" in the logical data frame
significant_padj_logical[is.na(significant_padj_logical)] <- "NA"

# Define the colors for the annotations, including NA as grey
annotation_colors <- list(
  padj_patients = c("TRUE" = "black", "FALSE" = "gray", "NA" = "white"),
  padj_JQ1047 = c("TRUE" = "black", "FALSE" = "gray", "NA" = "white"),
  padj_JQ1090 = c("TRUE" = "black", "FALSE" = "gray", "NA" = "white")
)

# Create the annotation object for the left side
significant_annotation <- rowAnnotation(
  padj_patients = anno_simple(as.character(significant_padj_logical[, "padj_patients"]), col = annotation_colors$padj_patients),
  padj_JQ1047 = anno_simple(as.character(significant_padj_logical[, "padj_JQ1047"]), col = annotation_colors$padj_JQ1047),
  padj_JQ1090 = anno_simple(as.character(significant_padj_logical[, "padj_JQ1090"]), col = annotation_colors$padj_JQ1090),
  annotation_name_gp = gpar(fontsize = 10)
)

# Transpose the heatmap matrix to switch axes
heatmap_matrix_transposed <- t(heatmap_matrix)
significant_padj_transposed <- t(significant_padj)

# Replace NA values with 1 in the significant_padj matrix for comparison purposes
significant_padj_transposed[is.na(significant_padj_transposed)] <- 1

# Custom function to add stars for significant values
add_stars <- function(j, i, x, y, width, height, fill) {
  if (significant_padj_transposed[i, j] < 0.001) {
    grid.text("***", x = x, y = y, gp = gpar(fontsize = 10, col = "black"))
  } else if (significant_padj_transposed[i, j] < 0.01) {
    grid.text("**", x = x, y = y, gp = gpar(fontsize = 10, col = "black"))
  } else if (significant_padj_transposed[i, j] < 0.05) {
    grid.text("*", x = x, y = y, gp = gpar(fontsize = 10, col = "black"))
  }
}

# Create and plot the heatmap with stars for significant values
heatmap <- Heatmap(
  heatmap_matrix_transposed,
  name = "log2FC_patients",
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  cell_fun = add_stars,
  heatmap_legend_param = list(
    title = "log2FC",
    at = c(-3, 0, 3),
    labels = c("-3", "0", "3"),
    legend_direction = "horizontal",
    legend_width = unit(4, "cm")
  )
)

# Save the heatmap to a PDF file
pdf("TFs_heatmap_panorama.pdf", width = 12, height = 8)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()

