##### Fernando Zamuner
##### Figure 5AB: eRNA Analysis
##### A. Volcano plot of eRNA (SED) expression
##### B. Hexbin plot comparing eRNA vs mRNA expression

# Load necessary libraries
library(optparse)
library(dplyr)    
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(tibble)
library(ggthemes)
library(differential.coverage)
library(purrr)
library(rtracklayer)
library(plyranges)
library(tidyr)
library(EnhancedVolcano)

# Clear the environment
rm(list=ls())

# Set working directory
setwd("~/data/")

# 1) DATA PROCESSING AND PREPARATION

# Read the eRNA data
erna <- as.data.frame(readRDS("SE_domains-deseq-res.rds"))
SE_ranges <- as.data.frame(readRDS("SE_domains_ranges.rds")) %>%
  filter(seqnames != "chrY" & seqnames != "chrM")

# Add SED_ID to each dataset
erna <- tibble::rownames_to_column(erna, "SED_ID")
SE_ranges <- tibble::rownames_to_column(SE_ranges, "SED_ID")

# Join the eRNA and SE_ranges datasets by SED_ID
combined_df <- inner_join(erna, SE_ranges, by = "SED_ID")
colnames(combined_df)[colnames(combined_df) == "log2FoldChange"] <- "lfc_erna"

# Load the gene annotation data
the_genes <- differential.coverage::get.Known.Gene.List(genome.annotation.id = "gencode19")

# Exclude genes located on chromosome Y and mitochondria while keeping GRanges format
the_genes_no_y_no_m <- the_genes[!(seqnames(the_genes) %in% c("chrY", "chrM")), ]

# Load RNAseq data for 77 samples
mrna <- readRDS('T_to_N_DE_results_hg19_ensembl_77_samples.rds') %>%
  as.data.frame %>%
  rownames_to_column(var = "ensembl")

# Map ENSG IDs to gene names using the gene annotations and filter missing values
mrna <- inner_join(as.data.frame(the_genes_no_y_no_m), mrna, by = "ensembl") %>%
  filter(!is.na(gene_name) & !is.na(log2FoldChange) & !is.na(padj))

# Handle genes with multiple entries by keeping the one with the lowest p-adjusted value
gene_name_counts <- table(mrna$gene_name)
nonunique_genes <- names(gene_name_counts[gene_name_counts > 1])
for (dname in nonunique_genes) {
  myrows <- which(mrna$gene_name == dname)
  winner <- which.min(mrna[myrows, "padj"])
  nonwinnerrows <- myrows[-winner]
  mrna <- mrna[-nonwinnerrows, ]
}

# Filter the gene list to include only the genes that appear in the mRNA data
filtered_genes <- the_genes_no_y_no_m[mcols(the_genes_no_y_no_m)$gene_name %in% mrna$gene_name]

# Calculate the midpoint of each SED region
SED_gr <- GRanges(
  seqnames = combined_df$seqnames, 
  ranges = IRanges(start = combined_df$start, end = combined_df$end), 
  SED_ID = combined_df$SED_ID
)
SED_mids <- SED_gr
start(SED_mids) <- (start(SED_gr) + end(SED_gr)) / 2
width(SED_mids) <- 0

# Find the nearest gene for each SED using the filtered gene list
nearest_genes <- nearest(SED_mids, filtered_genes)
distances <- distance(SED_mids, filtered_genes[nearest_genes])
TSS_in_pairs <- filtered_genes[nearest_genes]
SED_in_pairs <- SED_mids

# Create a dataframe with the nearest gene information for each SED
nearest_gene_info <- data.frame(
  SED_ID = SED_in_pairs$SED_ID,
  nearest_gene = mcols(TSS_in_pairs)$gene_name,
  nearest_gene_start = start(TSS_in_pairs),
  nearest_gene_end = end(TSS_in_pairs),
  distance_to_nearest_gene = distances
)

# Combine the nearest gene information with the original eRNA data
eRNA_expression <- combined_df %>%
  left_join(nearest_gene_info, by = "SED_ID")

# Ensure combined_expression has the same rows as eRNA_expression by performing a left join with mRNA data
combined_expression <- eRNA_expression %>%
  left_join(mrna, by = c("nearest_gene" = "gene_name"))

# Rename the mRNA log2FoldChange column for clarity
colnames(combined_expression)[colnames(combined_expression) == "log2FoldChange"] <- "lfc_mrna"

# Save the combined dataframe for further analysis
saveRDS(combined_expression, "combined_eRNA_mRNA_expression.rds")

cat("Data processing completed. Combined dataset saved.\n")
cat("Dataset dimensions:", nrow(combined_expression), "rows,", ncol(combined_expression), "columns\n")

# 2) FIGURE 5A: VOLCANO PLOT - eRNA (SED) Expression

# Filter the erna data for significant changes
erna_subset <- combined_expression %>%
  filter(abs(lfc_erna) > 1.5 & padj.x < 0.05)

cat("Significant eRNA changes (|FC| > 1.5, padj < 0.05):", nrow(erna_subset), "\n")

# Save the subset
write.csv(erna_subset, "erna_subset.csv", row.names = FALSE)
saveRDS(erna_subset, "erna_subset.rds")

# Create volcano plot
volcano_plot <- EnhancedVolcano(
  combined_expression,
  lab = ifelse(combined_expression$padj.x <= 0.05 & abs(combined_expression$lfc_erna) >= 1.5, 
               combined_expression$nearest_gene, NA),
  x = "lfc_erna",
  y = "padj.x",
  selectLab = c(
    # Upregulated genes
    'EGFR', 'TP63', 'CLDN1', 'SLC2A1', 'ITGA2', 'PLAU', 'ABCC5', 'ALCAM',
    'EGLN3', 'FAM83B', 'SERPINB11', 'JUP', 'CAMK1G', 'FMNL2', 'RBBP8',
    # Downregulated genes
    'MEF2C', 'IL16', 'EMP1', 'POU2F2', 'ECM1', 'TCF7', 'SPINK5', 'KRT78', 'PAX5'
  ),
  title = "eRNA (SED)",
  subtitle = NULL,
  pCutoff = 0.05,
  pCutoffCol = "padj.x",
  FCcutoff = 1.5,
  pointSize = 2.3,
  labSize = 5.0,
  labCol = 'black',
  boxedLabels = TRUE,
  legendPosition = "none",
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = 'black',
  col = c("black", "black", "black", "red"),
  maxoverlapsConnectors = 20,
  max.overlaps = 15,
  xlab = bquote(~Log[2]~ "Fold Change")
)

print(volcano_plot)

# Save Figure 5A
ggsave("fig5A_volcano_eRNA_logFC.pdf", volcano_plot, width = 8, height = 8)
cat("Figure 5A (volcano plot) saved as: fig5A_volcano_eRNA_logFC.pdf\n")

# 3) FIGURE 5B: HEXBIN PLOT - eRNA vs mRNA Expression

# Calculate the correlation between eRNA and mRNA log2 fold changes
correlation_test <- cor.test(combined_expression$lfc_mrna, combined_expression$lfc_erna, 
                           method = "kendall", use = "complete.obs")

cat("Correlation analysis (Kendall's tau):\n")
cat("Correlation coefficient:", round(correlation_test$estimate, 4), "\n")
cat("P-value:", format(correlation_test$p.value, scientific = TRUE), "\n")

# Create hexbin plot
hexplot <- ggplot(combined_expression, aes(x = lfc_mrna, y = lfc_erna)) + 
  geom_hline(yintercept = 0, colour = 'grey15') +
  geom_vline(xintercept = 0, colour = 'grey15') +
  geom_hex(aes(fill = factor(source), colour = factor(source), alpha = ..count..), 
           bins = 50, linewidth = 0.25) +
  geom_smooth(method = "lm", colour = "black", linetype = "solid") +
  scale_colour_manual('source', values = c(
    'UPPP' = 'blue',
    'tumor' = 'red',
    'mixed' = 'gray20'
  ), aesthetics = c('colour', 'fill')) +
  scale_alpha_continuous(name = "Count", range = c(0.4, 1)) +
  xlab('mRNA log2 Fold Change') +
  ylab('eRNA log2 Fold Change') +
  ggtitle('eRNA vs mRNA') +
  theme_bw(base_size = 12) +
  labs(caption = paste0('n = ', nrow(combined_expression), 
                       ', τ = ', round(correlation_test$estimate, 3),
                       ', p = ', format(correlation_test$p.value, scientific = TRUE, digits = 3)))

print(hexplot)

# Save Figure 5B
ggsave("fig5B_hex_SED_erna_vs_mrna.pdf", hexplot, width = 8, height = 8)
cat("Figure 5B (hexbin plot) saved as: fig5B_hex_SED_erna_vs_mrna.pdf\n")

# 4) SUMMARY STATISTICS

cat("\n=== FIGURE 5AB SUMMARY ===\n")
cat("Total SEDs analyzed:", nrow(combined_expression), "\n")
cat("SEDs with significant eRNA changes:", nrow(erna_subset), "\n")
cat("Correlation (eRNA vs mRNA):", round(correlation_test$estimate, 4), "\n")
cat("P-value:", format(correlation_test$p.value, scientific = TRUE), "\n")

# Source distribution
source_summary <- table(combined_expression$source)
cat("Source distribution:\n")
print(source_summary)

cat("\nFiles generated:\n")
cat("- combined_eRNA_mRNA_expression.rds (processed data)\n")
cat("- erna_subset.csv and erna_subset.rds (significant changes)\n")
cat("- fig5A_volcano_eRNA_logFC.pdf (Figure 5A)\n")
cat("- fig5B_hex_SED_erna_vs_mrna.pdf (Figure 5B)\n")
cat("========================\n")
