####Fernando Zamuner

# Load required packages
library(rtracklayer)
library(dplyr)
library(plyranges)
library(tidyr)
library(EnhancedVolcano)

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/PROJECTS/Other_Projects/SE-HNSCC-Daria/Codes_Ridge_plots/hexagon_plot")

# Read the eRNA data
erna <- as.data.frame(readRDS("combined_eRNA_mRNA_expression.rds"))
head(erna)

# Filter the erna data for lfc_erna greater than 1.5 or less than -1.5, and padj.x < 0.05
erna_subset <- erna %>%
  filter(abs(lfc_erna) > 1.5 & padj.x < 0.05)

# View the subset
head(erna_subset)

# Optionally, save the subset to a CSV or RDS file
write.csv(erna_subset, "erna_subset.csv", row.names = FALSE)
saveRDS(erna_subset, "erna_subset.rds")

# Volcano plot with adjusted arrows and label size
plot <- EnhancedVolcano(
  erna,
  lab = ifelse(erna$padj.x <= 0.05 & abs(erna$lfc_erna) >= 1.5, erna$nearest_gene, NA), # Only label genes meeting conditions
  x = "lfc_erna",
  y = "padj.x",
  selectLab = c(
    # Upregulated genes
    'EGFR',
    'TP63',
    'CLDN1',
    'SLC2A1',
    'ITGA2',
    'PLAU',
    'ABCC5',
    'ALCAM',
    'EGLN3',
    'FAM83B',
    'SERPINB11',
    'JUP',
    'CAMK1G',
    'FMNL2',
    'RBBP8',
    # Downregulated genes
    'MEF2C',
    'IL16',
    'EMP1',
    'POU2F2',
    'ECM1',
    'TCF7',
    'SPINK5',
    'KRT78',
    'PAX5'),
  title = "eRNA (SED)",  # Title
  subtitle = NULL,
  pCutoff = 0.05,
  pCutoffCol = "padj.x",
  FCcutoff = 1.5,
  pointSize = 2.3,
  labSize = 5.0,                # Increased label size
  labCol = 'black',
  boxedLabels = TRUE,
  legendPosition = "none",
  drawConnectors = TRUE,         # Ensures connectors are drawn
  widthConnectors = 1.0,         # Increased width of connectors
  colConnectors = 'black',
  col = c("black", "black", "black", "red"),
  maxoverlapsConnectors = 20,    # Increased this to reduce overlap
  max.overlaps = 15,             # Reduce overlap in labels
  xlab = bquote(~Log[2]~ "Fold Change")
)

plot

# Save the plot
ggsave("volcano_eRNA_logFC.pdf", plot, width = 8, height = 8)
