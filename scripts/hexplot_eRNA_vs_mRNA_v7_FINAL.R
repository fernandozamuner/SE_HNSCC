##Fernando Zamuner
##Sep_2024

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

# Clear the environment
rm(list=ls())  # Removes all objects from the current environment

# Set working directory
setwd("~/data/")

# Read the eRNA data
erna <- as.data.frame(readRDS("SE_domains-deseq-res.rds"))  # Load eRNA data from RDS file and convert to data frame
SE_ranges <- as.data.frame(readRDS("SE_domains_ranges.rds")) %>%
  filter(seqnames != "chrY" & seqnames != "chrM")  # Exclude data from chromosome Y and mitochondria

# Add SED_ID to each dataset
erna <- tibble::rownames_to_column(erna, "SED_ID")  # Add rownames as a new column (SED_ID) in the erna data frame
SE_ranges <- tibble::rownames_to_column(SE_ranges, "SED_ID")  # Same for SE_ranges

# Join the eRNA and SE_ranges datasets by SED_ID
combined_df <- inner_join(erna, SE_ranges, by = "SED_ID")  # Merge the two datasets based on SED_ID
colnames(combined_df)[colnames(combined_df) == "log2FoldChange"] <- "lfc_erna"  # Rename the log2FoldChange column

# Load the gene annotation data
the_genes <- differential.coverage::get.Known.Gene.List(genome.annotation.id = "gencode19")  # Get gene annotation data

# Exclude genes located on chromosome Y and mitochondria while keeping GRanges format
the_genes_no_y_no_m <- the_genes[!(seqnames(the_genes) %in% c("chrY", "chrM")), ]  # Filter out chrY and chrM genes

# Load RNAseq data for 77 samples
mrna <- readRDS('T_to_N_DE_results_hg19_ensembl_77_samples.rds') %>%
  as.data.frame %>%
  rownames_to_column(var = "ensembl")  # Convert RNAseq data into a data frame and add ensembl IDs as a column

# Map ENSG IDs to gene names using the gene annotations and filter missing values
mrna <- inner_join(as.data.frame(the_genes_no_y_no_m), mrna, by = "ensembl") %>%
  filter(!is.na(gene_name) & !is.na(log2FoldChange) & !is.na(padj))  # Remove rows with missing gene names, log2FoldChange, or padj

# Handle genes with multiple entries by keeping the one with the lowest p-adjusted value
gene_name_counts <- table(mrna$gene_name)  # Count occurrences of each gene name
nonunique_genes <- names(gene_name_counts[gene_name_counts > 1])  # Identify genes with multiple entries
for (dname in nonunique_genes) {
  myrows <- which(mrna$gene_name == dname)  # Get the row indices for each duplicated gene
  winner <- which.min(mrna[myrows, "padj"])  # Find the row with the smallest padj value
  nonwinnerrows <- myrows[-winner]  # Get the non-winning rows
  mrna <- mrna[-nonwinnerrows, ]  # Remove the non-winning rows
}

# Filter the gene list to include only the genes that appear in the mRNA data
filtered_genes <- the_genes_no_y_no_m[mcols(the_genes_no_y_no_m)$gene_name %in% mrna$gene_name]

# Calculate the midpoint of each SED region
SED_gr <- GRanges(
  seqnames = combined_df$seqnames, 
  ranges = IRanges(start = combined_df$start, end = combined_df$end), 
  SED_ID = combined_df$SED_ID  # Create GRanges object for SEDs
)
SED_mids <- SED_gr
start(SED_mids) <- (start(SED_gr) + end(SED_gr)) / 2  # Calculate the midpoint
width(SED_mids) <- 0  # Set width to 0, since we're only interested in the midpoint

# Find the nearest gene for each SED using the filtered gene list
nearest_genes <- nearest(SED_mids, filtered_genes)  # Get the indices of the nearest genes
distances <- distance(SED_mids, filtered_genes[nearest_genes])  # Calculate the distances to the nearest genes
TSS_in_pairs <- filtered_genes[nearest_genes]  # Get the nearest gene objects
SED_in_pairs <- SED_mids  # Keep the corresponding SEDs

# Create a dataframe with the nearest gene information for each SED
nearest_gene_info <- data.frame(
  SED_ID = SED_in_pairs$SED_ID,  # SED ID
  nearest_gene = mcols(TSS_in_pairs)$gene_name,  # Nearest gene name
  nearest_gene_start = start(TSS_in_pairs),  # Start position of the nearest gene
  nearest_gene_end = end(TSS_in_pairs),  # End position of the nearest gene
  distance_to_nearest_gene = distances  # Distance to the nearest gene
)

# Combine the nearest gene information with the original eRNA data
eRNA_expression <- combined_df %>%
  left_join(nearest_gene_info, by = "SED_ID")  # Merge eRNA and nearest gene information

# Ensure combined_expression has the same rows as eRNA_expression by performing a left join with mRNA data
combined_expression <- eRNA_expression %>%
  left_join(mrna, by = c("nearest_gene" = "gene_name"))  # Join with mRNA data using the gene names

# Rename the mRNA log2FoldChange column for clarity
colnames(combined_expression)[colnames(combined_expression) == "log2FoldChange"] <- "lfc_mrna"  # Rename mRNA fold change column

# Save the combined dataframe for further analysis
saveRDS(combined_expression, "combined_eRNA_mRNA_expression.rds")  # Save the final dataset

# Calculate the correlation between eRNA and mRNA log2 fold changes
correlation_test <- cor.test(combined_expression$lfc_mrna, combined_expression$lfc_erna, method = "kendall")  # Perform Kendall correlation test

# Print correlation test results
print(correlation_test)  # Output the result of the correlation test

# Create a hexbin plot to visualize the relationship between eRNA and mRNA fold changes
plot_linear <- ggplot(combined_expression, aes(x = lfc_mrna, y = lfc_erna)) + 
  geom_hline(yintercept = 0, colour = 'grey15') +  # Add horizontal line at y=0
  geom_vline(xintercept = 0, colour = 'grey15') +  # Add vertical line at x=0
  geom_hex(aes(fill = factor(source), colour = factor(source), alpha = ..count..), bins = 50, linewidth = 0.25) +  # Create hexbin plot
  geom_smooth(method = "lm", colour = "black", linetype = "solid") +  # Add a linear regression line
  scale_colour_manual('source', values = c(
    'UPPP' = 'blue',  # Color UPPP samples in blue
    'tumor' = 'red',  # Color tumor samples in red
    'mixed' = 'gray20'  # Color mixed samples in gray
  ), aesthetics = c('colour', 'fill')) +  # Customize color and fill scale
  scale_alpha_continuous(name = "Count", range = c(0.4, 1)) +  # Adjust transparency based on count
  xlab('mRNA log2 Fold Change') +  # Label x-axis
  ylab('eRNA log2 Fold Change') +  # Label y-axis
  ggtitle('eRNA vs mRNA') +  # Add plot title
  theme_bw(base_size = 12) +  # Use a clean theme with base size 12
  labs(caption = paste0('n = ', nrow(combined_expression)))  # Add a caption showing the sample size

# Display the plot
plot_linear  # Render the plot

# Save the plot as a PDF file
ggsave("hex-SED-erna-vs-mrna3.pdf", plot_linear, width = 8, height = 8)  # Save the plot to a PDF file
