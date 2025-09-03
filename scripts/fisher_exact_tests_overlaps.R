##### Fisher's Exact Tests for Gene Set Overlaps
##### Generates fisher_results.tsv

# 1) LIBRARIES & GLOBAL SETTINGS
library(dplyr)
library(tibble)
library(differential.coverage)
library(GenomicRanges)

# Set working directory
setwd("~/data/")

# 2) GENE ANNOTATION: load & filter
genes_all <- differential.coverage::get.Known.Gene.List("gencode19")
dups      <- names(which(table(genes_all$ensembl) > 1))
genes     <- genes_all[
  !(genes_all$ensembl %in% dups) &
    seqnames(genes_all) != "chrY"
]
cat("Kept", nrow(genes@elementMetadata), "unique genes.\n")

ensembl_to_name <- setNames(genes$gene_name, genes$ensembl)

attach_gene_names <- function(df) {
  df %>%
    as.data.frame() %>%
    rownames_to_column("ensembl") %>%
    left_join(
      tibble(ensembl = names(ensembl_to_name),
             gene_name = ensembl_to_name),
      by = "ensembl"
    )
}

# 3) READ & ANNOTATE DE RESULTS
mrna_O47 <- attach_gene_names(
  readRDS("JQ1_vs_DMSO_DE_results_O47_JQ1.rds")
)
mrna_O90 <- attach_gene_names(
  readRDS("JQ1_vs_DMSO_DE_results_O90_JQ1.rds")
)
mrna_pts <- attach_gene_names(
  readRDS("T_to_N_DE_results_hg19_ensembl_77_samples.rds")
)

# 4) PATIENT "SED" DATA: load & filter
genes_pat <- readRDS("fold_and_distance_77_samples_SED.rds") %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(SED_type == "tumor", padj < 0.05, log2FoldChange > 0)
cat("Tumor-SED upregulated genes:", nrow(genes_pat), "\n")

# 5) CELL-LINE LISTS: load
g047 <- readRDS("down_047.rds") %>%
  { df <- as.data.frame(.);
  if (!"gene_id" %in% names(df)) df$gene_id <- rownames(df);
  rownames(df) <- NULL; df
  }

g090 <- readRDS("down_090.rds") %>%
  { df <- as.data.frame(.);
  if (!"gene_id" %in% names(df)) df$gene_id <- rownames(df);
  rownames(df) <- NULL; df
  }

cat("UM-SCC-047 down-regulated genes:", nrow(g047), "\n")
cat("UPCI-SCC-090 down-regulated genes:", nrow(g090), "\n")

# 6) OVERLAP CALCULATIONS

# Calculate overlaps between cell lines and Tumor-SED
overlap_TumorSED_vs_047 <- genes_pat %>%
  dplyr::select(gene_id, log2FC_patients = log2FoldChange) %>%
  dplyr::inner_join(
    g047 %>%
      dplyr::select(gene_id, log2FC_047 = log2FoldChange),
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    gene_name = ensembl_to_name[gene_id]
  )

overlap_TumorSED_vs_090 <- genes_pat %>%
  dplyr::select(gene_id, log2FC_patients = log2FoldChange) %>%
  dplyr::inner_join(
    g090 %>%
      dplyr::select(gene_id, log2FC_090 = log2FoldChange),
    by = "gene_id"
  ) %>%
  dplyr::mutate(
    gene_name = ensembl_to_name[gene_id]
  )

# 7) FISHER'S EXACT TESTS FOR OVERLAPS

# Define parameters for Fisher's exact tests
total_genes <- length(unique(mrna_pts$ensembl))
n_red       <- nrow(genes_pat)                    # Tumor-SED upregulated
n_green_047 <- nrow(g047)                        # UM-SCC-047 downregulated
n_green_090 <- nrow(g090)                        # UPCI-SCC-090 downregulated

# Calculate overlap sizes
ov_047    <- nrow(overlap_TumorSED_vs_047)       # UM-SCC-047 vs Tumor-SED overlap
ov_090    <- nrow(overlap_TumorSED_vs_090)       # UPCI-SCC-090 vs Tumor-SED overlap
ov_common <- length(                             # Common genes (047 & 090) vs Tumor-SED
  intersect(intersect(g047$gene_id, g090$gene_id), genes_pat$gene_id)
)

cat("\n=== OVERLAP SUMMARY ===\n")
cat("Total genes in background:", total_genes, "\n")
cat("Tumor-SED upregulated genes:", n_red, "\n")
cat("UM-SCC-047 downregulated genes:", n_green_047, "\n")
cat("UPCI-SCC-090 downregulated genes:", n_green_090, "\n")
cat("UM-SCC-047 vs Tumor-SED overlap:", ov_047, "\n")
cat("UPCI-SCC-090 vs Tumor-SED overlap:", ov_090, "\n")
cat("Common genes (047&090) vs Tumor-SED overlap:", ov_common, "\n")
cat("========================\n")

# Function to perform Fisher's exact test
run_fisher_table <- function(overlap, down_size, label) {
  # Create 2x2 contingency table
  # Rows: in/not in downregulated set
  # Columns: in/not in upregulated set
  m <- matrix(c(
    overlap,                                    # overlap
    n_red        - overlap,                     # upregulated but not downregulated
    down_size    - overlap,                     # downregulated but not upregulated
    total_genes - n_red - down_size + overlap   # neither upregulated nor downregulated
  ), nrow = 2)
  
  # Perform Fisher's exact test
  res <- fisher.test(m, conf.level = 0.95)
  
  # Extract results
  ci   <- signif(res$conf.int, 3)
  or   <- signif(res$estimate, 3)
  pval <- signif(res$p.value, 3)
  
  # Print contingency table for verification
  cat("\n", label, ":\n")
  cat("Contingency table:\n")
  rownames(m) <- c("In downregulated", "Not in downregulated")
  colnames(m) <- c("In upregulated", "Not in upregulated")
  print(m)
  cat("Odds Ratio:", or, "\n")
  cat("95% CI: [", ci[1], ",", ci[2], "]\n")
  cat("P-value:", pval, "\n")
  
  # Return results as data frame
  data.frame(
    Comparison = label,
    Overlap = overlap,
    Down_Size = down_size,
    Up_Size = n_red,
    Total_Background = total_genes,
    Odds_Ratio = or,
    CI_Lower = ci[1],
    CI_Upper = ci[2],
    P_Value = pval,
    stringsAsFactors = FALSE
  )
}

# Perform Fisher's exact tests for each comparison
cat("\n=== FISHER'S EXACT TESTS ===\n")

tab1 <- run_fisher_table(ov_047, n_green_047, "UM-SCC-047 vs Tumor-SED")
tab2 <- run_fisher_table(ov_090, n_green_090, "UPCI-SCC-090 vs Tumor-SED")
tab3 <- run_fisher_table(
  ov_common, 
  length(intersect(g047$gene_id, g090$gene_id)),
  "91-gene core set vs Tumor-SED"
)

# Combine all results
results_table <- rbind(tab1, tab2, tab3)

# Display results
cat("\n=== FINAL RESULTS ===\n")
print(results_table)

# 8) EXPORT RESULTS

# Export to TSV file
write.table(
  results_table, 
  file = "fisher_results.tsv", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Also export to CSV for easier viewing
write.csv(
  results_table,
  file = "fisher_results.csv",
  row.names = FALSE
)

cat("\n=== FILES EXPORTED ===\n")
cat("Fisher's exact test results saved to:\n")
cat("- fisher_results.tsv\n")
cat("- fisher_results.csv\n")

# 9) SUMMARY STATISTICS

cat("\n=== ENRICHMENT SUMMARY ===\n")
for(i in 1:nrow(results_table)) {
  row <- results_table[i, ]
  enrichment <- ifelse(row$Odds_Ratio > 1, "ENRICHED", "DEPLETED")
  significance <- ifelse(row$P_Value < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")
  
  cat(sprintf("%s: %s (OR=%.2f, p=%.3e) - %s\n", 
              row$Comparison, 
              enrichment, 
              row$Odds_Ratio, 
              row$P_Value, 
              significance))
}
cat("============================\n")
