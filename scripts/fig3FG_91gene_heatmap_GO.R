##### Figure 3FG: 91-gene core set analysis
##### F. Heatmap of 91-gene core set vs Tumor-SED
##### G. 91-gene core set: GO_BP Enrichment

# 1) LIBRARIES & GLOBAL SETTINGS
library(dplyr)
library(tibble)
library(differential.coverage)
library(GenomicRanges)
library(writexl)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Set working directory and paths
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

# 6) OVERLAP TABLES (patients vs each cell line)
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

# FIGURE 3F: HEATMAP FOR 91-GENE CORE SET vs Tumor-SED

# 1) Identify 91-gene core set (common down-regulated genes in both 047 & 090)
common_genes <- intersect(
  overlap_TumorSED_vs_047$gene_id,
  overlap_TumorSED_vs_090$gene_id
)

cat("91-gene core set size:", length(common_genes), "genes\n")

# 2) Build metadata data.frame with all log2FCs + gene_name
heat_df_common <- tibble(gene_id = common_genes) %>%
  left_join(
    overlap_TumorSED_vs_047 %>% select(gene_id, log2FC_047),
    by = "gene_id"
  ) %>%
  left_join(
    overlap_TumorSED_vs_090 %>% select(gene_id, log2FC_090),
    by = "gene_id"
  ) %>%
  left_join(
    genes_pat %>% select(gene_id, log2FC_patients = log2FoldChange),
    by = "gene_id"
  ) %>%
  mutate(
    gene_name = ensembl_to_name[gene_id]
  ) %>%
  column_to_rownames("gene_name")

# 3) Create matrix for heatmap, dropping any rows with NA
logfc_cols <- c("log2FC_patients", "log2FC_047", "log2FC_090")
mat_common <- as.matrix(na.omit(heat_df_common[, logfc_cols]))

# 4) Define diverging color palette
max_fc <- max(abs(mat_common), na.rm = TRUE)
col_fun <- colorRamp2(
  breaks = c(-max_fc, 0, max_fc),
  colors = c("green4", "white", "red4")
)

# 5) Plot & save heatmap PDF - Figure 3F
pdf("heatmap_91gene_core_set_vs_TumorSED.pdf", width = 6, height = 8)
Heatmap(
  mat_common,
  name            = "log2FC",
  col             = col_fun,
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  show_row_names  = TRUE,
  row_names_side  = "right",
  row_dend_side   = "left",
  row_names_gp    = gpar(fontsize = 6),
  column_title    = "91-gene core set vs Tumor-SED",
  column_names_gp = gpar(fontsize = 8)
)
dev.off()
cat("Saved Figure 3F heatmap → heatmap_91gene_core_set_vs_TumorSED.pdf\n")

# 6) Export metadata list with all columns to Excel & CSV
export_df <- heat_df_common %>%
  rownames_to_column("gene_name") %>%
  select(
    gene_name,
    gene_id,
    log2FC_patients,
    log2FC_047,
    log2FC_090
  )

# Excel
write_xlsx(
  export_df,
  path = "91gene_core_set_metadata.xlsx"
)
# CSV
write.csv(
  export_df,
  "91gene_core_set_metadata.csv",
  row.names = FALSE
)
cat("Exported 91-gene core set metadata →\n",
    "- 91gene_core_set_metadata.xlsx\n",
    "- 91gene_core_set_metadata.csv\n")

# FIGURE 3G: GO_BP ENRICHMENT FOR 91-GENE CORE SET

# 1) Define 91-gene core set (Ensembl IDs)
core_set_ensembl <- intersect(
  overlap_TumorSED_vs_047$gene_id,
  overlap_TumorSED_vs_090$gene_id
)
cat("91-gene core set for GO enrichment:", length(core_set_ensembl), "genes\n")

# 2) Clean Ensembl IDs (remove version suffix) & map to Entrez
clean_ids <- sub("\\..*$", "", core_set_ensembl)
map_df   <- bitr(
  clean_ids,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
entrez_ids <- unique(map_df$ENTREZID)
cat("Mapped", length(entrez_ids), "genes to Entrez IDs\n")

# Identify unmapped Ensembl IDs for manual inspection
unmapped <- setdiff(clean_ids, map_df$ENSEMBL)
cat("Number of unmapped Ensembl IDs:", length(unmapped), "\n")
if(length(unmapped) > 0) {
  write.table(
    unmapped,
    file = "unmapped_ensembl_ids_91gene_core_set.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = "Ensembl_ID"
  )
}

# 3) Run GO–BP enrichment
ego_core_set_bp <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",    # input type
  ont           = "BP",          # Biological Process
  pvalueCutoff  = 0.05,          # raw p-value ≤ 0.05
  pAdjustMethod = "BH",          # Benjamini–Hochberg correction
  qvalueCutoff  = 0.05,          # FDR ≤ 0.05
  readable      = TRUE           # convert Entrez → gene symbols
)

# Inspect top terms
cat("Top GO_BP terms for 91-gene core set:\n")
head(ego_core_set_bp, 10)

# 4) Visualize as dot plot - Figure 3G
pdf("GO_BP_91gene_core_set.pdf", width = 8, height = 6)
print(
  dotplot(
    ego_core_set_bp,
    showCategory = 10,
    title = "91-gene core set: GO_BP Enrichment"
  )
)
dev.off()
cat("Saved Figure 3G dotplot → GO_BP_91gene_core_set.pdf\n")

# 5) Export results to CSV
write.csv(
  as.data.frame(ego_core_set_bp),
  "GO_BP_91gene_core_set.csv",
  row.names = FALSE
)
cat("Exported GO_BP results → GO_BP_91gene_core_set.csv\n")

# SUMMARY
cat("\n=== FIGURE 3FG SUMMARY ===\n")
cat("Figure 3F: Heatmap of 91-gene core set vs Tumor-SED\n")
cat("  - File: heatmap_91gene_core_set_vs_TumorSED.pdf\n")
cat("  - Genes:", nrow(mat_common), "\n")
cat("Figure 3G: GO_BP Enrichment for 91-gene core set\n")
cat("  - File: GO_BP_91gene_core_set.pdf\n")
cat("  - Enriched terms:", nrow(ego_core_set_bp), "\n")
cat("========================\n")
