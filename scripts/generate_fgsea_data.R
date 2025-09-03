##### FGSEA Analysis for Pathway Enrichment
##### Generates pathway enrichment files for Figure 6

## Load required libraries
if (!require(differential.coverage)) {
    devtools::install_github("favorov/differential.coverage") 
    library(differential.coverage)
}
library(rtracklayer)
library(dplyr)
library(plyranges)
library(data.table)
library(tidyr)
library(fgsea)
library(msigdbr)
library(ggplot2)

## Set seed for determinism
set.seed(42)

# Set working directory
setwd("~/data/")

## Get gene annotations
the_genes <- differential.coverage::get.Known.Gene.List(genome.annotation.id = "gencode19")

## Filter genes - remove duplicates on chrY, keep chrX
names_abundance_table <- table(the_genes$ensembl)
nonunique_names <- names(which(names_abundance_table > 1))
the_genes <- the_genes %>% filter((!ensembl %in% nonunique_names) | seqnames != "chrY")

## Create ensembl to gene name mapping
ensembl_to_genes <- setNames(the_genes$gene_name, the_genes$ensembl)

cat("Loaded", nrow(the_genes@elementMetadata), "unique genes for analysis\n")

## Load differential expression results
cat("Loading differential expression results...\n")
mrna_047 <- readRDS("JQ1_vs_DMSO_DE_results_O47_JQ1.rds")
mrna_090 <- readRDS("JQ1_vs_DMSO_DE_results_O90_JQ1.rds")
mrna_pts <- readRDS("T_to_N_DE_results_hg19_ensembl_77_samples.rds")

cat("Loaded DE results for:")
cat("\n- 047:", nrow(mrna_047), "genes")
cat("\n- 090:", nrow(mrna_090), "genes")
cat("\n- Patients:", nrow(mrna_pts), "genes\n")

## Retrieve Hallmark gene sets from MSigDB
cat("Loading Hallmark gene sets from MSigDB...\n")
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", collection = "H")
msigdbr_pathways <- split(x = hallmark_gene_sets$ensembl_gene, f = hallmark_gene_sets$gs_name)

cat("Loaded", length(msigdbr_pathways), "Hallmark pathways\n")

## FGSEA analysis function
fgsea_analysis <- function(rds, sample_name) {
    cat("Running FGSEA for", sample_name, "...\n")
    
    ## Create ranks vector (genes as names, log2FoldChange as values)
    ranks <- setNames(rds$log2FoldChange, rownames(rds))
    ranks <- na.omit(ranks)
    
    cat("  Using", length(ranks), "genes for ranking\n")
    
    ## Run FGSEA
    fgseaRes <- fgsea(pathways = msigdbr_pathways,
                      stats    = ranks,
                      nproc    = 1,     # Use single core to avoid crashing
                      eps      = 0.0)   # Estimate p-value accurately
    
    ## Add significance column
    fgseaRes$adjPvalue <- ifelse(fgseaRes$padj <= 0.05, "significant", "non-significant")
    
    ## Create output directory
    output_dir <- "mrna-analysis"
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }
    
    ## Generate bar plot
    cols <- c("non-significant" = "grey", "significant" = "red")
    p <- ggplot(fgseaRes, aes(reorder(.data$pathway, .data$NES), .data$NES, fill = .data$adjPvalue)) +
        geom_col() +
        scale_fill_manual(values = cols) +
        theme_bw() +
        theme(
            axis.text.y = element_text(size = 5),
            plot.title.position = "plot", 
            plot.title = element_text(size = 8, face = "bold")
        ) +
        coord_flip() +
        labs(x = "Pathway", y = "Normalized Enrichment Score", title = sample_name)
    
    ## Save plot
    pdfpath <- file.path(output_dir, paste0("fgsea-", sample_name, ".pdf"))
    ggsave(pdfpath, plot = p, width = 10, height = 8)
    cat("  Saved plot:", pdfpath, "\n")
    
    ## Save RDS file
    rdspath <- file.path(output_dir, paste0("fgsea-", sample_name, ".rds"))
    saveRDS(fgseaRes, rdspath)
    cat("  Saved RDS:", rdspath, "\n")
    
    ## Prepare and save CSV file with required columns
    fgseaRes_export <- fgseaRes %>%
        select(
            "pathway",
            "pval",
            "padj",
            "log2err",
            "ES",
            "NES",
            "size",
            "adjPvalue"
        ) %>%
        arrange(.data[["padj"]])  # Sort by adjusted p-value
    
    csvpath <- file.path(output_dir, paste0("fgsea-", sample_name, ".csv"))
    write.csv(fgseaRes_export, file = csvpath, row.names = FALSE, quote = TRUE)
    cat("  Saved CSV:", csvpath, "\n")
    
    ## Print summary
    significant_pathways <- sum(fgseaRes$padj <= 0.05, na.rm = TRUE)
    cat("  Found", significant_pathways, "significant pathways (padj <= 0.05)\n")
    
    return(fgseaRes)
}

## Run FGSEA analysis for all datasets
cat("\n=== Running FGSEA Analysis ===\n")

# Analysis for 047
result_047 <- fgsea_analysis(mrna_047, "mrna_047")

# Analysis for 090  
result_090 <- fgsea_analysis(mrna_090, "mrna_090")

# Analysis for patients
result_pts <- fgsea_analysis(mrna_pts, "mrna_pts")

## Summary of results
cat("\n=== FGSEA ANALYSIS SUMMARY ===\n")
cat("Analysis completed for 3 datasets:\n")
cat("- 047: ", sum(result_047$padj <= 0.05, na.rm = TRUE), " significant pathways\n")
cat("- 090: ", sum(result_090$padj <= 0.05, na.rm = TRUE), " significant pathways\n")
cat("- Patients: ", sum(result_pts$padj <= 0.05, na.rm = TRUE), " significant pathways\n")

cat("\nFiles generated in 'mrna-analysis/' directory:\n")
cat("- fgsea-mrna_047.csv (for Figure 6 dotplot)\n")
cat("- fgsea-mrna_090.csv (for Figure 6 dotplot)\n")
cat("- fgsea-mrna_pts.csv (for Figure 6 dotplot)\n")
cat("- Corresponding .pdf and .rds files\n")

cat("\nTo create Figure 6 dotplot, use these CSV files with fig6_pathway_dotplot.R\n")
cat("================================\n")
