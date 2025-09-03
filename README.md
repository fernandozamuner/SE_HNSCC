# SE_HNSCC

Analysis of super-enhancer domains in head and neck squamous cell carcinoma (HNSCC).

## Overview

This repository contains R scripts for analyzing super-enhancer (SE) domains and their associated gene expression changes in HNSCC samples. The analysis includes differential expression, pathway enrichment, and visualization of SE-regulated genes.

## Key Features
- **Super-Enhancer Analysis**: Identify and characterize SE domains and their target genes
- **Differential Expression Analysis**: Compare gene expression between treatment conditions (JQ1 vs DMSO) and tumor vs normal samples
- **Transcription Factor Analysis**: Identify key transcription factors regulating SE-associated genes
- **Pathway Enrichment**: GSEA analysis using Hallmark gene sets
- **Data Visualization**: Generate publication-ready figures including volcano plots, heatmaps, and pathway dotplots
- **Statistical Analysis**: Fisher's exact tests, contingency analyses, and transcription factor enrichment
- **Multi-sample Integration**: Analysis of 77 HNSCC samples plus cell line experiments

## Requirements

- R (≥4.0)
- Required packages: DESeq2, fgsea, ggplot2, dplyr, ComplexHeatmap, clusterProfiler, msigdbr

## Script Execution Order

### Phase 1: Data Preparation and Processing
1. `1_prepare_and_organise_SEs.R` - Process and organize super-enhancer (SE) data from H3K27ac scores
2. `2_prepare_domain_regions.R` - Prepare SE domain regions for analysis
3. `3_salmon2ensembletxi.r` - Process RNA-seq data using Salmon and tximport (77 samples)
4. `1_1_salmon2ensembletxiJQ1.r` - Process JQ1 treatment RNA-seq data
5. `4_deseq_gencode19_ensembl_77_samples.R` - Differential expression analysis (tumor vs normal, 77 samples)
6. `1_2_deseq_gencode19_JQ1_O47.R` - Differential expression analysis (JQ1 vs DMSO, O47 cell line)
7. `1_3_deseq_gencode19_JQ1_O90.R` - Differential expression analysis (JQ1 vs DMSO, O90 cell line)

### Phase 2: SE Domain Analysis
8. `5_prepare_fold_and_distance_object_77_samples_SED.r` - Prepare fold change and distance objects for SE domains
9. `6_SE_domains_txi.R` - Process SE domain transcriptomic data
10. `7_SE_domains_deseq2_res.R` - DESeq2 analysis for SE domains
11. `make_domains_fasta.r` - Generate FASTA files for SE domains

### Phase 3: Transcription Factor Analysis
12. `ParseCistromeCounts.R` - Parse Cistrome database counts
13. `SElogodds.R` - Calculate super-enhancer log odds

### Phase 4: Pathway and Functional Analysis
14. `generate_fgsea_data.R` - Generate FGSEA pathway enrichment results
15. `fisher_exact_tests_overlaps.R` - Statistical analysis of gene set overlaps
16. `CompContingencies.R` - Contingency table analysis
17. `ContingenciesResults.R` - Process contingency analysis results

### Phase 5: Figure Generation
18. `fig2AB_TF_odds_pval.R` - Generate Figure 2A-B (transcription factor analysis)
19. `fig2C_TFs_heatmap_panorama.R` - Generate Figure 2C (TF heatmap)
20. `fig3ABCDE_volcano_plots.R` - Generate Figure 3A-E (volcano plots)
21. `fig3FG_91gene_heatmap_GO.R` - Generate Figure 3F-G (heatmap and GO analysis)
22. `fig4AB_ridge_plots.R` - Generate Figure 4A-B (ridge plots)
23. `fig5AB_volcanoSED-eRNA_hexplot_eRNAvsmRNA.R` - Generate Figure 5A-B (volcano and correlation plots)
24. `fig6_pathway_dotplot.R` - Generate Figure 6 (pathway dotplots)
25. `figS1_scatterplot_cell_lines.R` - Generate Supplementary Figure 1
