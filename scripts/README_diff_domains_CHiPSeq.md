# Domain-Level ChIP-seq Pipeline

This script (`diff_domains_CHiPSeq_EDGER.R`) quantifies ChIP-seq signal across pre-defined promoter, enhancer, and super-enhancer domains and compares tumor versus matching normal samples. It produces per-domain coverage tables, model outputs, and annotated CSV files that can be fed into the transcription-factor enrichment workflow.

## What the script does

1. **Loads domain geometries** from `rds/domains/*_domains_ranges.rds`. If the ranges are missing, it invokes `Scripts/domains/prepare_domain_regions.R` to recreate them.
2. **Counts aligned reads** by running `featureCounts` on the four H3K27ac replicates (two normal, two tumor) together with their matching input libraries.
3. **Fits a generalized linear model** that contrasts tumor and normal signal while accounting for input/background tracks.
4. **Exports summaries** including:
   - full model tables with coordinates (`results/domains/*_domains_*_on_CHiPSeq.csv`),
   - significant domains (`results/domains/*_diff_domains_*_on_CHiPSeq.csv`),
   - raw counts and normalized CPM matrices (`results/domains/*_domains_*_coverage_and_CPM.csv`),
   - serialized R objects in `rds/domains/` for downstream reuse.

## Running the workflow

1. Ensure the domain range files exist under `rds/domains/`. If they do not, provide the inputs required by `Scripts/domains/prepare_domain_regions.R`.
2. Update `bam_folder` inside the script to point to the directory containing the four H3K27ac BAM files and their 2% input controls. The current path is a placeholder.
3. Adjust the `samples` vector if your replicate names differ from `3/4/5/6`.
4. From the repository root, run:

   ```bash
   Rscript Scripts/diff_domains_CHiPSeq_EDGER.R
   ```

5. Inspect the `results/domains/` directory for the refreshed CSV outputs. These are the tables expected by the TF overlap pipeline (e.g., `tf_analysis/source_data/*_domains_*_coverage_and_CPM.csv`).

## Notes

- The script expects eight BAM files: four H3K27ac enrichments followed by their four input controls, all aligned to hg19.
- Intermediate count objects are cached in `rds/domains/*_domains_counts_for_CHiPSeq.rds` to avoid re-running `featureCounts` unnecessarily.
- The contrast vector (`c(-1, 1, 1, -1)`) represents a tumor-minus-normal comparison while subtracting the corresponding input difference; modify it if your experimental design changes.

## Cistrome peak set

The TF overlap scripts expect a set of hg19 peak calls under `tf_analysis/source_data/hg19_cismotifs/`. The upstream dataset contains genome-wide TF binding regions as BED files (`<TF>.<reliability>.bed`), organized by transcription factor (UniProt identifier) and reliability class (A–D). The motif subset (`hg19_cismotifs`) adds a fourth column with the maximal HOCOMOCO motif score for each segment (−log10 scale; score ≥ 4 corresponds to a motif P-value ≤ 0.0001). Keep that folder accessible (or symlinked) so the pipeline can merge TF peaks before computing enrichment statistics.
