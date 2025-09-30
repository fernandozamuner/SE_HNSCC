#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
  library(Rsubread)
  library(dplyr)
  library(plyranges)
  library(stringr)
  library(edgeR)
  library(statmod)
})

# Ensure domain geometries are prepared
domains_ready <- file.exists('rds/domains/SE_domains_ranges.rds') &
  file.exists('rds/domains/E_domains_ranges.rds') &
  file.exists('rds/domains/P_domains_ranges.rds')

if (!domains_ready) {
  message('Preparing domain ranges...')
  source('Scripts/domains/prepare_domain_regions.R')
}
rm(list = c('domains_ready'))

run_domain_analysis <- function(me) {
  message('Processing ', me, ' domains...')
  domain_ranges <- readRDS(str_c('rds/domains/', me, '_domains_ranges.rds'))
  domains_counts_file <- str_c('rds/domains/', me, '_domains_counts_for_CHiPSeq.rds')

  if (!file.exists(domains_counts_file)) {
    domain_ranges_as_SAF <- data.frame(
      GeneID = names(domain_ranges),
      Chr = as.character(seqnames(domain_ranges)),
      Start = start(ranges(domain_ranges)),
      End = end(ranges(domain_ranges)),
      Strand = as.character(strand(domain_ranges))
    )

    bam_folder <- '/data/ChipSeqHPVHNSCC/JHUJC01008_000_analysis/hg19_alignments/bam_files/'
    samples <- c('3', '4', '5', '6')
    bams <- c()
    inputs <- c()

    for (sample in samples) {
      dirbams <- dir(path = bam_folder, pattern = str_c(sample, '_h3k27ac*'), full.names = TRUE)
      bam <- (dirbams[grep('bam.bai', dirbams, invert = TRUE)])[1]
      bams <- c(bams, bam)

      dirinputs <- dir(path = bam_folder, pattern = str_c(sample, '_2percent*'), full.names = TRUE)
      input <- (dirinputs[grep('bam.bai', dirinputs, invert = TRUE)])[1]
      inputs <- c(inputs, input)
    }

    counts <- featureCounts(c(bams, inputs), annot.ext = domain_ranges_as_SAF)
    saveRDS(counts, domains_counts_file)
  } else {
    counts <- readRDS(domains_counts_file)
  }

  coverage <- counts$counts[, 1:8]

  sample.status <- factor(c('N', 'N', 'T', 'T', 'N', 'N', 'T', 'T'))
  track.type <- factor(c('CS', 'CS', 'CS', 'CS', 'IN', 'IN', 'IN', 'IN'))
  groups <- sample.status:track.type

  design <- model.matrix(~0 + groups)
  rownames(design) <- colnames(coverage)

  y <- DGEList(coverage, group = groups)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design = design, robust = TRUE)

  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, contrast = c(-1, 1, 1, -1))

  cpm_vals <- cpm(y, log = FALSE, normalized.lib.sizes = TRUE)

  saveRDS(y, str_c('rds/domains/', me, '_EDGER_y_object.rds'))
  saveRDS(cpm_vals, str_c('rds/domains/', me, '_EDGER_cpm.rds'))

  diff_domains <- topTags(qlf, adjust.method = 'fdr', p.value = 0.05, n = Inf)[[1]]
  all_domains <- topTags(qlf, adjust.method = 'fdr', p.value = 1, n = Inf)[[1]]

  coords_frame <- data.frame(domain_ranges)
  coords_frame$ID <- names(domain_ranges)

  all_domains$ID <- rownames(all_domains)
  all_domains <- inner_join(all_domains, coords_frame, by = 'ID')
  saveRDS(all_domains, str_c('rds/domains/', me, '_domains_EDGER_on_CHiPSeq.rds'))
  write.csv(all_domains, file = str_c('results/domains/', me, '_domains_EDGER_on_CHiPSeq.csv'), quote = FALSE, row.names = FALSE)

  diff_domains$ID <- rownames(diff_domains)
  diff_domains <- inner_join(diff_domains, coords_frame, by = 'ID')
  write.csv(diff_domains, file = str_c('results/domains/', me, '_diff_domains_EDGER_on_CHiPSeq.csv'), quote = FALSE, row.names = FALSE)

  cpm_frame <- data.frame(ID = rownames(cpm_vals), cpm_vals)
  colnames(cpm_frame) <- str_replace(colnames(cpm_frame), 'bam', 'bam.cpm')

  coverage_frame <- data.frame(ID = rownames(coverage), coverage)

  all_domains_complete <- inner_join(all_domains, coverage_frame, by = 'ID')
  all_domains_complete <- inner_join(all_domains_complete, cpm_frame, by = 'ID')

  write.csv(all_domains_complete, file = str_c('results/domains/', me, '_domains_EDGER_coverage_and_CPM.csv'), quote = FALSE, row.names = FALSE)
}

invisible(lapply(c('P', 'E', 'SE'), run_domain_analysis))
