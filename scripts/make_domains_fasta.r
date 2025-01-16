# Load necessary libraries
library(rtracklayer)
library(dplyr)
library(plyranges)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)

# Define input and output paths
domain <- "/data/SE_domains_ranges.rds"
outputDir <- "/data/fasta"

# Create the output folder if it doesn't already exist
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Process and write SE domains to FASTA
writeXStringSet(
  getSeq(
    BSgenome.Hsapiens.UCSC.hg19,  # Get sequences from hg19 genome
    readRDS(domain)  # Load SE ranges from RDS file
  ),
  filepath = file.path(outputDir, 'SE_domains.fasta'),  # Write to FASTA file
  format = 'fasta'
)