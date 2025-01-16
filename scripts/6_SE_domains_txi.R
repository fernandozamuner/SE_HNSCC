
# Load necessary libraries
library(tximport)
library(readr)
library(data.table)

# Define file paths (Modify these paths as needed)
domain <- "SE"  # Domain identifier (e.g., "SE", "P", "E")
domainFile <- "/data/fasta/SE_domains.fasta"
outputDir <- "/results"
samplesFile <- "/data/samples-HPV-RNA.txt"
bedFile <- "/SE_domains.bed"

# Log file paths for confirmation
cat("Domain: ", domain, "\n")
cat("Domain file: ", domainFile, "\n")
cat("Output directory: ", outputDir, "\n")
cat("Samples file: ", samplesFile, "\n")
cat("BED file: ", bedFile, "\n")

# Read samples file
cat("Reading samples file...\n")
samples <- fread(samplesFile, header = FALSE)
colnames(samples) <- c("code")  # Rename column
rownames(samples) <- samples$code  # Set rownames

# Assign conditions based on sample code
cat("Assigning sample conditions...\n")
samples$condition <- ifelse(
  samples$code < 100, 'Tumor',
  ifelse(samples$code < 200, 'Normal',
         ifelse(samples$code < 300, 'Xenograph', 'Cell line'))
)

# Construct file paths for Salmon quantification results
cat("Constructing file paths for Salmon quantification...\n")
files <- file.path(outputDir, "salmon_result_files", samples$code, domain, "quant.sf.gz")
names(files) <- samples$code  # Assign sample names

# Check if all quant.sf.gz files exist
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  cat("Warning: The following files are missing:\n")
  print(missing_files)
}

# Read BED file to create tx2gene mapping
cat("Reading BED file...\n")
bed <- fread(bedFile, header = FALSE)
if (ncol(bed) < 4) {
  stop("Error: BED file must have at least 4 columns (chrom, start, end, name).")
}
tx2gene <- data.frame('Transcript_ID' = bed$V4, 'Gene_ID' = bed$V4)

# Perform tximport
cat("Running tximport...\n")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

# Create output directory if it does not exist
saveDir <- file.path(outputDir, 'diffexp')
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}

# Save tximport results as an RDS file
output_rds <- file.path(saveDir, paste0(domain, "-txi.rds"))
cat("Saving tximport results to:", output_rds, "\n")
saveRDS(txi, output_rds)

cat("Process completed successfully!\n")

