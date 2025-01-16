# Load necessary libraries
library(DESeq2)
library(apeglm)
library(data.table)

# Define file paths (Modify these paths as needed)
domain <- "SE"  # Domain identifier (e.g., "SE", "P", "E")
outputDir <- "/results"
samplesFile <- "/data/samples-HPV-RNA.txt"
tximportRDS <- "/data/SE_domains-txi.rds"

# Log file paths for confirmation
cat("Domain: ", domain, "\n")
cat("Output directory: ", outputDir, "\n")
cat("Samples file: ", samplesFile, "\n")
cat("Tximport RDS file: ", tximportRDS, "\n")

# Read the samples file
cat("Reading samples file...\n")
samples <- fread(samplesFile, header = FALSE)
colnames(samples) <- c("code")  # Rename column
rownames(samples) <- samples$code  # Set rownames to match sample codes

# Assign conditions based on sample code
cat("Assigning sample conditions...\n")
samples$condition <- ifelse(
  samples$code < 100, 'Tumor',
  ifelse(samples$code < 200, 'Normal',
         ifelse(samples$code < 300, 'Xenograph', 'Cell line'))
)

# Load tximport RDS file
cat("Loading tximport RDS file...\n")
txi <- readRDS(tximportRDS)

# Create output directory if it does not exist
saveDir <- file.path(outputDir, "diffexp")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}

# Construct DESeqDataSet from Tximport and sample information
cat("Constructing DESeqDataSet...\n")
ddsTxi <- DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ condition
)

# Perform differential expression analysis
cat("Running DESeq2 for differential expression analysis...\n")
dds <- DESeq(ddsTxi)

# Shrink log fold changes using apeglm
cat("Applying log fold change shrinkage using apeglm...\n")
resLFC <- lfcShrink(dds, coef = 2, type = "apeglm")

# Save results to RDS and CSV
deseqrdspath <- file.path(saveDir, paste0(domain, "-deseq-res.rds"))
cat("Saving DESeq results to RDS file: ", deseqrdspath, "\n")
saveRDS(resLFC, deseqrdspath)

deseqcsvpath <- file.path(saveDir, paste0(domain, "-deseq-res.csv"))
cat("Saving DESeq results to CSV file: ", deseqcsvpath, "\n")
write.csv(as.data.frame(resLFC), file = deseqcsvpath)

cat("Differential expression analysis complete!\n")
