# Load necessary libraries
library(ggplot2)
library(ggridges)
library(viridis)
library(reshape2)

#setdir
setwd("/data") 

# Define p-value cutoffs for analysis
pval_cutoffs <- c(0.05)

# Loop through each p-value cutoff
for (pval_cutoff in pval_cutoffs) {
  # Define input and output directories
  inDir <- "/data"
  outDir <- "/data"
  
  # Read data
  df <- as.data.frame(readRDS(file.path(inDir, "fold_and_distance_77_samples_SED.rds")))
  
  # Set Tumor vs Normal as factor and rename levels
  df$SED_type <- ifelse(df$SED_type == "UPPP", "normal", df$SED_type)
  df$SED_type <- factor(df$SED_type, levels = c("tumor", "normal"))
  levels(df$SED_type) <- c("tumor", "normal")
  
  # changing "specific" to TRUE and "non_specific" to FALSE in df$SED_is_specific
  df$SED_is_specific <- ifelse(df$SED_is_specific == "specific", TRUE, FALSE)
  
  # Filter data based on specified conditions
  df.sub <- subset(df, SED_is_specific == TRUE & !is.na(log2FoldChange))
  df.sub$distance <- abs(df.sub$distance) / 1000000
  df.sub.filt <- subset(df.sub, distance < 2 & padj <= pval_cutoff)
  df.sub.filt$distbin <- cut(df.sub.filt$distance, seq(0, 2, 0.1), right = FALSE, labels = as.character(seq(0.1, 2, 0.1)))
  
  # Save filtered dataset
  out_filt <- paste0("fold_and_distance_77_samples_SED.filt.pval", pval_cutoff, ".rds")
  saveRDS(df.sub.filt, file.path(outDir, out_filt))
  
  # Plot log2FoldChange as ggridges plot
  pdf(paste0("lfc_dist_SED.pval", pval_cutoff, ".pdf"), height = 10, width = 10)
  p1 <- ggplot(df.sub.filt, aes(y = distbin, x = log2FoldChange, fill = SED_type)) +
    geom_density_ridges(alpha = 0.7, color = "white") +
    scale_fill_manual(values = c("normal" = "lightskyblue3", "tumor" = "coral3")) +  # Set colors for tumor and normal)
    labs(
      x = "mRNA log2 Fold Change",
      y = "Distance (100 kilobases)",
      title = "Differencial expression of SEDs by distance to the nearest gene (withing 2Mb)"
    ) +
    scale_y_discrete(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    theme_classic() # U
  
  print(p1)
  dev.off()
  
  # Perform statistical tests
  distbin <- c()
  wilc_pval <- c()
  
  for (each in levels(df.sub.filt$distbin)) {
    testdf <- subset(df.sub.filt, distbin == each)
    
    normal.dist <- subset(testdf, SED_type == "normal")$log2FoldChange
    tumor.dist <- subset(testdf, SED_type == "tumor")$log2FoldChange
    wilc <- wilcox.test(normal.dist, tumor.dist)
    distbin <- append(distbin, each)
    wilc_pval <- append(wilc_pval, wilc$p.value)
  }
  
  
  # Write statistical results to CSV
  outdf <- as.data.frame(cbind(distbin, wilc_pval))
  outdf <- cbind(outdf, FDR = p.adjust(outdf$wilc_pval, method = "BH"))
  write.csv(outdf, file = paste0(outDir, "/lfc_dist_SED.ks_stats.pval", pval_cutoff, ".csv"))
  cols <- viridis(7)
}


#################Genes with p<0.05
library(dplyr)

df <- readRDS("fold_and_distance_77_samples_SED.filt.pval0.05.rds")

result <- df %>%
  group_by(distbin, SED_type) %>%
  summarise(count = dplyr::n())

# Printing the result
print(result)

# Exporting the result as a CSV file
write.csv(result, "genes_per_bin_pval0.05.csv", row.names = FALSE)

# Your data
data <- data.frame(
  distbin = rep(seq(0.1, 2, by = 0.1), each = 2),
  SED_type = rep(c("tumor", "normal"), times = 20),
  count = result$count
)

# Creating the plot
plot <- ggplot(data, aes(y = as.factor(distbin), x = count, fill = SED_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Counts of 'Normal' and 'Tumor' within Each Distbin",
       x = "Count", y = "Distbin") +
  scale_fill_manual(values = c("normal" = "lightskyblue3", "tumor" = "coral3")) +  # Set consistent colors
  theme_minimal() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(title = "SED_type"))

# Saving the plot as a PDF file
filename <- paste0("counts_each_distbin_pval", 0.05, ".pdf")
ggsave(file.path(outDir, filename), plot, width = 10, height = 6, units = "in")
