#' Generate Scatterplots for Cell Lines
#' This script reads data for two cell lines (JQ1_047 and JQ1_090), filters genes with negative log2FoldChange, and generates scatterplots comparing these cell lines with patient data.

########## Fernando ######## Feb 2024 ################  

# Install and Load Libraries 
library(scales) 
library(ggpubr) 
library(dplyr) 
library(tidyverse) 

setwd("/data") 

# Read data for JQ1_047 and JQ1_090 
fold_and_distance_JQ1_047 <- readRDS("fold_and_distance_JQ1_047.rds") 
fold_and_distance_JQ1_090 <- readRDS("fold_and_distance_JQ1_090.rds") 

# Remove NA values from relevant columns 
fold_and_distance_JQ1_047 <- na.omit(fold_and_distance_JQ1_047) 
fold_and_distance_JQ1_090 <- na.omit(fold_and_distance_JQ1_090) 

# Calculate counts of log2FoldChange negative for both datasets 
log2FC_negative_047 <- sum(   
  fold_and_distance_JQ1_047$CLSE_is_confirmed &     
    fold_and_distance_JQ1_047$padj < 0.05 &     
    fold_and_distance_JQ1_047$lfcSE * 3 < abs(fold_and_distance_JQ1_047$log2FoldChange) &     
    fold_and_distance_JQ1_047$log2FoldChange < 0 
)  

log2FC_negative_090 <- sum(   
  fold_and_distance_JQ1_090$CLSE_is_confirmed &     
    fold_and_distance_JQ1_090$padj < 0.05 &     
    fold_and_distance_JQ1_090$lfcSE * 3 < abs(fold_and_distance_JQ1_090$log2FoldChange) &     
    fold_and_distance_JQ1_090$log2FoldChange < 0 
)  

# Filter genes with log2FoldChange negative for both datasets
down_047 <- fold_and_distance_JQ1_047 %>%   
  filter(     
    CLSE_is_confirmed,     
    padj < 0.05,     
    lfcSE * 3 < abs(log2FoldChange),     
    log2FoldChange < 0   
  ) %>%   
  mutate(gene_id = row.names(.))  # Create a new column for gene_id  

down_090 <- fold_and_distance_JQ1_090 %>%   
  filter(     
    CLSE_is_confirmed,     
    padj < 0.05,     
    lfcSE * 3 < abs(log2FoldChange),     
    log2FoldChange < 0   
  ) %>%   
  mutate(gene_id = row.names(.))  # Create a new column for gene_id  

# Save the subsetted overlapping genes to a new RDS file 
saveRDS(down_047, file = "down_047.rds") 
saveRDS(down_090, file = "down_090.rds")  

generate_scatterplot <- function(cell_lines = c("090", "047")) {   
  # Load data for patients   
  genes_patients <- readRDS('fold_and_distance_77_samples_SED.rds')      
  
  # Load cell line data for each specified cell line   
  genes_cell_lines <- lapply(cell_lines, function(cell_line) {     
    file_name <- paste0('down_', cell_line, '.rds')
    readRDS(file_name)   
  })      
  
  # Add a new column for gene_id in patients data   
  genes_patients <- genes_patients %>%     
    mutate(gene_id = row.names(.))      
  
  # Add a new column for gene_id in each cell line data   
  genes_cell_lines <- lapply(genes_cell_lines, function(genes_cell_line) {     
    genes_cell_line %>%       
      mutate(gene_id = row.names(.))   
  })      
  
  # Filter genes_patients based on conditions   
  genes_patients_filtered <- genes_patients %>%     
    filter(       
      SED_type == "tumor",       
      padj < 0.05,       
      log2FoldChange > 0     
    )      
  
  # Subset genes from fold_and_distance_77_samples_SED using genes_cell_lines   
  genes_patients_filtered <- genes_patients_filtered %>%     
    filter(gene_id %in% unlist(lapply(genes_cell_lines, function(x) x$gene_id)))      
  
  # Generate scatterplots for each cell line   
  for (i in seq_along(cell_lines)) {     
    cell_line <- cell_lines[i]     
    cell_line_data <- genes_cell_lines[[i]]          
    
    # Merge the data frames based on gene_id     
    merged_data <- merge(cell_line_data, genes_patients_filtered, by = "gene_id", suffixes = c(paste0("_", cell_line), "_patients"))          
    
    # Scatterplot with p-values     
    title <- paste("Scatterplot of log2FoldChange - Cell Line", cell_line)          
    
    scatterplot <- ggplot(merged_data, aes(x = !!sym(paste0("log2FoldChange_", cell_line)), y = log2FoldChange_patients)) +       
      geom_point() +       
      labs(title = title,            
           x = paste("log2FoldChange - Cell Line", cell_line),            
           y = "log2FoldChange - Patients") +       
      theme_minimal() +       
      geom_smooth(method = "lm", se = FALSE) +       
      annotate("text", x = Inf, y = -Inf, label = paste("n =", nrow(merged_data)), hjust = 1, vjust = 0) +       
      stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, size = 4, vjust = 1, hjust = 0)          
    
    # Display the scatterplot     
    print(scatterplot)          
    
    # Define file name for saving scatterplot     
    file_name <- paste0('scatterplot_', cell_line, '.pdf')          
    
    # Save the scatterplot as a PDF file     
    ggsave(file_name, scatterplot, width = 8, height = 6)   
  } 
}  

#fig3_DE
# Generate scatterplots for both cell lines 090 and 047 
generate_scatterplot()
