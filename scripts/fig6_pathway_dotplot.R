##### Figure 6: Pathway Enrichment Dotplot
##### Gene set enrichment analysis of Hallmark pathways

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(cluster)
library(RColorBrewer)
library(viridis)
library(scales)

# Set working directory
setwd("~/data/")

# Function to read and process pathway enrichment files
read_pathway_files <- function(input_dir, file_pattern = "fgsea-mrna_.*\\.csv$") {
  # Get CSV files matching the pattern (updated to match mrna naming with underscore)
  files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
  
  # Debug information
  cat("Looking in directory:", input_dir, "\n")
  cat("Pattern:", file_pattern, "\n")
  all_files <- list.files(input_dir)
  cat("All files in directory:", paste(all_files, collapse = ", "), "\n")
  cat("Matching files:", paste(files, collapse = ", "), "\n")
  
  if (length(files) == 0) {
    stop("No CSV files found matching pattern 'fgsea-mrna_.*\\.csv$' in the specified directory")
  }
  
  # Required columns
  mandatory_columns <- c('pathway', 'NES', 'padj')
  
  # Initialize list to store data
  data_list <- list()
  
  # Read each file
  for (file_path in files) {
    # Extract sample name from filename (remove .csv extension)
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Read the data (CSV format with comma separators)
    df <- read_csv(file_path, show_col_types = FALSE)
    
    # Check for mandatory columns
    missing_cols <- setdiff(mandatory_columns, colnames(df))
    if (length(missing_cols) > 0) {
      warning(paste("Missing columns in", sample_name, ":", paste(missing_cols, collapse = ", ")))
      next
    }
    
    cat("Min padj for", sample_name, "is", min(df$padj, na.rm = TRUE), 
        ", max padj is", max(df$padj, na.rm = TRUE), "\n")
    
    # Add sample name column
    df$sample <- sample_name
    
    # Store in list
    data_list[[sample_name]] <- df
  }
  
  # Combine all data
  combined_data <- bind_rows(data_list)
  
  return(combined_data)
}

# Function to perform hierarchical clustering of pathways
cluster_pathways <- function(data, cluster_cols = c('NES'), max_clusters = 10) {
  # Create a wide format for clustering
  wide_data <- data %>%
    select(pathway, sample, all_of(cluster_cols)) %>%
    pivot_wider(names_from = sample, 
                values_from = all_of(cluster_cols),
                names_sep = "_")
  
  # Remove pathways with missing data
  complete_data <- wide_data[complete.cases(wide_data), ]
  
  # Prepare matrix for clustering (exclude pathway column)
  cluster_matrix <- as.matrix(complete_data[, -1])
  rownames(cluster_matrix) <- complete_data$pathway
  
  # Perform hierarchical clustering
  dist_matrix <- dist(cluster_matrix, method = "euclidean")
  hc <- hclust(dist_matrix, method = "ward.D2")
  
  # Cut tree to get clusters
  clusters <- cutree(hc, k = min(max_clusters, nrow(complete_data)))
  
  # Create pathway order based on clustering
  pathway_order <- complete_data$pathway[hc$order]
  
  # Add cluster information
  cluster_df <- data.frame(
    pathway = names(clusters),
    cluster = clusters,
    stringsAsFactors = FALSE
  )
  
  return(list(
    pathway_order = pathway_order,
    clusters = cluster_df,
    hclust = hc
  ))
}

# Function to sort pathways based on different criteria
sort_pathways <- function(data, pathway_sort = "cluster", pvalue_threshold = 0.05) {
  
  # Get sample names in reverse order (like Python version)
  sample_names <- sort(unique(data$sample), decreasing = TRUE)
  
  if (pathway_sort == "cluster") {
    # Use hierarchical clustering (only NES for clustering like updated Python)
    clustering_result <- cluster_pathways(data, cluster_cols = c('NES'))
    pathway_order <- clustering_result$pathway_order
    clusters <- clustering_result$clusters
    
  } else if (pathway_sort == "first_enrich") {
    # Sort by NES of first sample
    first_sample <- sample_names[1]
    sort_order_data <- data %>% 
      filter(sample == first_sample) %>%
      arrange(NES)
    pathway_order <- sort_order_data$pathway
    clusters <- NULL
    
  } else if (pathway_sort == "enrich_threshold") {
    # Sort by significance threshold then NES of first sample
    first_sample <- sample_names[1]
    sort_order_data <- data %>% 
      filter(sample == first_sample) %>%
      mutate(pval_thres = padj <= pvalue_threshold) %>%
      arrange(pval_thres, NES)
    pathway_order <- sort_order_data$pathway
    clusters <- NULL
    
  } else {
    stop("Unknown pathway_sort option. Use 'cluster', 'first_enrich', or 'enrich_threshold'")
  }
  
  return(list(
    pathway_order = pathway_order,
    clusters = clusters,
    sample_order = sample_names
  ))
}

# Function to filter pathways (only keep those with at least one significant result)
filter_pathways <- function(data, pvalue_threshold = 0.05) {
  # Find pathways where at least one sample has padj < threshold
  significant_pathways <- data %>%
    group_by(pathway) %>%
    summarise(has_significant = any(padj < pvalue_threshold, na.rm = TRUE), .groups = "drop") %>%
    filter(has_significant) %>%
    pull(pathway)
  
  # Filter data to only include these pathways
  filtered_data <- data %>%
    filter(pathway %in% significant_pathways)
  
  cat("Filtered", length(unique(data$pathway)), "pathways to", 
      length(significant_pathways), "pathways with significant results\n")
  
  return(filtered_data)
}

# Function to create pathway dotplot
create_pathway_dotplot <- function(data, 
                                   plot_type = "size", 
                                   pathway_sort = "cluster",
                                   pvalue_threshold = 0.05,
                                   title = "Gene set enrichment analysis of Hallmark pathways") {
  
  # Filter pathways to only include those with significant results
  filtered_data <- filter_pathways(data, pvalue_threshold)
  
  # Sort pathways and samples
  sorting_result <- sort_pathways(filtered_data, pathway_sort, pvalue_threshold)
  
  # Set sample order (reverse sorted like Python)
  filtered_data$sample <- factor(filtered_data$sample, levels = sorting_result$sample_order)
  
  # Set pathway order
  filtered_data$pathway <- factor(filtered_data$pathway, levels = sorting_result$pathway_order)
  
  # Add cluster information if clustering was used
  if (!is.null(sorting_result$clusters)) {
    filtered_data <- filtered_data %>%
      left_join(sorting_result$clusters, by = "pathway") %>%
      filter(!is.na(cluster))
  }
  
  # Calculate dot size based on plot type
  if (plot_type == "size") {
    # Use the same formula as Python: 50 * (-log10(padj)) / 1.6
    filtered_data$dot_size <- 50 * (-log10(filtered_data$padj)) / 1.6
  } else {
    filtered_data$dot_size <- 25 * 5  # Fixed size like Python
  }
  
  # Determine color based on significance for gray plot
  if (plot_type == "gray") {
    filtered_data$significant <- filtered_data$padj <= pvalue_threshold
  }
  
  # Create the base plot
  p <- ggplot(filtered_data, aes(x = sample, y = pathway))
  
  # Add dots based on plot type
  if (plot_type == "size") {
    p <- p + 
      geom_point(aes(size = dot_size, color = NES), alpha = 1) +
      scale_size_identity() +  # Use actual values for size
      scale_color_gradient2(name = "NES", 
                            low = "blue", mid = "white", high = "red",
                            midpoint = 0,
                            guide = guide_colorbar(order = 1)) +
      guides(size = "none")  # Don't show size legend
    
  } else if (plot_type == "gray") {
    p <- p + 
      geom_point(aes(color = ifelse(significant, NES, NA)), 
                 size = 3, alpha = 1) +
      geom_point(data = filter(filtered_data, !significant),
                 color = "gray", size = 3, alpha = 0.6) +
      scale_color_gradient2(name = "NES", 
                            low = "blue", mid = "white", high = "red",
                            midpoint = 0,
                            na.value = "gray")
    
  } else if (plot_type == "white") {
    # Only show significant pathways
    significant_data <- filter(filtered_data, padj <= pvalue_threshold)
    p <- ggplot(significant_data, aes(x = sample, y = pathway)) +
      geom_point(aes(color = NES), size = 3, alpha = 1) +
      scale_color_gradient2(name = "NES", 
                            low = "blue", mid = "white", high = "red",
                            midpoint = 0)
  }
  
  # Customize the plot
  p <- p +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray", linetype = "dashed", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.ontop = FALSE
    ) +
    labs(title = title)
  
  return(p)
}

# Main analysis function
create_figure6_dotplot <- function(input_dir, 
                                   output_file = "fig6_pathway_dotplot.pdf",
                                   plot_type = "size",
                                   pathway_sort = "cluster",
                                   pvalue_threshold = 0.05,
                                   width = 11, 
                                   height = 12) {
  
  cat("Reading pathway enrichment files from:", input_dir, "\n")
  
  # Read and combine data
  pathway_data <- read_pathway_files(input_dir)
  
  cat("Loaded data for", length(unique(pathway_data$sample)), "samples\n")
  cat("Found", length(unique(pathway_data$pathway)), "unique pathways\n")
  
  # Create the dotplot
  p <- create_pathway_dotplot(pathway_data, 
                              plot_type = plot_type,
                              pathway_sort = pathway_sort,
                              pvalue_threshold = pvalue_threshold)
  
  # Save the plot
  ggsave(output_file, plot = p, width = width, height = height, units = "in")
  
  cat("Figure 6 saved as:", output_file, "\n")
  
  # Print summary statistics
  cat("\n=== PATHWAY ANALYSIS SUMMARY ===\n")
  
  summary_stats <- pathway_data %>%
    group_by(sample) %>%
    summarise(
      total_pathways = dplyr::n(),
      significant_pathways = sum(padj <= pvalue_threshold, na.rm = TRUE),
      mean_NES = mean(NES, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(summary_stats)
  
  return(list(plot = p, data = pathway_data, summary = summary_stats))
}

# Example usage - uncomment and modify paths as needed:
# 
# # For size-based dotplot with clustering (default)
# result <- create_figure6_dotplot(
#   input_dir = "pathway_enrichment_results/",
#   output_file = "fig6A_pathway_dotplot_size.pdf",
#   plot_type = "size",
#   pathway_sort = "cluster"
# )
# 
# For gray-based dotplot with enrichment threshold sorting
result_gray <- create_figure6_dotplot(
  input_dir = "mrna-analysis/",
  output_file = "fig6B_pathway_dotplot_gray.pdf",
  plot_type = "gray",
  pathway_sort = "enrich_threshold",
  pvalue_threshold = 0.05
)
# 
# # For white-based dotplot with first sample enrichment sorting
# result_white <- create_figure6_dotplot(
#   input_dir = "pathway_enrichment_results/",
#   output_file = "fig6C_pathway_dotplot_white.pdf",
#   plot_type = "white",
#   pathway_sort = "first_enrich",
#   pvalue_threshold = 0.05
# )

cat("Figure 6 dotplot functions loaded successfully!\n")
cat("Use create_figure6_dotplot() function to generate your pathway dotplots.\n")
cat("\nRequired data format: CSV files matching 'fgsea-mrna*.csv' with columns 'pathway', 'NES', 'padj'\n")
cat("Available plot types: 'size', 'gray', 'white'\n")
cat("Available pathway sorts: 'cluster', 'first_enrich', 'enrich_threshold'\n")
