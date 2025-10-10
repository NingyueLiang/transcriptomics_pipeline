
#####
# Principal Component Analysis Functions

# Required libraries
library(edgeR)      # For CPM normalization
library(ggplot2)    # For creating plots
library(tidyr)      # For data tidying
library(dplyr)      # For data manipulation
library(stats)      # For prcomp() and other statistical operations
library(scales)     # For axis formatting


perform_pca_analysis <- function(data, groups_to_include = NULL, min_cpm = 1, 
                                 n_top_features = 50) {
  # Input validation
  if (!("gene_symbol" %in% colnames(data))) {
    stop("Data must contain a 'gene_symbol' column")
  }
  
  # Parse sample information from column names
  sample_cols <- colnames(data)[!colnames(data) %in% "gene_symbol"]
  sample_info <- data.frame(
    sample_name = sample_cols,
    group = gsub("_[0-9]+$", "", sample_cols)  # Remove replicate numbers
  )
  
  # Filter groups if specified
  if (!is.null(groups_to_include)) {
    included_samples <- sample_info$sample_name[sample_info$group %in% groups_to_include]
    data <- data[, c("gene_symbol", included_samples)]
    sample_info <- sample_info[sample_info$sample_name %in% included_samples, ]
  }
  
  # Prepare counts matrix
  counts_matrix <- as.matrix(data[, -1])  # Remove gene_symbol column
  rownames(counts_matrix) <- data$gene_symbol
  
  # Calculate CPM normalization
  cpm_normalized <- cpm(counts_matrix)
  
  # Filter out genes with low counts or zero variance
  keep_genes <- rowSums(cpm_normalized >= min_cpm) >= 2
  gene_vars <- apply(cpm_normalized[keep_genes, ], 1, var)
  keep_genes_var <- gene_vars > 0
  
  # Apply filters
  cpm_filtered <- cpm_normalized[keep_genes, ][keep_genes_var, ]
  
  message(sprintf("Original number of genes: %d", nrow(cpm_normalized)))
  message(sprintf("Number of genes after filtering: %d", nrow(cpm_filtered)))
  
  if(nrow(cpm_filtered) == 0) {
    stop("No genes remaining after filtering. Consider adjusting the min_cpm threshold.")
  }
  
  # Perform PCA
  pca_result <- prcomp(t(cpm_filtered), 
                       center = TRUE,
                       scale. = TRUE)
  
  # Calculate variance explained
  var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
  
  # Create scree plot data
  scree_data <- data.frame(
    PC = 1:length(var_explained),
    Variance = var_explained,
    Cumulative = cumsum(var_explained)
  )
  
  # Create scree plot
  scree_plot <- ggplot(scree_data, aes(x = PC)) +
    geom_bar(aes(y = Variance), stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 3) +
    scale_y_continuous(
      name = "Variance Explained (%)",
      sec.axis = sec_axis(~., name = "Cumulative Variance (%)")
    ) +
    xlab("Principal Component") +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle("Scree Plot of Principal Components")
  
  # Extract top contributing features for each PC
  loadings <- pca_result$rotation
  top_features <- lapply(1:ncol(loadings), function(pc) {
    pc_loadings <- loadings[, pc]
    top_indices <- order(abs(pc_loadings), decreasing = TRUE)[1:n_top_features]
    data.frame(
      Feature = rownames(loadings)[top_indices],
      Loading = pc_loadings[top_indices],
      PC = paste0("PC", pc)
    )
  })
  top_features_df <- do.call(rbind, top_features)
  
  # Create PCA plot
  pc_scores <- as.data.frame(pca_result$x)
  plot_data <- cbind(pc_scores, 
                     treatment = sample_info$group,
                     sample = sample_info$sample_name)
  
  pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = treatment)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = sample), hjust = -0.2, vjust = -0.2, size = 3) +
    theme_bw() +
    xlab(sprintf("PC1 (%.1f%%)", var_explained[1])) +
    ylab(sprintf("PC2 (%.1f%%)", var_explained[2])) +
    ggtitle("Principal Components Analysis") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  # Return results as a list
  return(list(
    normalized_counts = cpm_filtered,
    pca_result = pca_result,
    pca_plot = pca_plot,
    scree_plot = scree_plot,
    variance_explained = var_explained,
    sample_info = sample_info,
    genes_retained = rownames(cpm_filtered),
    top_features = top_features_df,
    pc_scores = pc_scores
  ))
}

# Run analysis
results <- perform_pca_analysis(
  data = xen_tran_2024_06_assay_data,
  groups_to_include = c("Veh_SF", "WB3_SF")
)

# If needed, adjust the CPM threshold
# results_sf <- perform_pca_analysis(
#   data = xen_tran_2024_06_assay_data,
#   groups_to_include = c("Veh_SF", "WB3_SF"),
#   min_cpm = 0.5  # Lower threshold if too many genes are being filtered out
# )

# View scree plot to see variance explained
print(results$scree_plot)

# Look at top contributing features for PC1
head(subset(results$top_features, PC == "PC1"))

# Get PC scores for downstream modeling
pc_scores_for_modeling <- results$pc_scores

print(results_sf$pca_plot)


