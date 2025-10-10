#' Principal Component Analysis
#' 
#' PCA analysis and visualization for transcriptomics data.

# Required libraries
suppressPackageStartupMessages({
    library(edgeR)      # For CPM normalization
    library(ggplot2)    # For creating plots
    library(tidyr)      # For data tidying
    library(dplyr)      # For data manipulation
    library(stats)      # For prcomp() and other statistical operations
    library(scales)     # For axis formatting
    library(here)       # For path management
})

# Optional logger package
if (requireNamespace("logger", quietly = TRUE)) {
    library(logger)
} else {
    # Fallback logging function
    log_info <- function(...) message(...)
    log_error <- function(...) stop(...)
    log_warn <- function(...) warning(...)
}

# Configuration
pca_config <- list(
    default_min_cpm = 1,
    default_min_samples = 2,
    default_n_top_features = 50,
    default_n_pcs_to_plot = 10
)

#' Validate PCA parameters
#' @param data Input data frame
#' @param groups_to_include Groups to include in analysis
#' @param min_cpm Minimum CPM threshold
#' @param n_top_features Number of top features to extract
validate_pca_params <- function(data, groups_to_include, min_cpm, n_top_features) {
    if (!is.data.frame(data) || nrow(data) == 0) {
        stop("data must be a non-empty data frame")
    }
    
    if (!"gene_symbol" %in% colnames(data)) {
        stop("Data must contain a 'gene_symbol' column")
    }
    
    if (!is.numeric(min_cpm) || min_cpm < 0) {
        stop("min_cpm must be a non-negative number")
    }
    
    if (!is.numeric(n_top_features) || n_top_features < 1) {
        stop("n_top_features must be a positive integer")
    }
    
    if (!is.null(groups_to_include) && !is.character(groups_to_include)) {
        stop("groups_to_include must be a character vector or NULL")
    }
}

#' Parse sample information from column names
#' @param data Input data frame
#' @param groups_to_include Groups to include
#' @param sample_metadata Optional sample metadata data frame
#' @return List containing filtered data and sample info
parse_sample_info <- function(data, groups_to_include = NULL, sample_metadata = NULL) {
    # Get sample columns
    sample_cols <- colnames(data)[!colnames(data) %in% "gene_symbol"]
    
    if (length(sample_cols) == 0) {
        stop("No sample columns found in data")
    }
    
    # Use provided sample metadata if available
    if (!is.null(sample_metadata)) {
        # Match sample columns to metadata
        sample_info <- sample_metadata[sample_metadata$sample_id %in% sample_cols, ]
        
        if (nrow(sample_info) == 0) {
            stop("No matching samples found in metadata")
        }
        
        # Use condition as group if available
        if ("condition" %in% colnames(sample_info)) {
            sample_info$group <- sample_info$condition
        } else {
            # Fallback to parsing from sample names
            sample_info$group <- gsub("_[0-9]+$", "", sample_info$sample_id)
        }
        
        # Rename sample_id to sample_name for consistency
        sample_info$sample_name <- sample_info$sample_id
    } else {
        # Create sample info data frame from column names
        sample_info <- data.frame(
            sample_name = sample_cols,
            group = gsub("_[0-9]+$", "", sample_cols),  # Remove replicate numbers
            stringsAsFactors = FALSE
        )
    }
    
    # Filter groups if specified
    if (!is.null(groups_to_include)) {
        included_samples <- sample_info$sample_name[sample_info$group %in% groups_to_include]
        
        if (length(included_samples) == 0) {
            stop("No samples found for specified groups: ", paste(groups_to_include, collapse = ", "))
        }
        
        data <- data[, c("gene_symbol", included_samples)]
        sample_info <- sample_info[sample_info$sample_name %in% included_samples, ]
        
        message("Filtered to ", length(included_samples), " samples from groups: ", 
                paste(groups_to_include, collapse = ", "))
    }
    
    return(list(data = data, sample_info = sample_info))
}

#' Prepare and normalize counts matrix
#' @param data Input data frame
#' @param min_cpm Minimum CPM threshold
#' @param min_samples Minimum number of samples
#' @return Normalized and filtered counts matrix
prepare_counts_matrix <- function(data, min_cpm = pca_config$default_min_cpm, 
                                 min_samples = pca_config$default_min_samples) {
    # Prepare counts matrix
    counts_matrix <- as.matrix(data[, -1])  # Remove gene_symbol column
    rownames(counts_matrix) <- data$gene_symbol
    
    # Validate counts matrix
    if (any(is.na(counts_matrix)) || any(is.infinite(counts_matrix))) {
        warning("Counts matrix contains NA or infinite values")
    }
    
    # Calculate CPM normalization
    message("Calculating CPM normalization...")
    cpm_normalized <- cpm(counts_matrix)
    
    # Filter out genes with low counts
    keep_genes <- rowSums(cpm_normalized >= min_cpm) >= min_samples
    message("Genes passing CPM filter: ", sum(keep_genes), " out of ", nrow(cpm_normalized))
    
    # Filter out genes with zero variance
    gene_vars <- apply(cpm_normalized[keep_genes, ], 1, var)
    keep_genes_var <- gene_vars > 0
    message("Genes with non-zero variance: ", sum(keep_genes_var), " out of ", sum(keep_genes))
    
    # Apply filters
    cpm_filtered <- cpm_normalized[keep_genes, ][keep_genes_var, ]
    
    if (nrow(cpm_filtered) == 0) {
        stop("No genes remaining after filtering. Consider adjusting the min_cpm threshold.")
    }
    
    message("Final filtered matrix: ", nrow(cpm_filtered), " genes x ", ncol(cpm_filtered), " samples")
    
    return(cpm_filtered)
}

#' Perform PCA with optimization
#' @param counts_matrix Normalized counts matrix
#' @param center Whether to center the data
#' @param scale Whether to scale the data
#' @return PCA result object
perform_pca <- function(counts_matrix, center = TRUE, scale = TRUE) {
    message("Performing PCA...")
    
    # Validate input
    if (nrow(counts_matrix) < 2 || ncol(counts_matrix) < 2) {
        stop("Insufficient data for PCA (need at least 2 genes and 2 samples)")
    }
    
    # Check for constant genes (should be filtered out already)
    gene_vars <- apply(counts_matrix, 1, var)
    if (any(gene_vars == 0)) {
        warning("Some genes have zero variance and will be excluded from PCA")
        counts_matrix <- counts_matrix[gene_vars > 0, ]
    }
    
    # Perform PCA
    pca_result <- prcomp(t(counts_matrix), center = center, scale. = scale)
    
    # Calculate variance explained
    var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
    
    message("PCA completed. First 5 PCs explain ", 
            round(sum(var_explained[1:min(5, length(var_explained))]), 1), "% of variance")
    
    return(list(
        pca_result = pca_result,
        variance_explained = var_explained
    ))
}

#' Create scree plot
#' @param variance_explained Vector of variance explained by each PC
#' @param n_pcs Number of PCs to plot
#' @return ggplot object
create_scree_plot <- function(variance_explained, n_pcs = pca_config$default_n_pcs_to_plot) {
    n_pcs <- min(n_pcs, length(variance_explained))
    
    # Create scree plot data
    scree_data <- data.frame(
        PC = 1:n_pcs,
        Variance = variance_explained[1:n_pcs],
        Cumulative = cumsum(variance_explained[1:n_pcs])
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
    
    return(scree_plot)
}

#' Extract top contributing features for each PC
#' @param pca_result PCA result object
#' @param n_top_features Number of top features to extract
#' @return Data frame with top features
extract_top_features <- function(pca_result, n_top_features = pca_config$default_n_top_features) {
    loadings <- pca_result$rotation
    n_pcs <- min(ncol(loadings), 10)  # Limit to first 10 PCs
    
    top_features <- lapply(1:n_pcs, function(pc) {
        pc_loadings <- loadings[, pc]
        top_indices <- order(abs(pc_loadings), decreasing = TRUE)[1:n_top_features]
        data.frame(
            Feature = rownames(loadings)[top_indices],
            Loading = pc_loadings[top_indices],
            PC = paste0("PC", pc),
            stringsAsFactors = FALSE
        )
    })
    
    top_features_df <- do.call(rbind, top_features)
    return(top_features_df)
}

#' Create PCA plot
#' @param pca_result PCA result object
#' @param sample_info Sample information data frame
#' @param variance_explained Vector of variance explained
#' @param pc_x Principal component for x-axis (default: 1)
#' @param pc_y Principal component for y-axis (default: 2)
#' @return ggplot object
create_pca_plot <- function(pca_result, sample_info, variance_explained, 
                           pc_x = 1, pc_y = 2) {
    # Extract PC scores
    pc_scores <- as.data.frame(pca_result$x)
    
    # Validate PC indices
    if (pc_x > ncol(pc_scores) || pc_y > ncol(pc_scores)) {
        stop("Requested PC indices exceed available components")
    }
    
    # Create plot data
    plot_data <- cbind(
        pc_scores[, c(pc_x, pc_y)],
        treatment = sample_info$group,
        sample = sample_info$sample_name
    )
    
    # Create plot
    pca_plot <- ggplot(plot_data, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y), color = "treatment")) +
        geom_point(size = 3, alpha = 0.8) +
        geom_text(aes(label = sample), hjust = -0.2, vjust = -0.2, size = 3) +
        theme_bw() +
        xlab(sprintf("PC%d (%.1f%%)", pc_x, variance_explained[pc_x])) +
        ylab(sprintf("PC%d (%.1f%%)", pc_y, variance_explained[pc_y])) +
        ggtitle("Principal Components Analysis") +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "right",
            legend.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )
    
    return(pca_plot)
}

#' Save PCA results
#' @param results PCA results list
#' @param output_dir Output directory
#' @param experiment_name Experiment name
save_pca_results <- function(results, output_dir, experiment_name) {
    if (!dir.exists(output_dir)) {
        if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
            stop("Failed to create output directory: ", output_dir)
        }
    }
    
    tryCatch({
        # Save plots
        plot_file <- file.path(output_dir, paste0(experiment_name, "_pca_plot.png"))
        ggsave(plot_file, results$pca_plot, width = 10, height = 8, dpi = 300)
        message("PCA plot saved to: ", plot_file)
        
        scree_file <- file.path(output_dir, paste0(experiment_name, "_scree_plot.png"))
        ggsave(scree_file, results$scree_plot, width = 10, height = 6, dpi = 300)
        message("Scree plot saved to: ", scree_file)
        
        # Save results
        results_file <- file.path(output_dir, paste0(experiment_name, "_pca_results.rds"))
        saveRDS(results, results_file)
        message("PCA results saved to: ", results_file)
        
        # Save top features
        features_file <- file.path(output_dir, paste0(experiment_name, "_top_features.csv"))
        write.csv(results$top_features, features_file, row.names = FALSE)
        message("Top features saved to: ", features_file)
        
    }, error = function(e) {
        stop("Failed to save PCA results: ", e$message)
    })
}

#' Main PCA analysis function
#' @param data Input data frame with gene expression data
#' @param groups_to_include Groups to include in analysis (optional)
#' @param min_cpm Minimum CPM threshold for gene filtering
#' @param n_top_features Number of top features to extract for each PC
#' @param experiment_name Name of experiment (for saving results)
#' @param save_results Whether to save results to files
#' @param output_dir Output directory for saving results
#' @return List containing PCA results
perform_pca_analysis <- function(data, 
                                groups_to_include = NULL, 
                                min_cpm = pca_config$default_min_cpm, 
                                n_top_features = pca_config$default_n_top_features,
                                experiment_name = NULL,
                                save_results = FALSE,
                                output_dir = NULL,
                                sample_metadata = NULL) {
    
    # Validate parameters
    validate_pca_params(data, groups_to_include, min_cpm, n_top_features)
    
    message("Starting PCA analysis...")
    if (!is.null(experiment_name)) {
        message("Experiment: ", experiment_name)
    }
    
    # Parse sample information
    parsed_data <- parse_sample_info(data, groups_to_include, sample_metadata)
    filtered_data <- parsed_data$data
    sample_info <- parsed_data$sample_info
    
    # Prepare and normalize counts matrix
    counts_matrix <- prepare_counts_matrix(filtered_data, min_cpm)
    
    # Perform PCA
    pca_output <- perform_pca(counts_matrix)
    pca_result <- pca_output$pca_result
    variance_explained <- pca_output$variance_explained
    
    # Create plots
    scree_plot <- create_scree_plot(variance_explained)
    pca_plot <- create_pca_plot(pca_result, sample_info, variance_explained)
    
    # Extract top features
    top_features <- extract_top_features(pca_result, n_top_features)
    
    # Create results list
    results <- list(
        normalized_counts = counts_matrix,
        pca_result = pca_result,
        pca_plot = pca_plot,
        scree_plot = scree_plot,
        variance_explained = variance_explained,
        sample_info = sample_info,
        genes_retained = rownames(counts_matrix),
        top_features = top_features,
        pc_scores = as.data.frame(pca_result$x)
    )
    
    # Save results if requested
    if (save_results && !is.null(experiment_name)) {
        if (is.null(output_dir)) {
            output_dir <- here("results", experiment_name, "pca")
        }
        save_pca_results(results, output_dir, experiment_name)
    }
    
    message("PCA analysis completed successfully")
    return(results)
}

# Additional functions from original implementation

#' Parse and filter sample information
#' @param data Input data frame
#' @param groups_to_include Groups to include in analysis
#' @return List containing filtered data and sample info
parse_and_filter_samples <- function(data, groups_to_include = NULL) {
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
  
  return(list(data = data, sample_info = sample_info))
}

#' Prepare and normalize counts matrix with filtering
#' @param data Input data frame
#' @param min_cpm Minimum CPM threshold
#' @return Normalized and filtered counts matrix
prepare_and_normalize_counts <- function(data, min_cpm = 1) {
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
  
  return(cpm_filtered)
}

#' Perform PCA and calculate variance explained
#' @param cpm_filtered Normalized and filtered counts matrix
#' @return List containing PCA result and variance explained
perform_pca_calculation <- function(cpm_filtered) {
  # Perform PCA
  pca_result <- prcomp(t(cpm_filtered), 
                       center = TRUE,
                       scale. = TRUE)
  
  # Calculate variance explained
  var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
  
  return(list(
    pca_result = pca_result,
    variance_explained = var_explained
  ))
}

#' Create comprehensive scree plot
#' @param variance_explained Vector of variance explained by each PC
#' @return ggplot object
create_comprehensive_scree_plot_internal <- function(variance_explained) {
  # Create scree plot data
  scree_data <- data.frame(
    PC = 1:length(variance_explained),
    Variance = variance_explained,
    Cumulative = cumsum(variance_explained)
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
  
  return(scree_plot)
}

#' Extract top contributing features for each PC
#' @param pca_result PCA result object
#' @param n_top_features Number of top features to extract
#' @return Data frame with top features
extract_comprehensive_top_features_internal <- function(pca_result, n_top_features = 50) {
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
  
  return(top_features_df)
}

#' Create comprehensive PCA plot
#' @param pca_result PCA result object
#' @param sample_info Sample information data frame
#' @param variance_explained Vector of variance explained
#' @return ggplot object
create_comprehensive_pca_plot_internal <- function(pca_result, sample_info, variance_explained) {
  # Create PCA plot
  pc_scores <- as.data.frame(pca_result$x)
  plot_data <- cbind(pc_scores, 
                     treatment = sample_info$group,
                     sample = sample_info$sample_name)
  
  pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = treatment)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = sample), hjust = -0.2, vjust = -0.2, size = 3) +
    theme_bw() +
    xlab(sprintf("PC1 (%.1f%%)", variance_explained[1])) +
    ylab(sprintf("PC2 (%.1f%%)", variance_explained[2])) +
    ggtitle("Principal Components Analysis") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  return(pca_plot)
}

#' Perform comprehensive PCA analysis with full workflow
#' @param data Input data frame with gene expression data
#' @param groups_to_include Groups to include in analysis (optional)
#' @param min_cpm Minimum CPM threshold for gene filtering
#' @param n_top_features Number of top features to extract for each PC
#' @return List containing comprehensive PCA results
perform_pca_analysis_comprehensive <- function(data, groups_to_include = NULL, min_cpm = 1, 
                                 n_top_features = 50) {
  # Input validation
  if (!("gene_symbol" %in% colnames(data))) {
    stop("Data must contain a 'gene_symbol' column")
  }
  
  # 1. Parse and filter sample information
  sample_data <- parse_and_filter_samples(data, groups_to_include)
  
  # 2. Prepare and normalize counts matrix
  cpm_filtered <- prepare_and_normalize_counts(sample_data$data, min_cpm)
  
  # 3. Perform PCA and calculate variance explained
  pca_output <- perform_pca_calculation(cpm_filtered)
  
  # 4. Create scree plot
  scree_plot <- create_comprehensive_scree_plot_internal(pca_output$variance_explained)
  
  # 5. Extract top features
  top_features_df <- extract_comprehensive_top_features_internal(pca_output$pca_result, n_top_features)
  
  # 6. Create PCA plot
  pca_plot <- create_comprehensive_pca_plot_internal(pca_output$pca_result, sample_data$sample_info, pca_output$variance_explained)
  
  # Return results as a list
  return(list(
    normalized_counts = cpm_filtered,
    pca_result = pca_output$pca_result,
    pca_plot = pca_plot,
    scree_plot = scree_plot,
    variance_explained = pca_output$variance_explained,
    sample_info = sample_data$sample_info,
    genes_retained = rownames(cpm_filtered),
    top_features = top_features_df,
    pc_scores = as.data.frame(pca_output$pca_result$x)
  ))
}

#' Create comprehensive scree plot with detailed formatting
#' @param variance_explained Vector of variance explained by each PC
#' @param n_pcs Number of PCs to plot
#' @return ggplot object
create_comprehensive_scree_plot <- function(variance_explained, n_pcs = 10) {
  n_pcs <- min(n_pcs, length(variance_explained))
  
  # Create scree plot data
  scree_data <- data.frame(
    PC = 1:n_pcs,
    Variance = variance_explained[1:n_pcs],
    Cumulative = cumsum(variance_explained[1:n_pcs])
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
  
  return(scree_plot)
}

#' Create comprehensive PCA plot with detailed formatting
#' @param pca_result PCA result object
#' @param sample_info Sample information data frame
#' @param variance_explained Vector of variance explained
#' @param pc_x Principal component for x-axis (default: 1)
#' @param pc_y Principal component for y-axis (default: 2)
#' @return ggplot object
create_comprehensive_pca_plot <- function(pca_result, sample_info, variance_explained, 
                           pc_x = 1, pc_y = 2) {
  # Extract PC scores
  pc_scores <- as.data.frame(pca_result$x)
  
  # Validate PC indices
  if (pc_x > ncol(pc_scores) || pc_y > ncol(pc_scores)) {
    stop("Requested PC indices exceed available components")
  }
  
  # Create plot data
  plot_data <- cbind(
    pc_scores[, c(pc_x, pc_y)],
    treatment = sample_info$group,
    sample = sample_info$sample_name
  )
  
  # Create plot
  pca_plot <- ggplot(plot_data, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y), color = "treatment")) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = sample), hjust = -0.2, vjust = -0.2, size = 3) +
    theme_bw() +
    xlab(sprintf("PC%d (%.1f%%)", pc_x, variance_explained[pc_x])) +
    ylab(sprintf("PC%d (%.1f%%)", pc_y, variance_explained[pc_y])) +
    ggtitle("Principal Components Analysis") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  return(pca_plot)
}

#' Extract comprehensive top contributing features for each PC
#' @param pca_result PCA result object
#' @param n_top_features Number of top features to extract
#' @return Data frame with top features
extract_comprehensive_top_features <- function(pca_result, n_top_features = 50) {
  loadings <- pca_result$rotation
  n_pcs <- min(ncol(loadings), 10)  # Limit to first 10 PCs
  
  top_features <- lapply(1:n_pcs, function(pc) {
    pc_loadings <- loadings[, pc]
    top_indices <- order(abs(pc_loadings), decreasing = TRUE)[1:n_top_features]
    data.frame(
      Feature = rownames(loadings)[top_indices],
      Loading = pc_loadings[top_indices],
      PC = paste0("PC", pc),
      stringsAsFactors = FALSE
    )
  })
  
  top_features_df <- do.call(rbind, top_features)
  return(top_features_df)
}

#' Save comprehensive PCA results
#' @param results PCA results list
#' @param output_dir Output directory
#' @param experiment_name Experiment name
save_comprehensive_pca_results <- function(results, output_dir, experiment_name) {
  if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Failed to create output directory: ", output_dir)
  }
  
  tryCatch({
    # Save plots
    plot_file <- file.path(output_dir, paste0(experiment_name, "_pca_plot.png"))
    ggsave(plot_file, results$pca_plot, width = 10, height = 8, dpi = 300)
    message("PCA plot saved to: ", plot_file)
    
    scree_file <- file.path(output_dir, paste0(experiment_name, "_scree_plot.png"))
    ggsave(scree_file, results$scree_plot, width = 10, height = 6, dpi = 300)
    message("Scree plot saved to: ", scree_file)
    
    # Save results
    results_file <- file.path(output_dir, paste0(experiment_name, "_pca_results.rds"))
    saveRDS(results, results_file)
    message("PCA results saved to: ", results_file)
    
    # Save top features
    features_file <- file.path(output_dir, paste0(experiment_name, "_top_features.csv"))
    write.csv(results$top_features, features_file, row.names = FALSE)
    message("Top features saved to: ", features_file)
    
  }, error = function(e) {
    stop("Failed to save PCA results: ", e$message)
  })
}

message("PCA functions loaded successfully. Use perform_pca_analysis() or perform_pca_analysis_comprehensive() to perform principal component analysis.")
