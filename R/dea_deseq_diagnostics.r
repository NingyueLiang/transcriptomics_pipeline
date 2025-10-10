#' Diagnostic Functions
#' 
#' Generate quality control plots and reports for DESeq2 analysis.

# Required libraries
suppressPackageStartupMessages({
    library(DESeq2)         # For DESeq2 analysis
    library(ggplot2)        # For creating plots
    library(dplyr)          # For data manipulation
    library(tidyr)          # For data tidying
    library(here)           # For path management
    library(flextable)      # For creating formatted tables
    library(officer)        # For saving Word documents
    library(vsn)            # For meanSdPlot
    library(pheatmap)       # For heatmaps
    library(hexbin)         # For hexbin plots
    library(S4Vectors)      # For mcols function
    library(BiocGenerics)   # For basic Bioconductor functions
    library(SummarizedExperiment) # For assays function
    library(GenomicRanges)  # For genomic ranges functionality
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
diagnostics_config <- list(
    default_plot_params = list(
        width = 10,
        height = 8,
        dpi = 300
    ),
    valid_analysis_types = c("simple", "complex"),
    default_plot_types = c("dispersion", "ma", "pvalue_histogram", "cooks_distance", "pca", "mean_variance", "size_factors", "sample_distances")
)

#' Validate diagnostics parameters
#' @param experiment_name Name of experiment
#' @param analysis_type Type of analysis
#' @param plot_types Types of plots to generate
validate_diagnostics_params <- function(experiment_name, analysis_type, plot_types) {
    if (is.null(experiment_name) || nchar(experiment_name) == 0) {
        stop("experiment_name cannot be empty")
    }
    
    if (!analysis_type %in% diagnostics_config$valid_analysis_types) {
        stop("Invalid analysis_type. Must be one of: ", paste(diagnostics_config$valid_analysis_types, collapse = ", "))
    }
    
    if (!all(plot_types %in% diagnostics_config$default_plot_types)) {
        stop("Invalid plot_types. Must be one of: ", paste(diagnostics_config$default_plot_types, collapse = ", "))
    }
}

#' Load experiment data with validation
#' @param experiment_name Name of experiment
#' @return Experiment object
load_experiment_data <- function(experiment_name) {
    experiment_path <- here("data", experiment_name, paste0(experiment_name, ".rds"))
    
    if (!file.exists(experiment_path)) {
        stop("Experiment data not found at: ", experiment_path)
    }
    
    tryCatch({
        experiment_obj <- readRDS(experiment_path)
        
        if (!is.list(experiment_obj) || !"metadata" %in% names(experiment_obj)) {
            stop("Invalid experiment data structure")
        }
        
        return(experiment_obj)
    }, error = function(e) {
        stop("Failed to load experiment data: ", e$message)
    })
}

#' Load DESeq2 results with validation
#' @param experiment_name Name of experiment
#' @param analysis_type Type of analysis
#' @return DESeq2 results
load_deseq2_results <- function(experiment_name, analysis_type) {
    results_filename <- if (analysis_type == "complex") {
        paste0(experiment_name, "_de_complex_results.rds")
    } else {
        paste0(experiment_name, "_de_results.rds")
    }
    
    results_path <- here("results", experiment_name, "de", "deseq2", analysis_type, results_filename)
    
    if (!file.exists(results_path)) {
        stop("DESeq2 results not found at: ", results_path)
    }
    
    tryCatch({
        dea_results <- readRDS(results_path)
        
        if (!is.list(dea_results)) {
            stop("Invalid DESeq2 results structure")
        }
        
        return(dea_results)
    }, error = function(e) {
        stop("Failed to load DESeq2 results: ", e$message)
    })
}

#' Create output directories for diagnostics
#' @param experiment_name Name of experiment
#' @param analysis_type Type of analysis
#' @return List of directory paths
create_diagnostics_output_directories <- function(experiment_name, analysis_type) {
    base_dir <- here("results", experiment_name, "de", "deseq2", analysis_type)
    plots_dir <- file.path(base_dir, "plots", "diagnostics")
    
    if (!dir.exists(plots_dir)) {
        if (!dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)) {
            stop("Failed to create plots directory: ", plots_dir)
        }
    }
    
    return(list(plots_dir = plots_dir))
}

#' Create dispersion plot data
#' @param dds DESeq2 dataset
#' @return Dispersion plot data
create_dispersion_data <- function(dds) {
    # Get dispersion estimates
    dispersions <- mcols(dds)$dispGeneEst
    dispersions_fit <- mcols(dds)$dispFit
    dispersions_final <- dispersions(dds)
    
    # Create data frame
    plot_data <- data.frame(
        mean = mcols(dds)$baseMean,
        gene_est = dispersions,
        fit = dispersions_fit,
        final = dispersions_final,
        stringsAsFactors = FALSE
    )
    
    # Remove rows with missing values
    plot_data <- plot_data[complete.cases(plot_data), ]
    
    return(plot_data)
}

#' Create dispersion plot
#' @param plot_data Dispersion plot data
#' @param title Plot title
#' @return ggplot object
create_dispersion_plot <- function(plot_data, title) {
    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$mean)) +
        ggplot2::geom_point(ggplot2::aes(y = .data$gene_est), color = "gray60", alpha = 0.6, size = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = .data$fit), color = "red", size = 1) +
        ggplot2::geom_point(ggplot2::aes(y = .data$final), color = "dodgerblue", alpha = 0.6, size = 0.5) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
            title = title,
            x = "Mean of Normalized Counts",
            y = "Dispersion"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = ggplot2::element_text(size = 12, face = "bold"),
            axis.text = ggplot2::element_text(size = 10)
        )
    
    return(p)
}

#' Create MA plot data
#' @param deseq_result DESeq2 result object
#' @return MA plot data
create_ma_data <- function(deseq_result) {
    # Create data frame
    plot_data <- data.frame(
        mean = deseq_result$baseMean,
        log2fc = deseq_result$log2FoldChange,
        padj = deseq_result$padj,
        stringsAsFactors = FALSE
    )
    
    # Remove rows with missing values
    plot_data <- plot_data[complete.cases(plot_data), ]
    
    # Add significance information
    plot_data$significant <- ifelse(plot_data$padj < 0.05, "Significant", "Not Significant")
    
    return(plot_data)
}

#' Create MA plot
#' @param plot_data MA plot data
#' @param title Plot title
#' @return ggplot object
create_ma_plot <- function(plot_data, title) {
    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$mean, y = .data$log2fc)) +
        ggplot2::geom_point(ggplot2::aes(color = .data$significant), alpha = 0.6, size = 0.5) +
        ggplot2::scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray60")) +
        ggplot2::scale_x_log10() +
        ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        ggplot2::labs(
            title = title,
            x = "Mean of Normalized Counts",
            y = "Log2 Fold Change"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = ggplot2::element_text(size = 12, face = "bold"),
            axis.text = ggplot2::element_text(size = 10),
            legend.title = ggplot2::element_text(size = 10, face = "bold"),
            legend.text = ggplot2::element_text(size = 9)
        )
    
    return(p)
}

#' Create p-value histogram data
#' @param deseq_result DESeq2 result object
#' @return P-value histogram data
create_pvalue_histogram_data <- function(deseq_result) {
    # Get p-values
    pvalues <- deseq_result$pvalue
    
    # Remove missing values
    pvalues <- pvalues[!is.na(pvalues)]
    
    return(pvalues)
}

#' Create p-value histogram
#' @param pvalues P-values
#' @param title Plot title
#' @return ggplot object
create_pvalue_histogram <- function(pvalues, title) {
    # Create data frame
    plot_data <- data.frame(pvalue = pvalues)
    
    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$pvalue)) +
        ggplot2::geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
        ggplot2::labs(
            title = title,
            x = "P-value",
            y = "Frequency"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = ggplot2::element_text(size = 12, face = "bold"),
            axis.text = ggplot2::element_text(size = 10)
        )
    
    return(p)
}

#' Create Cook's distance plot data
#' @param dds DESeq2 dataset
#' @return Cook's distance plot data
create_cooks_distance_data <- function(dds) {
    # Get Cook's distances
    cooks_distances <- assays(dds)[["cooks"]]
    
    # Calculate mean Cook's distance for each gene
    mean_cooks <- rowMeans(cooks_distances, na.rm = TRUE)
    
    # Create data frame
    plot_data <- data.frame(
        gene = rownames(dds),
        cooks_distance = mean_cooks,
        stringsAsFactors = FALSE
    )
    
    # Remove rows with missing values
    plot_data <- plot_data[complete.cases(plot_data), ]
    
    return(plot_data)
}

#' Create Cook's distance plot
#' @param plot_data Cook's distance plot data
#' @param title Plot title
#' @return ggplot object
create_cooks_distance_plot <- function(plot_data, title) {
    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = seq_along(.data$cooks_distance), y = .data$cooks_distance)) +
        ggplot2::geom_point(color = "steelblue", alpha = 0.6, size = 0.5) +
        ggplot2::labs(
            title = title,
            x = "Gene Index",
            y = "Mean Cook's Distance"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = ggplot2::element_text(size = 12, face = "bold"),
            axis.text = ggplot2::element_text(size = 10)
        )
    
    return(p)
}

#' Create PCA plot data
#' @param dds DESeq2 dataset
#' @return PCA plot data
create_pca_data <- function(dds) {
    # Check if we have enough samples for PCA
    if (ncol(dds) < 3) {
        stop("Not enough samples for PCA (need at least 3 samples)")
    }
    
    # Get VST transformed counts
    vst_counts <- DESeq2::vst(dds, blind = FALSE)
    
    # Perform PCA
    pca_result <- prcomp(t(assay(vst_counts)))
    
    # Check if we have enough components
    if (ncol(pca_result$x) < 2) {
        stop("Not enough principal components for PCA plot (need at least 2)")
    }
    
    # Get variance explained
    var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
    
    # Create data frame
    plot_data <- data.frame(
        sample = rownames(pca_result$x),
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        stringsAsFactors = FALSE
    )
    
    # Add sample information
    sample_info <- colData(dds)
    plot_data <- cbind(plot_data, as.data.frame(sample_info))
    
    return(list(
        plot_data = plot_data,
        var_explained = var_explained,
        dds = dds
    ))
}

#' Create PCA plot
#' @param pca_data PCA plot data
#' @param title Plot title
#' @return ggplot object
create_pca_plot <- function(pca_data, title) {
    # Get the first variable from the design formula
    design_formula <- DESeq2::design(pca_data$dds)
    condition <- all.vars(design_formula)[1]
    
    # Create the plot
    p <- ggplot2::ggplot(pca_data$plot_data, ggplot2::aes_string(x = "PC1", y = "PC2", color = condition)) +
        ggplot2::geom_point(size = 3, alpha = 0.8) +
        ggplot2::labs(
            title = title,
            x = paste0("PC1 (", round(pca_data$var_explained[1], 1), "%)"),
            y = paste0("PC2 (", round(pca_data$var_explained[2], 1), "%)")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = ggplot2::element_text(size = 12, face = "bold"),
            axis.text = ggplot2::element_text(size = 10),
            legend.title = ggplot2::element_text(size = 10, face = "bold"),
            legend.text = ggplot2::element_text(size = 9)
        )
    
    return(p)
}

#' Create mean-variance relationship plot
#' @param dds DESeq2 dataset
#' @param title Plot title
#' @return ggplot object
create_mean_variance_plot <- function(dds, title) {
    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)
    
    # Calculate mean and variance
    mean_counts <- rowMeans(norm_counts)
    var_counts <- apply(norm_counts, 1, var)
    
    # Create data frame
    plot_data <- data.frame(
        mean = mean_counts,
        variance = var_counts,
        stringsAsFactors = FALSE
    )
    
    # Remove rows with missing values
    plot_data <- plot_data[complete.cases(plot_data), ]
    
    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$mean, y = .data$variance)) +
        ggplot2::geom_point(alpha = 0.6, size = 0.5, color = "steelblue") +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
            title = title,
            x = "Mean of Normalized Counts",
            y = "Variance"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = ggplot2::element_text(size = 12, face = "bold"),
            axis.text = ggplot2::element_text(size = 10)
        )
    
    return(p)
}

#' Create size factors distribution plot
#' @param dds DESeq2 dataset
#' @param title Plot title
#' @return ggplot object
create_size_factors_plot <- function(dds, title) {
    # Get size factors
    size_factors <- sizeFactors(dds)
    
    # Create data frame
    plot_data <- data.frame(
        sample = names(size_factors),
        size_factor = size_factors,
        stringsAsFactors = FALSE
    )
    
    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$sample, y = .data$size_factor)) +
        ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
        ggplot2::labs(
            title = title,
            x = "Sample",
            y = "Size Factor"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = ggplot2::element_text(size = 12, face = "bold"),
            axis.text = ggplot2::element_text(size = 10),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    
    return(p)
}

#' Create sample-to-sample distances heatmap
#' @param dds DESeq2 dataset
#' @param title Plot title
#' @return ggplot object
create_sample_distances_plot <- function(dds, title) {
    # Get VST transformed counts
    vst_counts <- DESeq2::vst(dds, blind = FALSE)
    
    # Calculate sample distances
    sample_dists <- dist(t(assay(vst_counts)))
    sample_dist_matrix <- as.matrix(sample_dists)
    
    # Convert to long format for ggplot
    dist_df <- expand.grid(
        sample1 = rownames(sample_dist_matrix),
        sample2 = colnames(sample_dist_matrix)
    )
    dist_df$distance <- as.vector(sample_dist_matrix)
    
    # Create the plot
    p <- ggplot2::ggplot(dist_df, ggplot2::aes(x = .data$sample1, y = .data$sample2, fill = .data$distance)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient(low = "white", high = "red") +
        ggplot2::labs(
            title = title,
            x = "Sample",
            y = "Sample",
            fill = "Distance"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = ggplot2::element_text(size = 12, face = "bold"),
            axis.text = ggplot2::element_text(size = 10),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    
    return(p)
}

#' Generate comprehensive diagnostic report
#' @param dds DESeq2 dataset
#' @param contrast_name Contrast name
#' @param experiment_name Experiment name
#' @param analysis_type Analysis type
#' @param plots_dir Output directory for plots
#' @return Diagnostic statistics
generate_diagnostic_report <- function(dds, contrast_name, experiment_name, analysis_type, plots_dir) {
    # Get statistics first for diagnostic checks
    dispersions <- mcols(dds)$dispersion
    cooks <- assays(dds)[["cooks"]]
    size_factors <- sizeFactors(dds)
    
    # Diagnostic Checks
    diagnostic_messages <- c()
    
    # 1. Dispersion Checks
    dispersion_mean <- mean(dispersions, na.rm = TRUE)
    dispersion_cv <- sd(dispersions, na.rm = TRUE) / dispersion_mean
    
    if (dispersion_mean > 4) {
        diagnostic_messages <- c(
            diagnostic_messages,
            "WARNING: High mean dispersion (>4). This suggests high variability between replicates."
        )
    }
    if (dispersion_cv > 0.8) {
        diagnostic_messages <- c(
            diagnostic_messages,
            "WARNING: Large variation in dispersions. Some genes may be poorly modeled."
        )
    }
    
    # 2. Cook's Distance Checks
    outlier_threshold <- 4 / ncol(dds) # DESeq2's default
    outlier_proportion <- mean(cooks > outlier_threshold, na.rm = TRUE)
    
    if (outlier_proportion > 0.1) {
        diagnostic_messages <- c(
            diagnostic_messages,
            sprintf(
                "WARNING: High proportion of outliers (%.1f%%). Consider checking sample quality.",
                outlier_proportion * 100
            )
        )
    }
    
    # 3. Size Factor Checks
    size_factor_cv <- sd(size_factors) / mean(size_factors)
    
    if (size_factor_cv > 0.7) {
        diagnostic_messages <- c(
            diagnostic_messages,
            "WARNING: Large variation in size factors. Check for library size differences."
        )
    }
    if (any(size_factors < 0.1) || any(size_factors > 10)) {
        diagnostic_messages <- c(
            diagnostic_messages,
            "WARNING: Extreme size factors detected. Some samples may have very different sequencing depth."
        )
    }
    
    # 4. Sample Count Check
    if (ncol(dds) < 6) {
        diagnostic_messages <- c(
            diagnostic_messages,
            "CAUTION: Low sample count (<6). Statistical power may be limited."
        )
    }
    
    # 5. Zero Count Check
    zero_proportions <- rowMeans(counts(dds) == 0)
    high_zero_genes <- mean(zero_proportions > 0.75)
    
    if (high_zero_genes > 0.30) {
        diagnostic_messages <- c(
            diagnostic_messages,
            sprintf(
                "WARNING: %.1f%% of genes have zeros in >75%% of samples. Consider filtering.",
                high_zero_genes * 100
            )
        )
    }
    
    # Print diagnostic messages
    if (length(diagnostic_messages) > 0) {
        message("\nDiagnostic Results for ", contrast_name, ":")
        message(paste(diagnostic_messages, collapse = "\n"))
    } else {
        message("\nAll diagnostics look good for ", contrast_name)
    }
    
    # Get results for additional statistics
    res <- results(dds)
    pvals <- res$pvalue[!is.na(res$pvalue)]
    
    # Check p-value distribution
    null_pvals <- pvals[pvals >= 0.5 & pvals <= 1]
    ks_test <- ks.test(null_pvals, "punif", min = 0.5, max = 1)
    
    # Calculate R-squared using fitted values and observed counts
    calculate_rsquared <- function(dds) {
        # Get raw counts and fitted values
        raw_counts <- counts(dds, normalized = FALSE)
        fit <- fitted(dds)
        size_factors <- sizeFactors(dds)
        
        # Calculate R-squared for each gene
        rsq <- sapply(seq_len(nrow(raw_counts)), function(i) {
            observed <- raw_counts[i, ]
            predicted <- fit[i, ]
            
            # Only calculate for genes with sufficient expression and variation
            if (mean(observed) > 10 && var(observed) > 0) {
                # Calculate pseudo R-squared using deviance residuals
                null_model <- mean(observed / size_factors)
                null_dev <- sum((observed - null_model * size_factors)^2)
                model_dev <- sum((observed - predicted)^2)
                
                if (null_dev > 0) {
                    r2 <- 1 - (model_dev / null_dev)
                    # Bound R-squared between 0 and 1
                    return(max(0, min(1, r2)))
                }
            }
            return(NA)
        })
        
        return(rsq[!is.infinite(rsq) & !is.nan(rsq)])
    }
    
    # Check dispersion trend
    check_dispersion_trend <- function(dds) {
        # Get mean counts and dispersions
        mean_counts <- rowMeans(counts(dds, normalized = TRUE))
        dispersions <- dispersions(dds)
        
        # Calculate correlation (should be negative for proper trend)
        cor_coef <- cor(log10(mean_counts), log10(dispersions), method = "spearman")
        
        # Check if trend is decreasing
        is_decreasing <- cor_coef < 0
        
        return(list(
            correlation = cor_coef,
            is_decreasing = is_decreasing
        ))
    }
    
    # Get statistics
    r_squared_values <- calculate_rsquared(dds)
    dispersion_trend <- check_dispersion_trend(dds)
    
    stats <- list(
        dispersion = list(
            mean = mean(dispersions, na.rm = TRUE),
            median = median(dispersions, na.rm = TRUE),
            sd = sd(dispersions, na.rm = TRUE),
            min = min(dispersions, na.rm = TRUE),
            max = max(dispersions, na.rm = TRUE)
        ),
        cooks = list(
            mean = mean(cooks, na.rm = TRUE),
            median = median(cooks, na.rm = TRUE),
            sd = sd(cooks, na.rm = TRUE),
            min = min(cooks, na.rm = TRUE),
            max = max(cooks, na.rm = TRUE),
            outliers = sum(cooks > 4 / ncol(dds), na.rm = TRUE)
        ),
        size_factors = list(
            values = size_factors,
            mean = mean(size_factors),
            sd = sd(size_factors),
            cv = sd(size_factors) / mean(size_factors)
        ),
        independent_filtering = if (!is.null(metadata(dds)$filterThreshold)) {
            list(
                threshold = metadata(dds)$filterThreshold,
                filtered_fraction = mean(mcols(dds)$baseMean < metadata(dds)$filterThreshold)
            )
        } else {
            NULL
        },
        diagnostic_checks = list(
            dispersion_mean = dispersion_mean,
            dispersion_cv = dispersion_cv,
            outlier_proportion = outlier_proportion,
            size_factor_cv = size_factor_cv,
            high_zero_proportion = high_zero_genes,
            warnings = diagnostic_messages
        ),
        pvalue_distribution = list(
            n_tests = length(pvals),
            prop_significant = mean(pvals < 0.05, na.rm = TRUE),
            ks_test_pvalue = ks_test$p.value,
            uniform_under_null = ks_test$p.value >= 0.05
        ),
        model_fit = list(
            median_rsquared = median(r_squared_values, na.rm = TRUE),
            mean_rsquared = mean(r_squared_values, na.rm = TRUE),
            prop_high_rsquared = mean(r_squared_values > 0.5, na.rm = TRUE),
            zero_count_proportion = mean(counts(dds) == 0),
            extreme_cooks_count = sum(apply(assays(dds)[["cooks"]], 1, function(x) any(x > 4 / ncol(dds))))
        ),
        dispersion_trend = dispersion_trend
    )
    
    # Save statistics
    saveRDS(
        stats,
        file.path(plots_dir, paste0(experiment_name, "_", contrast_name, "_stats.rds"))
    )
    
    # Create summary document
    summary_text <- c(
        "=== DESeq2 Model Diagnostics ===\n",
        sprintf("Experiment: %s", experiment_name),
        sprintf("Contrast: %s\n", contrast_name),
        "Model Description:",
        "DESeq2 uses a negative binomial model for differential expression analysis.",
        "Input Requirements:",
        "- Raw (unnormalized) integer counts from RNA-seq",
        "- No prior normalization (e.g., TPM, RPKM, or FPKM) should be applied",
        "- Counts should represent the number of reads/fragments mapped to each gene",
        "Key model characteristics:",
        "- Variance is modeled using gene-wise dispersion estimates",
        "- Empirical Bayes shrinkage is applied to dispersion estimates",
        "- Log fold changes are shrunk using empirical Bayes to reduce noise",
        "- Size factors normalize for sequencing depth differences",
        "- Independent filtering optimizes power by removing low-count genes",
        "- Maximum likelihood estimation (MLE) is used for initial parameter estimates",
        "- Wald test is used for significance testing of coefficients",
        "- Cook's distance identifies outlier samples",
        "- Adjusted p-values control for multiple testing using Benjamini-Hochberg method\n",
        if (analysis_type == "complex") {
            c(
                "Complex Model Additional Characteristics:",
                "- Includes interaction terms between conditions and timepoints",
                "- Models time-dependent drug effects",
                "- Accounts for dose-response relationships",
                "- Enables testing of complex hypotheses about treatment effects",
                "- Allows for nested comparisons within treatment groups",
                "- Provides more sophisticated variance modeling across conditions",
                sprintf("- Full model formula: %s", deparse(formula(design(dds)))),
                "- Interaction contrasts use custom coefficient combinations",
                "- Handles multi-level factorial designs\n"
            )
        } else {
            NULL
        }
    )
    
    if (!is.null(metadata(dds)$filterThreshold)) {
        summary_text <- c(
            summary_text,
            "\nIndependent Filtering:",
            sprintf("- Threshold: %.2f", stats$independent_filtering$threshold),
            sprintf(
                "- Fraction filtered: %.2f%%",
                stats$independent_filtering$filtered_fraction * 100
            )
        )
    }
    
    # Add warnings to summary document
    if (length(diagnostic_messages) > 0) {
        summary_text <- c(
            summary_text,
            "\nDiagnostic Warnings:",
            diagnostic_messages
        )
    } else {
        summary_text <- c(
            summary_text,
            "\nNo diagnostic warnings - all checks passed."
        )
    }
    
    # Add p-value distribution check to summary text
    summary_text <- c(
        summary_text,
        "\nP-value Distribution Check:",
        sprintf("- Number of tests: %d", length(pvals)),
        sprintf(
            "- Proportion of p-values < 0.05: %.1f%%",
            mean(pvals < 0.05, na.rm = TRUE) * 100
        ),
        sprintf(
            "- Uniformity test of null p-values (0.5-1.0): p = %.3e",
            ks_test$p.value
        ),
        if (ks_test$p.value < 0.05) {
            "WARNING: P-value distribution may not be uniform under the null"
        } else {
            "- P-value distribution appears uniform under the null (as expected)"
        }
    )
    
    # Add dispersion trend check to summary text
    summary_text <- c(
        summary_text,
        "\nDispersion Trend Check:",
        if (dispersion_trend$is_decreasing) {
            "- Dispersion estimates follow expected decreasing trend with mean"
        } else {
            "WARNING: Dispersion estimates do not follow expected decreasing trend with mean"
        },
        sprintf("- Correlation coefficient: %.3f", dispersion_trend$correlation)
    )
    
    # Add model fit statistics to summary text
    summary_text <- c(
        summary_text,
        "\nModel Fit Statistics:",
        sprintf(
            "- Median R-squared of variance fit: %.3f",
            median(r_squared_values, na.rm = TRUE)
        ),
        sprintf(
            "- Mean R-squared of variance fit: %.3f",
            mean(r_squared_values, na.rm = TRUE)
        ),
        sprintf(
            "- Proportion of genes with R-squared > 0.5: %.1f%%",
            mean(r_squared_values > 0.5, na.rm = TRUE) * 100
        ),
        sprintf(
            "- Proportion of zero counts: %.1f%%",
            mean(counts(dds) == 0) * 100
        ),
        sprintf(
            "- Number of genes with extreme Cook's distances: %d",
            sum(apply(assays(dds)[["cooks"]], 1, function(x) any(x > 4 / ncol(dds))))
        )
    )
    
    # Create and save summary document
    ft <- flextable::flextable(data.frame(Summary = summary_text)) %>%
        flextable::delete_part(part = "header") %>%
        flextable::border_remove() %>%
        flextable::padding(padding = 0) %>%
        flextable::autofit()
    
    flextable::save_as_docx(ft,
        path = file.path(
            plots_dir,
            paste0(
                experiment_name, "_",
                contrast_name, "_diagnostics.docx"
            )
        )
    )
    
    return(stats)
}

#' Generate diagnostic plots for simple analysis
#' @param dea_results DESeq2 results
#' @param experiment_obj Experiment object
#' @param plot_types Types of plots to generate
#' @param plots_dir Output directory for plots
#' @param analysis_type Analysis type
#' @param experiment_name Experiment name
generate_simple_diagnostics <- function(dea_results, experiment_obj, plot_types, plots_dir, analysis_type, experiment_name) {
    for (contrast_name in names(dea_results)) {
        message("\nGenerating diagnostic plots for contrast: ", contrast_name)
        
        # Get the result object
        result <- dea_results[[contrast_name]]
        
        # Generate plots
        generate_contrast_diagnostics(result, contrast_name, plot_types, plots_dir, analysis_type, experiment_name)
    }
}

#' Generate diagnostic plots for complex analysis
#' @param dea_results DESeq2 results
#' @param experiment_obj Experiment object
#' @param plot_types Types of plots to generate
#' @param plots_dir Output directory for plots
#' @param analysis_type Analysis type
#' @param experiment_name Experiment name
generate_complex_diagnostics <- function(dea_results, experiment_obj, plot_types, plots_dir, analysis_type, experiment_name) {
    message("Generating diagnostic plots for complex analysis...")
    
    # Get the dds object
    dds <- dea_results$dds
    
    # Generate plots for the overall analysis
    generate_contrast_diagnostics(list(dds = dds), "overall", plot_types, plots_dir, analysis_type, experiment_name)
    
    # Generate plots for each contrast
    for (contrast_name in names(dea_results$deseq_results)) {
        message("\nGenerating diagnostic plots for contrast: ", contrast_name)
        
        # Create a result object with the dds and the specific contrast result
        result <- list(
            dds = dds,
            deseq_result = dea_results$deseq_results[[contrast_name]]
        )
        
        # Generate plots
        generate_contrast_diagnostics(result, contrast_name, plot_types, plots_dir, analysis_type, experiment_name)
    }
}

#' Generate diagnostic plots for a specific contrast
#' @param result Result object
#' @param contrast_name Contrast name
#' @param plot_types Types of plots to generate
#' @param plots_dir Output directory for plots
#' @param analysis_type Analysis type
#' @param experiment_name Experiment name
generate_contrast_diagnostics <- function(result, contrast_name, plot_types, plots_dir, analysis_type, experiment_name) {
    # Generate each requested plot type
    for (plot_type in plot_types) {
        tryCatch({
            generate_single_diagnostic_plot(result, contrast_name, plot_type, plots_dir, analysis_type)
        }, error = function(e) {
            message("Error generating ", plot_type, " plot for ", contrast_name, ": ", e$message)
        })
    }
    
    # Generate comprehensive diagnostic report
    tryCatch({
        generate_diagnostic_report(result$dds, contrast_name, experiment_name, analysis_type, plots_dir)
    }, error = function(e) {
        message("Error generating diagnostic report for ", contrast_name, ": ", e$message)
    })
}

#' Generate a single diagnostic plot
#' @param result Result object
#' @param contrast_name Contrast name
#' @param plot_type Type of plot to generate
#' @param plots_dir Output directory for plots
#' @param analysis_type Analysis type
generate_single_diagnostic_plot <- function(result, contrast_name, plot_type, plots_dir, analysis_type) {
    # Create filename
    filename <- file.path(plots_dir, paste0(contrast_name, "_", plot_type, ".png"))
    
    # Debug: Check if filename is valid
    if (is.null(filename) || length(filename) == 0 || filename == "") {
        stop("Invalid filename generated: '", filename, "'")
    }
    
    # Generate plot based on type
    plot_obj <- switch(plot_type,
        "dispersion" = {
            plot_data <- create_dispersion_data(result$dds)
            create_dispersion_plot(plot_data, paste("Dispersion Plot:", contrast_name))
        },
        "ma" = {
            plot_data <- create_ma_data(result$deseq_result)
            create_ma_plot(plot_data, paste("MA Plot:", contrast_name))
        },
        "pvalue_histogram" = {
            pvalues <- create_pvalue_histogram_data(result$deseq_result)
            create_pvalue_histogram(pvalues, paste("P-value Histogram:", contrast_name))
        },
        "cooks_distance" = {
            plot_data <- create_cooks_distance_data(result$dds)
            create_cooks_distance_plot(plot_data, paste("Cook's Distance Plot:", contrast_name))
        },
        "pca" = {
            pca_data <- create_pca_data(result$dds)
            create_pca_plot(pca_data, paste("PCA Plot:", contrast_name))
        },
        "mean_variance" = {
            create_mean_variance_plot(result$dds, paste("Mean-Variance Plot:", contrast_name))
        },
        "size_factors" = {
            create_size_factors_plot(result$dds, paste("Size Factors Plot:", contrast_name))
        },
        "sample_distances" = {
            create_sample_distances_plot(result$dds, paste("Sample Distances Plot:", contrast_name))
        },
        stop("Unknown plot type: ", plot_type)
    )
    
    # Save plot as PNG
    ggplot2::ggsave(filename, plot_obj, 
                   width = diagnostics_config$default_plot_params$width,
                   height = diagnostics_config$default_plot_params$height,
                   dpi = diagnostics_config$default_plot_params$dpi)
    
    message("Saved ", plot_type, " plot: ", filename)
}

#' Main diagnostics function with comprehensive workflow
#' @param experiment_name Name of experiment
#' @param analysis_type Type of analysis
#' @param plot_types Types of plots to generate
#' @return Experiment object
generate_deseq2_diagnostics <- function(experiment_name,
                                      analysis_type = "simple",
                                      plot_types = diagnostics_config$default_plot_types) {
    message("Starting DESeq2 diagnostics generation...")
    
    # Validate parameters
    validate_diagnostics_params(experiment_name, analysis_type, plot_types)
    
    # Load experiment data
    experiment_obj <- load_experiment_data(experiment_name)
    
    # Load DESeq2 results
    dea_results <- load_deseq2_results(experiment_name, analysis_type)
    
    # Create output directories
    dirs <- create_diagnostics_output_directories(experiment_name, analysis_type)
    
    # Generate diagnostic plots based on analysis type
    if (analysis_type == "complex") {
        generate_complex_diagnostics(dea_results, experiment_obj, plot_types, dirs$plots_dir, analysis_type, experiment_name)
    } else {
        generate_simple_diagnostics(dea_results, experiment_obj, plot_types, dirs$plots_dir, analysis_type, experiment_name)
    }
    
    message("\nDiagnostics generation completed!")
    return(experiment_obj)
}

message("DESeq2 diagnostics functions loaded successfully. Use generate_deseq2_diagnostics() to generate diagnostic plots.")
