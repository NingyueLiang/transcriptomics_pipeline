#' DESeq2 Diagnostics Functions
#' 
#' This module provides functions for generating comprehensive diagnostic plots
#' and quality control metrics for DESeq2 analysis with error handling and validation.

# Required libraries
suppressPackageStartupMessages({
    library(DESeq2)         # For DESeq2 analysis
    library(ggplot2)        # For creating plots
    library(dplyr)          # For data manipulation
    library(tidyr)          # For data tidying
    library(here)           # For path management
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
    default_plot_types = c("dispersion", "ma", "pvalue_histogram", "cooks_distance", "pca")
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

#' Create output directories
#' @param experiment_name Name of experiment
#' @param analysis_type Type of analysis
#' @return List of directory paths
create_output_directories <- function(experiment_name, analysis_type) {
    base_dir <- here("results", experiment_name, "de", "deseq2", analysis_type)
    plots_dir <- file.path(base_dir, "plots", "diagnostics")
    
    if (!dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Failed to create plots directory: ", plots_dir)
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
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mean)) +
        ggplot2::geom_point(ggplot2::aes(y = gene_est), color = "gray60", alpha = 0.6, size = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = fit), color = "red", size = 1) +
        ggplot2::geom_point(ggplot2::aes(y = final), color = "dodgerblue", alpha = 0.6, size = 0.5) +
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
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mean, y = log2fc)) +
        ggplot2::geom_point(ggplot2::aes(color = significant), alpha = 0.6, size = 0.5) +
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
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = pvalue)) +
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
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = seq_along(cooks_distance), y = cooks_distance)) +
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
    # Get VST transformed counts
    vst_counts <- DESeq2::vst(dds, blind = FALSE)
    
    # Perform PCA
    pca_result <- prcomp(t(assay(vst_counts)))
    
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
        var_explained = var_explained
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

#' Generate diagnostic plots for simple analysis
#' @param dea_results DESeq2 results
#' @param experiment_obj Experiment object
#' @param plot_types Types of plots to generate
#' @param plots_dir Output directory for plots
#' @param analysis_type Analysis type
generate_simple_diagnostics <- function(dea_results, experiment_obj, plot_types, plots_dir, analysis_type) {
    for (contrast_name in names(dea_results)) {
        message("\nGenerating diagnostic plots for contrast: ", contrast_name)
        
        # Get the result object
        result <- dea_results[[contrast_name]]
        
        # Generate plots
        generate_contrast_diagnostics(result, contrast_name, plot_types, plots_dir, analysis_type)
    }
}

#' Generate diagnostic plots for complex analysis
#' @param dea_results DESeq2 results
#' @param experiment_obj Experiment object
#' @param plot_types Types of plots to generate
#' @param plots_dir Output directory for plots
#' @param analysis_type Analysis type
generate_complex_diagnostics <- function(dea_results, experiment_obj, plot_types, plots_dir, analysis_type) {
    message("Generating diagnostic plots for complex analysis...")
    
    # Get the dds object
    dds <- dea_results$dds
    
    # Generate plots for the overall analysis
    generate_contrast_diagnostics(list(dds = dds), "overall", plot_types, plots_dir, analysis_type)
    
    # Generate plots for each contrast
    for (contrast_name in names(dea_results$deseq_results)) {
        message("\nGenerating diagnostic plots for contrast: ", contrast_name)
        
        # Create a result object with the dds and the specific contrast result
        result <- list(
            dds = dds,
            deseq_result = dea_results$deseq_results[[contrast_name]]
        )
        
        # Generate plots
        generate_contrast_diagnostics(result, contrast_name, plot_types, plots_dir, analysis_type)
    }
}

#' Generate diagnostic plots for a specific contrast
#' @param result Result object
#' @param contrast_name Contrast name
#' @param plot_types Types of plots to generate
#' @param plots_dir Output directory for plots
#' @param analysis_type Analysis type
generate_contrast_diagnostics <- function(result, contrast_name, plot_types, plots_dir, analysis_type) {
    # Generate each requested plot type
    for (plot_type in plot_types) {
        tryCatch({
            generate_single_diagnostic_plot(result, contrast_name, plot_type, plots_dir, analysis_type)
        }, error = function(e) {
            message("Error generating ", plot_type, " plot for ", contrast_name, ": ", e$message)
        })
    }
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
        stop("Unknown plot type: ", plot_type)
    )
    
    # Save plot
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
    dirs <- create_output_directories(experiment_name, analysis_type)
    
    # Generate diagnostic plots based on analysis type
    if (analysis_type == "complex") {
        generate_complex_diagnostics(dea_results, experiment_obj, plot_types, dirs$plots_dir, analysis_type)
    } else {
        generate_simple_diagnostics(dea_results, experiment_obj, plot_types, dirs$plots_dir, analysis_type)
    }
    
    message("\nDiagnostics generation completed!")
    return(experiment_obj)
}

message("DESeq2 diagnostics functions loaded successfully. Use generate_deseq2_diagnostics() to generate diagnostic plots.")
