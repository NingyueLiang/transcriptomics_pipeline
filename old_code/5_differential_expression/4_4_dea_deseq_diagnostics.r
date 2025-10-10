#####
# DESeq2 Model Diagnostics Functions

# Required libraries
library(DESeq2) # For accessing DESeq2 objects and diagnostics
library(dplyr) # For data manipulation
library(magrittr) # For pipe operations
library(flextable) # For creating formatted tables
library(officer) # For saving Word documents
library(grid) # For graphics
library(vsn) # For meanSdPlot
library(pheatmap) # For heatmaps
library(hexbin)
library(ggplot2) # For ggplot2
library(S4Vectors) # For mcols function
library(BiocGenerics) # For basic Bioconductor functions
library(SummarizedExperiment) # For assays function
library(GenomicRanges) # For genomic ranges functionality

# Look at deseq2 dispersion and cooks distance
plot_deseq2_diagnostics <- function(dds, experiment_name, contrast_name, analysis_type = "simple") {
    # Check if dds object is valid
    if (is.null(dds)) {
        stop("DESeq2 object is NULL")
    }

    if (!is(dds, "DESeqDataSet")) {
        stop("Object is not a DESeqDataSet. Got class: ", class(dds))
    }

    # Check if the object has been run through DESeq
    if (is.null(mcols(dds)$dispersion)) {
        message("No dispersion values found. Running DESeq...")
        dds <- DESeq(dds)
    }

    # Create results directory under analysis_type/diagnostics
    plot_dir <- file.path("results", experiment_name, "de", "deseq2", analysis_type, "diagnostics")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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

    # 1. Dispersion Plot
    pdf(file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_dispersion.pdf")))
    plotDispEsts(dds, main = "Dispersion Estimates")
    dev.off()

    # 2. Cook's distance plot
    pdf(file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_cooks.pdf")))
    par(mar = c(8, 5, 2, 2))
    boxplot(log10(assays(dds)[["cooks"]]),
        range = 0,
        las = 2,
        main = "Cook's distance per sample",
        ylab = "log10(Cook's distance)"
    )
    dev.off()

    # 3. Mean-Variance relationship plot
    pdf(file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_mean_variance.pdf")))
    vsn::meanSdPlot(log2(counts(dds, normalized = TRUE) + 1),
        ranks = FALSE,
        main = "Mean-Variance Relationship"
    )
    dev.off()

    # 4. NEW: Size factors distribution plot
    pdf(file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_size_factors.pdf")))
    barplot(sizeFactors(dds),
        main = "Size Factors Distribution",
        las = 2,
        ylab = "Size Factor"
    )
    abline(h = 1, col = "red", lty = 2)
    dev.off()

    # 5. PCA of variance stabilized data
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

    # Get the first column that isn't 'sample_id'
    group_col <- colnames(colData(dds))[!colnames(colData(dds)) %in% "sample_id"][1]

    pdf(file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_pca.pdf")))

    # Create PCA plot with proper grouping
    pcaData <- plotPCA(vsd, intgroup = group_col, returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    pca_plot <- ggplot(pcaData, aes_string(x = "PC1", y = "PC2", color = group_col)) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_bw() +
        ggtitle(paste("PCA Plot -", contrast_name))

    print(pca_plot)
    dev.off()

    # 6. NEW: Sample-to-sample distances heatmap
    pdf(file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_sample_distances.pdf")))
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    pheatmap::pheatmap(sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        main = "Sample Distances"
    )
    dev.off()

    # Create MA plot
    pdf(file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_MA_plot.pdf")))
    res <- results(dds)
    plotMA(res, ylim = c(-5, 5))
    dev.off()

    # Create P-value distribution plot
    pdf(file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_pvalue_dist.pdf")))
    hist(res$pvalue[!is.na(res$pvalue)],
        breaks = 50,
        main = "P-value Distribution",
        xlab = "P-value",
        ylab = "Frequency"
    )
    dev.off()

    # Check p-value distribution
    pvals <- res$pvalue[!is.na(res$pvalue)]
    # Check uniformity of p-values between 0.5 and 1 (where we expect null hypothesis to be true)
    null_pvals <- pvals[pvals >= 0.5 & pvals <= 1]
    ks_test <- ks.test(null_pvals, "punif", min = 0.5, max = 1)

    # 3. Get summary statistics
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
            outliers = sum(cooks > 4 / ncol(dds), na.rm = TRUE) # DESeq2's default threshold
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
        )
    )

    # Save statistics
    saveRDS(
        stats,
        file.path(plot_dir, paste0(experiment_name, "_", contrast_name, "_stats.rds"))
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

    # Calculate R-squared using fitted values and observed counts
    calculate_rsquared <- function(dds) {
        # Get raw counts and fitted values
        raw_counts <- counts(dds, normalized = FALSE)
        fit <- fitted(dds)
        size_factors <- sizeFactors(dds)

        # Calculate R-squared for each gene
        rsq <- sapply(1:nrow(raw_counts), function(i) {
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
        dispersions <- dispersions(dds) # Using proper DESeq2 accessor function

        # Calculate correlation (should be negative for proper trend)
        cor_coef <- cor(log10(mean_counts), log10(dispersions), method = "spearman")

        # Check if trend is decreasing
        is_decreasing <- cor_coef < 0

        return(list(
            correlation = cor_coef,
            is_decreasing = is_decreasing
        ))
    }

    # Add dispersion trend check to summary text
    dispersion_trend <- check_dispersion_trend(dds)

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
    r_squared_values <- calculate_rsquared(dds)

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
            plot_dir,
            paste0(
                experiment_name, "_",
                contrast_name, "_diagnostics.docx"
            )
        )
    )

    return(stats)
}

# Function to run diagnostics for both models
run_model_diagnostics <- function(experiment_name, simple_results = NULL, complex_results = NULL) {
    message("Starting DESeq2 model diagnostics...")

    # Load experiment data
    experiment_obj <- readRDS(file.path(
        "data", experiment_name,
        paste0(experiment_name, ".rds")
    ))

    # Load results if not provided
    if (is.null(simple_results)) {
        simple_results_path <- file.path(
            "results", experiment_name, "de", "deseq2",
            "simple", paste0(experiment_name, "_de_results.rds")
        )
        if (file.exists(simple_results_path)) {
            message("Loading simple model results from: ", simple_results_path)
            simple_results <- readRDS(simple_results_path)
        }
    }

    if (is.null(complex_results)) {
        complex_results_path <- file.path(
            "results", experiment_name, "de", "deseq2",
            "complex", paste0(experiment_name, "_de_complex_results.rds")
        )
        if (file.exists(complex_results_path)) {
            message("Loading complex model results from: ", complex_results_path)
            complex_results <- readRDS(complex_results_path)
        }
    }

    # Process simple model results
    if (!is.null(simple_results)) {
        message("\nProcessing simple model contrasts:")
        for (contrast_name in names(simple_results)) {
            message("- Processing contrast: ", contrast_name)
            if (!is.null(simple_results[[contrast_name]]$dds)) {
                dds <- simple_results[[contrast_name]]$dds
                diagnostics <- plot_deseq2_diagnostics(dds, experiment_name, contrast_name, "simple")
            } else {
                warning("No DESeq object found for simple model contrast: ", contrast_name)
            }
        }
    }

    # Process complex model results
    if (!is.null(complex_results)) {
        message("\nProcessing complex model:")
        if (!is.null(complex_results$dds)) {
            dds <- complex_results$dds

            # Run overall model diagnostics
            diagnostics <- plot_deseq2_diagnostics(dds, experiment_name, "complex_model", "complex")

            # Run diagnostics for each contrast
            message("Processing individual contrasts in complex model:")
            for (contrast_name in names(complex_results$deseq_results)) {
                message("- Processing contrast: ", contrast_name)
                diagnostics <- plot_deseq2_diagnostics(
                    dds, experiment_name,
                    paste0("complex_", contrast_name), "complex"
                )
            }
        } else {
            warning("No DESeq object found in complex model results")
        }
    }

    message("\nDiagnostics complete!")
    return(experiment_obj)
}

# Example usage of the combined function:
xen_tran_2024_03 <- run_model_diagnostics(
    experiment_name = "xen_tran_2024_03"
)

# Example usage for xen_tran_2024_12:
xen_tran_2024_12 <- run_model_diagnostics(
    experiment_name = "xen_tran_2024_12"
)
