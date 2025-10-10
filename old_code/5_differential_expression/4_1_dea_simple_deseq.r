#####
# Simple DESeq2 Analysis Functions

# Required libraries
library(DESeq2) # For differential expression analysis
library(dplyr) # For data manipulation
library(magrittr) # For pipe operations
library(flextable) # For creating formatted tables
library(officer) # For saving Word documents
library(stringr) # For string manipulation

# Run differential expression analysis using DESeq2
prepare_comparison_data <- function(counts_matrix, sample_info, comparison) {
    cat("\nPreparing comparison for variable:", comparison$variable, "\n")

    # Initial droplevels
    sample_info <- droplevels(sample_info)

    # Apply pattern-based filtering if specified
    if (!is.null(comparison$filter_patterns)) {
        cat("Applying pattern filters:\n")
        keep_samples <- rep(TRUE, nrow(sample_info))

        for (pattern in comparison$filter_patterns) {
            if (pattern$action == "remove") {
                remove_idx <- grepl(pattern$pattern, sample_info$sample_id,
                    ignore.case = pattern$ignore_case
                )
                keep_samples <- keep_samples & !remove_idx
                cat("- Removing samples where ID matches", pattern$pattern, "\n")
            }
        }

        sample_info <- sample_info[keep_samples, ]
        counts_matrix <- counts_matrix[, sample_info$sample_id]
        sample_info <- droplevels(sample_info)
    }

    # Clean any potential whitespace in relevant columns
    if (comparison$variable %in% names(sample_info)) {
        sample_info[[comparison$variable]] <- trimws(sample_info[[comparison$variable]])
    }

    # If the variable contains multiple components (e.g., drug_treatment_sample_prep)
    if (grepl("_", comparison$variable)) {
        possible_cols <- names(sample_info)
        cols <- c()

        # Find the matching columns from sample_info
        for (col in possible_cols) {
            if (grepl(col, comparison$variable, fixed = TRUE)) {
                cols <- c(cols, col)
            }
        }

        if (length(cols) == 0) {
            stop("Could not find matching columns for: ", comparison$variable)
        }

        # Clean whitespace in component columns
        for (col in cols) {
            sample_info[[col]] <- trimws(sample_info[[col]])
        }

        # Create combined variable
        sample_info[[comparison$variable]] <- apply(sample_info[cols], 1, paste, collapse = "_")
    } else {
        # For single column comparisons, verify the column exists
        if (!comparison$variable %in% names(sample_info)) {
            stop("Column not found: ", comparison$variable)
        }
    }

    # Keep only experimental and control groups - do this before factor conversion
    keep_groups <- sample_info[[comparison$variable]] %in%
        c(comparison$experimental, comparison$control)
    sample_info <- sample_info[keep_groups, ]
    sample_info <- droplevels(sample_info) # Additional droplevels after filtering

    # Set factor levels using simple factor conversion like the old script
    sample_info[[comparison$variable]] <- factor(
        sample_info[[comparison$variable]],
        levels = c(comparison$control, comparison$experimental)
    )

    # Sort samples by ID to ensure consistent order
    sample_ids <- sort(sample_info$sample_id)
    sample_info <- sample_info[match(sample_ids, sample_info$sample_id), ]

    # Get filtered counts in the same order as sample_info
    counts_filtered <- counts_matrix[, sample_ids]

    # Print summary
    cat(
        "\nGroups being compared:",
        paste(levels(sample_info[[comparison$variable]]), collapse = " vs "), "\n"
    )
    cat("Number of samples per group:\n")
    print(table(sample_info[[comparison$variable]]))

    return(list(
        counts = counts_filtered,
        sample_info = sample_info
    ))
}

# Function to automatically detect and create metadata mappings
create_metadata_mappings <- function(metadata_df) {
    # Get all unique columns except sample_id and any numeric columns
    metadata_cols <- names(metadata_df)[!names(metadata_df) %in% c("sample_id") &
        !sapply(metadata_df, is.numeric)]

    # Initialize mappings list
    mappings <- list()

    # For each metadata column
    for (col in metadata_cols) {
        # Get unique values
        unique_values <- unique(metadata_df[[col]])

        # Create mapping for each unique value
        col_mapping <- list()
        for (val in unique_values) {
            # Create potential variations of the value
            variations <- c(
                val,
                gsub(" ", "_", val),
                gsub(" ", ".", val),
                gsub(" ", "", val),
                substr(val, 1, 1) # First letter abbreviation
            )

            # Add to mapping
            col_mapping[[val]] <- unique(variations)
        }

        # Add to main mappings
        mappings[[col]] <- col_mapping
    }

    return(mappings)
}

# Filter genes based on count threshold
filter_low_counts <- function(counts_matrix,
                              min_count = 3,
                              min_sample_percent = 0.2) {
    sample_threshold <- ncol(counts_matrix) * min_sample_percent
    genes_to_keep <- apply(counts_matrix, 1, function(x) {
        sum(x >= min_count) >= sample_threshold
    })

    filtered_counts <- counts_matrix[genes_to_keep, ]

    filter_stats <- list(
        original_genes = nrow(counts_matrix),
        retained_genes = nrow(filtered_counts),
        removed_genes = nrow(counts_matrix) - nrow(filtered_counts),
        parameters = list(
            min_count = min_count,
            min_sample_percent = min_sample_percent
        )
    )

    return(list(
        counts = filtered_counts,
        stats = filter_stats
    ))
}

save_analysis_summary <- function(all_results, experiment_name) {
    # Create summary text
    summary_text <- c(
        "=== Differential Expression Analysis Summary ===\n",
        sprintf("Experiment: %s", experiment_name),
        sprintf("Total contrasts analyzed: %d\n", length(all_results)),
        "Contrasts analyzed:"
    )

    for (contrast_name in names(all_results)) {
        result <- all_results[[contrast_name]]
        n_sig <- sum(!is.na(result$deseq_results$padj) & result$deseq_results$padj < 0.05)
        summary_text <- c(
            summary_text,
            sprintf("\n- %s", contrast_name),
            sprintf("  Control: %s", result$comparison_info$control),
            sprintf("  Experimental: %s", result$comparison_info$experimental),
            sprintf("  Significant genes (padj < 0.05): %d", n_sig)
        )
    }

    # Create a flextable
    ft <- flextable::flextable(data.frame(Summary = summary_text)) %>%
        flextable::delete_part(part = "header") %>%
        flextable::border_remove() %>%
        flextable::padding(padding = 0) %>%
        flextable::autofit()

    # Create simple directory if it doesn't exist
    simple_dir <- file.path("results", experiment_name, "de", "deseq2", "simple")
    dir.create(simple_dir, recursive = TRUE, showWarnings = FALSE)

    # Save to docx with updated path
    summary_file <- file.path(simple_dir, paste0(experiment_name, "_de_summary.docx"))
    flextable::save_as_docx(ft, path = summary_file)

    return(summary_file)
}

# 1. Modified helper function for contrasts with interactive selection and reference levels
generate_all_contrasts <- function(metadata_df, comparison_cols = NULL,
                                   reference_levels = NULL,
                                   experimental_levels = NULL,
                                   contrast_filters = NULL,
                                   interactive = TRUE) {
    all_contrasts <- list()

    # Clean strings helper function
    clean_string <- function(x) {
        x <- tolower(trimws(x))
        x <- gsub("\\s+", "_", x) # replace spaces with underscores
        return(x)
    }

    # Helper function to check if contrast meets filter criteria
    should_include_contrast <- function(variables_present, val1, val2) {
        if (!is.null(contrast_filters) && !is.null(contrast_filters$valid_combinations)) {
            required_vars <- contrast_filters$valid_combinations[[1]]

            # For single-variable contrasts
            if (length(variables_present) == 1) {
                cat(
                    "\n[Automatically skipped - requires both",
                    paste(required_vars, collapse = " and "), "]\n"
                )
                return(FALSE)
            }

            # For combined contrasts, check if it has all required variables
            missing_vars <- setdiff(required_vars, variables_present)
            if (length(missing_vars) > 0) {
                cat(
                    "\n[Automatically skipped - missing:",
                    paste(missing_vars, collapse = ", "), "]\n"
                )
                return(FALSE)
            }
        }

        # Additional checks for experimental levels if provided
        if (!is.null(experimental_levels)) {
            for (var in names(experimental_levels)) {
                if (var %in% variables_present) {
                    pattern <- experimental_levels[[var]]
                    if (!any(grepl(pattern, c(val1, val2)))) {
                        cat("\n[Automatically skipped - no matching experimental level for", var, "]\n")
                        return(FALSE)
                    }
                }
            }
        }

        return(TRUE)
    }

    # If no specific columns provided, use all columns except sample_id
    if (is.null(comparison_cols)) {
        comparison_cols <- names(metadata_df)[!names(metadata_df) %in% c("sample_id")]
    }

    # Clean the metadata values
    for (col in comparison_cols) {
        metadata_df[[col]] <- clean_string(metadata_df[[col]])
    }

    # Clean reference levels if provided
    if (!is.null(reference_levels)) {
        reference_levels <- lapply(reference_levels, clean_string)
    }

    # Generate contrasts for individual columns
    for (col in comparison_cols) {
        values <- unique(as.character(metadata_df[[col]]))
        cat(sprintf("\nProcessing contrasts for column: %s\n", col))
        cat("Available values:", paste(values, collapse = ", "), "\n")

        if (!is.null(reference_levels) && !is.null(reference_levels[[col]])) {
            # Get primary and secondary reference levels
            primary_ref <- reference_levels[[col]][1]
            secondary_ref <- if (length(reference_levels[[col]]) > 1) {
                reference_levels[[col]][2]
            } else {
                NULL
            }

            cat(sprintf("\nPrimary reference level: %s\n", primary_ref))
            if (!is.null(secondary_ref)) {
                cat(sprintf("Secondary reference level: %s\n", secondary_ref))
            }

            # Generate all possible pairs
            for (i in 1:(length(values) - 1)) {
                for (j in (i + 1):length(values)) {
                    val1 <- values[i]
                    val2 <- values[j]

                    cat("\nContrast:", val1, "vs", val2, "\n")
                    cat("Variables involved:", col, "\n")

                    # Check if contrast should be included based on filters
                    if (!should_include_contrast(c(col), val1, val2)) {
                        next
                    }

                    # Determine if we should use automatic reference levels
                    use_primary <- primary_ref %in% c(val1, val2)
                    use_secondary <- !use_primary && !is.null(secondary_ref) &&
                        secondary_ref %in% c(val1, val2)

                    if (interactive) {
                        if (use_primary) {
                            cat("Primary reference level (", primary_ref, ") found\n")
                        } else if (use_secondary) {
                            cat("Secondary reference level (", secondary_ref, ") found\n")
                        }

                        cat("Include this contrast? (y/n): ")
                        if (tolower(readline()) != "y") next

                        if (use_primary) {
                            control <- primary_ref
                            experimental <- if (val1 == primary_ref) val2 else val1
                        } else if (use_secondary) {
                            control <- secondary_ref
                            experimental <- if (val1 == secondary_ref) val2 else val1
                        } else {
                            cat("Values: 1 =", val1, ", 2 =", val2, "\n")
                            cat("Choice (1=first as control, 2=second as control): ")
                            choice <- as.numeric(readline())

                            if (!choice %in% c(1, 2)) next

                            control <- if (choice == 1) val1 else val2
                            experimental <- if (choice == 1) val2 else val1
                        }
                    } else {
                        # Non-interactive mode
                        if (use_primary) {
                            control <- primary_ref
                            experimental <- if (val1 == primary_ref) val2 else val1
                        } else if (use_secondary) {
                            control <- secondary_ref
                            experimental <- if (val1 == secondary_ref) val2 else val1
                        } else {
                            control <- val1
                            experimental <- val2
                        }
                    }

                    contrast_name <- paste0(col, "_", experimental, "_vs_", control)
                    all_contrasts[[contrast_name]] <- list(
                        variable = col,
                        control = control,
                        experimental = experimental
                    )
                }
            }
        }
    }

    # For combined contrasts
    if (length(comparison_cols) > 1 && interactive) {
        cat("\nWould you like to generate contrasts holding other variables constant? (y/n): ")
        if (tolower(readline()) == "y") {
            # Create a temporary combined column
            temp_col <- paste(comparison_cols, collapse = "_")
            metadata_df[[temp_col]] <- apply(metadata_df[comparison_cols], 1, paste, collapse = "_")

            # Get actual existing combinations
            existing_combinations <- unique(metadata_df[[temp_col]])

            cat("\nProcessing combined contrasts for:", paste(comparison_cols, collapse = " + "), "\n")
            cat("Available combinations:", paste(existing_combinations, collapse = ", "), "\n")

            # Generate contrasts for existing combinations
            if (length(existing_combinations) > 1) {
                for (i in 1:(length(existing_combinations) - 1)) {
                    for (j in (i + 1):length(existing_combinations)) {
                        val1 <- existing_combinations[i]
                        val2 <- existing_combinations[j]

                        cat("\nCombined contrast:", val1, "vs", val2, "\n")
                        cat("Variables involved:", paste(comparison_cols, collapse = ", "), "\n")

                        # Check if combined contrast should be included
                        if (!should_include_contrast(comparison_cols, val1, val2)) {
                            next
                        }

                        # For combined contrasts, check reference levels in components
                        val1_parts <- strsplit(val1, "_")[[1]]
                        val2_parts <- strsplit(val2, "_")[[1]]
                        use_as_control <- NULL

                        # Check if any part matches primary reference levels
                        for (col in names(reference_levels)) {
                            if (!is.null(reference_levels[[col]])) {
                                primary_ref <- reference_levels[[col]][1]
                                if (primary_ref %in% val1_parts) {
                                    use_as_control <- val1
                                    break
                                } else if (primary_ref %in% val2_parts) {
                                    use_as_control <- val2
                                    break
                                }
                            }
                        }

                        if (interactive) {
                            if (!is.null(use_as_control)) {
                                cat("Reference level found in:", use_as_control, "\n")
                            }

                            cat("Include this contrast? (y/n): ")
                            if (tolower(readline()) != "y") next

                            if (is.null(use_as_control)) {
                                cat("Values: 1 =", val1, ", 2 =", val2, "\n")
                                cat("Choice (1=first as control, 2=second as control): ")
                                choice <- as.numeric(readline())

                                if (!choice %in% c(1, 2)) next

                                control <- if (choice == 1) val1 else val2
                                experimental <- if (choice == 1) val2 else val1
                            } else {
                                control <- use_as_control
                                experimental <- if (val1 == control) val2 else val1
                            }
                        } else {
                            # Non-interactive mode
                            if (!is.null(use_as_control)) {
                                control <- use_as_control
                                experimental <- if (val1 == control) val2 else val1
                            } else {
                                control <- val1
                                experimental <- val2
                            }
                        }

                        contrast_name <- paste0(temp_col, "_", experimental, "_vs_", control)
                        all_contrasts[[contrast_name]] <- list(
                            variable = temp_col,
                            control = control,
                            experimental = experimental
                        )
                    }
                }
            }
        }
    }

    # Print summary of contrasts
    if (length(all_contrasts) > 0) {
        message("\nGenerated ", length(all_contrasts), " contrasts:")
        for (name in names(all_contrasts)) {
            message("- ", name)
            message("  Control: ", all_contrasts[[name]]$control)
            message("  Experimental: ", all_contrasts[[name]]$experimental)
        }
    } else {
        message("\nNo contrasts were generated.")
    }

    return(all_contrasts)
}

# Flexible DESeq2 analysis function
run_flexible_deseq2 <- function(counts_matrix, sample_info, comparison,
                                min_count = 3, min_sample_percent = 0.2) {
    # Prepare data for comparison
    prepared_data <- prepare_comparison_data(counts_matrix, sample_info, comparison)

    # Filter low count genes
    filtered_data <- filter_low_counts(prepared_data$counts,
        min_count = min_count,
        min_sample_percent = min_sample_percent
    )

    # Create DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = filtered_data$counts,
        colData = prepared_data$sample_info,
        design = as.formula(paste0("~", comparison$variable))
    )

    # Run DESeq2
    dds <- DESeq2::DESeq(dds)

    # Get results with the exact parameters 
    res <- DESeq2::results(dds,
        contrast = c(
            comparison$variable,
            comparison$experimental,
            comparison$control
        ),
        pAdjustMethod = "fdr",
        cooksCutoff = FALSE,
        independentFiltering = FALSE,
        alpha = 0.05
    )

    # Add debugging
    cat("With alpha=0.05:\n")
    cat("Total results:", nrow(res), "\n")
    cat("NA p-values:", sum(is.na(res$pvalue)), "\n")
    cat("NA padj values:", sum(is.na(res$padj)), "\n")
    cat("Significant genes:", sum(res$padj < 0.05, na.rm = TRUE), "\n")

    return(list(
        dds = dds,
        deseq_results = res,
        comparison_info = comparison,
        filter_stats = filtered_data$stats
    ))
}

# Main DESeq2 analysis function
run_deseq2_analysis <- function(experiment_name,
                                comparison_cols = NULL,
                                reference_levels = NULL,
                                experimental_levels = NULL,
                                contrast_filters = NULL,
                                interactive = TRUE,
                                batch_mode = FALSE,
                                contrasts_file = NULL,
                                min_count = 3,
                                min_sample_percent = 0.2) {
    # Load experiment data with updated path
    experiment_obj <- readRDS(file.path(
        "data", experiment_name,
        paste0(experiment_name, ".rds")
    ))

    # Load cleaned data with updated path
    cleaned_data <- readRDS(file.path(
        "data", experiment_name, "cleaned_data",
        paste0(experiment_name, "_cleaned.rds")
    ))

    # Generate or load contrasts
    if (batch_mode && !is.null(contrasts_file)) {
        message("Loading contrasts from file:", contrasts_file)
        all_contrasts <- readRDS(contrasts_file)
    } else {
        all_contrasts <- generate_all_contrasts(
            metadata_df = cleaned_data$sample_info,
            comparison_cols = comparison_cols,
            reference_levels = reference_levels,
            experimental_levels = experimental_levels,
            contrast_filters = contrast_filters,
            interactive = interactive
        )

        # Save contrasts with updated path if file specified
        if (!is.null(contrasts_file)) {
            # Create deseq2/simple directory if it doesn't exist
            simple_dir <- file.path("results", experiment_name, "de", "deseq2", "simple")
            dir.create(simple_dir, recursive = TRUE, showWarnings = FALSE)

            # Update contrasts file path
            contrasts_file <- file.path(simple_dir, "contrasts.rds")
            message("Saving contrasts to file:", contrasts_file)
            saveRDS(all_contrasts, contrasts_file)
        }
    }

    # Run DESeq2 for each contrast
    all_results <- list()
    for (contrast_name in names(all_contrasts)) {
        message("\nProcessing contrast:", contrast_name)
        comparison <- all_contrasts[[contrast_name]]

        # Run DESeq2 analysis
        result <- run_flexible_deseq2(
            counts_matrix = cleaned_data$counts,
            sample_info = cleaned_data$sample_info,
            comparison = comparison,
            min_count = min_count,
            min_sample_percent = min_sample_percent
        )

        all_results[[contrast_name]] <- result
    }

    # Create deseq2/simple directory if it doesn't exist
    deseq2_dir <- file.path("results", experiment_name, "de", "deseq2", "simple")
    dir.create(deseq2_dir, recursive = TRUE, showWarnings = FALSE)

    # Save results to RDS file with updated path
    results_file <- file.path(deseq2_dir, paste0(experiment_name, "_de_results.rds"))
    saveRDS(all_results, results_file)
    message("\nResults saved to:", results_file)

    # Print and save analysis summary
    message("\n=== Differential Expression Analysis Complete ===")
    message(sprintf("Analyzed %d contrasts for experiment: %s", length(all_results), experiment_name))
    message("\nContrasts analyzed:")
    for (contrast_name in names(all_results)) {
        result <- all_results[[contrast_name]]
        n_sig <- sum(!is.na(result$deseq_results$padj) & result$deseq_results$padj < 0.05)
        message(sprintf("\n- %s", contrast_name))
        message(sprintf("  Control: %s", result$comparison_info$control))
        message(sprintf("  Experimental: %s", result$comparison_info$experimental))
        message(sprintf("  Significant genes (padj < 0.05): %d", n_sig))
    }

    # Save summary to Word document
    summary_file <- save_analysis_summary(all_results, experiment_name)
    message("\nSummary saved to:", summary_file)
    message("\nResults directory:", deseq2_dir)
    message("===================================================\n")

    return(all_results)
}

xen_tran_2024_03_de_results <- run_deseq2_analysis(
    experiment_name = "xen_tran_2024_03",
    comparison_cols = c("drug_treatment", "sample_prep"),
    reference_levels = list(
        drug_treatment = c("vehicle"),
        sample_prep = c("snap_frozen")
    ),
    experimental_levels = list(
        drug_treatment = c("wb3")
    ),
    interactive = TRUE,
    contrasts_file = file.path("results", "xen_tran_2024_03", "de", "deseq2", "simple", "contrasts.rds")
)

# Later - batch mode using saved selections
xen_tran_2024_03_de_results <- run_deseq2_analysis(
    experiment_name = "xen_tran_2024_03",
    batch_mode = TRUE,
    contrasts_file = file.path("results", "xen_tran_2024_03", "de", "deseq2", "simple", "contrasts.rds")
)

xen_tran_2024_03_de_results <- readRDS(file.path("results", "xen_tran_2024_03", "de", "deseq2", "simple", "xen_tran_2024_03_de_results.rds"))

# Check what contrasts are available
names(xen_tran_2024_03_de_results)
head(str(xen_tran_2024_03_de_results))
str(xen_tran_2024_03_de_results[1])
str(xen_tran_2024_03_de_results)

# Example usage:
xen_tran_2024_12_de_results <- run_deseq2_analysis(
    experiment_name = "xen_tran_2024_12",
    comparison_cols = c("condition", "timepoint"),
    reference_levels = list(
        condition = c("vehicle"),
        timepoint = c("r", "t")
    ),
    # experimental_levels = list(
    #    condition = "^wb3"
    # ),
    contrast_filters = list(
        valid_combinations = list(
            c("condition", "timepoint") # This should PREVENT single-variable contrasts from being shown
        ),
        min_samples = 3
    ),
    interactive = TRUE,
    contrasts_file = file.path("results", "xen_tran_2024_12", "de", "deseq2", "simple", "contrasts.rds")
)

# Later - batch mode using saved selections
xen_tran_2024_12_de_results <- run_deseq2_analysis(
    experiment_name = "xen_tran_2024_12",
    batch_mode = TRUE,
    contrasts_file = file.path("results", "xen_tran_2024_12", "de", "deseq2", "simple", "contrasts.rds")
)

xen_tran_2024_12_de_results <- readRDS(file.path("results", "xen_tran_2024_12", "de", "deseq2", "simple", "xen_tran_2024_12_de_results.rds"))
names(xen_tran_2024_12_de_results)
str(xen_tran_2024_12_de_results$dds)


######ME/CFS Project

# Run DESeq2 analysis for muscle dataset
hum_tran_2024_03_deep_phenotype_muscle_de_results <- run_deseq2_analysis(
    experiment_name = "hum_tran_2024_03_deep_phenotype_muscle",
    comparison_cols = c("group", "sex"),
    reference_levels = list(
        group = c("hv"),  # lowercase in the data
        sex = c("female")  # lowercase in the data
    ),
    experimental_levels = list(
        group = c("pi-me/cfs")
    ),
    interactive = TRUE,
    contrasts_file = file.path("results", "hum_tran_2024_03_deep_phenotype_muscle", "de", "deseq2", "simple", "contrasts.rds")
)

# Later - batch mode using saved selections for muscle dataset
hum_tran_2024_03_deep_phenotype_muscle_de_results <- run_deseq2_analysis(
    experiment_name = "hum_tran_2024_03_deep_phenotype_muscle",
    batch_mode = TRUE,
    contrasts_file = file.path("results", "hum_tran_2024_03_deep_phenotype_muscle", "de", "deseq2", "simple", "contrasts.rds")
)

hum_tran_2024_03_deep_phenotype_muscle_de_results <- readRDS(file.path("results", "hum_tran_2024_03_deep_phenotype_muscle", "de", "deseq2", "simple", "hum_tran_2024_03_deep_phenotype_muscle_de_results.rds"))
names(hum_tran_2024_03_deep_phenotype_muscle_de_results)
str(hum_tran_2024_03_deep_phenotype_muscle_de_results)

# Run DESeq2 analysis for muscle dataset
hum_tran_2024_03_deep_phenotype_pbmcs_de_results <- run_deseq2_analysis(
    experiment_name = "hum_tran_2024_03_deep_phenotype_pbmcs",
    comparison_cols = c("group", "sex"),
    reference_levels = list(
        group = c("hv"),  # lowercase in the data
        sex = c("female")  # lowercase in the data
    ),
    experimental_levels = list(
        group = c("pi-me/cfs")
    ),
    interactive = TRUE,
    contrasts_file = file.path("results", "hum_tran_2024_03_deep_phenotype_pbmcs", "de", "deseq2", "simple", "contrasts.rds")
)

# Later - batch mode using saved selections for muscle dataset
hum_tran_2024_03_deep_phenotype_pbmcs_de_results <- run_deseq2_analysis(
    experiment_name = "hum_tran_2024_03_deep_phenotype_pbmcs",
    batch_mode = TRUE,
    contrasts_file = file.path("results", "hum_tran_2024_03_deep_phenotype_pbmcs", "de", "deseq2", "simple", "contrasts.rds")
)

hum_tran_2024_03_deep_phenotype_pbmcs_de_results <- readRDS(file.path("results", "hum_tran_2024_03_deep_phenotype_pbmcs", "de", "deseq2", "simple", "hum_tran_2024_03_deep_phenotype_pbmcs.rds"))
names(hum_tran_2024_03_deep_phenotype_pbmcs_de_results)
str(hum_tran_2024_03_deep_phenotype_pbmcs_de_results)



# Run DESeq2 analysis for xen_tran_2025_02 dataset
xen_tran_2025_02_de_results <- run_deseq2_analysis(
    experiment_name = "xen_tran_2025_02",
    comparison_cols = c("condition"),  # Focus on condition comparisons
    # No reference levels filter - allow all contrasts
    contrast_filters = list(
        min_samples = 3
    ),
    interactive = TRUE,
    contrasts_file = file.path("results", "xen_tran_2025_02", "de", "deseq2", "simple", "contrasts.rds")
)

# Later - batch mode using saved selections
xen_tran_2025_02_de_results <- run_deseq2_analysis(
    experiment_name = "xen_tran_2025_02",
    batch_mode = TRUE,
    contrasts_file = file.path("results", "xen_tran_2025_02", "de", "deseq2", "simple", "contrasts.rds")
)

# Load and examine results
xen_tran_2025_02_de_results <- readRDS(file.path("results", "xen_tran_2025_02", "de", "deseq2", "simple", "xen_tran_2025_02_de_results.rds"))
names(xen_tran_2025_02_de_results)
str(xen_tran_2025_02_de_results)


