#' Simple DESeq2 Analysis Functions
#' 
#' This module provides functions for differential expression analysis using DESeq2
#' with comprehensive error handling, validation, and modular design.

# Required libraries
suppressPackageStartupMessages({
    library(DESeq2)    # For differential expression analysis
    library(dplyr)     # For data manipulation
    library(magrittr)  # For pipe operations
    library(flextable) # For creating formatted tables
    library(officer)   # For saving Word documents
    library(stringr)   # For string manipulation
    library(here)      # For path management
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
dea_config <- list(
    default_min_samples = 3,
    default_alpha = 0.05,
    default_lfc_threshold = 0,
    valid_analysis_types = c("simple", "complex")
)

#' Validate DEA parameters
#' @param experiment_name Name of experiment
#' @param comparison_cols Columns to use for comparisons
#' @param reference_levels Reference levels for factors
#' @param experimental_levels Experimental levels for factors
#' @param min_samples Minimum samples per group
validate_dea_params <- function(experiment_name, comparison_cols, reference_levels, 
                               experimental_levels, min_samples) {
    if (is.null(experiment_name) || nchar(experiment_name) == 0) {
        stop("experiment_name cannot be empty")
    }
    
    if (!is.character(comparison_cols) || length(comparison_cols) == 0) {
        stop("comparison_cols must be a non-empty character vector")
    }
    
    if (!is.numeric(min_samples) || min_samples < 2) {
        stop("min_samples must be a numeric value >= 2")
    }
    
    if (!is.null(reference_levels) && !is.list(reference_levels)) {
        stop("reference_levels must be a list or NULL")
    }
    
    if (!is.null(experimental_levels) && !is.list(experimental_levels)) {
        stop("experimental_levels must be a list or NULL")
    }
}

#' Load cleaned data for analysis
#' @param experiment_name Name of experiment
#' @return List containing counts matrix and sample info
load_cleaned_data <- function(experiment_name) {
    # Try to load from file first
    cleaned_data_path <- here("data", experiment_name, "cleaned_data", paste0(experiment_name, "_cleaned.rds"))
    
    if (!file.exists(cleaned_data_path)) {
        stop("Cleaned data file not found: ", cleaned_data_path)
    }
    
    tryCatch({
        cleaned_data <- readRDS(cleaned_data_path)
        
        if (!is.list(cleaned_data) || !"counts" %in% names(cleaned_data) || !"sample_info" %in% names(cleaned_data)) {
            stop("Invalid cleaned data structure")
        }
        
        return(cleaned_data)
    }, error = function(e) {
        stop("Failed to load cleaned data: ", e$message)
    })
}

#' Apply pattern-based filtering to sample data
#' @param sample_info Sample metadata data frame
#' @param filter_patterns List of filter patterns
#' @return Filtered sample info data frame
apply_pattern_filters <- function(sample_info, filter_patterns) {
    if (is.null(filter_patterns)) {
        return(sample_info)
    }
    
    keep_samples <- rep(TRUE, nrow(sample_info))
    
    for (pattern in filter_patterns) {
        if (pattern$action == "remove") {
            remove_idx <- grepl(pattern$pattern, sample_info$sample_id,
                              ignore.case = pattern$ignore_case)
            keep_samples <- keep_samples & !remove_idx
            message("Removing samples where ID matches: ", pattern$pattern)
        }
    }
    
    return(sample_info[keep_samples, ])
}

#' Prepare comparison data with validation
#' @param counts_matrix Counts matrix
#' @param sample_info Sample metadata
#' @param comparison Comparison configuration
#' @return List containing prepared data
prepare_comparison_data <- function(counts_matrix, sample_info, comparison) {
    message("Preparing comparison for variable: ", comparison$variable)
    
    # Validate inputs
    if (!is.matrix(counts_matrix) && !is.data.frame(counts_matrix)) {
        stop("counts_matrix must be a matrix or data frame")
    }
    
    if (!is.data.frame(sample_info)) {
        stop("sample_info must be a data frame")
    }
    
    # Initial droplevels
    sample_info <- droplevels(sample_info)
    
    # Apply pattern-based filtering if specified
    if (!is.null(comparison$filter_patterns)) {
        sample_info <- apply_pattern_filters(sample_info, comparison$filter_patterns)
        counts_matrix <- counts_matrix[, sample_info$sample_id]
        sample_info <- droplevels(sample_info)
    }
    
    # Clean any potential whitespace in relevant columns
    if (comparison$variable %in% names(sample_info)) {
        sample_info[[comparison$variable]] <- trimws(sample_info[[comparison$variable]])
    }
    
    # Handle composite variables
    if (grepl("_", comparison$variable)) {
        sample_info <- handle_composite_variable(sample_info, comparison$variable)
    }
    
    # Validate variable exists
    if (!comparison$variable %in% names(sample_info)) {
        stop("Variable '", comparison$variable, "' not found in sample_info")
    }
    
    # Check for sufficient samples
    if (nrow(sample_info) < 2) {
        stop("Insufficient samples after filtering")
    }
    
    return(list(
        counts = counts_matrix,
        sample_info = sample_info
    ))
}

#' Create metadata mappings for flexible matching
#' @param metadata_df Sample metadata data frame
#' @return List of mappings for each metadata column
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

#' Handle composite variables (e.g., drug_treatment_sample_prep)
#' @param sample_info Sample metadata
#' @param variable_name Variable name
#' @return Updated sample info with composite variable
handle_composite_variable <- function(sample_info, variable_name) {
    possible_cols <- names(sample_info)
    cols <- c()
    
    # Find matching columns
    for (col in possible_cols) {
        if (grepl(col, variable_name, ignore.case = TRUE)) {
            cols <- c(cols, col)
        }
    }
    
    if (length(cols) == 0) {
        stop("No matching columns found for composite variable: ", variable_name)
    }
    
    # Create composite variable
    sample_info[[variable_name]] <- do.call(paste, c(sample_info[cols], sep = "_"))
    
    return(sample_info)
}

#' Filter low-count genes with validation
#' @param counts_matrix Counts matrix
#' @param min_count Minimum count threshold
#' @param min_samples Minimum number of samples
#' @return Filtered counts matrix
filter_low_count_genes <- function(counts_matrix, min_count = 10, min_samples = 3) {
    if (!is.matrix(counts_matrix) && !is.data.frame(counts_matrix)) {
        stop("counts_matrix must be a matrix or data frame")
    }
    
    if (min_count < 0 || min_samples < 1) {
        stop("min_count must be >= 0 and min_samples must be >= 1")
    }
    
    # Count genes with sufficient counts
    keep_genes <- rowSums(counts_matrix >= min_count) >= min_samples
    
    if (sum(keep_genes) == 0) {
        stop("No genes pass the filtering criteria")
    }
    
    message("Keeping ", sum(keep_genes), " out of ", nrow(counts_matrix), " genes")
    
    return(counts_matrix[keep_genes, ])
}

#' Generate all possible contrasts with validation
#' @param sample_info Sample metadata
#' @param comparison_cols Columns to use for comparisons
#' @param reference_levels Reference levels
#' @param experimental_levels Experimental levels
#' @param min_samples Minimum samples per group
#' @param interactive Whether to run in interactive mode
#' @return List of selected contrasts
generate_contrasts <- function(sample_info, comparison_cols, reference_levels, 
                              experimental_levels, min_samples, interactive = TRUE) {
    contrasts <- list()
    
    for (col in comparison_cols) {
        if (!col %in% names(sample_info)) {
            warning("Column '", col, "' not found in sample_info, skipping")
            next
        }
        
        col_contrasts <- generate_column_contrasts(
            sample_info, col, reference_levels, experimental_levels, 
            min_samples, interactive
        )
        contrasts <- c(contrasts, col_contrasts)
    }
    
    return(contrasts)
}

#' Check if contrast should be included based on filters
#' @param variables_present Variables present in contrast
#' @param val1 First value
#' @param val2 Second value
#' @param contrast_filters Contrast filters
#' @param experimental_levels Experimental levels
#' @return Logical indicating if contrast should be included
should_include_contrast_filtered <- function(variables_present, val1, val2, contrast_filters, experimental_levels) {
    if (!is.null(contrast_filters) && !is.null(contrast_filters$valid_combinations)) {
        required_vars <- contrast_filters$valid_combinations[[1]]

        # For single-variable contrasts
        if (length(variables_present) == 1) {
            message("\n[Automatically skipped - requires both",
                   paste(required_vars, collapse = " and "), "]")
            return(FALSE)
        }

        # For combined contrasts, check if it has all required variables
        missing_vars <- setdiff(required_vars, variables_present)
        if (length(missing_vars) > 0) {
            message("\n[Automatically skipped - missing:",
                   paste(missing_vars, collapse = ", "), "]")
            return(FALSE)
        }
    }

    # Additional checks for experimental levels if provided
    if (!is.null(experimental_levels)) {
        for (var in names(experimental_levels)) {
            if (var %in% variables_present) {
                pattern <- experimental_levels[[var]]
                if (!any(grepl(pattern, c(val1, val2)))) {
                    message("\n[Automatically skipped - no matching experimental level for", var, "]")
                    return(FALSE)
                }
            }
        }
    }

    return(TRUE)
}

#' Clean string values for consistency
#' @param x String to clean
#' @return Cleaned string
clean_string <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("\\s+", "_", x) # replace spaces with underscores
    return(x)
}

#' Process single variable contrasts
#' @param metadata_df Sample metadata
#' @param col Column name
#' @param reference_levels Reference levels
#' @param interactive Whether to run interactively
#' @param contrast_filters Contrast filters
#' @param experimental_levels Experimental levels
#' @return List of contrasts for the column
process_single_variable_contrasts <- function(metadata_df, col, reference_levels, 
                                            interactive, contrast_filters, experimental_levels) {
    contrasts <- list()
    values <- unique(as.character(metadata_df[[col]]))
    message(sprintf("\nProcessing contrasts for column: %s", col))
    message("Available values:", paste(values, collapse = ", "))

    if (!is.null(reference_levels) && !is.null(reference_levels[[col]])) {
        # Get primary and secondary reference levels
        primary_ref <- reference_levels[[col]][1]
        secondary_ref <- if (length(reference_levels[[col]]) > 1) {
            reference_levels[[col]][2]
        } else {
            NULL
        }

        message(sprintf("\nPrimary reference level: %s", primary_ref))
        if (!is.null(secondary_ref)) {
            message(sprintf("Secondary reference level: %s", secondary_ref))
        }

        # Generate all possible pairs
        for (i in 1:(length(values) - 1)) {
            for (j in (i + 1):length(values)) {
                val1 <- values[i]
                val2 <- values[j]

                message("\nContrast:", val1, "vs", val2)
                message("Variables involved:", col)

                # Check if contrast should be included based on filters
                if (!should_include_contrast_filtered(c(col), val1, val2, contrast_filters, experimental_levels)) {
                    next
                }

                # Determine control and experimental levels
                control_exp <- determine_control_experimental_levels(
                    val1, val2, primary_ref, secondary_ref, interactive
                )
                
                if (is.null(control_exp)) next

                contrast_name <- paste0(col, "_", control_exp$experimental, "_vs_", control_exp$control)
                contrasts[[contrast_name]] <- list(
                    variable = col,
                    control = control_exp$control,
                    experimental = control_exp$experimental
                )
            }
        }
    }
    
    return(contrasts)
}

#' Determine control and experimental levels
#' @param val1 First value
#' @param val2 Second value
#' @param primary_ref Primary reference level
#' @param secondary_ref Secondary reference level
#' @param interactive Whether to run interactively
#' @return List with control and experimental levels
determine_control_experimental_levels <- function(val1, val2, primary_ref, secondary_ref, interactive) {
    # Determine if we should use automatic reference levels
    use_primary <- primary_ref %in% c(val1, val2)
    use_secondary <- !use_primary && !is.null(secondary_ref) &&
        secondary_ref %in% c(val1, val2)

    if (interactive) {
        if (use_primary) {
            message("Primary reference level (", primary_ref, ") found")
        } else if (use_secondary) {
            message("Secondary reference level (", secondary_ref, ") found")
        }

        message("Include this contrast? (y/n): ")
        if (tolower(readline()) != "y") return(NULL)

        if (use_primary) {
            control <- primary_ref
            experimental <- if (val1 == primary_ref) val2 else val1
        } else if (use_secondary) {
            control <- secondary_ref
            experimental <- if (val1 == secondary_ref) val2 else val1
        } else {
            message("Values: 1 =", val1, ", 2 =", val2)
            message("Choice (1=first as control, 2=second as control): ")
            choice <- as.numeric(readline())

            if (!choice %in% c(1, 2)) return(NULL)

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

    return(list(control = control, experimental = experimental))
}

#' Process combined variable contrasts
#' @param metadata_df Sample metadata
#' @param comparison_cols Columns to compare
#' @param reference_levels Reference levels
#' @param interactive Whether to run interactively
#' @param contrast_filters Contrast filters
#' @param experimental_levels Experimental levels
#' @return List of combined contrasts
process_combined_variable_contrasts <- function(metadata_df, comparison_cols, reference_levels,
                                              interactive, contrast_filters, experimental_levels) {
    contrasts <- list()
    
    if (length(comparison_cols) <= 1 || !interactive) {
        return(contrasts)
    }
    
    message("\nWould you like to generate contrasts holding other variables constant? (y/n): ")
    if (tolower(readline()) != "y") {
        return(contrasts)
    }
    
    # Create a temporary combined column
    temp_col <- paste(comparison_cols, collapse = "_")
    metadata_df[[temp_col]] <- apply(metadata_df[comparison_cols], 1, paste, collapse = "_")

    # Get actual existing combinations
    existing_combinations <- unique(metadata_df[[temp_col]])

    message("\nProcessing combined contrasts for:", paste(comparison_cols, collapse = " + "))
    message("Available combinations:", paste(existing_combinations, collapse = ", "))

    # Generate contrasts for existing combinations
    if (length(existing_combinations) > 1) {
        for (i in 1:(length(existing_combinations) - 1)) {
            for (j in (i + 1):length(existing_combinations)) {
                val1 <- existing_combinations[i]
                val2 <- existing_combinations[j]

                message("\nCombined contrast:", val1, "vs", val2)
                message("Variables involved:", paste(comparison_cols, collapse = ", "))

                # Check if combined contrast should be included
                if (!should_include_contrast_filtered(comparison_cols, val1, val2, contrast_filters, experimental_levels)) {
                    next
                }

                # Determine control and experimental for combined contrasts
                control_exp <- determine_combined_control_experimental(
                    val1, val2, reference_levels, interactive
                )
                
                if (is.null(control_exp)) next

                contrast_name <- paste0(temp_col, "_", control_exp$experimental, "_vs_", control_exp$control)
                contrasts[[contrast_name]] <- list(
                    variable = temp_col,
                    control = control_exp$control,
                    experimental = control_exp$experimental
                )
            }
        }
    }
    
    return(contrasts)
}

#' Determine control and experimental for combined contrasts
#' @param val1 First value
#' @param val2 Second value
#' @param reference_levels Reference levels
#' @param interactive Whether to run interactively
#' @return List with control and experimental levels
determine_combined_control_experimental <- function(val1, val2, reference_levels, interactive) {
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
            message("Reference level found in:", use_as_control)
        }

        message("Include this contrast? (y/n): ")
        if (tolower(readline()) != "y") return(NULL)

        if (is.null(use_as_control)) {
            message("Values: 1 =", val1, ", 2 =", val2)
            message("Choice (1=first as control, 2=second as control): ")
            choice <- as.numeric(readline())

            if (!choice %in% c(1, 2)) return(NULL)

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

    return(list(control = control, experimental = experimental))
}

#' Generate all possible contrasts with advanced interactive selection
#' @param metadata_df Sample metadata data frame
#' @param comparison_cols Columns to use for comparisons
#' @param reference_levels Reference levels for factors
#' @param experimental_levels Experimental levels for factors
#' @param contrast_filters Additional contrast filters
#' @param interactive Whether to run in interactive mode
#' @return List of selected contrasts
generate_all_contrasts <- function(metadata_df, comparison_cols = NULL,
                                   reference_levels = NULL,
                                   experimental_levels = NULL,
                                   contrast_filters = NULL,
                                   interactive = TRUE) {
    message("Generating all possible contrasts...")
    
    # If no specific columns provided, use all columns except sample_id
    if (is.null(comparison_cols)) {
        comparison_cols <- names(metadata_df)[!names(metadata_df) %in% c("sample_id")]
    }
    
    message("Comparison columns: ", paste(comparison_cols, collapse = ", "))

    # Clean the metadata values
    for (col in comparison_cols) {
        metadata_df[[col]] <- clean_string(metadata_df[[col]])
    }

    # Clean reference levels if provided
    if (!is.null(reference_levels)) {
        reference_levels <- lapply(reference_levels, clean_string)
    }

    # Generate single variable contrasts
    all_contrasts <- list()
    for (col in comparison_cols) {
        single_contrasts <- process_single_variable_contrasts(
            metadata_df, col, reference_levels, interactive, 
            contrast_filters, experimental_levels
        )
        all_contrasts <- c(all_contrasts, single_contrasts)
    }

    # Generate combined variable contrasts
    combined_contrasts <- process_combined_variable_contrasts(
        metadata_df, comparison_cols, reference_levels, interactive,
        contrast_filters, experimental_levels
    )
    all_contrasts <- c(all_contrasts, combined_contrasts)

    # Print summary of contrasts
    message("\n=== Contrast Generation Summary ===")
    message("Total contrasts generated: ", length(all_contrasts))
    
    if (length(all_contrasts) > 0) {
        message("\nGenerated contrasts:")
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

#' Generate contrasts for a single column
#' @param sample_info Sample metadata
#' @param col Column name
#' @param reference_levels Reference levels
#' @param experimental_levels Experimental levels
#' @param min_samples Minimum samples per group
#' @param interactive Whether to run in interactive mode
#' @return List of contrasts for the column
generate_column_contrasts <- function(sample_info, col, reference_levels, 
                                     experimental_levels, min_samples, interactive) {
    available_values <- unique(sample_info[[col]])
    available_values <- available_values[!is.na(available_values)]
    
    if (length(available_values) < 2) {
        message("Column '", col, "' has insufficient levels for comparison")
        return(list())
    }
    
    message("Processing contrasts for column: ", col)
    message("Available values: ", paste(available_values, collapse = ", "))
    
    # Determine reference level
    ref_level <- determine_reference_level(col, available_values, reference_levels)
    message("Primary reference level: ", ref_level)
    
    contrasts <- list()
    
    # Generate pairwise contrasts
    for (i in seq_along(available_values)) {
        for (j in seq_along(available_values)) {
            if (i != j) {
                contrast <- create_contrast(
                    sample_info, col, available_values[i], available_values[j],
                    ref_level, min_samples, interactive
                )
                if (!is.null(contrast)) {
                    contrasts <- c(contrasts, list(contrast))
                }
            }
        }
    }
    
    return(contrasts)
}

#' Determine reference level for a column
#' @param col Column name
#' @param available_values Available values
#' @param reference_levels Reference levels configuration
#' @return Reference level
determine_reference_level <- function(col, available_values, reference_levels) {
    if (!is.null(reference_levels) && col %in% names(reference_levels)) {
        ref_levels <- reference_levels[[col]]
        for (ref in ref_levels) {
            if (ref %in% available_values) {
                return(ref)
            }
        }
    }
    
    # Default to first value
    return(available_values[1])
}

#' Create a single contrast with validation
#' @param sample_info Sample metadata
#' @param col Column name
#' @param level1 First level
#' @param level2 Second level
#' @param ref_level Reference level
#' @param min_samples Minimum samples per group
#' @param interactive Whether to run in interactive mode
#' @return Contrast object or NULL
create_contrast <- function(sample_info, col, level1, level2, ref_level, 
                           min_samples, interactive) {
    # Check sample sizes
    n1 <- sum(sample_info[[col]] == level1, na.rm = TRUE)
    n2 <- sum(sample_info[[col]] == level2, na.rm = TRUE)
    
    if (n1 < min_samples || n2 < min_samples) {
        message("Contrast ", level1, " vs ", level2, " skipped - insufficient samples")
        return(NULL)
    }
    
    contrast_name <- paste(level1, "vs", level2)
    message("Contrast: ", contrast_name)
    message("Variables involved: ", col)
    
    if (ref_level %in% c(level1, level2)) {
        message("Primary reference level (", ref_level, ") found")
        
        if (interactive) {
            response <- readline("Include this contrast? (y/n): ")
            if (tolower(response) %in% c("y", "yes")) {
                return(list(
                    name = contrast_name,
                    variable = col,
                    level1 = level1,
                    level2 = level2,
                    reference = ref_level,
                    n1 = n1,
                    n2 = n2
                ))
            }
        } else {
            return(list(
                name = contrast_name,
                variable = col,
                level1 = level1,
                level2 = level2,
                reference = ref_level,
                n1 = n1,
                n2 = n2
            ))
        }
    }
    
    return(NULL)
}

#' Run DESeq2 analysis with comprehensive error handling
#' @param experiment_name Name of experiment
#' @param comparison_cols Columns to use for comparisons
#' @param reference_levels Reference levels for factors
#' @param experimental_levels Experimental levels for factors
#' @param min_samples Minimum samples per group
#' @param interactive Whether to run in interactive mode
#' @param batch_mode Whether to run in batch mode
#' @param contrasts_file Path to contrasts file
#' @param filter_patterns Patterns for filtering samples
#' @param contrast_filters Additional contrast filters
#' @return Analysis results
run_deseq2_analysis <- function(experiment_name,
                               comparison_cols = c("condition"),
                               reference_levels = NULL,
                               experimental_levels = NULL,
                               min_samples = dea_config$default_min_samples,
                               interactive = TRUE,
                               batch_mode = FALSE,
                               contrasts_file = NULL,
                               filter_patterns = NULL,
                               contrast_filters = NULL) {
    
    # Validate parameters
    validate_dea_params(experiment_name, comparison_cols, reference_levels, 
                       experimental_levels, min_samples)
    
    # Load cleaned data
    cleaned_data <- load_cleaned_data(experiment_name)
    counts_matrix <- cleaned_data$counts
    sample_info <- cleaned_data$sample_info
    
    # Create results directory
    results_dir <- here("results", experiment_name, "de", "deseq2", "simple")
    if (!dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Failed to create results directory: ", results_dir)
    }
    
    # Handle batch mode
    if (batch_mode) {
        if (is.null(contrasts_file)) {
            stop("contrasts_file is required for batch mode")
        }
        
        if (!file.exists(contrasts_file)) {
            stop("Contrasts file not found: ", contrasts_file)
        }
        
        contrasts <- readRDS(contrasts_file)
        message("Loading contrasts from file: ", contrasts_file)
    } else {
        # Generate contrasts using advanced function
        contrasts <- generate_all_contrasts(
            metadata_df = sample_info,
            comparison_cols = comparison_cols,
            reference_levels = reference_levels,
            experimental_levels = experimental_levels,
            contrast_filters = contrast_filters,
            interactive = interactive
        )
        
        # Save contrasts
        if (!is.null(contrasts_file)) {
            saveRDS(contrasts, contrasts_file)
            message("Saving contrasts to file: ", contrasts_file)
        }
    }
    
    if (length(contrasts) == 0) {
        warning("No contrasts were generated")
        return(list())
    }
    
    # Run DESeq2 for each contrast using flexible function
    results <- list()
    for (contrast in contrasts) {
        tryCatch({
            # Prepare data for this contrast
            comparison_data <- prepare_comparison_data(
                counts_matrix, sample_info, 
                list(variable = contrast$variable, filter_patterns = filter_patterns)
            )
            
            # Create DESeqDataSet
            dds <- DESeqDataSetFromMatrix(
                countData = comparison_data$counts,
                colData = comparison_data$metadata,
                design = comparison_data$design
            )
            
            # Run flexible DESeq2 analysis
            result <- run_flexible_deseq2(
                dds = dds,
                contrast = contrast,
                gene_symbols = NULL,  # Will be extracted from dds if available
                filter_low_counts = TRUE,
                min_count = 10,
                min_samples = 3,
                alpha = 0.05,
                lfc_threshold = 0,
                independent_filtering = TRUE,
                cooks_cutoff = TRUE,
                parallel = FALSE,
                BPPARAM = NULL
            )
            
            if (!is.null(result)) {
                results[[contrast$name]] <- result
            }
        }, error = function(e) {
            warning("Failed to analyze contrast '", contrast$name, "': ", e$message)
        })
    }
    
    # Save results
    results_file <- file.path(results_dir, paste0(experiment_name, "_de_results.rds"))
    saveRDS(results, results_file)
    message("Results saved to: ", results_file)
    
    # Generate summary
    save_analysis_summary(results, experiment_name, results_dir)
    
    return(results)
}

#' Run DESeq2 analysis for a single contrast
#' @param counts_matrix Counts matrix
#' @param sample_info Sample metadata
#' @param contrast Contrast configuration
#' @param filter_patterns Filter patterns
#' @param contrast_filters Contrast filters
#' @return DESeq2 results or NULL
run_single_contrast <- function(counts_matrix, sample_info, contrast, 
                               filter_patterns, contrast_filters) {
    # Prepare data for this contrast
    comparison_data <- prepare_comparison_data(
        counts_matrix, sample_info, 
        list(variable = contrast$variable, filter_patterns = filter_patterns)
    )
    
    # Filter low-count genes
    filtered_counts <- filter_low_count_genes(
        comparison_data$counts, 
        min_count = 10, 
        min_samples = 3
    )
    
    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(
        countData = filtered_counts,
        colData = comparison_data$sample_info,
        design = as.formula(paste("~", contrast$variable))
    )
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds, contrast = c(contrast$variable, contrast$level1, contrast$level2))
    
    return(list(
        contrast = contrast,
        results = res,
        dds = dds
    ))
}

#' Create DESeq2 dataset from prepared data
#' @param filtered_counts Filtered counts matrix
#' @param sample_info Sample metadata
#' @param comparison Comparison configuration
#' @return DESeq2 dataset
create_deseq2_dataset <- function(filtered_counts, sample_info, comparison) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = filtered_counts,
        colData = sample_info,
        design = as.formula(paste0("~", comparison$variable))
    )
    return(dds)
}

#' Run DESeq2 analysis on dataset
#' @param dds DESeq2 dataset
#' @return Fitted DESeq2 dataset
run_deseq2_analysis <- function(dds) {
    message("Running DESeq2 analysis...")
    dds <- DESeq2::DESeq(dds)
    return(dds)
}

#' Extract DESeq2 results with specific parameters
#' @param dds Fitted DESeq2 dataset
#' @param comparison Comparison configuration
#' @param alpha Significance threshold
#' @return DESeq2 results
extract_deseq2_results <- function(dds, comparison, alpha = 0.05) {
    message("Extracting DESeq2 results...")
    
    res <- DESeq2::results(dds,
        contrast = c(
            comparison$variable,
            comparison$experimental,
            comparison$control
        ),
        pAdjustMethod = "fdr",
        cooksCutoff = FALSE,
        independentFiltering = FALSE,
        alpha = alpha
    )
    
    # Add debugging information
    message("Results summary:")
    message("- Total results: ", nrow(res))
    message("- NA p-values: ", sum(is.na(res$pvalue)))
    message("- NA padj values: ", sum(is.na(res$padj)))
    message("- Significant genes (padj < ", alpha, "): ", sum(res$padj < alpha, na.rm = TRUE))
    
    return(res)
}

#' Flexible DESeq2 analysis function
#' @param counts_matrix Counts matrix
#' @param sample_info Sample metadata
#' @param comparison Comparison configuration
#' @param min_count Minimum count threshold
#' @param min_sample_percent Minimum sample percentage
#' @return DESeq2 analysis results
run_flexible_deseq2 <- function(counts_matrix, sample_info, comparison,
                                min_count = 3, min_sample_percent = 0.2) {
    message("Starting flexible DESeq2 analysis...")
    message("Comparison: ", comparison$experimental, " vs ", comparison$control)
    
    # Prepare data for comparison
    prepared_data <- prepare_comparison_data(counts_matrix, sample_info, comparison)

    # Filter low count genes
    filtered_data <- filter_low_counts(prepared_data$counts,
        min_count = min_count,
        min_sample_percent = min_sample_percent
    )

    # Create DESeq2 object
    dds <- create_deseq2_dataset(filtered_data$counts, prepared_data$sample_info, comparison)

    # Run DESeq2
    dds <- run_deseq2_analysis(dds)

    # Get results
    res <- extract_deseq2_results(dds, comparison)

    message("DESeq2 analysis completed successfully")
    
    return(list(
        dds = dds,
        deseq_results = res,
        comparison_info = comparison,
        filter_stats = filtered_data$stats
    ))
}

#' Filter genes based on count threshold
#' @param counts_matrix Counts matrix
#' @param min_count Minimum count threshold
#' @param min_sample_percent Minimum sample percentage
#' @return Filtered counts and statistics
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

#' Generate analysis summary
#' @param results Analysis results
#' @param experiment_name Experiment name
#' @param results_dir Results directory
generate_analysis_summary <- function(results, experiment_name, results_dir) {
    message("\n=== Differential Expression Analysis Complete ===")
    message("Analyzed ", length(results), " contrasts for experiment: ", experiment_name)
    message("\nContrasts analyzed:")
    
    if (length(results) > 0) {
        for (name in names(results)) {
            message("- ", name)
        }
    } else {
        message("None")
    }
    
    # Create summary document
    summary_file <- file.path(results_dir, paste0(experiment_name, "_de_summary.docx"))
    create_summary_document(results, experiment_name, summary_file)
    
    message("\nSummary saved to: ", summary_file)
    message("Results directory: ", results_dir)
    message("===================================================")
}

#' Create summary document
#' @param results Analysis results
#' @param experiment_name Experiment name
#' @param output_file Output file path
create_summary_document <- function(results, experiment_name, output_file) {
    tryCatch({
        doc <- read_docx()
        doc <- doc %>%
            body_add_par(paste("Differential Expression Analysis Summary:", experiment_name), 
                        style = "heading 1") %>%
            body_add_par(paste("Analysis completed on:", Sys.time()), 
                        style = "Normal") %>%
            body_add_par(paste("Number of contrasts analyzed:", length(results)), 
                        style = "Normal")
        
        if (length(results) > 0) {
            doc <- doc %>% body_add_par("Contrasts analyzed:", style = "heading 2")
            for (name in names(results)) {
                doc <- doc %>% body_add_par(paste("-", name), style = "Normal")
            }
        }
        
        print(doc, target = output_file)
    }, error = function(e) {
        warning("Failed to create summary document: ", e$message)
    })
}

#' Create summary text for analysis results
#' @param all_results Analysis results
#' @param experiment_name Experiment name
#' @return Summary text vector
create_summary_text <- function(all_results, experiment_name) {
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
    
    return(summary_text)
}

#' Create formatted table from summary text
#' @param summary_text Summary text vector
#' @return Formatted flextable
create_summary_table <- function(summary_text) {
    ft <- flextable::flextable(data.frame(Summary = summary_text)) %>%
        flextable::delete_part(part = "header") %>%
        flextable::border_remove() %>%
        flextable::padding(padding = 0) %>%
        flextable::autofit()
    
    return(ft)
}

#' Save summary table to Word document
#' @param ft Formatted flextable
#' @param experiment_name Experiment name
#' @return Path to summary file
save_summary_document <- function(ft, experiment_name) {
    # Create simple directory if it doesn't exist
    simple_dir <- file.path("results", experiment_name, "de", "deseq2", "simple")
    dir.create(simple_dir, recursive = TRUE, showWarnings = FALSE)

    # Save to docx with updated path
    summary_file <- file.path(simple_dir, paste0(experiment_name, "_de_summary.docx"))
    flextable::save_as_docx(ft, path = summary_file)
    
    message("Analysis summary saved to: ", summary_file)
    return(summary_file)
}

#' Save analysis summary to Word document
#' @param all_results All analysis results
#' @param experiment_name Experiment name
#' @return Path to summary file
save_analysis_summary <- function(all_results, experiment_name) {
    message("Creating analysis summary...")
    
    # Create summary text
    summary_text <- create_summary_text(all_results, experiment_name)
    
    # Create formatted table
    ft <- create_summary_table(summary_text)
    
    # Save to document
    summary_file <- save_summary_document(ft, experiment_name)
    
    return(summary_file)
}

message("DEA functions loaded successfully. Use run_deseq2_analysis() to perform differential expression analysis.")
