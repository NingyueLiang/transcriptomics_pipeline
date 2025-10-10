#' Data Cleaning Functions for Drug Explorer Pipeline
#' 
#' This module provides functions for cleaning and preparing omics data
#' for differential expression analysis with proper validation and error handling.

# Required libraries
suppressPackageStartupMessages({
    library(dplyr)    # For data manipulation
    library(tidyr)    # For data tidying
    library(magrittr) # For pipe operations
    library(tibble)   # For tibble operations
    library(stringr)  # For string manipulation
    library(here)     # For path management
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
cleaning_config <- list(
    default_metadata_cols = list(
        sample_id = "sample_id",
        group_cols = c("sample_prep", "drug_treatment")
    ),
    default_factor_levels = list(),
    default_numeric_cols = character(0)
)

#' Clean string values and handle NAs properly
#' @param x Vector of strings to clean
#' @return Cleaned vector
clean_string <- function(x) {
    if (!is.character(x) && !is.factor(x)) {
        stop("Input must be character or factor vector")
    }
    
    # Convert to character if factor
    x <- as.character(x)
    
    # Convert "NA" strings to actual NA
    na_patterns <- c("NA", "na", "N/A", "n/a", "", " ")
    x[x %in% na_patterns] <- NA
    
    # Only process non-NA values
    non_na <- !is.na(x)
    if (any(non_na)) {
        x[non_na] <- trimws(x[non_na]) # Remove leading/trailing spaces
        x[non_na] <- tolower(x[non_na]) # Convert to lowercase
        x[non_na] <- gsub("\\s+", "_", x[non_na]) # Replace spaces with underscores
    }
    
    return(x)
}

#' Clean numeric columns with validation
#' @param x Vector to convert to numeric
#' @return Numeric vector
clean_numeric <- function(x) {
    if (is.numeric(x)) {
        return(x)
    }
    
    # Convert to numeric, handling NA strings automatically
    result <- suppressWarnings(as.numeric(as.character(x)))
    
    # Check for conversion issues
    if (any(is.na(result) & !is.na(x))) {
        warning("Some values could not be converted to numeric and were set to NA")
    }
    
    return(result)
}

#' Convert assay data to matrix format
#' @param assay_data Data frame containing assay data
#' @return Matrix with gene symbols as rownames
convert_assay_to_matrix <- function(assay_data) {
    if (!"gene_symbol" %in% colnames(assay_data)) {
        stop("Assay data must contain a 'gene_symbol' column")
    }
    
    # Check for duplicate gene symbols
    if (any(duplicated(assay_data$gene_symbol))) {
        warning("Duplicate gene symbols found in assay data. Using first occurrence.")
        assay_data <- assay_data[!duplicated(assay_data$gene_symbol), ]
    }
    
    rownames(assay_data) <- assay_data$gene_symbol
    assay_matrix <- assay_data[, -which(colnames(assay_data) == "gene_symbol")]
    assay_matrix <- as.matrix(assay_matrix)
    
    return(assay_matrix)
}

#' Clean and format sample metadata
#' @param sample_data Data frame containing sample metadata
#' @param metadata_cols Metadata column configuration
#' @param factor_levels Factor level configuration
#' @param numeric_cols Numeric column names
#' @return Cleaned sample data frame
clean_sample_metadata <- function(sample_data, metadata_cols, factor_levels, numeric_cols) {
    # Check for required sample_id column
    if (!metadata_cols$sample_id %in% colnames(sample_data)) {
        stop("Sample data missing required column: ", metadata_cols$sample_id)
    }
    
    # Process each grouping column
    for (col in metadata_cols$group_cols) {
        if (!col %in% colnames(sample_data)) {
            stop("Sample data missing group column: ", col)
        }
        
        if (col %in% numeric_cols) {
            # Handle numeric columns
            sample_data[[col]] <- clean_numeric(sample_data[[col]])
        } else {
            # Clean the column values
            sample_data[[col]] <- clean_string(sample_data[[col]])
            
            # Set factor levels if provided
            if (!is.null(factor_levels[[col]])) {
                # Clean the factor levels too
                factor_levels[[col]] <- clean_string(factor_levels[[col]])
                sample_data[[col]] <- factor(sample_data[[col]],
                                           levels = factor_levels[[col]])
            } else {
                sample_data[[col]] <- factor(sample_data[[col]])
            }
        }
    }
    
    return(sample_data)
}

#' Validate cleaning parameters
#' @param experiment_name Name of experiment
#' @param experiment_obj Experiment object
#' @param metadata_cols Metadata column configuration
#' @param factor_levels Factor level configuration
#' @param numeric_cols Numeric column names
validate_cleaning_params <- function(experiment_name, experiment_obj, metadata_cols, 
                                   factor_levels, numeric_cols) {
    if (is.null(experiment_obj) && is.null(experiment_name)) {
        stop("Must provide either experiment_name or experiment_obj")
    }
    
    if (!is.list(metadata_cols) || !"sample_id" %in% names(metadata_cols)) {
        stop("metadata_cols must be a list containing 'sample_id'")
    }
    
    if (!is.list(factor_levels)) {
        stop("factor_levels must be a list")
    }
    
    if (!is.character(numeric_cols)) {
        stop("numeric_cols must be a character vector")
    }
}

#' Load experiment object from file or environment
#' @param experiment_name Name of experiment
#' @param experiment_obj Experiment object (if already loaded)
#' @return Experiment object
load_experiment_object <- function(experiment_name, experiment_obj) {
    if (!is.null(experiment_obj)) {
        return(experiment_obj)
    }
    
    # Try to load from file first
    experiment_path <- here("data", experiment_name, paste0(experiment_name, ".rds"))
    if (file.exists(experiment_path)) {
        tryCatch({
            return(readRDS(experiment_path))
        }, error = function(e) {
            stop("Failed to load experiment from file: ", e$message)
        })
    }
    
    # Fall back to global environment
    if (exists(experiment_name, envir = .GlobalEnv)) {
        return(get(experiment_name, envir = .GlobalEnv))
    }
    
    stop("Experiment object not found in file or global environment")
}

#' Create final cleaned data object with consistent ordering
#' @param assay_matrix Assay data matrix
#' @param sample_data Sample metadata
#' @param metadata_cols Metadata column configuration
#' @param original_assay_data Original assay data for dimensions
#' @param cleaned_data_dir Directory to save cleaned data
#' @param prefix Experiment prefix
#' @return Name of cleaned data object
create_cleaned_data_object <- function(assay_matrix, sample_data, metadata_cols, 
                                      original_assay_data, cleaned_data_dir, prefix) {
    # Ensure consistent ordering
    header_order <- sample_data[[metadata_cols$sample_id]]
    
    # Verify all sample IDs exist in assay data
    missing_samples <- setdiff(header_order, colnames(assay_matrix))
    if (length(missing_samples) > 0) {
        stop("Some sample IDs from sample_data are missing in assay_data: ",
             paste(missing_samples, collapse = ", "))
    }
    
    # Reorder the assay matrix
    assay_matrix_reordered <- assay_matrix[, as.character(header_order)]
    
    # Create cleaned data object with proper naming
    cleaned_data_name <- paste0(prefix, "_cleaned")
    
    # Create list with cleaned data
    cleaned_data <- list(
        counts = assay_matrix_reordered,
        sample_info = sample_data,
        cleaning_info = list(
            original_dims = list(
                assay = dim(original_assay_data),
                sample = dim(sample_data)
            ),
            cleaned_dims = list(
                counts = dim(assay_matrix_reordered),
                sample_info = dim(sample_data)
            ),
            metadata_cols = metadata_cols,
            factor_levels = lapply(
                metadata_cols$group_cols,
                function(col) levels(sample_data[[col]])
            )
        )
    )
    
    # Save cleaned data to file instead of global environment
    tryCatch({
        save_path <- file.path(cleaned_data_dir, paste0(cleaned_data_name, ".rds"))
        saveRDS(cleaned_data, save_path)
        message("Cleaned data saved to: ", save_path)
    }, error = function(e) {
        stop("Failed to save cleaned data: ", e$message)
    })
    
    return(cleaned_data_name)
}

#' Main function to clean and prepare data for differential expression analysis
#' @param experiment_name Name of experiment
#' @param experiment_obj Experiment object (optional)
#' @param metadata_cols Metadata column configuration
#' @param factor_levels Factor level configuration
#' @param numeric_cols Numeric column names
#' @return Name of cleaned data object
clean_data_for_de <- function(
    experiment_name = NULL,
    experiment_obj = NULL,
    metadata_cols = cleaning_config$default_metadata_cols,
    factor_levels = cleaning_config$default_factor_levels,
    numeric_cols = cleaning_config$default_numeric_cols
) {
    # Validate parameters
    validate_cleaning_params(experiment_name, experiment_obj, metadata_cols, 
                           factor_levels, numeric_cols)
    
    # Load experiment object
    experiment_obj <- load_experiment_object(experiment_name, experiment_obj)

    # Get config from experiment object
    config <- experiment_obj$metadata$config
    if (is.null(config) || is.null(config$data_dir)) {
        stop("Invalid experiment object: missing config or data_dir")
    }
    
    # Create cleaned data directory
    cleaned_data_dir <- file.path(config$data_dir, "cleaned_data")
    tryCatch({
        if (!dir.create(cleaned_data_dir, recursive = TRUE, showWarnings = FALSE)) {
            stop("Failed to create cleaned data directory: ", cleaned_data_dir)
        }
    }, error = function(e) {
        stop("Error creating cleaned data directory: ", e$message)
    })

    # Get data using the experiment name as prefix
    prefix <- experiment_obj$metadata$prefix
    if (is.null(prefix)) {
        stop("Experiment object missing prefix in metadata")
    }
    
    sample_data_name <- paste0(prefix, "_sample_data")
    assay_data_name <- paste0(prefix, "_assay_data")

    # Extract data from experiment object with validation
    if (!sample_data_name %in% names(experiment_obj$sample_data)) {
        stop("Sample data not found: ", sample_data_name)
    }
    if (!assay_data_name %in% names(experiment_obj$assay_data)) {
        stop("Assay data not found: ", assay_data_name)
    }
    
    sample_data <- experiment_obj$sample_data[[sample_data_name]]
    assay_data <- experiment_obj$assay_data[[assay_data_name]]
    
    # Validate data
    if (is.null(sample_data) || nrow(sample_data) == 0) {
        stop("Sample data is empty or NULL")
    }
    if (is.null(assay_data) || nrow(assay_data) == 0) {
        stop("Assay data is empty or NULL")
    }

    # Step 1: Convert assay data to matrix
    assay_matrix <- convert_assay_to_matrix(assay_data)
    
    # Step 2: Clean and format sample metadata
    sample_data <- clean_sample_metadata(sample_data, metadata_cols, factor_levels, numeric_cols)

    # Step 3: Ensure consistent ordering and create final object
    cleaned_data_name <- create_cleaned_data_object(
        assay_matrix, sample_data, metadata_cols, 
        assay_data, cleaned_data_dir, prefix
    )
    
    message(sprintf("Created cleaned data object '%s' in %s", cleaned_data_name, cleaned_data_dir))
    
    return(cleaned_data_name)
}

message("Cleaning functions loaded successfully. Use clean_data_for_de() to clean experiment data.")
