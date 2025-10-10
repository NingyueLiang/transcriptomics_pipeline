#####
# Data Cleaning Functions

# Required libraries
library(dplyr) # For data manipulation
library(tidyr) # For data tidying
library(magrittr) # For pipe operations
library(tibble) # For tibble operations
library(stringr) # For string manipulation (used in clean_string function)

# Helper function to clean strings and handle NAs properly
clean_string <- function(x) {
    # Convert "NA" strings to actual NA
    x[x == "NA" | x == "na" | x == "N/A" | x == "n/a"] <- NA

    # Only process non-NA values
    non_na <- !is.na(x)
    if (any(non_na)) {
        x[non_na] <- trimws(x[non_na]) # Remove leading/trailing spaces
        x[non_na] <- tolower(x[non_na]) # Convert to lowercase
        x[non_na] <- gsub("\\s+", "_", x[non_na]) # Replace spaces with underscores
    }

    return(x)
}

# Helper function to clean numeric columns
clean_numeric <- function(x) {
    # Convert to numeric, handling NA strings automatically
    x <- as.numeric(as.character(x))
    return(x)
}

# Main function to clean and prepare data for differential expression analysis
clean_data_for_de <- function(
    experiment_name = NULL,
    experiment_obj = NULL,
    # Make column names configurable
    metadata_cols = list(
        sample_id = "sample_id", # Required
        group_cols = c("sample_prep", "drug_treatment") # Optional grouping columns
    ),
    # Make factor levels optional and flexible
    factor_levels = list(), # Named list matching group_cols
    numeric_cols = c() # Add parameter for numeric columns
    ) {
    # Get experiment object if name provided
    if (is.null(experiment_obj) && !is.null(experiment_name)) {
        experiment_obj <- get(experiment_name, envir = .GlobalEnv)
    }
    if (is.null(experiment_obj)) {
        stop("Must provide either experiment_name or experiment_obj")
    }

    # Get config from experiment object
    config <- experiment_obj$metadata$config
    # Check if config$data_dir exists and is not NULL
    if (is.null(config$data_dir)) {
        stop("config$data_dir is NULL in the experiment object")
    }
    
    # Update directory name from "cleaning" to "cleaned_data"
    cleaned_data_dir <- file.path(config$data_dir, "cleaned_data")
    dir.create(cleaned_data_dir, recursive = TRUE, showWarnings = FALSE)

    # Get data using the experiment name as prefix
    prefix <- experiment_obj$metadata$prefix
    sample_data_name <- paste0(prefix, "_sample_data")
    assay_data_name <- paste0(prefix, "_assay_data")

    # Extract data from experiment object
    sample_data <- experiment_obj$sample_data[[sample_data_name]]
    assay_data <- experiment_obj$assay_data[[assay_data_name]]

    # Step 1: Convert assay data to matrix ----------------------------------------
    if (!"gene_symbol" %in% colnames(assay_data)) {
        stop("Assay data must contain a 'gene_symbol' column")
    }

    rownames(assay_data) <- assay_data$gene_symbol
    assay_matrix <- assay_data[, -which(colnames(assay_data) == "gene_symbol")]
    assay_matrix <- as.matrix(assay_matrix)

    # Step 2: Clean and format sample metadata -----------------------------------
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
                    levels = factor_levels[[col]]
                )
            } else {
                sample_data[[col]] <- factor(sample_data[[col]])
            }
        }
    }

    # Step 3: Ensure consistent ordering ----------------------------------------
    header_order <- sample_data[[metadata_cols$sample_id]]

    # Verify all sample IDs exist in assay data
    missing_samples <- setdiff(header_order, colnames(assay_matrix))
    if (length(missing_samples) > 0) {
        stop(
            "Some sample IDs from sample_data are missing in assay_data: ",
            paste(missing_samples, collapse = ", ")
        )
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
                assay = dim(assay_data),
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

    # Save cleaned data to cleaned_data directory
    saveRDS(cleaned_data, file.path(cleaned_data_dir, paste0(cleaned_data_name, ".rds")))

    # Assign cleaned data to global environment
    assign(cleaned_data_name, cleaned_data, envir = .GlobalEnv)

    message(sprintf("Created cleaned data object '%s' in %s", cleaned_data_name, cleaned_data_dir))

    return(cleaned_data_name)
}

cleaned_data_name_1 <- clean_data_for_de(
    experiment_name = "xen_tran_2024_03",
    metadata_cols = list(
        sample_id = "sample_id",
        group_cols = c("sample_prep", "drug_treatment")
    ),
    factor_levels = list(
        sample_prep = c("snap frozen", "trizol"),
        drug_treatment = c("vehicle", "WB3")
    )
)

cleaned_data_name_2 <- clean_data_for_de(
    experiment_name = "xen_tran_2024_12",
    metadata_cols = list(
        sample_id = "sample_id",
        group_cols = c("condition", "timepoint", "dose", "drug", "replicate")
    ),
    factor_levels = list(
        condition = c("WB3_100", "Ketamine_500", "Vehicle", "WB3_200", "Ketamine_100", "Etomidate_100", "Propofol_100"),
        timepoint = c("T", "R"),
        dose = c("100", "200", "500"),
        drug = c("WB3", "Ketamine", "Vehicle", "Etomidate", "Propofol"),
        replicate = c("1", "2", "3", "4")
    ),
    numeric_cols = c("dose")
)

# Clean data for hum_tran_2024_03_deep_phenotype_muscle
cleaned_data_name_3 <- clean_data_for_de(
    experiment_name = "hum_tran_2024_03_deep_phenotype_muscle",
    metadata_cols = list(
        sample_id = "sample_id",
        group_cols = c("source_name", "group", "subject_status", "sex")
    ),
    factor_levels = list(
        group = c("HV", "PI-ME/CFS"),
        subject_status = c("healthy volunteer (HV)", 
                         "post-infectious Myalgic encephalomyelitis/chronic fatigue syndrome (PI-ME/CFS) volunteer"),
        sex = c("Male", "Female")
    )
)

# Clean data for hum_tran_2024_03_deep_phenotype_pbmcs
cleaned_data_name_4 <- clean_data_for_de(
    experiment_name = "hum_tran_2024_03_deep_phenotype_pbmcs",
    metadata_cols = list(
        sample_id = "sample_id",
        group_cols = c("group", "subject_status", "sex")
    ),
    factor_levels = list(
        group = c("HV", "PI-ME/CFS"),
        subject_status = c("healthy volunteer (HV)", 
                         "post-infectious Myalgic encephalomyelitis/chronic fatigue syndrome (PI-ME/CFS) volunteer"),
        sex = c("Male", "Female")
    )
)

# Load the cleaned data
cleaned_data <- readRDS("data/hum_tran_2024_03_deep_phenotype_pbmcs/cleaned_data/hum_tran_2024_03_deep_phenotype_pbmcs_cleaned.rds")

# Print first few rows
head(cleaned_data)

# Clean data for xen_tran_2025_02
cleaned_data_name_5 <- clean_data_for_de(
    experiment_name = "xen_tran_2025_02",
    metadata_cols = list(
        sample_id = "sample_id",
        group_cols = c("condition", "collection_day", "collection_time", "timepoint", "drug", "dose", "replicate")
    ),
    factor_levels = list(
        condition = c("DMSO_Wed_R", "DMSO_Wed_T", "DMSO_Mon_AM", "DMSO_Mon_PM", 
                     "Ket_500_Wed_R", "Ket_500_Wed_T", "MMR_Mon_AM", "MMR_Mon_PM", 
                     "MMR_Wed_R", "MMR_Wed_T", "WB3_100_Mon_AM", "WB3_100_Mon_PM", 
                     "WB3_100_Wed_R", "WB3_100_Wed_T", "WB3_200_Wed_R", "WB3_200_Wed_T", 
                     "WC22_50_Mon_AM", "WC22_50_Mon_PM", "WC22_50_Wed_R", "WC22_50_Wed_T"),
        collection_day = c("Mon", "Wed"),
        collection_time = c("AM", "PM", ""),
        timepoint = c("R", "T"),
        drug = c("DMSO", "Ketamine", "MMR", "WB3", "WC22"),
        replicate = c("1", "2", "3")
    ),
    numeric_cols = c("dose")
)

cleaned_data <- readRDS("data/xen_tran_2025_02/cleaned_data/xen_tran_2025_02_cleaned.rds")

