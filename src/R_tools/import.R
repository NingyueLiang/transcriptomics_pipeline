#' Data Import Functions
#' 
#' Import transcriptomics data from local files or Pluto API.

# Required libraries
suppressPackageStartupMessages({
    library(dplyr)    # For data manipulation
    library(tidyr)    # For data tidying
    library(readr)    # For reading CSV files
    library(magrittr) # For pipe operations
    library(pluto)    # For Pluto API interactions
    library(tibble)   # For tibble operations
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
import_config <- list(
    valid_species = c("xenopus", "human", "zebrafish", "mouse", "rat"),
    valid_data_types = c("transcriptomics", "genomics", "metabolomics", "proteomics", "scRNAseq"),
    default_paths = list(
        xenbase_file = "input_files/Xenbase_Xenopus_Orthology_Predictions.txt",
        hcop_file = "input_files/HCOP_Xenopus_Orthology_Predictions.txt"
    )
)

#' Validate input parameters
#' @param species Species name
#' @param data_type Type of omics data
#' @param data_source Source of data
validate_import_params <- function(species, data_type, data_source) {
    if (!species %in% import_config$valid_species) {
        stop("Invalid species. Must be one of: ", paste(import_config$valid_species, collapse = ", "))
    }
    
    if (!data_type %in% import_config$valid_data_types) {
        stop("Invalid data type. Must be one of: ", paste(import_config$valid_data_types, collapse = ", "))
    }
    
    if (!data_source %in% c("pluto", "local")) {
        stop("Invalid data source. Must be 'pluto' or 'local'")
    }
}

#' Create experiment configuration
#' @param experiment_prefix Prefix for experiment name
#' @param experiment_date Date of experiment
#' @param data_source Source of data ("pluto" or "local")
#' @param data_type Type of omics data
#' @param species Species name
#' @param metadata_columns Columns to use for metadata
#' @param additional_notes Additional notes about experiment
#' @return List containing experiment configuration
setup_experiment <- function(experiment_prefix, experiment_date,
                             data_source = "pluto",
                             data_type = "transcriptomics",
                             species = "xenopus",
                             metadata_columns = NULL,
                             additional_notes = NULL) {
    
    # Validate inputs
    validate_import_params(species, data_type, data_source)
    
    if (is.null(experiment_prefix) || nchar(experiment_prefix) == 0) {
        stop("experiment_prefix cannot be empty")
    }
    
    if (is.null(experiment_date) || nchar(experiment_date) == 0) {
        stop("experiment_date cannot be empty")
    }

    # Create base configuration
    config <- list(
        name = experiment_prefix,
        date = experiment_date,
        data_dir = here("data", experiment_prefix),
        results_dir = here("results", experiment_prefix),
        data_source = data_source,
        data_type = data_type,
        species = species,
        metadata_columns = metadata_columns,
        additional_notes = additional_notes %||% NA,
        xenbase_file = here(import_config$default_paths$xenbase_file),
        hcop_file = here(import_config$default_paths$hcop_file)
    )
    
    # Create necessary directories with error handling
    tryCatch({
        create_experiment_directories(config)
    }, error = function(e) {
        stop("Failed to create experiment directories: ", e$message)
    })
    
    return(config)
}

#' Create experiment directory structure
#' @param config Experiment configuration list
create_experiment_directories <- function(config) {
    # Main directories
    main_dirs <- c(config$data_dir, config$results_dir)
    for (dir in main_dirs) {
        if (!dir.exists(dir)) {
            if (!dir.create(dir, recursive = TRUE, showWarnings = FALSE)) {
                stop("Failed to create directory: ", dir)
            }
        }
    }
    
    # Data subdirectories
    data_subdirs <- c("cleaned_data")
    for (subdir in data_subdirs) {
        subdir_path <- file.path(config$data_dir, subdir)
        if (!dir.exists(subdir_path)) {
            if (!dir.create(subdir_path, showWarnings = FALSE)) {
                stop("Failed to create data subdirectory: ", subdir_path)
            }
        }
    }
    
    # Results subdirectories
    results_subdirs <- c("de", "mapping", "gsea")
    for (subdir in results_subdirs) {
        subdir_path <- file.path(config$results_dir, subdir)
        if (!dir.exists(subdir_path)) {
            if (!dir.create(subdir_path, showWarnings = FALSE)) {
                stop("Failed to create results subdirectory: ", subdir_path)
            }
        }
    }
}

#' Handle Pluto API credentials with validation
#' @param api_key API key (optional)
#' @param save_credentials Whether to save credentials
#' @return API key string
get_pluto_api_key <- function(api_key = NULL, save_credentials = FALSE) {
    # Check provided API key
    if (!is.null(api_key)) {
        if (!is.character(api_key) || nchar(api_key) == 0) {
            stop("API key must be a non-empty character string")
        }
        return(api_key)
    }
    
    # Check environment variable
    key <- Sys.getenv("PLUTO_API_KEY")
    if (key != "") {
        return(key)
    }
    
    # Interactive input
    if (interactive()) {
        message("No Pluto API key found. Please enter your API key:")
        key <- readline(prompt = "API key: ")
        
        if (save_credentials) {
            save_credentials_to_renviron(key)
        }
        return(key)
    } else {
        stop("No Pluto API key found and not in interactive mode. Please set PLUTO_API_KEY environment variable.")
    }
}

#' Save API key to .Renviron file
#' @param api_key API key to save
save_credentials_to_renviron <- function(api_key) {
    tryCatch({
        renviron_path <- file.path(Sys.getenv("HOME"), ".Renviron")
        
        # Read existing content
        existing_content <- if (file.exists(renviron_path)) {
            readLines(renviron_path)
        } else {
            character(0)
        }
        
        # Remove existing PLUTO_API_KEY entries
        existing_content <- existing_content[!grepl("^PLUTO_API_KEY=", existing_content)]
        
        # Add new key
        new_content <- c(existing_content, sprintf('PLUTO_API_KEY="%s"', api_key))
        
        # Write to file
        writeLines(new_content, renviron_path)
        message("API key saved to .Renviron. Please restart R for changes to take effect.")
        
    }, error = function(e) {
        warning("Failed to save API key: ", e$message)
    })
}

# Modified extract_conditions function to use sample data
extract_conditions <- function(sample_data, metadata_columns) {
  # Initialize list to store different condition types
  condition_types <- list()

  # Process each metadata column
  for (col in metadata_columns) {
    if (col %in% colnames(sample_data)) {
      # Get unique values for this column
      unique_values <- unique(sample_data[[col]])
      # Remove NA values
      unique_values <- unique_values[!is.na(unique_values)]
      # Store in condition_types list
      condition_types[[col]] <- as.character(unique_values)
    }
  }

  return(condition_types)
}

# Function to create experiment metadata
create_experiment_metadata <- function(data_source, data_type, species, experiment_date = NULL,
                                      description = NULL, experimenter = NULL, additional_notes = NULL) {
  metadata <- list(
    data_source = data_source,
    data_type = data_type,
    species = species,
    experiment_date = experiment_date,
    description = description,
    experimenter = experimenter,
    additional_notes = additional_notes,
    import_timestamp = Sys.time()
  )

  return(metadata)
}

# Function to make valid R object names
make_valid_name <- function(name) {
  # Replace invalid characters with underscores
  name <- gsub("[^A-Za-z0-9_]", "_", name)
  # Ensure it starts with a letter or underscore
  if (!grepl("^[A-Za-z_]", name)) {
    name <- paste0("exp_", name)
  }
  return(name)
}

# Function to create experiment object
create_experiment_object <- function(metadata, data_dictionary, sample_data, assay_data, experiment_name) {
  experiment_obj <- list(
    metadata = metadata,
    data_dictionary = data_dictionary,
    sample_data = setNames(list(sample_data), paste0(experiment_name, "_sample_data")),
    assay_data = setNames(list(assay_data), paste0(experiment_name, "_assay_data"))
  )

  class(experiment_obj) <- "omics_experiment"

  return(experiment_obj)
}

#' Suggest gene column from available column names
#' @param column_names Vector of column names
#' @return Name of suggested gene column or NULL
suggest_gene_column <- function(column_names) {
    common_gene_names <- c(
        "gene_symbol", "Gene", "GENE", "gene", "Gene_Symbol", "GENE_SYMBOL",
        "gene_name", "GeneName", "GENE_NAME", "symbol", "Symbol", "SYMBOL",
        "Protein_ID", "protein_id", "PROTEIN_ID"
    )
    
    # Find matching columns
    matches <- column_names[tolower(column_names) %in% tolower(common_gene_names)]
    
    if (length(matches) == 0) {
        return(NULL)
    }
    return(matches[1])
}

#' Standardize gene symbol column in assay data
#' @param assay_data Data frame containing assay data
#' @return Data frame with standardized gene_symbol column
standardize_gene_symbol_column <- function(assay_data) {
    if (!"gene_symbol" %in% colnames(assay_data)) {
        # Check if it might be named differently
        possible_gene_cols <- c("Gene", "GENE", "gene", "Gene_Symbol", "GENE_SYMBOL",
                                "Protein_ID", "protein_id", "PROTEIN_ID")
        gene_col <- which(colnames(assay_data) %in% possible_gene_cols)[1]
        
        if (!is.na(gene_col)) {
            # Rename the found column to gene_symbol
            colnames(assay_data)[gene_col] <- "gene_symbol"
        } else {
            stop("No gene symbol column found in assay data")
        }
    }
    
    # Move gene_symbol to first column if it's not already
    if (which(colnames(assay_data) == "gene_symbol") != 1) {
        assay_data <- assay_data %>%
            select(gene_symbol, everything())
    }
    
    return(assay_data)
}

#' Handle gene symbol column interactively for local imports
#' @param assay_data Data frame containing assay data
#' @return Data frame with standardized gene_symbol column
handle_gene_symbol_column_interactive <- function(assay_data) {
    if (!"gene_symbol" %in% colnames(assay_data)) {
        suggested_col <- suggest_gene_column(colnames(assay_data))
        
        if (!is.null(suggested_col)) {
            if (interactive()) {
                message("Found potential gene identifier column: '", suggested_col, "'")
                message("Would you like to use this as the gene symbol column? (yes/no)")
                answer <- readline(prompt = "Answer: ")
                
                if (tolower(answer) %in% c("yes", "y")) {
                    colnames(assay_data)[colnames(assay_data) == suggested_col] <- "gene_symbol"
                    message("Column '", suggested_col, "' renamed to 'gene_symbol'")
                } else {
                    assay_data <- select_gene_column_interactive(assay_data)
                }
            } else {
                # Non-interactive mode - use suggested column
                colnames(assay_data)[colnames(assay_data) == suggested_col] <- "gene_symbol"
                message("Using suggested gene column: '", suggested_col, "'")
            }
        } else {
            if (interactive()) {
                assay_data <- select_gene_column_interactive(assay_data)
            } else {
                stop("No gene identifier column found and not in interactive mode")
            }
        }
    }
    
    # Move gene_symbol to first column if it's not already
    if (which(colnames(assay_data) == "gene_symbol") != 1) {
        assay_data <- assay_data %>%
            select(gene_symbol, everything())
    }
    
    return(assay_data)
}

#' Interactive gene column selection
#' @param assay_data Data frame containing assay data
#' @return Data frame with gene_symbol column
select_gene_column_interactive <- function(assay_data) {
    message("No gene identifier column found. Available columns:")
    for (i in seq_along(colnames(assay_data))) {
        message(i, ": ", colnames(assay_data)[i])
    }
    message("\nWhich column contains gene identifiers? (Enter column number)")
    col_num <- as.numeric(readline(prompt = "Column number: "))
    
    if (!is.na(col_num) && col_num <= ncol(assay_data) && col_num > 0) {
        old_name <- colnames(assay_data)[col_num]
        colnames(assay_data)[col_num] <- "gene_symbol"
        message("Column '", old_name, "' renamed to 'gene_symbol'")
        return(assay_data)
    } else {
        stop("Invalid column selection")
    }
}

#' Main import function with comprehensive error handling
#' @param data_source Source of data ("pluto" or "local")
#' @param data_type Type of omics data
#' @param species Species name
#' @param experiment_id Experiment ID for Pluto
#' @param experiment_name Name of experiment
#' @param sample_path Path to sample data file
#' @param assay_path Path to assay data file
#' @param api_key Pluto API key
#' @param save_credentials Whether to save credentials
#' @param experiment_date Date of experiment
#' @param description Experiment description
#' @param experimenter Name of experimenter
#' @param additional_notes Additional notes
#' @param metadata_columns Metadata columns to use
#' @return Experiment object
import_omics_data <- function(
    data_source = c("pluto", "local"),
    data_type = c("transcriptomics", "genomics", "metabolomics", "proteomics", "scRNAseq"),
    species = c("xenopus", "human", "zebrafish", "mouse", "rat"),
    experiment_id = NULL,
    experiment_name = NULL,
    sample_path = NULL,
    assay_path = NULL,
    api_key = NULL,
    save_credentials = FALSE,
    experiment_date = NULL,
    description = NULL,
    experimenter = NULL,
    additional_notes = NULL,
    metadata_columns = NULL) {
    
    # Match and validate arguments
    data_source <- match.arg(data_source)
    data_type <- match.arg(data_type)
    species <- match.arg(species)
    
    # Validate required parameters
    if (data_source == "pluto" && is.null(experiment_id)) {
        stop("experiment_id is required for Pluto import")
    }
    
    if (data_source == "local" && is.null(assay_path)) {
        stop("assay_path is required for local import")
    }
    
    if (is.null(experiment_name)) {
        experiment_name <- experiment_id
    }
    
    if (is.null(experiment_name)) {
        stop("Either experiment_name or experiment_id must be provided")
    }

    # Handle Pluto authentication if needed
    if (data_source == "pluto") {
        tryCatch({
            key <- get_pluto_api_key(api_key, save_credentials)
            pluto::pluto_login(key)
        }, error = function(e) {
            stop("Pluto authentication failed: ", e$message)
        })
    }

    # Create initial metadata
    metadata <- create_experiment_metadata(
        data_source = data_source,
        data_type = data_type,
        species = species,
        experiment_date = experiment_date,
        description = description,
        experimenter = experimenter,
        additional_notes = additional_notes
    )
    metadata$metadata_columns <- metadata_columns
    
    # Create config for the experiment
    config <- list(
        name = experiment_name,
        date = experiment_date,
        data_dir = here("data", experiment_name),
        results_dir = here("results", experiment_name),
        data_source = data_source,
        data_type = data_type,
        species = species,
        metadata_columns = metadata_columns,
        additional_notes = additional_notes
    )
    metadata$config <- config
    metadata$prefix <- experiment_name

    # Create valid experiment name
    experiment_name <- make_valid_name(experiment_name %||% experiment_id)

    # Import data based on source
    if (data_source == "pluto") {
        tryCatch({
            message("Importing data from Pluto for experiment: ", experiment_id)
            
            sample_data <- pluto::pluto_read_data(
                experiment_id = experiment_id,
                table_type = "sample"
            )
            assay_data <- pluto::pluto_read_data(
                experiment_id = experiment_id,
                table_type = "assay"
            )
            
            # Validate imported data
            if (is.null(sample_data) || nrow(sample_data) == 0) {
                stop("No sample data found for experiment: ", experiment_id)
            }
            if (is.null(assay_data) || nrow(assay_data) == 0) {
                stop("No assay data found for experiment: ", experiment_id)
            }

            # Ensure gene_symbol column exists and is the first column
            assay_data <- standardize_gene_symbol_column(assay_data)
            
        }, error = function(e) {
            stop("Failed to read from Pluto: ", e$message)
        })
    } else {
        # Local import
        tryCatch({
            message("Importing data from local files...")
            
            # Validate file paths
            if (!file.exists(assay_path)) {
                stop("Assay file not found: ", assay_path)
            }
            
            # Read assay data
            assay_data <- readr::read_csv(assay_path, show_col_types = FALSE)
            # Convert tibble to data frame to allow rownames
            assay_data <- as.data.frame(assay_data)
            
            # Validate assay data
            if (is.null(assay_data) || nrow(assay_data) == 0) {
                stop("No data found in assay file: ", assay_path)
            }
            
            # Handle gene symbol column
            assay_data <- handle_gene_symbol_column_interactive(assay_data)
            
            # Handle sample data creation
            if (is.null(sample_path)) {
                sample_cols <- setdiff(colnames(assay_data), "gene_symbol")
                sample_data <- data.frame(
                    sample_name = sample_cols,
                    condition = "default",
                    stringsAsFactors = FALSE
                )
            } else {
                if (!file.exists(sample_path)) {
                    stop("Sample file not found: ", sample_path)
                }
                sample_data <- readr::read_csv(sample_path, show_col_types = FALSE)
                # Convert tibble to data frame
                sample_data <- as.data.frame(sample_data)
            }
            
        }, error = function(e) {
            stop("Failed to read local files: ", e$message)
        })
    }

    # Handle duplicate gene symbols by aggregating (summing) counts
    if ("gene_symbol" %in% colnames(assay_data)) {
        if (any(duplicated(assay_data$gene_symbol))) {
            n_dups <- sum(duplicated(assay_data$gene_symbol))
            message(sprintf("Found %d duplicate gene symbols - aggregating by sum", n_dups))
            numeric_cols <- setdiff(colnames(assay_data), "gene_symbol")
            assay_data <- assay_data %>%
                dplyr::group_by(gene_symbol) %>%
                dplyr::summarise(across(all_of(numeric_cols), sum, na.rm = TRUE)) %>%
                as.data.frame()
        }
        rownames(assay_data) <- assay_data$gene_symbol
    }
    
    # Create data dictionary
    data_dict <- create_data_dictionary(sample_data, assay_data)
    
    # Create experiment object
    experiment_obj <- list(
        metadata = metadata,
        data_dictionary = data_dict,
        sample_data = setNames(list(sample_data), paste0(experiment_name, "_sample_data")),
        assay_data = setNames(list(assay_data), paste0(experiment_name, "_assay_data"))
    )
    
    # Extract conditions only if metadata columns are provided and not NULL
    if (!is.null(metadata_columns) && length(metadata_columns) > 0) {
        experiment_obj$metadata$conditions <- extract_conditions(sample_data, metadata_columns)
    } else {
        experiment_obj$metadata$conditions <- list(condition = "default")
    }
    
    class(experiment_obj) <- "omics_experiment"
    
    # Save experiment object to file instead of global environment
    tryCatch({
        save_experiment_object(experiment_obj, experiment_name)
        message("Created and saved experiment object: ", experiment_name)
    }, error = function(e) {
        stop("Failed to save experiment object: ", e$message)
    })
    
    return(experiment_obj)
}

#' Create data dictionary for experiment
#' @param sample_data Sample metadata data frame
#' @param assay_data Assay data frame
#' @return List containing data dictionary
create_data_dictionary <- function(sample_data, assay_data) {
    list(
        sample_data_columns = data.frame(
            column_name = names(sample_data),
            description = NA,
            data_type = sapply(sample_data, class),
            stringsAsFactors = FALSE
        ),
        assay_data_columns = data.frame(
            column_name = names(assay_data),
            description = NA,
            data_type = sapply(assay_data, class),
            stringsAsFactors = FALSE
        )
    )
}

#' Save experiment object to file
#' @param experiment_obj Experiment object to save
#' @param experiment_name Name of experiment
save_experiment_object <- function(experiment_obj, experiment_name) {
    # Create data directory if it doesn't exist
    data_dir <- here("data", experiment_name)
    if (!dir.exists(data_dir)) {
        dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Save to RDS file
    save_path <- file.path(data_dir, paste0(experiment_name, ".rds"))
    saveRDS(experiment_obj, save_path)
    message("Experiment object saved to: ", save_path)
}

# Convenience function for importing with config
import_omics_data_with_config <- function(config, ...) {
  return(import_omics_data(
    data_source = config$data_source,
    data_type = config$data_type,
    species = config$species,
    experiment_date = config$date,
    metadata_columns = config$metadata_columns,
    additional_notes = config$additional_notes,
    ...
  ))
}

# Standalone TPP import function
import_tpp_data <- function(experiment_prefix,
                            experiment_date,
                            species = c("xenopus", "human", "zebrafish", "mouse", "rat"),
                            experiment_id = NULL,
                            api_key = NULL,
                            save_credentials = FALSE,
                            description = NULL,
                            experimenter = NULL,
                            additional_notes = NULL) {
  species <- match.arg(species)

  # Create experiment config
  config <- setup_experiment(
    experiment_prefix = experiment_prefix,
    experiment_date = experiment_date,
    data_source = "pluto",
    data_type = "proteomics",
    species = species,
    metadata_columns = NULL,
    additional_notes = additional_notes
  )

  # Import TPP data
  tpp_obj <- import_omics_data_with_config(
    config = config,
    experiment_id = experiment_id,
    api_key = api_key,
    save_credentials = save_credentials,
    description = description,
    experimenter = experimenter
  )

  return(tpp_obj)
}

# Standalone binding import function
import_binding_data <- function(experiment_prefix,
                                experiment_date,
                                species = c("xenopus", "human", "zebrafish", "mouse", "rat"),
                                experiment_id = NULL,
                                api_key = NULL,
                                save_credentials = FALSE,
                                description = NULL,
                                experimenter = NULL,
                                additional_notes = NULL) {
  species <- match.arg(species)

  # Create experiment config
  config <- setup_experiment(
    experiment_prefix = experiment_prefix,
    experiment_date = experiment_date,
    data_source = "pluto",
    data_type = "genomics",
    species = species,
    metadata_columns = NULL,
    additional_notes = additional_notes
  )

  # Import binding data
  binding_obj <- import_omics_data_with_config(
    config = config,
    experiment_id = experiment_id,
    api_key = api_key,
    save_credentials = save_credentials,
    description = description,
    experimenter = experimenter
  )

  return(binding_obj)
}

message("Import functions loaded successfully. Use import_omics_data() or import_omics_data_with_config() to import data.")
