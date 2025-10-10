#####
# Import Functions

# Required libraries
library(dplyr) # For data manipulation
library(tidyr) # For data tidying
library(readr) # For reading CSV files
library(magrittr) # For pipe operations
library(pluto) # For Pluto API interactions
library(tibble) # For tibble operations

# Configuration function
setup_experiment <- function(experiment_prefix, experiment_date,
                             data_source = "pluto",
                             data_type = "transcriptomics",
                             species = "xenopus",
                             metadata_columns = NULL,
                             additional_notes = NULL) {
  # Validate species
  valid_species <- c("xenopus", "human", "zebrafish", "mouse", "rat")
  if (!species %in% valid_species) {
    stop("Invalid species. Must be one of: ", paste(valid_species, collapse = ", "))
  }

  # Validate data type
  valid_data_types <- c("transcriptomics", "genomics", "metabolomics", "proteomics", "scRNAseq")
  if (!data_type %in% valid_data_types) {
    stop("Invalid data type. Must be one of: ", paste(valid_data_types, collapse = ", "))
  }

  # Create base configuration
  config <- list(
    name = experiment_prefix,
    date = experiment_date,
    data_dir = file.path("data", experiment_prefix),
    results_dir = file.path("results", experiment_prefix),
    data_source = data_source,
    data_type = data_type,
    species = species,
    metadata_columns = metadata_columns,
    additional_notes = additional_notes %||% NA,
    xenbase_file = "data_files/Xenbase_Xenopus_Orthology_Predictions.txt",
    hcop_file = "data_files/HCOP_Xenopus_Orthology_Predictions.txt"
  )

  # Create necessary directories
  for (dir in c(config$data_dir, config$results_dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Create data subdirectories (removed "raw")
  data_subdirs <- c("cleaned_data")
  for (subdir in data_subdirs) {
    dir.create(file.path(config$data_dir, subdir), showWarnings = FALSE)
  }

  # Create results subdirectories
  results_subdirs <- c("de", "mapping", "gsea")
  for (subdir in results_subdirs) {
    dir.create(file.path(config$results_dir, subdir), showWarnings = FALSE)
  }

  return(config)
}

# Function to handle Pluto API credentials
get_pluto_api_key <- function(api_key = NULL) {
  if (!is.null(api_key)) {
    return(api_key)
  }

  key <- Sys.getenv("PLUTO_API_KEY")
  if (key != "") {
    return(key)
  }

  message("No Pluto API key found. Please enter your API key:")
  key <- readline(prompt = "API key: ")

  save_key <- readline(prompt = "Save API key for future use? (yes/no): ")
  if (tolower(save_key) %in% c("yes", "y")) {
    renviron_path <- file.path(Sys.getenv("HOME"), ".Renviron")
    if (!file.exists(renviron_path)) file.create(renviron_path)
    writeLines(sprintf('PLUTO_API_KEY="%s"', key), renviron_path, append = TRUE)
    message("API key saved. Please restart R for changes to take effect.")
  }

  return(key)
}

# Modified extract_conditions function to use sample data
extract_conditions <- function(sample_data, metadata_columns) {
  # Initialize list to store different condition types
  condition_types <- list()

  # For each metadata column, collect all unique values
  for (col_name in metadata_columns) {
    if (col_name %in% colnames(sample_data)) {
      unique_values <- unique(sample_data[[col_name]])
      if (length(unique_values) > 1) {
        condition_types[[col_name]] <- unique_values
      }
    }
  }

  return(condition_types)
}

# Helper function for valid names
make_valid_name <- function(name) {
  if (is.null(name)) {
    name <- format(Sys.time(), "exp_%Y%m%d_%H%M%S")
  }
  return(make.names(name))
}

# Metadata creation function
create_experiment_metadata <- function(data_source, data_type, species,
                                       experiment_date = NULL,
                                       description = NULL,
                                       experimenter = NULL,
                                       additional_notes = NULL) {
  metadata <- list(
    data_source = data_source,
    data_type = data_type,
    species = species,
    experiment_date = experiment_date %||% Sys.Date(),
    description = description %||% "No description provided",
    experimenter = experimenter %||% Sys.info()["user"],
    additional_notes = additional_notes %||% NA,
    creation_timestamp = Sys.time()
  )
  return(metadata)
}

# Modified import function with config and automatic fixes
import_omics_data_with_config <- function(config, experiment_id = NULL,
                                          sample_path = NULL, assay_path = NULL,
                                          api_key = NULL, save_credentials = FALSE,
                                          description = NULL, experimenter = NULL,
                                          additional_notes = NULL) {
  experiment_obj <- import_omics_data(
    data_source = config$data_source,
    data_type = config$data_type,
    species = config$species,
    experiment_id = experiment_id,
    experiment_name = config$name,
    sample_path = sample_path,
    assay_path = assay_path,
    api_key = api_key,
    save_credentials = save_credentials,
    experiment_date = config$date,
    description = description,
    experimenter = experimenter,
    additional_notes = additional_notes,
    metadata_columns = config$metadata_columns
  )

  # Add config information to metadata
  experiment_obj$metadata$config <- config
  experiment_obj$metadata$prefix <- config$name
  experiment_obj$metadata$date_format <- "month"
  experiment_obj$metadata$creation_date <- Sys.Date()
  experiment_obj$metadata$experiment_date <- as.Date(paste0(config$date, "-01"))

  # Fix the assay and sample data naming
  names(experiment_obj$assay_data)[1] <- paste0(config$name, "_assay_data")
  names(experiment_obj$sample_data)[1] <- paste0(config$name, "_sample_data")

  # Set the class
  class(experiment_obj) <- "omics_experiment"

  # Save imported data to the experiment's data directory
  saveRDS(
    experiment_obj,
    file.path(config$data_dir, paste0(config$name, ".rds"))
  )

  return(experiment_obj)
}

# Add helper function for column mapping
suggest_gene_column <- function(column_names) {
  common_gene_names <- c(
    "gene_symbol", "Gene", "GENE", "gene", "Gene_Symbol", "GENE_SYMBOL",
    "gene_name", "GeneName", "GENE_NAME", "symbol", "Symbol", "SYMBOL"
  )

  # Find matching columns
  matches <- column_names[tolower(column_names) %in% tolower(common_gene_names)]

  if (length(matches) == 0) {
    return(NULL)
  }
  return(matches[1])
}

# Main import function
import_omics_data <- function(
    # Basic parameters
    data_source = c("pluto", "local"),
    data_type = c("transcriptomics", "genomics", "metabolomics", "proteomics", "scRNAseq"),
    species = c("xenopus", "human", "zebrafish", "mouse", "rat"),
    # Data source parameters
    experiment_id = NULL,
    experiment_name = NULL,
    sample_path = NULL,
    assay_path = NULL,
    api_key = NULL,
    save_credentials = FALSE,
    # Optional metadata
    experiment_date = NULL,
    description = NULL,
    experimenter = NULL,
    additional_notes = NULL,
    metadata_columns = NULL) {
  # Match and validate arguments
  data_source <- match.arg(data_source)
  data_type <- match.arg(data_type)
  species <- match.arg(species)

  # Handle Pluto authentication if needed
  if (data_source == "pluto") {
    tryCatch(
      {
        key <- get_pluto_api_key(api_key)

        # Save credentials if requested
        if (save_credentials && !is.null(api_key)) {
          renviron_path <- file.path(Sys.getenv("HOME"), ".Renviron")
          if (!file.exists(renviron_path)) {
            file.create(renviron_path)
          }
          current_content <- readLines(renviron_path)
          current_content <- current_content[!grepl("^PLUTO_API_KEY=", current_content)]
          writeLines(c(current_content, sprintf('PLUTO_API_KEY="%s"', api_key)), renviron_path)
          message("API key saved to .Renviron. Restart R for changes to take effect.")
        }

        pluto::pluto_login(key)
      },
      error = function(e) {
        stop("Pluto authentication failed: ", e$message)
      }
    )
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

  # Add metadata columns to metadata
  metadata$metadata_columns <- metadata_columns

  # Create valid experiment name
  experiment_name <- make_valid_name(experiment_name %||% experiment_id)

  # Import data based on source
  if (data_source == "pluto") {
    if (is.null(experiment_id)) {
      stop("experiment_id is required for Pluto import")
    }

    tryCatch(
      {
        sample_data <- pluto::pluto_read_data(
          experiment_id = experiment_id,
          table_type = "sample"
        )
        assay_data <- pluto::pluto_read_data(
          experiment_id = experiment_id,
          table_type = "assay"
        )

        # Ensure gene_symbol column exists and is the first column
        if (!"gene_symbol" %in% colnames(assay_data)) {
          # Check if it might be named differently
          possible_gene_cols <- c("Gene", "GENE", "gene", "Gene_Symbol", "GENE_SYMBOL")
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
      },
      error = function(e) {
        stop("Failed to read from Pluto: ", e$message)
      }
    )
  } else {
    if (is.null(assay_path)) {
      stop("assay_path is required for local import")
    }

    tryCatch(
      {
        # Read assay data
        assay_data <- readr::read_csv(assay_path)

        # Handle gene symbol column
        if (!"gene_symbol" %in% colnames(assay_data)) {
          suggested_col <- suggest_gene_column(colnames(assay_data))

          if (!is.null(suggested_col)) {
            message("Found potential gene identifier column: '", suggested_col, "'")
            message("Would you like to use this as the gene symbol column? (yes/no)")
            answer <- readline(prompt = "Answer: ")

            if (tolower(answer) %in% c("yes", "y")) {
              # Rename the column
              colnames(assay_data)[colnames(assay_data) == suggested_col] <- "gene_symbol"
              message("Column '", suggested_col, "' renamed to 'gene_symbol'")
            } else {
              # Let user choose column
              message("\nAvailable columns:")
              for (i in seq_along(colnames(assay_data))) {
                message(i, ": ", colnames(assay_data)[i])
              }
              message("\nWhich column contains gene identifiers? (Enter column number)")
              col_num <- as.numeric(readline(prompt = "Column number: "))

              if (!is.na(col_num) && col_num <= ncol(assay_data)) {
                old_name <- colnames(assay_data)[col_num]
                colnames(assay_data)[col_num] <- "gene_symbol"
                message("Column '", old_name, "' renamed to 'gene_symbol'")
              } else {
                stop("Invalid column selection")
              }
            }
          } else {
            # No suggestion found - ask user directly
            message("No gene identifier column found. Available columns:")
            for (i in seq_along(colnames(assay_data))) {
              message(i, ": ", colnames(assay_data)[i])
            }
            message("\nWhich column contains gene identifiers? (Enter column number)")
            col_num <- as.numeric(readline(prompt = "Column number: "))

            if (!is.na(col_num) && col_num <= ncol(assay_data)) {
              old_name <- colnames(assay_data)[col_num]
              colnames(assay_data)[col_num] <- "gene_symbol"
              message("Column '", old_name, "' renamed to 'gene_symbol'")
            } else {
              stop("Invalid column selection")
            }
          }
        }

        # Move gene_symbol to first column if it's not already
        if (which(colnames(assay_data) == "gene_symbol") != 1) {
          assay_data <- assay_data %>%
            select(gene_symbol, everything())
        }

        # Handle sample data creation
        if (is.null(sample_path)) {
          sample_cols <- setdiff(colnames(assay_data), "gene_symbol")
          sample_data <- data.frame(
            sample_name = sample_cols,
            condition = "default"
          )
        } else {
          sample_data <- readr::read_csv(sample_path)
        }
      },
      error = function(e) {
        stop("Failed to read local files: ", e$message)
      }
    )
  }

  # Create data dictionary
  data_dict <- list(
    sample_data_columns = data.frame(
      column_name = names(sample_data),
      description = NA,
      data_type = sapply(sample_data, class)
    ),
    assay_data_columns = data.frame(
      column_name = names(assay_data),
      description = NA,
      data_type = sapply(assay_data, class)
    )
  )

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

  # Safely assign to global environment
  tryCatch(
    {
      assign(experiment_name, experiment_obj, envir = .GlobalEnv)
      message("Created experiment object: ", experiment_name)
    },
    error = function(e) {
      stop("Failed to create experiment object: ", e$message)
    }
  )

  return(experiment_obj)
}

# Standalone TPP import function
import_tpp_data <- function(experiment_prefix,
                            experiment_date,
                            species = c("xenopus", "human", "zebrafish", "mouse", "rat"),
                            tpp_data_path,
                            protein_info_path = NULL,
                            description = NULL,
                            experimenter = NULL,
                            additional_notes = NULL) {
  species <- match.arg(species)

  # Create experiment name
  experiment_name <- paste(substr(species, 1, 3), "tpp", experiment_date, sep = "_")

  # Create directories
  data_dir <- file.path("data", experiment_name)
  results_dir <- file.path("results", experiment_name)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # Read TPP data
  if (!file.exists(tpp_data_path)) {
    stop("TPP data file not found: ", tpp_data_path)
  }
  tpp_data <- readr::read_csv(tpp_data_path)

  # Read protein info if provided
  protein_info <- if (!is.null(protein_info_path)) {
    if (!file.exists(protein_info_path)) {
      stop("Protein info file not found: ", protein_info_path)
    }
    readr::read_csv(protein_info_path)
  } else {
    NULL
  }

  # Create metadata
  metadata <- create_experiment_metadata(
    data_source = "local",
    data_type = "tpp",
    species = species,
    experiment_date = experiment_date,
    description = description,
    experimenter = experimenter,
    additional_notes = additional_notes
  )

  # Create TPP experiment object
  tpp_obj <- list(
    metadata = metadata,
    tpp_data = tpp_data,
    protein_info = protein_info,
    data_dictionary = list(
      tpp_columns = data.frame(
        column_name = names(tpp_data),
        description = NA,
        data_type = sapply(tpp_data, class)
      )
    )
  )

  class(tpp_obj) <- "tpp_experiment"

  # Save the object
  saveRDS(tpp_obj, file.path(data_dir, paste0(experiment_name, ".rds")))

  return(tpp_obj)
}

# Standalone binding targets import function
import_binding_data <- function(experiment_prefix,
                                experiment_date,
                                species = c("xenopus", "human", "zebrafish", "mouse", "rat"),
                                targets_path,
                                binding_data_path = NULL,
                                description = NULL,
                                experimenter = NULL,
                                additional_notes = NULL) {
  species <- match.arg(species)

  # Create experiment name
  experiment_name <- paste(substr(species, 1, 3), "drug_targets", experiment_date, sep = "_")

  # Create directories
  data_dir <- file.path("data", experiment_name)
  results_dir <- file.path("results", experiment_name)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # Read targets data
  if (!file.exists(targets_path)) {
    stop("Binding targets file not found: ", targets_path)
  }
  targets_data <- readr::read_csv(targets_path)

  # Read additional binding data if provided
  binding_data <- if (!is.null(binding_data_path)) {
    if (!file.exists(binding_data_path)) {
      stop("Binding data file not found: ", binding_data_path)
    }
    readr::read_csv(binding_data_path)
  } else {
    NULL
  }

  # Create metadata
  metadata <- create_experiment_metadata(
    data_source = "local",
    data_type = "binding_targets",
    species = species,
    experiment_date = experiment_date,
    description = description,
    experimenter = experimenter,
    additional_notes = additional_notes
  )

  # Create binding experiment object
  binding_obj <- list(
    metadata = metadata,
    targets_data = targets_data,
    binding_data = binding_data,
    data_dictionary = list(
      targets_columns = data.frame(
        column_name = names(targets_data),
        description = NA,
        data_type = sapply(targets_data, class)
      )
    )
  )

  class(binding_obj) <- "binding_experiment"

  # Save the object
  saveRDS(binding_obj, file.path(data_dir, paste0(experiment_name, ".rds")))

  return(binding_obj)
}

####DARPA-ABC Project 

# Setup first experiment
experiment_config_1 <- setup_experiment(
  experiment_prefix = "xen_tran_2024_03",
  experiment_date = "2024-03",
  data_source = "pluto",
  data_type = "transcriptomics",
  species = "xenopus",
  metadata_columns = c("drug_treatment", "sample_prep"),
  additional_notes = "Initial experiment"
)

# Import first experiment
xen_tran_2024_03 <- import_omics_data_with_config(
  config = experiment_config_1,
  experiment_id = "PLX119351",
  api_key = Sys.getenv("PLUTO_API_KEY"),
  save_credentials = TRUE,
  description = "Xenopus transcriptomics experiment March 2024",
  experimenter = "Veda Kim"
)

# Access metadata columns
xen_tran_2024_03$metadata$metadata_columns
xen_tran_2024_03$metadata$conditions$drug_treatment
xen_tran_2024_03$metadata$conditions$sample_prep
xen_tran_2024_03$metadata
xen_tran_2024_03$data_dictionary
xen_tran_2024_03$sample_data
head(xen_tran_2024_03$assay_data[[1]])

# Later, to run another experiment:
experiment_config_2 <- setup_experiment(
  experiment_prefix = "xen_tran_2024_12",
  experiment_date = "2024-12",
  data_source = "pluto",
  data_type = "transcriptomics",
  species = "xenopus",
  metadata_columns = c("condition", "drug", "dose", "timepoint", "replicate"),
  additional_notes = "Second experiment"
)

# Import second experiment
xen_tran_2024_12 <- import_omics_data_with_config(
  config = experiment_config_2,
  experiment_id = "PLX073248",
  api_key = Sys.getenv("PLUTO_API_KEY"),
  description = "Xenopus transcriptomics experiment December 2024",
  experimenter = "Allison Grossberg"
)

# Access the created object
xen_tran_2024_12$metadata$metadata_columns
xen_tran_2024_12$metadata$conditions$condition
xen_tran_2024_12$metadata$conditions$drug
xen_tran_2024_12$metadata$conditions$dose
xen_tran_2024_12$metadata$conditions$timepoint
xen_tran_2024_12$metadata$conditions$replicate
xen_tran_2024_12$metadata
xen_tran_2024_12$data_dictionary
xen_tran_2024_12$sample_data
head(xen_tran_2024_12$assay_data[[1]])


####Neurobots vs Xenobots Project 

# Setup experiment config
experiment_config_3 <- setup_experiment(
  experiment_prefix = "NBvsXB_xen_tran_2024_12",
  experiment_date = "2024-12",
  data_source = "local", 
  data_type = "transcriptomics",
  species = "xenopus",
  metadata_columns = NULL,
  additional_notes = "Neurobots vs Xenobots"
)

# Import experiment with local files
NBvsXB_xen_tran_2024_12 <- import_omics_data_with_config(
  config = experiment_config_3,
  assay_path = "/Users/wyssuser/Desktop/NBvsXB_Novogene_All_Genes.csv", # Add path to assay data file
  description = "Neurobots vs Xenobots transcriptomics experiment December 2024",
  experimenter = "Haleh Fotowat"
)

str(NBvsXB_xen_tran_2024_12)

# Access metadata columns
NBvsXB_xen_tran_2024_12$metadata$metadata_columns
NBvsXB_xen_tran_2024_12$metadata$conditions$drug_treatment
NBvsXB_xen_tran_2024_12$metadata$conditions$sample_prep
NBvsXB_xen_tran_2024_12$metadata
NBvsXB_xen_tran_2024_12$data_dictionary
NBvsXB_xen_tran_2024_12$sample_data
head(NBvsXB_xen_tran_2024_12$assay_data[[1]])


####ME/CFS Project

experiment_config_4 <- setup_experiment(
  experiment_prefix = "hum_tran_2024_03_deep_phenotype_muscle",
  experiment_date = "2024-03",
  data_source = "pluto",
  data_type = "transcriptomics",
  species = "human",
  metadata_columns = c("sample_id", "group", "subject_status", "sex"),
  additional_notes = "GEO Dataset: Deep phenotyping of Post-infectious Myalgic Encephalomyelitis/Chronic Fatigue Syndrome [RNA-Seq]"
)

hum_tran_2024_03_deep_phenotype_muscle <- import_omics_data_with_config(
  config = experiment_config_4,
  experiment_id = "PLX124291",
  api_key = Sys.getenv("PLUTO_API_KEY"),
  description = "ME/CFS GEO Dataset: Deep phenotyping of Post-infectious Myalgic Encephalomyelitis/Chronic Fatigue Syndrome [RNA-Seq]",
  experimenter = "Shreya Sharma"
)

# Access the created object
hum_tran_2024_03_deep_phenotype_muscle$metadata$metadata_columns
hum_tran_2024_03_deep_phenotype_muscle$metadata$conditions$sample_id
hum_tran_2024_03_deep_phenotype_muscle$metadata$conditions$group
hum_tran_2024_03_deep_phenotype_muscle$metadata$conditions$subject_status  
hum_tran_2024_03_deep_phenotype_muscle$metadata$conditions$sex
hum_tran_2024_03_deep_phenotype_muscle$metadata
hum_tran_2024_03_deep_phenotype_muscle$data_dictionary
hum_tran_2024_03_deep_phenotype_muscle$sample_data
head(hum_tran_2024_03_deep_phenotype_muscle$assay_data[[1]])

#Human PBMCs Experiment (GEO Dataset)

experiment_config_5 <- setup_experiment(
  experiment_prefix = "hum_tran_2024_03_deep_phenotype_pbmcs",
  experiment_date = "2024-03",
  data_source = "pluto",
  data_type = "transcriptomics",
  species = "human",
  metadata_columns = c("sample_id", "group", "subject_status", "sex"),
  additional_notes = "GEO Dataset: Deep phenotyping of Post-infectious Myalgic Encephalomyelitis/Chronic Fatigue Syndrome [PBMC RNA-Seq]"
)

hum_tran_2024_03_deep_phenotype_pbmcs <- import_omics_data_with_config(
  config = experiment_config_5,
  experiment_id = "PLX073352",
  api_key = Sys.getenv("PLUTO_API_KEY"),
  description = "ME/CFS GEO Dataset: Deep phenotyping of Post-infectious Myalgic Encephalomyelitis/Chronic Fatigue Syndrome [PBMC RNA-Seq]",
  experimenter = "Shreya Sharma"
)

# Access the created object
hum_tran_2024_03_deep_phenotype_pbmcs$metadata$metadata_columns
hum_tran_2024_03_deep_phenotype_pbmcs$metadata$conditions$sample_id
hum_tran_2024_03_deep_phenotype_pbmcs$metadata$conditions$group
hum_tran_2024_03_deep_phenotype_pbmcs$metadata$conditions$subject_status  
hum_tran_2024_03_deep_phenotype_pbmcs$metadata$conditions$sex
hum_tran_2024_03_deep_phenotype_pbmcs$metadata
hum_tran_2024_03_deep_phenotype_pbmcs$data_dictionary
hum_tran_2024_03_deep_phenotype_pbmcs$sample_data
head(hum_tran_2024_03_deep_phenotype_pbmcs$assay_data[[1]])


# Setup first experiment
experiment_config_6 <- setup_experiment(
  experiment_prefix = "xen_tran_2025_02",
  experiment_date = "2025-02",
  data_source = "pluto",
  data_type = "transcriptomics",
  species = "xenopus",
  metadata_columns = c("condition",	"collection_day",	"collection_time",	"timepoint",	"tube_label",	"drug",	"dose",	"replicate"),
  additional_notes = "Confounding variables experiment"
)

xen_tran_2025_02 <- import_omics_data_with_config(
  config = experiment_config_6,
  experiment_id = "PLX166849",
  api_key = Sys.getenv("PLUTO_API_KEY"),
  save_credentials = TRUE,
  description = "Xenopus transcriptomics experiment February 2025",
  experimenter = "Allison Grossberg"
)

# Access metadata columns
xen_tran_2025_02$metadata$metadata_columns
xen_tran_2025_02$metadata
xen_tran_2025_02$data_dictionary
xen_tran_2025_02$sample_data
head(xen_tran_2025_02$assay_data[[1]])

