#' Package Setup and Dependencies
#' 
#' Functions for installing and managing package dependencies.

# Configuration
setup_config <- list(
    cran_packages = c(
        # Core data manipulation and visualization
        "dplyr", "tidyr", "flextable", "officer", "ggplot2", "ggsignif", "rstatix",
        # Network analysis and visualization
        "ggraph", "tidygraph", "igraph", "viridis",
        # Additional utilities
        "purrr", "pluto", "tibble", "magrittr", "ggrepel",
        # Heatmap visualization
        "circlize", "grid", "RColorBrewer",
        # UpSet Plots
        "UpSetR",
        # R Environment set-up
        "jsonlite", "cli", "logger", "here"
    ),
    bioc_packages = c("DESeq2", "ComplexHeatmap", "limma", "TPP", "MSnbase", 
                     "msigdbr", "fgsea", "STRINGdb", "org.Hs.eg.db", "AnnotationDbi")
)

#' Install single package
#' @param package_name Package name
#' @param source Package source
#' @return Success status
install_single_package <- function(package_name, source = "cran") {
    # Validate inputs
    if (is.null(package_name)) {
        stop("package_name cannot be NULL")
    }
    if (!is.character(package_name) || nchar(package_name) == 0) {
        stop("package_name must be a non-empty character string")
    }
    if (!source %in% c("cran", "bioc")) {
        stop("Invalid source. Must be 'cran' or 'bioc'")
    }
    
    tryCatch({
        if (source == "cran") {
            if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
                message("Installing CRAN package: ", package_name)
                install.packages(package_name, dependencies = TRUE)
                return(require(package_name, character.only = TRUE, quietly = TRUE))
            }
        } else if (source == "bioc") {
            if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
                message("Installing Bioconductor package: ", package_name)
                BiocManager::install(package_name, dependencies = TRUE)
                return(require(package_name, character.only = TRUE, quietly = TRUE))
            }
        }
        return(TRUE)
    }, error = function(e) {
        warning("Failed to install ", package_name, ": ", e$message)
        return(FALSE)
    })
}

#' Install BiocManager if not available
#' 
#' @return Logical indicating success
#' @export
install_bioc_manager <- function() {
    tryCatch({
        if (!require("BiocManager", quietly = TRUE)) {
            message("Installing BiocManager...")
            install.packages("BiocManager", dependencies = TRUE)
        }
        return(require("BiocManager", quietly = TRUE))
    }, error = function(e) {
        stop("Failed to install BiocManager: ", e$message)
    })
}

#' Main setup function
#' 
#' @param config List containing package configuration
#' @return List with installation results
#' @export
setup_packages <- function(config = setup_config) {
    message("Starting package setup...")
    
    # Validate configuration
    if (!is.list(config) || !all(c("cran_packages", "bioc_packages") %in% names(config))) {
        stop("Invalid configuration provided")
    }
    
    results <- list(
        cran_success = character(0),
        cran_failed = character(0),
        bioc_success = character(0),
        bioc_failed = character(0)
    )
    
    # Install BiocManager first
    if (!install_bioc_manager()) {
        stop("Cannot proceed without BiocManager")
    }
    
    # Install CRAN packages
    message("Installing CRAN packages...")
    for (pkg in config$cran_packages) {
        if (install_single_package(pkg, "cran")) {
            results$cran_success <- c(results$cran_success, pkg)
        } else {
            results$cran_failed <- c(results$cran_failed, pkg)
        }
    }
    
    # Install Bioconductor packages
    message("Installing Bioconductor packages...")
    for (pkg in config$bioc_packages) {
        if (install_single_package(pkg, "bioc")) {
            results$bioc_success <- c(results$bioc_success, pkg)
        } else {
            results$bioc_failed <- c(results$bioc_failed, pkg)
        }
    }
    
    # Print summary
    message("\n=== Installation Summary ===")
    message("CRAN packages - Success: ", length(results$cran_success), 
            ", Failed: ", length(results$cran_failed))
    message("Bioconductor packages - Success: ", length(results$bioc_success), 
            ", Failed: ", length(results$bioc_failed))
    
    if (length(results$cran_failed) > 0) {
        warning("Failed CRAN packages: ", paste(results$cran_failed, collapse = ", "))
    }
    if (length(results$bioc_failed) > 0) {
        warning("Failed Bioconductor packages: ", paste(results$bioc_failed, collapse = ", "))
    }
    
    return(results)
}

#' Create a reusable R environment for transcriptomics analysis
#' 
#' This function sets up a complete R environment with all required packages
#' and creates necessary directory structures for transcriptomics analysis.
#' 
#' @param install_packages Whether to install packages (default: TRUE)
#' @param create_directories Whether to create standard directory structure (default: TRUE)
#' @param base_dir Base directory for the project (default: current working directory)
#' @return List with setup results and environment information
#' @export
create_transcriptomics_environment <- function(install_packages = TRUE, 
                                             create_directories = TRUE,
                                             base_dir = getwd()) {
    message("Setting up transcriptomics analysis environment...")
    
    # Initialize results list
    setup_results <- list(
        packages_installed = FALSE,
        directories_created = FALSE,
        environment_info = list(),
        errors = character(0)
    )
    
    # 1. Install packages if requested
    if (install_packages) {
        message("Installing required packages...")
        tryCatch({
            package_results <- setup_packages()
            setup_results$packages_installed <- TRUE
            setup_results$package_results <- package_results
            message("✓ Package installation completed")
        }, error = function(e) {
            setup_results$errors <- c(setup_results$errors, paste("Package installation failed:", e$message))
            message("✗ Package installation failed: ", e$message)
        })
    }
    
    # 2. Create directory structure if requested
    if (create_directories) {
        message("Creating directory structure...")
        tryCatch({
            create_standard_directories(base_dir)
            setup_results$directories_created <- TRUE
            message("✓ Directory structure created")
        }, error = function(e) {
            setup_results$errors <- c(setup_results$errors, paste("Directory creation failed:", e$message))
            message("✗ Directory creation failed: ", e$message)
        })
    }
    
    # 3. Set up environment information
    setup_results$environment_info <- list(
        r_version = R.version.string,
        platform = R.version$platform,
        working_directory = getwd(),
        base_directory = base_dir,
        setup_date = Sys.time(),
        loaded_packages = .packages()
    )
    
    # 4. Load the transcriptomics pipeline
    message("Loading transcriptomics pipeline...")
    tryCatch({
        if (requireNamespace("devtools", quietly = TRUE)) {
            devtools::load_all()
            setup_results$pipeline_loaded <- TRUE
            message("✓ Transcriptomics pipeline loaded successfully")
        } else {
            setup_results$errors <- c(setup_results$errors, "devtools not available for loading pipeline")
            message("✗ devtools not available")
        }
    }, error = function(e) {
        setup_results$errors <- c(setup_results$errors, paste("Pipeline loading failed:", e$message))
        message("✗ Pipeline loading failed: ", e$message)
    })
    
    # 5. Print summary
    message("\n=== Environment Setup Summary ===")
    message("Packages installed: ", setup_results$packages_installed)
    message("Directories created: ", setup_results$directories_created)
    message("Pipeline loaded: ", setup_results$pipeline_loaded %||% FALSE)
    
    if (length(setup_results$errors) > 0) {
        message("\nErrors encountered:")
        for (error in setup_results$errors) {
            message("- ", error)
        }
    } else {
        message("\n✓ Environment setup completed successfully!")
        message("You can now run transcriptomics analyses using the pipeline functions.")
    }
    
    return(setup_results)
}

#' Create standard directory structure for transcriptomics analysis
#' 
#' @param base_dir Base directory for the project
#' @return List of created directories
create_standard_directories <- function(base_dir = getwd()) {
    # Define standard directory structure
    directories <- c(
        "data",
        "results", 
        "input_data",
        "input_files",
        "scripts",
        "reports",
        "temp"
    )
    
    created_dirs <- character(0)
    
    for (dir in directories) {
        dir_path <- file.path(base_dir, dir)
        if (!dir.exists(dir_path)) {
            if (dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)) {
                created_dirs <- c(created_dirs, dir_path)
                message("Created directory: ", dir_path)
            } else {
                warning("Failed to create directory: ", dir_path)
            }
        } else {
            message("Directory already exists: ", dir_path)
        }
    }
    
    return(created_dirs)
}

#' Check if the transcriptomics environment is properly set up
#' 
#' @return List with environment status
#' @export
check_environment <- function() {
    message("Checking transcriptomics environment...")
    
    status <- list(
        packages_available = FALSE,
        directories_exist = FALSE,
        pipeline_loaded = FALSE,
        issues = character(0)
    )
    
    # Check if required packages are available
    required_packages <- c("DESeq2", "dplyr", "ggplot2", "here")
    missing_packages <- character(0)
    
    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            missing_packages <- c(missing_packages, pkg)
        }
    }
    
    if (length(missing_packages) == 0) {
        status$packages_available <- TRUE
        message("✓ All required packages are available")
    } else {
        status$issues <- c(status$issues, paste("Missing packages:", paste(missing_packages, collapse = ", ")))
        message("✗ Missing packages: ", paste(missing_packages, collapse = ", "))
    }
    
    # Check if standard directories exist
    required_dirs <- c("data", "results", "input_data", "input_files")
    missing_dirs <- character(0)
    
    for (dir in required_dirs) {
        if (!dir.exists(dir)) {
            missing_dirs <- c(missing_dirs, dir)
        }
    }
    
    if (length(missing_dirs) == 0) {
        status$directories_exist <- TRUE
        message("✓ Standard directories exist")
    } else {
        status$issues <- c(status$issues, paste("Missing directories:", paste(missing_dirs, collapse = ", ")))
        message("✗ Missing directories: ", paste(missing_dirs, collapse = ", "))
    }
    
    # Check if pipeline functions are available
    pipeline_functions <- c("import_omics_data", "run_deseq2_analysis", "perform_pca_analysis")
    missing_functions <- character(0)
    
    for (func in pipeline_functions) {
        if (!exists(func)) {
            missing_functions <- c(missing_functions, func)
        }
    }
    
    if (length(missing_functions) == 0) {
        status$pipeline_loaded <- TRUE
        message("✓ Pipeline functions are available")
    } else {
        status$issues <- c(status$issues, paste("Missing functions:", paste(missing_functions, collapse = ", ")))
        message("✗ Missing functions: ", paste(missing_functions, collapse = ", "))
    }
    
    # Overall status
    if (length(status$issues) == 0) {
        message("\n✓ Environment is properly set up!")
    } else {
        message("\n✗ Environment has issues:")
        for (issue in status$issues) {
            message("- ", issue)
        }
    }
    
    return(status)
}

