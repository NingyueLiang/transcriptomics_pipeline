#' Setup and Install Required Packages
#' 
#' This module provides functions for setting up and installing all required packages
#' for the transcriptomics pipeline, including both CRAN and Bioconductor packages
#' with proper error handling.

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

#' Install a single package with error handling
#' 
#' @param package_name Name of the package to install
#' @param source Source of the package ("cran" or "bioc")
#' @return Logical indicating success
#' @export
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

