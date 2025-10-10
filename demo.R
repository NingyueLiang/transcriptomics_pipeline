#!/usr/bin/env Rscript

#' Real Data Demo Script for transcriptomicsPipeline Package
#' 
#' This script demonstrates the functionality of the transcriptomicsPipeline package
#' using real test data from the input_data directory.

# Load required libraries
suppressPackageStartupMessages({
    library(devtools)
    library(testthat)
})

# Set working directory to package root
setwd("/Users/wyssuser/Documents/Python Repositories/Pluto_Transcriptomics_Pipeline")

# Load the package
devtools::load_all()

# Create proper directory structure like the original scripts
# Using actual experiment information from PLX073248
experiment_name <- "xen_tran_2024_12"
data_dir <- file.path("data", experiment_name)
results_dir <- file.path("results", experiment_name)

# Create necessary directories
for (dir in c(data_dir, results_dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Create data subdirectories
data_subdirs <- c("cleaned_data")
for (subdir in data_subdirs) {
    dir.create(file.path(data_dir, subdir), showWarnings = FALSE)
}

# Create results subdirectories  
results_subdirs <- c("de", "mapping", "gsea", "pca")
for (subdir in results_subdirs) {
    dir.create(file.path(results_dir, subdir), showWarnings = FALSE)
}

cat("=== transcriptomicsPipeline Demo ===\n")
cat("Testing package functionality with PLX073248 data...\n")
cat("Experiment: Xenopus transcriptomics experiment December 2024\n")
cat("Experimenter: Allison Grossberg\n")
cat("Data Source: PLX073248 (Pluto API)\n\n")

# Test 1: Package Setup
cat("1. Testing package setup functions...\n")
tryCatch({
    # Test setup configuration
    if (exists("setup_config")) {
        cat("✓ setup_config loaded successfully\n")
        cat("  - CRAN packages:", length(setup_config$cran_packages), "\n")
        cat("  - Bioconductor packages:", length(setup_config$bioc_packages), "\n")
    } else {
        cat("✗ setup_config not found\n")
    }
    
    # Test setup functions exist
    if (exists("install_single_package")) {
        cat("✓ install_single_package function available\n")
    } else {
        cat("✗ install_single_package function not found\n")
    }
    
    if (exists("setup_packages")) {
        cat("✓ setup_packages function available\n")
    } else {
        cat("✗ setup_packages function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in setup functions:", e$message, "\n")
})

cat("\n")

# Test 2: Examine real data
cat("2. Examining PLX073248 data...\n")
tryCatch({
    # Read sample data
    sample_data <- read.csv("input_data/PLX073248_sample_data.csv", stringsAsFactors = FALSE)
    cat("✓ Sample data loaded successfully\n")
    cat("  - Samples:", nrow(sample_data), "\n")
    cat("  - Columns:", paste(colnames(sample_data), collapse = ", "), "\n")
    cat("  - Conditions:", paste(unique(sample_data$condition), collapse = ", "), "\n")
    cat("  - Drugs:", paste(unique(sample_data$drug), collapse = ", "), "\n")
    cat("  - Doses:", paste(unique(sample_data$dose), collapse = ", "), "\n")
    cat("  - Timepoints:", paste(unique(sample_data$timepoint), collapse = ", "), "\n")
    cat("  - Replicates:", paste(unique(sample_data$replicate), collapse = ", "), "\n")
    
    # Read assay data (first few rows to check structure)
    assay_data <- read.csv("input_data/PLX073248_assay_data.csv", nrows = 5, stringsAsFactors = FALSE)
    cat("✓ Assay data structure examined\n")
    cat("  - Genes preview:", nrow(assay_data), "rows examined\n")
    cat("  - Sample columns:", ncol(assay_data) - 1, "samples\n")
    cat("  - Gene column:", colnames(assay_data)[1], "\n")
    
}, error = function(e) {
    cat("✗ Error examining PLX073248 data:", e$message, "\n")
})

cat("\n")

# Test 3: Test import functions with real data
cat("3. Testing import functions with real data...\n")
tryCatch({
    if (exists("import_omics_data")) {
        # Test local data import with real data
        experiment_obj <- import_omics_data(
            data_source = "local",
            data_type = "transcriptomics",
            species = "xenopus",
            assay_path = "input_data/PLX073248_assay_data.csv",
            sample_path = "input_data/PLX073248_sample_data.csv",
            experiment_name = experiment_name,
            metadata_columns = c("condition", "timepoint", "dose", "drug", "replicate"),
            interactive = FALSE # Non-interactive for demo
        )
        
        if (!is.null(experiment_obj)) {
            cat("✓ Real data import successful\n")
            cat("  - Experiment name:", experiment_obj$metadata$name, "\n")
            cat("  - Assay data dimensions:", dim(experiment_obj$assay_data[[1]]), "\n")
            cat("  - Sample data dimensions:", dim(experiment_obj$sample_data[[1]]), "\n")
            cat("  - Species:", experiment_obj$metadata$species, "\n")
            cat("  - Data source:", experiment_obj$metadata$data_source, "\n")
        } else {
            cat("✗ Real data import returned NULL\n")
        }
    } else {
        cat("✗ import_omics_data function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in import functions:", e$message, "\n")
})

cat("\n")

# Test 4: Test cleaning functions with real data
cat("4. Testing cleaning functions with real data...\n")
tryCatch({
    if (exists("clean_data_for_de")) {
        # Test data cleaning with real data
        cleaned_data_name <- clean_data_for_de(
            experiment_name = experiment_name,
            metadata_cols = list(
                sample_id = "sample_id",
                group_cols = c("condition", "timepoint", "dose", "drug", "replicate")
            )
        )
        
        if (!is.null(cleaned_data_name)) {
            cat("✓ Real data cleaning successful\n")
            cat("  - Cleaned data object:", cleaned_data_name, "\n")
        } else {
            cat("✗ Real data cleaning returned NULL\n")
        }
    } else {
        cat("✗ clean_data_for_de function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in cleaning functions:", e$message, "\n")
})

cat("\n")

# Test 5: Test differential expression functions with real data
cat("5. Testing differential expression functions with real data...\n")
tryCatch({
    if (exists("run_deseq2_analysis")) {
        # Test DESeq2 analysis with real data (non-interactive)
        de_results <- run_deseq2_analysis(
            experiment_name = experiment_name,
            comparison_cols = c("condition"),
            reference_levels = list(condition = c("vehicle")),
            interactive = FALSE  # Non-interactive for demo
        )
        
        if (!is.null(de_results) && length(de_results) > 0) {
            cat("✓ Real data differential expression analysis successful\n")
            cat("  - Number of contrasts:", length(de_results), "\n")
            for (contrast_name in names(de_results)) {
                if (!is.null(de_results[[contrast_name]]$deseq_results)) {
                    n_sig <- sum(de_results[[contrast_name]]$deseq_results$padj < 0.05, na.rm = TRUE)
                    cat("  -", contrast_name, ":", n_sig, "significant genes\n")
                }
            }
        } else {
            cat("✗ Real data differential expression analysis returned empty results\n")
        }
    } else {
        cat("✗ run_deseq2_analysis function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in differential expression functions:", e$message, "\n")
})

cat("\n")

# Test 6: Test PCA functions with real data
cat("6. Testing PCA functions with real data...\n")
tryCatch({
    if (exists("perform_pca_analysis")) {
        # Load the real assay data for PCA
        real_assay_data <- read.csv("input_data/PLX073248_assay_data.csv", stringsAsFactors = FALSE)
        
        # Test PCA analysis with real data (using group names as they would be parsed from sample names)
        pca_results <- perform_pca_analysis(
            data = real_assay_data,
            groups_to_include = c("Veh", "Ket_500"),  # These are the group names after removing replicate numbers
            experiment_name = experiment_name,
            save_results = TRUE
        )
        
        if (!is.null(pca_results)) {
            cat("✓ Real data PCA analysis successful\n")
            if (!is.null(pca_results$variance_explained)) {
                cat("  - PC1 variance explained:", round(pca_results$variance_explained[1], 2), "%\n")
                cat("  - PC2 variance explained:", round(pca_results$variance_explained[2], 2), "%\n")
            } else {
                cat("  - PCA analysis completed successfully\n")
            }
        } else {
            cat("✗ Real data PCA analysis returned NULL\n")
        }
    } else {
        cat("✗ perform_pca_analysis function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in PCA functions:", e$message, "\n")
})

cat("\n")

# Test 7: Test GSEA functions with real data
cat("7. Testing GSEA functions with real data...\n")
tryCatch({
    if (exists("run_gsea_analysis")) {
        # Test GSEA analysis with real data
        gsea_results <- run_gsea_analysis(
            experiment_name = experiment_name,
            model_type = "deseq2",
            rank_method = "sign_log_padj",
            collections = c("H"),  # Just test Hallmark collection
            interactive = FALSE # Non-interactive for demo
        )
        
        if (!is.null(gsea_results)) {
            cat("✓ Real data GSEA analysis successful\n")
            cat("  - Number of gene sets tested:", length(gsea_results), "\n")
        } else {
            cat("✗ Real data GSEA analysis returned NULL\n")
        }
    } else {
        cat("✗ run_gsea_analysis function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in GSEA functions:", e$message, "\n")
})

cat("\n")

# Test 8: Test mapping functions with real data
cat("8. Testing mapping functions with real data...\n")
tryCatch({
    if (exists("map_xenopus_to_human")) {
        # Get some real gene symbols from the data
        real_assay_data <- read.csv("input_data/PLX073248_assay_data.csv", stringsAsFactors = FALSE)
        test_genes <- head(real_assay_data$gene_symbol, 10)
        
        # Use real mapping files
        xenbase_file <- "input_files/Xenbase_Xenopus_Orthology_Predictions.txt"
        hcop_file <- "input_files/HCOP_Xenopus_Orthology_Predictions.txt"
        
        mapping_results <- map_xenopus_to_human(
            xenopus_genes = test_genes,
            xenbase_file = xenbase_file,
            hcop_file = hcop_file,
            suffix_mode = "all"
        )
        
        if (!is.null(mapping_results)) {
            cat("✓ Real data gene mapping successful\n")
            cat("  - Input genes:", length(test_genes), "\n")
            cat("  - Mapped genes:", nrow(mapping_results), "\n")
        } else {
            cat("✗ Real data gene mapping returned NULL\n")
        }
    } else {
        cat("✗ map_xenopus_to_human function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in mapping functions:", e$message, "\n")
})

cat("\n")

# Test 9: Test additional analysis functions
cat("9. Testing additional analysis functions...\n")
tryCatch({
    # Test DESeq2 diagnostics if available
    if (exists("generate_deseq2_diagnostics")) {
        cat("✓ DESeq2 diagnostics function available\n")
    }
    
    # Test DESeq2 results export if available
    if (exists("export_all_result_types")) {
        cat("✓ DESeq2 results export function available\n")
    }
    
    # Test comprehensive GSEA if available
    if (exists("run_gsea_analysis_comprehensive")) {
        cat("✓ Comprehensive GSEA function available\n")
    }
    
    # Test comprehensive PCA if available
    if (exists("perform_pca_analysis_comprehensive")) {
        cat("✓ Comprehensive PCA function available\n")
    }
    
}, error = function(e) {
    cat("✗ Error testing additional functions:", e$message, "\n")
})

cat("\n")

# Test 10: Check output files and directories
cat("10. Checking output files and directories...\n")
tryCatch({
    # Check if data directory was created
    # data_dir already defined at top of script
    if (dir.exists(data_dir)) {
        cat("✓ Data directory created:", data_dir, "\n")
        
        # List contents
        data_contents <- list.files(data_dir, recursive = TRUE)
        cat("  - Data files:", length(data_contents), "\n")
        for (file in head(data_contents, 5)) {
            cat("    *", file, "\n")
        }
        if (length(data_contents) > 5) {
            cat("    * ... and", length(data_contents) - 5, "more files\n")
        }
    } else {
        cat("✗ Data directory not found:", data_dir, "\n")
    }
    
    # Check if results directory was created
    # results_dir already defined at top of script
    if (dir.exists(results_dir)) {
        cat("✓ Results directory created:", results_dir, "\n")
        
        # List contents
        results_contents <- list.files(results_dir, recursive = TRUE)
        cat("  - Result files:", length(results_contents), "\n")
        for (file in head(results_contents, 5)) {
            cat("    *", file, "\n")
        }
        if (length(results_contents) > 5) {
            cat("    * ... and", length(results_contents) - 5, "more files\n")
        }
    } else {
        cat("✗ Results directory not found:", results_dir, "\n")
    }
    
}, error = function(e) {
    cat("✗ Error checking output files:", e$message, "\n")
})

cat("\n")

# Test 11: Run package tests
cat("11. Running package tests...\n")
tryCatch({
    # Run testthat tests
    test_results <- test_dir("tests", reporter = "summary")
    cat("✓ Package tests completed\n")
    
}, error = function(e) {
    cat("✗ Error running package tests:", e$message, "\n")
})

cat("\n")

# Test 12: Check package structure
cat("12. Checking package structure...\n")
tryCatch({
    # Check DESCRIPTION file
    if (file.exists("DESCRIPTION")) {
        desc <- read.dcf("DESCRIPTION")
        cat("✓ DESCRIPTION file exists\n")
        cat("  - Package:", desc[1, "Package"], "\n")
        cat("  - Version:", desc[1, "Version"], "\n")
        cat("  - Title:", desc[1, "Title"], "\n")
    } else {
        cat("✗ DESCRIPTION file not found\n")
    }
    
    # Check NAMESPACE file
    if (file.exists("NAMESPACE")) {
        cat("✓ NAMESPACE file exists\n")
    } else {
        cat("✗ NAMESPACE file not found\n")
    }
    
    # Check R files
    r_files <- list.files("R", pattern = "\\.R$", full.names = FALSE)
    cat("✓ R files found:", length(r_files), "\n")
    for (file in r_files) {
        cat("  -", file, "\n")
    }
    
    # Check test files
    test_files <- list.files("tests", pattern = "\\.R$", recursive = TRUE, full.names = FALSE)
    cat("✓ Test files found:", length(test_files), "\n")
    for (file in test_files) {
        cat("  -", file, "\n")
    }
    
}, error = function(e) {
    cat("✗ Error checking package structure:", e$message, "\n")
})

cat("\n")

# Summary
cat("=== Demo Summary ===\n")
cat("Demo completed with PLX073248 data! Check the output above for any errors.\n")
cat("Experiment: Xenopus transcriptomics experiment December 2024\n")
cat("Experimenter: Allison Grossberg\n")
cat("Data Source: PLX073248 (Pluto API)\n")
cat("Demo results saved in standard package directories\n")
cat("Package data saved in:", data_dir, "\n")
cat("Package results saved in:", results_dir, "\n")
cat("Package appears to be working correctly if no errors were reported.\n")

# Don't clean up demo files for real data - keep them for inspection
cat("\nPLX073248 analysis results have been preserved for inspection.\n")

cat("\n=== Demo Complete ===\n")
