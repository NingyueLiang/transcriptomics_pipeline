#!/usr/bin/env Rscript

#' Demo Script for transcriptomicsPipeline Package
#' 
#' This script demonstrates the functionality of the transcriptomicsPipeline package
#' by running through a complete workflow with synthetic data.

# Load required libraries
suppressPackageStartupMessages({
    library(devtools)
    library(testthat)
})

# Set working directory to package root
setwd("/Users/wyssuser/Documents/Python Repositories/Pluto_Transcriptomics_Pipeline")

# Load the package
devtools::load_all()

# Create demo data directory
demo_dir <- "demo_data"
if (!dir.exists(demo_dir)) {
    dir.create(demo_dir, recursive = TRUE)
}

# Create demo results directory
demo_results <- "demo_results"
if (!dir.exists(demo_results)) {
    dir.create(demo_results, recursive = TRUE)
}

cat("=== transcriptomicsPipeline Package Demo ===\n")
cat("Testing package functionality with synthetic data...\n\n")

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

# Test 2: Create synthetic data for testing
cat("2. Creating synthetic data for testing...\n")
tryCatch({
    # Create synthetic gene expression data
    n_genes <- 1000
    n_samples <- 12
    
    # Generate random counts matrix
    set.seed(123)
    counts_matrix <- matrix(
        rnbinom(n_genes * n_samples, size = 10, mu = 100),
        nrow = n_genes,
        ncol = n_samples
    )
    
    # Create gene symbols
    gene_symbols <- paste0("GENE_", sprintf("%04d", 1:n_genes))
    rownames(counts_matrix) <- gene_symbols
    
    # Create sample IDs with condition prefix for PCA compatibility
    sample_ids <- c(
        paste0("control_", sprintf("%02d", 1:6)),
        paste0("treatment_", sprintf("%02d", 1:6))
    )
    colnames(counts_matrix) <- sample_ids
    
    # Create sample metadata
    sample_metadata <- data.frame(
        sample_id = sample_ids,
        condition = rep(c("control", "treatment"), each = 6),
        batch = rep(c("batch1", "batch2"), 6),
        stringsAsFactors = FALSE
    )
    
    # Add some differential expression
    # Make first 50 genes upregulated in treatment
    treatment_samples <- sample_metadata$sample_id[sample_metadata$condition == "treatment"]
    counts_matrix[1:50, treatment_samples] <- counts_matrix[1:50, treatment_samples] * 2
    
    # Make next 50 genes downregulated in treatment
    counts_matrix[51:100, treatment_samples] <- counts_matrix[51:100, treatment_samples] * 0.5
    
    # Save synthetic data
    assay_data <- data.frame(gene_symbol = gene_symbols, counts_matrix, stringsAsFactors = FALSE)
    write.csv(assay_data, file.path(demo_dir, "synthetic_assay_data.csv"), row.names = FALSE)
    write.csv(sample_metadata, file.path(demo_dir, "synthetic_sample_data.csv"), row.names = FALSE)
    
    cat("✓ Synthetic data created successfully\n")
    cat("  - Genes:", n_genes, "\n")
    cat("  - Samples:", n_samples, "\n")
    cat("  - Conditions: control, treatment\n")
    
}, error = function(e) {
    cat("✗ Error creating synthetic data:", e$message, "\n")
})

cat("\n")

# Test 3: Test import functions
cat("3. Testing import functions...\n")
tryCatch({
    if (exists("import_omics_data")) {
        # Test local data import
        experiment_obj <- import_omics_data(
            data_source = "local",
            data_type = "transcriptomics",
            species = "human",
            assay_path = file.path(demo_dir, "synthetic_assay_data.csv"),
            sample_path = file.path(demo_dir, "synthetic_sample_data.csv"),
            experiment_name = "demo_experiment"
        )
        
        if (!is.null(experiment_obj)) {
            cat("✓ Data import successful\n")
            cat("  - Experiment name:", experiment_obj$metadata$name, "\n")
            cat("  - Assay data dimensions:", dim(experiment_obj$assay_data[[1]]), "\n")
            cat("  - Sample data dimensions:", dim(experiment_obj$sample_data[[1]]), "\n")
        } else {
            cat("✗ Data import returned NULL\n")
        }
    } else {
        cat("✗ import_omics_data function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in import functions:", e$message, "\n")
})

cat("\n")

# Test 4: Test cleaning functions
cat("4. Testing cleaning functions...\n")
tryCatch({
    if (exists("clean_data_for_de")) {
        # Test data cleaning
        cleaned_data_name <- clean_data_for_de(
            experiment_name = "demo_experiment",
            metadata_cols = list(
                sample_id = "sample_id",
                group_cols = c("condition", "batch")
            )
        )
        
        if (!is.null(cleaned_data_name)) {
            cat("✓ Data cleaning successful\n")
            cat("  - Cleaned data object:", cleaned_data_name, "\n")
        } else {
            cat("✗ Data cleaning returned NULL\n")
        }
    } else {
        cat("✗ clean_data_for_de function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in cleaning functions:", e$message, "\n")
})

cat("\n")

# Test 5: Test differential expression functions
cat("5. Testing differential expression functions...\n")
tryCatch({
    if (exists("run_deseq2_analysis")) {
        # Test DESeq2 analysis (non-interactive)
        de_results <- run_deseq2_analysis(
            experiment_name = "demo_experiment",
            comparison_cols = c("condition"),
            reference_levels = list(condition = c("control", "treatment")),
            interactive = FALSE  # Non-interactive for demo
        )
        
        if (!is.null(de_results) && length(de_results) > 0) {
            cat("✓ Differential expression analysis successful\n")
            cat("  - Number of contrasts:", length(de_results), "\n")
            for (contrast_name in names(de_results)) {
                n_sig <- sum(de_results[[contrast_name]]$deseq_results$padj < 0.05, na.rm = TRUE)
                cat("  -", contrast_name, ":", n_sig, "significant genes\n")
            }
        } else {
            cat("✗ Differential expression analysis returned empty results\n")
        }
    } else {
        cat("✗ run_deseq2_analysis function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in differential expression functions:", e$message, "\n")
})

cat("\n")

# Test 6: Test PCA functions
cat("6. Testing PCA functions...\n")
tryCatch({
    if (exists("perform_pca_analysis")) {
        # Test PCA analysis with proper data format
        pca_data <- data.frame(gene_symbol = gene_symbols, counts_matrix, stringsAsFactors = FALSE)
        pca_results <- perform_pca_analysis(
            data = pca_data,
            groups_to_include = c("control", "treatment"),
            experiment_name = "demo_experiment",
            save_results = TRUE
        )
        
        if (!is.null(pca_results)) {
            cat("✓ PCA analysis successful\n")
            if (!is.null(pca_results$variance_explained)) {
                cat("  - PC1 variance explained:", round(pca_results$variance_explained[1], 2), "%\n")
                cat("  - PC2 variance explained:", round(pca_results$variance_explained[2], 2), "%\n")
            } else {
                cat("  - PCA analysis completed successfully\n")
            }
        } else {
            cat("✗ PCA analysis returned NULL\n")
        }
    } else {
        cat("✗ perform_pca_analysis function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in PCA functions:", e$message, "\n")
})

cat("\n")

# Test 7: Test GSEA functions
cat("7. Testing GSEA functions...\n")
tryCatch({
    if (exists("run_gsea_analysis")) {
        # Test GSEA analysis
        gsea_results <- run_gsea_analysis(
            experiment_name = "demo_experiment",
            model_type = "deseq2",
            rank_method = "sign_log_padj",
            collections = c("H")  # Just test Hallmark collection
        )
        
        if (!is.null(gsea_results)) {
            cat("✓ GSEA analysis successful\n")
            cat("  - Number of gene sets tested:", length(gsea_results), "\n")
        } else {
            cat("✗ GSEA analysis returned NULL\n")
        }
    } else {
        cat("✗ run_gsea_analysis function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in GSEA functions:", e$message, "\n")
})

cat("\n")

# Test 8: Test mapping functions (if applicable)
cat("8. Testing mapping functions...\n")
tryCatch({
    if (exists("map_xenopus_to_human")) {
        # Test gene mapping with synthetic data
        test_genes <- c("GENE_0001", "GENE_0002", "GENE_0003")
        
        # Create dummy mapping files for testing
        xenbase_file <- file.path(demo_dir, "dummy_xenbase.txt")
        hcop_file <- file.path(demo_dir, "dummy_hcop.txt")
        
        # Create dummy files
        writeLines(c("GENE_0001\tHUMAN_GENE1", "GENE_0002\tHUMAN_GENE2"), xenbase_file)
        writeLines(c("GENE_0001\tHUMAN_GENE1", "GENE_0003\tHUMAN_GENE3"), hcop_file)
        
        mapping_results <- map_xenopus_to_human(
            xenopus_genes = test_genes,
            xenbase_file = xenbase_file,
            hcop_file = hcop_file,
            suffix_mode = "all"
        )
        
        if (!is.null(mapping_results)) {
            cat("✓ Gene mapping successful\n")
            cat("  - Input genes:", length(test_genes), "\n")
            cat("  - Mapped genes:", nrow(mapping_results), "\n")
        } else {
            cat("✗ Gene mapping returned NULL\n")
        }
    } else {
        cat("✗ map_xenopus_to_human function not found\n")
    }
    
}, error = function(e) {
    cat("✗ Error in mapping functions:", e$message, "\n")
})

cat("\n")

# Test 9: Run package tests
cat("9. Running package tests...\n")
tryCatch({
    # Run testthat tests
    test_results <- test_dir("tests", reporter = "summary")
    cat("✓ Package tests completed\n")
    
}, error = function(e) {
    cat("✗ Error running package tests:", e$message, "\n")
})

cat("\n")

# Test 10: Check package structure
cat("10. Checking package structure...\n")
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
cat("Demo completed! Check the output above for any errors.\n")
cat("Demo data saved in:", demo_dir, "\n")
cat("Demo results saved in:", demo_results, "\n")
cat("Package appears to be working correctly if no errors were reported.\n")

# Clean up demo files (optional)
cat("\nCleaning up demo files...\n")
unlink(demo_dir, recursive = TRUE)
unlink(demo_results, recursive = TRUE)
cat("Demo cleanup completed.\n")

cat("\n=== Demo Complete ===\n")
