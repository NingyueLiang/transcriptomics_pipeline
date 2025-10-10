#!/usr/bin/env Rscript

#' Example Usage of transcriptomicsPipeline Package
#' 
#' This script demonstrates how to use the transcriptomicsPipeline package
#' with real data (you'll need to provide your own data files).

# Load the package
library(transcriptomicsPipeline)

cat("=== transcriptomicsPipeline Usage Example ===\n\n")

# Example 1: Setup and Install Dependencies
cat("1. Setting up package dependencies...\n")
cat("   Run: setup_results <- setup_packages()\n")
cat("   This will install all required CRAN and Bioconductor packages.\n\n")

# Example 2: Import Data
cat("2. Importing data from local files...\n")
cat("   # For local data import:\n")
cat("   experiment_obj <- import_omics_data(\n")
cat("       data_source = 'local',\n")
cat("       data_type = 'transcriptomics',\n")
cat("       species = 'human',\n")
cat("       assay_path = 'path/to/your/assay_data.csv',\n")
cat("       sample_path = 'path/to/your/sample_data.csv',\n")
cat("       experiment_name = 'my_experiment'\n")
cat("   )\n\n")

cat("   # For Pluto API import:\n")
cat("   experiment_obj <- import_omics_data(\n")
cat("       data_source = 'pluto',\n")
cat("       data_type = 'transcriptomics',\n")
cat("       species = 'xenopus',\n")
cat("       experiment_id = 'your_experiment_id',\n")
cat("       api_key = 'your_api_key'\n")
cat("   )\n\n")

# Example 3: Clean Data
cat("3. Cleaning data for analysis...\n")
cat("   cleaned_data_name <- clean_data_for_de(\n")
cat("       experiment_name = 'my_experiment',\n")
cat("       metadata_cols = list(\n")
cat("           sample_id = 'sample_id',\n")
cat("           group_cols = c('condition', 'treatment')\n")
cat("       )\n")
cat("   )\n\n")

# Example 4: Differential Expression Analysis
cat("4. Running differential expression analysis...\n")
cat("   de_results <- run_deseq2_analysis(\n")
cat("       experiment_name = 'my_experiment',\n")
cat("       comparison_cols = c('condition'),\n")
cat("       reference_levels = list(condition = c('control', 'treatment')),\n")
cat("       interactive = TRUE\n")
cat("   )\n\n")

# Example 5: Principal Component Analysis
cat("5. Performing PCA analysis...\n")
cat("   pca_results <- perform_pca_analysis(\n")
cat("       data = your_expression_data,\n")
cat("       groups_to_include = c('control', 'treatment'),\n")
cat("       experiment_name = 'my_experiment',\n")
cat("       save_results = TRUE\n")
cat("   )\n\n")

# Example 6: Gene Set Enrichment Analysis
cat("6. Running GSEA analysis...\n")
cat("   gsea_results <- run_gsea_analysis(\n")
cat("       experiment_name = 'my_experiment',\n")
cat("       model_type = 'deseq2',\n")
cat("       rank_method = 'sign_log_padj',\n")
cat("       collections = c('H', 'C2', 'C5')\n")
cat("   )\n\n")

# Example 7: Gene Mapping (for Xenopus data)
cat("7. Mapping genes (for Xenopus to human)...\n")
cat("   mapping_results <- map_xenopus_to_human(\n")
cat("       xenopus_genes = rownames(counts_matrix),\n")
cat("       xenbase_file = 'path/to/xenbase_file.txt',\n")
cat("       hcop_file = 'path/to/hcop_file.txt',\n")
cat("       suffix_mode = 'all'\n")
cat("   )\n\n")

# Example 8: Complete Workflow
cat("8. Complete workflow example...\n")
cat("   # Step 1: Setup\n")
cat("   setup_packages()\n\n")
cat("   # Step 2: Import\n")
cat("   experiment_obj <- import_omics_data(...)\n\n")
cat("   # Step 3: Clean\n")
cat("   cleaned_data <- clean_data_for_de(...)\n\n")
cat("   # Step 4: Analyze\n")
cat("   de_results <- run_deseq2_analysis(...)\n")
cat("   pca_results <- perform_pca_analysis(...)\n")
cat("   gsea_results <- run_gsea_analysis(...)\n\n")

cat("=== Data Requirements ===\n")
cat("For local data import, you need:\n")
cat("- Assay data file (CSV): Gene symbols in first column, sample counts in other columns\n")
cat("- Sample metadata file (CSV): Sample information including condition/treatment columns\n")
cat("- Gene symbol column should be named 'gene_symbol' or similar\n\n")

cat("For Pluto API import, you need:\n")
cat("- Valid Pluto API key\n")
cat("- Experiment ID from Pluto database\n\n")

cat("=== Output Structure ===\n")
cat("The package creates organized output directories:\n")
cat("results/\n")
cat("├── experiment_name/\n")
cat("│   ├── de/                 # Differential expression results\n")
cat("│   ├── mapping/            # Gene mapping results\n")
cat("│   ├── pca/               # PCA results and plots\n")
cat("│   └── gsea/              # GSEA results and plots\n\n")

cat("=== Getting Help ===\n")
cat("For detailed documentation, see the README.md file or run:\n")
cat("?function_name  # for help on specific functions\n")
cat("vignette('transcriptomicsPipeline')  # for package vignette (if available)\n\n")

cat("=== Example Complete ===\n")
