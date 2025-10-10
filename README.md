# transcriptomicsPipeline

A comprehensive R package for transcriptomics data analysis including data import, cleaning, differential expression analysis, gene mapping, PCA, and gene set enrichment analysis. Supports both Pluto API and local data sources with extensive validation and error handling.

## Features

- **Data Import**: Import data from Pluto API or local files with comprehensive validation
- **Data Cleaning**: Clean and prepare omics data for analysis with flexible filtering options
- **Differential Expression Analysis**: Perform DEA using DESeq2 with interactive contrast selection
- **Gene Mapping**: Map Xenopus genes to human orthologs using Xenbase and HCOP databases
- **Principal Component Analysis**: Comprehensive PCA with visualization and feature extraction
- **Gene Set Enrichment Analysis**: GSEA using MSigDB gene sets with multiple ranking methods
- **Package Setup**: Automated installation of all required dependencies

## Installation

### From GitHub

```r
# Install devtools if not already installed
if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install the package
devtools::install_github("agrossberg/transcriptomics_pipeline")
```

### From Source

```r
# Clone the repository
git clone https://github.com/agrossberg/transcriptomics_pipeline.git
cd transcriptomics_pipeline

# Install dependencies and build package
R CMD INSTALL .
```

## Quick Start

### 1. Setup and Install Dependencies

```r
library(transcriptomicsPipeline)

# Install all required packages
setup_results <- setup_packages()
```

### 2. Import Data

```r
# Import from Pluto API
experiment_obj <- import_omics_data(
    data_source = "pluto",
    data_type = "transcriptomics",
    species = "xenopus",
    experiment_id = "your_experiment_id",
    api_key = "your_api_key"
)

# Or import from local files
experiment_obj <- import_omics_data(
    data_source = "local",
    data_type = "transcriptomics",
    species = "xenopus",
    assay_path = "path/to/assay_data.csv",
    sample_path = "path/to/sample_data.csv"
)
```

### 3. Clean Data

```r
# Clean data for differential expression analysis
cleaned_data_name <- clean_data_for_de(
    experiment_name = "your_experiment",
    metadata_cols = list(
        sample_id = "sample_id",
        group_cols = c("condition", "treatment")
    )
)
```

### 4. Differential Expression Analysis

```r
# Run DESeq2 analysis
de_results <- run_deseq2_analysis(
    experiment_name = "your_experiment",
    comparison_cols = c("condition"),
    reference_levels = list(condition = c("control", "treatment")),
    interactive = TRUE
)
```

### 5. Gene Mapping (for Xenopus data)

```r
# Map Xenopus genes to human orthologs
mapping_results <- map_xenopus_to_human(
    xenopus_genes = rownames(counts_matrix),
    xenbase_file = "path/to/xenbase_file.txt",
    hcop_file = "path/to/hcop_file.txt",
    suffix_mode = "all"
)
```

### 6. Principal Component Analysis

```r
# Perform PCA analysis
pca_results <- perform_pca_analysis(
    data = your_expression_data,
    groups_to_include = c("control", "treatment"),
    experiment_name = "your_experiment",
    save_results = TRUE
)
```

### 7. Gene Set Enrichment Analysis

```r
# Run GSEA analysis
gsea_results <- run_gsea_analysis(
    experiment_name = "your_experiment",
    model_type = "deseq2",
    rank_method = "sign_log_padj",
    collections = c("H", "C2", "C5")
)
```

## Package Structure

The package is organized into the following modules:

- **setup.R**: Package installation and dependency management
- **import.R**: Data import functions for various sources
- **cleaning.R**: Data cleaning and preparation functions
- **differential_expression.R**: DESeq2-based differential expression analysis
- **mapping.R**: Gene mapping between species (Xenopus to human)
- **pca.R**: Principal component analysis and visualization
- **gsea.R**: Gene set enrichment analysis using MSigDB

## Configuration

The package uses configuration objects to manage default parameters:

- `setup_config`: Package installation configuration
- `import_config`: Data import configuration
- `cleaning_config`: Data cleaning configuration
- `dea_config`: Differential expression analysis configuration
- `mapping_config`: Gene mapping configuration
- `pca_config`: PCA analysis configuration
- `gsea_config`: GSEA analysis configuration

## Data Requirements

### For Pluto API Import
- Valid Pluto API key
- Experiment ID from Pluto database

### For Local Import
- Assay data file (CSV format) with gene symbols in first column
- Sample metadata file (CSV format) with sample information
- Gene symbol column should be named "gene_symbol" or similar

### For Gene Mapping
- Xenbase orthology predictions file
- HCOP orthology predictions file

## Output Structure

The package creates organized output directories:

```
results/
├── experiment_name/
│   ├── de/                 # Differential expression results
│   ├── mapping/            # Gene mapping results
│   ├── pca/               # PCA results and plots
│   └── gsea/              # GSEA results and plots
```

## Error Handling

The package includes comprehensive error handling and validation:

- Input parameter validation
- File existence checks
- Data format validation
- Graceful error messages with suggestions
- Logging for debugging

## Dependencies

### CRAN Packages
- dplyr, tidyr, readr, magrittr, tibble
- ggplot2, ggsignif, rstatix, viridis
- flextable, officer
- ggraph, tidygraph, igraph
- circlize, grid, RColorBrewer
- UpSetR, purrr, pluto
- jsonlite, cli, logger, here

### Bioconductor Packages
- DESeq2, ComplexHeatmap, limma
- TPP, MSnbase
- msigdbr, fgsea
- STRINGdb, org.Hs.eg.db, AnnotationDbi

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Allison Grossberg**  
Wyss Institute for Biologically Inspired Engineering  
Harvard University  
Email: allison.grossberg@wyss.harvard.edu
