# transcriptomicsPipeline

An R package for Xenopus laevis transcriptomics analysis. Provides automated workflows for differential expression analysis, gene mapping to human orthologs, and pathway enrichment.

## Overview

This package was developed to analyze transcriptomics data from Xenopus laevis experiments, particularly focusing on anesthetic drug responses. It handles the complete analysis pipeline from raw data import through pathway enrichment analysis.

## Key Features

- **Differential Expression**: DESeq2-based analysis with automated contrast generation
- **Gene Mapping**: Xenopus to human ortholog mapping using Xenbase and HCOP databases  
- **Pathway Analysis**: GSEA using MSigDB gene sets with human gene symbols
- **Visualization**: Automated generation of volcano plots, heatmaps, and PCA plots
- **Quality Control**: Comprehensive diagnostic plots and statistical reports

## Installation

```r
# Install from GitHub
devtools::install_github("agrossberg/transcriptomics_pipeline")

# Or install from source
git clone https://github.com/agrossberg/transcriptomics_pipeline.git
R CMD INSTALL transcriptomics_pipeline
```

## Quick Start

```r
library(transcriptomicsPipeline)

# 1. Import data
experiment_obj <- import_omics_data(
    data_source = "local",
    data_type = "transcriptomics", 
    species = "xenopus",
    assay_path = "assay_data.csv",
    sample_path = "sample_data.csv",
    experiment_name = "my_experiment"
)

# 2. Clean data
clean_data_for_de(
    experiment_name = "my_experiment",
    metadata_cols = list(
        sample_id = "sample_id",
        group_cols = c("condition", "timepoint")
    )
)

# 3. Run differential expression analysis
de_results <- run_deseq2_analysis(
    experiment_name = "my_experiment",
    comparison_cols = c("condition", "timepoint"),
    reference_levels = list(condition = "Vehicle", timepoint = "T"),
    interactive = FALSE
)

# 4. Map genes to human orthologs
mapping_results <- map_xenopus_to_human(
    xenopus_genes = gene_list,
    xenbase_file = "xenbase_orthologs.txt",
    hcop_file = "hcop_orthologs.txt",
    experiment_name = "my_experiment"
)

# 5. Run pathway enrichment
gsea_results <- run_gsea_analysis(
    experiment_name = "my_experiment",
    model_type = "deseq2",
    suffix_mode = "all",
    mapping_type = "clean",
    collections = c("H", "C2", "C5")
)
```

## Data Requirements

### Input Files
- **Assay data**: CSV with gene symbols in first column, samples as columns
- **Sample metadata**: CSV with sample information and experimental conditions
- **Orthology files**: Xenbase and HCOP orthology prediction files

### Sample Metadata Columns
- `sample_id`: Unique sample identifiers
- `condition`: Treatment conditions (e.g., Vehicle, Ketamine_100, etc.)
- `timepoint`: Time points (e.g., T, R)
- `dose`: Drug concentrations
- `drug`: Drug names

## Output Structure

```
results/
└── experiment_name/
    ├── de/                    # Differential expression results
    │   ├── csv/              # CSV files with results
    │   ├── tables/           # Formatted tables
    │   └── plots/            # Volcano plots and heatmaps
    ├── mapping/              # Gene mapping results
    ├── pca/                  # PCA analysis and plots
    └── gsea/                 # Pathway enrichment results
```

## Analysis Workflow

1. **Data Import**: Load expression data and sample metadata
2. **Quality Control**: Generate diagnostic plots and reports
3. **Differential Expression**: Identify significantly changed genes
4. **Gene Mapping**: Convert Xenopus genes to human orthologs
5. **Pathway Analysis**: Enrichment analysis using human gene sets
6. **Visualization**: Generate publication-ready plots

## Configuration

The package uses sensible defaults but allows customization through configuration objects:

- `dea_config`: DESeq2 analysis parameters
- `mapping_config`: Gene mapping settings
- `gsea_config`: Pathway analysis parameters
- `pca_config`: PCA analysis options

## Dependencies

### Core R Packages
- dplyr, tidyr, ggplot2, readr
- DESeq2, edgeR, limma
- msigdbr, fgsea
- ComplexHeatmap, viridis

### Bioconductor
- DESeq2, ComplexHeatmap
- msigdbr, fgsea
- org.Hs.eg.db, AnnotationDbi


## Contact

**Allison Grossberg**  
Wyss Institute for Biologically Inspired Engineering  
Harvard University  
Email: allison.grossberg@wyss.harvard.edu

## License

MIT License - see LICENSE file for details.