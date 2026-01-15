#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(here)
})

get_repo_root <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- args[grep("^--file=", args)]
    if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg[1])
    } else {
        script_path <- file.path(getwd(), "src", "scripts", "run_R_pipeline.R")
    }
    normalizePath(file.path(dirname(script_path), "..", ".."))
}

repo_root <- get_repo_root()
setwd(repo_root)
here::i_am("src/scripts/run_R_pipeline.R")

# Source all R files from src/R_tools/ to load functions directly
load_pipeline_functions <- function() {
    r_dir <- file.path(repo_root, "src", "R_tools")
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE, ignore.case = TRUE)
    for (f in r_files) {
        source(f, local = FALSE)
    }
    return(invisible(TRUE))
}

option_list <- list(
    make_option("--experiment-name", type = "character", default = "xen_tran_2024_12",
        help = "Experiment identifier used for data/results folders"),
    make_option("--experiment-date", type = "character",
        default = format(Sys.Date(), "%Y_%m_%d"),
        help = "Date string stored in the experiment metadata"),
    make_option("--assay-file", type = "character",
        default = file.path("data", "PLX073248_assay_data.csv"),
        help = "Path to the assay count matrix CSV"),
    make_option("--sample-file", type = "character",
        default = file.path("data", "PLX073248_sample_data.csv"),
        help = "Path to the sample metadata CSV"),
    make_option("--xenbase-file", type = "character",
        default = file.path("data", "Xenbase_Xenopus_Orthology_Predictions.txt"),
        help = "Path to Xenbase orthology predictions"),
    make_option("--hcop-file", type = "character",
        default = file.path("data", "HCOP_Xenopus_Orthology_Predictions.txt"),
        help = "Path to HCOP orthology predictions"),
    make_option("--suffix-mode", type = "character", default = "all",
        help = "Suffix handling mode for gene mapping"),
    make_option("--collections", type = "character", default = "H,C2,C5",
        help = "Comma-separated MSigDB collections for GSEA"),
    make_option("--reference-condition", type = "character", default = "Vehicle",
        help = "Reference condition level for DESeq2"),
    make_option("--reference-timepoint", type = "character", default = "T",
        help = "Reference timepoint level for DESeq2"),
    make_option("--group-cols", type = "character", default = "condition,timepoint,dose,drug,replicate",
        help = "Comma-separated metadata columns to use as grouping variables"),
    make_option("--comparison-cols", type = "character", default = "condition,timepoint",
        help = "Comma-separated columns for DESeq2 comparison"),
    make_option("--skip-mapping", action = "store_true", default = FALSE,
        help = "Skip the gene mapping step (useful when orthology files unavailable)"),
    make_option("--skip-gsea", action = "store_true", default = FALSE,
        help = "Skip the GSEA step"),
    make_option("--significance", type = "double", default = 0.05,
        help = "Adjusted p-value cutoff used in the JSON summary"),
    make_option("--summary-file", type = "character",
        default = file.path("results", "pipeline_summary.json"),
        help = "Path to write the pipeline summary JSON (placed inside the experiment folder)")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

load_pipeline_functions()

experiment_name <- opt$`experiment-name`
summary_path <- opt$`summary-file`

group_cols_vec <- strsplit(opt$`group-cols`, ",", fixed = TRUE)[[1]]
group_cols_vec <- trimws(group_cols_vec)
group_cols_vec <- group_cols_vec[nzchar(group_cols_vec)]

comparison_cols_vec <- strsplit(opt$`comparison-cols`, ",", fixed = TRUE)[[1]]
comparison_cols_vec <- trimws(comparison_cols_vec)
comparison_cols_vec <- comparison_cols_vec[nzchar(comparison_cols_vec)]

metadata_cols <- list(
    sample_id = "sample_id",
    group_cols = group_cols_vec
)

assay_path <- normalizePath(opt$`assay-file`, mustWork = TRUE)
sample_path <- normalizePath(opt$`sample-file`, mustWork = TRUE)

skip_mapping <- isTRUE(opt$`skip-mapping`)
skip_gsea <- isTRUE(opt$`skip-gsea`)

if (!skip_mapping) {
    xenbase_path <- normalizePath(opt$`xenbase-file`, mustWork = TRUE)
    hcop_path <- normalizePath(opt$`hcop-file`, mustWork = TRUE)
} else {
    xenbase_path <- NULL
    hcop_path <- NULL
}

experiment_config <- setup_experiment(
    experiment_prefix = experiment_name,
    experiment_date = opt$`experiment-date`,
    data_source = "local",
    data_type = "transcriptomics",
    species = "xenopus",
    metadata_columns = metadata_cols$group_cols,
    additional_notes = NULL
)

experiment_obj <- import_omics_data(
    data_source = "local",
    data_type = "transcriptomics",
    species = "xenopus",
    assay_path = assay_path,
    sample_path = sample_path,
    experiment_name = experiment_name,
    metadata_columns = metadata_cols$group_cols
)

clean_data_for_de(
    experiment_name = experiment_name,
    metadata_cols = metadata_cols
)

reference_levels <- list(
    condition = tolower(opt$`reference-condition`)
)
if ("timepoint" %in% comparison_cols_vec) {
    reference_levels$timepoint <- tolower(opt$`reference-timepoint`)
}

de_results <- run_deseq2_analysis(
    experiment_name = experiment_name,
    comparison_cols = comparison_cols_vec,
    reference_levels = reference_levels,
    interactive = FALSE
)

assay_entry <- experiment_obj$assay_data[[1]]
xenopus_genes <- rownames(assay_entry)

if (!skip_mapping) {
    mapping_results <- map_xenopus_to_human(
        xenopus_genes = xenopus_genes,
        xenbase_file = xenbase_path,
        hcop_file = hcop_path,
        suffix_mode = opt$`suffix-mode`,
        experiment_name = experiment_name,
        model_type = "deseq2"
    )
} else {
    mapping_results <- list(data = NULL)
    message("Skipping gene mapping step (--skip-mapping flag set)")
}

collection_vec <- strsplit(opt$collections, ",", fixed = TRUE)[[1]]
collection_vec <- trimws(collection_vec)
collection_vec <- collection_vec[nzchar(collection_vec)]

if (!skip_gsea && !skip_mapping) {
    gsea_results <- run_gsea_analysis(
        experiment_name = experiment_name,
        model_type = "deseq2",
        suffix_mode = opt$`suffix-mode`,
        mapping_type = "clean",
        collections = collection_vec
    )
} else {
    gsea_results <- NULL
    if (skip_gsea) {
        message("Skipping GSEA step (--skip-gsea flag set)")
    } else {
        message("Skipping GSEA step (requires gene mapping)")
    }
}

summarize_de <- function(results, padj_cutoff = 0.05) {
    lapply(names(results), function(name) {
        result <- results[[name]]
        deseq_tbl <- as.data.frame(result$deseq_results)
        deseq_tbl$gene_symbol <- rownames(result$deseq_results)
        deseq_tbl <- deseq_tbl[order(deseq_tbl$padj), ]
        list(
            contrast = name,
            experimental = result$comparison_info$experimental,
            control = result$comparison_info$control,
            significant_genes = sum(deseq_tbl$padj < padj_cutoff, na.rm = TRUE),
            top_genes = head(deseq_tbl[!is.na(deseq_tbl$padj),
                c("gene_symbol", "log2FoldChange", "padj")], 10)
        )
    })
}

summarize_gsea <- function(gsea_df, padj_cutoff = 0.05) {
    if (is.null(gsea_df) || nrow(gsea_df) == 0) {
        return(list())
    }
    gsea_df <- gsea_df[order(gsea_df$padj), ]
    list(
        significant_pathways = sum(gsea_df$padj < padj_cutoff, na.rm = TRUE),
        top_pathways = head(gsea_df[, c("pathway", "NES", "padj", "size")], 15)
    )
}

mapping_summary <- list(
    mapped_genes = if (!is.null(mapping_results$data)) nrow(mapping_results$data) else 0,
    unique_hugo_symbols = if (!is.null(mapping_results$data)) {
        length(unique(mapping_results$data$hugo_symbol))
    } else {
        0
    }
)

data_files_list <- list(assay = assay_path, sample = sample_path)
if (!is.null(xenbase_path)) data_files_list$xenbase <- xenbase_path
if (!is.null(hcop_path)) data_files_list$hcop <- hcop_path

params_list <- list(
    reference_condition = reference_levels$condition,
    padj_cutoff = opt$significance,
    suffix_mode = opt$`suffix-mode`,
    collections = collection_vec,
    comparison_cols = comparison_cols_vec
)
if (!is.null(reference_levels$timepoint)) {
    params_list$reference_timepoint <- reference_levels$timepoint
}

summary_list <- list(
    experiment = list(
        name = experiment_name,
        date = opt$`experiment-date`,
        data_files = data_files_list
    ),
    timestamps = list(
        completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    differential_expression = summarize_de(de_results, opt$significance),
    mapping = mapping_summary,
    gsea = summarize_gsea(gsea_results, opt$significance),
    parameters = params_list
)

summary_path <- normalizePath(summary_path, mustWork = FALSE)
summary_dir <- dirname(summary_path)
if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
}
write_json(summary_list, summary_path, auto_unbox = TRUE, pretty = TRUE)
message("Pipeline summary written to: ", summary_path)
