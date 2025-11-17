#!/usr/bin/env Rscript

install_if_missing <- function(pkgs) {
    installable <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(installable) == 0) {
        return(invisible(TRUE))
    }
    repos <- getOption("repos")
    if (is.null(repos) || repos["CRAN"] == "@CRAN@") {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
    install.packages(installable, dependencies = TRUE)
}

install_if_missing(c("optparse", "jsonlite", "here", "devtools"))

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(here)
})

# Determine repository root so here::here() works even when invoked from anywhere
get_repo_root <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- args[grep("^--file=", args)]
    if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg[1])
    } else {
        script_path <- file.path(getwd(), "scripts", "setup_r_env.R")
    }
    normalizePath(file.path(dirname(script_path), ".."))
}

repo_root <- get_repo_root()
setwd(repo_root)
here::i_am("scripts/setup_r_env.R")

# Install pipeline dependencies (CRAN + Bioc) before loading package code
ensure_pipeline_dependencies <- function() {
    setup_file <- file.path(repo_root, "R", "setup.R")
    if (!file.exists(setup_file)) {
        stop("Cannot find setup.R at ", setup_file)
    }
    repos <- getOption("repos")
    if (is.null(repos) || repos["CRAN"] == "@CRAN@") {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
    setup_env <- new.env(parent = globalenv())
    sys.source(setup_file, envir = setup_env)
    if (!is.function(setup_env$setup_packages)) {
        stop("setup_packages() not defined in setup.R")
    }
    message("Ensuring transcriptomicsPipeline dependencies are installed...")
    results <- setup_env$setup_packages(setup_env$setup_config)
    failed <- c(results$cran_failed, results$bioc_failed, results$github_failed)
    if (length(failed) > 0) {
        stop("Failed to install the following packages: ", paste(failed, collapse = ", "))
    }
    message("Dependency installation complete.")
}

# Load the package so we can access create_transcriptomics_environment()
load_pipeline_package <- function() {
    if ("package:transcriptomicsPipeline" %in% search()) {
        return(invisible(TRUE))
    }
    if (!requireNamespace("transcriptomicsPipeline", quietly = TRUE)) {
        if (!requireNamespace("devtools", quietly = TRUE)) {
            stop("devtools is required to load the transcriptomicsPipeline package from source.")
        }
        tryCatch(
            devtools::load_all(repo_root, quiet = TRUE),
            error = function(e) {
                message("devtools::load_all() failed: ", e$message)
                message("Attempting to install missing dependencies and retry...")
                ensure_pipeline_dependencies()
                devtools::load_all(repo_root, quiet = TRUE)
            }
        )
    }
    suppressPackageStartupMessages(library(transcriptomicsPipeline))
    return(invisible(TRUE))
}

option_list <- list(
    make_option("--base-dir",
        type = "character",
        default = repo_root,
        help = "Base directory for creating standard data/results folders [default: %default]"
    ),
    make_option("--skip-packages",
        action = "store_true",
        default = FALSE,
        help = "Skip installing R dependencies"
    ),
    make_option("--skip-directories",
        action = "store_true",
        default = FALSE,
        help = "Skip creating the standard directory tree"
    ),
    make_option("--output",
        type = "character",
        default = file.path("results", "environment_setup.json"),
        help = "Path to write a JSON summary of the environment setup"
    )
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

install_packages <- !isTRUE(opt$`skip-packages`)
if (install_packages) {
    ensure_pipeline_dependencies()
} else {
    message("Skipping dependency installation (--skip-packages).")
}
load_pipeline_package()
create_directories <- !isTRUE(opt$`skip-directories`)

message("Creating transcriptomics analysis environment...")
results <- transcriptomicsPipeline::create_transcriptomics_environment(
    install_packages = install_packages,
    create_directories = create_directories,
    base_dir = opt$`base-dir`
)

output_path <- normalizePath(opt$output, mustWork = FALSE)
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

write_json(results, output_path, auto_unbox = TRUE, pretty = TRUE)
message("Environment setup summary written to: ", output_path)
