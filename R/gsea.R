#' Gene Set Enrichment Analysis
#' 
#' GSEA analysis using MSigDB gene sets and fgsea.

# Required libraries
suppressPackageStartupMessages({
    library(dplyr)         # For data manipulation
    library(tidyr)         # For data tidying
    library(msigdbr)       # For accessing MSigDB gene sets
    library(fgsea)         # For running GSEA
    library(ggplot2)       # For creating lollipop plots
    library(stringr)       # For string manipulation
    library(tibble)        # For tibble operations
    library(purrr)         # For map functions
    library(grid)          # For graphics
    library(RColorBrewer)  # For color palettes
    library(ComplexHeatmap) # For legend
    library(circlize)      # For circos plots
    library(here)          # For path management
})

# Optional logger package
if (requireNamespace("logger", quietly = TRUE)) {
    library(logger)
} else {
    # Fallback logging function
    log_info <- function(...) message(...)
    log_error <- function(...) stop(...)
    log_warn <- function(...) warning(...)
}

# Configuration
gsea_config <- list(
    valid_rank_methods = c("log2fc", "sign_log_padj", "sign_log_pval"),
    valid_suffix_modes = c("L", "S", "average", "highest_effect", "all"),
    valid_mapping_types = c("clean", "messy"),
    valid_model_types = c("deseq2", "edgeR"),
    default_min_size = 10,
    default_max_size = 1000,
    msigdb_collections = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
)

#' Validate GSEA parameters
#' @param experiment_name Name of experiment
#' @param model_type Type of model
#' @param rank_method Ranking method
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @param suffix_mode Suffix mode (optional)
#' @param mapping_type Mapping type (optional)
validate_gsea_params <- function(experiment_name, model_type, rank_method, 
                                min_size, max_size, suffix_mode = NULL, mapping_type = NULL) {
    if (is.null(experiment_name) || nchar(experiment_name) == 0) {
        stop("experiment_name cannot be empty")
    }
    
    if (!model_type %in% gsea_config$valid_model_types) {
        stop("Invalid model_type. Must be one of: ", paste(gsea_config$valid_model_types, collapse = ", "))
    }
    
    if (!rank_method %in% gsea_config$valid_rank_methods) {
        stop("Invalid rank_method. Must be one of: ", paste(gsea_config$valid_rank_methods, collapse = ", "))
    }
    
    if (!is.numeric(min_size) || min_size < 1) {
        stop("min_size must be a positive integer")
    }
    
    if (!is.numeric(max_size) || max_size < min_size) {
        stop("max_size must be >= min_size")
    }
    
    if (!is.null(suffix_mode) && !suffix_mode %in% gsea_config$valid_suffix_modes) {
        stop("Invalid suffix_mode. Must be one of: ", paste(gsea_config$valid_suffix_modes, collapse = ", "))
    }
    
    if (!is.null(mapping_type) && !mapping_type %in% gsea_config$valid_mapping_types) {
        stop("Invalid mapping_type. Must be one of: ", paste(gsea_config$valid_mapping_types, collapse = ", "))
    }
}

#' Calculate ranking statistic with validation
#' @param data Data frame with gene expression results
#' @param method Ranking method
#' @return Named vector of ranking statistics
calculate_rank_stat <- function(data, method = "log2fc") {
    if (!is.data.frame(data) || nrow(data) == 0) {
        stop("data must be a non-empty data frame")
    }
    
    if (!method %in% gsea_config$valid_rank_methods) {
        stop("Invalid method. Must be one of: ", paste(gsea_config$valid_rank_methods, collapse = ", "))
    }
    
    message("Calculating ranking statistic using method: ", method)
    message("Initial number of rows: ", nrow(data))
    message("Initial number of unique genes: ", n_distinct(data$hugo_symbol))
    
    # Validate required columns
    required_cols <- c("hugo_symbol", "log2_fc", "p_value")
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Convert character values to numeric
    data <- data %>%
        mutate(
            p_value = as.numeric(p_value),
            adj_p_value = as.numeric(adj_p_value),
            log2_fc = as.numeric(log2_fc)
        )
    
    # Check for problematic values
    message("Checking for problematic values:")
    message("- NA in log2_fc: ", sum(is.na(data$log2_fc)))
    message("- NA in p_value: ", sum(is.na(data$p_value)))
    message("- Infinite values in log2_fc: ", sum(is.infinite(data$log2_fc)))
    message("- Zero p-values: ", sum(data$p_value == 0, na.rm = TRUE))
    
    # Handle zero p-values
    min_non_zero_pval <- min(data$p_value[data$p_value > 0], na.rm = TRUE)
    if (is.finite(min_non_zero_pval)) {
        data <- data %>%
            mutate(p_value = ifelse(p_value == 0, min_non_zero_pval/10, p_value))
    }
    
    # Remove problematic values
    clean_data <- data %>%
        filter(!is.na(log2_fc), 
               !is.na(p_value),
               is.finite(log2_fc), 
               p_value > 0)
    
    message("After cleaning: ", nrow(clean_data), " rows")
    
    if (nrow(clean_data) == 0) {
        stop("No valid data remaining after cleaning")
    }
    
    # Calculate ranking statistic based on method
    rank_stat <- switch(method,
        "log2fc" = {
            clean_data %>%
                arrange(desc(log2_fc)) %>%
                pull(log2_fc, name = hugo_symbol)
        },
        "sign_log_padj" = {
            clean_data %>%
                mutate(rank_stat = sign(log2_fc) * -log10(adj_p_value)) %>%
                arrange(desc(rank_stat)) %>%
                pull(rank_stat, name = hugo_symbol)
        },
        "sign_log_pval" = {
            clean_data %>%
                mutate(rank_stat = sign(log2_fc) * -log10(p_value)) %>%
                arrange(desc(rank_stat)) %>%
                pull(rank_stat, name = hugo_symbol)
        }
    )
    
    # Deduplicate gene names for fgsea (requires unique stats)
    if (any(duplicated(names(rank_stat)))) {
        dup_count <- sum(duplicated(names(rank_stat)))
        message("Collapsing ", dup_count, " duplicated gene entries for GSEA statistics.")
        rank_stat_tbl <- tibble(
            gene = names(rank_stat),
            score = as.numeric(rank_stat)
        ) %>%
            group_by(gene) %>%
            summarize(score = mean(score, na.rm = TRUE), .groups = "drop")
        rank_stat <- setNames(rank_stat_tbl$score, rank_stat_tbl$gene)
    }
    
    message("Final ranking statistic: ", length(rank_stat), " genes")
    return(rank_stat)
}

#' Load MSigDB gene sets with validation
#' @param collections MSigDB collections to load
#' @param species Species for gene sets
#' @return List of gene sets
load_msigdb_genesets <- function(collections = gsea_config$msigdb_collections, 
                                species = "Homo sapiens") {
    if (!is.character(collections) || length(collections) == 0) {
        stop("collections must be a non-empty character vector")
    }
    
    message("Loading MSigDB gene sets for species: ", species)
    message("Collections: ", paste(collections, collapse = ", "))
    
    tryCatch({
        genesets <- list()
        
        for (collection in collections) {
            message("Loading collection: ", collection)
            collection_data <- msigdbr(species = species, category = collection)
            
            if (nrow(collection_data) > 0) {
                # Convert to list format for fgsea
                collection_genesets <- split(collection_data$gene_symbol, collection_data$gs_name)
                genesets <- c(genesets, collection_genesets)
                message("Loaded ", length(collection_genesets), " gene sets from ", collection)
            } else {
                warning("No gene sets found for collection: ", collection)
            }
        }
        
        message("Total gene sets loaded: ", length(genesets))
        return(genesets)
    }, error = function(e) {
        stop("Failed to load MSigDB gene sets: ", e$message)
    })
}

#' Run GSEA for multiple collections
#' @param rank_stat Named vector of ranking statistics
#' @param genesets List of gene sets
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @param nperm Number of permutations
#' @return GSEA results
run_multiple_gsea <- function(rank_stat, genesets, min_size = gsea_config$default_min_size, 
                             max_size = gsea_config$default_max_size, nperm = 10000) {
    if (is.null(rank_stat) || length(rank_stat) == 0) {
        stop("rank_stat cannot be empty")
    }
    
    if (is.null(genesets) || length(genesets) == 0) {
        stop("genesets cannot be empty")
    }
    
    message("Running GSEA with ", length(genesets), " gene sets")
    message("Gene set size range: ", min_size, " - ", max_size)
    message("Number of permutations: ", nperm)
    
    tryCatch({
        # Filter gene sets by size
        filtered_genesets <- genesets[sapply(genesets, function(x) {
            length(intersect(x, names(rank_stat))) >= min_size && 
            length(intersect(x, names(rank_stat))) <= max_size
        })]
        
        message("Gene sets after size filtering: ", length(filtered_genesets))
        
        if (length(filtered_genesets) == 0) {
            warning("No gene sets pass size filtering criteria")
            return(data.frame())
        }
        
        # Run GSEA
        gsea_results <- fgsea(
            pathways = filtered_genesets,
            stats = rank_stat,
            minSize = min_size,
            maxSize = max_size,
            nperm = nperm
        )
        
        message("GSEA completed. Found ", sum(gsea_results$padj < 0.05, na.rm = TRUE), 
                " significantly enriched pathways")
        
        return(gsea_results)
    }, error = function(e) {
        stop("GSEA analysis failed: ", e$message)
    })
}

#' Create lollipop plot for significant pathways
#' @param gsea_results GSEA results data frame
#' @param top_n Number of top pathways to plot
#' @param title Plot title
#' @return ggplot object
create_lollipop_plot <- function(gsea_results, top_n = 20, title = "Top Enriched Pathways") {
    if (is.null(gsea_results) || nrow(gsea_results) == 0) {
        warning("No GSEA results to plot")
        return(NULL)
    }
    
    # Filter significant results and sort by NES
    significant_results <- gsea_results %>%
        filter(padj < 0.05) %>%
        arrange(desc(abs(NES))) %>%
        head(top_n)
    
    if (nrow(significant_results) == 0) {
        warning("No significant pathways found for plotting")
        return(NULL)
    }
    
    # Create plot
    p <- ggplot(significant_results, aes(x = reorder(pathway, NES), y = NES)) +
        geom_segment(aes(xend = pathway, yend = 0), color = "grey50") +
        geom_point(aes(color = padj, size = size)) +
        scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
        scale_size_continuous(name = "Gene Set Size", range = c(2, 8)) +
        coord_flip() +
        labs(
            title = title,
            x = "Pathway",
            y = "Normalized Enrichment Score (NES)"
        ) +
        theme_minimal() +
        theme(
            axis.text.y = element_text(size = 8),
            legend.position = "right"
        )
    
    return(p)
}

#' Export GSEA results
#' @param gsea_results GSEA results
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @param model_type Model type
#' @param suffix_mode Suffix mode (optional)
#' @param mapping_type Mapping type (optional)
export_gsea_results <- function(gsea_results, output_dir, experiment_name, model_type, 
                               suffix_mode = NULL, mapping_type = NULL) {
    if (is.null(gsea_results) || nrow(gsea_results) == 0) {
        warning("No GSEA results to export")
        return()
    }
    
    # Create output directory
    if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Failed to create output directory: ", output_dir)
    }
    
    # Create filename components
    filename_parts <- c(experiment_name, model_type)
    if (!is.null(suffix_mode)) {
        filename_parts <- c(filename_parts, suffix_mode)
    }
    if (!is.null(mapping_type)) {
        filename_parts <- c(filename_parts, mapping_type)
    }
    
    base_filename <- paste(filename_parts, collapse = "_")
    
    tryCatch({
        # Export CSV
        csv_file <- file.path(output_dir, paste0(base_filename, "_gsea_results.csv"))
        write_csv(gsea_results, csv_file)
        message("GSEA results exported to: ", csv_file)
        
        # Create and save lollipop plot
        plot <- create_lollipop_plot(gsea_results, title = paste("GSEA Results:", experiment_name))
        if (!is.null(plot)) {
            plot_file <- file.path(output_dir, paste0(base_filename, "_lollipop_plot.png"))
            ggsave(plot_file, plot, width = 12, height = 8, dpi = 300)
            message("Lollipop plot saved to: ", plot_file)
        }
        
        # Save RDS file
        rds_file <- file.path(output_dir, paste0(base_filename, "_gsea_results.rds"))
        saveRDS(gsea_results, rds_file)
        message("GSEA results saved to: ", rds_file)
        
    }, error = function(e) {
        stop("Failed to export GSEA results: ", e$message)
    })
}

#' Load differential expression data for GSEA
#' @param experiment_name Name of experiment
#' @param model_type Type of model
#' @param suffix_mode Suffix mode (optional)
#' @param mapping_type Mapping type (optional)
#' @return Data frame with DE results
load_de_data_for_gsea <- function(experiment_name, model_type, suffix_mode = NULL, 
                                 mapping_type = NULL) {
    # Always use unmapped DE results and apply mapping on the fly
    data_file <- here("results", experiment_name, "de", model_type, "simple", 
                     paste0(experiment_name, "_de_results.rds"))
    
    if (!file.exists(data_file)) {
        stop("Data file not found: ", data_file)
    }
    
    tryCatch({
        data <- readRDS(data_file)
        
        # Extract DE results and apply mapping if needed
        if (is.list(data) && length(data) > 0) {
            # Get first contrast results
            first_result <- data[[1]]
            if (is.list(first_result) && "deseq_results" %in% names(first_result)) {
                de_data <- as.data.frame(first_result$deseq_results)
                # Handle case where rownames are missing
                if (length(rownames(de_data)) == 0) {
                    de_data$xenopus_symbol <- paste0("Gene_", 1:nrow(de_data))
                } else {
                    de_data$xenopus_symbol <- rownames(de_data)
                }
                
                # Apply mapping if suffix_mode and mapping_type are provided
                if (!is.null(suffix_mode) && !is.null(mapping_type)) {
                    # Load mapping dictionary to convert Xenopus to human genes
                    mapping_file <- here("results", experiment_name, "mapping", model_type, suffix_mode, 
                                       paste0(experiment_name, "_", suffix_mode, "_", model_type, "_mapping.csv"))
                    
                    if (file.exists(mapping_file)) {
                        mapping_dict <- read.csv(mapping_file, stringsAsFactors = FALSE)
                        
                        # Create clean gene symbols by removing suffixes for mapping
                        de_data$clean_xenopus_symbol <- gsub("\\.L$", "", de_data$xenopus_symbol)
                        de_data$clean_xenopus_symbol <- gsub("\\.S$", "", de_data$clean_xenopus_symbol)
                        
                        # Merge with mapping dictionary using clean symbols
                        de_data <- merge(de_data, mapping_dict, 
                                       by.x = "clean_xenopus_symbol", by.y = "gene", 
                                       all.x = TRUE)
                        
                        # Use human symbol as hugo_symbol, fallback to xenopus if no mapping
                        de_data$hugo_symbol <- ifelse(!is.na(de_data$hugo_symbol) & de_data$hugo_symbol != "", 
                                                     de_data$hugo_symbol, de_data$xenopus_symbol)
                        
                        message("Applied mapping: ", sum(!is.na(de_data$hugo_symbol) & de_data$hugo_symbol != ""), 
                               " genes mapped to human symbols")
                    } else {
                        warning("Mapping file not found: ", mapping_file, ". Using Xenopus symbols.")
                        de_data$hugo_symbol <- de_data$xenopus_symbol
                    }
                } else {
                    # No mapping requested, use xenopus symbols as hugo symbols
                    de_data$hugo_symbol <- de_data$xenopus_symbol
                }
            } else {
                stop("Unexpected DE results structure")
            }
        } else {
            stop("No DE results found")
        }
        
        # Standardize column names
        colnames(de_data)[colnames(de_data) == "log2FoldChange"] <- "log2_fc"
        colnames(de_data)[colnames(de_data) == "pvalue"] <- "p_value"
        colnames(de_data)[colnames(de_data) == "padj"] <- "adj_p_value"
        
        message("Loaded DE data: ", nrow(de_data), " genes")
        return(de_data)
        
    }, error = function(e) {
        stop("Failed to load DE data: ", e$message)
    })
}

#' Main GSEA analysis function
#' @param experiment_name Name of experiment
#' @param model_type Type of model
#' @param rank_method Ranking method
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @param suffix_mode Suffix mode (optional)
#' @param mapping_type Mapping type (optional)
#' @param collections MSigDB collections to use
#' @param nperm Number of permutations
#' @return GSEA results
run_gsea_analysis <- function(experiment_name,
                             model_type = "deseq2",
                             rank_method = "sign_log_padj",
                             min_size = gsea_config$default_min_size,
                             max_size = gsea_config$default_max_size,
                             suffix_mode = NULL,
                             mapping_type = NULL,
                             collections = gsea_config$msigdb_collections,
                             nperm = 10000) {
    
    # Validate parameters
    validate_gsea_params(experiment_name, model_type, rank_method, 
                        min_size, max_size, suffix_mode, mapping_type)
    
    message("Starting GSEA analysis for: ", experiment_name)
    message("Model type: ", model_type)
    message("Rank method: ", rank_method)
    if (!is.null(suffix_mode)) message("Suffix mode: ", suffix_mode)
    if (!is.null(mapping_type)) message("Mapping type: ", mapping_type)
    
    # Load differential expression data
    de_data <- load_de_data_for_gsea(experiment_name, model_type, suffix_mode, mapping_type)
    
    # Calculate ranking statistic
    rank_stat <- calculate_rank_stat(de_data, rank_method)
    
    # Load MSigDB gene sets
    genesets <- load_msigdb_genesets(collections)
    
    # Run GSEA
    gsea_results <- run_multiple_gsea(rank_stat, genesets, min_size, max_size, nperm)
    
    # Create output directory
    output_dir <- here("results", experiment_name, "gsea", model_type)
    if (!is.null(suffix_mode)) {
        output_dir <- file.path(output_dir, suffix_mode)
    }
    if (!is.null(mapping_type)) {
        output_dir <- file.path(output_dir, mapping_type)
    }
    
    # Export results
    export_gsea_results(gsea_results, output_dir, experiment_name, model_type, 
                       suffix_mode, mapping_type)
    
    message("GSEA analysis completed successfully")
    return(gsea_results)
}

# Additional functions from original implementation

#' Calculate ranking statistic with comprehensive validation
#' @param data Data frame with gene expression results
#' @param method Ranking method
#' @return Data frame with ranking statistics
calculate_rank_stat_comprehensive <- function(data, method = "log2fc") {
  # Initial data state
  message("\nData state tracking:")
  message("Initial number of rows: ", nrow(data))
  message("Initial number of unique genes: ", n_distinct(data$hugo_symbol))
  
  # Convert character p-values back to numeric
  data <- data %>%
    mutate(
      p_value = as.numeric(p_value),
      adj_p_value = as.numeric(adj_p_value),
      log2_fc = as.numeric(log2_fc)
    )
  
  message("\nChecking for problematic values:")
  message("- NA in log2_fc: ", sum(is.na(data$log2_fc)))
  message("- NA in p_value: ", sum(is.na(data$p_value)))
  message("- Infinite values in log2_fc: ", sum(is.infinite(as.numeric(data$log2_fc))))
  message("- Zero p-values: ", sum(data$p_value == 0))  # Add this check
  
  # Handle zero p-values before log transformation
  min_non_zero_pval <- min(data$p_value[data$p_value > 0], na.rm = TRUE)
  data <- data %>%
    mutate(p_value = ifelse(p_value == 0, min_non_zero_pval/10, p_value))
  
  # Remove NA and infinite values
  clean_data <- data %>%
    filter(!is.na(log2_fc), 
           !is.na(p_value),
           is.finite(log2_fc), 
           p_value > 0)  # Ensure no zero p-values
  
  message("\nAfter removing NA/Inf/zero values:")
  message("Number of rows: ", nrow(clean_data))
  
  # Remove duplicates
  final_data <- clean_data %>%
    distinct(hugo_symbol, .keep_all = TRUE)
  
  message("\nAfter removing duplicates:")
  message("Final number of rows: ", nrow(final_data))
  
  # Calculate ranking statistic with additional checks
  final_data <- final_data %>%
    mutate(rank_stat = case_when(
      method == "log2fc" ~ log2_fc,
      method == "pval" ~ -log10(p_value),
      method == "padj" ~ -log10(adj_p_value),
      method == "sign_log_pval" ~ sign(log2_fc) * -log10(p_value),
      method == "sign_log_padj" ~ sign(log2_fc) * -log10(adj_p_value),
      method == "log2fc_log_pval" ~ log2_fc * -log10(p_value),
      method == "log2fc_log_padj" ~ log2_fc * -log10(adj_p_value),
      TRUE ~ log2_fc
    ))
  
  # Final check for infinite values
  if(any(is.infinite(final_data$rank_stat))) {
    message("\nWarning: Infinite values in rank statistic after calculation")
    final_data <- final_data %>%
      filter(is.finite(rank_stat))
    message("Rows after removing infinite rank stats: ", nrow(final_data))
  }
  
  return(final_data)
}

#' Run multiple GSEA analyses with comprehensive collection handling
#' @param diff_expr_data Differential expression data
#' @param rank_method Ranking method
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @return List containing GSEA results and summary
run_multiple_gsea_comprehensive <- function(diff_expr_data, 
                            rank_method = "log2fc",
                            min_size = 10, 
                            max_size = 1000) {
  message("Starting GSEA analysis with ", rank_method, " ranking method")
  
  # Clean and prepare data with ranking statistic
  clean_data <- if("hugo_symbol" %in% colnames(diff_expr_data)) {
    # For mapped data
    diff_expr_data %>%
      mutate(hugo_symbol = toupper(hugo_symbol))
  } else if("gene_id" %in% colnames(diff_expr_data)) {
    # For unmapped data
    diff_expr_data %>%
      mutate(hugo_symbol = toupper(gene_id))  # Use gene_id as hugo_symbol
  } else {
    stop("Neither hugo_symbol nor gene_id column found in data")
  }
  
  clean_data <- calculate_rank_stat_comprehensive(clean_data, rank_method)
  
  # Define all collections including Hallmark
  collections <- list(
    list(category = "H"),
    list(category = "C1"),
    list(category = "C2", subcategory = "CGP"),
    list(category = "C2", subcategory = "CP"),
    list(category = "C2", subcategory = "CP:BIOCARTA"),
    list(category = "C2", subcategory = "CP:KEGG"),
    list(category = "C2", subcategory = "CP:PID"),
    list(category = "C2", subcategory = "CP:REACTOME"),
    list(category = "C2", subcategory = "CP:WIKIPATHWAYS"),
    list(category = "C3", subcategory = "MIR:MIR_Legacy"),  
    list(category = "C3", subcategory = "MIR:MIRDB"),  
    list(category = "C3", subcategory = "TFT:GTRD"),  
    list(category = "C3", subcategory = "TFT:TFT_Legacy"),  
    list(category = "C4", subcategory = "CGN"),
    list(category = "C4", subcategory = "CM"),
    list(category = "C5", subcategory = "GO:BP"),
    list(category = "C5", subcategory = "GO:CC"),
    list(category = "C5", subcategory = "GO:MF"),
    list(category = "C5", subcategory = "HPO"),
    list(category = "C6"),
    list(category = "C7", subcategory = "IMMUNESIGDB"),
    list(category = "C7", subcategory = "VAX"),
    list(category = "C8")
  )
  
  # Function to run single GSEA analysis
  run_single_gsea <- function(collection_info) {
    collection_name <- paste(c(collection_info$category, collection_info$subcategory), 
                           collapse = ":")
    message("\nProcessing collection: ", collection_name)
    
    tryCatch({
      # Get gene sets
      gene_sets <- msigdbr(species = "Homo sapiens", 
                          category = collection_info$category,
                          subcategory = collection_info$subcategory)
      
      if(nrow(gene_sets) == 0) {
        message("Warning: No gene sets found for ", collection_name)
        return(tibble())
      }
      
      gene_sets <- gene_sets %>%
        dplyr::select(gs_name, gene_symbol) %>%
        split(x = .$gene_symbol, f = .$gs_name)
      
      # Create ranked list
      ranks <- clean_data %>%
        dplyr::select(hugo_symbol, rank_stat) %>%
        deframe()
      
      # Run GSEA
      results <- fgseaMultilevel(pathways = gene_sets,
                                stats = ranks,
                                minSize = min_size,
                                maxSize = max_size) %>%
        as_tibble() %>%
        arrange(padj) %>%
        mutate(collection = collection_name)
      
      return(results)
    }, error = function(e) {
      message("Error processing ", collection_name, ": ", e$message)
      return(tibble())
    })
  }
  
  # Run GSEA for all collections
  all_results <- map_df(collections, run_single_gsea)
  
  # Create summary
  summary <- all_results %>%
    group_by(collection) %>%
    summarise(
      total_pathways = n(),
      significant_pathways = sum(padj < 0.05),
      min_padj = min(padj, na.rm = TRUE)
    )
  
  return(list(
    results = all_results,
    summary = summary,
    ranking_method = rank_method
  ))
}

#' Create lollipop plot for GSEA results
#' @param gsea_results GSEA results
#' @param collection_name Collection name
#' @param padj_cutoff P-value cutoff
#' @return ggplot object
create_gsea_lollipop <- function(gsea_results, collection_name, padj_cutoff = 0.05) {
  # Filter for significant pathways
  plot_data <- gsea_results %>%
    filter(collection == collection_name, padj < padj_cutoff) %>%
    arrange(padj)
  
  if(nrow(plot_data) == 0) {
    message("No significant pathways found for ", collection_name)
    return(NULL)
  }
  
  max_log_padj <- -log10(min(plot_data$padj))
  
  p <- ggplot(plot_data, aes(x = NES, y = reorder(pathway, NES))) +
    geom_segment(aes(x = 0, xend = NES, 
                     y = reorder(pathway, NES), 
                     yend = reorder(pathway, NES)),
                 color = "grey80") +
    geom_point(aes(color = -log10(padj), size = size), alpha = 0.8) +
    scale_color_gradient(
      low = "blue", high = "red",
      limits = c(0, max_log_padj),
      name = "Adjusted\np-value",
      labels = function(x) paste0("1e-", round(x))
    ) +
    scale_size_continuous(range = c(2, 8), name = "Gene Set Size") +
    labs(title = paste("Gene Set Enrichment Analysis:", collection_name),
         subtitle = paste0("Significant pathways (FDR < ", padj_cutoff, ")"),
         x = "Normalized Enrichment Score (NES)",
         y = NULL) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Create all lollipop plots
#' @param gsea_results GSEA results
#' @param padj_cutoff P-value cutoff
#' @return List of plots
create_all_lollipop_plots <- function(gsea_results, padj_cutoff = 0.05) {
  collections <- unique(gsea_results$collection)
  plots <- list()
  
  for(collection in collections) {
    plot <- create_gsea_lollipop(gsea_results, collection, padj_cutoff)
    if(!is.null(plot)) {
      plots[[collection]] <- plot
    }
  }
  
  return(plots)
}

#' Export GSEA results with comprehensive file handling
#' @param all_gsea_results GSEA results
#' @param all_plots Plots
#' @param experiment_name Experiment name
#' @param contrast_name Contrast name
#' @param model_type Model type
#' @param output_dir Output directory
#' @param suffix_mode Suffix mode
#' @param mapping_type Mapping type
export_gsea_results_comprehensive <- function(all_gsea_results, all_plots, experiment_name, contrast_name,
                              model_type,
                              output_dir = "gsea_results", suffix_mode = NULL, mapping_type = NULL) {
  # Create base directory
  base_dir <- file.path("results", experiment_name, "gsea", model_type)
  
  # If this is mapped data, add suffix_mode and mapping_type to the path
  if (!is.null(suffix_mode) && !is.null(mapping_type)) {
    output_dir <- file.path(base_dir, suffix_mode, mapping_type)
  } else {
    output_dir <- base_dir  # For unmapped data, use simpler path
  }
  
  # Create directories
  plots_dir <- file.path(output_dir, "plots")
  tables_dir <- file.path(output_dir, "tables")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Export results table
  significant_results <- all_gsea_results$results %>%
    filter(padj < 0.05) %>%
    arrange(collection, padj) %>%
    mutate(
      FDR = format(padj, scientific = FALSE, digits = 10),
      pval = format(pval, scientific = FALSE, digits = 10)
    ) %>%
    select(collection, pathway, NES, size, pval, FDR)
  
  write.csv(significant_results,
            file.path(tables_dir, 
                     paste0(experiment_name, "_", contrast_name, 
                           "_significant_gsea_results.csv")),
            row.names = FALSE)
  
  # Export collection-specific tables
  split_results <- split(significant_results, significant_results$collection)
  for(collection_name in names(split_results)) {
    safe_collection_name <- gsub(":", "_", collection_name)
    write.csv(split_results[[collection_name]],
              file.path(tables_dir,
                       paste0(experiment_name, "_", contrast_name, "_",
                             safe_collection_name, "_results.csv")),
              row.names = FALSE)
  }
  
  # Save plots
  for(name in names(all_plots)) {
    safe_name <- gsub(":", "_", name)
    ggsave(
      filename = file.path(plots_dir,
                          paste0(experiment_name, "_", contrast_name, "_",
                                safe_name, "_lollipop.png")),
      plot = all_plots[[name]],
      width = 12,
      height = 8,
      dpi = 300
    )
  }
  
  cat("GSEA results exported to:", output_dir, "\n")
}

#' Determine data type and set up directory paths
#' @param experiment_name Experiment name
#' @param model_type Model type
#' @param suffix_mode Suffix mode
#' @param mapping_type Mapping type
#' @return List containing data type info and paths
setup_gsea_paths <- function(experiment_name, model_type, suffix_mode, mapping_type) {
    # Determine if this is mapped or unmapped data
    is_mapped_data <- !is.null(suffix_mode) && !is.null(mapping_type)
    
    # Set the appropriate directory path based on data type
    if(is_mapped_data) {
        base_dir <- file.path("results", experiment_name, "mapping", model_type, suffix_mode)
        file_pattern <- paste0(".*_DEGs_", mapping_type, "\\.csv$")
    } else {
        base_dir <- file.path("results", experiment_name, "de", model_type, "simple")
        file_pattern <- ".*_de_results\\.rds$"
    }
    
    if(!dir.exists(base_dir)) {
        stop("Results directory not found at: ", base_dir)
    }
    
    return(list(
        is_mapped_data = is_mapped_data,
        base_dir = base_dir,
        file_pattern = file_pattern
    ))
}

#' Process mapped data (CSV files)
#' @param base_dir Base directory path
#' @param file_pattern File pattern to match
#' @param experiment_name Experiment name
#' @param rank_method Ranking method
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @param model_type Model type
#' @param suffix_mode Suffix mode
#' @param mapping_type Mapping type
#' @param debug Debug flag
#' @return List of GSEA results
process_mapped_data <- function(base_dir, file_pattern, experiment_name, rank_method, 
                               min_size, max_size, model_type, suffix_mode, mapping_type, debug) {
    # For mapped data, process CSV files
    deg_files <- list.files(base_dir, 
                           pattern = file_pattern, 
                           full.names = TRUE)
    if(length(deg_files) == 0) {
        stop("No DEG files found in ", base_dir)
    }
    
    gsea_results <- list()
    
    # Process each contrast
    for(deg_file in deg_files) {
        # Extract contrast name from filename
        contrast_name <- sub(paste0(experiment_name, "_"), "", basename(deg_file))
        contrast_name <- sub("_DEGs\\.csv$", "", contrast_name)
        
        if(debug) cat("\nProcessing contrast:", contrast_name, "\n")
        
        # Read the CSV file
        deg_data <- read.csv(deg_file)
        
        # Run GSEA for this contrast
        contrast_results <- run_multiple_gsea_comprehensive(
            diff_expr_data = deg_data,
            rank_method = rank_method,
            min_size = min_size,
            max_size = max_size
        )
        
        # Store results and create plots
        gsea_results[[contrast_name]] <- contrast_results
        lollipop_plots <- create_all_lollipop_plots(contrast_results$results)
        
        # Export results
        export_gsea_results_comprehensive(
            contrast_results,
            lollipop_plots,
            experiment_name = experiment_name,
            contrast_name = contrast_name,
            model_type = model_type,
            suffix_mode = suffix_mode,
            mapping_type = mapping_type
        )
    }
    
    return(gsea_results)
}

#' Process unmapped data (RDS files)
#' @param base_dir Base directory path
#' @param file_pattern File pattern to match
#' @param experiment_name Experiment name
#' @param rank_method Ranking method
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @param model_type Model type
#' @param suffix_mode Suffix mode
#' @param mapping_type Mapping type
#' @param debug Debug flag
#' @return List of GSEA results
process_unmapped_data <- function(base_dir, file_pattern, experiment_name, rank_method, 
                                 min_size, max_size, model_type, suffix_mode, mapping_type, debug) {
    # For unmapped data, process RDS file
    rds_files <- list.files(base_dir, 
                           pattern = file_pattern, 
                           full.names = TRUE)
    if(length(rds_files) == 0) {
        stop("No DESeq2 results files found in ", base_dir)
    }
    
    # Load the RDS file and check contrast names
    de_results <- readRDS(rds_files[1])
    cat("Contrast names in DESeq2 results:\n")
    print(names(de_results))
    
    gsea_results <- list()
    
    # Process each contrast in the DESeq2 results
    for(contrast_name in names(de_results)) {
        # Sanitize contrast name for file paths
        safe_contrast_name <- gsub("[^[:alnum:]]", "_", contrast_name)
        if(debug) cat("\nProcessing contrast:", contrast_name, 
                      "\nSanitized name:", safe_contrast_name, "\n")
        
        # Extract DESeq2 results for this contrast
        contrast_data <- de_results[[contrast_name]]$deseq_results
        
        if(debug) {
            cat("\nContrast data structure:\n")
            print(str(contrast_data))
        }
        
        # Create data frame with required columns
        deg_data <- data.frame(
            gene_id = rownames(contrast_data),
            log2_fc = contrast_data$log2FoldChange,
            p_value = contrast_data$pvalue,
            adj_p_value = contrast_data$padj,
            base_mean = contrast_data$baseMean
        )
        
        if(debug) {
            cat("\nDEG data dimensions:", dim(deg_data), "\n")
            cat("Column names:", colnames(deg_data), "\n")
        }
        
        # Run GSEA for this contrast
        contrast_results <- run_multiple_gsea_comprehensive(
            diff_expr_data = deg_data,
            rank_method = rank_method,
            min_size = min_size,
            max_size = max_size
        )
        
        # Store results and create plots
        gsea_results[[contrast_name]] <- contrast_results
        lollipop_plots <- create_all_lollipop_plots(contrast_results$results)
        
        # Export results with sanitized name
        export_gsea_results_comprehensive(
            contrast_results,
            lollipop_plots,
            experiment_name = experiment_name,
            contrast_name = safe_contrast_name,  # Use sanitized name
            model_type = model_type,
            suffix_mode = suffix_mode,
            mapping_type = mapping_type
        )
    }
    
    return(gsea_results)
}

#' Save GSEA results to RDS file
#' @param gsea_results GSEA results
#' @param experiment_name Experiment name
#' @param model_type Model type
#' @param suffix_mode Suffix mode
#' @param mapping_type Mapping type
#' @param debug Debug flag
#' @return Path to saved results file
save_gsea_results <- function(gsea_results, experiment_name, model_type, suffix_mode, mapping_type, debug) {
    # Save GSEA results to RDS file - simplified path handling
    gsea_dir <- if (!is.null(suffix_mode) && !is.null(mapping_type)) {
        # For mapped data
        file.path("results", experiment_name, "gsea", model_type, suffix_mode, mapping_type)
    } else {
        # For unmapped data - use simple path
        file.path("results", experiment_name, "gsea", model_type)
    }
    
    # Create directory and save results
    dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Simplified file naming
    results_file <- if (!is.null(suffix_mode) && !is.null(mapping_type)) {
        # For mapped data
        paste0(experiment_name, "_", model_type, "_", suffix_mode, "_", mapping_type, "_gsea_results.rds")
    } else {
        # For unmapped data
        paste0(experiment_name, "_", model_type, "_gsea_results.rds")
    }
    
    # Save results
    saveRDS(gsea_results, file.path(gsea_dir, results_file))
    
    if(debug) {
        cat("\nGSEA directory:", gsea_dir)
        cat("\nResults file:", results_file)
        cat("\nFull path:", file.path(gsea_dir, results_file), "\n")
    }
    
    return(file.path(gsea_dir, results_file))
}

#' Main GSEA wrapper function with comprehensive workflow
#' @param experiment_name Experiment name
#' @param model_type Model type
#' @param rank_method Ranking method
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @param suffix_mode Suffix mode
#' @param mapping_type Mapping type
#' @param debug Debug flag
#' @return GSEA results
run_gsea_analysis_comprehensive <- function(experiment_name,
                            model_type = "deseq2",
                            rank_method = "log2fc_log_pval",
                            min_size = 10,
                            max_size = 1000,
                            suffix_mode = NULL,
                            mapping_type = NULL,
                            debug = TRUE) {
    
    if(debug) cat("Starting GSEA analysis for", experiment_name, "using", model_type, "\n")
    
    # 1. Set up paths and determine data type
    path_info <- setup_gsea_paths(experiment_name, model_type, suffix_mode, mapping_type)
    
    # 2. Process data based on type
    if(path_info$is_mapped_data) {
        gsea_results <- process_mapped_data(
            path_info$base_dir, path_info$file_pattern, experiment_name, 
            rank_method, min_size, max_size, model_type, suffix_mode, 
            mapping_type, debug
        )
    } else {
        gsea_results <- process_unmapped_data(
            path_info$base_dir, path_info$file_pattern, experiment_name, 
            rank_method, min_size, max_size, model_type, suffix_mode, 
            mapping_type, debug
        )
    }
    
    # 3. Save results
    save_gsea_results(gsea_results, experiment_name, model_type, suffix_mode, mapping_type, debug)
    
    return(gsea_results)
}

message("GSEA functions loaded successfully. Use run_gsea_analysis() or run_gsea_analysis_comprehensive() to perform gene set enrichment analysis.")
