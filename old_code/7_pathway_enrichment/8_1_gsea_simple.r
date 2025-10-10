#####
# Gene Set Enrichment Analysis Functions

# Required libraries
library(dplyr)      # For data manipulation
library(tidyr)      # For data tidying
library(msigdbr)    # For accessing MSigDB gene sets
library(fgsea)      # For running GSEA
library(ggplot2)    # For creating lollipop plots
library(stringr)    # For string manipulation
library(tibble)     # For tibble operations
library(purrr)      # For map functions
library(grid)       # For graphics
library(RColorBrewer) # For color palettes
library(ComplexHeatmap) # For legend
library(circlize)     # For circos plots

# Calculation function for ranking statistic
calculate_rank_stat <- function(data, method = "log2fc") {
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

# Main GSEA function
run_multiple_gsea <- function(diff_expr_data, 
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
  
  clean_data <- calculate_rank_stat(clean_data, rank_method)
  
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

# Lollipop plot functions
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

# Export function
export_gsea_results <- function(all_gsea_results, all_plots, experiment_name, contrast_name,
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

# Main GSEA wrapper function
run_gsea_analysis <- function(experiment_name,
                            model_type = "deseq2",
                            rank_method = "log2fc_log_pval",
                            min_size = 10,
                            max_size = 1000,
                            suffix_mode = NULL,
                            mapping_type = NULL,
                            debug = TRUE) {
    
    if(debug) cat("Starting GSEA analysis for", experiment_name, "using", model_type, "\n")
    
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
    
    # Initialize results list
    gsea_results <- list()
    
    if(is_mapped_data) {
        # For mapped data, process CSV files
        deg_files <- list.files(base_dir, 
                               pattern = file_pattern, 
                               full.names = TRUE)
        if(length(deg_files) == 0) {
            stop("No DEG files found in ", base_dir)
        }
        
        # Process each contrast
        for(deg_file in deg_files) {
            # Extract contrast name from filename
            contrast_name <- sub(paste0(experiment_name, "_"), "", basename(deg_file))
            contrast_name <- sub("_DEGs\\.csv$", "", contrast_name)
            
            if(debug) cat("\nProcessing contrast:", contrast_name, "\n")
            
            # Read the CSV file
            deg_data <- read.csv(deg_file)
            
            # Run GSEA for this contrast
            contrast_results <- run_multiple_gsea(
                diff_expr_data = deg_data,
                rank_method = rank_method,
                min_size = min_size,
                max_size = max_size
            )
            
            # Store results and create plots
            gsea_results[[contrast_name]] <- contrast_results
            lollipop_plots <- create_all_lollipop_plots(contrast_results$results)
            
            # Export results
            export_gsea_results(
                contrast_results,
                lollipop_plots,
                experiment_name = experiment_name,
                contrast_name = contrast_name,
                model_type = model_type,
                suffix_mode = suffix_mode,
                mapping_type = mapping_type
            )
        }
    } else {
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
            contrast_results <- run_multiple_gsea(
                diff_expr_data = deg_data,
                rank_method = rank_method,
                min_size = min_size,
                max_size = max_size
            )
            
            # Store results and create plots
            gsea_results[[contrast_name]] <- contrast_results
            lollipop_plots <- create_all_lollipop_plots(contrast_results$results)
            
            # Export results with sanitized name
            export_gsea_results(
                contrast_results,
                lollipop_plots,
                experiment_name = experiment_name,
                contrast_name = safe_contrast_name,  # Use sanitized name
                model_type = model_type,
                suffix_mode = suffix_mode,
                mapping_type = mapping_type
            )
        }
    }
    
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
    
    return(gsea_results)
}

# Example usage:
# For DESeq2 results with "highest_effect" suffix mode and clean mapping
xen_tran_2024_03_deseq2_gsea_results_he_clean <- run_gsea_analysis(
  experiment_name = "xen_tran_2024_03",
  model_type = "deseq2",
  rank_method = "sign_log_padj",  
  min_size = 5,
  max_size = 1000,
  suffix_mode = "highest_effect",
  mapping_type = "clean"
)

# For DESeq2 results with "highest_effect" suffix mode and messy mapping
xen_tran_2024_03_deseq2_gsea_results_he_messy <- run_gsea_analysis(
  experiment_name = "xen_tran_2024_03",
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000,
  suffix_mode = "highest_effect",
  mapping_type = "messy"  
)

# For DESeq2 results with "highest_effect" suffix mode and messy mapping
xen_tran_2024_03_deseq2_gsea_results_all_clean <- run_gsea_analysis(
  experiment_name = "xen_tran_2024_03",
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000,
  suffix_mode = "all",
  mapping_type = "clean"  
)

# For DESeq2 results with "highest_effect" suffix mode and messy mapping
xen_tran_2024_03_deseq2_gsea_results_all_messy <- run_gsea_analysis(
  experiment_name = "xen_tran_2024_03",
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000,
  suffix_mode = "all",
  mapping_type = "messy"  
)

gsea_results <- readRDS("results/xen_tran_2024_03/gsea/deseq2/all/messy/xen_tran_2024_03_deseq2_all_messy_gsea_results.rds")
str(gsea_results)




# Example usage:
# For DESeq2 results with "highest_effect" suffix mode and clean mapping
xen_tran_2024_12_deseq2_gsea_results_he_clean <- run_gsea_analysis(
  experiment_name = "xen_tran_2024_12",
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000,
  suffix_mode = "highest_effect",
  mapping_type = "clean" 
)

# For DESeq2 results with "highest_effect" suffix mode and messy mapping
xen_tran_2024_12_deseq2_gsea_results_he_messy <- run_gsea_analysis(
  experiment_name = "xen_tran_2024_12",
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000,
  suffix_mode = "highest_effect",
  mapping_type = "messy"  
)

# For DESeq2 results with "highest_effect" suffix mode and messy mapping
xen_tran_2024_12_deseq2_gsea_results_all_clean <- run_gsea_analysis(
  experiment_name = "xen_tran_2024_12",
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000,
  suffix_mode = "all",
  mapping_type = "clean"  
)

# For DESeq2 results with "highest_effect" suffix mode and messy mapping
xen_tran_2024_12_deseq2_gsea_results_all_messy <- run_gsea_analysis(
  experiment_name = "xen_tran_2024_12",
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000,
  suffix_mode = "all",
  mapping_type = "messy"  
)




# Example usage for unmapped human data:
hum_tran_2024_03_deep_phenotype_muscle_gsea_results <- run_gsea_analysis(
  experiment_name = "hum_tran_2024_03_deep_phenotype_muscle",
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000
  # No suffix_mode or mapping_type needed for human data
)

# Example usage for unmapped human PBMCs data:
hum_tran_2024_03_deep_phenotype_pbmcs_gsea_results <- run_gsea_analysis(
  experiment_name = "hum_tran_2024_03_deep_phenotype_pbmcs", 
  model_type = "deseq2",
  rank_method = "sign_log_padj",
  min_size = 10,
  max_size = 1000
  # No suffix_mode or mapping_type needed for human data
)

#gsea_results <- readRDS("results/xen_tran_2024_03/gsea/deseq2/highest_effect/clean/xen_tran_2024_03_deseq2_highest_effect_clean_gsea_results.rds")
#str(gsea_results)



