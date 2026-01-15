#' Gene Mapping Functions
#' 
#' Map Xenopus genes to human orthologs using Xenbase and HCOP databases.

# Required libraries
suppressPackageStartupMessages({
    library(dplyr)         # For data manipulation
    library(tidyr)         # For data tidying
    library(rentrez)       # For NCBI Entrez queries
    library(stringr)       # For string manipulation
    library(readr)         # For reading files
    library(org.Hs.eg.db)  # For human gene mapping
    library(AnnotationDbi) # For human gene mapping
    library(officer)       # For creating documents
    library(fs)            # For file system operations
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
mapping_config <- list(
    valid_suffix_modes = c("L", "S", "average", "highest_effect", "all"),
    valid_model_types = c("deseq2", "edgeR"),
    default_paths = list(
        xenbase_file = "pipeline_input_data_files/Xenbase_Xenopus_Orthology_Predictions.txt",
        hcop_file = "pipeline_input_data_files/HCOP_Xenopus_Orthology_Predictions.txt"
    )
)

#' Validate mapping parameters
#' @param experiment_name Name of experiment
#' @param model_type Type of model
#' @param suffix_mode Suffix mode for processing
validate_mapping_params <- function(experiment_name, model_type, suffix_mode) {
    if (is.null(experiment_name) || nchar(experiment_name) == 0) {
        stop("experiment_name cannot be empty")
    }
    
    if (!model_type %in% mapping_config$valid_model_types) {
        stop("Invalid model_type. Must be one of: ", paste(mapping_config$valid_model_types, collapse = ", "))
    }
    
    if (!suffix_mode %in% mapping_config$valid_suffix_modes) {
        stop("Invalid suffix_mode. Must be one of: ", paste(mapping_config$valid_suffix_modes, collapse = ", "))
    }
}

#' Load orthology mapping files with validation
#' @param xenbase_file Path to Xenbase file
#' @param hcop_file Path to HCOP file
#' @return List containing mapping data
load_orthology_files <- function(xenbase_file = NULL, hcop_file = NULL) {
    # Use default paths if not provided
    if (is.null(xenbase_file)) {
        xenbase_file <- here(mapping_config$default_paths$xenbase_file)
    }
    if (is.null(hcop_file)) {
        hcop_file <- here(mapping_config$default_paths$hcop_file)
    }
    
    # Validate file existence
    if (!file.exists(xenbase_file)) {
        stop("Xenbase file not found: ", xenbase_file)
    }
    if (!file.exists(hcop_file)) {
        stop("HCOP file not found: ", hcop_file)
    }
    
    tryCatch({
        # Load Xenbase data
        xenbase_data <- read_tsv(xenbase_file, show_col_types = FALSE)
        message("Loaded Xenbase data: ", nrow(xenbase_data), " rows")
        
        # Load HCOP data
        hcop_data <- read_tsv(hcop_file, show_col_types = FALSE)
        message("Loaded HCOP data: ", nrow(hcop_data), " rows")
        
        return(list(
            xenbase = xenbase_data,
            hcop = hcop_data
        ))
    }, error = function(e) {
        stop("Failed to load orthology files: ", e$message)
    })
}

#' Prepare initial gene data structures
#' @param xenopus_genes Vector of Xenopus gene symbols
#' @return List containing original symbols and unique gene list
prepare_gene_data <- function(xenopus_genes) {
    # Create initial dataframes and lookups
    original_symbols_df <- data.frame(
        original_symbol = xenopus_genes,
        gene = remove_suffix(xenopus_genes)
    ) %>% distinct()
    
    # Create unique gene list for mapping
    gene_list_no_suffix <- remove_suffix(xenopus_genes)
    Xenopus_gene_list_unique <- data.frame(gene = unique(gene_list_no_suffix))
    
    return(list(
        original_symbols = original_symbols_df,
        unique_genes = Xenopus_gene_list_unique
    ))
}

#' Load and standardize orthology reference data
#' @param xenbase_file Path to Xenbase file
#' @param hcop_file Path to HCOP file
#' @return List containing standardized Xenbase and HCOP data
load_and_standardize_data <- function(xenbase_file, hcop_file) {
    # Load reference data
    xenbase_data <- read.delim(xenbase_file,
        header = TRUE, sep = "\t",
        quote = "", fill = TRUE
    )

    hcop_data <- read.delim(hcop_file,
        header = TRUE, sep = "\t",
        quote = "", fill = TRUE
    )

    # Standardize gene names (lowercase and trim)
    xenbase_data$Xenbase_gene_symbol <- tolower(trimws(xenbase_data$Xenbase_gene_symbol))
    hcop_data$xenopus_symbol <- tolower(trimws(hcop_data$xenopus_symbol))

    # Select relevant columns
    xenbase_data <- xenbase_data %>%
        dplyr::select(Xenbase_gene_symbol, Human_NCBI_Entrez_ID)

    hcop_data <- hcop_data %>%
        dplyr::select(xenopus_symbol, human_entrez_gene, human_ensembl_gene)
    
    return(list(
        xenbase = xenbase_data,
        hcop = hcop_data
    ))
}

#' Perform initial mapping to Xenbase and HCOP databases
#' @param unique_genes Unique gene list
#' @param xenbase_data Standardized Xenbase data
#' @param hcop_data Standardized HCOP data
#' @return List containing initial mappings and ID mappings
perform_initial_mappings <- function(unique_genes, xenbase_data, hcop_data) {
    # Standardize gene names
    unique_genes$gene <- tolower(trimws(unique_genes$gene))
    
    # Xenbase matching sequence
    xenbase_initial <- unique_genes %>%
        left_join(xenbase_data, by = c("gene" = "Xenbase_gene_symbol")) %>%
        filter(!is.na(Human_NCBI_Entrez_ID)) %>%
        mutate(
            human_entrez_id = as.character(Human_NCBI_Entrez_ID),
            source = "xenbase"
        ) %>%
        dplyr::select(-Human_NCBI_Entrez_ID) %>%
        distinct()

    # HCOP matching sequence
    hcop_initial <- unique_genes %>%
        left_join(hcop_data, by = c("gene" = "xenopus_symbol")) %>%
        filter(
            (!is.na(human_entrez_gene) & human_entrez_gene != "-") |
                (!is.na(human_ensembl_gene) & human_ensembl_gene != "-")
        ) %>%
        mutate(
            human_entrez_id = ifelse(human_entrez_gene == "-", NA, human_entrez_gene),
            human_ensembl_id = ifelse(human_ensembl_gene == "-", NA, human_ensembl_gene),
            source = "hcop"
        ) %>%
        dplyr::select(-human_entrez_gene, -human_ensembl_gene) %>%
        distinct()

    # Create ID mappings for duplicate analysis
    xenbase_id_mappings <- xenbase_initial %>%
        group_by(human_entrez_id) %>%
        summarise(
            xenbase_genes = paste(unique(gene), collapse = ";"),
            n_xenopus_genes = n_distinct(gene)
        )

    hcop_id_mappings <- hcop_initial %>%
        group_by(human_entrez_id, human_ensembl_id) %>%
        summarise(
            hcop_genes = paste(unique(gene), collapse = ";"),
            n_xenopus_genes = n_distinct(gene)
        )
    
    return(list(
        xenbase_initial = xenbase_initial,
        hcop_initial = hcop_initial,
        xenbase_id_mappings = xenbase_id_mappings,
        hcop_id_mappings = hcop_id_mappings
    ))
}

#' Print diagnostic information about gene mapping
#' @param original_symbols Original symbols data frame
#' @param unique_genes Unique genes data frame
#' @param xenopus_genes Original Xenopus genes vector
#' @param mappings Initial mapping results
#' @param debug Debug flag
print_mapping_diagnostics <- function(original_symbols, unique_genes, xenopus_genes, mappings, debug) {
    if (!debug) return()
    
    # Count genes with different suffix patterns
    suffix_patterns <- data.frame(
        original_symbol = original_symbols$original_symbol,
        stringsAsFactors = FALSE
    ) %>%
        mutate(
            has_L = grepl("\\.L$", original_symbol),
            has_S = grepl("\\.S$", original_symbol),
            has_neither = !has_L & !has_S
        ) %>%
        summarise(
            L_count = sum(has_L, na.rm = TRUE),
            S_count = sum(has_S, na.rm = TRUE),
            neither_count = sum(has_neither, na.rm = TRUE)
        )

    # Calculate pairs by looking at base genes
    pair_analysis <- data.frame(
        original_symbol = original_symbols$original_symbol,
        stringsAsFactors = FALSE
    ) %>%
        mutate(
            base_gene = remove_suffix(original_symbol),
            has_L = grepl("\\.L$", original_symbol),
            has_S = grepl("\\.S$", original_symbol)
        ) %>%
        group_by(base_gene) %>%
        summarise(
            has_L = any(has_L, na.rm = TRUE),
            has_S = any(has_S, na.rm = TRUE),
            is_pair = has_L & has_S,
            .groups = "drop"
        ) %>%
        summarise(
            total_pairs = sum(is_pair, na.rm = TRUE),
            L_only = sum(has_L & !has_S, na.rm = TRUE),
            S_only = sum(!has_L & has_S, na.rm = TRUE),
            neither = sum(!has_L & !has_S, na.rm = TRUE)
        )

    cat("\nOriginal Symbols Diagnostics:\n")
    cat("\nGene Counts:")
    cat("\n  All Initial Genes:", length(xenopus_genes))
    cat("\n  All Initial Unique Genes:", nrow(unique_genes))
    cat("\nSuffix counts:")
    cat("\n.L genes:", suffix_patterns$L_count)
    cat("\n.S genes:", suffix_patterns$S_count)
    cat("\nNeither .L nor .S:", suffix_patterns$neither_count)
    cat("\nComplete L/S pairs:", pair_analysis$total_pairs)
    cat("\nL only (no S partner):", pair_analysis$L_only)
    cat("\nS only (no L partner):", pair_analysis$S_only)
    cat("\nID Mapping Analysis:")
    cat("\nXenbase:")
    cat("\n  Total Unique Entrez IDs:", n_distinct(na.omit(mappings$xenbase_initial$human_entrez_id)))
    cat("\n  Total IDs:", length(na.omit(mappings$xenbase_initial$human_entrez_id)))
    cat("\n  Entrez IDs mapping to multiple Xenopus genes (Case 2 duplicates):", sum(mappings$xenbase_id_mappings$n_xenopus_genes > 1))
    cat("\n  Max Xenopus genes per Entrez ID:", max(mappings$xenbase_id_mappings$n_xenopus_genes))
    cat("\n\nHCOP:")
    cat("\n  Total Unique Entrez IDs:", n_distinct(na.omit(mappings$hcop_initial$human_entrez_id)))
    cat("\n  Total Entrez IDs:", length(na.omit(mappings$hcop_initial$human_entrez_id)))
    cat("\n  Total Unique Ensembl IDs:", n_distinct(na.omit(mappings$hcop_initial$human_ensembl_id)))
    cat("\n  Total Ensembl IDs:", length(na.omit(mappings$hcop_initial$human_ensembl_id)))
    cat("\n  IDs mapping to multiple Xenopus genes (Case 2 duplicates):", sum(mappings$hcop_id_mappings$n_xenopus_genes > 1))
    cat("\n  Max Xenopus genes per ID:", max(mappings$hcop_id_mappings$n_xenopus_genes))
}

#' Get HUGO symbols for human gene IDs
#' @param xenbase_initial Xenbase initial mappings
#' @param hcop_initial HCOP initial mappings
#' @return List containing HUGO symbol mappings
get_hugo_symbols <- function(xenbase_initial, hcop_initial) {
    # Get HUGO symbols for all potential human IDs
    hugo_from_entrez <- tryCatch(
        {
            suppressWarnings(
                AnnotationDbi::select(org.Hs.eg.db,
                    keys = na.omit(c(xenbase_initial$human_entrez_id, hcop_initial$human_entrez_id)),
                    columns = c("SYMBOL", "ENTREZID"),
                    keytype = "ENTREZID"
                )
            ) %>%
                dplyr::rename(
                    hugo_symbol_entrez = SYMBOL,
                    human_entrez_id = ENTREZID
                )
        },
        error = function(e) {
            cat("\nError mapping Entrez IDs:", e$message)
            data.frame(
                hugo_symbol_entrez = character(0),
                human_entrez_id = character(0)
            )
        }
    )

    hugo_from_ensembl <- tryCatch(
        {
            suppressWarnings(
                AnnotationDbi::select(org.Hs.eg.db,
                    keys = na.omit(hcop_initial$human_ensembl_id),
                    columns = c("SYMBOL", "ENSEMBL"),
                    keytype = "ENSEMBL"
                )
            ) %>%
                dplyr::rename(
                    hugo_symbol_ensembl = SYMBOL,
                    human_ensembl_id = ENSEMBL
                )
        },
        error = function(e) {
            cat("\nError mapping Ensembl IDs:", e$message)
            data.frame(
                hugo_symbol_ensembl = character(0),
                human_ensembl_id = character(0)
            )
        }
    )
    
    return(list(
        hugo_from_entrez = hugo_from_entrez,
        hugo_from_ensembl = hugo_from_ensembl
    ))
}

#' Merge and process gene mappings with HUGO symbols
#' @param xenbase_initial Xenbase initial mappings
#' @param hcop_initial HCOP initial mappings
#' @param hugo_symbols HUGO symbol mappings
#' @return Processed gene mapping data frame
merge_and_process_mappings <- function(xenbase_initial, hcop_initial, hugo_symbols) {
    # Merge Xenbase and HCOP data
    # Ensure consistent data types before combining
    xenbase_initial$source <- "xenbase"
    hcop_initial$source <- "hcop"
    
    # Ensure human_entrez_id is character in both data frames
    xenbase_initial$human_entrez_id <- as.character(xenbase_initial$human_entrez_id)
    hcop_initial$human_entrez_id <- as.character(hcop_initial$human_entrez_id)
    
    xenbase_hcop_merge <- bind_rows(xenbase_initial, hcop_initial) %>%
        distinct()

    # Add HUGO symbols
    entrez_xenbase_hcop_hugo_symbol_merge <- xenbase_hcop_merge %>%
        left_join(hugo_symbols$hugo_from_entrez,
            by = "human_entrez_id",
            relationship = "many-to-many"
        ) %>%
        distinct()

    entrez_ensemble_xenbase_hcop_hugo_symbol_merge <- entrez_xenbase_hcop_hugo_symbol_merge %>%
        left_join(hugo_symbols$hugo_from_ensembl,
            by = "human_ensembl_id",
            relationship = "many-to-many"
        ) %>%
        distinct()

    # Create final gene mapping data frame
    gene_entrez_hugo_df <- entrez_ensemble_xenbase_hcop_hugo_symbol_merge %>%
        # create hugo_symbol column with values from hugo_symbol_entrez or hugo_symbol_ensembl if entrez is NA
        mutate(hugo_symbol = ifelse(is.na(hugo_symbol_entrez), hugo_symbol_ensembl, hugo_symbol_entrez)) %>%
        dplyr::select(gene, human_entrez_id, source, hugo_symbol) %>%
        distinct()

    # Determine source overlap
    xenbase_list <- gene_entrez_hugo_df %>%
        dplyr::filter(source == "xenbase") %>%
        dplyr::select(gene) %>%
        distinct() %>%
        pull(gene)

    hcop_list <- gene_entrez_hugo_df %>%
        dplyr::filter(source == "hcop") %>%
        dplyr::select(gene) %>%
        distinct() %>%
        pull(gene)

    # Update source information
    gene_entrez_hugo_df <- gene_entrez_hugo_df %>%
        mutate(
            xenbase_gene = gene %in% xenbase_list,
            hcop_gene = gene %in% hcop_list
        ) %>%
        mutate(
            source = case_when(
                xenbase_gene & hcop_gene ~ "both",
                xenbase_gene ~ "xenbase",
                hcop_gene ~ "hcop",
                TRUE ~ "none"
            )
        ) %>%
        dplyr::select(-xenbase_gene, -hcop_gene)
    
    return(gene_entrez_hugo_df)
}

#' Handle duplicate mappings and create clean data
#' @param gene_entrez_hugo_df Gene mapping data frame
#' @return Clean mapping data frame
handle_duplicates_and_clean <- function(gene_entrez_hugo_df) {
    # Identify duplicate types
    gene_entrez_hugo_df_dupes <- gene_entrez_hugo_df %>%
        group_by(gene) %>%
        mutate(
            case_1_1 = n() > 1 & n_distinct(hugo_symbol) == 1,
            case_1_2 = n() > 1 & n_distinct(hugo_symbol) > 1
        ) %>%
        ungroup()
    
    # Report duplicate info
    cat(sprintf("Total rows: %d", nrow(gene_entrez_hugo_df)))
    cat(sprintf("\nTotal duplicate rows: %d", sum(gene_entrez_hugo_df_dupes$case_1_1 | gene_entrez_hugo_df_dupes$case_1_2)))
    cat(sprintf("\nCase 1.1 (multiple rows, same ID & HUGO): %d rows", sum(gene_entrez_hugo_df_dupes$case_1_1)))
    cat(sprintf("\nCase 1.2 (different IDs): %d rows", sum(gene_entrez_hugo_df_dupes$case_1_2)))

    # Clean duplicates
    merged_data_clean <- gene_entrez_hugo_df_dupes %>%
        group_by(gene) %>%
        # for each gene, if case 1.1, keep the first row
        filter(!case_1_2) %>%
        filter((!case_1_1) |
            (case_1_1 & row_number() == 1)) %>%
        # if case 1.2, drop all rows for that gene
        ungroup() %>%
        mutate(gene_symbol = gene)

    cat(sprintf("Total rows after handling Case 1 duplicates: %d", nrow(merged_data_clean)))
    
    return(merged_data_clean)
}

#' Enhance original symbols with duplicate type information
#' @param original_symbols Original symbols data frame
#' @param merged_data_clean Clean mapping data
#' @param gene_entrez_hugo_df Original gene mapping data
#' @return Enhanced original symbols data frame
enhance_original_symbols <- function(original_symbols, merged_data_clean, gene_entrez_hugo_df) {
    # Identify duplicate types
    gene_entrez_hugo_df_dupes <- gene_entrez_hugo_df %>%
        group_by(gene) %>%
        mutate(
            case_1_1 = n() > 1 & n_distinct(hugo_symbol) == 1,
            case_1_2 = n() > 1 & n_distinct(hugo_symbol) > 1
        ) %>%
        ungroup()
        
    # Identify case 2 duplicates (multiple Xenopus genes to same human gene)
    case_2_duplicates <- merged_data_clean %>%
        group_by(hugo_symbol) %>%
        filter(n() > 1) %>%
        pull(hugo_symbol) %>%
        unique()
    
    # Enhance original symbols with duplicate type information
    enhanced_symbols <- original_symbols %>%
        left_join(
            merged_data_clean %>% dplyr::select(gene, hugo_symbol),
            by = "gene"
        ) %>%
        # Left join with the duplicate information
        left_join(
            gene_entrez_hugo_df_dupes %>% 
                dplyr::select(gene, case_1_1, case_1_2) %>%
                distinct(gene, .keep_all = TRUE),
            by = "gene"
        ) %>%
        # Add case 2 duplicate information
        mutate(
            case_2 = !is.na(hugo_symbol) & hugo_symbol %in% case_2_duplicates,
            # Fill NA values with FALSE
            case_1_1 = ifelse(is.na(case_1_1), FALSE, case_1_1),
            case_1_2 = ifelse(is.na(case_1_2), FALSE, case_1_2),
            # Add a simplified duplicate_type column
            duplicate_type = case_when(
                case_1_2 ~ "case_1_2", # Conflicting human mappings
                case_1_1 ~ "case_1_1", # Same human mapping
                case_2 ~ "case_2",     # Multiple Xenopus genes to same human
                !is.na(hugo_symbol) ~ "clean", # Has mapping but no duplicates
                TRUE ~ "unmapped"     # No human mapping
            )
        )
    
    return(enhanced_symbols)
}

#' Print overlap analysis
#' @param merged_data_clean Clean mapping data
#' @param debug Debug flag
print_overlap_analysis <- function(merged_data_clean, debug) {
    if (!debug) return()
    
    cat("\nOverlap Analysis:")
    cat("\nTotal overlap between databases:", nrow(merged_data_clean %>% filter(source == "both")), "genes")
    cat("\n\nUnique matches:")
    cat("\nXenbase only:", nrow(merged_data_clean %>% filter(source == "xenbase")), "genes")
    cat("\nHCOP only:", nrow(merged_data_clean %>% filter(source == "hcop")), "genes")
}

#' Map Xenopus genes to human orthologs
#' @param xenopus_genes Vector of Xenopus gene symbols
#' @param xenbase_file Path to Xenbase file
#' @param hcop_file Path to HCOP file
#' @param suffix_mode Suffix mode for processing
#' @param de_results DE results (optional)
#' @param experiment_name Experiment name (optional)
#' @param model_type Model type (optional)
#' @param debug Debug flag
#' @return List containing mapping results and original symbols
map_xenopus_to_human <- function(xenopus_genes,
                                 xenbase_file,
                                 hcop_file,
                                 suffix_mode = "all",
                                 de_results = NULL,
                                 experiment_name = NULL,
                                 model_type = NULL,
                                 debug = TRUE) {
    # Validate suffix_mode
    valid_modes <- c("all", "L", "S", "average", "highest_effect")
    if (!suffix_mode %in% valid_modes) {
        stop("Invalid suffix_mode. Must be one of: ", paste(valid_modes, collapse = ", "))
    }

    # 1. Prepare initial gene data structures
    gene_data <- prepare_gene_data(xenopus_genes)
    
    # 2. Load and standardize reference data
    ref_data <- load_and_standardize_data(xenbase_file, hcop_file)
    
    # 3. Perform initial mappings
    mappings <- perform_initial_mappings(gene_data$unique_genes, ref_data$xenbase, ref_data$hcop)
    
    # 4. Print diagnostics
    print_mapping_diagnostics(gene_data$original_symbols, gene_data$unique_genes, 
                             xenopus_genes, mappings, debug)
    
    # 5. Get HUGO symbols
    hugo_symbols <- get_hugo_symbols(mappings$xenbase_initial, mappings$hcop_initial)
    
    # 6. Merge and process mappings
    gene_entrez_hugo_df <- merge_and_process_mappings(mappings$xenbase_initial, 
                                                     mappings$hcop_initial, 
                                                     hugo_symbols)
    
    # 7. Handle duplicates and create clean data
    merged_data_clean <- handle_duplicates_and_clean(gene_entrez_hugo_df)
    
    # 8. Print overlap analysis
    print_overlap_analysis(merged_data_clean, debug)
    
    # 9. Enhance original symbols with duplicate type information
    enhanced_symbols <- enhance_original_symbols(gene_data$original_symbols, 
                                                merged_data_clean, 
                                                gene_entrez_hugo_df)
    
    # Save results if experiment_name and model_type are provided
    if (!is.null(experiment_name) && !is.null(model_type)) {
        tryCatch({
            results_dir <- here("results", experiment_name, "mapping", model_type, suffix_mode)
            if (!dir.exists(results_dir)) {
                dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
            }
            
            # Save the complete results
            results_file <- file.path(results_dir, paste0(experiment_name, "_", suffix_mode, "_", model_type, "_results.rds"))
            saveRDS(list(
                data = merged_data_clean,
                original_symbols = enhanced_symbols
            ), results_file)
            
            # Save just the mapping data as CSV for easy viewing
            if (nrow(merged_data_clean) > 0) {
                csv_file <- file.path(results_dir, paste0(experiment_name, "_", suffix_mode, "_", model_type, "_mapping.csv"))
                write.csv(merged_data_clean, csv_file, row.names = FALSE)
            }
            
            message("Mapping results saved to: ", results_dir)
        }, error = function(e) {
            warning("Failed to save mapping results: ", e$message)
        })
    }
    
    # Return results
    return(list(
        data = merged_data_clean,
        original_symbols = enhanced_symbols
    ))
}

#' Process gene symbols based on suffix mode
#' @param gene_symbols Vector of gene symbols
#' @param suffix_mode Suffix mode
#' @return Processed gene symbols
process_gene_symbols <- function(gene_symbols, suffix_mode) {
    switch(suffix_mode,
        "L" = gene_symbols[grepl("\\.L$", gene_symbols)],
        "S" = gene_symbols[grepl("\\.S$", gene_symbols)],
        "average" = {
            # Remove suffixes and get unique genes
            base_genes <- gsub("\\.(L|S)$", "", gene_symbols)
            unique(base_genes)
        },
        "highest_effect" = {
            # This would need additional logic based on fold changes
            # For now, return all genes
            gene_symbols
        },
        "all" = gene_symbols,
        stop("Invalid suffix_mode: ", suffix_mode)
    )
}

#' Map genes using Xenbase database
#' @param gene_symbols Vector of gene symbols
#' @param xenbase_data Xenbase data frame
#' @return Mapping results
map_with_xenbase <- function(gene_symbols, xenbase_data) {
    tryCatch({
        # This is a simplified version - actual implementation would depend on Xenbase format
        # Assuming columns: Xenopus_gene, Human_gene, Confidence
        if (ncol(xenbase_data) < 2) {
            warning("Xenbase data format unexpected")
            return(data.frame())
        }
        
        # Perform mapping (simplified)
        mapping <- xenbase_data %>%
            filter(Xenopus_gene %in% gene_symbols) %>%
            select(gene = Xenopus_gene, hugo_symbol = Human_gene, 
                   confidence = if("Confidence" %in% names(.)) Confidence else NA) %>%
            mutate(source = "xenbase")
        
        return(mapping)
    }, error = function(e) {
        warning("Failed to map with Xenbase: ", e$message)
        return(data.frame())
    })
}

#' Map genes using HCOP database
#' @param gene_symbols Vector of gene symbols
#' @param hcop_data HCOP data frame
#' @return Mapping results
map_with_hcop <- function(gene_symbols, hcop_data) {
    tryCatch({
        # This is a simplified version - actual implementation would depend on HCOP format
        if (ncol(hcop_data) < 2) {
            warning("HCOP data format unexpected")
            return(data.frame())
        }
        
        # Perform mapping (simplified)
        mapping <- hcop_data %>%
            filter(Xenopus_gene %in% gene_symbols) %>%
            select(gene = Xenopus_gene, hugo_symbol = Human_gene,
                   confidence = if("Confidence" %in% names(.)) Confidence else NA) %>%
            mutate(source = "hcop")
        
        return(mapping)
    }, error = function(e) {
        warning("Failed to map with HCOP: ", e$message)
        return(data.frame())
    })
}

#' Combine mappings from different sources
#' @param xenbase_mapping Xenbase mapping results
#' @param hcop_mapping HCOP mapping results
#' @return Combined mapping results
combine_mappings <- function(xenbase_mapping, hcop_mapping) {
    tryCatch({
        # Combine mappings
        combined <- bind_rows(xenbase_mapping, hcop_mapping)
        
        if (nrow(combined) == 0) {
            return(data.frame())
        }
        
        # Remove duplicates and prioritize by source
        combined <- combined %>%
            group_by(gene, hugo_symbol) %>%
            summarise(
                source = paste(unique(source), collapse = ","),
                confidence = max(confidence, na.rm = TRUE),
                .groups = "drop"
            )
        
        return(combined)
    }, error = function(e) {
        stop("Failed to combine mappings: ", e$message)
    })
}

#' Analyze L/S divergence in gene expression
#' @param data Mapping data
#' @param de_results Differential expression results
#' @param debug Debug flag
#' @return Divergence analysis results
analyze_ls_divergence <- function(data, de_results, debug = TRUE) {
    # Process each contrast
    ls_divergence <- lapply(de_results, function(result) {
        de_result <- result$deseq_results
        de_df <- data.frame(
            gene_symbol = rownames(de_result),
            log2FoldChange = de_result$log2FoldChange,
            padj = de_result$padj,
            stringsAsFactors = FALSE
        ) %>%
            mutate(base_gene = remove_suffix(gene_symbol))

        # Find pairs with both L and S
        ls_pairs <- de_df %>%
            group_by(base_gene) %>%
            filter(
                n() == 2,
                any(grepl("\\.L$", gene_symbol)),
                any(grepl("\\.S$", gene_symbol))
            ) %>%
            summarise(
                gene_symbols = paste(sort(gene_symbol), collapse = "/"),
                fold_changes = paste(round(log2FoldChange, 2), collapse = "/"),
                padj_values = paste(round(padj, 4), collapse = "/"),
                direction_conflict = sign(log2FoldChange[1]) != sign(log2FoldChange[2]),
                fc_diff = abs(diff(log2FoldChange)),
                padj_diff = abs(diff(padj)),
                both_significant = all(padj < 0.05),
                .groups = "drop"
            )
        
        # Calculate divergence stats
        stats <- list(
            total_pairs = nrow(ls_pairs),
            direction_conflicts = sum(ls_pairs$direction_conflict),
            large_fc_diff = sum(ls_pairs$fc_diff > 1),
            large_padj_diff = sum(ls_pairs$padj_diff > 0.1),
            both_large_diff = sum(ls_pairs$fc_diff > 1 & ls_pairs$padj_diff > 0.1),
            sig_direction_conflicts = sum(ls_pairs$direction_conflict & ls_pairs$both_significant)
        )
        
        list(pairs = ls_pairs, stats = stats)
    })

    return(ls_divergence)
}

#' Process raw counts based on mode
#' @param counts_df Counts data frame
#' @param mapping_output Mapping results
#' @param suffix_mode Suffix mode
#' @return Processed counts data
process_raw_counts <- function(counts_df, mapping_output, suffix_mode) {
    if (is.null(counts_df) || nrow(counts_df) == 0) {
        stop("counts_df cannot be empty")
    }
    
    if (is.null(mapping_output) || is.null(mapping_output$data)) {
        stop("mapping_output must contain data")
    }
    
    tryCatch({
        # Create basic messy version with human mapping
        counts_messy <- counts_df %>%
            left_join(
                mapping_output$data,
                by = c("gene_symbol" = "gene")
            ) %>%
            filter(!is.na(hugo_symbol))
        
        # Create basic clean version
        counts_clean <- counts_messy %>%
            group_by(hugo_symbol) %>%
            filter(n_distinct(gene_symbol) == 1) %>%
            ungroup()
        
        # Process based on mode
        processed_counts <- switch(suffix_mode,
            "L" = filter(counts_messy, grepl("\\.L$", gene_symbol)),
            "S" = filter(counts_messy, grepl("\\.S$", gene_symbol)),
            "average" = {
                counts_messy %>%
                    group_by(hugo_symbol) %>%
                    summarise(across(where(is.numeric), mean), .groups = "drop") %>%
                    mutate(gene_symbol = hugo_symbol)
            },
            "highest_effect" = {
                # This would need additional logic based on effect sizes
                counts_clean
            },
            "all" = counts_messy,
            stop("Invalid suffix_mode: ", suffix_mode)
        )
        
        return(processed_counts)
    }, error = function(e) {
        stop("Failed to process raw counts: ", e$message)
    })
}


#' Load differential expression results
#' @param experiment_name Name of experiment
#' @param model_type Type of model
#' @return DE results data frame
load_de_results <- function(experiment_name, model_type) {
    results_file <- here("results", experiment_name, "de", model_type, "simple", 
                        paste0(experiment_name, "_de_results.rds"))
    
    if (!file.exists(results_file)) {
        stop("DE results file not found: ", results_file)
    }
    
    tryCatch({
        de_results <- readRDS(results_file)
        
        # Extract gene symbols from results
        if (is.list(de_results) && length(de_results) > 0) {
            # Get first contrast results
            first_result <- de_results[[1]]
            if (!is.null(first_result$results)) {
                gene_symbols <- rownames(first_result$results)
                return(data.frame(gene_symbol = gene_symbols))
            }
        }
        
        stop("Invalid DE results format")
    }, error = function(e) {
        stop("Failed to load DE results: ", e$message)
    })
}

#' Load raw counts data
#' @param experiment_name Name of experiment
#' @return Raw counts data frame
load_raw_counts <- function(experiment_name) {
    counts_file <- here("data", experiment_name, "cleaned_data", 
                       paste0(experiment_name, "_cleaned.rds"))
    
    if (!file.exists(counts_file)) {
        stop("Cleaned data file not found: ", counts_file)
    }
    
    tryCatch({
        cleaned_data <- readRDS(counts_file)
        
        if (!is.list(cleaned_data) || !"counts" %in% names(cleaned_data)) {
            stop("Invalid cleaned data format")
        }
        
        # Convert counts matrix to data frame
        counts_df <- as.data.frame(cleaned_data$counts)
        counts_df$gene_symbol <- rownames(counts_df)
        
        return(counts_df)
    }, error = function(e) {
        stop("Failed to load raw counts: ", e$message)
    })
}

#' Save mapping results
#' @param output Output data
#' @param experiment_name Experiment name
#' @param model_type Model type
#' @param suffix_mode Suffix mode
save_mapping_results <- function(output, experiment_name, model_type, suffix_mode) {
    results_dir <- here("results", experiment_name, "mapping", model_type, suffix_mode)
    
    tryCatch({
        if (!dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)) {
            stop("Failed to create results directory: ", results_dir)
        }
        
        # Save main results
        results_file <- file.path(results_dir, paste0(experiment_name, "_", suffix_mode, "_results.rds"))
        saveRDS(output, results_file)
        
        # Save mapping data separately
        if (!is.null(output$mapping_results)) {
            mapping_file <- file.path(results_dir, paste0(experiment_name, "_mapping_data.rds"))
            saveRDS(output$mapping_results, mapping_file)
        }
        
        message("Results saved to: ", results_dir)
    }, error = function(e) {
        stop("Failed to save mapping results: ", e$message)
    })
}

#' Run and capture output for a specific mode
#' @param experiment_name Name of experiment
#' @param model_type Type of model
#' @param suffix_mode Suffix mode
#' @param output_file Output file path
#' @return Analysis results
run_and_capture_mode <- function(experiment_name, model_type, suffix_mode, output_file) {
    validate_mapping_params(experiment_name, model_type, suffix_mode)
    
    message(paste0("\n\n========== STARTING ", suffix_mode, " MODE FOR ", experiment_name, " ==========\n"))
    
    # Create output directory
    output_dir <- dirname(output_file)
    if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Failed to create output directory: ", output_dir)
    }
    
    # Create a temporary text file for capturing output
    temp_file <- tempfile(fileext = ".txt")
    
    # Start capturing console output
    sink(temp_file)
    
    # Run the analysis
    result <- tryCatch({
        format_for_drughub(
            experiment_name = experiment_name,
            model_type = model_type,
            suffix_mode = suffix_mode
        )
    }, finally = {
        # Stop capturing
        sink(NULL)
    })
    
    # Read the captured output
    output_text <- suppressWarnings(readLines(temp_file))
    
    # Create a Word document
    tryCatch({
        doc <- officer::read_docx()
        doc <- officer::body_add_par(doc, paste0("Analysis for ", suffix_mode, " mode"), style = "heading 1")
        
        # Add the captured text
        for (line in output_text) {
            doc <- officer::body_add_par(doc, line, style = "Normal")
        }
        
        # Save the document
        print(doc, target = output_file)
        message("Analysis document saved to: ", output_file)
    }, error = function(e) {
        warning("Failed to create analysis document: ", e$message)
    })
    
    # Clean up
    unlink(temp_file)
    
    return(result)
}

# Additional functions from original implementation

#' Remove suffix function
remove_suffix <- function(x) {
    sub("\\.[LS]$", "", x)
}

#' Process DE results for DrugHub formatting
#' @param de_result DE result object
#' @param model_type Type of model
#' @param mapping_output Mapping results
#' @param suffix_mode Suffix mode
#' @param debug Debug flag
#' @return Processed DE results
process_de_results <- function(de_result, model_type, mapping_output, suffix_mode, debug = TRUE) {
    # Get base results based on model type
    df <- switch(model_type,
        "deseq2" = {
            data.frame(
                gene_symbol = remove_suffix(rownames(de_result$deseq_results)),
                log2FoldChange = de_result$deseq_results$log2FoldChange,
                pvalue = de_result$deseq_results$pvalue,
                padj = de_result$deseq_results$padj,
                stringsAsFactors = FALSE
            ) %>%
                left_join(
                    mapping_output$data %>%
                        dplyr::select(gene_symbol, hugo_symbol),
                    by = "gene_symbol",
                    relationship = "many-to-many"
                ) %>%
                filter(!is.na(hugo_symbol)) %>%
                mutate(
                    log2FoldChange = round(log2FoldChange, digits = 3)
                )
        }
    )
    
    # Count significant DEGs in original results
    orig_sig_count <- sum(de_result$deseq_results$padj < 0.05, na.rm = TRUE)
    
    # Count significant DEGs in mapped results
    mapped_sig_count <- sum(df$padj < 0.05, na.rm = TRUE)

    # Create clean version - properly removing all Case 2 duplicates
    # First identify which hugo symbols appear multiple times
    duplicate_hugos <- df %>%
        group_by(hugo_symbol) %>%
        summarise(count = n()) %>%
        filter(count > 1) %>%
        pull(hugo_symbol)

    # Then create clean version by excluding all rows with those hugo symbols
    df_clean <- df %>%
        filter(!hugo_symbol %in% duplicate_hugos)
        
    # Count significant DEGs in clean version
    clean_sig_count <- sum(df_clean$padj < 0.05, na.rm = TRUE)

    if (debug) {
        cat("\nDEA Results Analysis:")
        cat("\n--------------------")
        
        # Add DEG retention stats
        cat("\nSignificant DEG Retention:")
        cat(sprintf("\nOriginal significant DEGs: %d", orig_sig_count))
        cat(sprintf("\nMapped significant DEGs (messy): %d (%.1f%%)", 
                   mapped_sig_count, 100 * mapped_sig_count / orig_sig_count))
        cat(sprintf("\nMapped significant DEGs (clean): %d (%.1f%%)", 
                   clean_sig_count, 100 * clean_sig_count / orig_sig_count))
        
        # Analyze messy version
        cat("\nMESSY VERSION:")
        cat(sprintf("\nTotal rows: %d", nrow(df)))
        cat(sprintf("\nUnique hugo symbols: %d", n_distinct(df$hugo_symbol)))
        
        # Analyze clean version
        cat("\n\nCLEAN VERSION:")
        cat(sprintf("\nTotal rows: %d", nrow(df_clean)))
        cat(sprintf("\nUnique hugo symbols: %d", n_distinct(df_clean$hugo_symbol)))
        
        # Print Case 2 duplicate stats
        cat("\n\nCase 2 Duplicate Analysis:")
        cat(sprintf("\nHUGO symbols appearing multiple times (removed from clean): %d", length(duplicate_hugos)))
        cat(sprintf("\nTotal rows removed due to Case 2 duplicates: %d", nrow(df) - nrow(df_clean)))
    }

    # Handle different suffix modes
    if (suffix_mode == "average") {
        processed <- process_ls_pairs_average(
            all_mode_results = list(messy = df, clean = df_clean),
            mapping_output = mapping_output,
            debug = debug
        )
        return(processed)
    } else if (suffix_mode == "highest_effect") {
        processed <- process_ls_pairs_highest_effect(
            all_mode_results = list(messy = df, clean = df_clean),
            mapping_output = mapping_output,
            debug = debug
        )
        return(processed)
    } else if (suffix_mode %in% c("L", "S")) {
        if (debug) cat(sprintf("\nProcessing %s mode...\n", suffix_mode))

        # Re-map original symbols with suffixes
        df_with_suffixes <- df %>%
            left_join(
                mapping_output$original_symbols %>%
                    dplyr::select(gene, original_symbol),
                by = c("gene_symbol" = "gene")
            )

        # Filter for L or S genes using original symbols
        suffix_pattern <- paste0("\\.", suffix_mode, "$")
        df_messy <- df_with_suffixes %>%
            filter(grepl(suffix_pattern, original_symbol)) %>%
            dplyr::select(-original_symbol) # Remove the temporary column

        df_clean <- df_clean %>%
            left_join(
                mapping_output$original_symbols %>%
                    dplyr::select(gene, original_symbol),
                by = c("gene_symbol" = "gene")
            ) %>%
            filter(grepl(suffix_pattern, original_symbol)) %>%
            dplyr::select(-original_symbol)

        if (debug) {
            cat("\nProcessed Results Analysis:")
            cat(sprintf("\nTotal rows (messy): %d", nrow(df_messy)))
            cat(sprintf("\nTotal rows (clean): %d", nrow(df_clean)))
            cat(sprintf("\nUnique hugo symbols (messy): %d", n_distinct(df_messy$hugo_symbol)))
            cat(sprintf("\nUnique hugo symbols (clean): %d", n_distinct(df_clean$hugo_symbol)))
        }

        return(list(
            messy = df_messy %>% mutate(
                pvalue = sprintf("%.4f", as.numeric(pvalue)),
                padj = sprintf("%.4f", as.numeric(padj)),
                log2FoldChange = sprintf("%.4f", as.numeric(log2FoldChange))
            ),
            clean = df_clean %>% mutate(
                pvalue = sprintf("%.4f", as.numeric(pvalue)),
                padj = sprintf("%.4f", as.numeric(padj)),
                log2FoldChange = sprintf("%.4f", as.numeric(log2FoldChange))
            )
        ))
    }

    # Default return for "all" mode
    list(
        messy = df %>% mutate(
            pvalue = sprintf("%.4f", as.numeric(pvalue)),
            padj = sprintf("%.4f", as.numeric(padj)),
            log2FoldChange = sprintf("%.4f", as.numeric(log2FoldChange))
        ),
        clean = df_clean %>% mutate(
            pvalue = sprintf("%.4f", as.numeric(pvalue)),
            padj = sprintf("%.4f", as.numeric(padj)),
            log2FoldChange = sprintf("%.4f", as.numeric(log2FoldChange))
        )
    )
}

#' Combine p-values using Fisher's method
#' @param pvals Vector of p-values
#' @return Combined p-value
combine_pvalues_fisher <- function(pvals) {
    # Remove any NA values
    pvals <- pvals[!is.na(pvals)]

    # Handle edge cases
    if (length(pvals) == 0) {
        return(NA)
    }
    if (length(pvals) == 1) {
        return(pvals)
    }

    # Replace 0 with smallest possible number to avoid -Inf
    pvals[pvals == 0] <- .Machine$double.xmin

    # Fisher's method: -2 * sum(log(p-values)) follows chi-square with 2k df
    statistic <- -2 * sum(log(pvals))
    df <- 2 * length(pvals)

    # Calculate combined p-value
    combined_p <- pchisq(statistic, df, lower.tail = FALSE)

    return(combined_p)
}

#' Process L/S pairs for average mode
#' @param all_mode_results Results from all mode
#' @param mapping_output Mapping results
#' @param debug Debug flag
#' @return Processed results for average mode
process_ls_pairs_average <- function(all_mode_results, mapping_output, debug = TRUE) {
    # Get complete L/S pairs from original_symbols
    ls_pairs <- mapping_output$original_symbols %>%
        mutate(
            base_name = sub("\\.[LS]$", "", original_symbol)
        ) %>%
        group_by(base_name) %>%
        filter(
            n() == 2, # Only complete pairs
            all(grepl("\\.[LS]$", original_symbol))
        ) %>% # Both must have .L or .S
        ungroup()

    # Function to process messy/clean versions
    process_version <- function(df) {
        # Split the dataframe into pairs and non-pairs
        df_with_base <- df %>%
            mutate(
                base_gene = remove_suffix(gene_symbol)
            )

        pairs_data <- df_with_base %>%
            filter(base_gene %in% ls_pairs$base_name)

        non_pairs_data <- df_with_base %>%
            filter(!base_gene %in% ls_pairs$base_name)

        # Process pairs by averaging values and combining p-values
        processed_pairs <- pairs_data %>%
            group_by(base_gene, hugo_symbol) %>%
            summarise(
                gene_symbol = dplyr::first(base_gene),
                log2FoldChange = mean(log2FoldChange),
                pvalue = combine_pvalues_fisher(pvalue),
                padj = combine_pvalues_fisher(padj),
                .groups = "drop"
            )

        # Combine processed pairs with non-pairs
        bind_rows(
            processed_pairs %>% dplyr::select(-base_gene),
            non_pairs_data %>% dplyr::select(-base_gene)
        )
    }

    # Process both messy and clean versions
    df_messy <- process_version(all_mode_results$messy)
    df_clean <- process_version(all_mode_results$clean)
    
    # Count significant DEGs in average mode
    avg_messy_sig_count <- sum(df_messy$padj < 0.05, na.rm = TRUE)
    avg_clean_sig_count <- sum(df_clean$padj < 0.05, na.rm = TRUE)

    if (debug) {
        cat("\nProcessed Results Analysis (Average Mode):")
        cat(sprintf("\nTotal rows (messy): %d", nrow(df_messy)))
        cat(sprintf("\nTotal rows (clean): %d", nrow(df_clean)))
        cat(sprintf("\nUnique hugo symbols (messy): %d", n_distinct(df_messy$hugo_symbol)))
        cat(sprintf("\nUnique hugo symbols (clean): %d", n_distinct(df_clean$hugo_symbol)))
        cat(sprintf("\nNumber of .L/.S pairs processed: %d", n_distinct(ls_pairs$base_name)))
        cat(sprintf("\nSignificant DEGs (messy): %d", avg_messy_sig_count))
        cat(sprintf("\nSignificant DEGs (clean): %d", avg_clean_sig_count))
    }

    list(
        messy = df_messy,
        clean = df_clean
    )
}

#' Process L/S pairs for highest effect mode
#' @param all_mode_results Results from all mode
#' @param mapping_output Mapping results
#' @param debug Debug flag
#' @return Processed results for highest effect mode
process_ls_pairs_highest_effect <- function(all_mode_results, mapping_output, debug = TRUE) {
    # Get complete L/S pairs from original_symbols
    ls_pairs <- mapping_output$original_symbols %>%
        mutate(
            base_name = sub("\\.[LS]$", "", original_symbol)
        ) %>%
        group_by(base_name) %>%
        filter(
            n() == 2, # Only complete pairs
            all(grepl("\\.[LS]$", original_symbol))
        ) %>% # Both must have .L or .S
        ungroup()

    # Function to process messy/clean versions with new ranking logic
    process_version <- function(df) {
        # Split the dataframe into pairs and non-pairs
        df_with_base <- df %>%
            mutate(
                base_gene = remove_suffix(gene_symbol)
            )

        pairs_data <- df_with_base %>%
            filter(base_gene %in% ls_pairs$base_name)

        non_pairs_data <- df_with_base %>%
            filter(!base_gene %in% ls_pairs$base_name)

        # Process pairs using new ranking logic
        processed_pairs <- pairs_data %>%
            group_by(base_gene) %>%
            mutate(
                abs_log2FC = abs(log2FoldChange),
                # Create composite score based on fold change and -log10(padj)
                composite_score = abs_log2FC * (-log10(padj)),
                rank = dense_rank(dplyr::desc(composite_score)),
                is_L = grepl("\\.L$", gene_symbol),
                # If ranks are equal, prefer S over L
                selected = case_when(
                    n() == 1 ~ TRUE, # Keep single genes
                    rank == 1 & rank == lag(rank, default = 0) ~ !is_L, # If tied, choose S
                    rank == 1 ~ TRUE, # Choose highest rank
                    TRUE ~ FALSE
                )
            ) %>%
            ungroup() %>%
            filter(selected) %>%
            dplyr::select(-abs_log2FC, -composite_score, -rank, -is_L, -selected)

        # Combine processed pairs with non-pairs
        bind_rows(
            processed_pairs,
            non_pairs_data
        ) %>%
            dplyr::select(-base_gene)
    }

    # Process both messy and clean versions
    df_messy <- process_version(all_mode_results$messy)
    df_clean <- process_version(all_mode_results$clean)
    
    # Count significant DEGs in highest effect mode
    he_messy_sig_count <- sum(df_messy$padj < 0.05, na.rm = TRUE)
    he_clean_sig_count <- sum(df_clean$padj < 0.05, na.rm = TRUE)

    if (debug) {
        cat("\nProcessed Results Analysis (Highest Effect Mode):")
        cat(sprintf("\nTotal rows (messy): %d", nrow(df_messy)))
        cat(sprintf("\nTotal rows (clean): %d", nrow(df_clean)))
        cat(sprintf("\nUnique hugo symbols (messy): %d", n_distinct(df_messy$hugo_symbol)))
        cat(sprintf("\nUnique hugo symbols (clean): %d", n_distinct(df_clean$hugo_symbol)))
        cat(sprintf("\nNumber of .L/.S pairs processed: %d", n_distinct(ls_pairs$base_name)))
        cat(sprintf("\nSignificant DEGs (messy): %d", he_messy_sig_count))
        cat(sprintf("\nSignificant DEGs (clean): %d", he_clean_sig_count))
    }

    list(
        messy = df_messy,
        clean = df_clean
    )
}

#' Process data based on L or S suffix
#' @param data Data frame
#' @param suffix_type Suffix type (L or S)
#' @param original_symbols Original symbols lookup
#' @return Processed data
process_ls_data <- function(data, suffix_type, original_symbols) {
  # Validate suffix type
  if (!suffix_type %in% c("L", "S")) {
    stop("Invalid suffix_type. Must be either 'L' or 'S'")
  }
  
  # Create the suffix pattern
  suffix_pattern <- paste0("\\.", suffix_type, "$")
  
  # Join with original symbols to get the full symbols with suffixes
  data_with_suffixes <- data %>%
    left_join(
      original_symbols,
      by = c("gene_symbol" = "gene")
    )
  
  # Filter for genes with the specified suffix
  filtered_data <- data_with_suffixes %>%
    filter(grepl(suffix_pattern, original_symbol)) %>%
    # Keep only the necessary columns
    dplyr::select(-original_symbol)
  
  # Log processing info
  cat(sprintf("\nProcessed %s suffix data:", suffix_type))
  cat(sprintf("\n  Original rows: %d", nrow(data)))
  cat(sprintf("\n  Filtered rows: %d", nrow(filtered_data)))
  cat(sprintf("\n  Unique genes: %d", n_distinct(filtered_data$gene_symbol)))
  cat(sprintf("\n  Unique HUGO symbols: %d", n_distinct(filtered_data$hugo_symbol)))
  
  return(filtered_data)
}

message("Mapping functions loaded successfully. Use run_and_capture_mode() to perform gene mapping.")
