# Required libraries
library(dplyr) # For data manipulation
library(tidyr) # For data tidying
library(rentrez) # For NCBI Entrez queries
library(stringr) # For string manipulation
library(readr) # For reading files
library(org.Hs.eg.db) # For human gene mapping
library(AnnotationDbi) # For human gene mapping
library(officer) # For capturing console output
library(fs) # For capturing console output

# Create a function to run and capture output for a specific mode
run_and_capture_mode <- function(experiment_name, model_type, suffix_mode, output_file) {
    # Print start message
    message(paste0("\n\n========== STARTING ", suffix_mode, " MODE FOR ", experiment_name, " ==========\n"))
    
    # Create a temporary text file
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
    doc <- officer::read_docx()
    doc <- officer::body_add_par(doc, paste0("Analysis for ", suffix_mode, " mode"), style = "heading 1")
    
    # Add the captured text
    for (line in output_text) {
        doc <- officer::body_add_par(doc, line, style = "Normal")
    }
    
    # Save the document
    print(doc, target = output_file)
    
    # Clean up
    unlink(temp_file)
    
    # Get paths to the saved RDS files
    drughub_dir <- file.path(
        "results", experiment_name, "mapping",
        model_type,
        suffix_mode
    )
    mapping_rds_file <- file.path(drughub_dir, paste0(experiment_name, "_mapping_data.rds"))
    original_symbols_file <- file.path(drughub_dir, paste0(experiment_name, "_original_symbols.rds"))
    results_rds_file <- file.path(drughub_dir, paste0(experiment_name, "_", suffix_mode, "_results.rds"))
    
    # Print completion message
    message(paste0("\n✓ COMPLETED ", suffix_mode, " MODE FOR ", experiment_name))
    message(paste0("✓ Output saved to: ", output_file)) #Contains the captured console output from the entire analysis run
    message(paste0("✓ Mapping data saved to: ", mapping_rds_file)) #Contains the complete mapping output from map_xenopus_to_human()
    message(paste0("✓ Original symbols saved to: ", original_symbols_file)) #Contains the original symbols dataframe from the initial mapping step
    message(paste0("✓ Results data saved to: ", results_rds_file, "\n")) #Contains the complete results from the analysis run
    
    # Return both the result and the file paths
    return(list(
        result = result,
        output_file = output_file,
        mapping_file = mapping_rds_file,
        original_symbols_file = original_symbols_file,
        results_file = results_rds_file
    ))
}

# Remove suffix function
remove_suffix <- function(x) {
    sub("\\.[LS]$", "", x)
}

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

# Process raw counts based on mode
process_raw_counts <- function(counts_df, mapping_output, suffix_mode) {
    # First create basic messy version with human mapping
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

    # Process both versions based on mode
    counts_messy_processed <- switch(suffix_mode,
        "L" = process_ls_data(counts_messy, "L", mapping_output$original_symbols_lookup),
        "S" = process_ls_data(counts_messy, "S", mapping_output$original_symbols_lookup),
        "average" = {
            counts_messy %>%
                mutate(base_gene = remove_suffix(gene_symbol)) %>%
                group_by(base_gene, hugo_symbol) %>%
                summarise(
                    gene_symbol = paste(sort(unique(gene_symbol)), collapse = "/"),
                    across(where(is.numeric), mean, na.rm = TRUE),
                    .groups = "drop"
                )
        },
        "highest_effect" = {
            counts_messy %>%
                mutate(base_gene = remove_suffix(gene_symbol)) %>%
                group_by(base_gene, hugo_symbol) %>%
                slice_max(across(where(is.numeric), sum), n = 1) %>%
                ungroup() %>%
                dplyr::select(-base_gene)
        },
        "all" = counts_messy
    )

    counts_clean_processed <- switch(suffix_mode,
        "L" = process_ls_data(counts_clean, "L", mapping_output$original_symbols_lookup),
        "S" = process_ls_data(counts_clean, "S", mapping_output$original_symbols_lookup),
        "average" = {
            counts_clean %>%
                mutate(base_gene = remove_suffix(gene_symbol)) %>%
                group_by(base_gene, hugo_symbol) %>%
                summarise(
                    gene_symbol = paste(sort(unique(gene_symbol)), collapse = "/"),
                    across(where(is.numeric), mean, na.rm = TRUE),
                    .groups = "drop"
                )
        },
        "highest_effect" = {
            counts_clean %>%
                mutate(base_gene = remove_suffix(gene_symbol)) %>%
                group_by(base_gene, hugo_symbol) %>%
                slice_max(across(where(is.numeric), sum), n = 1) %>%
                ungroup() %>%
                dplyr::select(-base_gene)
        },
        "all" = counts_clean
    )

    # Add diagnostic counts
    if (debug) {
        cat("\nRaw Counts Processing:")
        cat("\nMessy version:")
        cat("\n  Total rows:", nrow(counts_messy_processed))
        cat("\n  Unique human genes:", n_distinct(counts_messy_processed$hugo_symbol))
        cat("\n  Averaged pairs:", if (suffix_mode == "average") sum(grepl("/", counts_messy_processed$gene_symbol)) else "N/A")

        cat("\n\nClean version:")
        cat("\n  Total rows:", nrow(counts_clean_processed))
        cat("\n  Unique human genes:", n_distinct(counts_clean_processed$hugo_symbol))
        cat("\n  Averaged pairs:", if (suffix_mode == "average") sum(grepl("/", counts_clean_processed$gene_symbol)) else "N/A")
    }

    list(
        messy = counts_messy_processed,
        clean = counts_clean_processed
    )
}

# Map Xenopus to human genes
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

    # 1. Create initial dataframes and lookups
    original_symbols_df <- data.frame(
        original_symbol = xenopus_genes,
        gene = remove_suffix(xenopus_genes)
    ) %>% distinct()
    
    # 2. Create unique gene list for mapping
    gene_list_no_suffix <- remove_suffix(xenopus_genes)
    Xenopus_gene_list_unique <- data.frame(gene = unique(gene_list_no_suffix))

    # 3. Load and clean reference data
    xenbase_data <- read.delim(xenbase_file,
        header = TRUE, sep = "\t",
        quote = "", fill = TRUE
    )

    hcop_data <- read.delim(hcop_file,
        header = TRUE, sep = "\t",
        quote = "", fill = TRUE
    )

    # 4. Standardize gene names (lowercase and trim)
    Xenopus_gene_list_unique$gene <- tolower(trimws(Xenopus_gene_list_unique$gene))
    xenbase_data$Xenbase_gene_symbol <- tolower(trimws(xenbase_data$Xenbase_gene_symbol))
    hcop_data$xenopus_symbol <- tolower(trimws(hcop_data$xenopus_symbol))

    xenbase_data <- xenbase_data %>%
        dplyr::select(Xenbase_gene_symbol, Human_NCBI_Entrez_ID)

    hcop_data <- hcop_data %>%
        dplyr::select(xenopus_symbol, human_entrez_gene, human_ensembl_gene)

    # 5. Initial mappings and duplicate analysis
    # First get raw matches and analyze duplicates in each database
    # Xenbase matching sequence
    xenbase_initial <- Xenopus_gene_list_unique %>%
        left_join(xenbase_data, by = c("gene" = "Xenbase_gene_symbol")) %>%
        filter(!is.na(Human_NCBI_Entrez_ID)) %>%
        mutate(
            human_entrez_id = as.character(Human_NCBI_Entrez_ID),
            source = "xenbase"
        ) %>%
        dplyr::select(-Human_NCBI_Entrez_ID) %>%
        distinct()

    write.csv(xenbase_initial, "xenbase_initial.csv", row.names = FALSE)

    # HCOP matching sequence
    hcop_initial <- Xenopus_gene_list_unique %>%
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

    write.csv(hcop_initial, "hcop_initial.csv", row.names = FALSE)

    # Create ID mappings
    # Analyze initial mapping duplicates - accounting for case 2 duplicates
    xenbase_id_mappings <- xenbase_initial %>%
        group_by(human_entrez_id) %>%
        summarise(
            xenbase_genes = paste(unique(gene), collapse = ";"),
            n_xenopus_genes = n_distinct(gene)
        )

    # write.csv(xenbase_id_mappings, "xenbase_id_mappings.csv", row.names = FALSE)

    hcop_id_mappings <- hcop_initial %>%
        group_by(human_entrez_id, human_ensembl_id) %>%
        summarise(
            hcop_genes = paste(unique(gene), collapse = ";"),
            n_xenopus_genes = n_distinct(gene)
        )

    # write.csv(hcop_id_mappings, "hcop_id_mappings.csv", row.names = FALSE)


    if (debug) {
        # Count genes with different suffix patterns
        suffix_patterns <- data.frame(
            original_symbol = original_symbols_df$original_symbol, # Use original_symbol column
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
            original_symbol = original_symbols_df$original_symbol, # Use original_symbol column
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
    }

    # 6. Debug output
    if (debug) {
        cat("\nOriginal Symbols Diagnostics:\n")

        cat("\nGene Counts:")
        cat("\n  All Initial Genes:", length(xenopus_genes))
        cat("\n  All Initial Unique Genes:", nrow(Xenopus_gene_list_unique))

        cat("\nSuffix counts:")
        cat("\n.L genes:", suffix_patterns$L_count)
        cat("\n.S genes:", suffix_patterns$S_count)
        cat("\nNeither .L nor .S:", suffix_patterns$neither_count)
        cat("\nComplete L/S pairs:", pair_analysis$total_pairs)
        cat("\nL only (no S partner):", pair_analysis$L_only)
        cat("\nS only (no L partner):", pair_analysis$S_only)

        cat("\nID Mapping Analysis:")
        cat("\nXenbase:")
        cat("\n  Total Unique Entrez IDs:", n_distinct(na.omit(xenbase_initial$human_entrez_id)))
        cat("\n  Total IDs:", length(na.omit(xenbase_initial$human_entrez_id)))
        cat("\n  Entrez IDs mapping to multiple Xenopus genes (Case 2 duplicates):", sum(xenbase_id_mappings$n_xenopus_genes > 1))
        cat("\n  Max Xenopus genes per Entrez ID:", max(xenbase_id_mappings$n_xenopus_genes))

        cat("\n\nHCOP:")
        cat("\n  Total Unique Entrez IDs:", n_distinct(na.omit(hcop_initial$human_entrez_id)))
        cat("\n  Total Entrez IDs:", length(na.omit(hcop_initial$human_entrez_id)))
        cat("\n  Total Unique Ensembl IDs:", n_distinct(na.omit(hcop_initial$human_ensembl_id)))
        cat("\n  Total Ensembl IDs:", length(na.omit(hcop_initial$human_ensembl_id)))
        cat("\n  IDs mapping to multiple Xenopus genes (Case 2 duplicates):", sum(hcop_id_mappings$n_xenopus_genes > 1))
        cat("\n  Max Xenopus genes per ID:", max(hcop_id_mappings$n_xenopus_genes))
    }

    # 7. Get HUGO symbols for all potential human IDs
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
                human_entrez_id = character(0),
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
                human_ensembl_id = character(0),
            )
        }
    )

    xenbase_hcop_merge <- bind_rows(
        xenbase_initial %>% mutate(source = "xenbase"),
        hcop_initial %>% mutate(source = "hcop")
    ) %>%
        distinct()

    entrez_xenbase_hcop_hugo_symbol_merge <- xenbase_hcop_merge %>%
        left_join(hugo_from_entrez,
            by = "human_entrez_id",
            relationship = "many-to-many"
        ) %>%
        distinct()

    entrez_ensemble_xenbase_hcop_hugo_symbol_merge <- entrez_xenbase_hcop_hugo_symbol_merge %>%
        left_join(hugo_from_ensembl,
            by = "human_ensembl_id",
            relationship = "many-to-many"
        ) %>%
        distinct()

    # We are just going to use the hugo_symbol_entrez column as the HUGO symbol
    # because it has more mappings and better ones in the cases where there are conflicts
    gene_entrez_hugo_df <- entrez_ensemble_xenbase_hcop_hugo_symbol_merge %>%
        # create hugo_symbol column with values from hugo_symbol_entrez or hugo_symbol_ensembl if entrez is NA
        mutate(hugo_symbol = ifelse(is.na(hugo_symbol_entrez), hugo_symbol_ensembl, hugo_symbol_entrez)) %>%
        dplyr::select(gene, human_entrez_id, source, hugo_symbol) %>%
        distinct()

    # Create lists of genes for each source
    # Get xenbase list
    xenbase_list <- gene_entrez_hugo_df %>%
        dplyr::filter(source == "xenbase") %>%
        dplyr::select(gene) %>%
        distinct()
    # Convert to a list
    xenbase_list <- xenbase_list$gene

    # Get hcop list
    hcop_list <- gene_entrez_hugo_df %>%
        dplyr::filter(source == "hcop") %>%
        dplyr::select(gene) %>%
        distinct()
    # Convert to a list
    hcop_list <- hcop_list$gene

    # For each row in gene_entrez_hugo_df, determine if the gene is in xenbase_list, hcop_list, or both
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

    # case 1.1 is when for each gene, there are multiple rows that all have the same hugo symbol
    # case 1.2 is when for each gene, there are multiple rows which do not all have the same hugo symbol
    # lets identify each case with new boolean columns first
    # 9. Handle Case 1 duplicates (after getting HUGO symbols but before final processing)
    gene_entrez_hugo_df_dupes <- gene_entrez_hugo_df %>%
        group_by(gene) %>%
        mutate(
            case_1_1 = n() > 1 & n_distinct(hugo_symbol) == 1,
            case_1_2 = n() > 1 & n_distinct(hugo_symbol) > 1
        ) %>%
        ungroup()
    # Report out duplicate info
    cat(sprintf("Total rows: %d", nrow(gene_entrez_hugo_df)))
    cat(sprintf("\nTotal duplicate trows: %d", sum(gene_entrez_hugo_df_dupes$case_1_1 | gene_entrez_hugo_df_dupes$case_1_2)))
    cat(sprintf("\nCase 1.1 (multiple rows, same ID & HUGO): %d rows", sum(gene_entrez_hugo_df_dupes$case_1_1)))
    cat(sprintf("\nCase 1.2 (different IDs): %d rows", sum(gene_entrez_hugo_df_dupes$case_1_2)))

    write.csv(gene_entrez_hugo_df_dupes, "gene_entrez_hugo_df_dupes.csv", row.names = FALSE)

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
    write.csv(merged_data_clean, "merged_data_clean.csv", row.names = FALSE)

    if (debug) {
        cat("\nOverlap Analysis:")
        cat("\nTotal overlap between databases:", nrow(merged_data_clean %>% filter(source == "both")), "genes")
        cat("\n\nUnique matches:")
        cat("\nXenbase only:", nrow(merged_data_clean %>% filter(source == "xenbase")), "genes")
        cat("\nHCOP only:", nrow(merged_data_clean %>% filter(source == "hcop")), "genes")
    }

    # 12. Process based on suffix mode
    processed_matches <- list(
        data = merged_data_clean,
        de_results = de_results,
        experiment_name = experiment_name,
        model_type = model_type
    )
    # Print gene counts breakdown first
    gene_counts <- merged_data_clean %>%
        mutate(
            has_L = grepl("\\.L$", gene_symbol),
            has_S = grepl("\\.S$", gene_symbol),
            base_gene = sub("\\.[LS]$", "", gene_symbol)
        ) %>%
        summarise(
            total = n(),
            L_count = sum(has_L),
            S_count = sum(has_S),
            neither_count = sum(!has_L & !has_S),
            unique_base_genes = n_distinct(base_gene)
        )

    # After we've created merged_data_clean, add duplicate type information
    
    # First, identify the different duplicate types
    gene_entrez_hugo_df_dupes <- gene_entrez_hugo_df %>%
        group_by(gene) %>%
        mutate(
            case_1_1 = n() > 1 & n_distinct(hugo_symbol) == 1,
            case_1_2 = n() > 1 & n_distinct(hugo_symbol) > 1
        ) %>%
        ungroup()
        
    # Then identify case 2 duplicates (multiple Xenopus genes to same human gene)
    case_2_duplicates <- merged_data_clean %>%
        group_by(hugo_symbol) %>%
        filter(n() > 1) %>%
        pull(hugo_symbol) %>%
        unique()
    
    # Now enhance original_symbols_df with duplicate type information
    original_symbols_df <- original_symbols_df %>%
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
    
    # Now write the enhanced version with hugo_symbol and duplicate info included
    write.csv(original_symbols_df, "original_symbols_df.csv", row.names = FALSE)
    
    # Return results
    return(list(
        data = processed_matches$data,
        original_symbols = original_symbols_df
    ))
}

format_for_drughub <- function(experiment_name,
                               model_type = "deseq2",
                               samples_to_keep = NULL,
                               suffix_mode = "all",
                               de_results = NULL,
                               debug = TRUE) {
    # Load experiment data
    experiment_path <- file.path("data", experiment_name, paste0(experiment_name, ".rds"))
    if (!file.exists(experiment_path)) {
        stop("Experiment data not found at: ", experiment_path)
    }
    experiment_obj <- readRDS(experiment_path)

    # Load DE results if not provided
    if (is.null(de_results)) {
        de_results_path <- file.path(
            "results", experiment_name, "de", model_type, "simple",
            paste0(experiment_name, "_de_results.rds")
        )
        if (!file.exists(de_results_path)) {
            stop("DE results not found at: ", de_results_path)
        }
        de_results <- readRDS(de_results_path)
    }

    # Process raw counts
    assay_data <- experiment_obj$assay_data[[paste0(experiment_name, "_assay_data")]]
    raw_counts_df <- as.data.frame(assay_data)

    # Add gene_symbol column if it's in rownames
    if (!"gene_symbol" %in% colnames(raw_counts_df)) {
        raw_counts_df$gene_symbol <- rownames(raw_counts_df)
    }

    # Run mapping with raw data
    mapping_output <- map_xenopus_to_human(
        xenopus_genes = raw_counts_df$gene_symbol,
        xenbase_file = file.path("pipeline_input_data_files", "Xenbase_Xenopus_Orthology_Predictions.txt"),
        hcop_file = file.path("pipeline_input_data_files", "HCOP_Xenopus_Orthology_Predictions.txt"),
        suffix_mode = "all",
        debug = debug
    )

    processed_results <- switch(suffix_mode,
        "all" = mapping_output,
        "average" = mapping_output, # Just pass through mapping data
        "highest_effect" = mapping_output, # Just pass through mapping data
        "L" = list(data = mapping_output$data %>% filter(grepl("\\.L$", gene_symbol))),
        "S" = list(data = mapping_output$data %>% filter(grepl("\\.S$", gene_symbol)))
    )

    # Analyze L/S divergence after we have both mapping and DE results
    ls_divergence <- NULL # Initialize to NULL
    if (!is.null(de_results) && suffix_mode == "all") { # Only run for "all" mode
        # Calculate basic L/S pair statistics first
        gene_pair_stats <- mapping_output$data %>%
            mutate(
                has_L = grepl("\\.L$", gene),
                has_S = grepl("\\.S$", gene),
                base_gene = remove_suffix(gene)
            ) %>%
            group_by(base_gene) %>%
            summarise(
                has_both = any(has_L) & any(has_S),
                .groups = "drop"
            )

        ls_divergence <- analyze_ls_divergence(
            data = mapping_output$data,
            de_results = de_results,
            debug = TRUE
        )
        # Print summary statistics
        cat("\nL/S Divergence Analysis:")
        cat("\nTotal genes with both L/S:", sum(gene_pair_stats$has_both))

        for (contrast in names(ls_divergence)) {
            stats <- ls_divergence[[contrast]]$stats
            total_pairs <- stats$total_pairs

            cat(sprintf("\n\nContrast: %s", contrast))
            cat(sprintf("\nTotal L/S pairs: %d", total_pairs))
            if (total_pairs > 0) {
                cat(sprintf(
                    "\nOpposite directions: %d (%.1f%%)",
                    stats$direction_conflicts,
                    100 * stats$direction_conflicts / total_pairs
                ))
                cat(sprintf(
                    "\nOpposite directions (both significant): %d (%.1f%%)",
                    stats$sig_direction_conflicts,
                    100 * stats$sig_direction_conflicts / total_pairs
                ))
                cat(sprintf(
                    "\nLarge FC difference (>1): %d (%.1f%%)",
                    stats$large_fc_diff,
                    100 * stats$large_fc_diff / total_pairs
                ))
                cat(sprintf(
                    "\nLarge padj difference (>0.1): %d (%.1f%%)",
                    stats$large_padj_diff,
                    100 * stats$large_padj_diff / total_pairs
                ))
                cat(sprintf(
                    "\nBoth large FC and padj differences: %d (%.1f%%)",
                    stats$both_large_diff,
                    100 * stats$both_large_diff / total_pairs
                ))
            }
        }

        # Save conflicts to file if they exist
        if (length(ls_divergence) > 0) {
            all_conflicts <- bind_rows(lapply(names(ls_divergence), function(contrast) {
                ls_divergence[[contrast]]$pairs %>%
                    filter(direction_conflict | fc_diff > 1 | padj_diff > 0.1) %>%
                    mutate(contrast = contrast)
            }))

            if (nrow(all_conflicts) > 0) {
                conflict_file <- file.path(
                    "results", experiment_name, "mapping", model_type,
                    suffix_mode, # Changed from nested "all" folder
                    paste0(experiment_name, "_LS_conflicts.csv")
                )
                # Ensure directory exists before writing
                dir.create(dirname(conflict_file), recursive = TRUE, showWarnings = FALSE)
                write.csv(all_conflicts, conflict_file, row.names = FALSE)
                message("Conflicts written to: ", conflict_file)
            }
        }
    }

    # Check if mapping was successful
    if (nrow(mapping_output$data) == 0) {
        stop("No genes were successfully mapped to human orthologs")
    }

    # If samples_to_keep is NULL, use all samples
    if (is.null(samples_to_keep)) {
        samples_to_keep <- setdiff(colnames(raw_counts_df), "gene_symbol")
    }

    # Process raw counts with mapping
    raw_counts_mapped <- raw_counts_df %>%
        left_join(
            mapping_output$data,
            by = c("gene_symbol" = "gene")
        ) %>%
        filter(!is.na(hugo_symbol))

    # Create messy and clean versions of raw counts
    raw_counts_messy <- raw_counts_mapped
    raw_counts_clean <- raw_counts_mapped %>%
        group_by(hugo_symbol) %>%
        filter(n_distinct(gene_symbol) == 1) %>%
        ungroup()

    # Process raw counts - only for "all" mode
    raw_counts_results <- if (suffix_mode == "all") {
        list(
            messy = raw_counts_messy,
            clean = raw_counts_clean
        )
    } else {
        NULL # No raw counts processing for other modes
    }

    # Process DEA results for each contrast
    if (debug) cat("\nProcessing DEA results with detailed diagnostics...\n")

    dea_results_list <- lapply(names(de_results), function(contrast_name) {
        result <- de_results[[contrast_name]]
        contrast_info <- result$comparison_info
        contrast_label <- paste(contrast_info$variable,
            contrast_info$experimental,
            "vs",
            contrast_info$control,
            sep = "_"
        )

        if (debug) cat(sprintf("\n\nProcessing contrast: %s\n", contrast_label))

        # Process DEA results
        processed_results <- process_de_results(result, model_type, mapping_output, suffix_mode, debug)

        list(
            name = contrast_label,
            data = processed_results
        )
    })

    # Create output directory
    drughub_dir <- file.path(
        "results", experiment_name, "mapping",
        model_type,
        suffix_mode
    )

    dir.create(drughub_dir, recursive = TRUE, showWarnings = FALSE)

    # Save both versions for each contrast
    for (contrast in dea_results_list) {
        # Format p-values for messy version
        messy_formatted <- contrast$data$messy %>%
            mutate(
                pvalue = sprintf("%.4f", as.numeric(pvalue)),
                padj = sprintf("%.4f", as.numeric(padj)),
                log2FoldChange = sprintf("%.4f", as.numeric(log2FoldChange))
            )

        # Format p-values for clean version
        clean_formatted <- contrast$data$clean %>%
            mutate(
                pvalue = sprintf("%.4f", as.numeric(pvalue)),
                padj = sprintf("%.4f", as.numeric(padj)),
                log2FoldChange = sprintf("%.4f", as.numeric(log2FoldChange))
            )


        # Only do the final column selection right before writing to file
        messy_final <- messy_formatted %>%
            dplyr::select(
                hugo_symbol,
                log2_fc = log2FoldChange,
                p_value = pvalue,
                adj_p_value = padj
            )

        clean_final <- clean_formatted %>%
            dplyr::select(
                hugo_symbol,
                log2_fc = log2FoldChange,
                p_value = pvalue,
                adj_p_value = padj
            )

        # Save messy version
        write.csv(
            messy_final,
            file.path(drughub_dir, paste0(
                experiment_name, "_",
                contrast$name, "_DEGs_messy.csv" # Added "_messy" suffix
            )),
            row.names = FALSE,
            quote = TRUE
        )

        # Save clean version
        write.csv(
            clean_final,
            file.path(drughub_dir, paste0(
                experiment_name, "_",
                contrast$name, "_DEGs_clean.csv"
            )),
            row.names = FALSE,
            quote = TRUE
        )
    }

    # Save raw counts if in "all" mode
    if (suffix_mode == "all") {
        write.csv(
            raw_counts_mapped, # Use the mapped counts directly
            file.path(drughub_dir, paste0(experiment_name, "_raw_counts.csv")), # Removed suffix
            row.names = FALSE,
            quote = TRUE
        )
    }
    
    # Save mapping data to RDS files
    mapping_rds_file <- file.path(drughub_dir, paste0(experiment_name, "_mapping_data.rds"))
    saveRDS(mapping_output, mapping_rds_file)
    message("Mapping data saved to: ", mapping_rds_file)
    
    # Save original symbols dataframe separately
    original_symbols_file <- file.path(drughub_dir, paste0(experiment_name, "_original_symbols.rds"))
    saveRDS(mapping_output$original_symbols, original_symbols_file)
    message("Original symbols data saved to: ", original_symbols_file)
    
    # Save processed results to RDS files
    results_rds_file <- file.path(drughub_dir, paste0(experiment_name, "_", suffix_mode, "_results.rds"))
    results_to_save <- list(
        raw_counts = raw_counts_results,
        dea_results = dea_results_list,
        mapping_output = mapping_output,
        original_symbols = mapping_output$original_symbols,
        divergence_analysis = ls_divergence,
        suffix_mode = suffix_mode
    )
    saveRDS(results_to_save, results_rds_file)
    message("Processed results saved to: ", results_rds_file)

    # Return results structure
    list(
        raw_counts = raw_counts_results,
        dea_results = dea_results_list,
        mapping_output = mapping_output,
        original_symbols = mapping_output$original_symbols,
        divergence_analysis = ls_divergence
    )
}

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
        
        # Count duplicate xenopus genes with distinct values
        duplicates_messy <- df %>%
            group_by(gene_symbol) %>%
            summarise(
                distinct_values = n_distinct(paste(log2FoldChange, pvalue, padj)),
                .groups = "drop"
            ) %>%
            filter(distinct_values > 1)
        
        cat(sprintf("\nDuplicate xenopus genes with distinct values: %d", nrow(duplicates_messy)))
        cat(sprintf("\nUnique hugo symbols: %d", n_distinct(df$hugo_symbol)))
        
        # Analyze clean version
        cat("\n\nCLEAN VERSION:")
        cat(sprintf("\nTotal rows: %d", nrow(df_clean)))
        
        # Count duplicate xenopus genes with distinct values in clean version
        duplicates_clean <- df_clean %>%
            group_by(gene_symbol) %>%
            summarise(
                distinct_values = n_distinct(paste(log2FoldChange, pvalue, padj)),
                .groups = "drop"
            ) %>%
            filter(distinct_values > 1)
        
        cat(sprintf("\nDuplicate xenopus genes with distinct values: %d", nrow(duplicates_clean)))
        cat(sprintf("\nUnique hugo symbols: %d", n_distinct(df_clean$hugo_symbol)))
        
        # Print Case 2 duplicate stats
        cat("\n\nCase 2 Duplicate Analysis:")
        cat(sprintf("\nHUGO symbols appearing multiple times (removed from clean): %d", length(duplicate_hugos)))
        cat(sprintf("\nTotal rows removed due to Case 2 duplicates: %d", nrow(df) - nrow(df_clean)))

        # Analyze Case 2 duplicates
        case2_analysis <- df %>%
            group_by(hugo_symbol) %>%
            summarise(
                xenopus_genes = list(sort(unique(gene_symbol))),
                n_xenopus = n_distinct(gene_symbol),
                # Remove .L/.S suffixes and compare base names
                base_names = list(sort(unique(sub("\\.[LS]$", "", gene_symbol)))),
                n_unique_base_names = n_distinct(sub("\\.[LS]$", "", gene_symbol)),
                # Calculate if genes are likely from same family
                is_same_family = length(unique(sub("[0-9]+[LS]?$", "", gene_symbol))) == 1,
                .groups = "drop"
            ) %>%
            filter(n_xenopus > 1) # Only look at multi-mapping cases

        # Calculate stats
        multi_map_groups <- nrow(case2_analysis)
        total_rows_to_drop <- sum(case2_analysis$n_xenopus) - multi_map_groups
        avg_xenopus_per_group <- mean(case2_analysis$n_xenopus)

        # Analyze gene family patterns
        same_family_count <- sum(case2_analysis$is_same_family)
        ls_pairs_count <- sum(case2_analysis$n_xenopus > case2_analysis$n_unique_base_names)

        cat("\nCleaning Impact Analysis:")
        cat(sprintf("\nHuman genes with multiple Xenopus mappings: %d", multi_map_groups))
        cat(sprintf("\nTotal rows to be dropped in clean mode: %d", total_rows_to_drop))
        cat(sprintf("\nAverage Xenopus genes per human gene (in multi-map groups): %.2f", avg_xenopus_per_group))

        cat("\n\nGene Family Analysis:")
        cat(sprintf(
            "\nMulti-mappings from same gene family: %d (%.1f%%)",
            same_family_count, 100 * same_family_count / multi_map_groups
        ))
        cat(sprintf(
            "\nMulti-mappings that are L/S pairs: %d (%.1f%%)",
            ls_pairs_count, 100 * ls_pairs_count / multi_map_groups
        ))
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

# Combine p-values using Fisher's method
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

# Process L/S pairs for average mode
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

# Process L/S pairs for highest effect mode
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

# Process data based on L or S suffix
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

# Example usage with different modes:
xen_tran_2024_03_deseq2_drughub_data_all <- run_and_capture_mode(
    "xen_tran_2024_03", "deseq2", "all",
    "results/xen_tran_2024_03/mapping/deseq2/all/xen_tran_2024_03_all_mode_analysis.docx"
)

# Load original symbols dataframe
original_symbols_df <- readRDS("results/xen_tran_2024_03/mapping/deseq2/all/xen_tran_2024_03_original_symbols.rds")

# View first few rows and summary
print("First few rows of original_symbols_df:")
print(head(original_symbols_df))
print("\nSummary of original_symbols_df:")
print(summary(original_symbols_df))

# Run the mode analysis
xen_tran_2024_03_deseq2_drughub_data_all <- run_and_capture_mode(
    "xen_tran_2024_03", "deseq2", "all",
    "results/xen_tran_2024_03/mapping/deseq2/all/xen_tran_2024_03_all_mode_analysis.docx"
)

xen_tran_2024_03_deseq2_drughub_data_avg <- run_and_capture_mode(
    "xen_tran_2024_03", "deseq2", "average",
    "results/xen_tran_2024_03/mapping/deseq2/average/xen_tran_2024_03_avg_mode_analysis.docx"
)

xen_tran_2024_03_deseq2_drughub_data_highest_effect <- run_and_capture_mode(
    "xen_tran_2024_03", "deseq2", "highest_effect",
    "results/xen_tran_2024_03/mapping/deseq2/highest_effect/xen_tran_2024_03_highest_effect_mode_analysis.docx"
)

xen_tran_2024_03_deseq2_drughub_data_L <- run_and_capture_mode(
    "xen_tran_2024_03", "deseq2", "L",
    "results/xen_tran_2024_03/mapping/deseq2/L/xen_tran_2024_03_L_mode_analysis.docx"
)

xen_tran_2024_03_deseq2_drughub_data_S <- run_and_capture_mode(
    "xen_tran_2024_03", "deseq2", "S",
    "results/xen_tran_2024_03/mapping/deseq2/S/xen_tran_2024_03_S_mode_analysis.docx"
)




# Example usage with different modes:
xen_tran_2024_12_deseq2_drughub_data_all <- run_and_capture_mode(
    "xen_tran_2024_12", "deseq2", "all",
    "results/xen_tran_2024_12/mapping/deseq2/all/xen_tran_2024_12_all_mode_analysis.docx"
)

xen_tran_2024_12_deseq2_drughub_data_avg <- run_and_capture_mode(
    "xen_tran_2024_12", "deseq2", "average",
    "results/xen_tran_2024_12/mapping/deseq2/average/xen_tran_2024_12_avg_mode_analysis.docx"
)

xen_tran_2024_12_deseq2_drughub_data_highest_effect <- run_and_capture_mode(
    "xen_tran_2024_12", "deseq2", "highest_effect",
    "results/xen_tran_2024_12/mapping/deseq2/highest_effect/xen_tran_2024_12_highest_effect_mode_analysis.docx"
)

xen_tran_2024_12_deseq2_drughub_data_L <- run_and_capture_mode(
    "xen_tran_2024_12", "deseq2", "L",
    "results/xen_tran_2024_12/mapping/deseq2/L/xen_tran_2024_12_L_mode_analysis.docx"
)

xen_tran_2024_12_deseq2_drughub_data_S <- run_and_capture_mode(
    "xen_tran_2024_12", "deseq2", "S",
    "results/xen_tran_2024_12/mapping/deseq2/S/xen_tran_2024_12_S_mode_analysis.docx"
)
