#####
# DESeq2 Results Export Functions

# Required libraries
library(DESeq2) # For accessing DESeq2 results
library(dplyr) # For data manipulation
library(magrittr) # For pipe operations
library(tibble) # For tibble operations
library(flextable) # For creating formatted tables
library(officer) # For saving Word documents
library(ggplot2) # For volcano plots
library(ggrepel) # For label repulsion in volcano plots
library(ComplexHeatmap) # For creating heatmaps
library(circlize) # For heatmap color functions
library(grid) # For heatmap graphics
library(RColorBrewer) # For color palettes

# Main export function with memory optimization
export_all_result_types <- function(experiment_name,
                                    de_params = list(min_padj = 0.05, min_log2fc = 0),
                                    table_params = list(direction = "both", n_genes = 25),
                                    volcano_params = list(n_labels = 15),
                                    heatmap_params = list(
                                      direction = "both", n_genes = 25,
                                      scale = TRUE, n_clusters = 3
                                    ),
                                    analysis_type = "simple") {
  message("Starting export of all result types...")

  # Load experiment data (existing code)
  experiment_obj <- readRDS(file.path(
    "data", experiment_name,
    paste0(experiment_name, ".rds")
  ))

  # Load DESeq2 results with correct filename pattern
  results_filename <- if (analysis_type == "complex") {
    paste0(experiment_name, "_de_complex_results.rds")
  } else {
    paste0(experiment_name, "_de_results.rds")
  }

  dea_results <- readRDS(file.path(
    "results", experiment_name, "de", "deseq2",
    analysis_type, results_filename
  ))

  # Base directories setup
  base_dir <- file.path("results", experiment_obj$metadata$config$name, "de", "deseq2", analysis_type)
  dirs <- c("csv", "tables", "plots/volcano", "plots/heatmap")
  for (d in dirs) {
    dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
  }

  # Process results based on analysis type
  if (analysis_type == "complex") {
    # For complex analysis, iterate through each contrast
    message("Processing complex analysis results...")

    # Extract design formula from dds object
    design <- as.character(dea_results$dds@design)[2]
    # Get the main variable (assuming it's the first term in the formula)
    main_variable <- strsplit(design, " *\\+ *")[[1]][1]
    main_variable <- gsub("^~", "", main_variable)

    # Process each contrast
    for (contrast_name in names(dea_results$deseq_results)) {
      message("\nProcessing contrast: ", contrast_name)

      # Create a result object with all necessary components
      result <- list(
        deseq_result = dea_results$deseq_results[[contrast_name]],
        dds = dea_results$dds,
        comparison_info = list(
          variable = main_variable,
          experimental = strsplit(contrast_name, "_vs_")[[1]][1],
          control = strsplit(contrast_name, "_vs_")[[1]][2],
          contrasts = dea_results$comparison_info$contrasts # Add this line
        )
      )

      # Use process_contrast helper function
      process_contrast(
        result, experiment_obj, contrast_name,
        de_params, table_params, volcano_params, heatmap_params,
        analysis_type
      )
    }
  } else {
    # Simple analysis remains unchanged
    for (contrast_name in names(dea_results)) {
      message("\nProcessing contrast: ", contrast_name)

      # Extract design formula from dds object
      design <- as.character(dea_results[[contrast_name]]$dds@design)[2]
      # Get the main variable (assuming it's the first term in the formula)
      main_variable <- strsplit(design, " *\\+ *")[[1]][1]
      main_variable <- gsub("^~", "", main_variable)

      # Add comparison info to the result
      dea_results[[contrast_name]]$comparison_info <- list(
        variable = main_variable,
        experimental = strsplit(contrast_name, "_vs_")[[1]][1],
        control = strsplit(contrast_name, "_vs_")[[1]][2]
      )

      # Use process_contrast helper function
      process_contrast(
        dea_results[[contrast_name]], experiment_obj, contrast_name,
        de_params, table_params, volcano_params, heatmap_params,
        analysis_type
      )
    }
  }

  message("\nExport completed!")
  return(experiment_obj)
}

# New helper function to process individual contrasts
process_contrast <- function(result, experiment_obj, contrast_name,
                             de_params, table_params, volcano_params, heatmap_params,
                             analysis_type) {
  # 1. CSV Export
  message("Exporting CSV results...")
  tryCatch(
    {
      for (dir in c("both", "up", "down")) {
        export_de_results(result, experiment_obj,
          direction = dir,
          min_padj = de_params$min_padj, min_log2fc = de_params$min_log2fc,
          contrast_name = contrast_name, analysis_type = analysis_type
        )
      }
    },
    error = function(e) message("Error in CSV export: ", e$message)
  )
  gc()

  # 2. Table Export
  message("Generating tables...")
  tryCatch(
    {
      generate_de_tables(result, experiment_obj,
        direction = table_params$direction,
        n_genes = table_params$n_genes,
        min_padj = de_params$min_padj,
        min_log2fc = de_params$min_log2fc,
        contrast_name = contrast_name,
        analysis_type = analysis_type
      )
    },
    error = function(e) message("Error in table generation: ", e$message)
  )
  gc()

  # 3. Volcano Plot
  message("Generating volcano plot...")
  tryCatch(
    {
      generate_volcano_plot(result, experiment_obj,
        n_labels = volcano_params$n_labels,
        contrast_name = contrast_name,
        analysis_type = analysis_type
      )
    },
    error = function(e) message("Error in volcano plot generation: ", e$message)
  )
  gc()

  # 4. Heatmap
  message("Generating heatmap...")
  tryCatch(
    {
      generate_complex_heatmap(
        result = result,
        experiment_obj = experiment_obj,
        direction = heatmap_params$direction,
        n_genes = heatmap_params$n_genes,
        scale = heatmap_params$scale,
        n_clusters = heatmap_params$n_clusters,
        min_padj = de_params$min_padj,
        min_log2fc = de_params$min_log2fc,
        contrast_name = contrast_name,
        analysis_type = analysis_type
      )
    },
    error = function(e) message("Error in heatmap generation: ", e$message)
  )
  gc()
}

# Modified export functions to use updated directory structure
export_de_results <- function(result, experiment_obj, direction = "both",
                              min_padj = 0.05, min_log2fc = 0, contrast_name,
                              analysis_type = "simple") {
  # Create DE results directory with updated path
  de_dir <- file.path("results", experiment_obj$metadata$config$name, "de", "deseq2", analysis_type)
  csv_dir <- file.path(de_dir, "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

  # Convert DESeqResults to data frame
  res_df <- as.data.frame(result$deseq_result) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::filter(!is.na(padj), padj < min_padj, abs(log2FoldChange) >= min_log2fc)

  if (direction == "up") {
    res_df <- dplyr::filter(res_df, log2FoldChange > 0)
  } else if (direction == "down") {
    res_df <- dplyr::filter(res_df, log2FoldChange < 0)
  }

  filename <- file.path(
    csv_dir,
    paste0(
      contrast_name, "_", direction,
      "_padj", min_padj, "_log2fc", min_log2fc, ".csv"
    )
  )

  write.csv(res_df, filename, row.names = FALSE)
}

# Keep the original create_de_table exactly as is for simple analysis
create_de_table <- function(result, direction = "both", n_genes = NULL,
                            min_padj = 0.05, min_log2fc = 0) {
  # Extract the DESeq result and dds object
  deseq_result <- result$deseq_result
  dds <- result$dds
  comparison_info <- result$comparison_info

  # Get normalized counts using VST
  vst_counts <- DESeq2::vst(dds, blind = FALSE)
  normalized_counts <- SummarizedExperiment::assay(vst_counts)

  # Get experimental and control samples
  experimental <- gsub(paste0(comparison_info$variable, "_"), "", comparison_info$experimental)
  control <- gsub(paste0(comparison_info$variable, "_"), "", comparison_info$control)

  experimental_samples <- rownames(colData(dds))[colData(dds)[[comparison_info$variable]] == experimental]
  control_samples <- rownames(colData(dds))[colData(dds)[[comparison_info$variable]] == control]

  # Convert DESeqResults to data frame
  res_df <- as.data.frame(deseq_result) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::filter(!is.na(padj), padj < min_padj, abs(log2FoldChange) >= min_log2fc) %>%
    dplyr::select(Gene, log2FoldChange, pvalue, padj)

  if (direction == "up") {
    res_df <- dplyr::filter(res_df, log2FoldChange > 0)
  } else if (direction == "down") {
    res_df <- dplyr::filter(res_df, log2FoldChange < 0)
  }

  res_df <- dplyr::arrange(res_df, padj)

  if (!is.null(n_genes)) {
    res_df <- head(res_df, n_genes)
  }

  if (nrow(res_df) > 0) {
    exp_means <- rowMeans(normalized_counts[res_df$Gene, experimental_samples, drop = FALSE])
    ctrl_means <- rowMeans(normalized_counts[res_df$Gene, control_samples, drop = FALSE])

    res_df[[paste0("Mean (", experimental, ")")]] <- exp_means
    res_df[[paste0("Mean (", control, ")")]] <- ctrl_means

    res_df <- res_df %>%
      dplyr::select(
        Gene,
        paste0("Mean (", experimental, ")"),
        paste0("Mean (", control, ")"),
        everything()
      )

    names(res_df) <- gsub("log2FoldChange", "Log2 FC", names(res_df))
    names(res_df) <- gsub("pvalue", "P-value", names(res_df))
    names(res_df) <- gsub("padj", "Adj. P-value", names(res_df))

    numeric_cols <- sapply(res_df, is.numeric)
    res_df[numeric_cols] <- lapply(res_df[numeric_cols], round, digits = 3)

    return(res_df)
  }

  return(NULL)
}

create_complex_de_table <- function(result, direction = "both", n_genes = NULL,
                                    min_padj = 0.05, min_log2fc = 0, contrast_name) {
  if (!is(result$dds, "DESeqDataSet")) {
    stop("Expected a DESeqDataSet object")
  }

  # Extract the DESeqDataSet and its metadata
  dds <- result$dds
  sample_info <- colData(dds)

  # Get the groups directly from the contrast object
  groups <- extract_contrast_groups(result$comparison_info$contrasts[[contrast_name]])
  if (length(groups) < 2) {
    stop("Could not identify comparison groups")
  }

  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)

  # Get the results
  res <- result$deseq_result

  # Convert to data frame
  res_df <- data.frame(
    Gene = rownames(res),
    Log2_FC = res$log2FoldChange,
    P_value = res$pvalue,
    Adj_P_value = res$padj,
    stringsAsFactors = FALSE
  )

  # Add filtering and sorting
  res_df <- res_df %>%
    dplyr::filter(!is.na(Log2_FC)) %>%
    dplyr::filter(abs(Log2_FC) >= min_log2fc) %>%
    dplyr::filter(!is.na(Adj_P_value) & Adj_P_value < min_padj)

  if (direction == "up") {
    res_df <- dplyr::filter(res_df, Log2_FC > 0)
  } else if (direction == "down") {
    res_df <- dplyr::filter(res_df, Log2_FC < 0)
  }

  res_df <- res_df %>%
    dplyr::arrange(Adj_P_value, desc(abs(Log2_FC)))

  if (!is.null(n_genes)) {
    res_df <- head(res_df, n_genes)
  }

  if (nrow(res_df) > 0) {
    # Calculate means for both groups using the extracted group names
    for (group in groups) {
      # Split group name back into condition and timepoint
      group_parts <- strsplit(group, "_")[[1]]
      condition <- paste(head(group_parts, -1), collapse = "_")
      timepoint <- tail(group_parts, 1)

      # Get samples for this group
      samples <- rownames(sample_info)[
        sample_info$condition == condition &
          sample_info$timepoint == timepoint
      ]

      if (length(samples) > 0) {
        means <- rowMeans(norm_counts[res_df$Gene, samples, drop = FALSE])
        res_df[[paste0("Mean (", group, ")")]] <- means
      }
    }

    # Reorder columns
    res_df <- res_df %>%
      dplyr::select(
        Gene,
        Log2_FC,
        starts_with("Mean"),
        P_value,
        Adj_P_value
      )

    # Rename columns for clarity
    names(res_df) <- gsub("Log2_FC", "Log2 FC", names(res_df))
    names(res_df) <- gsub("P_value", "P-value", names(res_df))
    names(res_df) <- gsub("Adj_P_value", "Adj. P-value", names(res_df))

    # Round numeric columns
    numeric_cols <- sapply(res_df, is.numeric)
    res_df[numeric_cols] <- lapply(res_df[numeric_cols], round, digits = 3)

    return(res_df)
  }

  return(NULL)
}

generate_de_tables <- function(result, experiment_obj, direction = "both",
                               n_genes = 25, min_padj = 0.05, min_log2fc = 0,
                               contrast_name, analysis_type = "simple") {
  tables_dir <- file.path(
    "results", experiment_obj$metadata$config$name, "de", "deseq2",
    analysis_type, "tables"
  )
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  # Use appropriate function based on analysis type
  table_data <- if (analysis_type == "complex") {
    create_complex_de_table(result, direction, n_genes, min_padj, min_log2fc, contrast_name) # Pass contrast_name
  } else {
    create_de_table(result, direction, n_genes, min_padj, min_log2fc)
  }

  if (!is.null(table_data)) {
    filename <- file.path(
      tables_dir,
      paste0(
        contrast_name, "_", direction, "_top", n_genes,
        "_padj", min_padj, "_log2fc", min_log2fc, ".docx"
      )
    )

    ft <- save_table_to_docx(table_data, filename)
    return(filename)
  } else {
    warning("No table data generated for ", contrast_name)
    return(NULL)
  }
}

# Create volcano plot data
create_volcano_data <- function(result) {
  # Extract results
  res <- result$deseq_result

  # Create data frame
  plot_data <- data.frame(
    gene = rownames(res),
    log2FoldChange = res$log2FoldChange,
    padj = res$padj,
    stringsAsFactors = FALSE
  )

  # Calculate -log10(padj)
  plot_data$neg_log10_padj <- -log10(plot_data$padj)

  # Create significance and direction columns
  plot_data$significant <- ifelse(plot_data$padj < 0.05, "Significant", "Not Significant")
  plot_data$direction <- paste0(
    ifelse(plot_data$log2FoldChange > 0, "Increased", "Decreased"),
    " & ",
    plot_data$significant
  )

  return(plot_data)
}

create_volcano_plot <- function(plot_data, title, n_labels = 15, analysis_type = "simple") {
  # Set colors
  increased_color <- "#E16A54" # salmon/coral
  decreased_color <- "#7C444F" # burgundy

  # Calculate plot dimensions
  width <- 10
  height <- 8
  label_size <- 3.5

  # Calculate axis limits
  x_limit <- max(abs(plot_data$log2FoldChange), na.rm = TRUE) * 1.1
  y_limit <- max(plot_data$neg_log10_padj, na.rm = TRUE) * 1.1

  # Create direction and significance categories
  plot_data$category <- ifelse(plot_data$padj < 0.05,
    ifelse(plot_data$log2FoldChange > 0,
      "Significantly Increased",
      "Significantly Decreased"
    ),
    ifelse(plot_data$log2FoldChange > 0,
      "Increased (Not Significant)",
      "Decreased (Not Significant)"
    )
  )

  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = log2FoldChange,
    y = neg_log10_padj
  )) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = category,
        color = category,
        shape = category
      ),
      size = 2, alpha = 0.7
    ) +
    ggplot2::scale_shape_manual(
      name = "Expression Change",
      values = c(
        "Significantly Increased" = 21, # filled circle
        "Significantly Decreased" = 21, # filled circle
        "Increased (Not Significant)" = 1, # open circle
        "Decreased (Not Significant)" = 1 # open circle
      )
    ) +
    ggplot2::scale_fill_manual(
      name = "Expression Change",
      values = c(
        "Significantly Increased" = increased_color,
        "Significantly Decreased" = decreased_color,
        "Increased (Not Significant)" = "white",
        "Decreased (Not Significant)" = "white"
      )
    ) +
    ggplot2::scale_color_manual(
      name = "Expression Change",
      values = c(
        "Significantly Increased" = increased_color,
        "Significantly Decreased" = decreased_color,
        "Increased (Not Significant)" = increased_color,
        "Decreased (Not Significant)" = decreased_color
      )
    ) +
    # Dynamic axis limits and theme
    ggplot2::scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = seq(-floor(x_limit), floor(x_limit), by = 2)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_limit),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] ~ "(adjusted p-value)")
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans"),
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = 16 * sqrt(width / 10),
        face = "bold",
        margin = ggplot2::margin(b = 20)
      ),
      plot.margin = ggplot2::margin(20, 20, 20, 20),
      axis.text = ggplot2::element_text(
        size = 12 * sqrt(width / 10),
        face = "bold", color = "black"
      ),
      axis.title.x = ggplot2::element_text(
        size = 14 * sqrt(width / 10),
        face = "bold", color = "black",
        margin = ggplot2::margin(t = 30)
      ),
      axis.title.y = ggplot2::element_text(
        size = 14 * sqrt(width / 10),
        face = "bold", color = "black",
        margin = ggplot2::margin(r = 30)
      ),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = 12 * sqrt(width / 10),
        margin = ggplot2::margin(b = 15)
      ),
      legend.title = ggplot2::element_text(
        size = 14 * sqrt(width / 10),
        face = "bold"
      ),
      legend.position = "right",
      legend.spacing.y = unit(1, "cm"),
      legend.key.size = unit(1, "cm"),
      legend.box.margin = ggplot2::margin(l = 10)
    )

  # Get significant genes
  sig_genes <- plot_data[plot_data$padj < 0.05, ]

  # Get top 15 up-regulated genes (sorted by significance then fold change)
  up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
  up_genes <- up_genes[order(up_genes$padj, -up_genes$log2FoldChange)[1:min(15, nrow(up_genes))], ]

  # Get top 15 down-regulated genes (sorted by significance then fold change)
  down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]
  down_genes <- down_genes[order(down_genes$padj, down_genes$log2FoldChange)[1:min(15, nrow(down_genes))], ]

  # Combine top up and down genes
  top_genes <- rbind(up_genes, down_genes)

  # Add gene labels
  if (nrow(top_genes) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(label = gene),
      size = label_size,
      fontface = "bold",
      box.padding = 1.5,
      point.padding = 0.5,
      min.segment.length = 0.1,
      force = 3,
      max.overlaps = Inf, # Allow all labels to be shown
      segment.color = "grey50",
      segment.size = 0.3,
      nudge_y = ifelse(top_genes$neg_log10_padj < 2, 1.5,
        ifelse(top_genes$neg_log10_padj < 4, 0.5, 0.1)
      ),
      nudge_x = ifelse(top_genes$log2FoldChange > 0, 1.5, -1.5),
      ylim = c(NA, NA),
      hjust = 0.5,
      vjust = 0.5
    )
  }

  return(list(plot = p, width = width, height = height))
}

# Update generate_volcano_plot to use the new dimensions
generate_volcano_plot <- function(result, experiment_obj, n_labels = 15,
                                  contrast_name, analysis_type = "simple") {
  # Create volcano plots directory with updated path
  plots_dir <- file.path(
    "results", experiment_obj$metadata$config$name, "de", "deseq2",
    analysis_type, "plots", "volcano"
  )
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  plot_data <- create_volcano_data(result)
  plot_result <- create_volcano_plot(
    plot_data, paste("Volcano Plot:", contrast_name),
    n_labels, analysis_type
  )

  filename <- file.path(
    plots_dir,
    paste0(contrast_name, "_volcano_labels", n_labels, ".pdf")
  )

  ggplot2::ggsave(filename, plot_result$plot,
    width = plot_result$width, height = plot_result$height
  )
  return(plot_result$plot)
}

create_heatmap_data <- function(deseq_results, min_padj = 0.05, direction = "both",
                                n_genes = NULL, scale = TRUE) {
  vst_counts <- DESeq2::vst(deseq_results$dds, blind = FALSE)
  norm_counts <- SummarizedExperiment::assay(vst_counts)
  res_df <- as.data.frame(deseq_results$deseq_result)

  # Filter by direction and significance
  if (direction == "up") {
    sig_genes <- rownames(res_df)[!is.na(res_df$padj) & res_df$padj < min_padj &
      res_df$log2FoldChange > 0]
  } else if (direction == "down") {
    sig_genes <- rownames(res_df)[!is.na(res_df$padj) & res_df$padj < min_padj &
      res_df$log2FoldChange < 0]
  } else {
    sig_genes <- rownames(res_df)[!is.na(res_df$padj) & res_df$padj < min_padj]
  }

  # Order by padj
  sig_genes <- sig_genes[order(res_df[sig_genes, "padj"])]

  # Take top n genes if specified
  if (!is.null(n_genes)) {
    sig_genes <- head(sig_genes, n_genes)
  }

  if (length(sig_genes) == 0) {
    stop("No genes meet the specified criteria")
  }

  sig_counts <- norm_counts[sig_genes, , drop = FALSE]

  if (scale) {
    scaled_counts <- t(scale(t(sig_counts), center = TRUE, scale = TRUE))
  } else {
    scaled_counts <- sig_counts
  }

  return(list(
    counts = scaled_counts,
    n_genes = length(sig_genes),
    gene_names = sig_genes
  ))
}

create_complex_heatmap <- function(heatmap_data, result, title, n_clusters = 3, analysis_type = "simple") {
  sample_info <- SummarizedExperiment::colData(result$dds)

  # Get condition from design formula more safely
  design_formula <- DESeq2::design(result$dds)
  condition <- all.vars(design_formula)[1] # Get first variable from design formula

  # Create a larger color palette that can handle up to 12 conditions
  condition_levels <- levels(sample_info[[condition]])
  n_colors_needed <- length(condition_levels)

  # Combine multiple color palettes for more unique colors
  condition_colors <- c(
    RColorBrewer::brewer.pal(8, "Set2"),
    RColorBrewer::brewer.pal(8, "Set1"),
    RColorBrewer::brewer.pal(8, "Dark2")
  )[1:n_colors_needed]

  names(condition_colors) <- condition_levels

  # Clean up column names by replacing underscores with spaces
  colnames(heatmap_data$counts) <- gsub("_", " ", colnames(heatmap_data$counts))

  # Calculate dimensions based on data size
  n_rows <- nrow(heatmap_data$counts)
  n_cols <- ncol(heatmap_data$counts)

  # Calculate base dimensions that will determine cell size
  base_width <- 12
  base_height <- 14

  # Adjust based on number of samples and genes using log scaling
  width_factor <- 1 + log10(n_cols) / 2
  height_factor <- 1 + log10(n_rows) / 2

  # Calculate final dimensions
  width <- min(max(base_width * width_factor, 8), 20) # Between 8 and 20 inches
  height <- min(max(base_height * height_factor, 10), 24) # Between 10 and 24 inches

  # Calculate font sizes based on number of rows/columns
  base_font_size <- 14 / (1 + log10(max(n_rows, n_cols)))

  # Scale different text elements relative to base font size
  title_font <- max(12, base_font_size * 2.0)
  gene_font_size <- max(8, base_font_size * 1.2)
  sample_font_size <- max(8, base_font_size * 1.2)
  legend_font_size <- max(8, base_font_size * 1.1)
  annotation_font_size <- max(8, base_font_size * 1.1)

  # Create annotation with scaled styling
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    Condition = sample_info[[condition]],
    col = list(Condition = condition_colors),
    annotation_name_gp = grid::gpar(
      fontsize = annotation_font_size,
      fontface = "bold",
      family = "sans"
    ),
    annotation_legend_param = list(
      Condition = list(
        title_gp = grid::gpar(fontsize = legend_font_size, fontface = "bold", family = "sans"),
        labels_gp = grid::gpar(fontsize = legend_font_size, family = "sans"),
        ncol = ceiling(length(condition_levels) / 12)
      )
    )
  )

  use_raster <- nrow(heatmap_data$counts) > 2000

  # Create the heatmap with adjusted parameters
  hm <- ComplexHeatmap::Heatmap(
    heatmap_data$counts,
    name = "Expression",

    # Title settings
    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = title_font,
      fontface = "bold",
      family = "sans"
    ),

    # Remove row title
    row_title = NULL,

    # Clustering settings
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",

    # Text settings with dynamic sizing
    row_names_gp = grid::gpar(
      fontsize = gene_font_size,
      family = "sans"
    ),
    column_names_gp = grid::gpar(
      fontsize = sample_font_size,
      family = "sans"
    ),

    # Layout settings
    row_names_max_width = grid::unit(8, "cm"),
    show_row_names = TRUE,
    row_names_side = "right",
    show_row_dend = FALSE,
    show_column_dend = FALSE,

    # Column names rotation
    column_names_rot = 45,

    # Raster settings for large matrices
    use_raster = use_raster,
    raster_quality = 2,

    # Annotation
    top_annotation = column_ha,

    # Legend settings
    heatmap_legend_param = list(
      title_gp = grid::gpar(
        fontsize = legend_font_size,
        fontface = "bold",
        family = "sans"
      ),
      labels_gp = grid::gpar(
        fontsize = legend_font_size,
        family = "sans"
      ),
      legend_height = grid::unit(4, "cm"),
      legend_width = grid::unit(1.5, "cm"),
      title_position = "leftcenter-rot"
    ),

    # Color scheme - black and white gradient
    col = circlize::colorRamp2(
      c(
        min(heatmap_data$counts),
        mean(range(heatmap_data$counts)),
        max(heatmap_data$counts)
      ),
      c("white", "gray50", "black")
    ),

    # Set width and height to control cell size
    width = grid::unit(width * 50, "points"),
    height = grid::unit(height * 50, "points")
  )

  return(list(heatmap = hm, width = width, height = height))
}

create_complex_deseq_heatmap <- function(heatmap_data, result, title, n_clusters = 3) {
  # Get the groups using extract_contrast_groups
  groups <- extract_contrast_groups(result$comparison_info$contrasts[[title]])
  if (length(groups) < 2) {
    stop("Could not identify comparison groups")
  }

  # Get sample information
  sample_info <- SummarizedExperiment::colData(result$dds)

  # Split groups to get condition and timepoint
  sample_matches <- lapply(groups, function(group) {
    group_parts <- strsplit(group, "_")[[1]]
    condition <- paste(head(group_parts, -1), collapse = "_")
    timepoint <- tail(group_parts, 1)
    sample_info$condition == condition & sample_info$timepoint == timepoint
  })

  # Combine all relevant samples
  relevant_samples <- Reduce(`|`, sample_matches)
  filtered_sample_info <- sample_info[relevant_samples, ]
  filtered_counts <- heatmap_data$counts[, relevant_samples]

  # Get the factor name from the design formula
  design_formula <- DESeq2::design(result$dds)
  factor_name <- all.vars(design_formula)[1]

  # Create a larger color palette that can handle up to 12 conditions
  condition_levels <- levels(filtered_sample_info[[factor_name]])
  n_colors_needed <- length(condition_levels)

  # Combine multiple color palettes for more unique colors
  condition_colors <- c(
    RColorBrewer::brewer.pal(8, "Set2"),
    RColorBrewer::brewer.pal(8, "Set1"),
    RColorBrewer::brewer.pal(8, "Dark2")
  )[1:n_colors_needed]

  names(condition_colors) <- condition_levels

  # Clean up column names by replacing underscores with spaces
  colnames(filtered_counts) <- gsub("_", " ", colnames(filtered_counts))

  # Calculate dimensions based on data size
  n_rows <- nrow(filtered_counts)
  n_cols <- ncol(filtered_counts)

  # Calculate base dimensions that will determine cell size
  base_width <- 12
  base_height <- 14

  # Adjust based on number of samples and genes using log scaling
  width_factor <- 1 + log10(n_cols) / 2
  height_factor <- 1 + log10(n_rows) / 2

  # Calculate final dimensions
  width <- min(max(base_width * width_factor, 8), 20) # Between 8 and 20 inches
  height <- min(max(base_height * height_factor, 10), 24) # Between 10 and 24 inches

  # Calculate font sizes based on number of rows/columns
  base_font_size <- 14 / (1 + log10(max(n_rows, n_cols)))

  # Scale different text elements relative to base font size
  title_font <- max(12, base_font_size * 2.0)
  gene_font_size <- max(8, base_font_size * 1.2)
  sample_font_size <- max(8, base_font_size * 1.2)
  legend_font_size <- max(8, base_font_size * 1.1)
  annotation_font_size <- max(8, base_font_size * 1.1)

  # Create annotation with scaled styling
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    Group = factor(sapply(seq_len(ncol(filtered_counts)), function(i) {
      # Find which group this sample belongs to
      for (group in groups) {
        group_parts <- strsplit(group, "_")[[1]]
        condition <- paste(head(group_parts, -1), collapse = "_")
        timepoint <- tail(group_parts, 1)
        if (filtered_sample_info$condition[i] == condition &&
          filtered_sample_info$timepoint[i] == timepoint) {
          return(group)
        }
      }
      return(NA)
    }), levels = groups),
    col = list(Group = setNames(condition_colors[1:length(groups)], groups)),
    annotation_name_gp = grid::gpar(
      fontsize = annotation_font_size,
      fontface = "bold",
      family = "sans"
    ),
    annotation_legend_param = list(
      Group = list(
        title_gp = grid::gpar(fontsize = legend_font_size, fontface = "bold", family = "sans"),
        labels_gp = grid::gpar(fontsize = legend_font_size, family = "sans"),
        ncol = ceiling(length(groups) / 12)
      )
    )
  )

  use_raster <- nrow(filtered_counts) > 2000

  # Create the heatmap with adjusted parameters
  hm <- ComplexHeatmap::Heatmap(
    filtered_counts,
    name = "Expression",

    # Title settings
    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = title_font,
      fontface = "bold",
      family = "sans"
    ),

    # Remove row title
    row_title = NULL,

    # Clustering settings
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",

    # Text settings with dynamic sizing
    row_names_gp = grid::gpar(
      fontsize = gene_font_size,
      family = "sans"
    ),
    column_names_gp = grid::gpar(
      fontsize = sample_font_size,
      family = "sans"
    ),

    # Layout settings
    row_names_max_width = grid::unit(8, "cm"),
    show_row_names = TRUE,
    row_names_side = "right",
    show_row_dend = FALSE,
    show_column_dend = FALSE,

    # Column names rotation
    column_names_rot = 45,

    # Raster settings for large matrices
    use_raster = use_raster,
    raster_quality = 2,

    # Annotation
    top_annotation = column_ha,

    # Legend settings
    heatmap_legend_param = list(
      title_gp = grid::gpar(
        fontsize = legend_font_size,
        fontface = "bold",
        family = "sans"
      ),
      labels_gp = grid::gpar(
        fontsize = legend_font_size,
        family = "sans"
      ),
      legend_height = grid::unit(4, "cm"),
      legend_width = grid::unit(1.5, "cm"),
      title_position = "leftcenter-rot"
    ),

    # Color scheme - black and white gradient
    col = circlize::colorRamp2(
      c(
        min(filtered_counts),
        mean(range(filtered_counts)),
        max(filtered_counts)
      ),
      c("white", "gray50", "black")
    ),

    # Set width and height to control cell size
    width = grid::unit(width * 50, "points"),
    height = grid::unit(height * 50, "points")
  )

  return(list(heatmap = hm, width = width, height = height))
}

generate_complex_heatmap <- function(result, experiment_obj, direction = "both",
                                     n_genes = 25, scale = TRUE, n_clusters = 3,
                                     min_padj = 0.05, min_log2fc = 0,
                                     contrast_name, analysis_type = "simple") {
  # Create base directory path
  base_dir <- file.path("results", experiment_obj$metadata$config$name)
  plots_dir <- file.path(base_dir, "de", "deseq2", analysis_type, "plots", "heatmap")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Create heatmap data
  heatmap_data <- create_heatmap_data(result, min_padj, direction, n_genes, scale)

  if (analysis_type == "complex") {
    # Get the groups for this contrast
    groups <- extract_contrast_groups(result$comparison_info$contrasts[[contrast_name]])

    # Create result object with the groups
    result_for_heatmap <- list(
      dds = result$dds,
      deseq_result = result$deseq_results[[contrast_name]],
      comparison_info = list(
        contrasts = list(),
        groups = groups
      )
    )
    result_for_heatmap$comparison_info$contrasts[[contrast_name]] <- result$comparison_info$contrasts[[contrast_name]]

    hm_result <- create_complex_deseq_heatmap(
      heatmap_data = heatmap_data,
      result = result_for_heatmap,
      title = contrast_name,
      n_clusters = n_clusters
    )
  } else {
    # Simple analysis remains unchanged
    hm_result <- create_complex_heatmap(
      heatmap_data = heatmap_data,
      result = result,
      title = contrast_name,
      n_clusters = n_clusters,
      analysis_type = analysis_type
    )
  }

  # Create filename and save
  filename <- file.path(
    plots_dir,
    paste0(
      contrast_name,
      "_",
      direction,
      "_top",
      n_genes,
      "_padj",
      min_padj,
      "_log2fc",
      min_log2fc,
      "_heatmap.pdf"
    )
  )

  # Save the heatmap
  tryCatch(
    {
      pdf(filename, width = hm_result$width, height = hm_result$height)
      ComplexHeatmap::draw(hm_result$heatmap)
      dev.off()
    },
    error = function(e) {
      if (dev.cur() > 1) dev.off()
      stop("Error in heatmap generation: ", e$message)
    }
  )

  return(hm_result$heatmap)
}

save_table_to_docx <- function(table_data, filename) {
  if (is.null(table_data)) {
    warning("No table data provided")
    return(NULL)
  }

  # Sort by absolute fold change and adjusted p-value
  if (all(c("Log2 FC", "Adj. P-value") %in% colnames(table_data))) {
    table_data <- table_data %>%
      dplyr::arrange(dplyr::desc(abs(`Log2 FC`)), `Adj. P-value`)
  } else if (all(c("log2FoldChange", "padj") %in% colnames(table_data))) {
    table_data <- table_data %>%
      dplyr::arrange(dplyr::desc(abs(log2FoldChange)), padj)
  }

  # Calculate dimensions
  n_rows <- nrow(table_data)
  n_cols <- ncol(table_data)

  # Limit table size if too large
  max_rows <- 1000 # Maximum rows to prevent Word crashes
  if (n_rows > max_rows) {
    warning(sprintf("Table truncated to %d rows", max_rows))
    table_data <- head(table_data, max_rows)
    n_rows <- max_rows
  }

  # Calculate optimal widths based on content
  col_widths <- sapply(table_data, function(col) {
    if (is.numeric(col)) {
      return(1) # Numeric columns get base width
    } else {
      max_length <- max(nchar(as.character(col)), nchar(names(table_data)[which(table_data == col)[1]]))
      return(min(3, max(1, max_length / 10))) # Scale width based on content length
    }
  })

  # Calculate font size based on number of rows (between 8 and 11)
  font_size <- max(8, min(11, 14 - log2(n_rows)))

  # Create flextable with smart sizing
  ft <- flextable::flextable(table_data) %>%
    flextable::theme_vanilla() %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = font_size, part = "all") %>%
    flextable::padding(padding = 2, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(align = "left", j = 1) # Left align first column (Gene names)

  # Format numeric columns
  numeric_cols <- sapply(table_data, is.numeric)
  if (any(numeric_cols)) {
    ft <- flextable::colformat_double(ft,
      j = which(numeric_cols),
      digits = 3,
      big.mark = "",
      na_str = "NA"
    )
  }

  # Set optimal widths
  for (i in seq_along(col_widths)) {
    ft <- flextable::width(ft, j = i, width = col_widths[i])
  }

  # Add borders
  ft <- flextable::border(ft,
    border.top = officer::fp_border(color = "gray80"),
    border.bottom = officer::fp_border(color = "gray80"),
    part = "all"
  )

  # Set header style
  ft <- flextable::bold(ft, part = "header") %>%
    flextable::bg(bg = "gray95", part = "header")

  # Calculate page dimensions
  n_cols_total <- sum(col_widths)
  page_width <- min(8.5, max(6, n_cols_total * 1.2)) # in inches
  page_height <- min(11, max(8, n_rows * 0.2)) # in inches

  # Save to docx with calculated dimensions
  tryCatch(
    {
      flextable::save_as_docx(
        ft,
        path = filename,
        pr_section = officer::prop_section(
          page_size = officer::page_size(
            orient = ifelse(page_width > page_height, "landscape", "portrait"),
            width = page_width,
            height = page_height
          ),
          type = "continuous",
          page_margins = officer::page_mar(
            bottom = 0.5,
            top = 0.5,
            right = 0.5,
            left = 0.5
          )
        )
      )
    },
    error = function(e) {
      warning(sprintf("Failed to save docx: %s. Trying with simplified formatting...", e$message))

      # Fallback with minimal formatting
      ft_simple <- flextable::flextable(table_data) %>%
        flextable::theme_vanilla() %>%
        flextable::fontsize(size = 9, part = "all") %>%
        flextable::autofit()

      flextable::save_as_docx(ft_simple, path = filename)
    }
  )

  return(ft)
}

# 3. Set up export parameters
export_params <- list(
  de_params = list(min_padj = 0.05, min_log2fc = 0),
  table_params = list(direction = "both", n_genes = 25),
  volcano_params = list(n_labels = 15),
  heatmap_params = list(direction = "both", n_genes = 25, scale = TRUE, n_clusters = 3)
)

# 4. Run the export
xen_tran_2024_03 <- export_all_result_types(
  # dea_results = xen_tran_2024_03_de_results,
  experiment_name = "xen_tran_2024_03",
  de_params = export_params$de_params,
  table_params = export_params$table_params,
  volcano_params = export_params$volcano_params,
  heatmap_params = export_params$heatmap_params,
  analysis_type = "simple"
)

# 4. Run the export
xen_tran_2024_12 <- export_all_result_types(
  # dea_results = xen_tran_2024_12_de_results,
  experiment_name = "xen_tran_2024_12",
  de_params = export_params$de_params,
  table_params = export_params$table_params,
  volcano_params = export_params$volcano_params,
  heatmap_params = export_params$heatmap_params,
  analysis_type = "simple"
)

# 4. Run the export
xen_tran_2024_12 <- export_all_result_types(
  # dea_results = xen_tran_2024_12_de_results_complex,
  experiment_name = "xen_tran_2024_12",
  de_params = export_params$de_params,
  table_params = export_params$table_params,
  volcano_params = export_params$volcano_params,
  heatmap_params = export_params$heatmap_params,
  analysis_type = "complex"
)
# str(xen_tran_2024_12_de_results_complex[1])
#
# # Add this diagnostic code before running export_all_result_types
# print("Detailed structure of complex results:")
# if (exists("xen_tran_2024_12_de_results_complex")) {
#     # Print overall structure
#     print("Top level structure:")
#     str(xen_tran_2024_12_de_results_complex, max.level = 1)
#
#     # Check if it's a list with the required components
#     if (is.list(xen_tran_2024_12_de_results_complex)) {
#         print("\nContains deseq_results?")
#         print("deseq_results" %in% names(xen_tran_2024_12_de_results_complex))
#
#         if ("deseq_results" %in% names(xen_tran_2024_12_de_results_complex)) {
#             print("\nStructure of deseq_results:")
#             str(xen_tran_2024_12_de_results_complex$deseq_results, max.level = 1)
#         }
#
#         print("\nContains dds?")
#         print("dds" %in% names(xen_tran_2024_12_de_results_complex))
#
#         if ("dds" %in% names(xen_tran_2024_12_de_results_complex)) {
#             print("\nClass of dds object:")
#             print(class(xen_tran_2024_12_de_results_complex$dds))
#         }
#     }
# } else {
#     print("xen_tran_2024_12_de_results_complex not found!")
# }
#
# colnames(colData(xen_tran_2024_12_de_results[[1]]$dds))
#
#
# # Add this diagnostic code before running the simple analysis
# print("Checking structure of simple results:")
# if (exists("xen_tran_2024_03_de_results")) {
#     print("Top level structure:")
#     str(xen_tran_2024_03_de_results, max.level = 1)
#
#     # Check first contrast
#     if (length(xen_tran_2024_03_de_results) > 0) {
#         print("\nStructure of first contrast:")
#         str(xen_tran_2024_03_de_results[[1]], max.level = 1)
#
#         print("\nComponents in first contrast:")
#         print(names(xen_tran_2024_03_de_results[[1]]))
#
#         if ("deseq_result" %in% names(xen_tran_2024_03_de_results[[1]])) {
#             print("\nClass of deseq_result:")
#             print(class(xen_tran_2024_03_de_results[[1]]$deseq_result))
#         }
#
#         if ("dds" %in% names(xen_tran_2024_03_de_results[[1]])) {
#             print("\nClass of dds:")
#             print(class(xen_tran_2024_03_de_results[[1]]$dds))
#
#             print("\nDesign formula:")
#             print(design(xen_tran_2024_03_de_results[[1]]$dds))
#
#             print("\nColumn names in colData:")
#             print(colnames(colData(xen_tran_2024_03_de_results[[1]]$dds)))
#         }
#     }
# } else {
#     print("xen_tran_2024_03_de_results not found!")
# }
#
# # Add this diagnostic line
# print("Column names in colData:")
# print(colnames(colData(xen_tran_2024_03_de_results[[1]]$dds)))
