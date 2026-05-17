#!/usr/bin/env Rscript

# =============================================================================
# Step 02: Plot Marker Gene Expression
# Generates UMAP feature plots for a user-specified set of marker genes.
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript 02_plot_marker_genes.R <input_rds> <gene_string> <output_pdf>")
}

input_rds   <- args[1]
genes_raw   <- args[2]
output_plot <- args[3]

message(glue("
========== Marker Gene Feature Plots (Step 02) ==========
  Input  : {input_rds}
  Genes  : {genes_raw}
  Output : {output_plot}
=========================================================
"))

# =============================================================================
# 2. Parse and Clean Gene List
# =============================================================================

features_requested <- genes_raw |>
    str_split_1(",") |>
    str_trim() |>
    keep(\(x) x != "")

message(glue("Requested {length(features_requested)} genes: {paste(features_requested, collapse = ', ')}"))

# =============================================================================
# 3. Load Object
# =============================================================================

seurat_obj <- readRDS(input_rds)
message(glue("Loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes"))

# =============================================================================
# 4. Validate Genes Against Dataset
# =============================================================================

final_genes   <- features_requested |> keep(\(g) g %in% rownames(seurat_obj))
skipped_genes <- setdiff(features_requested, final_genes)

if (length(skipped_genes) > 0) {
    warning(glue("Skipped (not in dataset): {paste(skipped_genes, collapse = ', ')}"))
}

if (length(final_genes) == 0) {
    stop("None of the requested genes were found in the dataset. Check your `marker_genes` parameter.")
}

message(glue("Plotting {length(final_genes)} of {length(features_requested)} requested genes."))

# =============================================================================
# 5. Generate Feature Plots
# =============================================================================

n_cols <- min(3, length(final_genes))
n_rows <- ceiling(length(final_genes) / n_cols)

p <- FeaturePlot(
    object   = seurat_obj,
    features = final_genes,
    pt.size  = 0.5,
    order    = TRUE,        # Sparse markers like Mesp1 — plot "hot" cells on top
    ncol     = n_cols
) &
    theme_minimal() &
    theme(plot.title = element_text(face = "bold"))

# =============================================================================
# 6. Save
# =============================================================================

ggsave(
    filename   = output_plot,
    plot       = p,
    width      = 5 * n_cols,
    height     = 5 * n_rows,
    limitsize  = FALSE
)

message(glue("Feature plots saved to: {output_plot}"))
message("Step 02 complete.")
