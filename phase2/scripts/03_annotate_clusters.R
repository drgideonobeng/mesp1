#!/usr/bin/env Rscript

# =============================================================================
# Step 03: Annotate Clusters
# Assigns biological identities to Seurat clusters using a CSV label file.
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args         <- commandArgs(trailingOnly = TRUE)
input_rds    <- args[1]   # Clustered Seurat object
project_name <- args[2]   # Project title for plot
labels_file  <- args[3]   # Path to cluster_labels.csv
output_rds   <- args[4]   # Annotated Seurat object
output_pdf   <- args[5]   # Annotated UMAP plot

message(glue("
========== Cluster Annotation (Step 03) ==========
  Input  : {input_rds}
  Labels : {labels_file}
  RDS    : {output_rds}
  Plot   : {output_pdf}
===================================================
"))

# =============================================================================
# 2. Load Cluster Labels
# =============================================================================

labels_df <- read_csv(labels_file, col_types = cols(.default = "c"))

cluster_identities <- labels_df |>
    deframe()   # Converts a two-column dataframe into a named vector

message(glue("Loaded {nrow(labels_df)} cluster labels from: {labels_file}"))

# =============================================================================
# 3. Load Seurat Object
# =============================================================================

message(glue("Loading Seurat object from: {input_rds}"))
seurat_obj <- readRDS(input_rds)

# =============================================================================
# 4. Safety Checks
# =============================================================================

data_clusters   <- levels(seurat_obj)
labeled_ids     <- names(cluster_identities)

missing_from_data <- setdiff(labeled_ids, data_clusters)
unlabeled         <- setdiff(data_clusters, labeled_ids)

if (length(missing_from_data) > 0) {
    warning(glue(
        "Labels defined in CSV but not found in data: ",
        "{paste(missing_from_data, collapse = ', ')}"
    ))
}

if (length(unlabeled) > 0) {
    warning(glue(
        "Clusters in data but missing from CSV (will keep numeric IDs): ",
        "{paste(unlabeled, collapse = ', ')}"
    ))
}

# =============================================================================
# 5. Apply Annotations
# =============================================================================

seurat_obj <- RenameIdents(seurat_obj, cluster_identities)
seurat_obj$cell_type <- Idents(seurat_obj)

# Sanity check: at least one cluster should have been renamed
if (all(levels(seurat_obj) %in% data_clusters)) {
    warning("No cluster identities were applied — check that cluster IDs in CSV match data clusters.")
}

message("Cluster identities assigned:")
seurat_obj$cell_type |>
    table() |>
    enframe(name = "cell_type", value = "n_cells") |>
    mutate(n_cells = as.integer(n_cells)) |>
    arrange(desc(n_cells)) |>
    print()

# =============================================================================
# 6. UMAP Visualization
# =============================================================================

annotated_umap <- DimPlot(
    seurat_obj,
    reduction = "umap",
    label     = TRUE,
    repel     = TRUE,
    pt.size   = 0.5
) +
    labs(
        title    = glue("{project_name}: Annotated Cell Types"),
        subtitle = glue("{ncol(seurat_obj)} cells | {length(cluster_identities)} populations")
    ) +
    theme_minimal() +
    theme(legend.position = "right")

ggsave(output_pdf, plot = annotated_umap, width = 10, height = 7)
message(glue("UMAP saved to: {output_pdf}"))

# =============================================================================
# 7. Save Annotated Object
# =============================================================================

saveRDS(seurat_obj, file = output_rds)
message(glue("Annotated Seurat object saved to: {output_rds}"))
message("Step 03 complete.")
