#!/usr/bin/env Rscript

# =============================================================================
# Step 07: Clustering Resolution Stability (Clustree)
# Tests multiple resolutions and visualizes cluster splitting/merging behavior
# to help justify the final resolution choice.
# =============================================================================

library(Seurat)
library(clustree)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args        <- commandArgs(trailingOnly = TRUE)
input_rds   <- args[1]
output_plot <- args[2]

SEED <- 42
resolutions <- seq(0.2, 1.2, by = 0.2)

message(glue("
========== Resolution Stability (Step 07) ==========
  Input       : {input_rds}
  Output      : {output_plot}
  Resolutions : {paste(resolutions, collapse = ', ')}
  Seed        : {SEED}
====================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

seurat_obj <- readRDS(input_rds)
message(glue("Loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes"))

# =============================================================================
# 3. Cluster at Multiple Resolutions
# =============================================================================

message("Clustering across resolutions...")
seurat_obj <- FindClusters(
    object      = seurat_obj,
    resolution  = resolutions,
    random.seed = SEED,
    verbose     = FALSE
)

# =============================================================================
# 4. Summarize Cluster Counts per Resolution
# =============================================================================

cluster_cols <- glue("SCT_snn_res.{resolutions}")

cluster_summary <- tibble(
    resolution = resolutions,
    n_clusters = map_int(cluster_cols, ~ length(unique(seurat_obj@meta.data[[.x]])))
)

message("Cluster count per resolution:")
print(cluster_summary)

# =============================================================================
# 5. Generate Clustree
# =============================================================================

message("Building cluster stability tree...")
tree_plot <- clustree(seurat_obj, prefix = "SCT_snn_res.") +
    ggtitle("Clustering Resolution Stability Tree") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# =============================================================================
# 6. Save
# =============================================================================

ggsave(output_plot, plot = tree_plot, width = 10, height = 12)
message(glue("Clustree saved to: {output_plot}"))
message("Step 07 complete.")