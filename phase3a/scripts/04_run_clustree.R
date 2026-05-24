#!/usr/bin/env Rscript

# =============================================================================
# Step 04: Clustering Resolution Stability (Clustree)
# Input : 03_seurat_clustered.rds — clustered integrated object (RNA_snn graph
#         already built in step 03, so only FindClusters is re-run here)
# Output: 04_clustree_resolutions.pdf — resolution stability tree
# =============================================================================

library(Seurat)
library(clustree)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript 04_run_clustree.R <input_rds> <output_pdf>")
}

input_rds   <- args[1]
output_plot <- args[2]

SEED        <- 42
resolutions <- seq(0.2, 1.2, by = 0.2)

message(glue("
========== Resolution Stability (Step 04) ==========
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
message(glue("Loaded: {ncol(seurat_obj)} cells"))

# =============================================================================
# 3. Cluster at Multiple Resolutions
# =============================================================================

message("Clustering across resolutions (RNA_snn graph from step 03)...")
seurat_obj <- FindClusters(
    object      = seurat_obj,
    resolution  = resolutions,
    random.seed = SEED,
    verbose     = FALSE
)

# =============================================================================
# 4. Summarize Cluster Counts per Resolution
# =============================================================================

cluster_cols <- glue("RNA_snn_res.{resolutions}")

cluster_summary <- tibble(
    resolution = resolutions,
    n_clusters = map_int(cluster_cols, \(col) length(unique(seurat_obj@meta.data[[col]])))
)

message("Cluster count per resolution:")
print(cluster_summary)

# =============================================================================
# 5. Generate Clustree
# =============================================================================

message("Building cluster stability tree...")
tree_plot <- clustree(seurat_obj, prefix = "RNA_snn_res.") +
    ggtitle("Integrated Atlas: Clustering Resolution Stability") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(output_plot, plot = tree_plot, width = 10, height = 12)
message(glue("Clustree saved to: {output_plot}"))
message("Step 04 complete.")
