#!/usr/bin/env Rscript

# =============================================================================
# Step 06: Find All Cluster Markers
# Input : 05_atlas_clean.rds — cleaned atlas (cluster 3 recovered)
# Output: 06_all_markers.csv              — full FindAllMarkers table
#         06_top5_markers_per_cluster.csv — top 5 by avg_log2FC per cluster
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript 06_find_markers.R <input_rds> <out_all_markers> <out_top5>")
}

input_rds   <- args[1]
output_all  <- args[2]
output_top5 <- args[3]

SEED <- 42

message(glue("
========== Find Cluster Markers (Step 06) ==========
  Input    : {input_rds}
  All mrkrs: {output_all}
  Top 5    : {output_top5}
  Seed     : {SEED}
=====================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

atlas_clean <- readRDS(input_rds)
atlas_clean <- JoinLayers(atlas_clean)

n_cells    <- ncol(atlas_clean)
n_clusters <- nlevels(Idents(atlas_clean))

message(glue("Loaded: {n_cells} cells, {n_clusters} clusters"))

message("\n--- Cluster sizes ---")
print(table(cluster = Idents(atlas_clean)))

# =============================================================================
# 3. FindAllMarkers
# =============================================================================

message("Running FindAllMarkers...")
set.seed(SEED)
markers <- FindAllMarkers(
    atlas_clean,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    verbose         = FALSE
)

# =============================================================================
# 4. Save Full Marker Table
# =============================================================================

write.csv(markers, output_all, row.names = FALSE)
message(glue("Full marker table saved: {nrow(markers)} rows → {output_all}"))

# =============================================================================
# 5. Compute Top 5 per Cluster
# =============================================================================

top5 <- markers |>
    group_by(cluster) |>
    slice_max(avg_log2FC, n = 5, with_ties = FALSE) |>
    select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

# =============================================================================
# 6. Save Top 5
# =============================================================================

write.csv(top5, output_top5, row.names = FALSE)
message(glue("Top-5 table saved: {nrow(top5)} rows → {output_top5}"))

# =============================================================================
# 7. Print Top 5 to Log
# =============================================================================

message("\n--- Top 5 markers per cluster ---")
print(top5, n = Inf)

# =============================================================================
# 8. Summary
# =============================================================================

n_clusters_with_markers <- n_distinct(markers$cluster)
n_sig_markers           <- nrow(markers)

message(glue("
========== Marker Summary ==========
  Clusters with markers : {n_clusters_with_markers} / {n_clusters}
  Total significant     : {n_sig_markers}
=====================================
"))

message("Step 06 complete.")
