#!/usr/bin/env Rscript

# =============================================================================
# Step 01: Find Marker Genes
# Identifies differentially expressed genes for each cluster using
# Wilcoxon Rank-Sum tests on the SCT assay.
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)
library(future)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript 01_find_markers.R <input_rds> <output_csv> [cpus] [top_n_csv]")
}

input_rds       <- args[1]
output_csv      <- args[2]
cpus_to_use     <- if (length(args) >= 3) as.numeric(args[3]) else 1
top_markers_csv <- if (length(args) >= 4) args[4] else "01_top_5_markers_per_cluster.csv"

TOP_N <- 5
SEED  <- 42

message(glue("
========== Marker Discovery (Step 01) ==========
  Input  : {input_rds}
  All    : {output_csv}
  Top-{TOP_N}  : {top_markers_csv}
  Cores  : {cpus_to_use}
================================================
"))

# =============================================================================
# 2. Configure Parallelism
# =============================================================================

plan(multisession, workers = cpus_to_use)
options(future.globals.maxSize = 10 * 1024^3)   # 10 GB

# =============================================================================
# 3. Load Object
# =============================================================================

seurat_obj <- readRDS(input_rds)
DefaultAssay(seurat_obj) <- "SCT"
message(glue("Loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes | {length(levels(seurat_obj))} clusters"))

# =============================================================================
# 4. Find All Markers
# =============================================================================

message("Running FindAllMarkers (Wilcoxon, positive markers only)...")
set.seed(SEED)

markers <- FindAllMarkers(
    object          = seurat_obj,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    verbose         = FALSE
)

# =============================================================================
# 5. Summary Diagnostics
# =============================================================================

per_cluster_counts <- markers |>
    as_tibble() |>
    count(cluster, name = "n_markers")

message(glue(
    "Total markers found  : {nrow(markers)}\n",
    "Clusters with markers: {nrow(per_cluster_counts)} / {length(levels(seurat_obj))}\n",
    "Median per cluster   : {median(per_cluster_counts$n_markers)}"
))

# =============================================================================
# 6. Save Full Marker List
# =============================================================================

write_csv(markers, output_csv)
message(glue("All markers saved to: {output_csv}"))

# =============================================================================
# 7. Top-N Cheat Sheet
# =============================================================================

top_genes <- markers |>
    as_tibble() |>
    group_by(cluster) |>
    slice_max(n = TOP_N, order_by = avg_log2FC, with_ties = FALSE) |>
    ungroup()

write_csv(top_genes, top_markers_csv)
message(glue("Top-{TOP_N} markers saved to: {top_markers_csv}"))

# =============================================================================
# 8. Display Signature Markers
# =============================================================================

message(glue("\n{str_dup('-', 50)}"))
message(glue("SIGNATURE MARKERS (TOP {TOP_N} PER CLUSTER)"))
message(str_dup("-", 50))

top_genes |>
    select(cluster, gene, avg_log2FC, p_val_adj) |>
    print(n = Inf)

message("\nStep 01 complete.")
