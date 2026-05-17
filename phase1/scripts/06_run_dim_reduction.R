#!/usr/bin/env Rscript

# =============================================================================
# Step 06: Dimensionality Reduction and Clustering
# Runs PCA, finds clusters (Louvain), and embeds cells in UMAP space.
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args        <- commandArgs(trailingOnly = TRUE)
input_rds   <- args[1]
pca_dims    <- as.numeric(args[2])
cluster_res <- as.numeric(args[3])
umap_title  <- args[4]
output_rds  <- args[5]
elbow_plot  <- args[6]
output_plot <- args[7]

SEED <- 42

message(glue("
========== Dim Reduction & Clustering (Step 06) ==========
  Input       : {input_rds}
  pca_dims    : {pca_dims}
  cluster_res : {cluster_res}
  Seed        : {SEED}
==========================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

seurat_obj <- readRDS(input_rds)
message(glue("Loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes"))

# =============================================================================
# 3. PCA + Elbow Plot
# =============================================================================

message("Running PCA...")
n_pcs       <- 50
seurat_obj  <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)

elbow_p <- ElbowPlot(seurat_obj, ndims = n_pcs) +
    ggtitle(glue("PCA Elbow Plot: {umap_title}")) +
    theme_minimal()

ggsave(elbow_plot, plot = elbow_p, width = 8, height = 6)
message(glue("Elbow plot saved to: {elbow_plot}"))

# =============================================================================
# 4. Neighbors + Clustering (Louvain)
# =============================================================================

message(glue("Finding neighbors using top {pca_dims} PCs..."))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pca_dims)

message(glue("Clustering at resolution {cluster_res}..."))
seurat_obj <- FindClusters(seurat_obj, resolution = cluster_res, random.seed = SEED)

# =============================================================================
# 5. UMAP
# =============================================================================

message("Running UMAP...")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:pca_dims, seed.use = SEED)

umap_p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
    ggtitle(glue("UMAP: {umap_title}")) +
    theme_minimal()

ggsave(output_plot, plot = umap_p, width = 8, height = 6)
message(glue("UMAP plot saved to: {output_plot}"))

# =============================================================================
# 6. Cluster Summary
# =============================================================================

message("Cells per cluster:")
table(Idents(seurat_obj)) |>
    enframe(name = "cluster", value = "n_cells") |>
    mutate(n_cells = as.integer(n_cells)) |>
    arrange(desc(n_cells)) |>
    print()

# =============================================================================
# 7. Save
# =============================================================================

saveRDS(seurat_obj, file = output_rds)
message(glue("Clustered object saved to: {output_rds}"))
message("Step 06 complete.")