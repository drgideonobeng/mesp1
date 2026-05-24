#!/usr/bin/env Rscript

# =============================================================================
# Step 03: Dimensionality Reduction and Clustering (integrated atlas)
# Input : 02_seurat_integrated.rds — fastMNN-corrected Seurat object
# Output: 03_seurat_clustered.rds  — UMAP-embedded, Louvain-clustered object
#         03_umap_clusters.pdf     — UMAP coloured by Louvain cluster
#         03_umap_samples.pdf      — UMAP coloured by timepoint (sample_id)
# Note  : FindNeighbors uses the "mnn" reduction so cluster structure reflects
#         batch-corrected space, not raw gene expression.
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop(paste(
        "Usage: Rscript 03_cluster_integrated.R",
        "<input_rds> <mnn_dims> <cluster_res> <project_name>",
        "<output_rds> <umap_clusters_pdf> <umap_samples_pdf>"
    ))
}

input_rds        <- args[1]
mnn_dims         <- as.numeric(args[2])
cluster_res      <- as.numeric(args[3])
project_name     <- args[4]
output_rds       <- args[5]
umap_cluster_pdf <- args[6]
umap_sample_pdf  <- args[7]

SEED <- 42

message(glue("
========== Dim Reduction & Clustering (Step 03) ==========
  Input       : {input_rds}
  mnn_dims    : {mnn_dims}
  cluster_res : {cluster_res}
  Seed        : {SEED}
==========================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

seurat_obj <- readRDS(input_rds)
message(glue(
    "Loaded: {ncol(seurat_obj)} cells | ",
    "{length(unique(seurat_obj$sample_id))} samples"
))

# =============================================================================
# 3. Graph Construction on MNN Embedding
# =============================================================================

# Set RNA as default assay so graph and cluster column names use "RNA_snn_res."
# prefix — compatible with clustree (step 04).
DefaultAssay(seurat_obj) <- "RNA"

message(glue("Building neighbor graph on MNN dims 1:{mnn_dims}..."))
seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = "mnn",
    dims      = 1:mnn_dims,
    verbose   = FALSE
)

message(glue("Clustering at resolution {cluster_res}..."))
seurat_obj <- FindClusters(
    seurat_obj,
    resolution  = cluster_res,
    random.seed = SEED
)

# =============================================================================
# 4. UMAP
# =============================================================================

message("Running UMAP on MNN embedding...")
seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "mnn",
    dims      = 1:mnn_dims,
    seed.use  = SEED
)

# =============================================================================
# 5. Cluster UMAP
# =============================================================================

n_clusters <- length(levels(seurat_obj))

p_cluster <- DimPlot(
    seurat_obj,
    reduction = "umap",
    label     = TRUE,
    repel     = TRUE,
    pt.size   = 0.3
) +
    labs(
        title    = glue("{project_name}: Louvain Clusters"),
        subtitle = glue("{ncol(seurat_obj)} cells | {n_clusters} clusters")
    ) +
    theme_minimal() +
    theme(legend.position = "right")

ggsave(umap_cluster_pdf, plot = p_cluster, width = 9, height = 7)
message(glue("Cluster UMAP saved to: {umap_cluster_pdf}"))

# =============================================================================
# 6. Sample / Timepoint UMAP
# =============================================================================

p_sample <- DimPlot(
    seurat_obj,
    reduction  = "umap",
    group.by   = "sample_id",
    pt.size    = 0.3,
    shuffle    = TRUE,
    seed       = SEED
) +
    labs(
        title    = glue("{project_name}: Sample / Timepoint"),
        subtitle = glue("{ncol(seurat_obj)} cells")
    ) +
    theme_minimal() +
    theme(legend.position = "right")

ggsave(umap_sample_pdf, plot = p_sample, width = 9, height = 7)
message(glue("Sample UMAP saved to: {umap_sample_pdf}"))

# =============================================================================
# 7. Cluster Summary
# =============================================================================

message("Cells per cluster:")
table(Idents(seurat_obj)) |>
    enframe(name = "cluster", value = "n_cells") |>
    mutate(n_cells = as.integer(n_cells)) |>
    arrange(desc(n_cells)) |>
    print()

composition <- seurat_obj@meta.data |>
    as_tibble() |>
    count(cluster = Idents(seurat_obj), sample_id) |>
    pivot_wider(names_from = sample_id, values_from = n, values_fill = 0L)

message("Timepoint composition per cluster:")
print(composition)

write_csv(composition, "03_timepoint_composition.csv")
message("Composition table saved to: 03_timepoint_composition.csv")

# =============================================================================
# 8. Save
# =============================================================================

saveRDS(seurat_obj, file = output_rds)
message(glue("Clustered object saved to: {output_rds}"))
message("Step 03 complete.")
