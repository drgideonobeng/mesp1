#!/usr/bin/env Rscript

# =============================================================================
# Step 07: Annotate Atlas
# Applies biological labels to jointly-clustered integrated object, then finds
# marker genes per annotated cell type to document and validate identities.
# Input : 05_atlas_clean.rds           — cleaned atlas (cluster 3 recovered)
#         project_name                 — title string for plots
#         cluster_labels_atlas.csv     — two-column CSV: cluster_id, label
#         cpus                         — worker count
# Output: 07_seurat_annotated.rds      — annotated Seurat object
#         07_umap_annotated.pdf        — UMAP coloured by cell type
#         07_umap_timepoints.pdf       — UMAP faceted by timepoint
#         07_all_markers.csv           — full Wilcoxon marker table
#         07_top_5_markers_per_celltype.csv
# Note  : Re-run only this step after updating cluster_labels_atlas.csv
#         (nextflow run main.nf ... -resume).
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)
library(future)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
    stop(paste(
        "Usage: Rscript 05_annotate_atlas.R",
        "<input_rds> <project_name> <labels_csv> <cpus>",
        "<output_rds> <umap_annotated_pdf> <umap_timepoints_pdf>",
        "<all_markers_csv> <top_markers_csv>"
    ))
}

input_rds         <- args[1]
project_name      <- args[2]
labels_file       <- args[3]
cpus              <- as.numeric(args[4])
output_rds        <- args[5]
umap_annotated    <- args[6]
umap_timepoints   <- args[7]
all_markers_csv   <- args[8]
top_markers_csv   <- args[9]

TOP_N <- 5
SEED  <- 42

message(glue("
========== Atlas Annotation (Step 05) ==========
  Input   : {input_rds}
  Labels  : {labels_file}
  Output  : {output_rds}
  Cores   : {cpus}
  Seed    : {SEED}
================================================
"))

# =============================================================================
# 2. Configure Parallelism
# =============================================================================

plan(multisession, workers = cpus)
options(future.globals.maxSize = 10 * 1024^3)   # 10 GB

# =============================================================================
# 3. Load Labels
# =============================================================================

labels_df <- read_csv(labels_file, col_types = cols(.default = "c"))

if (any(labels_df$label == "" | is.na(labels_df$label))) {
    stop(glue(
        "cluster_labels_atlas.csv has empty labels. ",
        "Fill in all labels before running step 05."
    ))
}

cluster_identities <- labels_df |> deframe()
message(glue("Loaded {nrow(labels_df)} cluster labels from: {labels_file}"))

# =============================================================================
# 4. Load Object
# =============================================================================

seurat_obj <- readRDS(input_rds)
message(glue(
    "Loaded: {ncol(seurat_obj)} cells | ",
    "{length(levels(seurat_obj))} clusters | ",
    "{length(unique(seurat_obj$sample_id))} samples"
))

# =============================================================================
# 5. Safety Checks
# =============================================================================

data_clusters     <- levels(seurat_obj)
labeled_ids       <- names(cluster_identities)
missing_from_data <- setdiff(labeled_ids, data_clusters)
unlabeled         <- setdiff(data_clusters, labeled_ids)

if (length(missing_from_data) > 0) {
    warning(glue(
        "Labels in CSV but not in data: ",
        "{paste(missing_from_data, collapse = ', ')}"
    ))
}

if (length(unlabeled) > 0) {
    warning(glue(
        "Clusters in data but missing from CSV (numeric IDs kept): ",
        "{paste(unlabeled, collapse = ', ')}"
    ))
}

# =============================================================================
# 6. Apply Annotations
# =============================================================================

seurat_obj <- RenameIdents(seurat_obj, cluster_identities)
seurat_obj$cell_type <- Idents(seurat_obj)

message("Cell counts per annotated type:")
seurat_obj$cell_type |>
    table() |>
    enframe(name = "cell_type", value = "n_cells") |>
    mutate(n_cells = as.integer(n_cells)) |>
    arrange(desc(n_cells)) |>
    print()

# =============================================================================
# 7. UMAP: Annotated Cell Types
# =============================================================================

p_annotated <- DimPlot(
    seurat_obj,
    reduction = "umap",
    label     = TRUE,
    repel     = TRUE,
    pt.size   = 0.3
) +
    labs(
        title    = glue("{project_name}: Annotated Cell Types"),
        subtitle = glue("{ncol(seurat_obj)} cells | {length(levels(seurat_obj))} cell types")
    ) +
    theme_minimal() +
    theme(legend.position = "right")

ggsave(umap_annotated, plot = p_annotated, width = 11, height = 8)
message(glue("Annotated UMAP saved to: {umap_annotated}"))

# =============================================================================
# 8. UMAP: Faceted by Timepoint
# =============================================================================

p_timepoints <- DimPlot(
    seurat_obj,
    reduction = "umap",
    split.by  = "sample_id",
    label     = FALSE,
    pt.size   = 0.2,
    ncol      = 2
) +
    labs(title = glue("{project_name}: Cell Types by Timepoint")) +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(umap_timepoints, plot = p_timepoints, width = 12, height = 10)
message(glue("Timepoint UMAP saved to: {umap_timepoints}"))

# =============================================================================
# 9. Find Markers (validates annotations)
# =============================================================================

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- JoinLayers(seurat_obj)

message("Running FindAllMarkers (Wilcoxon, positive markers only)...")
set.seed(SEED)

markers <- FindAllMarkers(
    object          = seurat_obj,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    verbose         = FALSE
)

per_type_counts <- markers |>
    as_tibble() |>
    count(cluster, name = "n_markers")

message(glue(
    "Total markers: {nrow(markers)} | ",
    "Median per cell type: {median(per_type_counts$n_markers)}"
))

write_csv(markers, all_markers_csv)
message(glue("All markers saved to: {all_markers_csv}"))

# =============================================================================
# 10. Top-N Cheat Sheet
# =============================================================================

top_genes <- markers |>
    as_tibble() |>
    group_by(cluster) |>
    slice_max(n = TOP_N, order_by = avg_log2FC, with_ties = FALSE) |>
    ungroup()

write_csv(top_genes, top_markers_csv)
message(glue("Top-{TOP_N} markers saved to: {top_markers_csv}"))

message(glue("\n{str_dup('-', 50)}"))
message(glue("SIGNATURE MARKERS (TOP {TOP_N} PER CELL TYPE)"))
message(str_dup("-", 50))
top_genes |> select(cluster, gene, avg_log2FC, p_val_adj) |> print(n = Inf)

# =============================================================================
# 11. Save Annotated Object
# =============================================================================

saveRDS(seurat_obj, file = output_rds)
message(glue("Annotated object saved to: {output_rds}"))
message("Step 05 complete.")
