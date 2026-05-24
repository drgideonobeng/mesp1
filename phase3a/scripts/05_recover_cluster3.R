#!/usr/bin/env Rscript

# =============================================================================
# Step 05: Cluster-3 Contamination Recovery
# Input : 03_seurat_clustered.rds — clustered integrated object
# Output: 05_atlas_clean.rds        — atlas with cluster 3 replaced by
#                                     cleaned sub-clusters
#         05_recovery_summary.txt   — cells removed / recovered / final counts
#         05_umap_atlas_clean.pdf   — UMAP coloured by seurat_clusters
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript 05_recover_cluster3.R <input_rds> <out_rds> <out_summary> <out_umap>")
}

input_rds   <- args[1]
output_rds  <- args[2]
output_txt  <- args[3]
output_umap <- args[4]

SEED <- 42

message(glue("
========== Cluster-3 Recovery (Step 05) ==========
  Input   : {input_rds}
  RDS out : {output_rds}
  Summary : {output_txt}
  UMAP    : {output_umap}
  Seed    : {SEED}
===================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

seu <- readRDS(input_rds)
message(glue("Loaded: {ncol(seu)} cells, {nlevels(Idents(seu))} clusters"))

# =============================================================================
# 3. Cluster Sizes and Composition Table
# =============================================================================

message("\n--- Cluster sizes ---")
print(table(Idents(seu)))

message("\n--- Cluster × timepoint composition ---")
print(table(
    cluster   = as.character(Idents(seu)),
    timepoint = seu$orig.ident
))

# =============================================================================
# 4. Subset Cluster 3
# =============================================================================

cluster3_cells <- WhichCells(seu, idents = "3")
message(glue("\nCluster 3: {length(cluster3_cells)} cells — subsetting for re-processing"))

cluster3 <- subset(seu, cells = cluster3_cells)

# =============================================================================
# 5. Re-process Cluster 3 in Isolation
# =============================================================================

message("Re-processing cluster 3 in isolation...")

DefaultAssay(cluster3) <- "RNA"
cluster3 <- JoinLayers(cluster3)

set.seed(SEED); cluster3 <- NormalizeData(cluster3, verbose = FALSE)
set.seed(SEED); cluster3 <- FindVariableFeatures(cluster3, nfeatures = 2000, verbose = FALSE)
set.seed(SEED); cluster3 <- ScaleData(cluster3, verbose = FALSE)
set.seed(SEED); cluster3 <- RunPCA(cluster3, npcs = 20, seed.use = SEED, verbose = FALSE)
set.seed(SEED); cluster3 <- FindNeighbors(cluster3, dims = 1:15, verbose = FALSE)

# Resolution sweep; use 0.4 for final sub-clustering
set.seed(SEED)
cluster3 <- FindClusters(
    cluster3,
    resolution  = c(0.2, 0.4, 0.6, 0.8),
    random.seed = SEED,
    verbose     = FALSE
)
Idents(cluster3) <- "RNA_snn_res.0.4"

set.seed(SEED); cluster3 <- RunUMAP(cluster3, dims = 1:15, seed.use = SEED, verbose = FALSE)
cluster3 <- JoinLayers(cluster3)

message(glue("Sub-clusters at res 0.4: {nlevels(Idents(cluster3))}")))
print(table(sub_cluster = Idents(cluster3)))

# =============================================================================
# 6. Find Sub-cluster Markers
# =============================================================================

message("Finding sub-cluster markers...")
set.seed(SEED)
sub_markers <- FindAllMarkers(
    cluster3,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    verbose         = FALSE
)

message("\n--- Top 3 markers per sub-cluster ---")
sub_markers |>
    group_by(cluster) |>
    slice_max(order_by = avg_log2FC, n = 3) |>
    ungroup() |>
    print()

# =============================================================================
# 7. Identify Contamination Sub-cluster
# =============================================================================

contam_genes <- c("Lbx1", "Erbb3")

top5 <- sub_markers |>
    group_by(cluster) |>
    slice_max(order_by = avg_log2FC, n = 5) |>
    ungroup()

contam_subclusters <- top5 |>
    filter(gene %in% contam_genes) |>
    pull(cluster) |>
    unique() |>
    as.character()

if (length(contam_subclusters) == 0) {
    warning(
        "No sub-cluster has Lbx1 or Erbb3 in its top 5 markers. ",
        "Proceeding without exclusion — review sub_markers manually."
    )
    contam_cells <- character(0)
} else {
    hit_genes <- intersect(
        top5$gene[top5$cluster %in% contam_subclusters],
        contam_genes
    )
    message(glue(
        "Contamination sub-cluster(s): {paste(contam_subclusters, collapse = ', ')} ",
        "(markers: {paste(hit_genes, collapse = ', ')})"
    ))
    contam_cells <- WhichCells(cluster3, idents = contam_subclusters)
    message(glue("Cells to remove: {length(contam_cells)}"))
}

# =============================================================================
# 8. Build atlas_clean
# =============================================================================

message("Building atlas_clean...")

kept_cells    <- setdiff(colnames(cluster3), contam_cells)
non_cl3_cells <- setdiff(colnames(seu), cluster3_cells)

# Subset original object — preserves MNN embedding and UMAP from step 03
atlas_clean <- subset(seu, cells = c(non_cl3_cells, kept_cells))

# Base labels from seurat_clusters; overwrite recovered cl3 cells
atlas_cluster_labels        <- as.character(atlas_clean$seurat_clusters)
names(atlas_cluster_labels) <- colnames(atlas_clean)

kept_sub_ids                     <- as.character(Idents(cluster3)[kept_cells])
atlas_cluster_labels[kept_cells] <- paste0("cl3_sub_", kept_sub_ids)

atlas_clean$atlas_cluster <- atlas_cluster_labels
Idents(atlas_clean) <- "atlas_cluster"

message(glue("atlas_clean assembled: {ncol(atlas_clean)} cells"))

# =============================================================================
# 9. Recovery Summary
# =============================================================================

n_original_cl3   <- length(cluster3_cells)
n_removed        <- length(contam_cells)
n_recovered      <- length(kept_cells)
n_final_cells    <- ncol(atlas_clean)
n_final_clusters <- length(unique(atlas_clean$atlas_cluster))

summary_lines <- c(
    "========== Cluster-3 Recovery Summary ==========",
    glue("  Original cluster 3 cells     : {n_original_cl3}"),
    glue("  Contamination cells removed  : {n_removed}"),
    glue("  Cluster 3 cells recovered    : {n_recovered}"),
    glue("  Final atlas cell count       : {n_final_cells}"),
    glue("  Final atlas cluster count    : {n_final_clusters}"),
    "================================================="
)

message(paste(summary_lines, collapse = "\n"))
writeLines(summary_lines, output_txt)

# =============================================================================
# 10. UMAP
# =============================================================================

message("Producing UMAP of atlas_clean...")

umap_plot <- DimPlot(
    atlas_clean,
    reduction = "umap",
    group.by  = "seurat_clusters",
    label     = TRUE,
    repel     = TRUE
) +
    ggtitle("Atlas Clean: Clusters (cluster 3 contamination removed)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

umap_plot2 <- DimPlot(
    atlas_clean,
    reduction = "umap",
    group.by  = "atlas_cluster",
    label     = TRUE,
    repel     = TRUE
) +
    ggtitle("Atlas Clean: Atlas Clusters (cl3 sub-clusters labeled)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(output_umap, width = 10, height = 8)
print(umap_plot)
print(umap_plot2)
dev.off()

message(glue("UMAP saved to: {output_umap}"))

# =============================================================================
# 11. Save
# =============================================================================

saveRDS(atlas_clean, output_rds)
message(glue("Saved: {output_rds}"))
message("Step 05 complete.")
