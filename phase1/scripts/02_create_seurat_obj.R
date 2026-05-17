#!/usr/bin/env Rscript

# =============================================================================
# Step 02: Create Seurat Object
# Reads 10x Genomics sparse matrix files and constructs a Seurat object
# with initial gene/cell filtering and mitochondrial percentage annotation.
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
    stop("Usage: Rscript 02_create_seurat_obj.R <barcodes> <features> <matrix> <min_cells> <min_features> <mt_pattern> <project_name> <output_rds>")
}

barcodes_path <- args[1]
features_path <- args[2]
matrix_path   <- args[3]
min_cells     <- as.numeric(args[4])
min_features  <- as.numeric(args[5])
mt_pattern    <- args[6]
project_name  <- args[7]
output_rds    <- args[8]

message(glue("
========== Create Seurat Object (Step 02) ==========
  Project        : {project_name}
  min_cells      : {min_cells}
  min_features   : {min_features}
  mt_pattern     : {mt_pattern}
  Output         : {output_rds}
====================================================
"))

# =============================================================================
# 2. Load 10x Genomics Matrix
# =============================================================================

# Nextflow stages all three input files in the current process directory,
# so Read10X reads from "." rather than from named paths.
message("Loading 10x sparse matrix from staged files...")
counts <- Read10X(data.dir = ".")
message(glue("Raw matrix: {nrow(counts)} genes x {ncol(counts)} cells"))

# =============================================================================
# 3. Create Seurat Object
# =============================================================================

message(glue("Building Seurat object with min.cells={min_cells}, min.features={min_features}..."))
seurat_obj <- CreateSeuratObject(
    counts       = counts,
    project      = project_name,
    min.cells    = min_cells,
    min.features = min_features
)

message(glue(
    "Filtered matrix: {nrow(seurat_obj)} genes x {ncol(seurat_obj)} cells ",
    "(removed {nrow(counts) - nrow(seurat_obj)} genes, {ncol(counts) - ncol(seurat_obj)} cells)"
))

# =============================================================================
# 4. Mitochondrial Percentage
# =============================================================================

message(glue("Calculating mitochondrial percentage using pattern: '{mt_pattern}'"))
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)

# Quick diagnostic summary of QC metrics
qc_summary <- tibble(
    metric  = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    median  = c(
        median(seurat_obj$nFeature_RNA),
        median(seurat_obj$nCount_RNA),
        round(median(seurat_obj$percent.mt), 2)
    ),
    min     = c(
        min(seurat_obj$nFeature_RNA),
        min(seurat_obj$nCount_RNA),
        round(min(seurat_obj$percent.mt), 2)
    ),
    max     = c(
        max(seurat_obj$nFeature_RNA),
        max(seurat_obj$nCount_RNA),
        round(max(seurat_obj$percent.mt), 2)
    )
)

message("QC metric distribution:")
print(qc_summary)

# =============================================================================
# 5. Save
# =============================================================================

saveRDS(seurat_obj, file = output_rds)
message(glue("Unfiltered Seurat object saved to: {output_rds}"))
message("Step 02 complete.")