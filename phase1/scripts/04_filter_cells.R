#!/usr/bin/env Rscript

# =============================================================================
# Step 04: Filter Cells
# Removes low-quality cells based on gene count and mitochondrial thresholds.
# =============================================================================

library(Seurat)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args           <- commandArgs(trailingOnly = TRUE)
input_rds      <- args[1]
min_genes      <- as.numeric(args[2])
max_genes      <- as.numeric(args[3])
max_mt_percent <- as.numeric(args[4])
output_rds     <- args[5]

message(glue("
========== Cell Filtering (Step 04) ==========
  Input          : {input_rds}
  Output         : {output_rds}
  min_genes      : {min_genes}
  max_genes      : {max_genes}
  max_mt_percent : {max_mt_percent}%
==============================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

seurat_obj      <- readRDS(input_rds)
pre_filter_n   <- ncol(seurat_obj)
message(glue("Loaded {pre_filter_n} cells from: {input_rds}"))

# =============================================================================
# 3. Apply Filters and Tally What Each Criterion Removes
# =============================================================================

meta <- seurat_obj@meta.data

filter_summary <- tibble(
    Criterion           = c("nFeature_RNA <= min", "nFeature_RNA >= max", "percent.mt >= max"),
    Threshold           = c(min_genes, max_genes, max_mt_percent),
    Cells_Failing       = c(
        sum(meta$nFeature_RNA <= min_genes),
        sum(meta$nFeature_RNA >= max_genes),
        sum(meta$percent.mt   >= max_mt_percent)
    )
)

message("Cells failing each criterion (overlap possible):")
print(filter_summary)

seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > min_genes &
             nFeature_RNA < max_genes &
             percent.mt   < max_mt_percent
)

# =============================================================================
# 4. Report and Save
# =============================================================================

post_filter_n <- ncol(seurat_obj)
removed_n     <- pre_filter_n - post_filter_n
pct_retained  <- round(100 * post_filter_n / pre_filter_n, 1)

message(glue(
    "Retained {post_filter_n} of {pre_filter_n} cells ",
    "({pct_retained}%) | Removed: {removed_n}"
))

saveRDS(seurat_obj, file = output_rds)
message(glue("Filtered object saved to: {output_rds}"))
message("Step 04 complete.")