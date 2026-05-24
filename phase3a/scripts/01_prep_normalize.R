#!/usr/bin/env Rscript

# =============================================================================
# Step 01: Prep / Normalize (LogNormalize + FindVariableFeatures per sample)
# Input : 04_seurat_filtered.rds  — Phase 1 filtered Seurat object
#         sample_id               — timepoint label (e.g. "e95")
# Output: 01_<sample_id>_normalized.rds — log-normalized object with HVGs and
#         sample_id stored in metadata; ready for multi-sample integration
# =============================================================================

library(Seurat)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript 01_prep_normalize.R <input_rds> <sample_id> <output_rds> <n_features>")
}

input_rds  <- args[1]
sample_id  <- args[2]
output_rds <- args[3]
n_features <- as.numeric(args[4])
SEED <- 42

message(glue("
========== Per-Sample Normalization (Step 01) ==========
  Sample     : {sample_id}
  Input      : {input_rds}
  Output     : {output_rds}
  n_features : {n_features}
  Seed       : {SEED}
=========================================================
"))

# =============================================================================
# 3. Load Object
# =============================================================================

message(glue("Loading filtered object from: {input_rds}"))
seurat_obj <- readRDS(input_rds)
message(glue("Loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes"))

# =============================================================================
# 4. Tag Cells with Timepoint
# =============================================================================

seurat_obj$sample_id <- sample_id

# =============================================================================
# 5. Log-Normalize + Find Variable Features
# =============================================================================

message("Running NormalizeData (LogNormalize, scale.factor = 10000)...")
seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor         = 10000,
    verbose              = FALSE
)

message(glue("Finding variable features (vst, nfeatures = {n_features})..."))
seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures        = n_features,
    verbose          = FALSE
)

n_var <- length(VariableFeatures(seurat_obj))
top10 <- head(VariableFeatures(seurat_obj), 10)
message(glue("Variable features identified: {n_var}"))
message(glue("Top 10: {paste(top10, collapse = ', ')}"))

# =============================================================================
# 6. Save
# =============================================================================

saveRDS(seurat_obj, file = output_rds)
message(glue("Normalized object saved to: {output_rds}"))
message("Step 01 complete.")
