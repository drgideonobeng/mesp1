#!/usr/bin/env Rscript

# =============================================================================
# Step 05: Normalize Data with SCTransform
# Performs regularized negative binomial normalization with mt% regression.
# =============================================================================

library(Seurat)
library(glue)
library(future)

# =============================================================================
# 1. Arguments
# =============================================================================

args        <- commandArgs(trailingOnly = TRUE)
input_rds   <- args[1]
output_rds  <- args[2]
cpus_to_use <- if (length(args) >= 3) as.numeric(args[3]) else 1

message(glue("
========== SCTransform Normalization (Step 05) ==========
  Input  : {input_rds}
  Output : {output_rds}
  Cores  : {cpus_to_use}
=========================================================
"))

# =============================================================================
# 2. Configure Parallelism
# =============================================================================

# multisession works reliably on macOS (multicore is fragile with fork())
plan(multisession, workers = cpus_to_use)

# Raise the per-worker object size cap (SCTransform sends the count matrix to workers)
options(future.globals.maxSize = 8 * 1024^3)   # 8 GB

# =============================================================================
# 3. Load Data
# =============================================================================

message(glue("Loading filtered object from: {input_rds}"))
seurat_obj <- readRDS(input_rds)
message(glue("Object loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes"))

# =============================================================================
# 4. SCTransform
# =============================================================================

message("Running SCTransform with percent.mt regression...")
seurat_obj <- SCTransform(
    object          = seurat_obj,
    vars.to.regress = "percent.mt",
    verbose         = TRUE
)

# =============================================================================
# 5. Summary
# =============================================================================

n_variable <- length(VariableFeatures(seurat_obj))
top10      <- head(VariableFeatures(seurat_obj), 10)

message(glue("Identified {n_variable} variable features."))
message(glue("Top 10: {paste(top10, collapse = ', ')}"))

# =============================================================================
# 6. Save
# =============================================================================

saveRDS(seurat_obj, file = output_rds)
message(glue("Normalized object saved to: {output_rds}"))
message("Step 05 complete.")