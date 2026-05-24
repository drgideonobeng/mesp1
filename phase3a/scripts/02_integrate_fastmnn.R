#!/usr/bin/env Rscript

# =============================================================================
# Step 02: Multi-Sample Integration (fastMNN via batchelor)
# Input : output_rds    — path for the integrated object
#         cpus          — worker count
#         n_features    — number of shared variable features for fastMNN
#         <rds1 ...>    — one or more 01_*_normalized.rds files (variable-length)
# Output: 02_seurat_integrated.rds — merged Seurat object with "mnn" DimReduc
#         (50 corrected dimensions); sample_id column retained in metadata
# Note  : Fixed args (output, cpus, n_features) come FIRST so variable-length
#         input file list can safely occupy all remaining positional arguments.
#         Uses RNA assay (LogNormalize + HVGs from step 01).
# =============================================================================

library(Seurat)
library(batchelor)
library(BiocParallel)
library(tidyverse)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: Rscript 02_integrate_fastmnn.R <output_rds> <cpus> <n_features> <mnn_k> <rds1> [rds2 ...]")
}

output_rds  <- args[1]
cpus        <- as.numeric(args[2])
n_features  <- as.numeric(args[3])
mnn_k       <- as.numeric(args[4])
input_files <- args[-(1:4)]

MNN_DIMS <- 50   # internal fastMNN PCA dims; downstream uses params$mnn_dims (≤50)
SEED     <- 42

message(glue("
========== Multi-Sample Integration (Step 02) ==========
  Samples    : {length(input_files)}
  Features   : {n_features}
  MNN dims   : {MNN_DIMS}
  k          : {mnn_k}
  Output     : {output_rds}
  Cores      : {cpus}
  Seed       : {SEED}
=========================================================
"))

# =============================================================================
# 2. Load Objects and Extract Sample IDs
# =============================================================================

# Filename convention from step 01: "01_<sample_id>_normalized.rds"
sample_ids <- basename(input_files) |>
    str_remove("^01_") |>
    str_remove("_normalized\\.rds$")

message(glue("Loading samples: {paste(sample_ids, collapse = ', ')}"))

objs <- map(input_files, readRDS)
names(objs) <- sample_ids

walk2(objs, sample_ids, \(obj, sid) {
    message(glue("  {sid}: {ncol(obj)} cells x {nrow(obj)} genes"))
})

# =============================================================================
# 3. Prefix Cell Barcodes with Sample ID
# =============================================================================

# Rename before merging so barcodes are unambiguous and consistent with
# the fastMNN output row order throughout this script.
objs <- imap(objs, \(obj, sid) {
    RenameCells(obj, new.names = paste0(sid, "_", colnames(obj)))
})

# =============================================================================
# 4. Select Shared Variable Features
# =============================================================================

message(glue("Selecting {n_features} integration features across all samples..."))
features <- SelectIntegrationFeatures(object.list = objs, nfeatures = n_features)

# Keep only features present in every sample's RNA data matrix
message("Extracting RNA data matrices...")
batches <- map(objs, \(obj) {
    GetAssayData(obj, assay = "RNA", layer = "data")
})

common_features <- features[features %in% Reduce(intersect, map(batches, rownames))]
message(glue(
    "Features retained after intersection across all batches: ",
    "{length(common_features)} of {length(features)}"
))

batches <- map(batches, \(mat) mat[common_features, ])

# =============================================================================
# 5. fastMNN Batch Correction
# =============================================================================

message("Running fastMNN batch correction...")

bpparam <- if (cpus > 1) {
    MulticoreParam(workers = cpus, RNGseed = SEED)
} else {
    SerialParam()
}

set.seed(SEED)
mnn_out <- do.call(
    fastMNN,
    c(
        unname(batches),
        list(d = MNN_DIMS, k = mnn_k, BPPARAM = bpparam)
    )
)

mnn_embeddings <- reducedDim(mnn_out, "corrected")
message(glue(
    "fastMNN complete: {nrow(mnn_embeddings)} cells x ",
    "{ncol(mnn_embeddings)} corrected dims"
))

# =============================================================================
# 5b. Elbow Plot — SD per corrected MNN dimension
# Use to validate the mnn_dims cutoff used downstream in step 03.
# fastMNN does not return percent-variance-explained; per-dim SD is the standard proxy.
# =============================================================================

dim_sd   <- apply(mnn_embeddings, 2, sd)
elbow_df <- data.frame(dim = seq_along(dim_sd), sd = dim_sd)

p_elbow <- ggplot(elbow_df, aes(x = dim, y = sd)) +
    geom_line(colour = "steelblue") +
    geom_point(size = 1.5, colour = "steelblue") +
    labs(
        title    = "fastMNN: SD per corrected dimension",
        subtitle = glue("Dashed line = mnn_dims cutoff used downstream"),
        x        = "MNN dimension",
        y        = "Std dev across cells"
    ) +
    theme_minimal()

ggsave("02_mnn_elbow.pdf", plot = p_elbow, width = 7, height = 4)
message("MNN elbow plot saved to: 02_mnn_elbow.pdf")

# =============================================================================
# 6. Merge Seurat Objects
# =============================================================================

message("Merging Seurat objects...")
merged <- merge(objs[[1]], y = objs[-1])

# Confirm cell order matches between fastMNN output and merged object.
# Both are produced from the same prefixed-barcode objects in the same order.
stopifnot(
    "Cell order mismatch between fastMNN output and merged Seurat object" =
        all(rownames(mnn_embeddings) == colnames(merged))
)

# =============================================================================
# 7. Store MNN Embedding as a Seurat DimReduc
# =============================================================================

merged[["mnn"]] <- CreateDimReducObject(
    embeddings = mnn_embeddings,
    key        = "mnn_",
    assay      = "RNA"
)

# =============================================================================
# 8. Summary Diagnostics
# =============================================================================

message(glue(
    "Integration summary: {ncol(merged)} total cells | ",
    "{length(unique(merged$sample_id))} samples"
))

message("Cells per sample:")
merged@meta.data |>
    as_tibble() |>
    count(sample_id, name = "n_cells") |>
    arrange(sample_id) |>
    print()

# =============================================================================
# 9. Save
# =============================================================================

saveRDS(merged, file = output_rds)
message(glue("Integrated object saved to: {output_rds}"))
message("Step 02 complete.")
