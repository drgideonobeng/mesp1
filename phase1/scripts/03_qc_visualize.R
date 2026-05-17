#!/usr/bin/env Rscript

# =============================================================================
# Step 03: Visualize Quality Control Metrics
# Generates violin and scatter plots of QC metrics (nFeature, nCount, percent.mt)
# to inform appropriate filtering thresholds.
# =============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)
library(glue)

# =============================================================================
# 1. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript 03_qc_visualize.R <input_rds> <output_rds> <output_pdf>")
}

input_rds   <- args[1]
output_rds  <- args[2]
output_plot <- args[3]

message(glue("
========== QC Visualization (Step 03) ==========
  Input  : {input_rds}
  RDS    : {output_rds}
  Plot   : {output_plot}
================================================
"))

# =============================================================================
# 2. Load Object
# =============================================================================

seurat_obj <- readRDS(input_rds)
message(glue("Loaded: {ncol(seurat_obj)} cells x {nrow(seurat_obj)} genes"))

# =============================================================================
# 3. Generate QC Plots
# =============================================================================

message("Generating violin and scatter plots...")

# Violin plots — distribution per metric
violin_plot <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol     = 3,
    pt.size  = 0.1
) &
    theme_minimal() &
    theme(legend.position = "none")

# Scatter plots — detect technical/biological outliers
scatter_count_vs_feature <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
) +
    theme_minimal()

scatter_count_vs_mt <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
) +
    theme_minimal()

# Combine into a 2-row layout: violins on top, scatters below
combined_plot <- (violin_plot) /
                 (scatter_count_vs_feature | scatter_count_vs_mt) +
    plot_annotation(
        title = glue("QC Metrics: {seurat_obj@project.name}"),
        subtitle = glue("{ncol(seurat_obj)} cells before filtering"),
        theme = theme(plot.title = element_text(face = "bold"))
    )

# =============================================================================
# 4. Save Plot
# =============================================================================

ggsave(output_plot, plot = combined_plot, width = 12, height = 10)
message(glue("QC plots saved to: {output_plot}"))

# =============================================================================
# 5. Print Summary Stats to Terminal
# =============================================================================

qc_summary <- tibble(
    metric = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    median = c(
        median(seurat_obj$nFeature_RNA),
        median(seurat_obj$nCount_RNA),
        round(median(seurat_obj$percent.mt), 2)
    ),
    Q1     = c(
        quantile(seurat_obj$nFeature_RNA, 0.25),
        quantile(seurat_obj$nCount_RNA,   0.25),
        round(quantile(seurat_obj$percent.mt, 0.25), 2)
    ),
    Q3     = c(
        quantile(seurat_obj$nFeature_RNA, 0.75),
        quantile(seurat_obj$nCount_RNA,   0.75),
        round(quantile(seurat_obj$percent.mt, 0.75), 2)
    )
)

message("QC distribution summary (use to guide threshold choices):")
print(qc_summary)

# =============================================================================
# 6. Save Object (Unchanged — Passes Through to Next Step)
# =============================================================================

saveRDS(seurat_obj, file = output_rds)
message(glue("Seurat object saved to: {output_rds}"))
message("Step 03 complete.")