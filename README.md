# Mesp1 E9.5 Single-Cell RNA-seq Pipeline

![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A526.04.0-brightgreen)
![R](https://img.shields.io/badge/R-%E2%89%A54.3.0-blue)
![Seurat](https://img.shields.io/badge/Seurat-5.x-orange)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

A reproducible Nextflow pipeline for single-cell RNA-seq analysis of **Mesp1-lineage cells** 
isolated from E9.5 mouse cardiopharyngeal region. The pipeline processes raw 10x Genomics count matrices from GEO accession [GSM4816085](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4816085), 
performs quality control, normalization, clustering, marker gene discovery, and cluster annotation 
— producing publication-ready outputs at each stage.

---

## Table of Contents

- [Background](#background)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Overview](#pipeline-overview)
- [Output Structure](#output-structure)
- [Parameters](#parameters)
- [Execution Profiles](#execution-profiles)
- [Cluster Annotations](#cluster-annotations)
- [Citation](#citation)

---

## Background

**Mesp1** is a transcription factor and early mesodermal marker expressed transiently at E6.5–E7.5 
in mice. Mesp1-expressing progenitors give rise to the majority of cardiovascular cell types 
including cardiomyocytes, endothelial cells, and smooth muscle cells, as well as hematopoietic 
progenitors. This pipeline characterizes the transcriptional diversity of Mesp1-lineage cells 
at E9.5 — a critical stage of heart morphogenesis.

---

## Requirements

| Tool      | Version       | Purpose                        |
|-----------|---------------|--------------------------------|
| Nextflow  | >= 26.04.0    | Pipeline orchestration         |
| conda     | >= 23.0       | Environment management         |
| R         | >= 4.3.0      | Analysis runtime               |
| Seurat    | 5.x           | Single-cell analysis framework |

**Hardware:**

| Profile          | Minimum RAM | Recommended         |
|------------------|-------------|---------------------|
| `low_mem`        | 8 GB        | Standard laptop/PC  |
| `apple_silicon`  | 32 GB       | M-series Mac        |
| `hpc`            | 64 GB       | SLURM cluster       |

**Disk space:** ~50 GB (raw data + work directory + results)

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/drgideonobeng/bioinformatics/scRNAseq/mesp1-E95.git
cd mesp1-E95

# 2. Create and activate the conda environment
conda env create -f env.yml
conda activate multiome

# 3. Verify Nextflow is available
nextflow -version
```

---

## Usage

### Basic run

```bash
# Apple Silicon Mac
nextflow run main.nf -profile conda,apple_silicon

# Standard laptop (8–16 GB RAM)
nextflow run main.nf -profile conda,low_mem

# HPC cluster (SLURM)
nextflow run main.nf -profile conda,hpc
```

### Resume after a failure

```bash
nextflow run main.nf -profile conda,apple_silicon -resume
```

### Override parameters at runtime

```bash
nextflow run main.nf -profile conda,apple_silicon \
    --cluster_res 0.8 \
    --max_mt_percent 10 \
    --marker_genes "Mesp1,Hand2,Nkx2-5,Tbx5"
```

---

## Pipeline Overview

```
Raw GEO Data (GSM4816085)
        │
        ▼
01  DOWNLOAD_DATA              Download barcodes, features, matrix from GEO
        │
        ▼
02  CREATE_SEURAT_OBJECT       Build Seurat object from 10x sparse matrix
        │
        ▼
03  VISUALIZE_QC               Compute QC metrics; generate violin and scatter plots
        │
        ▼
04  FILTER_CELLS               Remove low-quality cells based on QC thresholds
        │
        ▼
05  NORMALIZE_DATA             SCTransform normalization; identify variable features
        │
        ▼
06  DIM_REDUCTION_AND_CLUSTER  PCA → UMAP; Louvain clustering
        │
        ├──────────────────────────────────┐
        ▼                                  ▼
07  RUN_CLUSTREE               08  FIND_MARKERS
    Evaluate resolution             Wilcoxon rank-sum test
    stability                       across all clusters
                                           │
                                           ▼
                               09  PLOT_MARKER_GENES
                                    Feature plots for
                                    selected genes
                                           │
                                           ▼
                               10  ANNOTATE_CLUSTERS
                                    Assign biological
                                    identities from
                                    cluster_labels.csv
```

| Step | Process                   | Output                                |
|------|---------------------------|---------------------------------------|
| 01   | DOWNLOAD_DATA             | Raw `.tsv.gz` and `.mtx.gz` files     |
| 02   | CREATE_SEURAT_OBJECT      | `02_seurat_unfiltered.rds`            |
| 03   | VISUALIZE_QC              | `03_qc_plots.pdf`                     |
| 04   | FILTER_CELLS              | `04_seurat_filtered.rds`              |
| 05   | NORMALIZE_DATA            | `05_seurat_normalized.rds`            |
| 06   | DIM_REDUCTION_AND_CLUSTER | `06_seurat_clustered.rds`, UMAP PDF   |
| 07   | RUN_CLUSTREE              | `07_clustree_resolutions.pdf`         |
| 08   | FIND_MARKERS              | `08_all_markers.csv`, top 5 CSV       |
| 09   | PLOT_MARKER_GENES         | `09_marker_feature_plots.pdf`         |
| 10   | ANNOTATE_CLUSTERS         | `10_seurat_annotated.rds`, UMAP PDF   |

---

## Output Structure

```
results/
├── objects/                        # Seurat RDS objects at each stage
│   ├── 02_seurat_unfiltered.rds
│   ├── 04_seurat_filtered.rds
│   ├── 05_seurat_normalized.rds
│   ├── 06_seurat_clustered.rds
│   └── 10_seurat_annotated.rds
│
├── plots/                          # All visualizations
│   ├── 03_qc_plots.pdf
│   ├── 06_elbow.pdf
│   ├── 06_umap.pdf
│   ├── 07_clustree_resolutions.pdf
│   ├── 09_marker_feature_plots.pdf
│   └── 10_umap_annotated.pdf
│
├── tables/                         # Marker gene tables
│   ├── 08_all_markers.csv
│   └── 08_top5_markers_per_cluster.csv
│
└── pipeline_info/                  # Execution diagnostics
    ├── report.html
    ├── timeline.html
    ├── trace.txt
    └── flowchart.html
```

---

## Parameters

All parameters can be set in `nextflow.config` or overridden at runtime with `--param value`.

### QC Thresholds

| Parameter       | Default | Description                                        |
|-----------------|---------|----------------------------------------------------|
| `min_cells`     | 3       | Minimum cells a gene must appear in                |
| `min_genes`     | 200     | Minimum genes per cell (removes empty droplets)    |
| `max_genes`     | 6000    | Maximum genes per cell (removes doublets)          |
| `max_mt_percent`| 3       | Maximum mitochondrial % (removes dying cells)      |

### Clustering

| Parameter     | Default | Description                                             |
|---------------|---------|---------------------------------------------------------|
| `pca_dims`    | 20      | PCs used for UMAP and clustering (validate with elbow plot) |
| `cluster_res` | 0.5     | Louvain resolution (validate with clustree output)      |

### Biological

| Parameter        | Default                        | Description                        |
|------------------|--------------------------------|------------------------------------|
| `species`        | `mouse`                        | Used for gene name formatting      |
| `mt_pattern`     | `^mt-`                         | Mitochondrial gene prefix (`^MT-` for human) |
| `marker_genes`   | `Mesp1,Prdm1,Hand2,...`        | Genes visualized in Step 09        |
| `cluster_labels` | `cluster_labels.csv`           | Path to cluster annotation file    |

---

## Execution Profiles

Profiles are combined with the `-profile` flag. Always pair a compute profile with `conda`:

```bash
-profile conda,apple_silicon
```

| Profile          | Use case                                      |
|------------------|-----------------------------------------------|
| `standard`       | Uses your active terminal environment         |
| `conda`          | Reproducible run via `env.yml` (recommended)  |
| `apple_silicon`  | M1/M2/M3 Mac with 32 GB+ unified memory       |
| `low_mem`        | Standard laptop or desktop (8–16 GB RAM)      |
| `hpc`            | SLURM cluster with high-memory queue          |
| `docker`         | Containerized run via Seurat Docker image     |

---

## Cluster Annotations

After reviewing marker gene outputs (`08_all_markers.csv`, `09_marker_feature_plots.pdf`), 
update `cluster_labels.csv` in the project root to assign biological identities:

```csv
cluster_id,label
0,Lateral plate mesoderm
1,Somitic mesoderm
2,Proepicardium
3,Skeletal muscle progenitors
4,Surface ectoderm
5,Progenitor
6,Cardiomyocytes
7,Endothelial cells
8,Pharyngeal arch
9,SHF progenitors
10,Hematopoietic cells
```

Then resume the pipeline — only `ANNOTATE_CLUSTERS` will re-run:

```bash
nextflow run main.nf -profile conda,apple_silicon -resume
```

---

## Citation

If you use this pipeline in your research, please cite the following:

**Dataset:**
> Nomaru, Hiroko et al. (2021). Title. *Journal*. DOI.
> GEO Accession: [GSM4816085](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4816085)

**Seurat / SCTransform:**
> Hao et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573–3587.
> https://doi.org/10.1016/j.cell.2021.04.048

> Choudhary & Satija (2022). Comparison and evaluation of statistical error models for scRNA-seq. 
> *Genome Biology*, 23, 27. https://doi.org/10.1186/s13059-021-02584-9

**Nextflow:**
> Di Tommaso et al. (2017). Nextflow enables reproducible computational workflows. 
> *Nature Biotechnology*, 35, 316–319. https://doi.org/10.1038/nbt.3820

---

*Pipeline developed as part of a single-cell bioinformatics project.*
