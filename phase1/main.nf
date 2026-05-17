nextflow.enable.dsl=2

// =============================================================================
// main.nf — Phase 1: Single-Cell RNA-seq Pipeline (Mesp1-cardiopharyngeal lineage)
// 7-step workflow: download → object creation → QC → filter → normalize →
// dim reduction → clustree
// =============================================================================

// --- PROCESS 1: Download raw data from GEO ---
process DOWNLOAD_DATA {
    publishDir "${params.datadir}", mode: 'copy'
    input:
        val bc_url
        val feat_url
        val mtx_url
    output:
        path "barcodes.tsv.gz", emit: barcodes
        path "features.tsv.gz", emit: features
        path "matrix.mtx.gz",   emit: matrix
    script:
        """
        wget -O barcodes.tsv.gz "${bc_url}"
        wget -O features.tsv.gz "${feat_url}"
        wget -O matrix.mtx.gz   "${mtx_url}"
        """
}

// --- PROCESS 2: Create Seurat Object ---
process CREATE_SEURAT_OBJECT {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    input:
        path bc
        path ft
        path mx
        val  min_c
        val  min_f
        val  mt_p
        val  proj
    output:
        path "02_seurat_unfiltered.rds", emit: rds
    script:
        """
        Rscript ${projectDir}/scripts/02_create_seurat_obj.R \
            $bc $ft $mx $min_c $min_f "$mt_p" "$proj" "02_seurat_unfiltered.rds"
        """
}

// --- PROCESS 3: QC Visualization ---
process VISUALIZE_QC {
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    input:
        path rds
    output:
        path "03_seurat_with_qc.rds", emit: rds
        path "03_qc_plots.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/03_qc_visualize.R \
            $rds "03_seurat_with_qc.rds" "03_qc_plots.pdf"
        """
}

// --- PROCESS 4: Filter Cells ---
process FILTER_CELLS {
    publishDir "${params.resultsdir}/objects", mode: 'copy'
    input:
        path rds
        val  min_g
        val  max_g
        val  max_mt
    output:
        path "04_seurat_filtered.rds", emit: rds
    script:
        """
        Rscript ${projectDir}/scripts/04_filter_cells.R \
            $rds $min_g $max_g $max_mt "04_seurat_filtered.rds"
        """
}

// --- PROCESS 5: Normalize (SCTransform) ---
process NORMALIZE_DATA {
    publishDir "${params.resultsdir}/objects", mode: 'copy'
    input:
        path rds
    output:
        path "05_seurat_normalized.rds", emit: rds
    script:
        """
        Rscript ${projectDir}/scripts/05_normalize_data.R \
            $rds "05_seurat_normalized.rds" ${task.cpus}
        """
}

// --- PROCESS 6: Dimensionality Reduction & Clustering ---
process DIM_REDUCTION_AND_CLUSTER {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    input:
        path rds
        val  dims
        val  res
        val  title
    output:
        path "06_seurat_clustered.rds", emit: rds
        path "06_elbow.pdf"
        path "06_umap.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/06_run_dim_reduction.R \
            $rds $dims $res "$title" \
            "06_seurat_clustered.rds" "06_elbow.pdf" "06_umap.pdf"
        """
}

// --- PROCESS 7: Clustree Resolution Stability ---
process RUN_CLUSTREE {
    publishDir "${params.resultsdir}/plots", mode: 'copy'
    input:
        path rds
    output:
        path "07_clustree_resolutions.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/07_run_clustree.R \
            $rds "07_clustree_resolutions.pdf"
        """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    DOWNLOAD_DATA(
        params.barcodes_url,
        params.features_url,
        params.matrix_url
    )

    CREATE_SEURAT_OBJECT(
        DOWNLOAD_DATA.out.barcodes,
        DOWNLOAD_DATA.out.features,
        DOWNLOAD_DATA.out.matrix,
        params.min_cells,
        params.min_features,
        params.mt_pattern,
        params.project_name
    )

    VISUALIZE_QC(CREATE_SEURAT_OBJECT.out.rds)

    FILTER_CELLS(
        VISUALIZE_QC.out.rds,
        params.min_genes,
        params.max_genes,
        params.max_mt_percent
    )

    NORMALIZE_DATA(FILTER_CELLS.out.rds)

    DIM_REDUCTION_AND_CLUSTER(
        NORMALIZE_DATA.out.rds,
        params.pca_dims,
        params.cluster_res,
        params.project_name
    )

    RUN_CLUSTREE(DIM_REDUCTION_AND_CLUSTER.out.rds)
}
