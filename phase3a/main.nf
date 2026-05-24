nextflow.enable.dsl=2

// =============================================================================
// phase3a/main.nf — Integration Atlas Pipeline (Mesp1-cardiopharyngeal lineage)
// FAN-IN topology: PREP_NORMALIZE runs once per timepoint in parallel;
// INTEGRATE_FASTMNN collects all four outputs before proceeding.
// Steps: prep_normalize (×4) → integrate_fastmnn → cluster_integrated →
//        run_clustree | recover_cluster3 → find_markers → annotate_atlas
// =============================================================================

// --- PROCESS 1: Per-Sample Normalize (runs 4× in parallel) ---
process PREP_NORMALIZE {
    tag "$sample_id"
    publishDir "${params.resultsdir}/objects", mode: 'copy'
    input:
        tuple val(sample_id), path(filtered_rds)
    output:
        path "01_${sample_id}_normalized.rds", emit: rds
    script:
        """
        Rscript ${projectDir}/scripts/01_prep_normalize.R \
            $filtered_rds "$sample_id" \
            "01_${sample_id}_normalized.rds" ${params.n_features}
        """
}

// --- PROCESS 2: fastMNN Batch-Correction (FAN-IN: waits for all 4 samples) ---
process INTEGRATE_FASTMNN {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    input:
        path normalized_rds   // collected list of all normalized objects
    output:
        path "02_seurat_integrated.rds", emit: rds
        path "02_mnn_elbow.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/02_integrate_fastmnn.R \
            "02_seurat_integrated.rds" ${task.cpus} ${params.n_features} ${params.mnn_k} \
            ${normalized_rds.join(' ')}
        """
}

// --- PROCESS 3: UMAP + Louvain Clustering on MNN Embedding ---
process CLUSTER_INTEGRATED {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    publishDir "${params.resultsdir}/tables",  mode: 'copy', pattern: "*.csv"
    input:
        path rds
        val  mnn_dims
        val  cluster_res
        val  project_name
    output:
        path "03_seurat_clustered.rds", emit: rds
        path "03_umap_clusters.pdf"
        path "03_umap_samples.pdf"
        path "03_timepoint_composition.csv"
    script:
        """
        Rscript ${projectDir}/scripts/03_cluster_integrated.R \
            $rds $mnn_dims $cluster_res "$project_name" \
            "03_seurat_clustered.rds" "03_umap_clusters.pdf" "03_umap_samples.pdf"
        """
}

// --- PROCESS 4: Clustering Resolution Stability (Clustree) ---
process RUN_CLUSTREE {
    publishDir "${params.resultsdir}/plots", mode: 'copy'
    input:
        path rds
    output:
        path "04_clustree_resolutions.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/04_run_clustree.R \
            $rds "04_clustree_resolutions.pdf"
        """
}

// --- PROCESS 5: Cluster-3 Contamination Recovery ---
process RECOVER_CLUSTER3 {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    publishDir "${params.resultsdir}/tables",  mode: 'copy', pattern: "*.txt"
    input:
        path rds
    output:
        path "05_atlas_clean.rds",        emit: rds
        path "05_recovery_summary.txt"
        path "05_umap_atlas_clean.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/05_recover_cluster3.R \
            $rds \
            "05_atlas_clean.rds" "05_recovery_summary.txt" "05_umap_atlas_clean.pdf"
        """
}

// --- PROCESS 6: Find All Cluster Markers ---
process FIND_MARKERS {
    publishDir "${params.resultsdir}/tables", mode: 'copy'
    input:
        path rds
    output:
        path "06_all_markers.csv"
        path "06_top5_markers_per_cluster.csv"
    script:
        """
        Rscript ${projectDir}/scripts/06_find_markers.R \
            $rds \
            "06_all_markers.csv" "06_top5_markers_per_cluster.csv"
        """
}

// --- PROCESS 7: Annotate Clusters (atlas final output) ---
process ANNOTATE_ATLAS {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    publishDir "${params.resultsdir}/tables",  mode: 'copy', pattern: "*.csv"
    input:
        path rds
        val  project_name
        path labels
    output:
        path "07_seurat_annotated.rds"
        path "07_umap_annotated.pdf"
        path "07_umap_timepoints.pdf"
        path "07_all_markers.csv"
        path "07_top_5_markers_per_celltype.csv"
    script:
        """
        Rscript ${projectDir}/scripts/07_annotate_atlas.R \
            $rds "$project_name" "$labels" ${task.cpus} \
            "07_seurat_annotated.rds" "07_umap_annotated.pdf" \
            "07_umap_timepoints.pdf" "07_all_markers.csv" \
            "07_top_5_markers_per_celltype.csv"
        """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    // --- Validate inputs ---
    [
        ["filtered_e80",  params.filtered_e80],
        ["filtered_e825", params.filtered_e825],
        ["filtered_e95",  params.filtered_e95],
        ["filtered_e105", params.filtered_e105]
    ].each { name, path ->
        if (!path) error "params.${name} is required."
        def f = file(path)
        if (!f.exists()) error "Phase 1 filtered object not found for ${name}: ${path}\nRun Phase 1 first."
    }

    // --- FAN-OUT: build per-sample channel, normalize in parallel ---
    samples_ch = Channel.fromList([
        ["e80",  file(params.filtered_e80)],
        ["e825", file(params.filtered_e825)],
        ["e95",  file(params.filtered_e95)],
        ["e105", file(params.filtered_e105)]
    ])

    PREP_NORMALIZE(samples_ch)

    // --- FAN-IN: collect all 4 normalized objects, then integrate once ---
    INTEGRATE_FASTMNN(PREP_NORMALIZE.out.rds.collect())

    // --- Cluster on integrated embedding ---
    CLUSTER_INTEGRATED(
        INTEGRATE_FASTMNN.out.rds,
        params.mnn_dims,
        params.cluster_res,
        params.project_name
    )

    // --- Broadcast clustered object to two parallel downstream processes ---
    CLUSTER_INTEGRATED.out.rds
        .multiMap { rds ->
            clustree:  rds
            recover:   rds
        }
        .set { clustered }

    RUN_CLUSTREE(clustered.clustree)

    RECOVER_CLUSTER3(clustered.recover)

    FIND_MARKERS(RECOVER_CLUSTER3.out.rds)

    ANNOTATE_ATLAS(
        RECOVER_CLUSTER3.out.rds,
        params.project_name,
        params.cluster_labels
    )
}
