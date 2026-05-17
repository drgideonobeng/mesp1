nextflow.enable.dsl=2

// =============================================================================
// phase2/main.nf — Biological Insight Pipeline (Mesp1-cardiopharyngeal lineage)
// Consumes: 06_seurat_clustered.rds (Phase 1 output)
// Steps: find_markers → plot_marker_genes → annotate_clusters
// =============================================================================

// --- PROCESS 1: Find Marker Genes ---
process FIND_MARKERS {
    publishDir "${params.resultsdir}/tables", mode: 'copy'
    input:
        path rds
    output:
        path "01_all_markers.csv",              emit: csv
        path "01_top_5_markers_per_cluster.csv"
    script:
        """
        Rscript ${projectDir}/scripts/01_find_markers.R \
            $rds "01_all_markers.csv" ${task.cpus} "01_top_5_markers_per_cluster.csv"
        """
}

// --- PROCESS 2: Plot Marker Gene Feature Maps ---
process PLOT_MARKER_GENES {
    publishDir "${params.resultsdir}/plots", mode: 'copy'
    input:
        path rds
        val  genes
    output:
        path "02_marker_feature_plots.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/02_plot_marker_genes.R \
            $rds "$genes" "02_marker_feature_plots.pdf"
        """
}

// --- PROCESS 3: Annotate Clusters ---
process ANNOTATE_CLUSTERS {
    publishDir "${params.resultsdir}/objects", mode: 'copy', pattern: "*.rds"
    publishDir "${params.resultsdir}/plots",   mode: 'copy', pattern: "*.pdf"
    input:
        path rds
        val  proj
        path labels
    output:
        path "03_seurat_annotated.rds", emit: rds
        path "03_umap_annotated.pdf"
    script:
        """
        Rscript ${projectDir}/scripts/03_annotate_clusters.R \
            $rds "$proj" "$labels" \
            "03_seurat_annotated.rds" "03_umap_annotated.pdf"
        """
}

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    if (!params.clustered_rds) {
        error "params.clustered_rds is required. Provide it in your params file."
    }

    def clustered_file = file(params.clustered_rds)
    if (!clustered_file.exists()) {
        error "Phase 1 output not found: ${params.clustered_rds}\nRun Phase 1 first."
    }

    rds_ch = Channel.value(clustered_file)

    FIND_MARKERS(rds_ch)
    PLOT_MARKER_GENES(rds_ch, params.marker_genes)
    ANNOTATE_CLUSTERS(
        rds_ch,
        params.project_name,
        params.cluster_labels
    )
}
