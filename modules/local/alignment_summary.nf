process ALIGNMENT_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pysam=0.20.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.20.0--py310hff46b53_0' :
        'quay.io/biocontainers/pysam:0.20.0--py310hff46b53_0' }"

    input:
    tuple val(meta), path(reads), path(summary)

    output:
    tuple val(meta), path(summary), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'alignment_summary.py'
}
