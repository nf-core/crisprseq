process MERGING_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::biopython=1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(raw_reads), path(assembled_reads), path(trimmed_reads), path(trimmed_adapters)


    output:
    tuple val(meta), path("*_summary.csv"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'merging_summary.py'
}
