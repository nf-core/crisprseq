process UMI_TO_SEQUENCE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::p7zip==16.02"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'quay.io/biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), path(cluster)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'umi_to_sequence.py'
}
