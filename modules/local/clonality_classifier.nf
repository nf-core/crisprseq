process CLONALITY_CLASSIFIER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::biopython=1.78"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(indels), path(edition)


    output:
    tuple val(meta), path("*_edits_classified.csv"), emit: classified


    when:
    task.ext.when == null || task.ext.when

    script:
    template 'clonality_classifier.py $indels, $edition'
}
