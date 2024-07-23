process CLONALITY_CLASSIFIER {
    tag "$meta.id"
    label 'process_single'

    conda "pandas=2.2.0,numpy=1.26.3,statsmodels=0.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d836da785124bb367cbe6fbfc00dddd2107a4da:b033d6a4ea3a42a6f5121a82b262800f1219b382-0' :
        'biocontainers/mulled-v2-9d836da785124bb367cbe6fbfc00dddd2107a4da:b033d6a4ea3a42a6f5121a82b262800f1219b382-0' }"

    input:
    tuple val(meta), path(indels), path(edition)


    output:
    tuple val(meta), path("*_edits_classified.csv"), emit: classified
    path "versions.yml",                             emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    template 'clonality_classifier.py'
}
