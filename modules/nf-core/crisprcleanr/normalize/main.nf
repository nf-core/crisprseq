process CRISPRCLEANR_NORMALIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-crisprcleanr:3.0.0--r42hdfd78af_1':
        'biocontainers/r-crisprcleanr:3.0.0--r42hdfd78af_1' }"

    input:
    tuple val(meta), path(count_file), val(library_value), path(library_file)
    val(min_reads)
    val(min_targeted_genes)

    output:
    tuple val(meta), path("*_norm_table.tsv"), emit: norm_count_file
    tuple val(meta), path("*.RData"),          emit: counts_rdata
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(CRISPRcleanR)
    library(dplyr)
    
    


    """
}
