
process MAGECK_FLUTEMLE {
    tag "$prefix"
    label 'process_high'

    conda "bioconda::bioconductor-mageckflute==2.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-mageckflute:2.6.0--r43hdfd78af_0':
        'biocontainers/bioconductor-mageckflute:2.6.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(gene_summary)

    output:
    tuple val(meta), path("MAGeCKFlute_*/Enrichment/*") , emit: enrich
    tuple val(meta), path("MAGeCKFlute_*/QC/*")         , emit: qc
    tuple val(meta), path("MAGeCKFlute_*/Selection/*")  , emit: select
    tuple val(meta), path("MAGeCKFlute_*/PathwayView/*"), emit: pathwayview
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ' '
    prefix = meta.id ?: "${meta.treatment}_vs_${meta.reference}"
    template 'template_fluteMLE.R'
}
