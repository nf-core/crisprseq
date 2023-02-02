process TEMPLATE_REFERENCE {
    tag "$meta.id"
    label 'process_medium'

    conda 'r-seqinr=4.2_16 r-optparse=1.7.3 bioconductor-rsamtools=2.10.0 r-dplyr=1.0.10 r-plyr=1.8.7 r-stringr=1.4.1 bioconductor-shortread=1.52.0 bioconductor-genomicalignments=1.30.0 r-data.table=1.14.4 bioconductor-biomart=2.50.0 r-plotly=4.10.0 bioconductor-decipher=2.22.0 bioconductor-biostrings=2.62.0 r-parallelly=1.32.1 r-jsonlite=1.8.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a1551cc3e69024a9d829751175a0fe30dd5035e8:666109dfd8ca00d6bbfab5188ff2735bfcc400d7-0' :
        'quay.io/biocontainers/mulled-v2-a1551cc3e69024a9d829751175a0fe30dd5035e8:666109dfd8ca00d6bbfab5188ff2735bfcc400d7-0' }"

    input:
    tuple val(meta), path(reference), path(template)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    template_reference.R \\
        $args \\
        --reference=$reference \\
        --template=$template \\
        --prefix=$prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    """
}
