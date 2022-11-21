process ORIENT_REFERENCE {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "r-seqinr=4.2_16 bioconductor-biostrings=2.62.0 bioconductor-shortread=1.52.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-dd29a726a518ab18d2113433052947ad2389bdcd:04842a4c570550c67fec3f957a1aae3bd0771965-0' :
        'quay.io/biocontainers/mulled-v2-dd29a726a518ab18d2113433052947ad2389bdcd:04842a4c570550c67fec3f957a1aae3bd0771965-0' }"

    input:
    tuple val(meta), file(reference), val(protospacer)

    output:
    tuple val(meta), path('*_reference-correctOrient.fasta') , emit: reference
    path "versions.yml"                                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    revComp_reference.R \\
        $reference \\
        ${meta.id}_reference-correctOrient.fasta \\
        $protospacer;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(Rscript --version)
    END_VERSIONS
    """
}
