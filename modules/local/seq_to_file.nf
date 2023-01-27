process SEQ_TO_FILE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::p7zip==16.02"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'quay.io/biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), val(sequence)
    val(type)

    output:
    tuple val(meta), path('*.fasta') , emit: file
    path "versions.yml"              , emit: versions

    when:
    sequence != null

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ">${meta.id}" > ${meta.id}_${type}.fasta
    echo $sequence >> ${meta.id}_${type}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
