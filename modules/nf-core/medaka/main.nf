process MEDAKA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.3--py38h2e44183_0' :
        'biocontainers/medaka:1.11.3--py38h2e44183_0' }"

    input:
    tuple val(meta), path(reads), path(assembly)
    path model

    output:
    tuple val(meta), path("*_medakaConsensus.fasta"), emit: assembly
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ "${reads.extension}" == "gz" ]]; then
        gzip -df $reads
        reads=$reads.baseName
    else
        reads=$reads
    fi
    if [[ "${assembly.extension}" == "gz" ]]; then
        gzip -df $assembly
        assembly=$assembly.baseName
    else
        assembly=$assembly
    fi

    medaka_consensus \\
        -t $task.cpus \\
        $args \\
        -i \$reads \\
        -d \$assembly \\
        -o ./

    mv consensus.fasta ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
