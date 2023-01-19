process EXTRACT_UMIS {
    tag "$meta.id"
    label 'process_low'

    conda "pysam=0.19.1 edlib=1.3.9 tqdm=4.64.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-76fd6f1f9ee3b11637ca84d7dc19fbaae82046a1:95aaafb5dc1049377454ee94b01f176b26ec0d57-0' :
        'quay.io/biocontainers/mulled-v2-76fd6f1f9ee3b11637ca84d7dc19fbaae82046a1:95aaafb5dc1049377454ee94b01f176b26ec0d57-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(fasta) , emit: fasta
    tuple val(meta), path(tsv)   , emit: tsv
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_umis.py \
        --threads $task.cpus \
        -o ${prefix}_extractedUMI.fasta \
        --tsv ${prefix}_extractedUMI.tsv \
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
