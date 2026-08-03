process CRISPRDECODE_ASSIGN {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11.4"
    container 'python:3.11.4-slim-bookworm@sha256:17d62d681d9ecef20aae6c6605e9cf83b0ba3dc247013e2f43e1b5a045ad4901'

    input:
    tuple val(meta), path(reads), path(library)
    val r1_anchor
    val r2_anchor
    val r1_offset
    val r2_offset
    val reverse_complement_r1
    val reverse_complement_r2

    output:
    tuple val(meta), path("*.crisprdecode.counts.tsv"),             emit: counts
    tuple val(meta), path("*.crisprdecode.assignment_summary.tsv"), emit: summary
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r1_anchor_arg = r1_anchor ? "--r1-anchor ${r1_anchor}" : ''
    def r2_anchor_arg = r2_anchor ? "--r2-anchor ${r2_anchor}" : ''
    def reverse_complement_r1_arg = reverse_complement_r1 ? '--reverse-complement-r1' : ''
    def reverse_complement_r2_arg = reverse_complement_r2 ? '--reverse-complement-r2' : ''
    """
    crisprdecode_assign.py \
        --r1 ${reads[0]} \
        --r2 ${reads[1]} \
        --library $library \
        --sample ${meta.id} \
        $r1_anchor_arg \
        $r2_anchor_arg \
        --r1-offset $r1_offset \
        --r2-offset $r2_offset \
        $reverse_complement_r1_arg \
        $reverse_complement_r2_arg \
        --counts ${prefix}.crisprdecode.counts.tsv \
        --summary ${prefix}.crisprdecode.assignment_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.crisprdecode.counts.tsv
    touch ${prefix}.crisprdecode.assignment_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
