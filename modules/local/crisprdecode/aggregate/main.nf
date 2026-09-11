process CRISPRDECODE_AGGREGATE {
    tag "aggregate"
    label 'process_single'

    conda "conda-forge::python=3.11.4"
    container 'docker.io/library/python:3.11.4-bookworm@sha256:d7df302a1bcf4db50650da79c174f5d8d973fa4753e0275696644c5bdb477c00'

    input:
    path sample_counts
    path sample_summaries
    path library

    output:
    path "count_table.count.txt",   emit: count_matrix
    path "assignment_summary.tsv",  emit: assignment_summary
    path "library_recovery.tsv",    emit: library_recovery
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    crisprdecode_aggregate.py \
        --library $library \
        --sample-counts ${sample_counts.join(' ')} \
        --sample-summaries ${sample_summaries.join(' ')} \
        --count-matrix count_table.count.txt \
        --assignment-summary assignment_summary.tsv \
        --library-recovery library_recovery.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch count_table.count.txt
    touch assignment_summary.tsv
    touch library_recovery.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
