process CRISPRDECODE_AGGREGATE {
    tag "aggregate"
    label 'process_single'

    conda "conda-forge::python=3.11.4"
    container 'python:3.11.4-slim-bookworm@sha256:17d62d681d9ecef20aae6c6605e9cf83b0ba3dc247013e2f43e1b5a045ad4901'

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
