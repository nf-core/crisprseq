process CRISPRDECODE_VALIDATE_LIBRARY {
    tag "$library"
    label 'process_single'

    conda "conda-forge::python=3.11.4"
    container 'docker.io/library/python:3.11.4-bookworm@sha256:d7df302a1bcf4db50650da79c174f5d8d973fa4753e0275696644c5bdb477c00'

    input:
    path library

    output:
    path "validated_construct_library.tsv", emit: library
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    crisprdecode_validate_library.py \
        --library $library \
        --output validated_construct_library.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch validated_construct_library.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
