process CRISPRDECODE_VALIDATE_LIBRARY {
    tag "$library"
    label 'process_single'

    conda "conda-forge::python=3.11.4"
    container 'python:3.11.4-slim-bookworm@sha256:17d62d681d9ecef20aae6c6605e9cf83b0ba3dc247013e2f43e1b5a045ad4901'

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
