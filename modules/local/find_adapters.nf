process FIND_ADAPTERS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::p7zip==15.09" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:15.09--h2d50403_4' :
        'quay.io/biocontainers/p7zip:15.09--h2d50403_4' }"

    input:
    tuple val(meta), path(zip)

    output:
    tuple val(meta), env(test_lines), env(adapter_seq), emit: adapters
    path "versions.yml",                                emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    7za \\
        e \\
        -o"${zip.baseName}"/ \\
        $args \\
        $zip
    test_lines=`grep -A 4 "Overrepresented sequences" $zip.baseName/fastqc_data.txt | grep -v "No Hit" | head -3 | tail -1 | awk '{print NF}'`
    adapter_seq=`grep -A 4 "Overrepresented sequences" $zip.baseName/fastqc_data.txt | grep -v "No Hit" | head -3 | tail -1 | awk -F "\t" '{ print \$1 }'`

    # TODO: one the bug mentioned in https://github.com/nextflow-io/nextflow/issues/2812 is fixed, remove following lines
    echo "test_lines=\$test_lines" > .command.env
    echo "adapter_seq=\$adapter_seq" >> .command.env
    source .command.env

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
