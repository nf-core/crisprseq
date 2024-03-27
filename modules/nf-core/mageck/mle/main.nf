process MAGECK_MLE {
    tag "$prefix"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageck:0.5.9.5--py39h1f90b4d_3':
        'biocontainers/mageck:0.5.9.5--py39h1f90b4d_3' }"

    input:
    tuple val(meta), path(design_matrix), path(count_table)

    output:
    tuple val(meta), path("*.gene_summary.txt") , emit: gene_summary
    tuple val(meta), path("*.sgrna_summary.txt"), emit: sgrna_summary
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = meta.id ?: "${meta.treatment}_vs_${meta.reference}"
    def design_command = design_matrix ? "-d $design_matrix" : ''

    """
    mageck \\
        mle \\
        $args \\
        --threads $task.cpus \\
        -k $count_table \\
        -n $prefix     \\
        $design_command 
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
}
