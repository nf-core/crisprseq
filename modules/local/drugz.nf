process DRUGZ {
    tag "${meta.treatment}_${meta.reference}"
    label 'process_single'

    
    conda "python=3.11.4 pandas=2.0.3 numpy=1.25.1 scipy=1.11.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1ec3f69e7819b1ab3e6f57d16594eb40ed7d6792:f94a27287a1921ce0dacd411d48acff738d3ca90-0':
        'biocontainers/mulled-v2-1ec3f69e7819b1ab3e6f57d16594eb40ed7d6792:f94a27287a1921ce0dacd411d48acff738d3ca90-0' }"

    input:
    tuple val(meta), path(count_table)

    output:
    tuple val(meta), path("*"), emit: per_gene_results
    path "versions.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    drugz.py -i $count_table -o ${meta.treatment}_vs_${meta.reference} -c $meta.reference -x $meta.treatment $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
