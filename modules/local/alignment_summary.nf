process ALIGNMENT_SUMMARY {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(reads)//, path(summary)

    output:
    //tuple val(meta), path(summary), emit: summary
    tuple val(meta), path("*_summary.csv"), emit: summary
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
	mapped_count=`samtools view -c -b -F 4 $reads`
	total=`samtools view -c -b $reads`
	perc=`echo "scale=1;(\$mapped_count/\$total)*100" | bc`
	echo "aligned-reads," \$mapped_count "("\$perc"%)" >> ${prefix}_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
