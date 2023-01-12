process CLUSTERING_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::p7zip==16.02" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'quay.io/biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), path(reads), path(summary)

    output:
    tuple val(meta), path(summary), emit: summary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    7za \\
        e \\
        $args \\
        $reads \\
        > NUL

    pre_final_count=`cat ${reads.baseName} | wc -l`
    final_count=`echo \$(( "\$pre_final_count/4" ))`
    echo "clustered-reads," \$final_count >> $summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}



process summary_reads{

    input:
    set sampleID, file(summary), file(finalReads) from summaryReport.join(readsSummary)

    output:
    set val(sampleID), file(summary) into summaryReportReads

    script:
    """
    final_count=`expr \$(cat $finalReads | wc -l) / 4`
    echo "clustered-reads," \$final_count >> ${sampleID}_summary.csv
    """

    }
