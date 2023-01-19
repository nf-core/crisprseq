process DUMMY_FINAL_UMI {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::p7zip==16.02"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'quay.io/biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_clustered.fastq.gz"), emit: dummy
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clustered.fastq

    7za \\
        a \\
        $args \\
        ${prefix}_clustered.fastq.gz \\
        ${prefix}_clustered.fastq

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
