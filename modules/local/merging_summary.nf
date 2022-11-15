process MERGING_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::p7zip==16.02" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'quay.io/biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), path(raw_reads), path(assembled_reads), path(trimmed_reads), path(trimmed_adapters)


    output:
    tuple val(meta), path("*_summary.csv"), emit: summary
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    7za \\
        e \\
        $args \\
        ${raw_reads[0]} \\
        > NUL

    7za \\
        e \\
        $args \\
        $assembled_reads \\
        > NUL

    7za \\
        e \\
        $args \\
        $trimmed_reads \\
        > NUL

	#### Counts
	pre_raw_count=`cat ${raw_reads[0].baseName} | wc -l`
	raw_count=`echo "\$pre_raw_count/4" | bc`
	pre_merged_count=`cat ${assembled_reads.baseName} | wc -l`
	merged_count=`echo "\$pre_merged_count/4" | bc`
	adapters=`grep "Reads with adapters" ${trimmed_adapters} | awk -F " " '{ print \$(NF-1)" "\$NF }'`
	pre_filt=`cat ${trimmed_reads.baseName} | wc -l`
	filt=`echo "\$pre_filt/4" | bc`

	#### Percentages
	merged_perc=`echo "scale=1;(\$merged_count/\$raw_count)*100" | bc`
	filt_perc=`echo "scale=1;(\$filt/\$merged_count)*100" | bc`

	#### Table
	echo "class, count" > ${prefix}_summary.csv
	echo "raw-reads," \$raw_count "(100.0%)" >> ${prefix}_summary.csv
	echo "merged-reads," \$merged_count "("\$merged_perc"%)" >> ${prefix}_summary.csv
	echo "reads-with-adapters," \${adapters//,} >> ${prefix}_summary.csv
	echo "quality-filtered-reads," \$filt "("\$filt_perc"%)" >> ${prefix}_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
