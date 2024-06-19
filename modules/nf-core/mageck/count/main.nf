process MAGECK_COUNT {
    tag "$meta.id"
    label 'process_high'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageck:0.5.9.5--py39h1f90b4d_3':
        'biocontainers/mageck:0.5.9.5--py39h1f90b4d_3' }"

    input:
    tuple val(meta), path(fastq1), path(fastq2)
    path(library)

    output:
    tuple val(meta), path("*count.txt")            , emit: count
    tuple val(meta), path("*.count_normalized.txt"), emit: norm
    tuple val(meta), path("*.countsummary.txt")    , emit: summary
    tuple val(meta), path("*.count_normalized.txt"), emit: normalized
    tuple val(meta), path("*.log")                 , emit: logs
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_label = ("$fastq1".endsWith(".fastq.gz") || "$fastq1".endsWith(".fq.gz")) ? "--sample-label ${meta.id}" : ''
    
    if (meta.single_end && ("$fastq1".endsWith(".fastq.gz") || "$fastq1".endsWith(".fq.gz")) || "$fastq1".endsWith(".bam")) {
        input = "--fastq $fastq1" 
    } else {
        input = "--fastq $fastq1 --fastq-2 $fastq2" 
    }

    """
    mageck \\
        count \\
        $args \\
        -l $library \\
        -n $prefix \\
        $sample_label \\
        $input \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
}