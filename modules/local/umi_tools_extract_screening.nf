process UMI_TOOLS_EXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::umi_tools=1.1.5'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.5--py39hf95cd2a_0' :
        'biocontainers/umi_tools:1.1.5--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_extracted*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")                , emit: log
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extract_method = params.umi_extract_method ?: 'string'
    def pattern = params.umi_pattern ? "--bc-pattern='${params.umi_pattern}'" : ''
    def pattern2 = params.umi_pattern2 ? "--bc-pattern2='${params.umi_pattern2}'" : ''
    
    if (meta.single_end) {
        """
        umi_tools extract \\
            --extract-method=${extract_method} \\
            ${pattern} \\
            ${args} \\
            -I ${reads[0]} \\
            -S ${prefix}_extracted.fastq.gz \\
            -L ${prefix}_umi_extract.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umi_tools: \$(umi_tools --version 2>&1 | grep -oP 'UMI-tools version: \\K[0-9.]+' || echo '1.1.5')
        END_VERSIONS
        """
    } else {
        """
        umi_tools extract \\
            --extract-method=${extract_method} \\
            ${pattern} \\
            ${pattern2} \\
            ${args} \\
            -I ${reads[0]} \\
            --read2-in=${reads[1]} \\
            -S ${prefix}_extracted_R1.fastq.gz \\
            --read2-out=${prefix}_extracted_R2.fastq.gz \\
            -L ${prefix}_umi_extract.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umi_tools: \$(umi_tools --version 2>&1 | grep -oP 'UMI-tools version: \\K[0-9.]+' || echo '1.1.5')
        END_VERSIONS
        """
    }
}
