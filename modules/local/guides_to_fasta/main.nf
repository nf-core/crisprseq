process GUIDES_TO_FASTA {
    tag "$table"
    label 'process_single'

    conda "bioconda::mawk=1.3.4 conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    publishDir "${params.outdir}/reference", enabled: false, mode:'copy'

    input:
    path table

    output:
    path "*.fasta",      emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def delim = table.toString().endsWith('.csv') ? ',' : '\t'
    """
    awk -F'$delim' '(NR > 1) {print ">"\$1"\\n"\$2}' ${table} > ${table.baseName}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(mawk -W version 2> /dev/null | sed -n 's/^mawk \\([^\\n]*\\).*/\\1/p')
        sed: \$(sed --version | sed -n 's/sed (GNU sed) \\([^\\n]*\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    """
    touch ${table.basename}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(mawk -W version 2> /dev/null | sed -n 's/^mawk \\([^\\n]*\\).*/\\1/p')
        sed: \$(sed --version | sed -n 's/sed (GNU sed) \\([^\\n]*\\).*/\\1/p')
    END_VERSIONS
    """
}
