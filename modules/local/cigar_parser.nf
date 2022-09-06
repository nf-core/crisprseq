process CIGAR_PARSER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::cutadapt=3.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
        'quay.io/biocontainers/cutadapt:3.4--py39h38f01e4_1' }"

    input:
    tuple val(meta), path(reads), path(reference), val(protospacer), path(template), path(summary)

    output:
	tuple val(meta), path("*indels.csv"), path("*_subs-perc.csv"), emit: indels
	tuple val(meta), path("*.html"), path("*edits.csv")          , emit: edition
    tuple val(meta), path("*cutSite.json")                       , emit: cutsite
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cigar_parser.R \\
        $args \\
        --input = $reads \\
        --output = $prefix \\
        --reference = $reference \\
        --gRNA_sequence = $protospacer \\
        --sample_name = $prefix \\
        --template = $template \\
        --summary_file = $summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version)
    END_VERSIONS
    """
}
