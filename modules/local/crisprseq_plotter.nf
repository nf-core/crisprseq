process CRISPRSEQ_PLOTTER {
    tag "$meta.id"
    label 'process_medium'

    conda 'r-ggplot2=3.4.3 bioconductor-shortread=1.58.0 r-ggpubr=0.6.0 r-ggmsa=1.0.2 r-seqmagick=0.1.6 r-tidyr=1.3.0 r-ggseqlogo=0.1 r-cowplot=1.1.1 r-seqinr=4.2_30 r-optparse=1.7.3 r-dplyr=1.1.2 r-plyr=1.8.8 r-stringr=1.5.0 r-plotly=4.10.2'
    // TODO: change containers
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-61c59287265e27f2c4589cfc90013ef6c2c6acf1:fb3e48060a8c0e5108b1b60a2ad7e090cfb9eee5-0' :
        'biocontainers/mulled-v2-61c59287265e27f2c4589cfc90013ef6c2c6acf1:fb3e48060a8c0e5108b1b60a2ad7e090cfb9eee5-0' }"

    input:
    tuple val(meta), path(indels), path(substitutions), path(reference), val(protospacer)

    output:
    tuple val(meta), path("*_Deletions.html"), path("*_delAlleles_plot.png")        , emit: deletions
    tuple val(meta), path("*_Insertions.html")                                      , emit: insertions
    tuple val(meta), path("*_accumulative.html")                                    , emit: accumulative
    tuple val(meta), path("*_subs-perc_plot.png"), path("*_subs-perc_plot_LOGO.png"), emit: substitutions
    tuple val(meta), path("*_top-alleles_LOGO.png"), path("*_top.html")             , emit: topalleles
    path "versions.yml"                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotter.R \\
        $args \\
        --indels_info=$indels \\
        --reference=$reference \\
        --gRNA_sequence=$protospacer \\
        --sample_name=$prefix \\
        --substitutions_info=$substitutions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(Rscript --version)
    END_VERSIONS
    """
}
