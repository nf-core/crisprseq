
process MAGECK_FLUTEMLE {
    tag "$prefix"
    label 'process_high'

    conda "bioconda::bioconductor-mageckflute=2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageckflute:2.2.0--r42hdfd78af_0':
        'biocontainers/bioconductor-mageckflute:2.2.0--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(gene_summary)

    output:
    tuple val(meta), path("MAGeCKFlute_*/Enrichment/*"), emit: enrich
    tuple val(meta), path("MAGeCKFlute_*/QC/*"), emit: qc
    tuple val(meta), path("MAGeCKFlute_*/Selection/*"), emit: select
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = meta.id ?: "${meta.treatment}_vs_${meta.reference}"

    """
    #!/usr/bin/env Rscript
    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    #### graphs mageck MLE

    library(MAGeCKFlute)
    library(clusterProfiler)
    library(ggplot2)

    #library(pathview)
    options(ggrepel.max.overlaps = Inf)
    mle <- read.table("${gene_summary}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    print(head(mle))
    colnames(mle)
    #stop(colnames(mle))
    FluteMLE(mle, treatname= "${prefix}", proj="${prefix}", pathview.top=0, $args)

    version_file_path <- "versions.yml"
    version_flute <- paste(unlist(packageVersion("MAGeCKFlute")), collapse = ".")
    version_ggplot <- paste(unlist(packageVersion("ggplot2")), collapse = ".")

    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    MAGeCKFlute: ", f, sep = "")
    writeLines(version_flute, f)
    writeLines("    ggplot2: ", f, sep = "")
    writeLines(version_ggplot, f)
    close(f)
    """
}
