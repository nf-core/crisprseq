process MAGECK_GRAPHRRA {
    tag "$meta.treatment"
    label 'process_single'

    conda "bioconda::bioconductor-mageckflute=2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageckflute:2.2.0--r42hdfd78af_0':
        'biocontainers/bioconductor-mageckflute:2.2.0--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(gene_summary)

    output:
    tuple val(meta), path("*.png"), emit: graphs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    #### Orient a reference sequence according to reads orientation.

    library(MAGeCKFlute)
    library(ggplot2)
    options(ggrepel.max.overlaps = Inf)

    gdata = ReadRRA("$gene_summary")
    gdata <- transform(gdata, LogFDR = -log10(FDR))
    png(filename = paste0("$meta.treatment","_vs_","$meta.reference","_scatterview.png"), width = 6, height = 4, units = "in", res = 300)
    p1 = ScatterView(gdata, x = "Score", y = "LogFDR", label = "id",
                model = "volcano", top = 5)
    print(p1)
    dev.off()

    gdata <- transform(gdata, Rank = rank(Score))
    png(filename = paste0("$meta.treatment","_vs_","$meta.reference","_rank.png"), width = 6, height = 4, units = "in", res = 300)
    p1 = ScatterView(gdata, x = "Rank", y = "Score", label = "id",
                top = 5, auto_cut_y = TRUE, ylab = "Log2FC",
                groups = c("top", "bottom"))
    print(p1)
    dev.off()

    version_file_path <- "versions.yml"
    version_flute <- paste(unlist(packageVersion("MAGeCKFlute")), collapse = ".")
    version_ggplot <- paste(unlist(packageVersion("MAGeCKFlute")), collapse = ".")

    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    MAGeCKFlute: ", f, sep = "")
    writeLines(version_flute, f)
    writeLines("    ggplot2: ", f, sep = "")
    writeLines(version_ggplot, f)
    close(f)




    """


}
