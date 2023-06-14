process CRISPRCLEANR_NORMALIZE {
    tag "$meta"
    label 'process_medium'

    conda "bioconda::r-crisprcleanr=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-crisprcleanr:3.0.0--r42hdfd78af_1':
        'quay.io/biocontainers/r-crisprcleanr:3.0.0--r42hdfd78af_1' }"

    input:
    tuple val(meta), path(count_file), val(library_file)
    val(min_reads)
    val(min_targeted_genes)

    output:
    tuple val(meta), path("*_norm_table.tsv"), emit: norm_count_file
    tuple val(meta), path("*.RData"),          emit: counts_rdata
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    #!/usr/bin/env Rscript
    library(CRISPRcleanR)
    library(dplyr)
    data('${library_file}')
    count_file <- read.delim('${count_file}',header=T,sep = "\t")
    count_file_to_normalize <- count_file  %>% dplyr::left_join(get('${library_file}'), by=c("sgRNA"="Target.Context.Sequence"),multiple = "all")

    count_file_to_normalize <- count_file_to_normalize %>% 
        dplyr::select(colnames(count_file),CODE,-sgRNA)

    names(count_file_to_normalize)[names(count_file_to_normalize) == 'Gene'] <- 'gene'
    names(count_file_to_normalize)[names(count_file_to_normalize) == 'CODE'] <- 'sgRNA'
    count_file_to_normalize <- count_file_to_normalize %>% dplyr::select(sgRNA, gene, everything())

    #crisprcleanr function
    normANDfcs <- ccr.NormfoldChanges(Dframe=count_file_to_normalize,saveToFig = FALSE,min_reads=${min_reads},EXPname="${prefix}", libraryAnnotation=get('${library_file}'),display=FALSE)
    gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs[["logFCs"]],get('${library_file}'))
    correctedFCs <- ccr.GWclean(gwSortedFCs,display=FALSE,label='${meta}')
    correctedCounts <- ccr.correctCounts('${meta}',
                            normANDfcs[["norm_counts"]],
                            correctedFCs,
                            get('${library_file}'),
                            minTargetedGenes=${min_targeted_genes},
                            OutDir='./')

    write.table(correctedCounts, file=paste0("${prefix}","_norm_table.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

    #version
    version_file_path <- "versions.yml"
    version_crisprcleanr <- paste(unlist(packageVersion("CRISPRcleanR")), collapse = ".")
    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    crisprcleanr: ", f, sep = "")
    writeLines(version_crisprcleanr, f)
    close(f)

    """
}
