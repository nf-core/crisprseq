process FIND_ADAPTERS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::r-fastqcr==0.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'quay.io/biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), path(zip)

    output:
    tuple val(meta), path("overrepresented.fasta"),     emit: adapters
    path "versions.yml",                                emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    library(fastqcr)
    library(magrittr)

    unzip("$zip", files = c("${zip.baseName}/fastqc_data.txt"))

    fastqcr::qc_read("${zip.baseName}/fastqc_data.txt")\$overrepresented_sequences %>%
    dplyr::mutate(name=paste(">",1:dplyr::n(),"-",Count,sep=""),fa=paste(name,Sequence,sep="\n")) %>%
    dplyr::pull(fa) %>%
    readr::write_lines("overrepresented.fasta")

    version_file_path <- "versions.yml"
    version_ampir <- paste(unlist(packageVersion("fastqcr")), collapse = ".")
    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    fastqcr: ", f, sep = "")
    writeLines(version_ampir, f)
    close(f)
    """
}
