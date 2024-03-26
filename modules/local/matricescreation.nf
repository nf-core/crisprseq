process MATRICESCREATION {
    label 'process_single'

    conda 'r-ggplot2=3.4.3 bioconductor-shortread=1.58.0 r-ggpubr=0.6.0 r-ggmsa=1.0.2 r-seqmagick=0.1.6 r-tidyr=1.3.0 r-ggseqlogo=0.1 r-cowplot=1.1.1 r-seqinr=4.2_30 r-optparse=1.7.3 r-dplyr=1.1.2 r-plyr=1.8.8 r-stringr=1.5.0 r-plotly=4.10.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-6de07928379e6eface08a0019c4a1d6b5192e805:0d77388f37ddd923a087f7792e30e83ab54c918c-0' :
        'biocontainers/mulled-v2-6de07928379e6eface08a0019c4a1d6b5192e805:0d77388f37ddd923a087f7792e30e83ab54c918c-0' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("*.txt"), emit: design_matrix

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    #!/usr/bin/env Rscript
    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    ####

    # Loop through each row in the data
    control_samples <- unlist(strsplit('${meta.reference}', ","))
    treatment_samples <- unlist(strsplit('$meta.treatment', ","))
    all_samples <- unique(c(control_samples, treatment_samples))
    design_matrix <- data.frame(matrix(0, nrow = length(all_samples), ncol = 3,
                                dimnames = list(all_samples,
                                                c("Samples", "baseline",
    paste0(gsub(',', '_', '$meta.treatment'),"_vs_",gsub(',','_','$meta.reference'))))))
    name = paste0(gsub(',', '_', '$meta.treatment' ),"_vs_", gsub(',', '_','$meta.reference'))
    # Set baseline and treatment values in the design matrix
    design_matrix[, "Samples"] <- rownames(design_matrix)
    design_matrix\$baseline <- 1
    design_matrix[treatment_samples, name] <- 1
    design_matrix[treatment_samples, paste0(gsub(',', '_', '$meta.treatment'),"_vs_",gsub(",","_",'$meta.reference'))] <- 1

    # Print the design matrix to a file
    output_file <- paste0(gsub(',', '_', '$meta.treatment' ),"_vs_",gsub(",","_",'$meta.reference'),".txt")
    write.table(design_matrix, output_file, sep = "\t", quote = FALSE, row.names=FALSE)

    """
}
