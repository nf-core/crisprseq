process MATRICESCREATION {
    label 'process_single'

    conda 'r-ggplot2=3.4.3 bioconductor-shortread=1.58.0 r-ggpubr=0.6.0 r-ggmsa=1.0.2 r-seqmagick=0.1.6 r-tidyr=1.3.0 r-ggseqlogo=0.1 r-cowplot=1.1.1 r-seqinr=4.2_30 r-optparse=1.7.3 r-dplyr=1.1.2 r-plyr=1.8.8 r-stringr=1.5.0 r-plotly=4.10.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-6de07928379e6eface08a0019c4a1d6b5192e805:0d77388f37ddd923a087f7792e30e83ab54c918c-0' :
        'biocontainers/mulled-v2-6de07928379e6eface08a0019c4a1d6b5192e805:0d77388f37ddd923a087f7792e30e83ab54c918c-0' }"

    input:
    path(contrasts)

    output:
    path("*.txt"), emit: design_matrix

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    #!/usr/bin/env Rscript
    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    #### Create design matrices

    data <- read.table("$contrasts", header = TRUE, sep = ";", stringsAsFactors = FALSE)
    # Loop through each row in the data
    for (i in 1:nrow(data)) {
        # Extract control and treatment samples for the current row
        control_samples <- unlist(strsplit(data\$reference[i], ","))
        treatment_samples <- unlist(strsplit(data\$treatment[i], ","))

        # Create a vector of all unique samples
        all_samples <- unique(c(control_samples, treatment_samples))

        # Initialize a matrix to store the design matrix
        design_matrix <- data.frame(matrix(0, nrow = length(all_samples), ncol = 3,
                        dimnames = list(all_samples, c("Samples", "baseline", paste0(gsub(',', '_', data\$treatment[i] ),"_vs_", data\$reference[i])))))

        # Set baseline and treatment values in the design matrix
        design_matrix[, "Samples"] <- rownames(design_matrix)
        design_matrix\$baseline <- 1
        design_matrix[treatment_samples, paste0(gsub(',', '_', data\$treatment[1] ),"_vs_",gsub(",","_",data\$reference[i]))] <- 1

        # Print the design matrix to a file
        output_file <- paste0(gsub(',', '_', data\$treatment[1] ),"_vs_",gsub(",","_",data\$reference[i]),".txt")
        write.table(design_matrix, output_file, sep = "\t", quote = FALSE, row.names=FALSE)
    }
    """
}
