process MATRICESCREATION {
    tag "${meta.id}"
    label 'process_single'

    conda 'conda-forge::r-base==4.0'
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
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####

    # Loop through each row in the data
    control_samples <- unlist(strsplit('${meta.reference}', ","))
    treatment_samples <- unlist(strsplit('${meta.treatment}', ","))
    all_samples <- unique(c(control_samples, treatment_samples))
    name = '${meta.id}'
    design_matrix <- data.frame(matrix(0, nrow = length(all_samples), ncol = 3,
                                dimnames = list(all_samples,
                                                c("Samples", "baseline",
                                                    name))))
    # R automatically converts "-" to "." in the column names
    # so here we re-assign the column names to keep the dashes defined by the user
    colnames(design_matrix) <- c("Samples", "baseline", name)

    # Set baseline and treatment values in the design matrix
    design_matrix[, "Samples"] <- rownames(design_matrix)
    design_matrix\$baseline <- 1
    design_matrix[treatment_samples, name] <- 1

    # Print the design matrix to a file
    output_file <- paste0('${prefix}',".txt")
    write.table(design_matrix, output_file, sep = "\t", quote = FALSE, row.names=FALSE)

    """
}
