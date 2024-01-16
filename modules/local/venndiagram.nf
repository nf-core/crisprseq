process VENNDIAGRAM {
   // tag "$meta.treatment"
    label 'process_low'


    conda "bioconda::r-venndiagram=1.6.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-venndiagram:1.6.16--r3.3.2_0':
        'biocontainers/r-venndiagram:1.6.16--r3.3.2_0' }"

    input:
    tuple val(meta), path(bagel_pr), path(gene_summary)

    output:
    tuple val(meta), path("*.txt"), emit: common_list
    tuple val(meta), path("*.png"), emit: venn_diagram

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.treatment}"

    """
    #!/usr/bin/env Rscript

    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    ####  produce a  venn diagram

    library(VennDiagram)
    mle = read.table('$gene_summary', sep = "\t",
                header=TRUE)
    bagel = read.table('$bagel_pr', sep = "\t",
                        header=TRUE)

    filtered_precision_recall <- subset(bagel, FDR > 0.1)
    name <- paste0('${meta.id}',".fdr")
    filtered_mageck_mle <- subset(mle,name > 0.1)
    common_genes <- intersect(filtered_mageck_mle\$Gene,
                        filtered_precision_recall\$Gene)

    write.table(common_genes, paste0('${meta.id}',"_common_genes_bagel_mle.txt"),sep = "\t", quote = FALSE, row.names=FALSE)

    venn.diagram(
        x = list(Set1 = filtered_mageck_mle\$Gene, Set2 = filtered_precision_recall\$Gene),
        category.names = c("MAGeCK MLE", "BAGEL2"),
        filename = '14_venn_diagramm.png',
        output=FALSE
        )
   # venn.plot <- venn.diagram(list(Set1 = filtered_mageck_mle\$Gene, Set2 = filtered_precision_recall\$Gene),
    #                    filename= paste0('${meta.id}',"_venn.svg"), category.names = c("MAGeCK MLE", "BAGEL2"),
     #                   # Output features
      #                  imagetype = "svg",
       #                 height = 480,
        #                width = 480,
         #               resolution = 300,
          #              compression = "lzw",

           #             # Circles
            #            lwd = 2,
             #           lty = 'blank',
              #          fill = c("#999999", "#E69F00"),

               #         # Numbers
                #        cex = .6,
                 #       fontface = "bold",
                  #      fontfamily = "sans",
                   #     cat.cex = 0.4,
                    #    cat.dist = c(-0.05, -0.02)
    #)


    """


}
