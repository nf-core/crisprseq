process VENNDIAGRAM {
    tag "${meta.treatment}_vs_${meta.reference}"
    label 'process_low'


    conda "bioconda::r-venndiagram=1.6.16"
    container "ghcr.io/qbic-pipelines/rnadeseq:dev"

    input:
    tuple val(meta), path(bagel_pr), path(gene_summary)

    output:
    tuple val(meta), path("*.txt"), emit: common_list
    tuple val(meta), path("*.png"), emit: venn_diagram

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.treatment}_vs_${meta.reference}"

    """
    #!/usr/bin/env Rscript

    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    ####  produce a  venn diagram
    library(ggvenn)
    mle = read.table('$gene_summary', sep = "\t",
                header=TRUE)
    bagel = read.table('$bagel_pr', sep = "\t",
                        header=TRUE)

    filtered_precision_recall <- subset(bagel, FDR > 0.1)
    name <- paste0('${prefix}',".fdr")
    filtered_mageck_mle <- subset(mle,name > 0.1)
    common_genes <- intersect(filtered_mageck_mle\$Gene,
                        filtered_precision_recall\$Gene)
    data <- list(Bagel2 = filtered_precision_recall\$Gene,
                        MAGeCK_MLE = filtered_mageck_mle\$Gene)

    plot_test <- ggvenn(data)
    ggsave("venn_bagel2_mageckmle.png",plot_test)
    write.table(common_genes, paste0('${prefix}',"_common_genes_bagel_mle.txt"),sep = "\t", quote = FALSE, row.names=FALSE)

    """


}
