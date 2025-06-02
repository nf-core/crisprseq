process VENNDIAGRAM {
    tag "${meta.treatment}_vs_${meta.reference}"
    label 'process_low'


    conda "conda-forge::r-ggvenn=0.1.10"
    container "ghcr.io/qbic-pipelines/rnadeseq:dev"

    input:
    tuple val(meta), path(bagel_drugz_genes), path(gene_summary)

    output:
    tuple val(meta), path("*.txt"), emit: common_list
    tuple val(meta), path("*.png"), emit: venn_diagram
    path "versions.yml",            emit: versions

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
    bagel_drugz_genes = read.table('$bagel_drugz_genes', sep = "\t",
                        header=TRUE)

    if (any(grepl("sumZ", colnames(bagel_drugz_genes)))) {
            filtered_precision_recall <- subset(bagel_drugz_genes, fdr_supp < 0.1)
            name_module = 'drugZ'
    } else {
            filtered_precision_recall <- subset(bagel_drugz_genes, FDR < 0.1)
            name_module = 'BAGEL2'
    }

    name <- gsub(",","_",paste0('${prefix}',".fdr"))
    filtered_mageck_mle <- mle[mle[, name] < 0.1, ]
    common_genes <- intersect(filtered_mageck_mle\$Gene,
                        filtered_precision_recall[[1]])
    data <- list(Bagel2 = filtered_precision_recall[[1]],
                        MAGeCK_MLE = filtered_mageck_mle\$Gene)

    plot_test <- ggvenn(data)
    ggsave(paste0("venn_",name_module,"_mageckmle_",name,".png"),plot_test)
    write.table(common_genes, paste0(name,'_',name_module,"_common_genes_mle.txt"),sep = "\t", quote = FALSE, row.names=FALSE)

    #version
    version_file_path <- "versions.yml"
    version_ggvenn <- paste(unlist(packageVersion("ggvenn")), collapse = ".")
    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    ggvenn: ", f, sep = "")
    writeLines(version_ggvenn, f)
    close(f)

    """


}
