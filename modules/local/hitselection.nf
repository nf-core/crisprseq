// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process HITSELECTION {
    tag "$meta.id"
    label 'process_high'

    conda "r-igraph=2.0.3 r-dplyr=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7452b49b3e42a6771eca2e3ab4ecae4e0fd990ee:0':
        'biocontainers/mulled-v2-7452b49b3e42a6771eca2e3ab4ecae4e0fd990ee:0' }"

    input:
    tuple val(meta), path(gene_summary)
    path(biogrid)

    output:
    //path("overrepresented.png")   , emit: png
    path("test.csv")   , emit: png

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    #### author: Metin Yazar
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####

    library(igraph)
    library(dplyr)
    #library(ggplot2)


    interactions_df <- read.csv('${biogrid}',sep=",")
    interactions_genes <- interactions_df\$hgnc_id
    g <- graph_from_data_frame(interactions_df[, c("hgnc_id_1", "hgnc_id_2")], directed=FALSE)


    screen_genes_df <- read.csv('${gene_summary}',sep="\t")

    # Define the prefix variable and get the beta column
    prefix <- '${meta.id}'
    beta_col <- paste0(prefix, ".beta")

    # Arrange the data frame by the beta column
    print(colnames(screen_genes_df))
    genes_ordered <- screen_genes_df %>% arrange(.data[[beta_col]])
    genes_ordered <- genes_ordered[[beta_col]][1:1000]


    #valid_hit_genes <- hit_genes[hit_genes %in% V(g)\$name]

    # Calculating degree for each gene (node)
    degree <- degree(g)
    names(degree) <- V(g)\$name

    set.seed(123)
    EdgePermutationDegreeConserved <- function(permutation,frequency,degree,hit.genes,graph) {
        edges.permutation.degree.conserved <- rep(NA,permutation)
        if(length(hit.genes) > 0) {
            for(j in 1:permutation) {
                set.seed(j) # seed has to be changed with each permutation so that they are permuted differently
                hit.genes.permuted.degree.conserved <- NULL
                if(sum(frequency[which(as.numeric(names(frequency)) > 500)]) > 0){ # Check if there is any gene with degree greater than 500
                    sample.temp <- which(degree > 500)
                    hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                        names(sample(x = sample.temp,size = sum(frequency[which(as.numeric(names(frequency)) > 500)]), replace = FALSE)))
            }
            if(sum(frequency[which(as.numeric(names(frequency)) <= 500 & as.numeric(names(frequency)) > 100)]) > 0){
                sample.temp <- which(degree <= 500 & degree > 100)
                hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                        names(sample(x = sample.temp,size = sum(frequency[which(as.numeric(names(frequency)) <= 500 & as.numeric(names(frequency)) > 100)]), replace = FALSE)))
            }
            if(sum(frequency[which(as.numeric(names(frequency)) <= 100 & as.numeric(names(frequency)) > 50)]) > 0){
                sample.temp <- which(degree <= 100 & degree > 50)
                hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                        names(sample(x = sample.temp,size = sum(frequency[which(as.numeric(names(frequency)) <= 100 & as.numeric(names(frequency)) > 50)]), replace = FALSE)))
            }
            if(sum(frequency[which(as.numeric(names(frequency)) <= 50 & as.numeric(names(frequency)) > 30)]) > 0){
                sample.temp <- which(degree <= 50 & degree > 30)
                hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                        names(sample(x = sample.temp,size = sum(frequency[which(as.numeric(names(frequency)) <= 50 & as.numeric(names(frequency)) > 30)]), replace = FALSE)))
            }
            if(sum(frequency[which(as.numeric(names(frequency)) <= 30 & as.numeric(names(frequency)) > 10)]) > 0){
                sample.temp <- which(degree <= 30 & degree > 10)
                hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                    names(sample(x = sample.temp,size = sum(frequency[which(as.numeric(names(frequency)) <= 30 & as.numeric(names(frequency)) > 10)]), replace = FALSE)))
            }
            if(length(which(as.numeric(names(frequency)) <= 10)) > 0) {
                temp.frequency <- frequency[which(as.numeric(names(frequency)) <= 10)]
                for(k in 1:length(temp.frequency)){
                    sample.temp <- which(degree == as.numeric(names(temp.frequency[k])))
                    hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,names(sample(x = sample.temp,size = temp.frequency[k], replace = FALSE)))
                }
            }
            if(length(hit.genes) != length(hit.genes.permuted.degree.conserved)) {
                print("we have a problem")
            }
            #      edges.permutation.degree.conserved[j] <- length(E(induced.subgraph(graph,intersect(hit.genes.permuted.degree.conserved,V(graph)))))
            edges.permutation.degree.conserved[j] <- length(E(induced.subgraph(graph,hit.genes.permuted.degree.conserved)))
            }
    } else {
        edges.permutation.degree.conserved[1:permutation] = 0
    }
    return(edges.permutation.degree.conserved)
    }

    Find.Significance <- function(graph, hit.genes, degree, permutation){
        subGraph <- induced.subgraph(graph, hit.genes)
        edges <- length(E(subGraph))	# Compute number of edges in the network
        frequency <- table(degree[which((names(degree) %in% hit.genes) == TRUE)])  # Find degrees in the original graph
        edges.permutation.degree.conserved <- EdgePermutationDegreeConserved(permutation,frequency,degree,hit.genes,graph)
        avg_edges_permutation_degree_conserved <- mean(edges.permutation.degree.conserved)
        logp.val.degree.conserved = - ppois(edges-1, avg_edges_permutation_degree_conserved,lower.tail = FALSE,log.p = TRUE)
        no.nodes <- length(V(subGraph))
        no.edges <- length(E(subGraph))

        return(data.frame(logp.val.degree.conserved,edges,avg_edges_permutation_degree_conserved))
    }

    min <- 0
    max <- 1000
    steps <- max - min
    hit.genes.last <- NULL
    final_results <- list() # Initialize as an empty list
    for(i in 1:steps) {
        print(i)
        final_hit.genes <- intersect(interactions_genes[c(1:i)], V(g)\$name)
        hit.genes.last <- final_hit.genes
        result <- Find.Significance(g, hit.genes.last, degree, permutation = 500)
        final_results[[i]]  <- result
    }

    print(final_hit.genes)
    final_dataframe <- bind_rows(final_results)
    #final_dataframe\$gene_symbols <- gene_symbols_1000

    #write.table(final_dataframe, file = "test.tsv",sep="\t")

    #png("plot.png")
    #plot(final_dataframe\$logp.val.degree.conserved)
    #dev.off()

    version_file_path <- "versions.yml"
    version_igraph <- paste(unlist(packageVersion("igraph")), collapse = ".")
    version_dplyr <- paste(unlist(packageVersion("dplyr")), collapse = ".")

    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    igraph: ", f, sep = "")
    writeLines(version_igraph, f)
    writeLines("    dplyr: ", f, sep = "")
    writeLines(version_dplyr, f)
    close(f)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hitselection: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
