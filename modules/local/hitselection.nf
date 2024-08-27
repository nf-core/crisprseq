process HITSELECTION {
    tag "${meta.treatment}_${meta.reference}"
    label 'process_high'

    conda "r-base=4.4.1,r-igraph=2.0.3,r-dplyr=1.1.4,r-tidyr=1.3.1,r-readr=2.1.5,r-ggplot2=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-dcfb6eba6adda57b9d4a990b9096cb040320914f:588e2c290fce5c5c11ef2340b6184370efd2c628-0':
        'biocontainers/mulled-v2-dcfb6eba6adda57b9d4a990b9096cb040320914f:588e2c290fce5c5c11ef2340b6184370efd2c628-0' }"

    input:
    tuple val(meta), path(per_gene_results)
    path(biogrid)
    path(hgnc)

    output:
    path("*_converted.txt"),        emit: png
    path("*_hitselection.tsv"),     emit: hitselection
    path "versions.yml",            emit: versions

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
    library(ggplot2)
    library(readr)

    screen <- read_delim("${per_gene_results}",
                    delim = '\t')

    hgnc <- read_delim("${hgnc}", delim = '\t')
    columns_to_include <- c('hgnc_id', 'symbol', 'prev_symbol', 'ensembl_gene_id', 'alias_symbol', 'entrez_id')
    hgnc <- hgnc[columns_to_include]
    hgnc_dict <- split(hgnc, seq(nrow(hgnc)))
    update_gene_columns <- function(gene_symbol) {
        gene_symbol_str <- as.character(gene_symbol)
        parts <- strsplit(gene_symbol_str, '-')[[1]]
        if (length(parts) > 1 && nchar(tail(parts, n=1)) >= 2) {
            return(list(paste(parts[-length(parts)], collapse = '-'), tail(parts, n = 1)))
        } else {
            return(list(gene_symbol_str, NA))
        }
    }

    updated_genes <- do.call(rbind, lapply(screen[[1]], update_gene_columns))
    screen[[1]] <- updated_genes[, 1]
    screen\$Gene_2 <- updated_genes[, 2]


    gene_lookup <- list()
    for (entry in hgnc_dict) {
        gene_lookup[[entry\$symbol]] <- list(entry\$hgnc_id, entry\$ensembl_gene_id)

        if (!is.na(entry\$prev_symbol)) {
            prev_symbols <- unlist(strsplit(entry\$prev_symbol, split = "|", fixed = TRUE))
            for (prev_symbol in prev_symbols) {
                gene_lookup[[prev_symbol]] <- list(entry\$hgnc_id, entry\$ensembl_gene_id)
            }
        }

        if (!is.na(entry\$alias_symbol)) {
            alias_symbols <- unlist(strsplit(entry\$alias_symbol, split = "|", fixed = TRUE))
            for (alias in alias_symbols) {
                gene_lookup[[alias]] <- list(entry\$hgnc_id, entry\$ensembl_gene_id)
            }
        }
    }

    lookup_gene_info <- function(gene_symbol, gene_symbol_2) {
        result <- gene_lookup[[gene_symbol]]
        if (is.null(result)) {
            result <- gene_lookup[[gene_symbol_2]]
            if (is.null(result)) {
                result <- list('No entry is found', 'No entry is found')
            }
        }
        return(result)
    }

    info <- do.call(rbind, apply(screen, 1, function(row) lookup_gene_info(as.character(row[1]), as.character(row\$Gene_2))))
    info <- as.data.frame(info)
    screen\$hgnc_id <- as.character(info[, 1])
    screen\$ensembl_gene_id <- as.character(info[, 2])
    screen\$GENE <- as.character(screen\$GENE)

    write_delim(as.data.frame(screen), 'treatment1_vs_control1_drugz_output_converted.txt', delim = '\t')
    screen_genes <- screen\$hgnc_id
    gene_symbols_1000 <- screen\$Gene[1:1000]

    interactions_df <- read.csv("${biogrid}")

    g <- graph_from_data_frame(interactions_df[, c("hgnc_id_1", "hgnc_id_2")], directed=FALSE)

    degree <- degree(g)
    names(degree) <- V(g)\$name

    EdgePermutationDegreeConserved <- function(permutation,frequency,degree,hit.genes,graph)
    {
        edges.permutation.degree.conserved <- rep(NA,permutation)
        if(length(hit.genes) > 0)
        {
            for(j in 1:permutation)
            {
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

    Find.Significance <- function(graph, hit.genes, degree, permutation)
    {
        subGraph <- induced.subgraph(graph, hit.genes)
        edges <- length(E(subGraph))	# Compute number of edges in the network
        frequency <- table(degree[which((names(degree) %in% hit.genes) == TRUE)])  # Find degrees in the original graph
        edges.permutation.degree.conserved <- EdgePermutationDegreeConserved(permutation,frequency,degree,hit.genes,graph)
        avg_edges_permutation_degree_conserved <- mean(edges.permutation.degree.conserved)
        logp.val.degree.conserved = - ppois(edges-1, avg_edges_permutation_degree_conserved,lower.tail = FALSE,log.p = TRUE)
        #logp.val.degree.conserved <-  -log10(p.val.degree.conserved)
        no.nodes <- length(V(subGraph))
        no.edges <- length(E(subGraph))
        #  connected.components <- length(connectedComp(subGraph))
        #  clustering.coefficient <- clusteringCoef(subGraph)
        #  return(data.frame(p.val.degree.conserved,p.val.random,no.nodes,no.edges,connected.components,clustering.coefficient))

        return(data.frame(logp.val.degree.conserved,edges,avg_edges_permutation_degree_conserved))
    }

    min <- 0
    max <- 1000
    steps <- max - min
    hit.genes.last <- NULL
    final_results <- list() # Initialize as an empty list
    for(i in 1:steps) {
        print(i)
        final_hit.genes <- intersect(screen_genes[c(1:i)], V(g)\$name)
        hit.genes.last <- final_hit.genes
        result <- Find.Significance(g, hit.genes.last, degree, permutation = 500)
        final_results[[i]]  <- result
    }

    final_dataframe <- bind_rows(final_results)
    final_dataframe\$gene_symbols <- gene_symbols_1000
    write.table(final_dataframe, file="${meta.treatment}_vs_${meta.reference}_hitselection.tsv",row.names=FALSE,quote=FALSE,sep="\t")

    version_file_path <- "versions.yml"
    version_igraph <- paste(unlist(packageVersion("igraph")), collapse = ".")
    version_dplyr <- paste(unlist(packageVersion("dplyr")), collapse = ".")
    version_ggplot2 <- paste(unlist(packageVersion("ggplot2")), collapse = ".")

    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    igraph: ", f, sep = "")
    writeLines(version_igraph, f)
    writeLines("    dplyr: ", f, sep = "")
    writeLines(version_dplyr, f)
    writeLines("    ggplot2: ", f, sep = "")
    writeLines(version_ggplot2, f)
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
