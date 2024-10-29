process HITSELECTION {
    tag "${meta.treatment}_${meta.reference}"
    label 'process_high'

    conda "r-base=4.4.1 r-igraph=2.0.3 r-dplyr=1.1.4 r-tidyr=1.3.1 r-readr=2.1.5 r-ggplot2=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-dcfb6eba6adda57b9d4a990b9096cb040320914f:588e2c290fce5c5c11ef2340b6184370efd2c628-0':
        'biocontainers/mulled-v2-dcfb6eba6adda57b9d4a990b9096cb040320914f:588e2c290fce5c5c11ef2340b6184370efd2c628-0' }"

    input:
    tuple val(meta), path(per_gene_results)
    path(biogrid)
    path(hgnc)
    val(hit_selection_iteration_nb)

    output:
    path("*_gene_conversion.txt"),        emit: gene_converted
    path("*_hitselection.tsv"),     emit: hitselection
    path("*_logpvalue_distribution.png"),         emit: plot
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

    # Function to process gene symbols from the screen data
    update_gene_columns <- function(gene_symbol) {
        gene_symbol_str <- as.character(gene_symbol)
        parts <- strsplit(gene_symbol_str, '-')[[1]]  # Split gene symbol by hyphen
        if (length(parts) > 1 && nchar(tail(parts, n=1)) >= 2) {
            # If there's a suffix with at least 2 characters, split it off as Gene_2
            return(list(paste(parts[-length(parts)], collapse = '-'), tail(parts, n = 1)))
        } else {
            # If no valid suffix, return the full symbol and NA for Gene_2
            return(list(gene_symbol_str, NA))
        }
    }

    # Function to look up gene information based on gene symbols
    lookup_gene_info <- function(gene_symbol, gene_symbol_2) {
        result <- gene_lookup[[gene_symbol]]  # Try to find information for the main symbol
        if (is.null(result)) {
            result <- gene_lookup[[gene_symbol_2]]  # Try the second symbol (suffix) if the first one fails
            if (is.null(result)) {
                # If neither symbol is found, return a "No entry is found" message
                result <- list('No entry is found', 'No entry is found')
            }
        }
        return(result)
    }

    # Function to perform degree-conserved edge permutations
    EdgePermutationDegreeConserved <- function(permutation, frequency, degree, hit.genes, graph) {
        edges.permutation.degree.conserved <- rep(NA, permutation)  # Initialize result vector
        if(length(hit.genes) > 0) {
            for(j in 1:permutation) {
                set.seed(j)  # Change the seed for each permutation
                hit.genes.permuted.degree.conserved <- NULL
                # Permute genes based on degree thresholds
                if(sum(frequency[which(as.numeric(names(frequency)) > 500)]) > 0){
                    sample.temp <- which(degree > 500)
                    hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                            names(sample(x = sample.temp, size = sum(frequency[which(as.numeric(names(frequency)) > 500)]), replace = FALSE)))
                }
                if(sum(frequency[which(as.numeric(names(frequency)) <= 500 & as.numeric(names(frequency)) > 100)]) > 0){
                    sample.temp <- which(degree <= 500 & degree > 100)
                    hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                            names(sample(x = sample.temp, size = sum(frequency[which(as.numeric(names(frequency)) <= 500 & as.numeric(names(frequency)) > 100)]), replace = FALSE)))
                }
                if(sum(frequency[which(as.numeric(names(frequency)) <= 100 & as.numeric(names(frequency)) > 50)]) > 0){
                    sample.temp <- which(degree <= 100 & degree > 50)
                    hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                            names(sample(x = sample.temp, size = sum(frequency[which(as.numeric(names(frequency)) <= 100 & as.numeric(names(frequency)) > 50)]), replace = FALSE)))
                }
                if(sum(frequency[which(as.numeric(names(frequency)) <= 50 & as.numeric(names(frequency)) > 30)]) > 0){
                    sample.temp <- which(degree <= 50 & degree > 30)
                    hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                            names(sample(x = sample.temp, size = sum(frequency[which(as.numeric(names(frequency)) <= 50 & as.numeric(names(frequency)) > 30)]), replace = FALSE)))
                }
                if(sum(frequency[which(as.numeric(names(frequency)) <= 30 & as.numeric(names(frequency)) > 10)]) > 0){
                    sample.temp <- which(degree <= 30 & degree > 10)
                    hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                            names(sample(x = sample.temp, size = sum(frequency[which(as.numeric(names(frequency)) <= 30 & as.numeric(names(frequency)) > 10)]), replace = FALSE)))
                }
                if(length(which(as.numeric(names(frequency)) <= 10)) > 0) {
                    temp.frequency <- frequency[which(as.numeric(names(frequency)) <= 10)]
                    for(k in 1:length(temp.frequency)){
                        sample.temp <- which(degree == as.numeric(names(temp.frequency[k])))
                        hit.genes.permuted.degree.conserved <- c(hit.genes.permuted.degree.conserved,
                                                                names(sample(x = sample.temp, size = temp.frequency[k], replace = FALSE)))
                    }
                }
                # Ensure the number of permuted genes matches the original hit genes
                if(length(hit.genes) != length(hit.genes.permuted.degree.conserved)) {
                    print("we have a problem")
                }
                # Count the number of edges in the permuted subgraph and store it
                edges.permutation.degree.conserved[j] <- length(E(induced.subgraph(graph, hit.genes.permuted.degree.conserved)))
            }
        } else {
            edges.permutation.degree.conserved[1:permutation] = 0  # If no hit genes, all edge counts are zero
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
        return(data.frame(logp.val.degree.conserved,edges,avg_edges_permutation_degree_conserved))
    }

    # Load necessary data
    screen <- read_delim("${per_gene_results}", delim = '\t')  # Load gene screening results

    beta_col <- grep("beta", colnames(screen), value = TRUE)
    bf_col <- grep("BF", colnames(screen), value = TRUE)
    rra_col <- grep("neg.*score", colnames(screen), value = TRUE)

    #For MAGeCK MLE output
    if(length(beta_col) >= 1) {
        screen\$Rank <- rank(screen[[beta_col]])
        screen <- screen[order(screen\$Rank), ]
    }

    #For RRA output
    if(length(rra_col) >= 1) {
        screen\$Rank <- rank(screen[[rra_col]])
        screen <- screen[order(screen\$Rank), ]
    }

    #For BAGEL2
    if(length(bf_col) >= 1) {
        screen\$Rank <- rank(screen[[bf_col]])
        screen <- screen[order(screen\$Rank), ]
    }

    hgnc <- read_delim("${hgnc}", delim = '\t')  # Load HGNC gene data

    # Select relevant columns from HGNC data
    columns_to_include <- c('hgnc_id', 'symbol', 'prev_symbol', 'ensembl_gene_id', 'alias_symbol', 'entrez_id')
    hgnc <- hgnc[columns_to_include]  # Keep only the necessary columns

    # Convert HGNC data to a list of dictionaries, where each entry corresponds to a gene
    hgnc_dict <- split(hgnc, seq(nrow(hgnc)))

    # Apply the update_gene_columns function to each gene symbol in the screen data
    updated_genes <- do.call(rbind, lapply(screen[[1]], update_gene_columns))
    screen[[1]] <- updated_genes[, 1]  # Update the main gene symbol column
    screen\$Gene_2 <- updated_genes[, 2]  # Add a new column for the suffixes (Gene_2)

    # Create a lookup dictionary for gene information based on HGNC data
    gene_lookup <- list()
    for (entry in hgnc_dict) {
        # Add main gene symbol to lookup table
        gene_lookup[[entry\$symbol]] <- list(entry\$hgnc_id, entry\$ensembl_gene_id)
        # Add previous symbols (if any) to lookup table
        if (!is.na(entry\$prev_symbol)) {
            prev_symbols <- unlist(strsplit(entry\$prev_symbol, split = "|", fixed = TRUE))
            for (prev_symbol in prev_symbols) {
                gene_lookup[[prev_symbol]] <- list(entry\$hgnc_id, entry\$ensembl_gene_id)
            }
        }
        # Add alias symbols (if any) to lookup table
        if (!is.na(entry\$alias_symbol)) {
            alias_symbols <- unlist(strsplit(entry\$alias_symbol, split = "|", fixed = TRUE))
            for (alias in alias_symbols) {
                gene_lookup[[alias]] <- list(entry\$hgnc_id, entry\$ensembl_gene_id)
            }
        }
    }


    # Apply the lookup_gene_info function to each row in the screen data
    info <- do.call(rbind, apply(screen, 1, function(row) lookup_gene_info(as.character(row[1]), as.character(row\$Gene_2))))
    info <- as.data.frame(info)

    # Add the HGNC and Ensembl gene IDs back into the screen data
    screen\$hgnc_id <- as.character(info[, 1])
    screen\$ensembl_gene_id <- as.character(info[, 2])

    if(length(beta_col) >= 1 || length(bf_col) >= 1) {
        screen\$Gene <- as.character(screen\$Gene)
    } else if(length(rra_col) >= 1) {
        screen\$id <- as.character(screen\$id)
    } else {
        screen\$GENE <- as.character(screen\$GENE)
    }

    # Save the updated screen data with HGNC and Ensembl IDs to a file
    if(length(beta_col) >= 1) {
        filename <- '${meta.treatment}_vs_${meta.reference}_mle'
    } else if(length(bf_col) >= 1) {
        filename <- '${meta.treatment}_vs_${meta.reference}_bagel2'
    } else if(length(rra_col) >= 1) {  # New else if condition
    filename <- '${meta.treatment}_vs_${meta.reference}_rra'
    } else {
        filename <- '${meta.treatment}_vs_${meta.reference}_drugz'
    }

    write_delim(as.data.frame(screen), paste0(filename, '_gene_conversion.txt'), delim = '\t')

    # Select the HGNC IDs from the screen data and the first 1000 gene symbols
    screen_genes <- screen\$hgnc_id
    gene_symbols_1000 <- screen[[1]][1:${hit_selection_iteration_nb}]
    # Load interaction data from BioGRID
    interactions_df <- read.csv("${biogrid}")

    # Create an undirected graph from the interaction data
    g <- graph_from_data_frame(interactions_df[, c("hgnc_id_1", "hgnc_id_2")], directed=FALSE)

    # Calculate the degree (number of connections) for each gene in the graph
    degree <- degree(g)
    names(degree) <- V(g)\$name

    min <- 0
    max <- ${hit_selection_iteration_nb}
    steps <- max - min
    hit.genes.last <- NULL
    final_results <- list() # Initialize as an empty list
    for(i in 1:steps) {
        final_hit.genes <- intersect(screen_genes[c(1:i)], V(g)\$name)
        hit.genes.last <- final_hit.genes
        result <- Find.Significance(g, hit.genes.last, degree, permutation = 500)
        final_results[[i]]  <- result
    }

    final_dataframe <- bind_rows(final_results)
    final_dataframe\$gene_symbols <- gene_symbols_1000
    write.table(final_dataframe, file=paste0(filename,"_hitselection.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

    max_value <- max(final_dataframe\$logp.val.degree.conserved)
    max_index <- which.max(final_dataframe\$logp.val.degree.conserved)

    # Add a column to the dataframe to indicate whether each point is the maximum or not
    final_dataframe <- final_dataframe %>%
    mutate(is_max = ifelse(logp.val.degree.conserved == max_value, "Maximum", "Other"))

    # Create the plot

    ggplot(final_dataframe, aes(x = seq_along(logp.val.degree.conserved), y = logp.val.degree.conserved)) +
        geom_point(aes(color = is_max, size = is_max)) +
        scale_color_manual(values = c("Maximum" = "red", "Other" = "black")) +
        scale_size_manual(values = c("Maximum" = 3, "Other" = 1)) +  # Adjust sizes as needed
        labs(title = "Log P-value (Degree Conserved) vs Index",
        x = "Index",
        y = "Log P-value (Degree Conserved)",
        color = "Point Type",
        size = "Point Size") +
    theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.position = "top",
        panel.grid.major = element_line(color = "gray80"),  # Major grid lines
        panel.grid.minor = element_line(color = "gray90")   # Minor grid lines
    )
    ggsave(paste0(filename, "_logpvalue_distribution.png"))

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
    touch "${meta.treatment}_vs_${meta.reference}_hitselection.tsv"
    touch "${meta.treatment}_vs_${meta.reference}_output_converted.txt"

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
}
