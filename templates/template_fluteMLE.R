    #!/usr/bin/env Rscript
    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    #### graphs mageck MLE

    # Required to fix corrupted cache from Singularity container
    library(BiocFileCache)
    bfc <- BiocFileCache("~/.cache/R/ExperimentHub")
    res <- bfcquery(bfc, "experimenthub.index.rds", field="rname", exact=TRUE)
    bfcremove(bfc, rids=res\$rid)
    library(ExperimentHub)
    eh = ExperimentHub()

    library(MAGeCKFlute)
    library(clusterProfiler)
    library(ggplot2)

    library(pathview)
    options(ggrepel.max.overlaps = Inf)
    mle <- read.table("${gene_summary}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    if("${prefix}" == "day0") {
        beta_strings <- grep("\\\\.beta", colnames(mle), value = TRUE)
        before_beta <- sub("\\\\.beta.*", "", beta_strings)
        unique_strings <- unique(before_beta)
        for(i in unique_strings) {
            tryCatch(
                {
                    FluteMLE(mle, treatname= i, proj=i, pathview.top=5)
                },
            error=function(e) {
                    print(paste("Could not run FluteMLE with project",i))
                }
            )
        }
    } else {
        beta_strings <- grep("\\\\.beta", colnames(mle), value = TRUE)
        before_beta <- sub("\\\\.beta.*", "", beta_strings)
        unique_strings <- unique(before_beta)
        for(i in unique_strings) {
            tryCatch(
                {
                FluteMLE(mle, treatname= i, proj=i, ${args}, pathview.top=5)
                },
            error=function(e) {
                    print(paste("Could not run FluteMLE with project",i))
                }
            )
        }
    }

    version_file_path <- "versions.yml"
    version_flute <- paste(unlist(packageVersion("MAGeCKFlute")), collapse = ".")
    version_ggplot <- paste(unlist(packageVersion("ggplot2")), collapse = ".")
    version_clusterprofiler <- paste(unlist(packageVersion("clusterProfiler")), collapse = ".")
    version_pathview <- paste(unlist(packageVersion("pathview")), collapse = ".")

    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    MAGeCKFlute: ", f, sep = "")
    writeLines(version_flute, f)
    writeLines("    ggplot2: ", f, sep = "")
    writeLines(version_ggplot, f)
    writeLines("    clusterProfiler: ", f, sep = "")
    writeLines(version_clusterprofiler, f)
    writeLines("    pathview: ", f, sep = "")
    writeLines(version_pathview, f)
    close(f)
