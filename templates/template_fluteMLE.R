    #!/usr/bin/env Rscript
    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    #### graphs mageck MLE

    library(MAGeCKFlute)
    library(clusterProfiler)
    library(ggplot2)

    library(pathview)
    options(ggrepel.max.overlaps = Inf)
    mle <- read.table("${gene_summary}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    #print(head(mle))
    #stop(colnames(mle))

    if("${prefix}" == "day0") {
        print("true")
        beta_strings <- grep("\\\\.beta", colnames(mle), value = TRUE)
        before_beta <- sub("\\\\.beta.*", "", beta_strings)
        unique_strings <- unique(before_beta)
        print(unique_strings)
        for(i in unique_strings) {
            FluteMLE(mle, treatname= i, proj="${prefix}", pathview.top=0)
            }
        #print(column_names)
    } else {
        FluteMLE(mle, treatname= "${prefix}", proj="${prefix}", ${args}, pathview.top=0)
    }

    version_file_path <- "versions.yml"
    version_flute <- paste(unlist(packageVersion("MAGeCKFlute")), collapse = ".")
    version_ggplot <- paste(unlist(packageVersion("ggplot2")), collapse = ".")

    f <- file(version_file_path, "w")
    writeLines('"${task.process}":', f)
    writeLines("    MAGeCKFlute: ", f, sep = "")
    writeLines(version_flute, f)
    writeLines("    ggplot2: ", f, sep = "")
    writeLines(version_ggplot, f)
    close(f)
