#!/usr/bin/env Rscript

############################
#### Gene editing variant calling --> Parser CIGAR
############################

Sys.setenv(HOME="/root")

args = commandArgs(trailingOnly=TRUE)

library(optparse)
library(Rsamtools)
library(plyr)
library(dplyr)
library(stringr)
library(ShortRead)
library(seqinr)
library(GenomicAlignments)
library(data.table)
library(biomaRt)
library(plotly)
### HDR
library(DECIPHER)
library(Biostrings)
library(parallel)
### GenomeBrowser
library(jsonlite)

#################
### Functions ###
#################

#######
#### Empty Plots
#######
empty_plot <- function(title = NULL){
    p <- plotly_empty(type = "scatter", mode = "markers") %>%
        config(
            displayModeBar = FALSE
        ) %>%
        layout(
            title = list(
                text = title,
                yref = "paper",
                y = 0.5
            )
        )
    return(p)
}

#######
#### Reverse complement
#######
strComp=function(X){
    return(c2s(rev(comp(s2c(X)))))
}

#######
### Trim masked nucleotides at the beginning and the end of the alignment
#######
ignore_masked <- function(df_aln){
    # Ignore S and H at the beginning or end of the sequence
    new_cigar <- lapply(df_aln$cigar, function(x){
            startS <- str_replace_all(x, "^\\d+[S|H]", "")
            endS <- str_replace_all(startS, "\\d+[S|H]$", "")
            return(endS)
    })
    df_aln$cigar <- unlist(new_cigar)

    return(df_aln)
}

########
### Remove not aligned reads from the data frame
########
removeNotAlign <- function(df.summary){
    incorrect_aligns <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(incorrect_aligns) <- c("class", "cigar", "start", "length", "count", "ids")

    not_aligned <- df.summary %>% filter(is.na(cigar))
    if (dim(not_aligned)[1] != 0){
        not_align_df <- data.frame("class" = rep("not_aligned", dim(not_aligned)[1]), "cigar" = not_aligned$cigar, "start" = rep("-", dim(not_aligned)[1]), "length" = rep("-", dim(not_aligned)[1]), "count" =    not_aligned$count, "ids" = not_aligned$ids)
        incorrect_aligns <- rbind(incorrect_aligns, not_align_df)
    }
    df.summary <- df.summary %>% filter(!is.na(cigar))
    results <- list(df.summary, incorrect_aligns)
    return(results)
}

########
### Classify genome edits into insertions and deletions (analysis of the CIGAR component of the BAM file)
########
variant_call <- function(df.summary, incorrect_aligns, perrMfilter=10){
    ## Get info (length and kind of alignment -match, insertion, deletio...-) from cigar
    cigar_lengths=explodeCigarOpLengths(df.summary$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
    cigar_types=explodeCigarOps(df.summary$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M

    ## Get indels and wt sequences
    wt <- 0
    bad_alignment <- 0
    bad_wt <- 0

    indels <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(indels) <- c("cigar", "class", "start", "length", "count", "ids")

    ## Required length for filtering wt.
    all_lengths <- lapply(cigar_lengths, sum)
    main_length <- as.numeric(names(which(table(unlist(all_lengths)) == max(table(unlist(all_lengths))))))

    for (i in 1:dim(df.summary)[1]){
        if (length(cigar_types[[i]]) == 1 && cigar_types[[i]] == "M") {
            if (cigar_lengths[[i]][1] > (main_length-perrMfilter)){
                wt <- wt + df.summary$count[i]
            } else {
                bad_wt <- bad_wt + df.summary$count[i]
            }
        } else if (length(cigar_types[[i]]) != 3) {
            bad_alignment <- bad_alignment + df.summary$count[i]
            incorrect_aligns[nrow(incorrect_aligns) + 1,] = c("bad_alignment", df.summary$cigar[i], "-", "-", df.summary$count[i], df.summary$ids[i])
        } else {
            if(cigar_types[[i]][2] == "D") {
                indel_type <- "del"
            } else if (cigar_types[[i]][2] == "I") {
                indel_type <- "ins"
            } else {
                indel_type <- "unknown"
            }
            position <- df.summary$pos[i] + cigar_lengths[[i]][1] #### Alert! I have removed the -1 of the positions, since modifications starts 1 nt after the matches
            indels[nrow(indels) + 1,] = c(df.summary$cigar[i], indel_type, position, cigar_lengths[[i]][2], df.summary$count[i], df.summary$ids[i])
        }
    }
    return(list(indels, wt, bad_alignment, bad_wt, incorrect_aligns))
}

########
### Re-classify genome edits that has been previously classificated as truncated (analysis of the CIGAR component of the BAM file)
########
variant_call_truncated <- function(df.summary, ins_len_mean, mocks){
    ####
    ## Get info (length and kind of alignment -match, insertion, deletion...-) from cigar
    cigar_lengths=explodeCigarOpLengths(df.summary$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
    cigar_types=explodeCigarOps(df.summary$cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M

    ## Select which alignments have been previously classified as truncated
    pre_truncated <- which(unlist(lapply(cigar_lengths, length)) > 3)
    df.summary_trunc <- df.summary[pre_truncated,]
    cigar_lengths <- cigar_lengths[pre_truncated]
    cigar_types <- cigar_types[pre_truncated]

    ## Get delins and indels that perhaps are together woith other events that can be classified as noise
    def_bad_alignment <- 0

    delins <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(delins) <- c("cigar", "class", "start", "length", "count", "ids")

    if(dim(df.summary_trunc)[1] > 0){
        for (i in 1:dim(df.summary_trunc)[1]){
            if (length(cigar_types[[i]]) > 7) {
                def_bad_alignment <- def_bad_alignment + df.summary_trunc$count[i]
                #incorrect_aligns[nrow(incorrect_aligns) + 1,] = c("bad_alignment", df.summary$cigar[i], "-", "-", df.summary$count[i], df.summary$ids[i])
            } else if(length(cigar_types[[i]]) == 5){
                if(cigar_lengths[[i]][2] + cigar_lengths[[i]][4] > cigar_lengths[[i]][3] && cigar_lengths[[i]][3] < ins_len_mean){
                    ### Clear delins
                    position <- df.summary_trunc$pos[i] + cigar_lengths[[i]][1] #### Alert! I have removed the -1 of the positions, since modifications starts 1 nt after the matches
                    delins[nrow(delins) + 1,] = c(df.summary_trunc$cigar[i], "delin", position, cigar_lengths[[i]][2]+cigar_lengths[[i]][3]+cigar_lengths[[i]][4], df.summary_trunc$count[i], df.summary_trunc$ids[i]) # as length we are saving the length of the deletion plus the insertion
                } else if (mocks == "true") {
                    ### Big insertions to be delins for sure (are they edits + noise)
                    ######## FIRST INDEL #########
                    if (cigar_types[[i]][2] == "D"){
                        first_type <- "del"
                        del_positions <- cigar_lengths[[i]][2]
                    } else {
                        first_type <- "ins"
                        del_positions <- 0
                    }
                    position <- df.summary_trunc$pos[i] + cigar_lengths[[i]][1]
                    delins[nrow(delins) + 1,] = c(df.summary_trunc$cigar[i], first_type, position, cigar_lengths[[i]][2], df.summary_trunc$count[i], df.summary_trunc$ids[i])
                    ######## SECOND INDEL #########
                    if (cigar_types[[i]][4] == "D"){
                        second_type <- "del"
                    } else {
                        second_type <- "ins"
                    }
                    position_2 <- position + del_positions + cigar_lengths[[i]][3]
                    delins[nrow(delins) + 1,] = c(df.summary_trunc$cigar[i], second_type, position_2, cigar_lengths[[i]][4], df.summary_trunc$count[i], df.summary_trunc$ids[i])
                } else {
                    def_bad_alignment <- def_bad_alignment + df.summary_trunc$count[i]
                }
            } else if((length(cigar_types[[i]]) == 7) && (mocks == "true")){
                ######## FIRST INDEL #########
                if (cigar_types[[i]][2] == "D"){
                    first_type <- "del"
                    del_positions <- cigar_lengths[[i]][2]
                } else {
                    first_type <- "ins"
                    del_positions <- 0
                }
                position <- df.summary_trunc$pos[i] + cigar_lengths[[i]][1]
                delins[nrow(delins) + 1,] = c(df.summary_trunc$cigar[i], first_type, position, cigar_lengths[[i]][2], df.summary_trunc$count[i], df.summary_trunc$ids[i])
                ######## SECOND INDEL #########
                if (cigar_types[[i]][4] == "D"){
                    second_type <- "del"
                    del_positions_2 <- cigar_lengths[[i]][4]
                } else {
                    second_type <- "ins"
                    del_positions_2 <- 0
                }
                position_2 <- position + del_positions + cigar_lengths[[i]][3]
                delins[nrow(delins) + 1,] = c(df.summary_trunc$cigar[i], second_type, position_2, cigar_lengths[[i]][4], df.summary_trunc$count[i], df.summary_trunc$ids[i])
                ######## THIRD INDEL #########
                if (cigar_types[[i]][6] == "D"){
                    third_type <- "del"
                } else {
                    third_type <- "ins"
                }
                position_3 <- position_2 + del_positions_2 + cigar_lengths[[i]][5]
                delins[nrow(delins) + 1,] = c(df.summary_trunc$cigar[i], third_type, position_3, cigar_lengths[[i]][6], df.summary_trunc$count[i], df.summary_trunc$ids[i])
            } else if(length(cigar_types[[i]]) == 4 && cigar_lengths[[i]][2] != "M" && cigar_lengths[[i]][3] != "M")    {
                position <- df.summary_trunc$pos[i] + cigar_lengths[[i]][1]
                delins[nrow(delins) + 1,] = c(df.summary_trunc$cigar[i], "delin", position, cigar_lengths[[i]][2]+cigar_lengths[[i]][3], df.summary_trunc$count[i], df.summary_trunc$ids[i])
            } else {
                def_bad_alignment <- def_bad_alignment + df.summary_trunc$count[i]
            }
        }
    }
    return(list(delins, def_bad_alignment))
}

#######
### Deconvolute indels, one row for read
#######
indels_per_read <- function(indels, reads_rm){
    indel_reads <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(indel_reads) <- c("Modification", "Start", "Length","Ids")
    indel_reads <- apply(as.matrix(indels %>% filter(count != 0)),1,function(X) {
        split_ids <- str_split(X[6], ",")
        ### Added for spikes
        if (reads_rm == TRUE){
            spl_list <- 1:X[5]
        } else {
            spl_list <- 1:length(split_ids[[1]])
        }
        ###
        for (sing_id in spl_list){
            indel_reads[nrow(indel_reads) + 1,] = c(X[[2]], X[[3]], X[[4]], split_ids[[1]][sing_id])
        }
        return(indel_reads)
    } )
    indel_reads <- do.call(rbind, indel_reads)
    indel_reads$Start <- as.numeric(indel_reads$Start)
    indel_reads$Length <- as.numeric(indel_reads$Length)
    return(indel_reads)
}

#######
### Filter by error rate
#######
error_rate_filter <- function(indel_count, indel_reads, wt){
    ins_grouped <- indel_count %>% filter(Modification == "ins")

    ## Insertions
    good_qual_ins <- c()
    if (dim(ins_grouped)[1] > 0){
        for (i in c(1:dim(ins_grouped)[1])){
            each_ins_group <- indel_reads %>% filter(Start == as.numeric(ins_grouped[i,]$Start)) %>% filter(Length == ins_grouped[i,]$Length) %>% filter(Modification == ins_grouped[i,]$Modification)

            ## Insertions
            cigar_test <- bam[[1]]$cigar[which(bam[[1]]$qname %in% each_ins_group$Ids)]
            qual_test <- bam[[1]]$qual[which(bam[[1]]$qname %in% each_ins_group$Ids)]
            c_l <- explodeCigarOpLengths(cigar_test)
            c_t <- explodeCigarOps(cigar_test)

            pos <- lapply(1:length(c_l), function(row){
                sum(c_l[[row]][1:( which(c_t[[row]] == "I")) - 1])
            })

            phred <- lapply(1:length(pos), function(each) {
                substr(as.character(qual_test[[each]]), pos[[each]], pos[[each]])
            })

            Q <-    lapply(1:length(phred), function(q) {
                as.numeric(charToRaw(phred[[q]]))-33
            })

            P <- lapply(1:length(Q), function(p) {
                10^(-Q[[p]]/10)
            })

            good_qual_ins <- c(good_qual_ins, dim(each_ins_group)[1]/(dim(indel_reads)[1]+wt) > mean(unlist(P)))
        }
        filt_ins <- ins_grouped[good_qual_ins,]
    } else {
        filt_ins <- c()
    }

    ## Deletions
    del_grouped <- indel_count %>% filter(Modification == "del")
    if(dim(del_grouped)[1] > 0){
        good_qual_del <- c()
        for (i in c(1:dim(del_grouped)[1])){
            each_del_group <- indel_reads %>% filter(Start == as.numeric(del_grouped[i,]$Start)) %>% filter(Length == del_grouped[i,]$Length) %>% filter(Modification == del_grouped[i,]$Modification)

            ## Insertions
            cigar_test <- bam[[1]]$cigar[which(trimws(bam[[1]]$qname) %in% trimws(each_del_group$Ids))] ### some problem here...
            qual_test <- bam[[1]]$qual[which(bam[[1]]$qname %in% each_del_group$Ids)]
            c_l <- explodeCigarOpLengths(cigar_test)
            c_t <- explodeCigarOps(cigar_test)

            ### Previous position
            pos <- lapply(1:length(c_l), function(row){
                sum(c_l[[row]][1:( which(c_t[[row]] == "D")) - 1])
            })
            phred <- lapply(1:length(pos), function(each) {
                substr(as.character(qual_test[[each]]), pos[[each]], pos[[each]])
            })
            Q <-    lapply(1:length(phred), function(q) {
                as.numeric(charToRaw(phred[[q]]))-33
            })
            P <- lapply(1:length(Q), function(p) {
                10^(-Q[[p]]/10)
            })

            ### Position after the deletion
            pos_2 <- lapply(1:length(c_l), function(row){
                sum(c_l[[row]][1: which(c_t[[row]] == "D") - 1])+1
            })
            phred_2 <- lapply(1:length(pos_2), function(each) {
                substr(as.character(qual_test[[each]]), pos_2[[each]], pos_2[[each]])
            })
            Q_2 <-    lapply(1:length(phred_2), function(q) {
                as.numeric(charToRaw(phred_2[[q]]))-33
            })
            P_2 <- lapply(1:length(Q_2), function(p) {
                10^(-Q_2[[p]]/10)
            })

            if (mean(unlist(P_2)) < mean(unlist(P))){
                good_qual_del <- c(good_qual_del, dim(each_del_group)[1]/(dim(indel_reads)[1]+wt) > mean(unlist(P_2)))
            } else {
                good_qual_del <- c(good_qual_del, dim(each_del_group)[1]/(dim(indel_reads)[1]+wt) > mean(unlist(P)))
            }
        }
        filt_del <- del_grouped[good_qual_del,]

    } else {
        filt_del <- c()
    }


    filt <- rbind(filt_del, filt_ins)
    return(filt)
}


#########
### Filter by pick
#########
pick_filter <- function(indels_dataset){
    dels <- indels_dataset %>% filter(Modification == "del")
    # Count accumulated insertions and deletions by position
    accum.dels=c()
    accum.dels=unlist(apply(dels,1,function(X) return( seq(as.integer(X[2]), as.integer(X[2]) + as.integer(X[3])) )))
    # Get intervals of consecutive edited positions
    intervals_list <- list("start")
    start <- 0
    for (j in as.numeric(names(table(accum.dels)))) {
        if (start == 0){
            start <- 1
        } else {
            if(i+1 == j){
                intervals_list[length(intervals_list)][[1]] <-    c(intervals_list[length(intervals_list)][[1]], j)
            } else {
                intervals_list[length(intervals_list)+1][[1]] <-    "start"
            }
        }
        i <- j
    }
    # Get the longer interval and filter all data that is not in this interval
    longer_interval <- intervals_list[which(lengths(intervals_list) == max(lengths(intervals_list)))]
    beg <- as.numeric(longer_interval[[1]][2])
    end <- as.numeric(longer_interval[[1]][length(longer_interval[[1]])])
    indels_dataset$Start <- as.numeric(indels_dataset$Start)
    data_filt <- indels_dataset %>% filter(Start > beg) %>% filter(Start < end)
    return(data_filt)
}

#######
### Table with percentage of unique genotypes
#######
indel_count_seq <- function(indel_data){
    # Get table with percentages of indels among all genotypes
    indel_count <- plyr::count(indel_data, vars = c("Start", "Length", "Modification"))
    indel_count <- indel_count[order(indel_count$freq, decreasing = TRUE),]
    indel_count$Perc <- (indel_count$freq/sum(indel_count$freq)) * 100
    return(indel_count)
}

#######
### Get sequences of genotypes
#######
get_sequences <- function(indels, ref_seq){
    edited_sequences <- c()
    seqs <- lapply(c(1:dim(indels)[1]), function(num) {
        ref_seq_test <- as.character(sread(ref_seq))
        if (indels$Modification[[num]] == "del"){
            substr(ref_seq_test, indels$Start[[num]], indels$Start[[num]]+indels$Length[[num]]-1) <- paste0(rep("-", indels$Length[[num]]), collapse = "")
            edited_sequences <- c(edited_sequences, ref_seq_test)
        } else {
            seq <- paste(substring(ref_seq_test, 1, indels$Start[[num]]), paste0(rep("N", indels$Length[[num]]), collapse = ""), substring(ref_seq_test, indels$Start[[num]]+1, nchar(ref_seq_test)), sep = "")
            edited_sequences <- c(edited_sequences, seq)
        }
    })
}

######
### Cut site
######
get_cutSite <- function(gRNA_seq, reference, rel_cut){
    ### gRNA seq to upper case
    gRNA_seq <- toupper(gRNA_seq)
    ### Check orientation of gRNA in reference sequence and cut site
    rvComp_seq <- strComp(as.character(sread(reference)[[1]]))
    align <- pairwiseAlignment(toupper(gRNA_seq), toupper(as.character(sread(reference)[[1]])), type="global-local")
    alignRevComp <- pairwiseAlignment(toupper(gRNA_seq), toupper(rvComp_seq), type="global-local")

    if (score(align) > score(alignRevComp)){
        if (rel_cut > 0){
            cut_site <- start(subject(align))+rel_cut-1
        } else {
            cut_site <- start(subject(align)) + ( nchar(gRNA_seq) + rel_cut ) -1
        }
    } else {
        if (rel_cut > 0){
            cut_site <- nchar(as.character(sread(reference)[[1]])) - end(subject(alignRevComp)) + (nchar(gRNA_seq)-rel_cut)
        } else {
            cut_site <- nchar(as.character(sread(reference)[[1]])) - end(subject(alignRevComp)) - rel_cut
        }
    }
    return(cut_site)
}

#######
### Classes of deletions
#######
compare_kmers_left <- function(n, junction, del_sequences, i){
    ### Compare k-mer of left site of the deletion with right side of the deleted region
    remaining <- substr(junction[1], nchar(junction[1])-n+1, nchar(junction[1]))
    deleted_region <- substr(del_sequences[[i]][[2]], nchar(del_sequences[[i]][[2]])-n+1, nchar(del_sequences[[i]][[2]]))
    return(remaining == deleted_region)
}

compare_kmers_right <- function(n, junction, del_sequences, i){
    ### Compare k-mer of rigth site of the deletion with left side of the deleted region
    remaining <- substr(junction[2], 1, n)
    deleted_region <- substr(del_sequences[[i]][[2]], 1, n)
    return(remaining == deleted_region)
}

classify_deletions <- function(dels_classes, reference){
    # Get the deletions sequences
    ref_seq <- as.character(sread(reference)[[1]])
    del_sequences <- lapply(c(1:dim(dels_classes)[1]), function(line) {
        ref_deleted <- paste0(substr(x = ref_seq, start = 1, stop = dels_classes[line,]$Start-1), paste(rep("-", dels_classes[line,]$Length), collapse = ""), substr(x = ref_seq, start = dels_classes[line,]$Start+dels_classes[line,]$Length, stop = nchar(ref_seq)))
        del_nts <- substr(x = ref_seq, start = dels_classes[line,]$Start, stop = dels_classes[line,]$Start+dels_classes[line,]$Length-1)
        return(list(ref_deleted, del_nts))
    })
    # Searching MH patterns
    mmej_list <- list()
    for (i in c(1:(length(del_sequences)))) {
        len_del <- nchar(del_sequences[[i]][[2]])
        junction <- strsplit(del_sequences[[i]][[1]],    paste(rep("-", len_del), collapse = ""))[[1]]

        check_left <- TRUE
        n <- 1
        while (check_left && nchar(junction[1]) >= n) {
            check_left <- compare_kmers_left(n, junction, del_sequences, i)
            n <- n+1
        }
        n <- n-2

        check_right <- TRUE
        m <- 1
        while (check_right && nchar(junction[2]) >= m ) { ## stop when known region after the deletion is not enough long to be compared
            check_right <- compare_kmers_right(m, junction, del_sequences, i)
            m <- m+1
        }
        m <- m-2

        if (n == 0 && m == 0) {
            mmej <- "NHEJ"
        } else if (m > n){
            mmej <- substr(junction[2], 1, m)
        } else {
            mmej <- substr(junction[1], nchar(junction[1])-n+1, nchar(junction[1]))
        }
        mmej_list <- c(mmej_list,mmej)
    }
    dels_classes$patterns <- unlist(mmej_list)
    return(dels_classes)
}

#######
### Get inserted nt
#######
get_nts_insertion <- function(bam_info, separated_ins, where_nts = "ins_pos", howmuch_nts = 3){
    ### Use a bam information list with the qname and seq to get the inserted nts from the insertions separated data frame
    ### You can get the inserted nts or those that are after or before the insertion
    inserted_nt <- lapply(c(1:dim(separated_ins)[1]), function(ind){

        if (is.na(separated_indels_ins$Ids[ind])){
            return("-")
        } else {
            bam_index <- which(bam[[1]]$qname == separated_indels_ins$Ids[ind])
            cigar_lengths <- explodeCigarOpLengths(bam[[1]]$cigar[bam_index])
            cigar_types <- explodeCigarOps(bam[[1]]$cigar[bam_index])
            i_pos <- which(cigar_types[[1]] == "I")
            nts_bef_nt <- sum(cigar_lengths[[1]][1:i_pos-1])

            if (where_nts == "ins_pos"){
                inserted_nt <- substr(as.character(bam[[1]]$seq[bam_index]), nts_bef_nt +1, nts_bef_nt + separated_indels_ins$Length[ind])
            } else if (where_nts == "pre_ins_pos"){
                inserted_nt <- substr(as.character(bam_info[[1]]$seq[bam_index]), nts_bef_nt+1-howmuch_nts, nts_bef_nt)
            } else if (where_nts == "post_ins_pos"){
                inserted_nt <- substr(as.character(bam_info[[1]]$seq[bam_index]), nts_bef_nt+1+separated_indels_ins$Length[ind], nts_bef_nt+separated_indels_ins$Length[ind]+howmuch_nts)
            } else {
                return(c())
            }
            ### Return list with the nucleotides insterted or before or after the insertion
            return(inserted_nt)
        }
    })
}

correct_insertion_vc <- function(cut_site_position, insertion_row, gRNA){
    if(insertion_row$Start <= cut_site_position && insertion_row$Start + insertion_row$Length >= cut_site_position){
        # Sequences to be able to compare the insertion called by cut position
        before_cut <- substr(gRNA, 15, 17)
        after_cut <- substr(gRNA, 18, 20)
        # Get difference between the reported position and the expected by cut
        diff <- cut_site_position - insertion_row$Start + 1

        cigar_lengths <- explodeCigarOpLengths(bam[[1]]$cigar[bam[[1]]$qname == insertion_row$Ids])
        cigar_types <- explodeCigarOps(bam[[1]]$cigar[bam[[1]]$qname == insertion_row$Ids])
        i_pos <- which(cigar_types[[1]] == "I")
        nts_bef_nt <- sum(cigar_lengths[[1]][1:i_pos-1])

        before_inserted <- substr(as.character(bam[[1]]$seq[bam[[1]]$qname == insertion_row$Ids]), nts_bef_nt+1-3 + diff, nts_bef_nt + diff)
        inserted <- substr(as.character(bam[[1]]$seq[bam[[1]]$qname == insertion_row$Ids]), nts_bef_nt + 1 + diff, nts_bef_nt + insertion_row$Length + diff)
        after_inserted <- substr(as.character(bam[[1]]$seq[bam[[1]]$qname == insertion_row$Ids]), nts_bef_nt+1+insertion_row$Length + diff, nts_bef_nt+insertion_row$Length+3 + diff)

        if(before_inserted == before_cut && after_inserted == after_cut){
            return(list(Modification="ins", Start=cut_site_position+1, Length=insertion_row$Length, Ids=insertion_row$Ids, above_error_rate=insertion_row$above_error_rate, in_pick = insertion_row$in_pick, freq = insertion_row$freq, Perc=insertion_row$Perc, patterns=insertion_row$patterns, pre_ins_nt=before_inserted, ins_nt=inserted, post_ins_nt=after_inserted))
        } else {
            return(insertion_row)
        }

    } else {
        return(insertion_row)
    }
}


#######
### Get SNPs of the reference region
#######
get_SNPs <- function(fasta_coords_name, genome_ref, filter=0.0001){
    if (genome_ref == "human" || genome_ref == "mouse") {
        chr <- sapply(strsplit(as.character(fasta_coords_name@id),"\\:"), `[`, 1)
        chr_num <- sapply(strsplit(chr,"\\chr"), `[`, 2)
        coords <- sapply(strsplit(as.character(fasta_coords_name@id),"\\:"), `[`, 2)
        start <- sapply(strsplit(coords,"\\-"), `[`, 1)
        end <- sapply(strsplit(coords,"\\-"), `[`, 2)

        ### Upload the required database from biomaRt
        if (genome_ref == "human"){
            snp_dataset = "hsapiens_snp"
        } else {
            snp_dataset = "mmusculus_snp"
        }
        snp.db <- useMart("ENSEMBL_MART_SNP", dataset=snp_dataset)
        snps <- getBM(c("refsnp_id", "allele", "chr_name", "chrom_start", "chrom_end", "minor_allele_freq"),
                                    filters = c("chr_name", "start", "end"),
                                    values = list(chr_num,as.numeric(start),as.numeric(end)), mart = snp.db)
    } else {
        chr <- "chr?"
        start <- 0
        snps <- data.frame(refsnp_id=c(), allele=c(), chr_name=c(), chrom_end=c(), minor_allele_freq=c())
    }
    return(snps %>% filter(minor_allele_freq >= filter))
    # snps %>% filter(minor_allele_freq >= 0.001) #With this filter, it just appear those SNPs that USCS genome shows you
    # snps_filt <- snps %>% filter(minor_allele_freq >= 0.0001) # 1 of each 10,000
}

#########
### Get substitutions (% by position) using pileup
#########
get_substitutions_pileup <- function(bam_file){
    #### Get bam file and generate pileup
    p_param <- PileupParam(min_mapq=13,
                                                 min_base_quality=10,
                                                 min_nucleotide_depth=4)
    res <- pileup(bam_file, pileupParam = p_param)
    res <- res %>% group_by(pos, nucleotide) %>% dplyr::summarise(count = sum(count)) ### to join counts of fw with rv alignment and avoid position and nucleotides duplicates
    ### Get percentages
    fperc <- res %>% group_by(pos) %>% transmute(percentage = round((count / sum(count))*100, 2))
    fperc$nucleotide <- res$nucleotide
    ### Return data frame with all positions and the percentage of each nucleotide by position
    return(fperc)
}

#########
### Get new reference for template-based edition
#########
newReference <- function(template_sequence, reference_sequence){
    ##########
    ##### This function generates a new reference with the change of the target integrated on it
    ##### Its based in the assumption that a template is composed by two homology arms and a modification between them
    ##########

    # Align the template in both directions to check which is the correct one
    fw_alignment <- pairwiseAlignment(reference_sequence, template_sequence, type="local")
    revComp_alignment <- pairwiseAlignment(reference_sequence, reverseComplement(template_sequence), type="local")

    # Check if the alignment is better in forward or in reverse complement direction and save alignment and template in right
    if (score(fw_alignment) > score(revComp_alignment)){
        rigth_seq <- template_sequence
        aln <- fw_alignment
    } else {
        rigth_seq <- reverseComplement(template_sequence)
        aln <- revComp_alignment
    }

    # Get template start (to see if left, right or both arms have aligned) and last position of template in reference
    temp_aln_start <- start(subject(aln))
    ref_aln_end <- start(pattern(aln)) + nmatch(aln) + nmismatch(aln)    #nchar(pattern(aln))

    # Which is the longer aligned arm? Template alignment starts at position 1 (--> left arm longer; aligns first)?
    if (temp_aln_start != 1 && start(pattern(aln)) != 1){ # when left arm has not aligned. We also check if the template starts before
        left_arm <- substr(rigth_seq, 1, temp_aln_start-1)
        # The alignment starts after the change, we have to re-align the first part (short left arm) to see where in the reference starts the template
        left_start <- start(pattern(pairwiseAlignment(reference_sequence, left_arm, type="local")))
        new_reference <- paste0(substr(reference_sequence,1,left_start-1), rigth_seq, substr(reference_sequence,ref_aln_end,nchar(reference_sequence)))
        # Start and end location of the modification in the reference sequence
        change_start <- end(subject(pairwiseAlignment(reference_sequence, left_arm, type="local")))
        change_end <- temp_aln_start
    } else {
        # Has the whole template aligned or we have to look for the end of the right template?
        if (nchar(aln) >= (nchar(template_sequence) - temp_aln_start + 1) || (ref_aln_end-1 == nchar(reference_sequence)) ){ # Whole template has aligned --> corrected in case template starts before reference. Or we have arrived to the end of reference
            # Check if there is a deletion in the template. Otherwise, we have to take into account the deletion to jump its length in the reference
            if (!identical(width(insertion(aln))[[1]], integer(0))) {    ## There is a deletion in the template
                ref_aln_end <- ref_aln_end + sum(width(insertion(aln))[[1]])
                change_start <- start(insertion(aln)[[1]][1])
                change_end <- change_start + sum(width(insertion(aln))[[1]])
            } else if (!identical(width(deletion(aln))[[1]], integer(0))){ ## There is an insertion in the template
                change_start <- start(deletion(aln)[[1]][1])
                change_end <- change_start + 1 #since this is a deletion the end should be the the same position, but this causes errors later...
            } else { ## It's just a substitution
                change_start <- mismatchTable(aln)$SubjectStart[1]
                change_end <- mismatchTable(aln)$SubjectEnd[dim(mismatchTable(aln))[1]]+1
                if (is.na(change_start)){ ### This means that there is no change in the template!!
                    change_start <- 0
                    change_end <- 0
                }
            }
            new_reference <- paste0(substr(reference_sequence,1,start(pattern(aln))-1), rigth_seq, substr(reference_sequence,ref_aln_end,nchar(reference_sequence)))
        } else { # We have to look for the end of the right template
            right_arm <- substr(rigth_seq, nchar(aln)+1, length(rigth_seq))
            right_end <- end(pattern(pairwiseAlignment(reference_sequence, right_arm, type="local")))
            new_reference <- paste0(substr(reference_sequence,1,start(pattern(aln))-1), rigth_seq, substr(reference_sequence,right_end+1,nchar(reference_sequence)))
            change_start <- nchar(subject(aln))
            change_end <- nchar(aln)+1
        }
    }

    NewRef <- paste(unlist(new_reference),collapse='')
    write.fasta(new_reference,"NewRef", file.out = "NewRef.fasta")

    return(c(new_reference, change_start, change_end))
}

#############
### Template-based quantification
#############
templateCount <- function(ori_ref_file, new_ref_temp, readsAlnInfo_collaps, sampleID, alnCompleteInfo_corrected, alnCompleteInfo_NOTcorrected){
    ########
    ### Function to get number of reads supporting template-based edits. The return will be a list with two values: count and kind of substitution
    ########
    #
    #### Align new reference with the change led by template with original reference
    system(paste0("minimap2 -d ", sampleID, "_reference.mmi ", ori_ref_file))
    system(paste0("minimap2 -a -A 29 -B 17 -O 25 -E 2 ", sampleID, "_reference.mmi ", new_ref_temp, "> template.sam;"))
    #### Get cigar
    param <- ScanBamParam(
        flag=scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isDuplicate = FALSE, isSupplementaryAlignment = FALSE), ## Get just primary alignments
        what=c("cigar"))
    sam_file_path = "template.sam"
    bam=scanBam(asBam(sam_file_path), param=param)
    temp_cigar=bam[[1]]$cigar[1]
    if(is.na(temp_cigar)){
        tempCount <- c(0, "no-changes")
    } else {
        #### Cigar length and type
        cigar_lengths=explodeCigarOpLengths(temp_cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> 50,2,70
        cigar_types=explodeCigarOps(temp_cigar) ##deconvolution of CIGAR field ie: M50I2M70 --> M, I, M
        #### Remove S from cigar (it happens when one of both arms of the template are longer than the reference)
        if (cigar_types[[1]][1] == "S"){
            s_slice <- cigar_lengths[[1]][1]
            cigar_lengths[[1]] <- cigar_lengths[[1]][2:length(cigar_lengths[[1]])]
            cigar_types[[1]] <- cigar_types[[1]][2:length(cigar_types[[1]])]
            temp_cigar <- strsplit(temp_cigar, "S")[[1]][2]
            if(cigar_types[[1]][length(cigar_types[[1]])] == "S"){
                temp_cigar <- substr(temp_cigar, 1, nchar(temp_cigar) - nchar(cigar_lengths[[1]][length(cigar_lengths[[1]])]))
                cigar_lengths[[1]] <- cigar_lengths[[1]][1:length(cigar_lengths[[1]])-1]
                cigar_types[[1]] <- cigar_types[[1]][1:length(cigar_types[[1]])-1]
            }
        } else if(cigar_types[[1]][length(cigar_types[[1]])] == "S"){
            s_slice <- 0
            temp_cigar <- strsplit(temp_cigar, "S")[[1]][1]
            temp_cigar <- substr(temp_cigar, 1, nchar(temp_cigar) - nchar(cigar_lengths[[1]][length(cigar_lengths[[1]])]))
            cigar_lengths[[1]] <- cigar_lengths[[1]][1:length(cigar_lengths[[1]])-1]
            cigar_types[[1]] <- cigar_types[[1]][1:length(cigar_types[[1]])-1]
        } else {
            s_slice <- 0
        }
        ##### Check if there is any read that looks like the template!
        if(dim(readsAlnInfo_collaps %>% filter(cigar == temp_cigar))[1] == 0){
            tempCount <- c(0, "no-changes")
        } else {
            ##### Let's check which is the kind of mutation that we are facing
            if ("I" %in% cigar_types[[1]] || "D" %in% cigar_types[[1]]){
                #### When it is not a simple substitution; not all nucleotides are matches of mismatches (or clipped nucleotides!!)
                rowReads <- readsAlnInfo_collaps %>% filter(cigar == temp_cigar)
                if ("I" %in% cigar_types[[1]] && "D" %in% cigar_types[[1]]){
                    tempCount <- c(rowReads$count, "delin", strsplit(rowReads$ids, ","))
                } else if ("I" %in% cigar_types[[1]]) {
                    if(cigar_lengths[[1]][which(cigar_types[[1]] == "I")] %% 3 != 0){
                        tempCount <- c(rowReads$count, "ins-out", strsplit(rowReads$ids, ","))
                    } else {
                        tempCount <- c(rowReads$count, "ins-in", strsplit(rowReads$ids, ","))
                    }
                } else if ("D" %in% cigar_types[[1]]) {
                    if(cigar_lengths[[1]][which(cigar_types[[1]] == "D")] %% 3 != 0){
                        tempCount <- c(rowReads$count, "dels-out", strsplit(rowReads$ids, ","))
                    } else {
                        tempCount <- c(rowReads$count, "dels-in", strsplit(rowReads$ids, ","))
                    }
                }
            } else { # there is a substitution. Alignment an pattern match of the substitution in that certain region
                ref_tempSeq <- sread(readFasta(new_ref_temp))[[1]]
                ref_oriSeq <- sread(readFasta(ori_ref_file))[[1]]
                aln <- pairwiseAlignment(ref_tempSeq, ref_oriSeq)
                if (nmismatch(aln) == 0){
                    tempCount <- c(0, "no-changes")
                } else {
                    s <- mismatchTable(aln)$PatternStart[1]
                    e <- mismatchTable(aln)$PatternStart[length(mismatchTable(aln)$PatternStart)]
                    ### Filter alignments without indels (candidates to find template-based substitutions and previously characterized as wt)
                    aln_cl=explodeCigarOpLengths(alnCompleteInfo_corrected$cigar) ##### Here is important to remove S!!! alnCompleteInfo <-- corrected
                    aln_m <- which(lapply(aln_cl, function(x){length(x)==1}) == TRUE)
                    not_indels_alins <- alnCompleteInfo_NOTcorrected[aln_m,] ####### without removing S previously!!
                    ### Get the correction factor for each read in function of it S
                    correction_factors <- lapply(not_indels_alins$cigar, function(x){ if(explodeCigarOps(x)[[1]][1] == "S"){ return(explodeCigarOpLengths(x)[[1]][1]) } else { return(0) } })
                    ### Check if the directed change is found
                    diffinTemp <- substr(ref_tempSeq, s, e)
                    matchChange <- lapply(c(1:length(not_indels_alins$seq)), function(x){ substr(not_indels_alins$seq[x],s-not_indels_alins$pos[x]+1+correction_factors[[x]]-s_slice,e-not_indels_alins$pos[x]+1+correction_factors[[x]]-s_slice) == as.character(diffinTemp) })
                    temp_ids <- not_indels_alins$id[unlist(matchChange)]
                    tempCount <- c(sum(unlist(matchChange)), "subs", temp_ids)
                }
            }
        }
    }

    return(tempCount)
}

############
### Spikes correction
############
spikes_correction <- function(vc_indels){
    ########
    ### Function to get number of reads correcting by insert length biases in amplification. The proportion of over estimated reads in function of the change length will be removed from dataframe
    ########
    #
    correction_factor <- 0.156
    #
    ### Original amount of molecules = Observed count - ( Observed count * ((deletion size * slop) / (deletion size * slop + 10)) ) ## n (origin) = 10
    vc_dels <- vc_indels %>% filter(class == "del" )
    dels_count_update <- round( as.numeric(vc_dels$count) - ( as.numeric(vc_dels$count) * (as.numeric(vc_dels$length) * correction_factor /    (as.numeric(vc_dels$length) * correction_factor + 10)) ))
    vc_dels$count <- dels_count_update
    #
    ### Original amount of molecules = Observed count + insertion size * slop
    vc_ins <- vc_indels %>% filter(class == "ins" )
    ins_count_update <- round( as.numeric(vc_ins$count) + ( as.numeric(vc_ins$count) * (as.numeric(vc_ins$length) * correction_factor /    (as.numeric(vc_ins$length) * correction_factor + 10)) ))
    vc_ins$count <- ins_count_update
    #
    return(rbind(vc_dels, vc_ins))
}

####################################
### Load command line arguments ####
####################################

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
        help="Input data path", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="out.txt",
        help="Output folder path", metavar="character"),
    make_option(c("-r", "--reference"), type="character", default=NULL,
        help="Reference fasta file", metavar="character"),
    make_option(c("-g", "--gRNA_sequence"), type="character", default=NULL,
        help="gRNA sequence", metavar="character"),
    make_option(c("-n", "--sample_name"), type="character", default=NULL,
        help="Sample ID", metavar="character"),
    make_option(c("-t", "--template"), type="character", default=NULL,
        help="Temporary folder", metavar="character"),
    make_option(c("-s", "--spikes"), type="character", default=NULL,
        help="If the sample is an spikes experiment [yes or no]", metavar="character"),
    make_option(c("-f", "--summary_file"), type="character", default=NULL,
        help="Output summary file name", metavar="character"),
    make_option(c("-c", "--cut_site"), type="numeric", default=NULL,
        help="Cut position", metavar="numeric"),
    make_option(c("-m", "--mock"), type="character", default=NULL,
        help="Mock data", metavar="character"),
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data_path = opt$input
results_path = opt$output
ref_fasta = opt$reference
gRNA_sequence = opt$gRNA_sequence
sample_id = opt$sample_name
temp = opt$template
spikes = opt$spikes ### yes or no
files_summary = opt$summary_file
cut_pos_prot = as.numeric(opt$cut_site)
mock = opt$mock

#### Reference sequence
ref <- read.fasta(ref_fasta, as.string = TRUE)
ref <- toupper(as.character(ref))

#### Get substitutions from pileup
subs_plup <- get_substitutions_pileup(data_path)
write.csv(subs_plup,file=paste0(results_path, "_subs-perc.csv"))

cut_site <- get_cutSite(gRNA_sequence, readFasta(ref_fasta), cut_pos_prot)
exportJson <- toJSON(cut_site)
write(exportJson, paste(sample_id,"_cutSite.json",sep = ""))


### Get alignments data
if (substr(data_path, nchar(data_path)-2, nchar(data_path)) == "bam"){

    ## Open bam file and filter to obtain just primary alignments
    param <- ScanBamParam(
        flag=scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isDuplicate = FALSE, isSupplementaryAlignment = FALSE), ## Get just primary alignments
        what=c("seq", "pos", "cigar", "qname", "flag", "qual"))
    bam_file_path = paste0(data_path)
    bam=scanBam(bam_file_path, param=param)

    ## Create data frame
    alignment_info <- data.frame(pos=bam[[1]]$pos,cigar=bam[[1]]$cigar, id=bam[[1]]$qname, seq=bam[[1]]$seq)

} else if(substr(data_path, nchar(data_path)-2, nchar(data_path)) == "csv"){ ### For aligners that don't generate a bam file (blat and pairwise alignment)
    alignment_info <- read.csv(data_path)
} else {
    alignment_info <- data.frame()
}

#If the alignment file is empty let's create empty documents to avoid the crush of the pipeline
if (dim(alignment_info)[1] != 0){
    #### Trim masked nucleotides from alignments and remove not aligned reads
    corrected_cigar <- ignore_masked(alignment_info)
    collapsed_df <- corrected_cigar %>% group_by(cigar, pos) %>% dplyr::summarise(count = n(), ids = paste0(id, collapse = ","))
    incorrect_aln <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(incorrect_aln) <- c("class", "cigar", "start", "length", "count", "ids")
    # collapsed_df_aligned <- removeNotAlign(collapsed_df) ### These lines are not necessary when the bam file is filtered by primary alignments
    # collapsed_df <- collapsed_df_aligned[[1]]
    # incorrect_aln <- collapsed_df_aligned[[2]]

    #### Classify the alignments into indels, wt, truncated alignments and not enough long alignment to be considered wt
    vc_result <- variant_call(collapsed_df, incorrect_aln)
    indels_sure <- vc_result[[1]]
    wt_reads <- vc_result[[2]]
    trunc_reads <- vc_result[[3]] ## This is the same as --> sum(as.numeric(vc_result[[5]]$count))
    incorrect_wt <- vc_result[[4]]
    incorrect_aln <- vc_result[[5]]

    #### Check if there is a germ line variant that is "truncating" our alignments


    #### Re-classify reads with more than one indel
    all_ins <- indels_sure %>% filter(class == "ins")
    if (dim(all_ins)[1] > 0){
        ins_len_average <- mean(as.numeric(rep(all_ins$length, all_ins$count)))
    } else {
        ins_len_average <- 1
    }
    delins_noisy_indels <- variant_call_truncated(collapsed_df, ins_len_average, mock)
    #### All together
    all_indels <- rbind(delins_noisy_indels[[1]], indels_sure)

    if ( dim(all_indels)[1] > 0 ){
        ##### Spikes correction
        if(spikes == "yes"){
            ##### Remove over-represented reads and get data frame with one read in each row
            vc_result_corrected <- spikes_correction(all_indels)
            separated_indels <- indels_per_read(vc_result_corrected, TRUE)
        } else {
            ##### Get data frame with one read in each row
            separated_indels <- indels_per_read(all_indels, FALSE)
        }

        ##### Count of unique genotypes
        unique_genotypes_count <- indel_count_seq(separated_indels)

        ##### Error rate filter & pick filter
        filter_er <- error_rate_filter(unique_genotypes_count, separated_indels, wt_reads)
        separated_indels$above_error_rate <- (paste0(separated_indels$Modification,separated_indels$Start,separated_indels$Length) %in% paste0(filter_er$Modification, filter_er$Start, filter_er$Length))
        filter_pick <- pick_filter(separated_indels)
        separated_indels$in_pick <- (paste0(separated_indels$Modification,separated_indels$Start,separated_indels$Length) %in% paste0(filter_pick$Modification, filter_pick$Start, filter_pick$Length))

        ##### Deletions and insertions classes: MMEJ detection
        ref_sread <- readFasta(ref_fasta)
        ## Deletions classes: MMEJ detection
        separated_indels_dels <- separated_indels %>% filter(Modification == "del")
        unique_genotypes_count_dels <- unique_genotypes_count %>% filter(Modification == "del")
        if( dim(unique_genotypes_count_dels)[1] > 0){
            del_classes <- classify_deletions(unique_genotypes_count_dels, ref_sread)
            separated_indels_dels <- merge(x = separated_indels_dels, y = del_classes, by = c("Modification", "Start", "Length"), all = TRUE)
        }

        ## Insertions classes
        separated_indels_ins <- separated_indels %>% filter(Modification == "ins")
        if( dim(separated_indels_ins)[1] > 0){
            pre_nts <- get_nts_insertion(bam, separated_indels_ins, where_nts = "pre_ins_pos")
            separated_indels_ins$pre_ins_nt <- unlist(pre_nts)
            ins_nts <- get_nts_insertion(bam, separated_indels_ins)
            separated_indels_ins$ins_nt <- unlist(ins_nts)
            post_nts <- get_nts_insertion(bam, separated_indels_ins, where_nts = "post_ins_pos")
            separated_indels_ins$post_ins_nt <- unlist(post_nts)
            unique_genotypes_count_ins <- unique_genotypes_count %>% filter(Modification == "ins")
            separated_indels_ins_all <- merge(x = separated_indels_ins, y = unique_genotypes_count_ins, by = c("Modification", "Start", "Length"), all = TRUE)
        } else {
            separated_indels_ins_all <- data.frame("Modification"=c(),"Start"=c(),"Length"=c(),"Ids"=c(),"above_error_rate"=c(),"in_pick"=c(),"pre_ins_nt"=c(),"ins_nt"=c(),"post_ins_nt"=c(),"freq"=c(), "Perc"=c())
        }

        # Join both data frames
        if(dim(separated_indels_dels)[1]>0 && dim(separated_indels_ins_all)[1]>0){
            delins <- separated_indels %>% filter(Modification == "delin")
            separated_indels <- merge(x = separated_indels_dels, y = separated_indels_ins_all, by = c("Modification", "Start", "Length", "Ids", "above_error_rate", "in_pick", "freq", "Perc"), all = TRUE)
            if(dim(delins)[1] >0 ){
                separated_indels <- merge(x = separated_indels, y = delins, by = c("Modification", "Start", "Length", "Ids", "above_error_rate", "in_pick"), all = TRUE)
            }
        } else if(dim(separated_indels_dels)[1]>0){
            separated_indels <- separated_indels_dels
            separated_indels["pre_ins_nt"]<-NA
            separated_indels["ins_nt"]<-NA
            separated_indels["post_ins_nt"]<-NA
        }else if(dim(separated_indels_ins_all)[1]>0){
            separated_indels <- separated_indels_ins_all
        }


        ## Correct insertions by cut site position
        cut_site <- get_cutSite(gRNA_sequence, ref_sread, cut_pos_prot)
        exportJson <- toJSON(cut_site)
        write(exportJson, paste(sample_id,"_cutSite.json",sep = "")) #CREATE JSON WITH CUT SITE FOR GENOME BROWSER
        surronding_ins <- separated_indels %>% filter(Start > cut_site-4) %>% filter(Start < cut_site+4) %>% filter(Modification == "ins") %>% filter(ins_nt != "-")
        corrected_df <- data.frame()
        if(dim(surronding_ins)[1] > 0){
            for(i in c(1:dim(surronding_ins)[1])){
                corrected_df <- rbind(corrected_df, as.data.frame(correct_insertion_vc(cut_site, surronding_ins[i,], gRNA_sequence)))
            }
            separated_indels <- anti_join(separated_indels, surronding_ins, by=names(separated_indels))
            separated_indels <- rbind(separated_indels, corrected_df)
        }


        ##### Save reads with indels and stimated percentage of edition

        #Now is placed down after Save edits count

        # edition <- (dim(separated_indels)[1]/(dim(separated_indels)[1]+wt_reads))*100
        # write.csv(edition,file=paste0(results_path, "_edit-perc.csv"))

    }
    #else {
    ##### Save reads with indels and stimated percentage of edition
    # write.csv("edits-not-found",file=paste0(results_path, "_indels.csv"))
    # edition <- 0
    # write.csv(edition,file=paste0(results_path, "_edit-perc.csv"))
    #}

    ########### Template-based edition
    ### Check if there is template sequence
    if(file.size(temp) == 0 | is.na(file.size(temp))){
        temp_seq_size <- 0
    } else {
        temp_seq_size <- length(sread(readFasta(temp))[[1]])
    }
    ### If so... Let's see how many reads are template based
    if (temp_seq_size > 0){
        ### Get template and reference sequences from fasta files
        temp_seq <- sread(readFasta(temp))[[1]]
        ref_seq <- sread(readFasta(ref_fasta))[[1]]
        ### Generate new reference (reference with change done by the template)
        newRef_Results <- newReference(temp_seq, ref_seq)
        NewRef <- newRef_Results[1]
        write.fasta(NewRef, "reference_template", "newRef.fa")
        ### Look for the kind of change and get the number of reads with that change
        # cigar_aligned_reads <- names(which(table(alignment_info$cigar) == max(table(alignment_info$cigar)))) #Check if there are S to slide subs search
        # if (explodeCigarOps(cigar_aligned_reads)[[1]][1] == "S"){
        #     cigar_s <- explodeCigarOpLengths(cigar_aligned_reads)[[1]][1]
        # } else {
        #     cigar_s <- 0
        # }
        t_info <- templateCount(ref_fasta, "newRef.fa", collapsed_df, sample_id, corrected_cigar, alignment_info) ### corrected_cigar instead of alignment_info
        t_reads <- as.numeric(t_info[1])
        t_type <- t_info[2]
        write.csv(t_info[3:length(t_info)],file=paste0(results_path, "_template-reads.csv"))
    } else {
        t_reads <- 0
        t_type <- "no_template"
    }

    ##### Get SNPs from ensembl in the amplicon sequence
    # organisme_genome <- "human" ## the other option is mouse
    # SNPs <- get_SNPs(ref_sread, organisme_genome)

    ########### Edition counts
    if (dim(all_indels)[1] > 0 ) {
        # Indels
        indels_count <- dim(separated_indels)[1]
        if ( t_type == "ins-out" || t_type == "dels-out" || t_type == "ins-in" || t_type == "dels-in"){
            indels_count <- indels_count - t_reads
        }
        dels <- separated_indels %>% filter(Modification == "del")
        dels_count <- dim(dels)[1]
        if ( t_type == "dels-out" || t_type == "dels-in"){
            dels_count <- dels_count - t_reads
        }
        ins <- separated_indels %>% filter(Modification == "ins")
        ins_count <- dim(ins)[1]
        if ( t_type == "ins-out" || t_type == "ins-in"){
            ins_count <- ins_count - t_reads
        }
        # Delins
        delin <- separated_indels %>% filter(Modification == "delin")
        delin_count <- dim(delin)[1]
        if ( t_type == "delin"){
            delin_count <- delin_count - t_reads
        }
        # Indels out and in frame
        out_frame_ins <- dim(ins[ins$Length %% 3 != 0,])[1]
        if ( t_type == "ins-out"){
            out_frame_ins <- out_frame_ins - t_reads
        }
        out_frame_dels <- dim(dels[dels$Length %% 3 != 0,])[1]
        if ( t_type == "dels-out"){
            out_frame_dels <- out_frame_dels - t_reads
        }
        in_frame_ins <- dim(ins[ins$Length %% 3 == 0,])[1]
        if ( t_type == "ins-in"){
            in_frame_ins <- in_frame_ins - t_reads
        }
        in_frame_dels <- dim(dels[dels$Length %% 3 == 0,])[1]
        if ( t_type == "dels-in"){
            in_frame_dels <- in_frame_dels - t_reads
        }

        ########## Read counts
        summary_counts <- read.csv(files_summary)
        ## Reads at the uploaded files
        pre_reads <- summary_counts %>% filter(class == "raw-reads")
        reads <- pre_reads$count
        # Reads after merge
        pre_merged_reads <- summary_counts %>% filter(class == "merged-reads")
        merged_reads <- pre_merged_reads$count
        ## Reads passing quality filter
        pre_trimmed_reads <- summary_counts %>% filter(class == "quality-filtered-reads")
        trimmed_reads <- pre_trimmed_reads$count
        ## Primary aligned reads
        pre_clustered_reads <- summary_counts %>% filter(class == "clustered-reads")
        clustered_reads <- pre_clustered_reads$count
        ## Primary aligned reads
        pre_aligned_reads <- summary_counts %>% filter(class == "aligned-reads")
        aligned_reads <- pre_aligned_reads$count


        reads_classes <- c("Raw reads", "Merged reads", "Quality filtered reads", "Clustered reads", "Aligned reads")
        reads_counts <- c(as.character(reads), as.character(merged_reads), as.character(trimmed_reads), as.character(clustered_reads), as.character(aligned_reads))
        reads_summary <- data.frame(classes = unlist(reads_classes), counts = unlist(reads_counts))

        ########### Pre-variant calling counts
        # Indel filters
        ep_pre <- separated_indels %>% filter(in_pick == TRUE) %>% filter(above_error_rate == TRUE)
        ep <- dim(ep_pre)[1]
        nep_pre <- separated_indels %>% filter(in_pick == TRUE) %>% filter(above_error_rate == FALSE)
        nep <- dim(nep_pre)[1]
        nenp_pre <- separated_indels %>% filter(in_pick == FALSE) %>% filter(above_error_rate == FALSE)
        nenp <- dim(nenp_pre)[1]
        enp_pre <- separated_indels %>% filter(in_pick == FALSE) %>% filter(above_error_rate == TRUE)
        enp <- dim(enp_pre)[1]

        prevc_classes <- c("Aligned reads", "Wt", "Indels", "Wt passing filter", "Wt NOT passing filter", "Indels passing filter", "Indels NOT passing filter",
                                             "Above error & in pick", "NOT above error & in pick", "NOT above error & NOT in pick", "Above error & NOT in pick")
        prevc_counts <- c(aligned_reads, wt_reads+incorrect_wt, dim(separated_indels)[1]+trunc_reads, wt_reads, incorrect_wt, dim(separated_indels)[1], trunc_reads,
                                            ep, nep, nenp, enp)
        prevc_summary <- data.frame(classes = unlist(prevc_classes), counts = unlist(prevc_counts))

    } else {
        ###########################In case we don't have indels but maybe there are template based editions
        ########## Read counts
        indels_count = 0
        delin_count = 0
        ins_count = 0
        dels_count = 0
        in_frame_dels = 0
        out_frame_dels = 0
        in_frame_ins = 0
        out_frame_ins = 0
        columns = c("Modification", "Start", "Length", "Ids", "above_error_rate", "in_pick", "freq", "Perc", "patterns", "pre_ins_nt","ins_nt","post_ins_nt")
        separated_indels = data.frame(matrix(nrow = 1, ncol = length(columns)))
        colnames(separated_indels) = columns
        ref_sread <- readFasta(ref_fasta)
        cut_site <- get_cutSite(gRNA_sequence, ref_sread, cut_pos_prot)
        summary_counts <- read.csv(files_summary)
        ## Reads at the uploaded files
        pre_reads <- summary_counts %>% filter(class == "raw-reads")
        reads <- pre_reads$count
        # Reads after merge
        pre_merged_reads <- summary_counts %>% filter(class == "merged-reads")
        merged_reads <- pre_merged_reads$count
        ## Reads passing quality filter
        pre_trimmed_reads <- summary_counts %>% filter(class == "quality-filtered-reads")
        trimmed_reads <- pre_merged_reads$count
        ## Primary aligned reads
        pre_clustered_reads <- summary_counts %>% filter(class == "clustered-reads")
        clustered_reads <- pre_clustered_reads$count
        ## Primary aligned reads
        pre_aligned_reads <- summary_counts %>% filter(class == "aligned-reads")
        aligned_reads <- pre_aligned_reads$count

        reads_classes <- c("Raw reads", "Merged reads", "Quality filtered reads", "Clustered reads", "Aligned reads")
        reads_counts <- c(as.character(reads), as.character(merged_reads), as.character(trimmed_reads), as.character(clustered_reads), as.character(aligned_reads))
        reads_summary <- data.frame(classes = unlist(reads_classes), counts = unlist(reads_counts))

        ########### Pre-variant calling counts
        # Indel filters
        ep = 0
        nep = 0
        nenp = 0
        enp = 0

        prevc_classes <- c("Aligned reads", "Wt", "Indels", "Wt passing filter", "Wt NOT passing filter", "Indels passing filter", "Indels NOT passing filter",
                                             "Above error & in pick", "NOT above error & in pick", "NOT above error & NOT in pick", "Above error & NOT in pick")
        prevc_counts <- c(aligned_reads, wt_reads+incorrect_wt, 0+trunc_reads, wt_reads, incorrect_wt, 0, trunc_reads,
                                            ep, nep, nenp, enp)
        prevc_summary <- data.frame(classes = unlist(prevc_classes), counts = unlist(prevc_counts))
        exportJson <- toJSON(cut_site)
        write(exportJson, paste(sample_id,"_cutSite.json",sep = ""))
    }

    ###########
    ##### Summary: edition
    ###########

    edit_classes <- c("Wt", "Template-based", "Indels", "Insertions", "Deletions", "Delins", "Dels inframe", "Dels outfarme", "Ins inframe", "Ins outfarme")

    ##### Update wt if template-based is a substitution
    if ( t_type == "subs"){
        wt_reads <- wt_reads - t_reads
    }

    edit_counts <- c(wt_reads, t_reads, indels_count, ins_count, dels_count, delin_count, in_frame_dels, out_frame_dels, in_frame_ins, out_frame_ins )
    edit_summary <- data.frame(classes = unlist(edit_classes), counts = unlist(edit_counts))

    total_reads <- wt_reads+t_reads+indels_count
    wt_perc <- (wt_reads/total_reads)*100
    temp_perc <- ((t_reads)/total_reads)*100
    indels_perc <- (indels_count/total_reads)*100
    ins_perc <- (ins_count/indels_count)*100
    dels_perc <- (dels_count/indels_count)*100
    delins_perc <- (delin_count/indels_count)*100
    in_frame_ins_perc <- (in_frame_ins/ins_count)*100
    out_frame_ins_perc <- (out_frame_ins/ins_count)*100
    in_frame_dels_perc <- (in_frame_dels/dels_count)*100
    out_frame_dels_perc <- (out_frame_dels/dels_count)*100

    edit_classes_perc <- c("Wt", "Template-based", "Indels", "Delins", "Insertions", "Ins inframe", "Ins outfarme", "Deletions", "Dels inframe", "Dels outfarme")
    edit_counts_perc <- c(round(wt_perc,2), round(temp_perc,2), round(indels_perc,2),    round(delins_perc,2), round(ins_perc,2), round(in_frame_ins_perc,2), round(out_frame_ins_perc,2), round(dels_perc,2), round(in_frame_dels_perc,2), round(out_frame_dels_perc,2 ))
    edit_summary_perc <- data.frame(classes = unlist(edit_classes_perc), counts = unlist(edit_counts_perc))

    ### Save edits count
    write.csv(edit_summary_perc,file=paste0(results_path, "_edits.csv"))

    ##### Save reads with indels and stimated percentage of edition
    separated_indels$sample <- sample_id
    separated_indels$cut_site <- cut_site
    separated_indels$aligned_reads <- aligned_reads
    separated_indels$wt_reads <- wt_reads
    separated_indels$t_reads <- t_reads

    if (dim(separated_indels)[1] == 1 && is.na(separated_indels$Modification[1])) {
        write.csv(separated_indels,file=paste0(results_path, "_Emptyindels.csv"))
    } else {
        write.csv(separated_indels,file=paste0(results_path, "_indels.csv"))
    }

    ########
    ### Interactive plots
    #########
    library(plotly)

    ### Processed reads plot
    reads_summary$counts <- unlist(lapply(1:length(reads_summary$counts), function(x){ as.numeric(strsplit(as.character(reads_summary$counts[x]), " ")[[1]][2]) }))
    reads_summary$parents = c("", "Raw reads", "Merged reads", "Quality filtered reads", "Clustered reads")
    fig <- plot_ly(reads_summary,
                                 labels = ~classes,
                                 parents = ~parents,
                                 values = ~counts,
                                 type = 'sunburst',
                                 branchvalues = 'total',
                                 textinfo = "label+percent entry",
                                 textfont = list(color = '#000000', size = 20),
                                 marker = list(colors = c("#f2f2f2", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7")))

    htmlwidgets::saveWidget(as_widget(fig), paste0(results_path,"_reads.html"))

    ### Indels quality
    prevc_summary$parents <- c("", "Aligned reads", "Aligned reads", "Wt", "Wt", "Indels", "Indels", "Indels passing filter", "Indels passing filter", "Indels passing filter", "Indels passing filter")
    fig <- plot_ly(prevc_summary,
                                 labels = ~classes,
                                 parents = ~parents,
                                 values = ~counts,
                                 type = 'sunburst',
                                 branchvalues = 'total',
                                 textinfo = "label+percent entry",
                                 textfont = list(color = '#000000', size = 20),
                                 marker = list(colors = c("#f2f2f2", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7")))

    htmlwidgets::saveWidget(as_widget(fig), paste0(results_path,"_QC-indels.html"))

    ### Kinds of edits plot
    edit_summary$parents = c("", "", "", "Indels", "Indels", "Indels", "Deletions", "Deletions", "Insertions", "Insertions")
    fig <- plot_ly(edit_summary,
                                 labels = ~classes,
                                 parents = ~parents,
                                 values = ~counts,
                                 type = 'sunburst',
                                 branchvalues = 'total',
                                 textinfo = "label+percent entry",
                                 textfont = list(color = '#000000', size = 20),
                                 marker = list(colors = c("#bebebe", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7", "#9394f7")))

    htmlwidgets::saveWidget(as_widget(fig), paste0(results_path,"_edition.html"))

}else{
    fig<-empty_plot("No alignments were produced.
    Please check your files and references")
    htmlwidgets::saveWidget(as_widget(fig), paste0(results_path,"_edition.html"))
    htmlwidgets::saveWidget(as_widget(fig), paste0(results_path,"_QC-indels.html"))
    htmlwidgets::saveWidget(as_widget(fig), paste0(results_path,"_reads.html"))
    columns = c("Modification", "Start", "Length", "Ids", "above_error_rate", "in_pick", "freq", "Perc", "patterns", "pre_ins_nt", "ins_nt", "post_ins_nt", "sample", "cut_site","aligned_reads","wt_reads", "t_reads")
    separated_indels = data.frame(matrix(nrow = 1, ncol = length(columns)))
    colnames(separated_indels) = columns
    write.csv(separated_indels,file=paste0(results_path, "_Badlyindels.csv"))
    edit_classes_perc <- c("Wt", "Template-based", "Indels", "Delins", "Insertions", "Ins inframe", "Ins outfarme", "Deletions", "Dels inframe", "Dels outfarme")
    edit_counts_perc <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    edit_summary_perc <- data.frame(classes = unlist(edit_classes_perc), counts = unlist(edit_counts_perc))
    write.csv(edit_summary_perc,file=paste0(results_path, "_edits.csv"))
}
