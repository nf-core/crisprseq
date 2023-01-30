#!/usr/bin/env Rscript

############################
#### Gene editing variant calling --> Parser CIGAR
#### author: Marta Sanvicente
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
        # Was the whole template aligned or do we have to look for the end of the right template?
        if (nchar(aln) >= (nchar(template_sequence) - temp_aln_start + 1) || (ref_aln_end-1 == nchar(reference_sequence)) ){ # Whole template was aligned --> corrected case template starts before reference. Or we have arrived at the end of reference
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

    return(new_reference)
}

####################################
### Load command line arguments ####
####################################

option_list = list(
    make_option(c("-r", "--reference"), type="character", default=NULL,
        help="Reference fasta file", metavar="character"),
    make_option(c("-t", "--template"), type="character", default=NULL,
        help="Temporary folder", metavar="character"),
    make_option(c("-p", "--prefix"), type="character", default=NULL,
        help="Sample prefix", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

ref_fasta = opt$reference
temp = opt$template
prefix = opt$prefix

### Get template and reference sequences from fasta files
temp_seq <- sread(readFasta(temp))[[1]]
ref_seq <- sread(readFasta(ref_fasta))[[1]]

### Generate new reference (reference with change done by the template)
NewRef <- newReference(temp_seq, ref_seq)
write.fasta(NewRef, "reference_template", file.out = paste(prefix, "NewReference.fasta", sep="_"))
