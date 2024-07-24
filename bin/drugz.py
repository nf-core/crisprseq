#!/usr/bin/env python
# VERSION = "1.1.0.2"
# BUILD   = 116

# ---------------------------------
# DRUGZ:  Identify drug-gene interactions in paired sample genomic perturbation screens
# Special thanks to Matej Usaj
# Last modified 20 April 2021
# Free to modify and redistribute according to the MIT License:

# MIT License

# Copyright (c) 2018 hart-lab

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ---------------------------------


import argparse
import logging as log
# ------------------------------------
# python modules
# ------------------------------------
import sys

import numpy as np
import pandas as pd
import scipy.stats as stats
import six

log.basicConfig(level=log.INFO)
log_ = log.getLogger(__name__)


pd.options.mode.chained_assignment = None
# default='warn' - printing a warning message
# None - ignoring the warning
# "raise" - raising an exception

# ------------------------------------
# constants
norm_value = 1e7
min_reads_thresh = 1
# ------------------------------------


def load_reads(filepath, index_column, genes_to_remove):
    """
    Load the input file (raw reads counts - guide level)
    and remove guides targerting control genes
    :param filepath: The path to the file to be loaded
    :param index_column: The column to use as an index (guide IDs)
    :param genes_to_remove: A string of comma separated control gene names
    :return: reads: A dataframe containing the read counts to be used in the further analysis
    """

    # Load the input file with guide ids as row ids (file is to be tab delimited)
    reads = pd.read_csv(filepath, index_col=index_column, delimiter="\t")

    # Remove guides targeting control genes
    if genes_to_remove is not None:
        gene_column = reads.columns.values[0]
        reads_wo_remove_genes = reads.loc[~reads[gene_column].isin(genes_to_remove), :]
        return reads_wo_remove_genes
    else:
        return reads


def normalize_readcounts(reads, treatment, control):
    """
    Normalise input read counts using the global variable norm_value
    :param reads: Dataframe containing reads counts (guide level)
    :param treatment: List of columns names for the samples in the treatment group
    :param control: List of column names for the samples in the control group
    :return: normalised_counts: A dataframe containing normalised read counts
    """

    reads_to_normalize = reads[control + treatment]

    # Normalize raw read counts using norm_value (1e7)

    normalized_counts = (
        norm_value * reads_to_normalize
    ) / reads_to_normalize.sum().values
    return normalized_counts


def calculate_fold_change(
    reads, normalized_counts, control_samples, treatment_samples, pseudocount, replicate
):
    """
    Create a dataframe with index as guide ids
    Calculate log2 ratio (foldchange) between treated and control reads
    :param reads: Dataframe containing read counts (guide level)
    :param normalized_counts: Dataframe containing normalized read counts
    :param control_samples: List of control sample names
    :param treatment_samples: List of treated sample names
    :param pseudocount: Constant value added to all reads (default 5) - prevents log(0) problems
    :return: A dataframe with calculated foldchange for each replicate and
    initialized columns for guide estimated variance and z-score
    """

    fold_change = pd.DataFrame(index=reads.index.values)
    fold_change["GENE"] = reads[reads.columns.values[0]]

    # Generate foldchange, estimated_variance, and foldchange zscore column ids for each replicate
    fc_replicate_id = "fc_{replicate}".format(replicate=replicate)
    fc_zscore_id = "zscore_" + fc_replicate_id
    empirical_bayes_id = "eb_std_{replicate}".format(replicate=replicate)
    one_based_idx = replicate + 1

    # Get the control and treatment sample ids for each replicate
    control_sample = control_samples[replicate]
    treatment_sample = treatment_samples[replicate]

    # Add the control sample column to the fold change dataframe and sort by this column
    fold_change[control_sample] = reads[control_sample]
    fold_change.sort_values(control_sample, ascending=False, inplace=True)

    # Extract the number of rows (number of guides) of the reads dataframe
    no_of_guides = reads.shape[0]

    # Fill in the estimated_variance and foldchange_zscore columns with 0s for now
    fold_change[empirical_bayes_id] = np.zeros(no_of_guides)
    fold_change[fc_zscore_id] = np.zeros(no_of_guides)

    # Calculate the log2 ratio of treatment normalised read counts to control - foldchange
    fold_change[fc_replicate_id] = np.log2(
        (normalized_counts[treatment_sample] + pseudocount)
        / (normalized_counts[control_sample] + pseudocount)
    )

    return fold_change


def empirical_bayes(
    fold_change,
    half_window_size,
    no_of_guides,
    fc_replicate_id,
    empirical_bayes_id,
    fc_zscore_id,
):
    """
    Calculate the variation present in foldchange between treatment and control read counts for bins of the data, smoothing
    the variation during this process thus ensuring it increases or remains the same as the estimate for the previous bin.
    The estimate of variation is first calculated for 2xhalf_window_size guides (sorted by reads in corresponding control sample).
    The default for half_window_size is 500. So, for the first 1000 values we calculate
    the estimate of variation and the set the 1st bin (0 to half_the_window_size - i.e. 0 to 500 by default) equal to this estimate.
    For each following bin between the half-window-size and n - (where n is the number of rows (guides)) we then calculate
    the estimate of variance for this bin:
        if the estimate is greater than for the previous bin, we keep this calculated estimated variance
        otherwise (etimate is less than for the previous bin) we use the estimate variance of the previous bin.
    This smooths the variance, i.e estimated variance only ever increases or remains flat between each bin.
    The final bin (n-half_window_size : n) is then set to the variance of the previous bin (e.g. the last estimate to be calculated).
    :param fold_change: A dataframe containing foldchange (log2 ratio of treatment read counts to control read counts)
    :param half_window_size: An integer value equal to the size of the first bin and half the size of the inital sample
    (window) to estimate StDev. Default is 500.
    :param empiric_bayes_id: Column header under which to store the variance estimates.
    :param no_of_guides: An integer value euqal to the number of rows (guides) of the data frame.
    :param fc_replicate_id: Column header under which the foldchange values are stored for the current comparison (replicate)
    :return: fold_change: Updated instance of the input dataframe
    """

    # Calculate the standard deviation of foldchange based on a 2 * define window size range
    std_dev = fold_change.iloc[0 : half_window_size * 2][fc_replicate_id].std()
    fold_change[empirical_bayes_id][0:half_window_size] = std_dev

    # Iterate in a range(half_window_size, n-half_window_size, 25) where and n is the number of guides
    for i in range(half_window_size, no_of_guides - half_window_size + 25, 25):
        # in every bin calculate stdev
        std_dev = fold_change.iloc[i - half_window_size : i + half_window_size][
            fc_replicate_id
        ].std()

        # If the current variation is greater than the one for previous bin then set variation equal to this
        if std_dev >= fold_change[empirical_bayes_id][i - 1]:
            fold_change[empirical_bayes_id][
                i : i + 25
            ] = std_dev  # set new std in whole step size (25)
        # Otherwise, set it equal to the variation of the previous bin
        # This allows variation estimate for each bin to only increase or stay the same as the previous
        else:
            fold_change[empirical_bayes_id][i : i + 25] = fold_change.iloc[i - 1][
                empirical_bayes_id
            ]

    # Get the variation estimate for the final bin and set the remaining values in the empirical bayes column
    # equal to this estimate
    results = fold_change.iloc[no_of_guides - (half_window_size + 1)][
        empirical_bayes_id
    ]
    fold_change[empirical_bayes_id][no_of_guides - half_window_size :] = results

    # Calculate the z_score for each guide (fc/eb_std)
    fold_change[fc_zscore_id] = (
        fold_change[fc_replicate_id] / fold_change[empirical_bayes_id]
    )

    return fold_change, fc_zscore_id


def calculate_drugz_score(fold_change, min_observations, columns):
    """
    Calculate per gene statistics for the zscores aggregated across all comparisons
    The summed zscores and the number of observations for each gene are first aggregated. These zscores are then
    normalised and pvalue estimates (assuming guassian distribution), rank position, and FDR are calculated
    The statistics are first (with normalised zscores ranked smallest to largest) to identify synthetic
    interactions and then (with normalised zscores now ranked largest to smallest) to identify suppressor interactions
    :param fold_change: Data frame containing calculated zscores per comparison
    :param min_observations: An integer value to act as a threshold for the minimum number observations to be included
    in the analysis (default=1)
    :param columns: The header (list) of the columns containing the zscores per replicate
    :return: per_gene_results: A dataframe of summary statistics for each gene
    """

    # This is going to produce a per gene summation of the zscores for each comparison. Missing values are converted
    # to zeros. The aggregate zscore will be stored in a column named sumZ.
    # The second column will be the number of all non_zero observations (guides/gene * replicates) for that gene
    per_gene_scores = fold_change.groupby("GENE")[columns].apply(
        lambda x: pd.Series([np.nansum(x.values), np.count_nonzero(x.values)])
    )
    # Set the header for the two columns
    per_gene_scores.columns = ["sumZ", "numObs"]

    # Get a dataframe of values for genes where the number of observations is greater than the minimum threshold
    per_gene_results = per_gene_scores.loc[
        per_gene_scores.numObs >= min_observations, :
    ]

    # Update the row number (number of genes) for this new dataframe.
    no_of_genes = per_gene_results.shape[0]

    # Calcualte normalized gene z-score by:
    # 1. normalizing the sumZ values by number of observations,
    # 2. renormalizing these values to fit uniform distribution of null p-vals
    normalised_z_scores = stats.zscore(
        per_gene_results["sumZ"] / np.sqrt(per_gene_results["numObs"])
    )
    per_gene_results["normZ"] = normalised_z_scores

    # Sort the data frame by normZ (ascending) to highlight synthetic interactions
    per_gene_results.sort_values("normZ", ascending=True, inplace=True)

    # Calculate pvals (from normal dist), and fdrs (by benjamini & hochberg correction)
    per_gene_results["pval_synth"] = stats.norm.sf(per_gene_results["normZ"] * -1)
    per_gene_results["rank_synth"] = np.arange(1, no_of_genes + 1)
    scale = per_gene_results["rank_synth"] / float(no_of_genes)
    per_gene_results["fdr_synth"] = per_gene_results["pval_synth"] / scale
    per_gene_results["fdr_synth"] = np.minimum.accumulate(
        per_gene_results["fdr_synth"][::-1]
    )[::-1]

    # Resort by normZ (descending) and recalculate above values to identify suppressor interactions
    per_gene_results = per_gene_results.sort_values("normZ", ascending=False)
    per_gene_results["pval_supp"] = stats.norm.sf(per_gene_results["normZ"])
    per_gene_results["rank_supp"] = np.arange(1, no_of_genes + 1)
    scale = per_gene_results["rank_supp"] / float(no_of_genes)
    per_gene_results["fdr_supp"] = per_gene_results["pval_supp"] / scale
    per_gene_results["fdr_supp"] = np.minimum.accumulate(
        per_gene_results["fdr_supp"][::-1]
    )[::-1]

    per_gene_results = per_gene_results.sort_values("normZ", ascending=True)

    return per_gene_results


def write_drugZ_output(outfile, output):
    """Write drugZ results to a file
    :param outfile: Output file
    :param output: Per gene calculated statistics
    """
    fout = outfile
    if not hasattr(fout, "write"):
        fout = open(fout, "w")
    fout.write("GENE")
    cols = output.columns.values
    for c in cols:
        fout.write("\t" + c)
    fout.write("\n")

    for i in output.index.values:
        fout.write(
            "{0:s}\t{1:3.2f}\t{2:d}\t{3:4.2f}\t{4:.3g}\t{5:d}\t{6:.3g}\t{7:.3g}\t{8:d}\t{9:.3g}\n".format(
                i,
                output.loc[i, "sumZ"],
                int(output.loc[i, "numObs"]),
                output.loc[i, "normZ"],
                output.loc[i, "pval_synth"],
                int(output.loc[i, "rank_synth"]),
                output.loc[i, "fdr_synth"],
                output.loc[i, "pval_supp"],
                int(output.loc[i, "rank_supp"]),
                output.loc[i, "fdr_supp"],
            )
        )
    fout.close()
    return fout


def get_args():
    """Parse user giver arguments"""

    parser = argparse.ArgumentParser(
        description="DrugZ for chemogenetic interaction screens",
        epilog="dependencies: pylab, pandas",
    )
    parser._optionals.title = "Options"
    parser.add_argument(
        "-i",
        dest="infile",
        type=argparse.FileType("r"),
        metavar="sgRNA_count.txt",
        help="sgRNA readcount file",
        default=sys.stdin,
    )
    parser.add_argument(
        "-o",
        dest="drugz_output_file",
        type=argparse.FileType("w"),
        metavar="drugz-output.txt",
        help="drugz output file",
        default=sys.stdout,
    )
    parser.add_argument(
        "-f",
        dest="fc_outfile",
        type=argparse.FileType("w"),
        metavar="drugz-foldchange.txt",
        help="drugz normalized foldchange file (optional",
    )
    parser.add_argument(
        "-c",
        dest="control_samples",
        metavar="control samples",
        required=True,
        help="control samples, comma delimited",
    )
    parser.add_argument(
        "-x",
        dest="drug_samples",
        metavar="drug samples",
        required=True,
        help="treatment samples, comma delimited",
    )
    parser.add_argument(
        "-r",
        dest="remove_genes",
        metavar="remove genes",
        help="genes to remove, comma delimited",
    )
    parser.add_argument(
        "-p",
        dest="pseudocount",
        type=int,
        metavar="pseudocount",
        help="pseudocount (default=5)",
        default=5,
    )
    parser.add_argument(
        "-I",
        dest="index_column",
        type=int,
        help="Index column in the input file (default=0; GENE_CLONE column)",
        default=0,
    )
    parser.add_argument(
        "--minobs",
        dest="minObs",
        type=int,
        metavar="minObs",
        help="min number of obs (default=1)",
        default=1,
    )
    parser.add_argument(
        "--half_window_size",
        dest="half_window_size",
        type=int,
        metavar="half_window_size",
        help="width of variance-estimation window",
        default=500,
    )
    parser.add_argument(
        "-q",
        dest="quiet",
        action="store_true",
        default=False,
        help="Be quiet, do not print log messages",
    )
    parser.add_argument(
        "-unpaired",
        dest="unpaired",
        action="store_true",
        default=False,
        help="comparison status, paired (default) or unpaired",
    )
    return parser.parse_args()


def calculate_unpaired_foldchange(
    reads, normalized_counts, control_samples, treatment_samples, pseudocount
):
    """
    Calculates unpaired foldchange for each guides
    """
    fold_change = pd.DataFrame(index=reads.index.values)
    fold_change["GENE"] = reads[reads.columns.values[0]]

    # Add the control samples mean readcounts column to the fold change dataframe and sort by this column
    fold_change["ctrl_mean_reads"] = reads[control_samples].mean(axis=1)
    fold_change.sort_values("ctrl_mean_reads", ascending=False, inplace=True)

    # Add the column for unpaired foldchange (i.e. mean foldchange)
    fold_change["mean_fc"] = np.log2(
        (normalized_counts[treatment_samples].mean(axis=1) + pseudocount)
        / (normalized_counts[control_samples].mean(axis=1) + pseudocount)
    )

    # set empty columns for eb variance and zscores
    fold_change["eb_std"] = np.zeros(normalized_counts.shape[0])
    fold_change["zscore"] = np.zeros(normalized_counts.shape[0])

    return fold_change


def calculate_drugz_score_unpaired(per_gene_matrix, min_observations):
    """
    Calculate per gene statistics for the zscores aggregated across comparisons from unpaired approach
    The summed zscores and the number of observations for each gene are first aggregated. These zscores are then
    normalised and pvalue estimates (assuming guassian distribution), rank position, and FDR are calculated
    The statistics are first (with normalised zscores ranked smallest to largest) to identify synthetic
    interactions and then (with normalised zscores now ranked largest to smallest) to identify suppressor interactions
    :param per_gene_matrix: Data frame containing per gene aggregated zscores
    :param min_observations: An integer value to act as a threshold for the minimum number observations to be included
    in the analysis (default=1)
    :return: per_gene_results: A dataframe of summary statistics for each gene
    """

    per_gene_results = per_gene_matrix.loc[
        per_gene_matrix.numObs >= min_observations, :
    ]

    # Update the row number (number of genes) for this new dataframe.
    no_of_genes = per_gene_results.shape[0]

    # Calcualte normalized gene z-score by:
    # 1. normalizing the sumZ values by number of observations,
    # 2. renormalizing these values to fit uniform distribution of null p-vals
    normalised_z_scores = stats.zscore(
        per_gene_results["sumZ"] / np.sqrt(per_gene_results["numObs"])
    )
    per_gene_results["normZ"] = normalised_z_scores

    # Sort the data frame by normZ (ascending) to highlight synthetic interactions
    per_gene_results.sort_values("normZ", ascending=True, inplace=True)

    # Calculate pvals (from normal dist), and fdrs (by benjamini & hochberg correction)
    per_gene_results["pval_synth"] = stats.norm.sf(per_gene_results["normZ"] * -1)
    per_gene_results["rank_synth"] = np.arange(1, no_of_genes + 1)
    scale = per_gene_results["rank_synth"] / float(no_of_genes)
    per_gene_results["fdr_synth"] = per_gene_results["pval_synth"] / scale
    per_gene_results["fdr_synth"] = np.minimum.accumulate(
        per_gene_results["fdr_synth"][::-1]
    )[::-1]

    # Resort by normZ (descending) and recalculate above values to identify suppressor interactions
    per_gene_results = per_gene_results.sort_values("normZ", ascending=False)
    per_gene_results["pval_supp"] = stats.norm.sf(per_gene_results["normZ"])
    per_gene_results["rank_supp"] = np.arange(1, no_of_genes + 1)
    scale = per_gene_results["rank_supp"] / float(no_of_genes)
    per_gene_results["fdr_supp"] = per_gene_results["pval_supp"] / scale
    per_gene_results["fdr_supp"] = np.minimum.accumulate(
        per_gene_results["fdr_supp"][::-1]
    )[::-1]

    per_gene_results = per_gene_results.sort_values("normZ", ascending=True)

    return per_gene_results


def drugZ_analysis(args):
    """Call all functions and run drugZ analysis
    :param args: User given arguments
    """

    log_.info("Initiating analysis")

    control_samples = args.control_samples.split(",")
    treatment_samples = args.drug_samples.split(",")

    if args.remove_genes == None:
        remove_genes = []
    else:
        remove_genes = args.remove_genes.split(",")

    log_.debug("Control samples:" + str(control_samples))
    log_.debug("Treated samples:" + str(treatment_samples))

    log_.info("Loading the read count matrix")
    reads = load_reads(
        filepath=args.infile, index_column=0, genes_to_remove=remove_genes
    )
    no_of_guides = reads.shape[0]

    normalized_counts = normalize_readcounts(
        reads=reads, control=control_samples, treatment=treatment_samples
    )
    log_.info("Normalizing read counts")
    num_replicates = len(control_samples)
    fc_zscore_ids = list()
    fold_changes = list()

    if (len(control_samples) != len(treatment_samples) and args.unpaired == True) or (
        len(control_samples) == len(treatment_samples) and args.unpaired == True
    ):
        log_.info("Calculating gene-level Zscores unpaired approach")
        fold_change2 = calculate_unpaired_foldchange(
            reads,
            normalized_counts,
            control_samples=control_samples,
            treatment_samples=treatment_samples,
            pseudocount=args.pseudocount,
        )

        fold_change2 = empirical_bayes(
            fold_change=fold_change2,
            half_window_size=args.half_window_size,
            no_of_guides=no_of_guides,
            fc_replicate_id="mean_fc",
            empirical_bayes_id="eb_std",
            fc_zscore_id="zscore",
        )[0]

        per_gene_scores2 = pd.DataFrame(
            fold_change2.groupby("GENE")["zscore"].apply(
                lambda x: pd.Series([np.nansum(x.values), np.count_nonzero(x.values)])
            )
        ).unstack()
        per_gene_scores2.columns = ["sumZ", "numObs"]
        gene_normZ2 = calculate_drugz_score_unpaired(
            per_gene_matrix=per_gene_scores2, min_observations=1
        )

        log_.info("Writing output file unpaired results")
        write_drugZ_output(outfile=args.drugz_output_file, output=gene_normZ2)

    elif len(control_samples) != len(treatment_samples) and args.unpaired == False:
        log_.error(
            "Must have the same number of control and drug samples to run the paired approach"
        )

    else:
        for i in range(num_replicates):
            fold_change = calculate_fold_change(
                reads,
                normalized_counts,
                control_samples=control_samples,
                treatment_samples=treatment_samples,
                pseudocount=args.pseudocount,
                replicate=i,
            )
            log_.info("Calculating raw fold change for replicate {0}".format(i + 1))

            fold_change, fc_zscore_id = empirical_bayes(
                fold_change=fold_change,
                half_window_size=args.half_window_size,
                no_of_guides=no_of_guides,
                fc_replicate_id="fc_{replicate}".format(replicate=i),
                empirical_bayes_id="eb_std_{replicate}".format(replicate=i),
                fc_zscore_id="zscore_fc_{replicate}".format(replicate=i),
            )

            log_.info(
                "Caculating smoothed Empirical Bayes estimates of stdev for replicate {0}".format(
                    i + 1
                )
            )

            fold_changes.append(fold_change)

            log_.info("Caculating guide-level Zscores for replicate {0}".format(i + 1))
            fc_zscore_ids.append(fc_zscore_id)
            fold_change = pd.concat(fold_changes, axis=1, sort=False)
            fold_change = fold_change.loc[:, ~fold_change.columns.duplicated()]

        if args.fc_outfile:
            with args.fc_outfile as fold_change_file:
                fold_change.to_csv(fold_change_file, sep="\t", float_format="%4.3f")

        log_.info("Caculating gene-level Zscores")
        gene_normZ = calculate_drugz_score(
            fold_change=fold_change, min_observations=1, columns=fc_zscore_ids
        )

        log_.info("Writing output file paired results")
        write_drugZ_output(outfile=args.drugz_output_file, output=gene_normZ)
    if args.unpaired == True:
        return gene_normZ2
    else:
        return gene_normZ


def main():
    args = get_args()

    drugZ_analysis(args)


if __name__ == "__main__":
    main()
