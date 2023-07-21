#!/usr/bin/env python

import click
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.linear_model import LinearRegression
import sys
import time

VERSION = 2.0

BUILD = 115



'''
Update history

Build 115
1. Add single pass without resampling

Build 114
1. Add option for normalizing rep counts

Build 113
1. Fixed bugs

Build 112
1. Add sgRNA filtering options
2. Implemented 'Click' library. Thanks to John McGonigle


Build 111
1. Add an option to equalize # of sgRNA per gene

Build 110
1. Enable multi control for fold-change calculation
2. Now user can input column names
3. Fix Threshold function

'''


class OptionRequiredIf(click.Option):

    def full_process_value(self, ctx, value):
        value = super(OptionRequiredIf, self).full_process_value(ctx, value)

        if value is None and ctx.params['filter_multi_target'] is True:
            msg = 'Error! Multi-target-filtering selected and not align-info provided.\n' \
                  'Please indicate align-info file.'
            raise click.MissingParameter(ctx=ctx, param=self, message=msg)
        return value


# ---------------------------------
# BAGEL:  Bayesian Analysis of Gene EssentaLity
# (c) Traver Hart <traver@hart-lab.org>, Eiru Kim <rooeikim@gmail.com> 2017.

# Acknowledgements: John McGonigle <j.e.mcgonigle@gmail.com>
# modified 10/2019
# Free to modify and redistribute with attribution
# ---------------------------------

# ------------------------------------
# constants



CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def round_to_hundredth(x):
    return np.around(x * 100) / 100.0


def func_linear(x, a, b):
    return (a * x) + b


class Training:
    def __init__(self, X, n=None, cvnum=10):

        if n == None:
            self._n = len(X)
        self._cvnum = cvnum
        self._bid = int(self._n / cvnum)
        self._bucket = np.arange(len(X))
        self._X = X
        self._step = 0


    def cross_validation(self):
        if self._bid < 1:  # bid check
            print("The number of genes is too small! n<" + str(self._cvnum))
            sys.exit(1)
        drawing = list()
        mask = np.array([True] * self._n)
        for j in range(self._bid):
            # drawing.append(delete(self._bucket, np.random.randrange(len(self._bucket))))
            select = np.random.randint(len(self._bucket))
            drawing.append(self._bucket[select])
            mask[self._bucket[select]] = False
            self._bucket = np.delete(self._bucket, select)
        if self._step < self._n % self._cvnum:  # for distribute remain..
            select = np.random.randint(len(self._bucket))
            drawing.append(self._bucket[select])
            mask[self._bucket[select]] = False
            self._bucket = np.delete(self._bucket, select)
        self._step += 1
        X_resample = self._X[mask]
        return X_resample, self._X[~mask]

    def get_cv_step(self):
        return self._step

    def bootstrap_resample(self):
        mask = np.array([False] * self._n)
        resample_i = np.floor(np.random.rand(self._n) * len(self._X)).astype(int)

        mask[resample_i] = True
        X_resample = self._X[mask]
        return X_resample, self._X[~mask]

    def get_data(self, method=0):
        if method == 0:
            train, test = self.bootstrap_resample()
        elif method == 1:
            train, test = self.cross_validation()
        else:
            print('Method passed as a value that was neither 0 nor 1.')
            sys.exit(1)
        return train, test


def fibo_weighted_sum(listofscore):
    value = p1 = p2 = 0.0
    c = 1.0  # current value
    for v in listofscore:
        value += v / c
        p2 = p1  # go one step
        p1 = c
        c = p1 + p2
    return value


# Note the \b in the doc string below is prevent click from wrapping the lines on the terminal.
@click.group(context_settings=CONTEXT_SETTINGS)
def bagel():
    """
    --------------------------------------------------------------------
    BAGEL.py
    --------------------------------------------------------------------
    A tool from the Bayesian Analysis of Gene EssentiaLity (BAGEL) suite.

    \b
    Calculate fold changes from read count data:

        \b
        BAGEL.py fc -i [read count file] -o [output label] -c [control column]

    Calculate Bayes Factors from foldchange data:

        \b
        BAGEL.py bf -i [fold change] -o [output file] -e [essentials genes] -n [nonessentials genes] -c [columns]


    Calculate precision-recall from Bayes Factors:

        \b
        BAGEL.py pr -i [Bayes Factor file] -o [output file] -e [essentials genes] -n [nonessentials genes]



    To print the current build and version use:

        \b
        BAGEL.py version
    """


@click.command(name='version')
def report_bagel_version():
    """
    Report the current build and version no.
    """
    print(
        'Bayesian Analysis of Gene EssentiaLity (BAGEL) suite:\n'
        'Version: {VERSION}\n'
        'Build: {BUILD}'.format(VERSION=VERSION, BUILD=BUILD)
    )


@click.command(name='fc')
@click.option('-i', '--read-count-file', required=True, type=click.Path(exists=True))
@click.option('-o', '--output-label', required=True)
@click.option('-c', '--control-columns', required=True)
@click.option('-m', '--min-reads', type=int, default=0)
@click.option('-Np', '--pseudo-count', type=int, default=5)
def calculate_fold_change(read_count_file, output_label, control_columns, min_reads, pseudo_count):
    """
    \b
    Calculate fold changes from read count data outputting a fold change column:

        \b
        BAGEL.py fc -i [read count file] -o [output label] -c [control column]

    \b
    Required options:
        -i --read-count-file       Tab-delimited file of reagents and fold changes.  See documentation for format.
        -o --output-label          Label for all output files
        -c --control-columns       A comma-delimited list of columns of control (T0 or plasmid) columns.
                                    Input can be either number or name.
    \b
    Other options:
        --min-reads=N              Discard gRNA with T0 counts < N (default 0)
        -Np, --pseudo-count=N	     Add a pseudocount of N to every readcount (default 5)
        -h, --help                 Show this help text


    \b
    Example:
        BAGEL.py fc -i readcount_file.txt -o experiment_name -c 1

    This command calculates fold change, and writes [output label].foldchange and [output label].normalized_reads

    """

    # ---------------------------------------------------------------- #
    #  Import raw read data, normalize, filter for T0 min read counts  #
    #  Output:   [output label].foldchange                             #
    # ---------------------------------------------------------------- #
    reads = pd.read_csv(read_count_file, sep='\t', index_col=0)

    reads[reads.columns.values[0]].fillna('NO_GENE_NAME', inplace=True)
    reads.fillna(0, inplace=True)
    control_columns = control_columns.split(",")
    #
    # check if controls are given as numeric or name
    #

    try:
        try:
            ctrl_columns = list(map(int, control_columns))
            ctrl_labels = reads.columns.values[ctrl_columns]
        except ValueError:
            ctrl_labels = control_columns

        ctrl_sum = reads[ctrl_labels].sum(axis=1)
        reads.drop(ctrl_labels, axis=1, inplace=True)
        ctrl_label_new = ';'.join(ctrl_labels)
        reads[ctrl_label_new] = ctrl_sum
    except:
        print(reads[ctrl_labels].sum(axis=1))
        print("Invalid input controls")
        sys.exit(1)

    numClones, numColumns = reads.shape
    print("Controls: " + ", ".join(ctrl_labels))

    #
    # Add pseudo count
    #

    reads.iloc[:, list(range(1, numColumns))] += pseudo_count

    #
    # normalize each sample to a fixed total readcount
    #
    sumReads = reads.iloc[:, list(range(1, numColumns))].sum(0)
    normed = pd.DataFrame(index=reads.index.values)
    normed['GENE'] = reads.iloc[:, 0]  # first column is gene name
    normed = reads.iloc[:, list(range(1, numColumns))] / np.tile(sumReads,
                                                                 [numClones, 1]) * 10000000  # normalize to 10M reads

    #
    # filter for minimum readcount
    #
    f = np.where(reads[ctrl_label_new] >= min_reads)[0]
    normed = normed.iloc[f, :]

    #
    # calculate fold change
    #
    foldchange = pd.DataFrame(index=normed.index.values)
    foldchange.index.name = 'REAGENT_ID'
    foldchange['GENE'] = reads.iloc[f, 0]  # dataframe 'normed' has no GENE column
    for i in range(numColumns - 1):
        foldchange[normed.columns.values[i]] = np.log2(
            (normed.loc[:, normed.columns.values[i]]) / normed[ctrl_label_new]
        )
    #
    # we have calculated a foldchange for the control column.  Drop it.
    #
    foldchange.drop(ctrl_label_new, axis=1, inplace=True)

    #
    # write normed readcount file
    # write foldchange file
    #

    foldchange_filename = output_label + '.foldchange'
    foldchange.to_csv(foldchange_filename, sep='\t', float_format='%4.3f')

    normedreads_filename = output_label + '.normed_readcount'
    normed.to_csv(normedreads_filename, sep='\t', float_format='%3.2f')


@click.command(name='bf')
@click.option('-i', '--fold-change', required=True, type=click.Path(exists=True))
@click.option('-o', '--output-file', required=True)
@click.option('-e', '--essential-genes', required=True, type=click.Path(exists=True))
@click.option('-n', '--non-essential-genes', required=True, type=click.Path(exists=True))
@click.option('-c', '--columns-to-test', required=True)
@click.option('-w', '--network-file', metavar='[network File]', default=None, type=click.Path(exists=True))
@click.option('-m', '--filter-multi-target', is_flag=True)
@click.option('-m0', '--loci-without-mismatch', type=int, default=10)
@click.option('-m1', '--loci-with-mismatch', type=int, default=10)
@click.option('--align-info', metavar='--align-info [File]', default=None,
              type=click.Path(exists=True), cls=OptionRequiredIf)
@click.option('-b', '--use-bootstrapping', is_flag=True)
@click.option('-NS', '--no-resampling', is_flag=True)
@click.option('-s', '--use-small-sample', is_flag=True)
@click.option('-N', '--no-of-cross-validations', type=int, default=10)
@click.option('-NB', '--bootstrap-iterations', type=int, default=1000)
@click.option('-r', '--sgrna-bayes-factors', is_flag=True)
@click.option('-f', '--equalise-sgrna-no', type=int)
@click.option('-s', '--seed', default=int(time.time() * 100000 % 100000), type=int)
@click.option('-t', '--run-test-mode', is_flag=True)
@click.option('-p', '--equalise-rep-no', type=int)
def calculate_bayes_factors(
        fold_change, output_file, essential_genes, non_essential_genes, columns_to_test, network_file, align_info,
        use_bootstrapping, no_resampling, use_small_sample, filter_multi_target, loci_without_mismatch, loci_with_mismatch,
        bootstrap_iterations, no_of_cross_validations, sgrna_bayes_factors, equalise_sgrna_no, seed, run_test_mode, equalise_rep_no
):
    """
    \b
    Calculate Bayes Factors from an input fold change file:

        \b
        BAGEL.py bf -i [fold change] -o [output file] -e [essentials genes] -n [nonessentials genes] -c [columns]


    \b
    Calculates a log2 Bayes Factor for each gene. Positive BFs indicate confidence that the gene is essential.
    Output written to the [output file] contains: gene name, mean Bayes Factor across all iterations, std deviation of
    BFs, and number of iterations in which the gene was part of the test set (and a BF was calculated[output file].


    \b
    Required options:
        -i --fold-change [fold change file]                     Tab-delimited file of reagents and fold changes
                                                                (see documentation for format).
        -o, --output-file [output file]                         Output filename
        -e, --essential-genes [reference essentials]            File with list of training set of essential genes
        -n, --non-essential-genes [reference nonessentials]     File with list of training set of nonessential genes
        -c [columns to test]                                    comma-delimited list of columns in input file to
                                                                include in analyisis

    \b
    Network options:
        -w  [network file]    Enable Network boosting. Tab-delmited file of edges. [GeneA (\\t) GeneB]\n'

    \b
    Multi-target guides filtering options:
        -m, --filter-multi-target     Enable filtering multi-targeting guide RNAs
        --align-info  [file]          Input precalculated align-info file
        -m0, --loci-without-mismatch  Filtering guide RNAs without mismatch targeting over than [N] loci, default = 10
        -m1, --loci-with-mismatch     Filtering guide RNAs with 1-bp mismatch targeting over than [N] loci, default = 10

    \b
    Other options:
        -b, --bootstrapping            Use bootstrapping instead of cross-validation (Slow)
        -NS, --no-resampling           Run BAGEL without resampling
        -s, --small-sample             Low-fat BAGEL, Only resampled training set (Bootstrapping, iteration = 100)
        -r  --sgrna-bayes-factors      Calculate sgRNA-wise Bayes Factor
        -f  --equalise-sgrna-no        Equalize the number of sgRNAs per gene to particular value [Number]
        -p  --equalise-rep-no          Equalize the number of repicates to particular value [Number]
        -N  --no-of-cross-validations  Number of sections for cross validation (default 10)
        -NB  --bootstraps-iterations   Number of bootstrap iterations (default 1000)
        -s, --seed=N                   Define random seed
        -h, --help                     Show this help text

    \b
    Example:

        \b
        BAGEL.py bf -i fc_file.txt -o results.bf -e ess_training_set.txt -n noness_training_set.txt -c 1,2,3

    """
    np.random.seed(seed)  # set random seed
    if network_file:
        network_boost = True
    else:
        network_boost = False

    if sgrna_bayes_factors:
        rna_level = True
    else:
        rna_level = False

    if network_file and sgrna_bayes_factors:
        network_boost = False

    if equalise_sgrna_no:
        flat_sgrna = True
    else:
        flat_sgrna = False

    if equalise_rep_no:
        flat_rep = True
    else:
        flat_rep = False

    if use_small_sample:
        train_method = 0
        bootstrap_iterations = 100

    elif use_bootstrapping:
        train_method = 0

    else:
        train_method = 1

    genes = {}
    fc = {}
    gene2rna = {}
    rna2gene = {}

    multi_targeting_sgrnas = dict()
    multi_targeting_sgrnas_info = dict()

    if filter_multi_target:

        try:
            aligninfo = pd.read_csv(align_info, header=None, index_col=0, sep="\t").fillna("")
            for seqid in aligninfo.index:
                perfectmatch = 0
                mismatch_1bp = 0
                perfectmatch_gene = 0
                mismatch_1bp_gene = 0
                if aligninfo[1][seqid] != "":
                    perfectmatch = len(aligninfo[1][seqid].split(","))
                if aligninfo[2][seqid] != "":
                    perfectmatch_gene = len(aligninfo[2][seqid].split(","))
                if aligninfo[3][seqid] != "":
                    mismatch_1bp = len(aligninfo[3][seqid].split(","))
                if aligninfo[4][seqid] != "":
                    mismatch_1bp_gene = len(aligninfo[4][seqid].split(","))
                if perfectmatch > loci_without_mismatch or mismatch_1bp > loci_with_mismatch:
                    multi_targeting_sgrnas[seqid] = True
                elif perfectmatch > 1 or mismatch_1bp > 0:
                    multi_targeting_sgrnas_info[seqid] = (
                        perfectmatch, mismatch_1bp, perfectmatch_gene, mismatch_1bp_gene
                    )

        except:
            print("Please check align-info file")
            sys.exit(1)

        print("Total %d multi-targeting gRNAs are discarded" % len(multi_targeting_sgrnas))

    #
    # LOAD FOLDCHANGES
    #
    rnatagset = set()
    with open(fold_change) as fin:
        fieldname = fin.readline().rstrip().split('\t')
        #
        # DEFINE CONTROLS
        #
        columns = columns_to_test.split(',')
        try:
            try:
                column_list = list(map(int, columns))
                column_labels = [fieldname[x + 1] for x in column_list]
            except ValueError:
                column_labels = columns
                column_list = [x for x in range(len(fieldname) - 1) if
                               fieldname[x + 1] in column_labels]  # +1 because of First column start 2
            print("Using column:  " + ", ".join(column_labels))
        # print "Using column:  " + ", ".join(map(str,column_list))

        except:
            print("Invalid columns")
            sys.exit(1)

        for line in fin:
            fields = line.rstrip().split('\t')
            rnatag = fields[0]
            if filter_multi_target is True:  # multitargeting sgrna filtering
                if rnatag in multi_targeting_sgrnas:
                    continue  # skip multitargeting sgrna.
            if rnatag in rnatagset:
                print("Error! sgRNA tag duplicates")
                sys.exit(1)
            rnatagset.add(rnatag)
            gsym = fields[1]

            genes[gsym] = 1
            if gsym not in gene2rna:
                gene2rna[gsym] = []
            gene2rna[gsym].append(rnatag)
            rna2gene[rnatag] = gsym
            fc[rnatag] = {}
            for i in column_list:
                fc[rnatag][i] = float(fields[i + 1])  # per user docs, GENE is column 0, first data column is col 1.

    genes_array = np.array(list(genes.keys()))
    gene_idx = np.arange(len(genes))
    print("Number of unique genes:  " + str(len(genes)))

    #
    # DEFINE REFERENCE SETS
    #
    coreEss = []

    with open(essential_genes) as fin:
        skip_header = fin.readline()
        for line in fin:
            coreEss.append(line.rstrip().split('\t')[0])
    coreEss = np.array(coreEss)
    print("Number of reference essentials: " + str(len(coreEss)))

    nonEss = []
    with open(non_essential_genes) as fin:
        skip_header = fin.readline()
        for line in fin:
            nonEss.append(line.rstrip().split('\t')[0])

    nonEss = np.array(nonEss)
    print("Number of reference nonessentials: " + str(len(nonEss)))

    #
    # LOAD NETWORK
    #

    if network_boost is True:
        network = {}
        edgecount = 0
        with open(network_file) as fin:
            for line in fin:
                linearray = line.rstrip().split('\t')  # GeneA \t GeneB format
                if linearray[0] in genes_array and linearray[1] in genes_array:
                    for i in [0, 1]:
                        if linearray[i] not in network:
                            network[linearray[i]] = {}
                        network[linearray[i]][linearray[-1 * (i - 1)]] = 1  # save edge information
                    edgecount += 1

        print("Number of network edges: " + str(edgecount))

    #
    # INITIALIZE BFS
    #

    # Define foldchange dynamic threshold. logarithm decay.
    # Parameters are defined by regression (achilles data)  2**-7 was used in previous version.

    FC_THRESH = 2 ** (-1.1535 * np.log(len(np.intersect1d(genes_array,
                                                          nonEss)) + 13.324) + 0.7728)
    bf = {}
    boostedbf = {}
    for g in genes_array:
        for rnatag in gene2rna[g]:
            bf[rnatag] = []

        boostedbf[g] = []  # boosted bf at gene level

    #
    # TRAINING
    #
    if use_small_sample:
        # declare training class
        # training_data = Training(setdiff1d(gene_idx,np.where(in1d(genes_array,coreEss))),cvnum=NUMCV)
        # declare training class (only for Gold-standard gene set)
        training_data = Training(np.where(np.in1d(genes_array, np.union1d(coreEss, nonEss)))[0],
                                 cvnum=no_of_cross_validations)
        # all non-goldstandards
        all_non_gs = np.where(np.logical_not(np.in1d(genes_array, np.union1d(coreEss, nonEss))))[0]
    else:
        training_data = Training(gene_idx, cvnum=no_of_cross_validations)  # declare training class

    if train_method == 0:
        LOOPCOUNT = bootstrap_iterations
    elif train_method == 1:
        LOOPCOUNT = no_of_cross_validations  # 10-folds

    if run_test_mode == True:
        fp = open(output_file + ".traininfo", "w")
        fp.write("#1: Loopcount\n#2: Training set\n#3: Testset\n")
    # No resampling option
    if no_resampling == True:
        print("# Caution: Resampling is disabled")
        LOOPCOUNT = 1

    print("Iter TrainEss TrainNon TestSet")
    sys.stdout.flush()
    for loop in range(LOOPCOUNT):
        currentbf = {}
        printstr = ""
        printstr += str(loop)

        #
        # bootstrap resample (10-folds cross-validation) from gene list to get the training set
        # test set for this iteration is everything not selected in bootstrap resampled (10-folds cross-validation)
        # training set
        # define essential and nonessential training sets:  arrays of indexes
        #
        if no_resampling == True:
            # no resampling
            gene_train_idx = gene_idx
            gene_test_idx = gene_idx
        else:
            # CV or bootstrapping
            gene_train_idx, gene_test_idx = training_data.get_data(train_method)
        if use_small_sample:
            # test set is union of rest of training set (gold-standard) and the other genes (all of non-gold-standard)
            gene_test_idx = np.union1d(gene_test_idx, all_non_gs)

        if run_test_mode:
            fp.write(
                "%d\n%s\n%s\n" % (loop, ",".join(genes_array[gene_train_idx]), ",".join(genes_array[gene_test_idx])))

        train_ess = np.where(np.in1d(genes_array[gene_train_idx], coreEss))[0]
        train_non = np.where(np.in1d(genes_array[gene_train_idx], nonEss))[0]
        printstr += " " + str(len(train_ess))
        printstr += " " + str(len(train_non))
        printstr += " " + str(len(gene_test_idx))
        print(printstr)
        sys.stdout.flush()
        #
        # define ess_train: vector of observed fold changes of essential genes in training set
        #
        ess_train_fc_list_of_lists = [fc[rnatag] for g in genes_array[gene_train_idx[train_ess]] for rnatag in
                                      gene2rna[g]]
        ess_train_fc_flat_list = [obs for sublist in ess_train_fc_list_of_lists for obs in list(sublist.values())]
        #
        # define non_train vector of observed fold changes of nonessential genes in training set
        #
        non_train_fc_list_of_lists = [fc[rnatag] for g in genes_array[gene_train_idx[train_non]] for rnatag in
                                      gene2rna[g]]
        non_train_fc_flat_list = [obs for sublist in non_train_fc_list_of_lists for obs in list(sublist.values())]
        #
        # calculate empirical fold change distributions for both
        #
        kess = stats.gaussian_kde(ess_train_fc_flat_list)
        knon = stats.gaussian_kde(non_train_fc_flat_list)
        #
        # define empirical upper and lower bounds within which to calculate BF = f(fold change)
        #
        x = np.arange(-10, 2, 0.01)
        nonfitx = knon.evaluate(x)
        # define lower bound empirical fold change threshold:  minimum FC np.where knon is above threshold
        f = np.where(nonfitx > FC_THRESH)
        xmin = round_to_hundredth(min(x[f]))
        # define upper bound empirical fold change threshold:  minimum value of log2(ess/non)
        subx = np.arange(xmin, max(x[f]), 0.01)
        logratio_sample = np.log2(kess.evaluate(subx) / knon.evaluate(subx))
        f = np.where(logratio_sample == logratio_sample.min())
        xmax = round_to_hundredth(subx[f])
        #
        # round foldchanges to nearest 0.01
        # precalculate logratios and build lookup table (for speed)
        #
        logratio_lookup = {}
        for i in np.arange(xmin, xmax + 0.01, 0.01):
            logratio_lookup[np.around(i * 100)] = np.log2(kess.evaluate(i) / knon.evaluate(i))
        #
        # calculate BFs from lookup table for withheld test set
        #

        # liner interpolation
        testx = list()
        testy = list()

        for g in genes_array[gene_train_idx]:
            for rnatag in gene2rna[g]:
                for foldchange in list(fc[rnatag].values()):
                    if foldchange >= xmin and foldchange <= xmax:
                        testx.append(np.around(foldchange * 100) / 100)
                        testy.append(logratio_lookup[np.around(foldchange * 100)][0])
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(np.array(testx), np.array(testy))
        except:
            print("Regression failed. Check quality of the screen")
            sys.exit(1)
        #
        # BF calculation
        #

        for g in genes_array[gene_test_idx]:
            for rnatag in gene2rna[g]:
                bayes_factor = []
                for rep in column_list:
                    bayes_factor.append(slope * fc[rnatag][rep] + intercept)
                bf[rnatag].append(bayes_factor)

    if run_test_mode == True:
        fp.close()

    num_obs = dict()
    if rna_level is False:
        bf_mean = dict()
        bf_std = dict()
        bf_norm = dict()  # sgRNA number complement
    if rna_level or filter_multi_target:
        bf_mean_rna_rep = dict()
        bf_std_rna_rep = dict()
    # bf_norm_rna_rep = dict()

    for g in gene2rna:
        num_obs[g] = len(bf[gene2rna[g][0]])
        if rna_level or filter_multi_target:
            for rnatag in gene2rna[g]:
                bf_mean_rna_rep[rnatag] = dict()
                bf_std_rna_rep[rnatag] = dict()
                t = list(zip(*bf[rnatag]))
                for rep in range(len(column_list)):
                    bf_mean_rna_rep[rnatag][column_list[rep]] = np.mean(t[rep])
                    bf_std_rna_rep[rnatag][column_list[rep]] = np.std(t[rep])

        if rna_level == False:
            sumofbf_list = list()
            for i in range(num_obs[g]):
                sumofbf = 0.0
                for rnatag in gene2rna[g]:
                    sumofbf += sum(bf[rnatag][i])
                sumofbf_list.append(sumofbf)  # append each iter
            bf_mean[g] = np.mean(sumofbf_list)
            bf_std[g] = np.std(sumofbf_list)

    #
    # BUILD MULTIPLE REGRESSION MODEL FOR MULTI TARGETING GUIDE RNAs
    #
    if filter_multi_target:
        count = 0
        trainset = dict()
        bf_multi_corrected_gene = dict()
        bf_multi_corrected_rna = dict()
        for gene in gene2rna:
            # multi_targeting_sgrnas_info[seqid] = (perfectmatch, mismatch_1bp, perfectmatch_gene, mismatch_1bp_gene)
            multitarget = list()
            onlytarget = list()
            for seqid in gene2rna[gene]:
                if seqid not in aligninfo.index:
                    continue
                if seqid in multi_targeting_sgrnas_info:
                    multitarget.append(seqid)
                else:
                    onlytarget.append(seqid)

            if len(onlytarget) > 0:  # comparsion between sgRNAs targeting one locus and multiple loci
                if len(multitarget) > 0:

                    bf_only = np.mean([sum(list(bf_mean_rna_rep[seqid].values())) for seqid in onlytarget])
                    for seqid in onlytarget:
                        trainset[seqid] = [1, 0, 0]

                    for seqid in multitarget:
                        if multi_targeting_sgrnas_info[seqid][2] > 1 or multi_targeting_sgrnas_info[seqid][
                            3] > 0:  # train model using multi-targeting only targeting one protein coding gene
                            continue

                        count += 1
                        increment = sum(list(bf_mean_rna_rep[seqid].values())) - bf_only

                        trainset[seqid] = [multi_targeting_sgrnas_info[seqid][0], multi_targeting_sgrnas_info[seqid][1],
                                           increment]

        if count < 10:
            print("Not enough train set for calculating multi-targeting effect.\n")
            print("It may cause due to unmatched gRNA names between the foldchange file and the align info file.\n")
            print("Filtering is not finished\n")
            filter_multi_target = False

        else:

            trainset = pd.DataFrame().from_dict(trainset).T
            X = trainset[[0, 1]]
            y = trainset[2]

            regressor = LinearRegression()
            regressor.fit(X, y)
            coeff_df = pd.DataFrame(regressor.coef_, X.columns, columns=['Coefficient'])
            for i in [0, 1]:
                if coeff_df['Coefficient'][i] < 0:
                    print("Regression coefficient is below than zero. Substituted to zero\n")
                    coeff_df['Coefficient'][i] = 0.0
            print("Multiple effects from perfect matched loci = %.3f and 1bp mis-matched loci = %.3f" % (
                coeff_df['Coefficient'][0], coeff_df['Coefficient'][1]))

            if rna_level == False:
                for g in gene2rna:
                    penalty = 0.0
                    for seqid in gene2rna[g]:
                        if seqid in multi_targeting_sgrnas_info:
                            penalty += float(multi_targeting_sgrnas_info[seqid][0] - 1) * coeff_df['Coefficient'][
                                0] + float(multi_targeting_sgrnas_info[seqid][1]) * coeff_df['Coefficient'][1]
                    bf_multi_corrected_gene[g] = bf_mean[g] - penalty
            else:
                for g in gene2rna:
                    for seqid in gene2rna[g]:
                        if seqid in multi_targeting_sgrnas_info:
                            penalty = float(multi_targeting_sgrnas_info[seqid][0] - 1) * coeff_df['Coefficient'][
                                0] + float(multi_targeting_sgrnas_info[seqid][1]) * coeff_df['Coefficient'][1]
                        else:
                            penalty = 0.0
                        bf_multi_corrected_rna[seqid] = sum(list(bf_mean_rna_rep[seqid].values())) - penalty

    #
    #  NORMALIZE sgRNA COUNT
    #
    if rna_level is False and flat_sgrna == True:
        if filter_multi_target == True:
            targetbf = bf_multi_corrected_gene
        else:
            targetbf = bf_mean

        for g in gene2rna:
            multiple_factor = equalise_sgrna_no / float(len(gene2rna[g]))
            bf_norm[g] = targetbf[g] * multiple_factor

    '''			
    if bf_std[rnatag] == 0.0:
        bf_norm[rnatag] = float('inf')
    else:
        bf_norm[g] = ( bf[rnatag] - bf_mean[rnatag] ) / bf_std[rnatag]
    '''
    training_data = Training(gene_idx)  # set training class reset

    #
    # calculate network scores
    #

    if network_boost == True and rna_level == False:  # Network boost is only working for gene level
        if run_test_mode == True:  # TEST MODE
            fp = open(output_file + ".netscore", "w")
        print("\nNetwork score calculation start\n")

        networkscores = {}
        for g in genes_array[gene_idx]:
            if g in network:
                templist = list()
                for neighbor in network[g]:
                    if neighbor in bf_mean:
                        templist.append(bf_mean[neighbor])

                templist.sort(reverse=True)

                networkscores[g] = fibo_weighted_sum(templist)
        #
        # start training
        #

        for loop in range(LOOPCOUNT):
            currentnbf = {}
            printstr = ""
            printstr += str(loop)

            #
            # draw train, test sets
            #
            gene_train_idx, gene_test_idx = training_data.get_data(train_method)
            #
            # define essential and nonessential training sets:  arrays of indexes
            #
            train_ess = np.where(np.in1d(genes_array[gene_train_idx], coreEss))[0]
            train_non = np.where(np.in1d(genes_array[gene_train_idx], nonEss))[0]
            printstr += " " + str(len(train_ess))
            printstr += " " + str(len(train_non))
            printstr += " " + str(len(gene_test_idx))

            sys.stdout.flush()
            #
            # calculate Network BF for test set
            #
            ess_ns_list = [networkscores[x] for x in genes_array[gene_train_idx[train_ess]] if x in networkscores]
            non_ns_list = [networkscores[x] for x in genes_array[gene_train_idx[train_non]] if x in networkscores]

            kess = stats.gaussian_kde(ess_ns_list)
            knon = stats.gaussian_kde(non_ns_list)
            #
            # set x boundary for liner regression
            #
            testx = list()
            testy = list()
            xmin = float(np.inf)
            xmax = float(-np.inf)

            for networkscore in np.arange(max(ess_ns_list), min(ess_ns_list), -0.01):
                density_ess = kess.evaluate(networkscore)[0]
                density_non = knon.evaluate(networkscore)[0]
                if density_ess == 0.0 or density_non == 0.0:
                    continue

                if np.log2(density_ess / density_non) > -5 and networkscore < np.array(ess_ns_list).mean():  # reverse
                    xmin = min(xmin, networkscore)

            for networkscore in np.arange(min(non_ns_list), max(non_ns_list), 0.01):
                density_ess = kess.evaluate(networkscore)[0]
                density_non = knon.evaluate(networkscore)[0]
                if density_ess == 0.0 or density_non == 0.0:
                    continue
                if np.log2(density_ess / density_non) < 5 and networkscore > np.array(non_ns_list).mean():  # reverse
                    xmax = max(xmax, networkscore)
            #
            # liner regression
            #
            testx = list()
            testy = list()
            for g in genes_array[gene_train_idx]:
                if g in networkscores:
                    if networkscores[g] >= xmin and networkscores[g] <= xmax:
                        testx.append(np.around(networkscores[g] * 100) / 100)
                        testy.append(np.log2(kess.evaluate(networkscores[g])[0] / knon.evaluate(networkscores[g])[0]))

            slope, intercept, r_value, p_value, std_err = stats.linregress(np.array(testx), np.array(testy))

            for g in genes_array[gene_test_idx]:
                if g in networkscores:
                    if run_test_mode == True:
                        fp.write("%s\t%f\t%f\n" % (g, networkscores[g], slope * networkscores[g] + intercept))
                    nbf = slope * networkscores[g] + intercept
                else:
                    nbf = 0.0

                boostedbf[g].append(bf_mean[g] + nbf)
                if flat_sgrna == True:
                    boostedbf[g].append(bf_norm[g] + nbf)

        if run_test_mode == True:
            fp.close()


    

    #
    # print out results
    #

    # Equalizing factor (Replicates)
    if flat_rep==True:
        eqf = equalise_rep_no/float(len(column_labels))
    else:
        eqf = 1

    # print out
    with open(output_file, 'w') as fout:

        if rna_level == True:
            fout.write('RNA\tGENE')
            for i in range(len(column_list)):
                fout.write('\t{0:s}'.format(column_labels[i]))
                if train_method == 0:
                    fout.write('\t{0:s}'.format(column_labels[i] + "_STD"))
            fout.write('\tBF')
            if train_method == 0:
                fout.write('\tNumObs')
            fout.write('\n')

            for rnatag in sorted(bf.keys()):
                # RNA tag
                fout.write('{0:s}\t'.format(rnatag))
                # Gene
                gene = rna2gene[rnatag]
                fout.write('{0:s}\t'.format(gene))

                # BF of replicates
                for rep in column_list:
                    fout.write('{0:4.3f}\t'.format(bf_mean_rna_rep[rnatag][rep]))
                    if train_method == 0:
                        fout.write('{0:4.3f}\t'.format(bf_std_rna_rep[rnatag][rep]))

                # Sum BF of replicates
                if filter_multi_target == True:
                    fout.write('{0:4.3f}'.format(float(bf_multi_corrected_rna[rnatag]) * eqf))  # eqf = equalizing factor for the number of replicates
                else:
                    fout.write('{0:4.3f}'.format(float(sum(list(bf_mean_rna_rep[rnatag].values()))) * eqf))

                # Num obs
                if train_method == 0:
                    fout.write('\t{0:d}'.format(num_obs[gene]))
                fout.write('\n')
        else:
            fout.write('GENE')
            if network_boost == True:
                fout.write('\tBoostedBF')
                if train_method == 0:
                    fout.write('\tSTD_BoostedBF')
            fout.write('\tBF')
            if train_method == 0:
                fout.write('\tSTD\tNumObs')
            if flat_sgrna == True:
                fout.write('\tNormBF')
            fout.write('\n')

            for g in sorted(genes.keys()):
                # Gene
                fout.write('{0:s}'.format(g))
                if network_boost == True:
                    boostedbf_mean = np.mean(boostedbf[g])
                    boostedbf_std = np.std(boostedbf[g])
                    fout.write('\t{0:4.3f}'.format(float(boostedbf_mean) * eqf))
                    if train_method == 0:
                        fout.write('\t{0:4.3f}'.format(float(boostedbf_std) * eqf))

                # BF
                if filter_multi_target == True:
                    fout.write('\t{0:4.3f}'.format(float(bf_multi_corrected_gene[g]) * eqf))  # eqf = equalizing factor for the number of replicates
                else:
                    fout.write('\t{0:4.3f}'.format(float(bf_mean[g]) * eqf ))
                # STD, Count
                if train_method == 0:
                    fout.write('\t{0:4.3f}\t{1:d}'.format(float(bf_std[g]), num_obs[g]))
                # Normalized BF
                if flat_sgrna == True:
                    fout.write('\t{0:4.3f}'.format(float(bf_norm[g])))

                fout.write('\n')


@click.command(name='pr')
@click.option('-i', '--bayes-factors', required=True,
              type=click.Path(exists=True))
@click.option('-o', '--output-file', required=True)
@click.option('-e', '--essential-genes', required=True,
              type=click.Path(exists=True))
@click.option('-n', '--non-essential-genes', required=True, type=click.Path(exists=True))
@click.option('-k', '--use-column', default=None)
def calculate_precision_recall(bayes_factors, output_file, essential_genes, non_essential_genes, use_column):
    """
    Calculate precision-recall from an input Bayes Factors file:

        \b
        BAGEL.py pr -i [Bayes Factor file] -o [output file] -e [essentials genes] -n [nonessentials genes]

    \b
    Required options:
        -i, --bayes-factors [Bayes factors]                     BAGEL output file.
        -o, --output-file [output file]                         Output filename
        -e, --essential-genes [reference essentials]            File with list of training set of essential genes
        -n, --non-essential-genes [reference nonessentials]     File with list of training set of nonessential genes

    \b
    Other options:
        -k [column name]    Use other column (default \BF\)

    \b
    Example:
         BAGEL.py pr  -i input.bf -o output.PR -e ref_essentials.txt -n ref_nonessentials.txt

    """
    #
    # test for availability of all files
    #
    essentials = pd.read_csv(essential_genes, index_col=0, sep="\t")
    nonessentials = pd.read_csv(non_essential_genes, index_col=0, sep="\t")
    bf = pd.read_csv(bayes_factors, index_col=0, sep="\t")

    if use_column is not None:
        bf_column = use_column
        if bf_column not in bf.dtypes.index:
            print("Error! the column name is not in the file")
            sys.exit(1)
    else:
        bf_column = 'BF'

    bf.sort_values(by=bf_column, ascending=False, inplace=True)

    cumulative_tp = 0.
    cumulative_fp = 0.
    precision = 1.
    recall = 0.
    # note float formats

    ess = essentials.index.values
    non = nonessentials.index.values
    totNumEssentials = len([x for x in bf.index.values if x in ess])

    with open(output_file, 'w') as fout:

        fout.write('Gene\t')
        fout.write(bf_column)
        fout.write('\tRecall\tPrecision\tFDR\n')

        for g in bf.index.values:
            if (g in ess):
                cumulative_tp += 1
            elif (g in non):
                cumulative_fp += 1
            recall = cumulative_tp / totNumEssentials
            if ((cumulative_tp > 0) | (cumulative_fp > 0)):
                precision = cumulative_tp / (cumulative_tp + cumulative_fp)
            fout.write('{0:s}\t{1:4.3f}\t{2:4.3f}\t{3:4.3f}\t{4:4.3f}\n'.format(g, bf.loc[g, bf_column], recall, precision, 1.0-precision))


if __name__ == '__main__':
    bagel.add_command(calculate_fold_change)
    bagel.add_command(calculate_bayes_factors)
    bagel.add_command(calculate_precision_recall)
    bagel.add_command(report_bagel_version)
    bagel()
