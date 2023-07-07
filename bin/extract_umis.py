#!/usr/bin/env python3
#
# This code is obtained from Oxford Nanopore Technologies pipeline-umi-amplicon
# Distributed under the Mozilla Public License, v. 2.0
#
# Original source code: https://github.com/nanoporetech/pipeline-umi-amplicon/blob/69ec3907879aea406b5fb02e3db83b579bfb8b45/lib/umi_amplicon_tools/extract_umis.py
#

import argparse
import logging
import sys

import edlib
import pysam
from tqdm import tqdm


def rev_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def str2bool(v):
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def parse_args(argv):
    """
    Commandline parser
    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to telemap"
    parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
            "CRITICAL",
            "debug",
            "info",
            "warning",
            "error",
            "critical",
        ],
        default="INFO",
        help="Print debug information",
    )
    parser.add_argument(
        "--max-error",
        dest="MAX_ERROR",
        type=int,
        default=2,
        help="Max edit distance for UMI",
    )
    parser.add_argument(
        "--adapter-length",
        dest="ADAPTER_LENGTH",
        type=int,
        default=250,
        help="Length of adapter",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="THREADS",
        type=int,
        default=1,
        help="Number of threads.",
    )
    parser.add_argument("--tsv", dest="TSV", type=str, required=False, help="TSV output file")
    parser.add_argument("-o", "--output", dest="OUT", type=str, required=False, help="FASTA output file")
    parser.add_argument(
        "--fwd-context",
        dest="FWD_CONTEXT",
        type=str,
        default="GTATCGTGTAGAGACTGCGTAGG",
        help="Forward upstream sequence",
    )
    parser.add_argument(
        "--rev-context",
        dest="REV_CONTEXT",
        type=str,
        default="AGTGATCGAGTCAGTGCGAGTG",
        help="Reverse upstream sequence",
    )
    parser.add_argument(
        "--fwd-umi",
        dest="FWD_UMI",
        type=str,
        default="TTTVVVVTTVVVVTTVVVVTTVVVVTTT",
        help="Forward UMI sequence",
    )
    parser.add_argument(
        "--rev-umi",
        dest="REV_UMI",
        type=str,
        default="AAABBBBAABBBBAABBBBAABBBBAAA",
        help="Reverse UMI sequence",
    )
    parser.add_argument("INPUT_FA", type=str, nargs="+", default="/dev/stdin", help="Detected UMIs")

    args = parser.parse_args(argv)

    return args


def align(query, pattern_info, max_ed, normalise=False):
    pattern, forward = pattern_info

    # move this somewhere not in a loop
    seq = pattern
    for c in "actgACTG":
        seq = seq.replace(c, "")
    wildcard = set("".join(seq))

    equalities = [
        ("M", "A"),
        ("M", "C"),
        ("R", "A"),
        ("R", "A"),
        ("W", "A"),
        ("W", "A"),
        ("S", "C"),
        ("S", "C"),
        ("Y", "C"),
        ("Y", "C"),
        ("K", "G"),
        ("K", "G"),
        ("V", "A"),
        ("V", "C"),
        ("V", "G"),
        ("H", "A"),
        ("H", "C"),
        ("H", "T"),
        ("D", "A"),
        ("D", "G"),
        ("D", "T"),
        ("B", "C"),
        ("B", "G"),
        ("B", "T"),
        ("N", "G"),
        ("N", "A"),
        ("N", "T"),
        ("N", "C"),
    ]

    result = edlib.align(
        pattern,
        query,
        task="path",
        mode="HW",
        k=max_ed,
        additionalEqualities=equalities,
    )
    if result["editDistance"] == -1:
        return None, None

    ed = result["editDistance"]
    if not normalise:
        locs = result["locations"][0]
        umi = query[locs[0] : locs[1] + 1]
        return ed, umi

    # Extract and normalise UMI
    umi = ""
    align = edlib.getNiceAlignment(result, pattern, query)
    for q, t in zip(align["query_aligned"], align["target_aligned"]):
        if q not in wildcard:
            continue
        if t == "-":
            umi += "N"
        else:
            umi += t

    if len(umi) != 16:
        raise RuntimeError("UMI length incorrect: {}".format(umi))

    return ed, umi


def count_reads_fastx(fasta_filename):
    n_read = 0
    logging.info("Counting reads in {}".format(fasta_filename))
    with pysam.FastxFile(fasta_filename) as fh:
        for entry in fh:
            n_read += 1
    return n_read


def extract_adapters(entry, max_adapter_length):
    read_5p_seq = None
    read_3p_seq = None
    if len(entry.sequence) > max_adapter_length:
        read_5p_seq = entry.sequence[:max_adapter_length]
        read_3p_seq = entry.sequence[-max_adapter_length:]

    return read_5p_seq, read_3p_seq


def detect_read_strand(read_5p_seq, upstream_context_fwd, upstream_context_rev):
    result_fwd = edlib.align(upstream_context_fwd, read_5p_seq, mode="HW", k=5)
    result_rev = edlib.align(upstream_context_rev, read_5p_seq, mode="HW", k=5)
    strand = "-"
    if result_fwd["editDistance"] < result_rev["editDistance"]:
        strand = "+"

    return strand


def combine_umis(seq_5p, seq_3p, strand):
    seq = seq_5p + seq_3p
    if strand != "+":
        seq = rev_comp(seq_3p) + rev_comp(seq_5p)
    return seq


def print_seq(
    entry,
    strand,
    result_5p_fwd_dist,
    result_3p_rev_dist,
    result_5p_fwd_seq,
    result_3p_rev_seq,
    out,
):
    seq = combine_umis(result_5p_fwd_seq, result_3p_rev_seq, strand)
    print(
        ">{};strand={};umi_fwd={};umi_rev={};umi_fwd_seq={};umi_rev_seq={};seq={}".format(
            entry.name,
            strand,
            result_5p_fwd_dist,
            result_3p_rev_dist,
            result_5p_fwd_seq,
            result_3p_rev_seq,
            entry.sequence,
        ),
        file=out,
    )
    print(seq, file=out)


def extract_umis(
    input_files,
    max_adapter_length,
    max_pattern_dist,
    upstream_context_fwd,
    upstream_context_rev,
    pattern_fwd,
    pattern_rev,
    output_file,
    out,
    tsv,
):
    n_both_umi = 0
    strand_stats = {"+": 0, "-": 0}
    for fasta_filename in input_files:
        n_read = count_reads_fastx(fasta_filename)
        with tqdm(total=n_read) as pbar:
            with pysam.FastxFile(fasta_filename) as fh:
                for entry in fh:
                    pbar.update(1)

                    read_5p_seq, read_3p_seq = extract_adapters(entry, max_adapter_length)
                    if not read_5p_seq or not read_3p_seq:
                        continue

                    strand = detect_read_strand(read_5p_seq, upstream_context_fwd, upstream_context_rev)
                    strand_stats[strand] += 1

                    # Extract fwd UMI
                    result_5p_fwd_dist, result_5p_fwd_seq = align(read_5p_seq, pattern_fwd, max_pattern_dist)
                    # Extract rev UMI
                    result_3p_rev_dist, result_3p_rev_seq = align(read_3p_seq, pattern_rev, max_pattern_dist)

                    if not result_5p_fwd_seq or not result_3p_rev_seq:
                        continue

                    n_both_umi += 1
                    print_seq(
                        entry,
                        strand,
                        result_5p_fwd_dist,
                        result_3p_rev_dist,
                        result_5p_fwd_seq,
                        result_3p_rev_seq,
                        out,
                    )

        fwd_rev_ratio = -1
        perc = -1.0
        if strand_stats["-"]:
            fwd_rev_ratio = strand_stats["+"] / strand_stats["-"]
        logging.info(
            "Found {} fwd and {} rev reads (ratio: {})".format(strand_stats["+"], strand_stats["-"], fwd_rev_ratio)
        )
        if n_read:
            perc = 100.0 * n_both_umi / n_read
            logging.info("{}% of reads contained both UMIs with max {} mismatches".format(perc, max_pattern_dist))
        if tsv:
            print(
                output_file,
                max_pattern_dist,
                strand_stats["+"],
                strand_stats["-"],
                fwd_rev_ratio,
                perc,
                file=tsv,
            )

    out.close()


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.
    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log.upper())
    logging.basicConfig(level=numeric_level, format="%(message)s")

    max_adapter_length = args.ADAPTER_LENGTH
    max_pattern_dist = args.MAX_ERROR
    upstream_context_fwd = args.FWD_CONTEXT
    upstream_context_rev = args.REV_CONTEXT
    output_file = args.OUT
    tsv_file = args.TSV
    input_files = args.INPUT_FA
    umi_fwd = args.FWD_UMI
    umi_rev = args.REV_UMI

    pattern_fwd = (
        umi_fwd,
        True,
    )
    pattern_rev = (
        umi_rev,
        False,
    )

    tsv = None
    if tsv_file:
        tsv = open(tsv_file, "w")

    with open(output_file, "w") as out:
        extract_umis(
            input_files,
            max_adapter_length,
            max_pattern_dist,
            upstream_context_fwd,
            upstream_context_rev,
            pattern_fwd,
            pattern_rev,
            output_file,
            out,
            tsv,
        )
    if tsv_file:
        tsv.close()


if __name__ == "__main__":
    main()
