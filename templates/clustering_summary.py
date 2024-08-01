#!/usr/bin/env python

############################
#### Summary of clustering
#### author: Júlia Mir @mirpedrol
#### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
############################

import gzip

import Bio
from Bio import SeqIO

with gzip.open("$reads", "rt") as handle:
    clusters_count = len(list(SeqIO.parse(handle, "fastq")))

with open("$summary", "r") as summary:
    summary_lines = summary.readlines()

add_line = True
outname = "$summary".replace("_preprocessing_summary.csv", "_clustering_summary.csv")
with open(outname, "w") as output_file:
    for line in summary_lines:
        if "clustered-reads" not in line:
            output_file.write(line)
        else:
            output_file.write(f"clustered-reads, {clusters_count}\\n")
            add_line = False
    if add_line:
        output_file.write(f"clustered-reads, {clusters_count}\\n")

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'  biopython: "{Bio.__version__}"\\n')
