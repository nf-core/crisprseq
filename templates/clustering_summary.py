#!/usr/bin/env python

import gzip

from Bio import SeqIO

with gzip.open("$reads", "rt") as handle:
    clusters_count = len(list(SeqIO.parse(handle, "fastq")))

with open("$summary", "a") as output_file:
    output_file.write(f"clustered-reads, {clusters_count}\\n")
