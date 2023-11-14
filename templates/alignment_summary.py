#!/usr/bin/env python

############################
#### Summary of alignment
#### author: Júlia Mir @mirpedrol
#### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
############################

import sys

import pysam

mapped_reads_count = int(pysam.view("-c", "-b", "-F", "4", "$reads"))
total_reads_count = int(pysam.view("-c", "-b", "$reads"))
mapped_reads_percentage = mapped_reads_count * 100 / total_reads_count

with open("$summary", "r") as summary:
    summary_lines = summary.readlines()

add_line = True
outname = "$summary".replace("_clustering_summary.csv", "_alignment_summary.csv")
with open(outname, "w") as output_file:
    for line in summary_lines:
        if "aligned-reads" not in line:
            output_file.write(line)
        else:
            output_file.write(f"aligned-reads, {mapped_reads_count} ({round(mapped_reads_percentage, 1)}%)\\n")
            add_line = False
    if add_line:
        output_file.write(f"aligned-reads, {mapped_reads_count} ({round(mapped_reads_percentage, 1)}%)\\n")

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'  pysam: "{pysam.__version__}"\\n')
