#!/usr/bin/env python

import sys

import pysam

mapped_reads_count = int(pysam.view("-c", "-b", "-F", "4", "$reads"))
total_reads_count = int(pysam.view("-c", "-b", "$reads"))
mapped_reads_percentage = mapped_reads_count * 100 / total_reads_count

with open("$summary", "r") as summary:
    summary_lines = summary.readlines()

add_line = True
with open("$summary", "w") as output_file:
    for line in summary_lines:
        if "aligned-reads" not in line:
            output_file.write(line)
        else:
            output_file.write(f"aligned-reads, {mapped_reads_count} ({round(mapped_reads_percentage, 1)}%)\\n")
            add_line = False
    if add_line:
        output_file.write(f"aligned-reads, {mapped_reads_count} ({round(mapped_reads_percentage, 1)}%)\\n")
