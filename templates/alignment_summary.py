
#!/usr/bin/env python

import pysam

mapped_reads_count = int(pysam.view("-c", "-b", "-F", "4", "$reads"))
total_reads_count = int(pysam.view("-c", "-b", "$reads"))
mapped_reads_percentage = mapped_reads_count * 100 / total_reads_count

with open("$summary", "a") as output_file:
    output_file.write(f"aligned-reads, {mapped_reads_count} ({round(mapped_reads_percentage, 1)}%)\\n")
