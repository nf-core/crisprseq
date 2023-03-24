#!/usr/bin/env python3

#### Extract the read sequence from a FASTA file of UMIs (output of extract_umis.py)
#### author: Júlia Mir Pedrol @mirpedrol
#### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.

with open("$cluster") as i:
    for line in i:
        split = line.split(";")
        with open("${cluster.baseName}_${task.ext.prefix}.fasta", "w") as o:
            o.write(f"{line.split(';')[0]}\n{line.split('=')[-1]}")
