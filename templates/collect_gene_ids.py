#!/usr/bin/env python3

data = "${gene_data}"
question = "${gpt_question}"

# Read input file
with open(data, 'r') as input_file:
    lines = input_file.readlines()

# Extract gene names while ignoring header row
gene_names = [line.split('\t')[0] for line in lines[1: ]]

# Write question and gene names to output file
with open("query.txt", 'w') as output_file:
    output_file.write(question + '''\n''')

    for gene in gene_names:
        output_file.write(gene + '''\n''')