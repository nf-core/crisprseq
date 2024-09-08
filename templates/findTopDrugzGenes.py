#!/usr/bin/env python

# Define variables
file_path = "${data}"
num_of_genes = "${amount}"
num_of_genes = int(num_of_genes)
question = "${question}"

# Read the file and split it into lines
with open(file_path, "r") as file:
    lines = file.readlines()

# Extract the header and the rows
header = lines[0].strip().split("\t")
rows = [line.strip().split("\t") for line in lines[1:]]

# Sort rows based on "rank_synth" (convert the value to int for proper sorting)
rows_sorted = sorted(rows, key=lambda row: int(row[header.index("rank_synth")]))

# Extract the top num_of_genes genes
top_genes = rows_sorted[:num_of_genes]

# Write the output to query.txt
with open("gpt_drugz_query.txt", "w") as output:
    output.write(question + """\n""")
    for gene in top_genes:
        output.write(gene[header.index("GENE")] + """\n""")
