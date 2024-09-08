#!/usr/bin/env python

# Define variables
file_path = "${data}"
num_of_genes = "${amount}"
num_of_genes = int(num_of_genes)
question = "${question}"

# Open and read the data file
with open(file_path, "r") as file:
    lines = file.readlines()

# Skip the header row and process the data
genes_data = []

for line in lines[1:]:
    gene_info = line.strip().split("\t")  # Split the line by tab
    gene = gene_info[0]  # Get the gene name
    bf_value = float(gene_info[1])  # Get the BF value as a float
    genes_data.append((gene, bf_value))  # Store gene and its BF value as tuple

# Sort the genes_data by BF values in descending order to get the most positive ones
for i in range(len(genes_data)):
    for j in range(i + 1, len(genes_data)):
        if genes_data[i][1] < genes_data[j][1]:
            # Swap the elements
            temp = genes_data[i]
            genes_data[i] = genes_data[j]
            genes_data[j] = temp

# Get the top genes with the most positive BF values
top_genes = genes_data[:num_of_genes]

# Write the question and top genes into a new output file
with open("gpt_bagel2_query.txt", "w") as outfile:
    outfile.write(question + """\n""")
    for gene, bf_value in top_genes:
        outfile.write(gene + """\n""")
