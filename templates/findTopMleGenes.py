#!/usr/bin/env python

# Define variables
file_path = "${data}"
num_of_genes = "${amount}"
num_of_genes = int(num_of_genes)
question = "${question}"


# Function to extract relevant data from the line (split by tabs)
def extract_gene_info(line):
    fields = line.strip().split("\t")
    gene = fields[0]
    p_value = float(
        fields[4]
    )  # control_vs_treatment|p-value is in the 5th column (index 4)
    return gene, p_value


# Read the data and store the genes with their p-values
genes_with_p_values = []
with open(file_path, "r") as file:
    # Skip the header line
    header = file.readline()

    # Process the rest of the lines
    for line in file:
        gene, p_value = extract_gene_info(line)
        genes_with_p_values.append((gene, p_value))

# Sort genes based on p-value (smallest first)
genes_with_p_values.sort(key=lambda item: item[1])

# Get the top genes with the smallest p-values
top_genes = genes_with_p_values[:num_of_genes]

# Write the question and top genes to the output file
with open("gpt_mle_query.txt", "w") as f:
    # Write the question
    f.write(question + """\n""")

    # Write the top genes (only their names)
    for gene, _ in top_genes:
        f.write(gene + """\n""")
