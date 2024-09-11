#!/usr/bin/env python

# Define process input variables
data_path = "${data}"
source = "${source}"
target_column = "${column}"
num_genes = "${count}"
num_genes = int(num_genes)
mode = "${mode}"
question = "${question}"

# Open data file
with open(data_path, "r") as file:
    # Read the header row and split into column names
    header = file.readline().strip().split("\t")

    # Check if specified column exists
    if target_column not in header:
        print(
            f"Error: The specified column '{target_column}' was not found in the data file!"
        )

    # Find the target column index
    target_index = header.index(target_column)

    # Initiate list to store gene id's with corresponding data values
    data = []

    # Prepare data file's rows
    for line in file:
        row = line.strip().split("\t")
        gene_id = row[0]
        value = float(row[target_index])
        data.append((gene_id, value))

    # Sort the data based on provided mode
    if mode == "low":
        sorted_data = sorted(data, key=lambda x: x[1])
    elif mode == "high":
        sorted_data = sorted(data, key=lambda x: x[1], reverse=True)
    else:
        print("Error: Please provide either 'low' or 'high' as mode.")

    # Extract num_genes many top genes
    top_gene_ids = [gene_id for gene_id, value in sorted_data[:num_genes]]

    # Write everything into a output file
    with open("gpt_${source}_query.txt", "w") as query_file:
        query_file.write(question + """\n""")
        for gene_id in top_gene_ids:
            query_file.write(gene_id + """\n""")
