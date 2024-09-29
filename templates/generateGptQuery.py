#!/usr/bin/env python

# Define process input variables
data_path = "${data}"
source = "${source}"
target_index = "${index}"
target_index = int(target_index)
num_genes = "${count}"
num_genes = int(num_genes)
mode = "${mode}"
question = "${question}"

# Open data file
with open(data_path, "r") as file:
    # Read the header row (though it’s not necessary for column names anymore)
    header = file.readline().strip().split("\t")
    
    # Ensure target_index is within the bounds of available columns
    if target_index >= len(header) or target_index < 0:
        print(f"Error: The specified column index {target_index} is out of range!")
    
    # Initiate list to store gene IDs with corresponding data values
    data = []

    # Prepare data file's rows
    for line in file:
        row = line.strip().split("\t")
        gene_id = row[0]  # Assume the first column is the gene ID
        value = float(row[target_index])  # Extract value using the provided index
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

    # Write everything into an output file
    with open(f"gpt_{source}_query.txt", "w") as query_file:
        query_file.write(question + """\n""")
        for gene_id in top_gene_ids:
            query_file.write(gene_id + """\n""")
