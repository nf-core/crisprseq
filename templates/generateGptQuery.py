#!/usr/bin/env python

# Open data file
with open("${data}", "r") as file:
    # Read the header row
    header = file.readline().strip().split("\t")

    # Ensure target_index is within the bounds of available columns
    if int("${index}") >= len(header) or int("${index}") < 0:
        raise ValueError(
            f"Error: The specified column index {int("${index}")} is out of range!"
        )

    # Initiate list to store gene IDs with corresponding data values
    data = []

    # Prepare data file's rows
    for line in file:
        row = line.strip().split("\t")
        # Assume the first column is the gene ID
        gene_id = row[0]
        # Extract value using the provided index
        value = float(row[int("${index}")])
        data.append((gene_id, value))

    # Sort the data based on provided mode
    if "${mode}" == "low":
        sorted_data = sorted(data, key=lambda x: x[1])
    elif "${mode}" == "high":
        sorted_data = sorted(data, key=lambda x: x[1], reverse=True)
    else:
        raise ValueError("Error: Please provide either 'low' or 'high' as mode.")

    # Extract num_genes many top genes
    top_gene_ids = [gene_id for gene_id, value in sorted_data[: int("${count}")]]

    # Write everything into an output file
    with open(f"gpt_{"${source}"}_query.txt", "w") as query_file:
        query_file.write("${question}" + """\n""")
        for gene_id in top_gene_ids:
            query_file.write(gene_id + """\n""")
