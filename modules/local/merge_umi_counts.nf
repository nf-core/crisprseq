process MERGE_UMI_COUNTS {
    tag "merge"
    label 'process_single'

    // Use umi_tools container which has Python - we'll write a simple merge without pandas
    conda 'bioconda::umi_tools=1.1.5'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.5--py39hf95cd2a_0' :
        'biocontainers/umi_tools:1.1.5--py39hf95cd2a_0' }"

    input:
    path(count_files)
    path(library)

    output:
    path("merged_counts.txt"), emit: counts
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import sys
    import os
    from collections import defaultdict

    count_files = "${count_files}".split()
    
    # Read library to get all sgRNAs and genes
    sgrna_to_gene = {}
    all_sgrnas = []
    with open("${library}", 'r') as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\\t')
            if len(parts) >= 3:
                sgrna, seq, gene = parts[0], parts[1], parts[2]
                sgrna_to_gene[sgrna] = gene
                all_sgrnas.append(sgrna)
    
    # Initialize counts dict
    counts = {sgrna: {} for sgrna in all_sgrnas}
    sample_names = []
    
    # Read each count file
    for count_file in count_files:
        if os.path.exists(count_file):
            with open(count_file, 'r') as f:
                header = f.readline().strip().split('\\t')
                sample_name = header[2] if len(header) >= 3 else os.path.basename(count_file).replace('.count.tsv', '')
                sample_names.append(sample_name)
                
                for line in f:
                    parts = line.strip().split('\\t')
                    if len(parts) >= 3:
                        sgrna = parts[0]
                        count = int(parts[2]) if parts[2].isdigit() else 0
                        if sgrna in counts:
                            counts[sgrna][sample_name] = count
    
    # Write merged output
    with open("merged_counts.txt", 'w') as f:
        # Header
        f.write('sgRNA\\tGene\\t' + '\\t'.join(sample_names) + '\\n')
        
        # Data rows
        for sgrna in all_sgrnas:
            gene = sgrna_to_gene.get(sgrna, 'Unknown')
            row_counts = [str(counts[sgrna].get(s, 0)) for s in sample_names]
            f.write(f'{sgrna}\\t{gene}\\t' + '\\t'.join(row_counts) + '\\n')
    
    # Version info
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
    """
}
