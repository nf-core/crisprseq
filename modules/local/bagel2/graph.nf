process BAGEL2_GRAPH {
    tag "${meta.treatment}_${meta.reference}"
    label 'process_single'

    conda "python=3.11.4 pandas=2.0.3 numpy=1.25.1 scikit-learn=1.3.0 click=8.1.6 matplotlib=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-54e0353146eca1531516863e8235bf7385d76663:c9ff1a9eec871c54cbea815eae778da702623978-0':
        'biocontainers/mulled-v2-54e0353146eca1531516863e8235bf7385d76663:c9ff1a9eec871c54cbea815eae778da702623978-0' }"

    input:
    tuple val(meta), path(pr)

    output:
    path("*.png")                   , emit: pictures
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3

    #### author: Laurence Kuhlburger
    #### Released under the MIT license. See git repository (https://github.com/nf-core/crisprseq) for full license text.
    ####
    #### Orient a reference sequence according to reads orientation.


    import pandas as pd
    import matplotlib.pyplot as plt

    pr_data = pd.read_table('${pr}', index_col=0)

    plt.plot(pr_data.Recall, pr_data.Precision, linewidth=2)
    plt.xlim(0, 1.01)
    plt.ylim(0, 1.02)
    plt.xlabel('Recall')
    plt.ylabel('Precision (1-FDR)')
    plt.title('Precision-Recall Plot')

    # Save the plot to a PNG file
    file_name = 'PR_plot_{}_vs_{}.png'.format('${meta.treatment}', '${meta.reference}')

    plt.savefig(file_name)

    # Show the plot (optional)
    plt.show()

    pr_data.hist('BF', bins=50, range=(-100,100))
    plt.xlabel('Bayes Factor')
    plt.ylabel('Number of Genes')

    file_name = 'barplot_{}_vs_{}.png'.format('${meta.treatment}', '${meta.reference}')

    plt.savefig(file_name)
    plt.show()

    # Output version information
    version = pd. __version__
    matplotlib_version = plt.matplotlib.__version__
    # alas, no `pyyaml` pre-installed in the cellranger container
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'  pandas: "{version}"\\n')
        f.write(f'  matplotlib.pyplot: "{matplotlib_version}"\\n')


    """


}
