// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process BAGEL2_GRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "python=3.11.4 pandas=2.0.3 numpy=1.25.1 scikit-learn=1.3.0 click=8.1.6 matplotlib=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-54e0353146eca1531516863e8235bf7385d76663:c9ff1a9eec871c54cbea815eae778da702623978-0':
        'biocontainers/mulled-v2-54e0353146eca1531516863e8235bf7385d76663:c9ff1a9eec871c54cbea815eae778da702623978-0' }"


    input:
    tuple val(meta), path(pr)

    output:
    //path "versions.yml"           , emit: versions
    path("*.png")                   , emit: pictures

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3

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

    """


}
