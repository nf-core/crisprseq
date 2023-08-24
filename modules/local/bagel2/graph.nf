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


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    /*
    # ---------------------------------
    # BAGEL:  Bayesian Analysis of Gene EssentaLity
    # (c) Traver Hart <traver@hart-lab.org>, Eiru Kim <rooeikim@gmail.com> 2017.

    # Acknowledgements: John McGonigle <j.e.mcgonigle@gmail.com>
    # modified 10/2019
    # Free to modify and redistribute with attribution
    # ---------------------------------

    # ------------------------------------
    # constants

"""
    MIT License

    Copyright (c) 2020 Hart Lab

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

        */
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
