process UMI_DEDUP_SCREENING {
    tag "$meta.id"
    label 'process_medium'

    // Uses library-guided matching approach - no need for Starcode
    // The script matches each read directly to the library (exact + fuzzy)
    // then clusters UMIs using umi_tools directional algorithm
    conda 'bioconda::umi_tools=1.1.5'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.5--py39hf95cd2a_0' :
        'biocontainers/umi_tools:1.1.5--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(reads)
    path(library)

    output:
    tuple val(meta), path("*.count.tsv"), emit: counts
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // UMI clustering parameters
    def umi_edit_dist = params.umi_edit_distance ?: 1
    def umi_method = params.umi_clustering_method ?: 'directional'
    
    // sgRNA matching parameters
    def sgrna_edit_dist = params.sgrna_edit_distance ?: 1
    def sgrna_exact = params.sgrna_exact_only ?: false
    def sgrna_start = params.sgRNA_start ?: 0
    def sgrna_length = params.sgRNA_length ?: 20
    
    def r1 = reads instanceof List ? reads[0] : reads
    def exact_flag = sgrna_exact ? '--sgrna-exact-only' : ''
    
    """
    umi_dedup_screening.py \\
        --r1 ${r1} \\
        --library ${library} \\
        --output ${prefix}.count.tsv \\
        --umi-edit-dist ${umi_edit_dist} \\
        --umi-method ${umi_method} \\
        --sgrna-edit-dist ${sgrna_edit_dist} \\
        --sgRNA-start ${sgrna_start} \\
        --sgRNA-length ${sgrna_length} \\
        --sample-name ${meta.id} \\
        ${exact_flag} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        umi_tools: \$(umi_tools --version 2>&1 | grep -oP 'UMI-tools version: \\K[0-9.]+' || echo '1.1.5')
        starcode: \$(starcode --version 2>&1 | grep -oP '[0-9.]+' || echo 'not available')
    END_VERSIONS
    """
}
