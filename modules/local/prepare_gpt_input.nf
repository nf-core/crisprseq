process PREPARE_GPT_INPUT {
    input:
    path gene_data
    val gpt_question

    output:
    path 'query.txt', emit: query

    script:
    template 'collect_gene_ids.py'
}
