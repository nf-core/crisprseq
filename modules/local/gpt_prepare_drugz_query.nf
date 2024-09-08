process GPT_PREPARE_DRUGZ_QUERY {
    input:
    path data
    val amount
    val question

    output:
    path 'gpt_drugz_query.txt', emit: query

    script:
    template 'findTopDrugzGenes.py'
}
