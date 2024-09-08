process GPT_PREPARE_BAGEL2_QUERY {
    input:
    path data
    val amount
    val question

    output:
    path 'gpt_bagel2_query.txt', emit: query

    script:
    template 'findTopBagel2Genes.py'
}
