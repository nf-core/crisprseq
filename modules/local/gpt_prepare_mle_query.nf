process GPT_PREPARE_MLE_QUERY {
    input:
    path data
    val amount
    val question

    output:
    path 'gpt_mle_query.txt', emit: query

    script:
    template 'findTopMleGenes.py'
}