process GPT_PREPARE_QUERY {
    input:
    path data
    val source
    val column
    val count
    val mode
    val question

    output:
    path "gpt_${source}_query.txt", emit: query

    script:
    template 'generateGptQuery.py'
}
