//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .multiMap {
            reads:       create_fastq_channel(it)
            reference:   create_reference_channel(it)
            template:    create_template_channel(it)
            protospacer: create_protospacer_channel(it)
        }
        .set { inputs }

    emit:
    reads = inputs.reads                      // channel: [ val(meta), [ reads ] ]
    reference = inputs.reference              // channel: [ val(meta), reference ]
    template = inputs.template                // channel: [ val(meta), template ]
    protospacer = inputs.protospacer          // channel: [ val(meta), protospacer ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.sample

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
        meta.single_end = true
    } else {
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
        meta.single_end = false
    }
    return fastq_meta
}

// Function to get a list of [ meta, reference ]
def create_reference_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.sample

    // add reference sequence to meta
    def reference_meta = []
    if (row.reference.length() <= 0) {
        meta.self_reference = true
        reference_meta = [ meta, null ]
    } else {
        meta.self_reference = false
        reference_meta = [ meta, row.reference ]
    }

    return reference_meta
}

// Function to  get a list of [ meta, template ]
def create_template_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.sample

    // add template sequence/path to meta
    def template_meta = []
    if (!row.template) {
        meta.template = false
        template_meta = [ meta, null ]
    } else if (!file(row.template).exists()) {
        meta.template = true
        template_meta = [ meta, row.template ]
    } else {
        meta.template = true
        template_meta = [ meta, file(row.template) ]
    }

    return template_meta
}

// Function to get a list of [ meta, protospacer ]
def create_protospacer_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.sample

    // add protospacer sequence to meta
    def protospacer_meta = []
    if (row.protospacer.length() <= 0) {
        exit 1, "ERROR: Please check input samplesheet -> Protospacer sequence is not provided!\n"
    } else {
        protospacer_meta = [ meta, row.protospacer ]
    }

    return protospacer_meta
}
