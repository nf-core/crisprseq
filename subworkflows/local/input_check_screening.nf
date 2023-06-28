//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK_SCREENING } from '../../modules/local/samplesheet_check_screening'

workflow INPUT_CHECK_SCREENING {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK_SCREENING ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                               // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK_SCREENING.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.condition  = row.condition

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        error("ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}")
    }
    if (!row.fastq_2) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            error("ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}")
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}
