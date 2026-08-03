include { CRISPRDECODE_VALIDATE_LIBRARY } from '../../../modules/local/crisprdecode/validate/main'
include { CRISPRDECODE_ASSIGN           } from '../../../modules/local/crisprdecode/assign/main'
include { CRISPRDECODE_AGGREGATE        } from '../../../modules/local/crisprdecode/aggregate/main'

workflow CRISPRDECODE_PAIRED_GUIDE {
    take:
    reads   // channel: [ val(meta), path([read_1, read_2]) ]
    library // channel: path(construct_library.tsv)

    main:
    CRISPRDECODE_VALIDATE_LIBRARY(library)

    reads
        .map { meta, read_files ->
            if (meta.single_end || read_files.size() != 2) {
                error "CRISPRDecode paired-guide counting requires exactly two FASTQ files per sample; sample '${meta.id}' is not paired-end."
            }
            [meta, read_files]
        }
        .combine(CRISPRDECODE_VALIDATE_LIBRARY.out.library)
        .set { ch_validate_for_assign }

    CRISPRDECODE_ASSIGN(
        ch_validate_for_assign,
        params.crisprdecode_r1_anchor ?: '',
        params.crisprdecode_r2_anchor ?: '',
        params.crisprdecode_r1_offset,
        params.crisprdecode_r2_offset,
        params.crisprdecode_reverse_complement_r1,
        params.crisprdecode_reverse_complement_r2
    )

    CRISPRDECODE_ASSIGN.out.counts
        .map { meta, counts -> counts }
        .collect()
        .set { ch_assign_counts_for_aggregate }

    CRISPRDECODE_ASSIGN.out.summary
        .map { meta, summary -> summary }
        .collect()
        .set { ch_assign_summaries_for_aggregate }

    CRISPRDECODE_AGGREGATE(
        ch_assign_counts_for_aggregate,
        ch_assign_summaries_for_aggregate,
        CRISPRDECODE_VALIDATE_LIBRARY.out.library
    )

    ch_versions = CRISPRDECODE_VALIDATE_LIBRARY.out.versions
        .mix(CRISPRDECODE_ASSIGN.out.versions)
        .mix(CRISPRDECODE_AGGREGATE.out.versions)

    emit:
    counts             = CRISPRDECODE_AGGREGATE.out.count_matrix
    assignment_summary = CRISPRDECODE_AGGREGATE.out.assignment_summary
    library_recovery   = CRISPRDECODE_AGGREGATE.out.library_recovery
    versions           = ch_versions
}
