#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/crisprseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/crisprseq
    Website: https://nf-co.re/crisprseq
    Slack  : https://nfcore.slack.com/channels/crisprseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_crisprseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_crisprseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_crisprseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.reference_fasta = params.reference_fasta ?: getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CRISPRSEQ_TARGETED  } from './workflows/crisprseq_targeted'
include { CRISPRSEQ_SCREENING } from './workflows/crisprseq_screening'

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CRISPRSEQ {

    take:
    reads_targeted // channel: fastqc files read in from --input
    reads_screening // channel: fastqc files read in from --input
    reference // channel: reference sequence read from --input
    protospacer // channel: protospacer sequence read from --input
    template // channel: template sequence read from --input

    main:
    //
    // WORKFLOW: Run pipeline
    //
    if ( params.analysis == "targeted" ) {
        CRISPRSEQ_TARGETED (
            reads_targeted,
            reference,
            template,
            protospacer
        )
        multiqc_report_ch = CRISPRSEQ_TARGETED.out.multiqc_report
    } else if ( params.analysis == "screening" ) {
        CRISPRSEQ_SCREENING (reads_screening)
        multiqc_report_ch = CRISPRSEQ_SCREENING.out.multiqc_report
    }

    emit:
    multiqc_report = multiqc_report_ch // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CRISPRSEQ (
        PIPELINE_INITIALISATION.out.reads_targeted,
        PIPELINE_INITIALISATION.out.fastqc_screening,
        PIPELINE_INITIALISATION.out.reference,
        PIPELINE_INITIALISATION.out.protospacer,
        PIPELINE_INITIALISATION.out.template
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_CRISPRSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
