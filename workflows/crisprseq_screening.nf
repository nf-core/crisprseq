/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCrisprseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.multiqc_config, params.fasta, params.library, params.mle_design_matrix ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (!params.count_table) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.library) { ch_library = file(params.library) }
if (params.crisprcleanr) { ch_crisprcleanr= Channel.value(params.crisprcleanr) }

if(params.mle_design_matrix) {
    Channel.fromPath(params.mle_design_matrix)
        .set { ch_design }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK_SCREENING } from '../subworkflows/local/input_check_screening'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { MAGECK_COUNT                } from '../modules/nf-core/mageck/count/main'
include { MAGECK_MLE                  } from '../modules/nf-core/mageck/mle/main'
include { MAGECK_TEST                 } from '../modules/nf-core/mageck/test/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CRISPRCLEANR_NORMALIZE      } from '../modules/nf-core/crisprcleanr/normalize/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CRISPRSEQ_SCREENING {

    ch_versions = Channel.empty()

    if(!params.count_table){
        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //
        INPUT_CHECK_SCREENING (
            ch_input
        )
        ch_versions = ch_versions.mix(INPUT_CHECK_SCREENING.out.versions)

        //
        // MODULE: Run FastQC
        //
        FASTQC (
            INPUT_CHECK_SCREENING.out.reads
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        INPUT_CHECK_SCREENING.out.reads
        .map { meta, fastq ->
            [meta.condition, fastq]
        }
        .reduce { a, b ->
            ["${a[0]},${b[0]}", a[1] + b[1]]
        }
        .map { condition, fastqs ->
            [[id: condition], fastqs]
        }
        .set { joined }

        //
        // MODULE: Run mageck count
        //
        MAGECK_COUNT (
            joined,
            ch_library
        )

        ch_versions = ch_versions.mix(MAGECK_COUNT.out.versions.first())


        MAGECK_COUNT.out.count.map {
        it -> it[1]
        }.set { ch_counts }

    } else {
        Channel.fromPath(params.count_table)
        .set { ch_counts }
    }


    if(params.crisprcleanr) {
        ch_crispr_normalize = Channel.of([id: "count_table_normalize"])
        CRISPRCLEANR_NORMALIZE(
            ch_crispr_normalize.concat(ch_counts,ch_crisprcleanr).collect(),
            params.min_reads,
            params.min_targeted_genes
        )

        CRISPRCLEANR_NORMALIZE.out.norm_count_file.map {
            it -> it[1]
        }.set { ch_counts }
    }

    if(params.rra_contrasts) {
        Channel.fromPath(params.rra_contrasts)
            .splitCsv(header:true, sep:',' )
            .set { ch_contrasts }
        counts = ch_contrasts.combine(ch_counts)

        MAGECK_TEST (
            counts
        )
    }

    if(params.mle_design_matrix) {
        ch_mle = ch_counts.combine(ch_design)
        ch_mle.map {
            it -> [[id: it[1].getBaseName()], it[0], it[1]]
        }.set { ch_designed_mle }

        MAGECK_MLE (
            ch_designed_mle
        )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCrisprseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCrisprseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    if(!params.count_table) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    } else {
        ch_multiqc_files = channel.empty()
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
