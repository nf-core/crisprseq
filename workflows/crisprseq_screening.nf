/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowCrisprseq.initialise(params, log)

// Set screening parameters and channels
if (params.library) { ch_library = file(params.library) }
if (params.crisprcleanr) {
    if(params.crisprcleanr.endsWith(".csv")) {
        ch_crisprcleanr_file = file(params.crisprcleanr)
    } else {
        ch_crisprcleanr_value = Channel.value(params.crisprcleanr)
        }
    }

if(params.mle_design_matrix) {
    Channel.fromPath(params.mle_design_matrix)
        .set { ch_design }
}

if(params.rra && params.mle_design_matrix) {
    warning "mle_design_matrix will only be used for the MAGeCK MLE computations"
    }

if(params.rra && !params.contrasts) {
    error "Please also provide the contrasts table to compare the samples for MAGeCK RRA"
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo )   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                    } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { MAGECK_COUNT                } from '../modules/nf-core/mageck/count/main'
include { MAGECK_MLE                  } from '../modules/nf-core/mageck/mle/main'
include { MAGECK_TEST                 } from '../modules/nf-core/mageck/test/main'
include { MAGECK_GRAPHRRA             } from '../modules/local/mageck/graphrra'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CRISPRCLEANR_NORMALIZE      } from '../modules/nf-core/crisprcleanr/normalize/main'
include { BAGEL2_FC                   } from '../modules/local/bagel2/fc'
include { BAGEL2_BF                   } from '../modules/local/bagel2/bf'
include { BAGEL2_PR                   } from '../modules/local/bagel2/pr'
include { BAGEL2_GRAPH                } from '../modules/local/bagel2/graph'
include { MATRICESCREATION            } from '../modules/local/matricescreation'

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
        // Create input channel from input file provided through params.input
        //
        Channel.fromSamplesheet("input")
        .map{ meta, fastq_1, fastq_2, x, y, z ->
            // x (reference), y (protospacer), and z (template) are part of the targeted workflows and we don't need them
        return   [ meta + [ single_end:fastq_2?false:true ], fastq_2?[ fastq_1, fastq_2 ]:[ fastq_1 ] ]        }
        .set { ch_input }


        //
        // MODULE: Run FastQC
        //
        FASTQC (
            ch_input
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())


        ch_input_cutadapt = ch_input.combine(Channel.value([[]]))

    if(params.cutadapt) {
        CUTADAPT(
            ch_input_cutadapt
        )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)

        CUTADAPT.out.reads
        .map{ meta, fastq  ->
            [meta, [fastq]]
        }
        .set { ch_input }
        }

        // this is to concatenate everything for mageck count

        ch_input
        .map { meta, fastqs  ->
            if(fastqs.size() == 1){
                [meta.condition, [fastqs[0]], meta.single_end, []]
            } else {
                [meta.condition, [fastqs[0]], meta.single_end, [fastqs[1]]]
            }
        }
        // if one element is paired-end and the other single-end throw an error
        // otherwise just concatenate the conditions and the fastqs
        .reduce { a, b ->
            if(a[2] != b[2] ) {
                error "Your samplesheet contains a mix of single-end and paired-end data. This is not supported."
            }
            return ["${a[0]},${b[0]}", a[1] + b[1], b[2] ,a[3] + b[3]]
        }
        .map { condition, fastqs_1, single_end, fastqs_2 ->
            [[id: condition, single_end: single_end], fastqs_1, fastqs_2]
        }
        .last()
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
        ch_crispr_normalize = Channel.of([id: "count_table_normalize"]).concat(ch_counts)

        if(params.crisprcleanr.endsWith(".csv")) {
            CRISPRCLEANR_NORMALIZE(
                ch_crispr_normalize.collect(),
                Channel.empty(),
                ch_crisprcleanr_file,
                params.min_reads,
                params.min_targeted_genes
        ) } else
        {
            ch_crispr_normalize = Channel.of([id: "count_table_normalize"]).concat(ch_counts)
            CRISPRCLEANR_NORMALIZE(
                ch_crispr_normalize.collect(),
                ch_crisprcleanr_value,
                Channel.empty(),
                params.min_reads,
                params.min_targeted_genes)
        }

        ch_versions = ch_versions.mix(CRISPRCLEANR_NORMALIZE.out.versions)


        CRISPRCLEANR_NORMALIZE.out.norm_count_file.map {
            it -> it[1]
        }.set { ch_counts }
    }

    if(params.rra) {
        Channel.fromPath(params.contrasts)
            .splitCsv(header:true, sep:';' )
            .set { ch_contrasts }
        counts = ch_contrasts.combine(ch_counts)

        MAGECK_TEST (
            counts
        )

        ch_versions = ch_versions.mix(MAGECK_TEST.out.versions)

        MAGECK_GRAPHRRA (
            MAGECK_TEST.out.gene_summary
        )
        ch_versions = ch_versions.mix(MAGECK_GRAPHRRA.out.versions)
    }

    if(params.contrasts) {
        Channel.fromPath(params.contrasts)
            .splitCsv(header:true, sep:';' )
            .set { ch_bagel }
    counts = ch_bagel.combine(ch_counts)

    //Define non essential and essential genes channels for bagel2
    ch_bagel_reference_essentials= Channel.value(params.bagel_reference_essentials)
    ch_bagel_reference_nonessentials= Channel.value(params.bagel_reference_nonessentials)

    BAGEL2_FC (
            counts
        )
    ch_versions = ch_versions.mix(BAGEL2_FC.out.versions)

    BAGEL2_BF (
        BAGEL2_FC.out.foldchange,
        ch_bagel_reference_essentials,
        ch_bagel_reference_nonessentials
    )

    ch_versions = ch_versions.mix(BAGEL2_BF.out.versions)


    ch_bagel_pr = BAGEL2_BF.out.bf.combine(ch_bagel_reference_essentials)
                                        .combine(ch_bagel_reference_nonessentials)

    BAGEL2_PR (
        ch_bagel_pr
    )

    ch_versions = ch_versions.mix(BAGEL2_PR.out.versions)

    BAGEL2_GRAPH (
        BAGEL2_PR.out.pr
    )

    ch_versions = ch_versions.mix(BAGEL2_GRAPH.out.versions)

    }

    if((params.mle_design_matrix) || (params.contrasts && !params.rra)) {
        if(params.mle_design_matrix) {
            ch_mle = ch_counts.combine(ch_design)
            }
        if(params.contrasts) {
            MATRICESCREATION(params.contrasts)
            ch_mle = ch_counts.combine(MATRICESCREATION.out.design_matrix)
        }
        ch_mle.map {
            it -> [[id: it[1].getBaseName()], it[0], it[1]]
        }.set { ch_designed_mle }

        MAGECK_MLE (
            ch_designed_mle
        )

        ch_versions = ch_versions.mix(MAGECK_MLE.out.versions)

    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCrisprseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCrisprseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
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
    NfcoreTemplate.dump_parameters(workflow, params)
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
