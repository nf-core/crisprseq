/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Local modules
include { BAGEL2_FC                                    } from '../modules/local/bagel2/fc'
include { BAGEL2_BF                                    } from '../modules/local/bagel2/bf'
include { BAGEL2_PR                                    } from '../modules/local/bagel2/pr'
include { BAGEL2_GRAPH                                 } from '../modules/local/bagel2/graph'
include { MATRICESCREATION                             } from '../modules/local/matricescreation'
include { MAGECK_FLUTEMLE                              } from '../modules/local/mageck/flutemle'
include { MAGECK_FLUTEMLE as MAGECK_FLUTEMLE_CONTRASTS } from '../modules/local/mageck/flutemle'
include { MAGECK_FLUTEMLE as MAGECK_FLUTEMLE_DAY0      } from '../modules/local/mageck/flutemle'
include { VENNDIAGRAM                                  } from '../modules/local/venndiagram'
include { GPT_PREPARE_BAGEL2_QUERY                     } from '../modules/local/gpt_prepare_bagel2_query'
include { GPT_PREPARE_DRUGZ_QUERY                      } from '../modules/local/gpt_prepare_drugz_query'
include { GPT_PREPARE_MLE_QUERY                        } from '../modules/local/gpt_prepare_mle_query'

// nf-core modules
include { FASTQC                                       } from '../modules/nf-core/fastqc/main'
include { CUTADAPT as CUTADAPT_THREE_PRIME             } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_FIVE_PRIME              } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                                      } from '../modules/nf-core/multiqc/main'
include { MAGECK_COUNT                                 } from '../modules/nf-core/mageck/count/main'
include { MAGECK_MLE                                   } from '../modules/nf-core/mageck/mle/main'
include { MAGECK_TEST                                  } from '../modules/nf-core/mageck/test/main'
include { MAGECK_GRAPHRRA                              } from '../modules/local/mageck/graphrra'
include { CRISPRCLEANR_NORMALIZE                       } from '../modules/nf-core/crisprcleanr/normalize/main'
include { MAGECK_MLE as MAGECK_MLE_MATRIX              } from '../modules/nf-core/mageck/mle/main'
include { MAGECK_MLE as MAGECK_MLE_DAY0                } from '../modules/nf-core/mageck/mle/main'
include { BOWTIE2_BUILD                                } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                                } from '../modules/nf-core/bowtie2/align/main'
// Local subworkflows
include { INITIALISATION_CHANNEL_CREATION_SCREENING    } from '../subworkflows/local/utils_nfcore_crisprseq_pipeline'
// Functions
include { paramsSummaryMap                             } from 'plugin/nf-validation'
include { gptPromptForText                             } from 'plugin/nf-gpt'
include { paramsSummaryMultiqc                         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                       } from '../subworkflows/local/utils_nfcore_crisprseq_pipeline'
include { validateParametersScreening                  } from '../subworkflows/local/utils_nfcore_crisprseq_pipeline'
include { DRUGZ                                        } from '../modules/local/drugz'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CRISPRSEQ_SCREENING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    // Set screening parameters and channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Validate parameters specific to the screening subworkflow
    validateParametersScreening()

    //
    // Initialise channels
    //
    INITIALISATION_CHANNEL_CREATION_SCREENING()

    if(!params.count_table){
        ch_input = ch_samplesheet

        //
        // MODULE: Run FastQC
        //
        FASTQC (
            ch_input
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        //set adapter seq to null to make it compatible with crispr targeted
        ch_cutadapt = ch_input.combine(Channel.value([[]]))
        if(params.five_prime_adapter) {
            CUTADAPT_FIVE_PRIME(
                ch_cutadapt
            )
            CUTADAPT_FIVE_PRIME.out.reads.combine(Channel.value([[]])).set { ch_cutadapt }
            ch_cutadapt.map{ meta, fastq, proto  ->
                meta.id = "${meta.id}_trim"
                [meta, fastq, proto]
            }.set { ch_cutadapt }

            ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_FIVE_PRIME.out.log.collect{it[1]})
            ch_versions = ch_versions.mix(CUTADAPT_FIVE_PRIME.out.versions)
        }

        if(params.three_prime_adapter) {
            CUTADAPT_THREE_PRIME(
                ch_cutadapt
            )
            ch_cutadapt = CUTADAPT_THREE_PRIME.out.reads.combine(Channel.value([[]]))
            ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_THREE_PRIME.out.log.collect{it[1]})
            ch_versions = ch_versions.mix(CUTADAPT_THREE_PRIME.out.versions)
        }


        if(params.five_prime_adapter || params.three_prime_adapter) {
            ch_cutadapt
            .map{ meta, fastq, empty  ->
                [meta, fastq]
            }
            .set { ch_input }
        }

        if(params.fasta){
            Channel.of("fasta")
                .combine(Channel.fromPath(params.fasta))
                .set{ ch_fasta }

            BOWTIE2_BUILD(ch_fasta)
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

            BOWTIE2_ALIGN (
            ch_input,
            BOWTIE2_BUILD.out.index,
            false,
            false
            )

            ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)


            BOWTIE2_ALIGN.out.aligned.map{ meta, bam ->
                [meta, [bam]]
            }.set{ch_input}
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
            INITIALISATION_CHANNEL_CREATION_SCREENING.out.library
        )

        ch_versions = ch_versions.mix(MAGECK_COUNT.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(MAGECK_COUNT.out.summary.collect{it[1]})

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
                '',
                INITIALISATION_CHANNEL_CREATION_SCREENING.out.crisprcleanr,
                params.min_reads,
                params.min_targeted_genes
        ) } else
        {
            ch_crispr_normalize = Channel.of([id: "count_table_normalize"]).concat(ch_counts)
            CRISPRCLEANR_NORMALIZE(
                ch_crispr_normalize.collect(),
                INITIALISATION_CHANNEL_CREATION_SCREENING.out.crisprcleanr,
                [],
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
            .set { ch_contrasts }
    counts = ch_contrasts.combine(ch_counts)


    //Define non essential and essential genes channels for bagel2
    ch_bagel_reference_essentials= Channel.fromPath(params.bagel_reference_essentials).first()
    ch_bagel_reference_nonessentials= Channel.fromPath(params.bagel_reference_nonessentials).first()

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

    if((params.mle_design_matrix) || (params.contrasts && !params.rra) || (params.day0_label)) {
        //if the user only wants to run mle through their own design matrices
        if(params.mle_design_matrix) {
            INITIALISATION_CHANNEL_CREATION_SCREENING.out.design.map {
                it -> [[id: it.getBaseName()], it]
                }.set { ch_designed_mle }

            ch_mle = ch_designed_mle.combine(ch_counts)
            MAGECK_MLE_MATRIX (ch_mle, INITIALISATION_CHANNEL_CREATION_SCREENING.out.mle_control_sgrna)
            ch_versions = ch_versions.mix(MAGECK_MLE_MATRIX.out.versions)
            MAGECK_FLUTEMLE(MAGECK_MLE_MATRIX.out.gene_summary)
            ch_versions = ch_versions.mix(MAGECK_FLUTEMLE.out.versions)
        }
        //if the user specified a contrast file
        if(params.contrasts) {
            MATRICESCREATION(ch_contrasts)
            ch_mle = MATRICESCREATION.out.design_matrix.combine(ch_counts)
            MAGECK_MLE (ch_mle, INITIALISATION_CHANNEL_CREATION_SCREENING.out.mle_control_sgrna)
            ch_versions = ch_versions.mix(MAGECK_MLE.out.versions)
            MAGECK_FLUTEMLE_CONTRASTS(MAGECK_MLE.out.gene_summary)
            ch_versions = ch_versions.mix(MAGECK_FLUTEMLE_CONTRASTS.out.versions)
            ch_venndiagram = BAGEL2_PR.out.pr.join(MAGECK_MLE.out.gene_summary)
            VENNDIAGRAM(ch_venndiagram)
            ch_versions = ch_versions.mix(VENNDIAGRAM.out.versions)
        }
        if(params.day0_label) {
            ch_mle = Channel.of([id: "day0"]).merge(Channel.of([[]])).merge(ch_counts)
            MAGECK_MLE_DAY0 (ch_mle, INITIALISATION_CHANNEL_CREATION_SCREENING.out.mle_control_sgrna)
            ch_versions = ch_versions.mix(MAGECK_MLE_DAY0.out.versions)
            MAGECK_FLUTEMLE_DAY0(MAGECK_MLE_DAY0.out.gene_summary)
            ch_versions = ch_versions.mix(MAGECK_FLUTEMLE_DAY0.out.versions)
        }
    }

    // Launch module drugZ
    if(params.drugz) {
        Channel.fromPath(params.drugz)
                .splitCsv(header:true, sep:';' )
                .set { ch_drugz }

        counts = ch_drugz.combine(ch_counts)
        DRUGZ (
            counts
            )
        ch_versions = ch_versions.mix(DRUGZ.out.versions)
    }

    //
    // Calling of nf-gpt plugin on drugZ or MAGeCK mle
    //
    if(params.gpt_interpretation.contains("drugz")) {
        def gpt_drugz_data = DRUGZ.out.per_gene_results.map { meta, genes -> genes }
        GPT_PREPARE_DRUGZ_QUERY(
            gpt_drugz_data,
            params.gpt_drugz_gene_amount,
            params.gpt_drugz_question
        )

        GPT_PREPARE_DRUGZ_QUERY.out.query.map {
            it -> it.text
        }
        .collect()
        .flatMap { it -> gptPromptForText(it[0]) }
        .collectFile( name: "${params.outdir}/gpt/gpt_drugz_output.txt", newLine: true, sort: false )
    }
    if(params.gpt_interpretation.contains("mle")) {
        def gpt_mle_data = MAGECK_MLE.out.gene_summary.map { meta, genes -> genes }
        GPT_PREPARE_MLE_QUERY(
            gpt_mle_data,
            params.gpt_mle_gene_amount,
            params.gpt_mle_question
        )

        GPT_PREPARE_MLE_QUERY.out.query.map {
            it -> it.text
        }
        .collect()
        .flatMap { it -> gptPromptForText(it[0]) }
        .collectFile( name: "${params.outdir}/gpt/gpt_mle_output.txt", newLine: true, sort: false )
    }
    if(params.gpt_interpretation.contains("bagel2")) {
        def gpt_bagel2_data = BAGEL2_BF.out.bf.map { meta, genes -> genes }
        GPT_PREPARE_BAGEL2_QUERY(
            gpt_bagel2_data,
            params.gpt_bagel2_gene_amount,
            params.gpt_bagel_question
        )

        GPT_PREPARE_BAGEL2_QUERY.out.query.map {
            it -> it.text
        }
        .collect()
        .flatMap { it -> gptPromptForText(it[0]) }
        .collectFile( name: "${params.outdir}/gpt/gpt_bagel2_output.txt", newLine: true, sort: false )
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
