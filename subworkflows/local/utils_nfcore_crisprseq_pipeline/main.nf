//
// Subworkflow with functionality specific to the nf-core/crisprseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    reads_targeted   = Channel.empty()
    reads_screening  = Channel.empty()
    fastqc_screening = Channel.empty()
    reference        = Channel.empty()
    protospacer      = Channel.empty()
    template         = Channel.empty()
    versions         = Channel.empty()

    //
    // Create channel from input file provided through params.input
    //
    if(params.input) {
        Channel
            .fromSamplesheet("input")
            .multiMap {
                meta, fastq_1, fastq_2, reference, protospacer, template ->
                    if (fastq_2) {
                        files = [ fastq_1, fastq_2 ]
                    } else {
                        files = [ fastq_1 ]
                    }
                    reads_targeted: [ meta.id, meta - meta.subMap('condition') + [ single_end : fastq_2 ? false : true, self_reference : reference ? false : true, template : template ? true : false ], files ]
                    reads_screening:[ meta + [ single_end:fastq_2?false:true ], files ]
                    reference:      [meta - meta.subMap('condition') + [ single_end : fastq_2 ? false : true, self_reference : reference ? false : true, template : template ? true : false ], reference]
                    protospacer:    [meta - meta.subMap('condition') + [ single_end : fastq_2 ? false : true, self_reference : reference ? false : true, template : template ? true : false ], protospacer]
                    template:       [meta - meta.subMap('condition') + [ single_end : fastq_2 ? false : true, self_reference : reference ? false : true, template : template ? true : false ], template]
            }
            .set { ch_input }

        //
        // Validate input samplesheet
        //
        ch_input.reads_targeted
            .groupTuple()
            .map {
                validateInputSamplesheet(it)
            }
            .set { reads_targeted }

        fastqc_screening = ch_input.reads_screening
        reference = ch_input.reference
        protospacer = ch_input.protospacer
        template = ch_input.template
    } else {
        ch_input = Channel.empty()
    }

    emit:
    reads_targeted
    fastqc_screening
    reference
    protospacer
    template
    versions = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE CHANNELS - SCREENING
========================================================================================
*/

workflow INITIALISATION_CHANNEL_CREATION_SCREENING {

    take:

    main:

    ch_library = Channel.empty()
    ch_crisprcleanr = Channel.empty()
    ch_design = Channel.empty()

    // Library
    if (params.library) {
        ch_library = Channel.fromPath(params.library)
    }

    // Crisprcleanr
    if (params.crisprcleanr) {
        if(params.crisprcleanr.endsWith(".csv")) {
            ch_crisprcleanr = Channel.fromPath(params.crisprcleanr)
        } else {
            ch_crisprcleanr = Channel.value(params.crisprcleanr)
        }
    }

    // MLE design matrix
    if(params.mle_design_matrix) {
        ch_design = Channel.fromPath(params.mle_design_matrix)
    }


    ch_biogrid = Channel.fromPath("$projectDir/assets/biogrid_hgncid_noduplicate_dropna.csv", checkIfExists: true)
    ch_hgnc = Channel.fromPath("$projectDir/assets/hgnc_complete_set.txt", checkIfExists: true)

    if(params.mle_control_sgrna) {
        ch_mle_control_sgrna = Channel.fromPath(params.mle_control_sgrna)
    } else {
        ch_mle_control_sgrna = []
    }

    emit:
    library = ch_library                     // channel: library file
    crisprcleanr = ch_crisprcleanr           // channel: crisprcleanr file or value
    design = ch_design                       // channel: design matrix file
    mle_control_sgrna = ch_mle_control_sgrna // channel: negative control sgRNA for MAGeCK MLE
    biogrid = ch_biogrid            // channel: biogrid
    hgnc = ch_hgnc                  // channel: hgnc

}

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE CHANNELS - TARGETED
========================================================================================
*/

workflow INITIALISATION_CHANNEL_CREATION_TARGETED {

    take:
    input_reads
    input_reference
    input_template
    input_protospacer

    main:

    //
    // Separate samples by the ones containing all reads in one file or the ones with many files to be concatenated
    //
    input_reads
    .groupTuple()
    .branch {
        meta, fastqs ->
            single  : fastqs.size() == 1
                return [ meta, fastqs.flatten() ]
            multiple: fastqs.size() > 1
                return [ meta, fastqs.flatten() ]
    }
    .set { ch_fastq }

    //
    // Add reference sequences to file
    //
    input_reference
    .tap{ meta_reference }
    .filter{ meta, sequence -> sequence instanceof String }
    .collectFile() { meta, reference ->
        [ "${meta.id}_reference.fasta", ">${meta.id}\n${reference}\n" ] // Write each reference sequence to a file
    }
    .map{ new_file ->
        [new_file.baseName.split("_reference")[0], new_file] // create a channel with the meta.id and the new file
    }
    .join(meta_reference
        .map{ meta, reference ->
            [meta.id, meta] // Join the channel by meta.id with the meta map
        }
    )
    .map{ metaid, new_file, meta ->
        [meta, new_file] // Obtain the final channel with meta map and the new file
    }
    .set{ ch_seq_reference }


    //
    // Add template sequences to file
    //
    input_template
    .tap{ meta_template }
    .filter{ meta, sequence -> sequence instanceof String }
    .collectFile() { meta, template ->
        [ "${meta.id}_template.fasta", ">${meta.id}\n${template}\n" ] // Write each template sequence to a file
    }
    .map{ new_file ->
        [new_file.baseName.split("_template")[0], new_file] // create a channel with the meta.id and the new file
    }
    .join(meta_template
        .map{ meta, template ->
            [meta.id, meta] // Join the channel by meta.id with the meta map
        }
    )
    .map{ metaid, new_file, meta ->
        [meta, new_file] // Obtain the final channel with meta map and the new file
    }
    .set{ ch_seq_template }


    // Join channels with reference and protospacer
    // to channel: [ meta, reference, protospacer]
    if (!params.reference_fasta && !params.protospacer) {
        ch_seq_reference
            .join(input_protospacer)
            .set{ reference_protospacer }
    } else if (!params.reference_fasta) {
        // If a protospacer was provided through the --protospacer param instead of the samplesheet
        ch_protospacer = Channel.of(params.protospacer)
        ch_seq_reference
            .combine(ch_protospacer)
            .set{ reference_protospacer }
    } else if (!params.protospacer) {
        // If a reference was provided through a fasta file or igenomes instead of the samplesheet
        ch_reference = Channel.fromPath(params.reference_fasta)
        input_protospacer
            .combine(ch_reference)
            .map{ meta, protospacer, reference -> [ meta, reference, protospacer ]} // Change the order of the channel
            .set{ reference_protospacer }
    } else {
        ch_reference = Channel.fromPath(params.reference_fasta)
        ch_protospacer = Channel.of(params.protospacer)
        input_reads
            .combine(ch_reference)
            .combine(ch_protospacer)
            .map{ meta, fastqs, reference, protospacer -> [ meta, reference, protospacer ]} // Don't add fastqs to the channel
            .set{ reference_protospacer }
    }

    emit:
    fastq_multiple = ch_fastq.multiple // [ meta, fastqs ] // Channel with the samples with multiple files
    fastq_single = ch_fastq.single // [ meta, fastqs ] // Channel with the samples with only one file
    template = ch_seq_template // [ meta, template ] // Channel with the template sequences
    reference_protospacer = reference_protospacer // [ meta, reference, protospacer] // Channel with the reference and protospacer sequences

}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample contain a reference or not
    def reference_ok = metas.collect{ it.self_reference }.unique().size == 1
    if (!reference_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must all contain a reference or not: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample contain a template or not
    def template_ok = metas.collect{ it.template }.unique().size == 1
    if (!template_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must all contain a template or not: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

def validateParametersScreening() {
    if(params.rra && params.mle_design_matrix) {
        warning "mle_design_matrix will only be used for the MAGeCK MLE computations"
    }

    if(params.fasta && params.count_table) {
        error "Please provide either a fasta file or a count_table"
    }

    if(params.fasta && !params.library) {
        error "Please provide a fasta file and the library file"
    }

    if(params.day0_label && params.mle_design_matrix) {
        warning "MAGeCK MLE module will be run twice, once with the design matrix and once with day0-label"
    }

    if(params.rra && params.mle_design_matrix) {
        warning "mle_design_matrix will only be used for the MAGeCK MLE computations"
    }

    if(params.rra && !params.contrasts) {
        error "Please also provide the contrasts table to compare the samples for MAGeCK RRA"
    }
}
