/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCrisprseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                      } from '../subworkflows/local/input_check'

//
// MODULE
//
include { FIND_ADAPTERS                    } from '../modules/local/find_adapters'
include { CUTADAPT                         } from '../modules/local/cutadapt_custom'
include { EXTRACT_UMIS                     } from '../modules/local/extract_umis'
include { SEQ_TO_FILE as SEQ_TO_FILE_REF   } from '../modules/local/seq_to_file'
include { SEQ_TO_FILE as SEQ_TO_FILE_TEMPL } from '../modules/local/seq_to_file'
include { ORIENT_REFERENCE                 } from '../modules/local/orient_reference'
include { CIGAR_PARSER                     } from '../modules/local/cigar_parser'
include { ALIGNMENT_SUMMARY                } from '../modules/local/alignment_summary'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { PEAR                        } from '../modules/nf-core/modules/pear/main'
include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { SEQTK_SEQ                   } from '../modules/nf-core/modules/seqtk/seq/main'
include { BOWTIE2_ALIGN               } from '../modules/nf-core/modules/bowtie2/align/main'
include { BOWTIE2_BUILD               } from '../modules/nf-core/modules/bowtie2/build/main'
include { BWA_MEM                     } from '../modules/nf-core/modules/bwa/mem/main'
include { BWA_INDEX                   } from '../modules/nf-core/modules/bwa/index/main'
include { MINIMAP2_ALIGN              } from '../modules/nf-core/modules/minimap2/align/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CRISPRSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            if (meta_clone.id.split('_').size() > 1) {
                meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            }
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    // Separate samples by the ones containing all reads in one file or the ones with many files to be concatenated
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    // Add reference sequences to file
    SEQ_TO_FILE_REF (
        INPUT_CHECK.out.reference,
        "reference"
    )

    //
    // Add template sequences to file
    //
    SEQ_TO_FILE_TEMPL (
        INPUT_CHECK.out.template,
        "template"
    )

    // Join channels with reference and protospacer
    // to channel: [ meta, reference, protospacer]
    SEQ_TO_FILE_REF.out.file
    .join(INPUT_CHECK.out.protospacer, by: 0)
    .set{ reference_protospacer }

    //
    // MODULE: Prepare reference sequence
    //
    ORIENT_REFERENCE (
        reference_protospacer
    )
    ch_versions = ch_versions.mix(ORIENT_REFERENCE.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .groupTuple(by: [0])
    .mix(ch_fastq.single)
    // Separate samples by paired-end or single-end
    .branch {
        meta, fastq ->
            single: fastq.size() == 1
                return [ meta, fastq ]
            paired: fastq.size() > 1
                return [ meta, fastq ]
    }
    .set { ch_cat_fastq }
    // comment version as there is an error (cat version is "cat: cat:")
    //ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // MODULE: Merge paired-end reads
    //
    PEAR (
        ch_cat_fastq.paired
    )
    .assembled
    .mix( ch_cat_fastq.single )
    .set { ch_pear_fastq }
    ch_versions = ch_versions.mix(PEAR.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_pear_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    FIND_ADAPTERS (
        FASTQC.out.zip
    )
    .adapters
    .join(ch_pear_fastq)
    .groupTuple(by: [0])
    // Separate samples by containing overrepresented sequences or not
    .branch {
        meta, adapter_lines, adapter_seqs, reads ->
            no_adapters: Integer.parseInt(adapter_lines.toString().substring(2,3)) < 6
                return [ meta, reads ]
            adapters   : Integer.parseInt(adapter_lines.toString().substring(2,3)) >= 6
                return [ meta, adapter_seqs, reads ]
    }
    .set { ch_adapter_seqs }

    //
    // MODULE: Trim adapter sequences
    //
    CUTADAPT (
        ch_adapter_seqs.adapters
    )

    ch_adapter_seqs.no_adapters
    .mix(CUTADAPT.out.reads)
    .groupTuple(by: [0])
    .map {
        meta, fastq ->
                return [ meta, fastq.flatten() ]
    }
    .set{ ch_trimmed }

    //
    // MODULE: Mask (convert to Ns) bases with quality lower than 20 and remove sequences shorter than 80
    //
    SEQTK_SEQ (
        ch_trimmed
    )

    /*
    Remove this step by now as I don't have test data
    //
    // MODULE: Extract UMI sequences
    //
    EXTRACT_UMIS (
        SEQTK_SEQ.out.fastx
    )

    To implement for UMIs

    vsearch
    get_ubs
    top_read
    polishing: minimap2, racon
    consensus
    join_reads
    fa2fq

    */

    //ch_mapped_bam = Channel.empty()

    if (params.aligner == "minimap2") {
        MINIMAP2_ALIGN (
            SEQTK_SEQ.out.fastx,
            ORIENT_REFERENCE.out.reference.map { it[1] },
            true,
            false,
            true
        )
        ch_mapped_bam = MINIMAP2_ALIGN.out.bam
        ch_mapped_bam_summary = MINIMAP2_ALIGN.out.bam
    }

    if (params.aligner == "bwa") {
        BWA_INDEX (
            ORIENT_REFERENCE.out.reference.map { it[1] }
        )
        BWA_MEM (
            SEQTK_SEQ.out.fastx,
            BWA_INDEX.out.index,
            true
        )
        ch_mapped_bam = BWA_MEM.out.bam
        ch_mapped_bam_summary = BWA_MEM.out.bam
    }

    if (params.aligner == "bowtie2") {
        BOWTIE2_BUILD (
            ORIENT_REFERENCE.out.reference.map { it[1] }
        )
        BOWTIE2_ALIGN (
            SEQTK_SEQ.out.fastx,
            BOWTIE2_BUILD.out.index,
            false,
            true
        )
        ch_mapped_bam = BOWTIE2_ALIGN.out.bam
        ch_mapped_bam_summary = BOWTIE2_ALIGN.out.bam
    }

    ALIGNMENT_SUMMARY (
        ch_mapped_bam_summary //.join(summary)
    )

    CIGAR_PARSER (
        ch_mapped_bam
        .join(ORIENT_REFERENCE.out.reference)
        .join(INPUT_CHECK.out.protospacer)
        .join(SEQ_TO_FILE_TEMPL.out.file)
        .join(ALIGNMENT_SUMMARY.out.summary)
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCrisprseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
