/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCrisprseq.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
include { INPUT_CHECK                      } from '../subworkflows/local/input_check'

//
// MODULE
//
include { FIND_ADAPTERS                    } from '../modules/local/find_adapters'
include { EXTRACT_UMIS                     } from '../modules/local/extract_umis'
include { SEQ_TO_FILE as SEQ_TO_FILE_REF   } from '../modules/local/seq_to_file'
include { SEQ_TO_FILE as SEQ_TO_FILE_TEMPL } from '../modules/local/seq_to_file'
include { ORIENT_REFERENCE                 } from '../modules/local/orient_reference'
include { CIGAR_PARSER                     } from '../modules/local/cigar_parser'
include { MERGING_SUMMARY                  } from '../modules/local/merging_summary'
include { CLUSTERING_SUMMARY               } from '../modules/local/clustering_summary'
include { ALIGNMENT_SUMMARY                } from '../modules/local/alignment_summary'
include { TEMPLATE_REFERENCE               } from '../modules/local/template_reference'
include { DUMMY_FINAL_UMI                  } from '../modules/local/dummy_final_umi'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PEAR                                          } from '../modules/nf-core/pear/main'
include { CAT_FASTQ                                     } from '../modules/nf-core/cat/fastq/main'
include { SEQTK_SEQ                                     } from '../modules/nf-core/seqtk/seq/main'
include { VSEARCH_CLUSTER                               } from '../modules/nf-core/vsearch/cluster/main'
include { BOWTIE2_ALIGN                                 } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD                                 } from '../modules/nf-core/bowtie2/build/main'
include { BWA_MEM                                       } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX                                     } from '../modules/nf-core/bwa/index/main'
include { MINIMAP2_ALIGN                                } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TEMPLATE     } from '../modules/nf-core/minimap2/align/main'
include { CUTADAPT                                      } from '../modules/nf-core/cutadapt/main'
include { SAMTOOLS_INDEX                                } from '../modules/nf-core/samtools/index/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE GROOVY FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def umi_to_sequence(cluster) {
    cluster.withReader { source ->
        String line
        while ( line=source.readLine() ) {
            if (line.startsWith(">")) {
                sequence = (line =~ /;seq=(.*$)/)[0][1]
                id = (line =~ /(>.*?);/)[0][1]
            }
        }
    }
    return id + "\n" + sequence
}


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
            [ meta - meta.subMap('id') + [id: meta.id.split('_')[0..-2].join('_')], fastq ]
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

    INPUT_CHECK.out.reference
    .map {
        meta, fastq ->
            [ meta - meta.subMap('id') + [id: meta.id.split('_')[0..-2].join('_')], fastq ]
    }

    //
    // MODULE: Add reference sequences to file
    //
    SEQ_TO_FILE_REF (
        INPUT_CHECK.out.reference
        .map {
            meta, fastq ->
                [ meta - meta.subMap('id') + [id: meta.id.split('_')[0..-2].join('_')], fastq ]
        },
        "reference"
    )
    ch_versions = ch_versions.mix(SEQ_TO_FILE_REF.out.versions)

    //
    // MODULE: Add template sequences to file
    //
    SEQ_TO_FILE_TEMPL (
        INPUT_CHECK.out.template
        .map {
            meta, fastq ->
                [ meta - meta.subMap('id') + [id: meta.id.split('_')[0..-2].join('_')], fastq ]
        },
        "template"
    )
    ch_versions = ch_versions.mix(SEQ_TO_FILE_TEMPL.out.versions)

    // Join channels with reference and protospacer
    // to channel: [ meta, reference, protospacer]
    SEQ_TO_FILE_REF.out.file
        .join(INPUT_CHECK.out.protospacer
            .map {
                meta, fastq ->
                    [ meta - meta.subMap('id') + [id: meta.id.split('_')[0..-2].join('_')], fastq ]
            },
            by: 0)
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
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

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

    //
    // MODULE: Find overrepresented sequences
    FIND_ADAPTERS (
        FASTQC.out.zip
    )
    .adapters
    .join(ch_pear_fastq)
    .groupTuple(by: [0])
    // Separate samples by containing overrepresented sequences or not
    .branch {
        meta, adapter_seqs, reads ->
            no_adapters: adapter_seqs[0].size() == 0
                return [ meta, reads[0] ]
            adapters   : adapter_seqs[0].size() > 0
                return [ meta, reads[0], adapter_seqs[0] ]
    }
    .set { ch_adapter_seqs }


    //
    // MODULE: Trim adapter sequences
    //
    CUTADAPT (
        ch_adapter_seqs.adapters
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

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
    ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions)


    //
    // MODULE: Summary of merged reads
    //
    MERGING_SUMMARY {
        ch_cat_fastq.paired
            .mix(ch_cat_fastq.single)
            .join(PEAR.out.assembled, remainder: true)
            .join(SEQTK_SEQ.out.fastx)
            .join(CUTADAPT.out.log)
            .map { meta, reads, assembled, masked, trimmed ->
                if (assembled == null) {
                    dummy = file('null')
                    return [ meta, reads, dummy, masked, trimmed ]
                } else {
                    return [ meta, reads, assembled, masked, trimmed ]
                }
            }
    }


    //
    // MODULE: Extract UMI sequences
    //
    EXTRACT_UMIS (
        SEQTK_SEQ.out.fastx
    )


    //
    // MODULE: Cluster UMIs
    //
    VSEARCH_CLUSTER (
        EXTRACT_UMIS.out.fasta
    )

    //  Obtain a file with UBS (UBI bin size) and UMI ID
    VSEARCH_CLUSTER.out.clusters
    .transpose()
    .collectFile( storeDir:params.outdir ) {
        it ->
            [ "${it[0].id}_ubs.txt", "${it[1].countFasta()}\t${it[1].baseName}\n" ]
    }

    // Branch the clusters into the ones containing only one sequence and the ones containing more than one sequences
    VSEARCH_CLUSTER.out.clusters
    .transpose()
    .branch{
        meta, cluster ->
            single: cluster.countFasta() == 1
                return [meta, cluster]
            cluster: cluster.countFasta() > 1
                return [meta, cluster]
    }
    .set{ ch_umi_bysize }

    // Get the correspondent fasta sequencences from single clusters
    ch_umi_bysize.single
    .map{ meta, cluster ->
        fasta_line = umi_to_sequence(cluster)
        [meta, cluster.baseName, fasta_line]
    }
    .collectFile() { meta, name, fasta ->
        [ "{$name}_consensus.fasta", fasta ]
    }
    .set{ ch_single_clusters_consensus }


    //
    // MODULE: Obtain longest read in cluster
    //
    //VSEARCH_SORT(

    //)

    /*
    The UMI clustering step is posponed until the next release, the steps to be implemented are listed below:


    Modules to implement:

    top_read
    polishing: minimap2, racon
    consensus
    join_reads
    fa2fq

    */

    // Dummy process simulating the last UMi clustering step to obtain clusters as fastq
    DUMMY_FINAL_UMI {
        SEQTK_SEQ.out.fastx
    }

    //
    // MODULE: Summary of clustered reads
    //
    CLUSTERING_SUMMARY (
        SEQTK_SEQ.out.fastx
            .join(MERGING_SUMMARY.out.summary)
    )


    //
    // MODULE: Mapping with Minimap2
    //
    if (params.aligner == "minimap2") {
        MINIMAP2_ALIGN (
            SEQTK_SEQ.out.fastx
                .join(ORIENT_REFERENCE.out.reference),
            true,
            false,
            true
        )
        ch_mapped_bam = MINIMAP2_ALIGN.out.bam
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    }

    //
    // MODULE: Mapping with BWA
    //
    if (params.aligner == "bwa") {
        BWA_INDEX (
            ORIENT_REFERENCE.out.reference.map { it[1] }
        )
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        BWA_MEM (
            SEQTK_SEQ.out.fastx,
            BWA_INDEX.out.index,
            true
        )
        ch_mapped_bam = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)
    }

    //
    // MODULE: Mapping with Bowtie2
    //
    if (params.aligner == "bowtie2") {
        BOWTIE2_BUILD (
            ORIENT_REFERENCE.out.reference.map { it[1] }
        )
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        BOWTIE2_ALIGN (
            SEQTK_SEQ.out.fastx,
            BOWTIE2_BUILD.out.index,
            false,
            true
        )
        ch_mapped_bam = BOWTIE2_ALIGN.out.bam
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
    }


    //
    // MODULE: Summary of mapped reads
    //
    ALIGNMENT_SUMMARY (
        ch_mapped_bam
            .join(CLUSTERING_SUMMARY.out.summary)
    )


    //
    // MODULE: Obtain .bam.bai files
    //
    SAMTOOLS_INDEX (
        ch_mapped_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //
    // MODULE: Obtain a new reference with the template modification
    //
    TEMPLATE_REFERENCE (
        ORIENT_REFERENCE.out.reference
            .join(SEQ_TO_FILE_TEMPL.out.file)
    )


    //
    // MODULE: Align new reference with the change led by template with original reference
    //
    MINIMAP2_ALIGN_TEMPLATE (
        TEMPLATE_REFERENCE.out.fasta
            .join(ORIENT_REFERENCE.out.reference),
        true,
        false,
        true
    )
    .bam
    .map {
        meta, bam ->
            if (bam.baseName.contains("template-align")) {
                return [ meta, bam ]
            } else {
                new_file = bam.parent / bam.baseName + "_template-align." + bam.extension
                bam.renameTo(new_file)
                return[ meta, new_file ]
            }
    }
    .set { ch_template_bam }
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_TEMPLATE.out.versions)

    ch_mapped_bam
        .join(SAMTOOLS_INDEX.out.bai)
        .join(ORIENT_REFERENCE.out.reference)
        .join(INPUT_CHECK.out.protospacer
            .map {
                meta, fastq ->
                    [ meta - meta.subMap('id') + [id: meta.id.split('_')[0..-2].join('_')], fastq ]
            }
        )
        .join(SEQ_TO_FILE_TEMPL.out.file, remainder: true)
        .join(ch_template_bam, remainder: true)
        .join(TEMPLATE_REFERENCE.out.fasta, remainder: true)
        .join(ALIGNMENT_SUMMARY.out.summary)
        .map { meta, reads, index, reference, protospacer, template, template_bam, reference_template, summary ->
            if (template == null) {
                template = file('null_t')
            }
            if (template_bam == null) {
                template_bam = file('null_b')
            }
            if (reference_template == null) {
                reference_template = file('null_r')
            }
            return [meta, reads, index, reference, protospacer, template, template_bam, reference_template, summary]
        }
        .set { ch_to_parse_cigar }


    //
    // MODULE: Parse cigar to find edits
    //
    CIGAR_PARSER (
        ch_to_parse_cigar
    )


    //
    // MODULE: Dump software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
