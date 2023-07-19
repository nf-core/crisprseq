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
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo )   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

//
// MODULE
//
include { FIND_ADAPTERS      } from '../modules/local/find_adapters'
include { EXTRACT_UMIS       } from '../modules/local/extract_umis'
include { ORIENT_REFERENCE   } from '../modules/local/orient_reference'
include { CIGAR_PARSER       } from '../modules/local/cigar_parser'
include { MERGING_SUMMARY    } from '../modules/local/merging_summary'
include { CLUSTERING_SUMMARY } from '../modules/local/clustering_summary'
include { ALIGNMENT_SUMMARY  } from '../modules/local/alignment_summary'
include { TEMPLATE_REFERENCE } from '../modules/local/template_reference'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS               } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PEAR                                      } from '../modules/nf-core/pear/main'
include { CAT_FASTQ                                 } from '../modules/nf-core/cat/fastq/main'
include { SEQTK_SEQ as SEQTK_SEQ_MASK               } from '../modules/nf-core/seqtk/seq/main'
include { SEQTK_SEQ as SEQTK_SEQ_FATOFQ             } from '../modules/nf-core/seqtk/seq/main'
include { VSEARCH_CLUSTER                           } from '../modules/nf-core/vsearch/cluster/main'
include { VSEARCH_SORT                              } from '../modules/nf-core/vsearch/sort/main'
include { RACON as RACON_1                          } from '../modules/nf-core/racon/main'
include { RACON as RACON_2                          } from '../modules/nf-core/racon/main'
include { BOWTIE2_ALIGN                             } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD                             } from '../modules/nf-core/bowtie2/build/main'
include { BWA_MEM                                   } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX                                 } from '../modules/nf-core/bwa/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ORIGINAL } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_UMI_1    } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_UMI_2    } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TEMPLATE } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_FAIDX                            } from '../modules/nf-core/samtools/faidx/main'
include { MINIMAP2_INDEX                            } from '../modules/nf-core/minimap2/index/main'
include { MEDAKA                                    } from '../modules/nf-core/medaka/main'
include { CUTADAPT                                  } from '../modules/nf-core/cutadapt/main'
include { SAMTOOLS_INDEX                            } from '../modules/nf-core/samtools/index/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE GROOVY FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def umi_to_sequence(cluster1) {
    String line1
    String sequences = ""
    String sequence1
    String id1
    cluster1.withReader {
        while ( line1=it.readLine() ) {
            if (line1.startsWith(">")) {
                sequence1 = (line1 =~ /;seq=(.*$)/)[0][1]
                id1 = (line1 =~ /(>.*?);/)[0][1]
                sequences = sequences + id1 + "\n" + sequence1 + "\n"
            }
        }
    }
    return sequences
}

def umi_to_sequence_centroid(cluster) {
    String line
    String sequence
    String id
    cluster.withReader {
        while ( line=it.readLine() ) {
            if (line.startsWith(">")) {
                sequence = (line =~ /;seq=(.*$)/)[0][1]
                id = (line =~ /(>.*?);/)[0][1]
                return id.replace(">", ">centroid_") + "\n" + sequence
            }
        }
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CRISPRSEQ_TARGETED {

    ch_versions = Channel.empty()

    //
    // Create input channel from input file provided through params.input
    //
    Channel.fromSamplesheet("input")
    .multiMap { meta, fastq_1, fastq_2, reference, protospacer, template ->
        // meta.condition is part of the screening workflow and we need to remove it
        reads:   [ meta.id, meta - meta.subMap('condition') + [ single_end:false, self_reference:reference?false:true, template:template?true:false ], fastq_2?[ fastq_1, fastq_2 ]:[ fastq_1 ] ]
        reference:   [meta - meta.subMap('condition') + [ single_end:true, self_reference:reference?false:true, template:template?true:false ], reference]
        protospacer: [meta - meta.subMap('condition') + [ single_end:true, self_reference:reference?false:true, template:template?true:false ], protospacer]
        template:    [meta - meta.subMap('condition') + [ single_end:true, self_reference:reference?false:true, template:template?true:false ], template]
    }
    .set { ch_input }

    ch_input
    .reads
    .groupTuple()
    .map {
        WorkflowCrisprseq.validateInput(it)
    }
    // Separate samples by the ones containing all reads in one file or the ones with many files to be concatenated
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
    ch_input.reference
    .tap{ meta_reference }
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
    ch_input.template
    .tap{ meta_template }
    .collectFile() { meta, template ->
        [ "${meta.if}_template.fasta", ">${meta.id}\n${template}\n" ] // Write each template sequence to a file
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
            .join(ch_input.protospacer)
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
        ch_input.protospacer
            .combine(ch_reference)
            .set{ reference_protospacer }
    } else {
        ch_reference = Channel.fromPath(params.reference_fasta)
        ch_protospacer = Channel.of(params.protospacer)
        ch_input.reads
            .combine(ch_reference)
            .combine(ch_protospacer)
            .set{ reference_protospacer }
    }


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

    ch_trimmed = Channel.empty()

    if (params.overrepresented) {
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
    } else {
        ch_trimmed = ch_pear_fastq
    }


    //
    // MODULE: Mask (convert to Ns) bases with quality lower than 20 and remove sequences shorter than 80
    //
    SEQTK_SEQ_MASK (
        ch_trimmed
    )
    ch_versions = ch_versions.mix(SEQTK_SEQ_MASK.out.versions)

    if (params.overrepresented) {
        ch_cat_fastq.paired
            .mix(ch_cat_fastq.single)
            .join(PEAR.out.assembled, remainder: true)
            .join(SEQTK_SEQ_MASK.out.fastx)
            .join(CUTADAPT.out.log)
            .set { ch_merging_summary_data }
    } else {
        ch_cat_fastq.paired
            .mix(ch_cat_fastq.single)
            .join(PEAR.out.assembled, remainder: true)
            .join(SEQTK_SEQ_MASK.out.fastx)
            .combine(Channel.value("null"))
            .map { meta, reads, assembled, masked, trimmed ->
                if (assembled == null) {
                    assembled = file('null_a')
                }
                if (trimmed == "null") {
                    trimmed = file('null_t')
                }
                return [ meta, reads, assembled, masked, trimmed ]
            }
            .set { ch_merging_summary_data }
    }

    //
    // MODULE: Summary of merged reads
    //
    MERGING_SUMMARY {
        ch_merging_summary_data
    }


    if (params.umi_clustering) {
        //
        // MODULE: Extract UMI sequences
        //
        EXTRACT_UMIS (
            SEQTK_SEQ_MASK.out.fastx
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
                single: cluster.countFasta() >= params.umi_bin_size && cluster.countFasta() == 1
                    return [meta, cluster]
                cluster: cluster.countFasta() >= params.umi_bin_size && cluster.countFasta() > 1
                    return [meta, cluster]
        }
        .set{ ch_umi_bysize }

        // Get the correspondent fasta sequencences from single clusters
        ch_umi_bysize.single
        .tap{ meta_channel }
        .map{ meta, cluster ->
            fasta_line = umi_to_sequence(cluster)
            [meta, cluster.baseName, fasta_line]
        }
        .collectFile() { meta, name, fasta ->
            [ "${name}_consensus.fasta", fasta ]
        }
        .map{ new_file ->
            [new_file.baseName, new_file]
        }
        .join(meta_channel
            .map { meta, original_file ->
                ["${original_file.baseName}_consensus", meta]
            })
        .map{ file_name, new_file, meta ->
            [meta, new_file]
        }
        .set{ ch_single_clusters_consensus }

        //
        // MODULE: Obtain most abundant UMI in cluster
        //
        VSEARCH_SORT(
            ch_umi_bysize.cluster,
            Channel.value("--sortbysize")
        )

        // Get the correspondent fasta sequencences from top cluster sequences
        // Replaces the sequence name adding the "centroid_" prefix to avoid having two sequences with the same name in following steps
        VSEARCH_SORT.out.fasta // [[id:sample_id, ...], sample_top.fasta]
        .tap{ meta_channel_2 } // [[id:sample_id, ...], sample_top.fasta]
        .map{ meta, fasta ->
            fasta_line = umi_to_sequence_centroid(fasta)
            [meta, fasta.baseName, fasta_line] // [[id:sample_id, ...], sample_top, >centroid_...]
        }
        .collectFile(cache:true,sort:true) { meta, name, fasta ->
            [ "${name}.fasta", fasta ] // >centroid_... -> sample_top.fasta
        }
        .map{ new_file ->
            [new_file.baseName, new_file] // [sample, sample_top.fasta]
        }
        .join(meta_channel_2
            .map { meta, original_file ->
                ["${original_file.baseName}", meta] // [sample, [id:sample_id, ...]]
            }) // [sample, sample_top.fasta, [id:sample_id, ...]]
        .map{ file_name, new_file, meta ->
            [meta + [cluster_id: file_name[0..-5]], new_file] // Add cluster ID to meta map // [[id:sample_id, ..., cluster_id:sample], sample_top.fasta]
        }
        .set{ ch_top_clusters_sequence } // [[id:sample_id, ..., cluster_id:sample], sample_top.fasta]


        // Get the correspondent fasta sequencences from UMI clusters
        ch_umi_bysize.cluster // [[id:sample_id, ...], sample]
        .tap{ meta_channel_3 } // [[id:sample_id, ...], sample]
        .map{ meta, cluster ->
            fasta_line = umi_to_sequence(cluster)
            [meta, cluster.baseName, fasta_line] // [[id:sample_id, ...], sample, >...]
        }
        .collectFile(cache:true,sort:true) { meta, name, fasta ->
            [ "${name}_sequences.fasta", fasta ] // >... -> sample_sequences.fasta
        }
        .map{ new_file ->
            [new_file.baseName[0..-11], new_file] // [sample, sample_sequences.fasta]
        }
        .join(meta_channel_3
            .map { meta, original_file ->
                ["${original_file.baseName}", meta] // [sample, [id:sample_id, ...]]
            }) // [sample, sample_sequences.fasta, [id:sample_id, ...]]
        .map{ file_name, new_file, meta ->
            [meta + [cluster_id: file_name], new_file] // Add cluster ID to meta map // [[id:sample_id, ..., cluster_id:sample], sample_sequences.fasta]
        }
        .set{ ch_clusters_sequence }


        // Cluster consensus & polishing
        // Two cycles of minimap2 + racon

        //
        // MODULE: Mapping with minimap2 - cycle 1
        //
        // Map each cluster against the top read (most abundant UMI) in the cluster
        MINIMAP2_ALIGN_UMI_1 (
            ch_clusters_sequence
                .join(ch_top_clusters_sequence),
            false, //output in paf format
            false,
            false
        )


        // Only continue with clusters that have aligned sequences
        MINIMAP2_ALIGN_UMI_1.out.paf
            .filter{ it[1].countLines() > 0 }
            .set{ ch_minimap_1 }

        //
        // MODULE: Improve top read from UMI cluster using cluster consensus - cycle 1
        //
        RACON_1 (
            ch_clusters_sequence
                .join(ch_top_clusters_sequence)
                .join(ch_minimap_1)
        )

        //
        // MODULE: Mapping with minimap2 - cycle 2
        //
        // Map each cluster against the top read (most abundant UMI) in the cluster
        MINIMAP2_ALIGN_UMI_2 (
            ch_clusters_sequence
                .join(RACON_1.out.improved_assembly),
            false, //output in paf format
            false,
            false
        )

        // Only continue with clusters that have aligned sequences
        MINIMAP2_ALIGN_UMI_2.out.paf
            .filter{ it[1].countLines() > 0 }
            .set{ ch_minimap_2 }

        //
        // MODULE: Improve top read from UMI cluster using cluster consensus - cycle 2
        //
        RACON_2 (
            ch_clusters_sequence
                .join(RACON_1.out.improved_assembly)
                .join(ch_minimap_2)
        )


        //
        // MODULE: Obtain a consensus sequence
        //
        MEDAKA (
            ch_clusters_sequence
                .join(RACON_2.out.improved_assembly)
        )

        // Collect all consensus UMI sequences into one single file per sample
        MEDAKA.out.assembly
        .tap{ meta_channel_4 }
        .map{ meta, file ->
            file_content = file.getText()
            [meta, file_content] // [[id:sample_id, ...], consensus_content]
        }
        .collectFile() { meta, file ->
            [ "${meta.id}_consensus.fasta", file ]
        }
        .map{ new_file ->
            [new_file.baseName, new_file]
        }
        .join(meta_channel_4
            .map{ meta, consensus ->
                ["${meta.id}_consensus", meta]
            }
        )
        .map{ name, file, meta ->
            [meta - meta.subMap('cluster_id'), file]
        }
        .set{ ch_umi_consensus }


        //
        // MODULE: Convert fasta to fastq
        //
        SEQTK_SEQ_FATOFQ (
            ch_umi_consensus
        )
    }

    ch_preprocess_reads = params.umi_clustering ? SEQTK_SEQ_FATOFQ.out.fastx : SEQTK_SEQ_MASK.out.fastx

    //
    // MODULE: Summary of clustered reads
    //
    CLUSTERING_SUMMARY (
        ch_preprocess_reads
            .join(MERGING_SUMMARY.out.summary)
    )


    //
    // MODULE: Mapping with Minimap2
    //
    if (params.aligner == "minimap2") {
        MINIMAP2_ALIGN_ORIGINAL (
            ch_preprocess_reads
                .join(ORIENT_REFERENCE.out.reference),
            true,
            false,
            true
        )
        ch_mapped_bam = MINIMAP2_ALIGN_ORIGINAL.out.bam
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_ORIGINAL.out.versions)
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
            ch_preprocess_reads,
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
            ch_preprocess_reads,
            BOWTIE2_BUILD.out.index,
            false,
            true
        )
        ch_mapped_bam = BOWTIE2_ALIGN.out.aligned
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
            .join(ch_seq_template)
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
        .join(ch_input.protospacer
            .map {
                meta, fastq ->
                    [ meta , fastq ]
            }
        )
        .join(ch_seq_template, remainder: true)
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

    methods_description    = WorkflowCrisprseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
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
