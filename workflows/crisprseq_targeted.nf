/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Local modules
include { FIND_ADAPTERS                             } from '../modules/local/find_adapters'
include { EXTRACT_UMIS                              } from '../modules/local/extract_umis'
include { ORIENT_REFERENCE                          } from '../modules/local/orient_reference'
include { CIGAR_PARSER                              } from '../modules/local/cigar_parser'
include { PREPROCESSING_SUMMARY                     } from '../modules/local/preprocessing_summary'
include { CLUSTERING_SUMMARY                        } from '../modules/local/clustering_summary'
include { ALIGNMENT_SUMMARY                         } from '../modules/local/alignment_summary'
include { TEMPLATE_REFERENCE                        } from '../modules/local/template_reference'
include { CRISPRSEQ_PLOTTER                         } from '../modules/local/crisprseq_plotter'
include { CLONALITY_CLASSIFIER                      } from '../modules/local/clonality_classifier'
// nf-core modules
include { FASTQC                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
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
include { MINIMAP2_INDEX                            } from '../modules/nf-core/minimap2/index/main'
include { MEDAKA                                    } from '../modules/nf-core/medaka/main'
include { CUTADAPT                                  } from '../modules/nf-core/cutadapt/main'
include { SAMTOOLS_INDEX                            } from '../modules/nf-core/samtools/index/main'
// Local subworkflows
include { INITIALISATION_CHANNEL_CREATION_TARGETED  } from '../subworkflows/local/utils_nfcore_crisprseq_pipeline'
// Functions
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_crisprseq_pipeline'

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

workflow CRISPRSEQ_TARGETED {

    take:
    ch_input_reads // channel: input reads read in from --input
    ch_input_reference // channel: reference sequence read in from --input
    ch_input_template // channel: template sequence read in from --input
    ch_input_protospacer // channel: protospacer sequence read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Initialise channels
    //
    INITIALISATION_CHANNEL_CREATION_TARGETED(
        ch_input_reads,
        ch_input_reference,
        ch_input_template,
        ch_input_protospacer
    )

    //
    // MODULE: Prepare reference sequence
    //
    ORIENT_REFERENCE (
        INITIALISATION_CHANNEL_CREATION_TARGETED.out.reference_protospacer
    )
    ch_versions = ch_versions.mix(ORIENT_REFERENCE.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        INITIALISATION_CHANNEL_CREATION_TARGETED.out.fastq_multiple
    )
    .reads
    .groupTuple(by: [0])
    .mix(INITIALISATION_CHANNEL_CREATION_TARGETED.out.fastq_single)
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
    .map {
        // Set single_end to true for the assembled reads
        meta, assembled ->
            return [ meta - meta.subMap('single_end') + [ single_end:true ], assembled ]
    }
    .mix( ch_cat_fastq.single )
    .set { ch_pear_fastq }
    ch_versions = ch_versions.mix(PEAR.out.versions)

    // Change reference, protospacer and template channels to have the same meta information as the reads
    ch_pear_fastq
        .map {meta, reads ->
            // save single_end value and remove the key from the meta map
            single_end = meta.single_end
            return [ meta - meta.subMap('single_end'), reads, single_end ]
        }
        .tap { no_single_end }
        .join(
            ORIENT_REFERENCE.out.reference
            .map {meta, reference ->
                // Remove single_end from the meta map to allow joining two channels with different single_end values
                return [ meta - meta.subMap('single_end'), reference ]
            }
        )
        .map {meta, reads, single_end, reference ->
            // Add the correct single_end value to the reference meta map.
            return [ meta + ["single_end": single_end], reference ]
        }
        .tap{ ch_oriented_reference }
    no_single_end
        .join(
            INITIALISATION_CHANNEL_CREATION_TARGETED.out.template
            .map {meta, template ->
                return [ meta - meta.subMap('single_end'), template ]
            }
        )
        .map {meta, reads, single_end, template ->
            return [ meta + ["single_end": single_end], template ]
        }
        .set{ ch_template }
    no_single_end
        .join(
            INITIALISATION_CHANNEL_CREATION_TARGETED.out.reference_protospacer
            .map {meta, reference, protospacer ->
                return [ meta - meta.subMap('single_end'), protospacer ]
            }
        )
        .map {meta, reads, single_end, protospacer ->
            return [ meta + ["single_end": single_end], protospacer ]
        }
        .set{ ch_protospacer }


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_pear_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
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
        ch_versions = ch_versions.mix(FIND_ADAPTERS.out.versions.first())


        //
        // MODULE: Trim adapter sequences
        //
        CUTADAPT (
            ch_adapter_seqs.adapters
        )
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]})
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
            .map { meta, reads, assembled ->
                // Remove the single_end key from the meta map to allow joining channels with different single_end values
                return [ meta - meta.subMap('single_end'), reads, assembled]
            }
            .join(
                SEQTK_SEQ_MASK.out.fastx
                .map { meta, masked ->
                    single_end = meta.single_end
                    return [ meta - meta.subMap('single_end'), masked, single_end]
                }
            )
            .join(
                CUTADAPT.out.log
                .map { meta, trimmed ->
                    return [ meta - meta.subMap('single_end'), trimmed]
                }
            )
            .map { meta, reads, assembled, masked, single_end, trimmed ->
                if (assembled == null) {
                    assembled = []
                }
                return [ meta + ["single_end": single_end], reads, assembled, masked, trimmed ]
            }
            .set { ch_preprocessing_summary_data }
    } else {
        ch_cat_fastq.paired
            .mix(ch_cat_fastq.single)
            .join(PEAR.out.assembled, remainder: true)
            .map { meta, reads, assembled ->
                // Remove the single_end key from the meta map to allow joining channels with different single_end values
                return [ meta - meta.subMap('single_end'), reads, assembled]
            }
            .join(
                SEQTK_SEQ_MASK.out.fastx
                .map { meta, masked ->
                    single_end = meta.single_end
                    return [ meta - meta.subMap('single_end'), masked, single_end]
                }
            )
            .map { meta, reads, assembled, masked, single_end ->
                if (assembled == null) {
                    assembled = []
                }
                return [ meta + ["single_end": single_end], reads, assembled, masked, [] ]
            }
            .set { ch_preprocessing_summary_data }
    }

    //
    // MODULE: Summary of merged reads
    //
    PREPROCESSING_SUMMARY {
        ch_preprocessing_summary_data
    }
    ch_versions = ch_versions.mix(PREPROCESSING_SUMMARY.out.versions)


    if (params.umi_clustering) {
        //
        // MODULE: Extract UMI sequences
        //
        EXTRACT_UMIS (
            SEQTK_SEQ_MASK.out.fastx
        )
        ch_versions = ch_versions.mix(EXTRACT_UMIS.out.versions.first())


        //
        // MODULE: Cluster UMIs
        //
        VSEARCH_CLUSTER (
            EXTRACT_UMIS.out.fasta
        )
        ch_versions = ch_versions.mix(VSEARCH_CLUSTER.out.versions.first())

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
        ch_versions = ch_versions.mix(VSEARCH_SORT.out.versions.first())

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
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_UMI_1.out.versions.first())


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
        ch_versions = ch_versions.mix(RACON_1.out.versions.first())

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
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_UMI_2.out.versions.first())

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
        ch_versions = ch_versions.mix(RACON_2.out.versions.first())


        //
        // MODULE: Obtain a consensus sequence
        //
        MEDAKA (
            ch_clusters_sequence
                .join(RACON_2.out.improved_assembly)
        )
        ch_versions = ch_versions.mix(MEDAKA.out.versions.first())

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
        ch_versions = ch_versions.mix(SEQTK_SEQ_FATOFQ.out.versions.first())
    }

    ch_preprocess_reads = params.umi_clustering ? SEQTK_SEQ_FATOFQ.out.fastx : SEQTK_SEQ_MASK.out.fastx

    //
    // MODULE: Summary of clustered reads
    //
    CLUSTERING_SUMMARY (
        ch_preprocess_reads
            .join(PREPROCESSING_SUMMARY.out.summary)
    )

    ch_versions = ch_versions.mix(CLUSTERING_SUMMARY.out.versions)



    //
    // MODULE: Mapping with Minimap2
    //
    if (params.aligner == "minimap2") {
        MINIMAP2_ALIGN_ORIGINAL (
            ch_preprocess_reads
                .join(ch_oriented_reference),
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
            ch_oriented_reference
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
            ch_oriented_reference
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
    ch_versions = ch_versions.mix(ALIGNMENT_SUMMARY.out.versions)


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
        ch_oriented_reference
            .join(ch_template)
    )
    ch_versions = ch_versions.mix(TEMPLATE_REFERENCE.out.versions.first())


    //
    // MODULE: Align new reference with the change led by template with original reference
    //
    MINIMAP2_ALIGN_TEMPLATE (
        TEMPLATE_REFERENCE.out.fasta
            .join(ch_oriented_reference),
        true,
        false,
        true
    )
    .bam
    .set { ch_template_bam }
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_TEMPLATE.out.versions)

    ch_mapped_bam
        .join(SAMTOOLS_INDEX.out.bai)
        .join(ch_oriented_reference)
        .join(ch_protospacer)
        .join(ch_template, remainder: true)
        .join(ch_template_bam, remainder: true)
        .join(TEMPLATE_REFERENCE.out.fasta, remainder: true)
        .join(ALIGNMENT_SUMMARY.out.summary)
        .map { meta, reads, index, reference, protospacer, template, template_bam, reference_template, summary ->
            if (template == null) {
                template = []
            }
            if (template_bam == null) {
                template_bam = []
            }
            if (reference_template == null) {
                reference_template = []
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
    ch_multiqc_files = ch_multiqc_files.mix(CIGAR_PARSER.out.processing.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(CIGAR_PARSER.out.edition.collect{it[2]})
    ch_multiqc_files = ch_multiqc_files.mix(CIGAR_PARSER.out.qcindels.collect{it[1]})
    ch_versions = ch_versions.mix(CIGAR_PARSER.out.versions.first())


    //
    // MODULE: Apply clonality classification
    //
    if (!params.skip_clonality) {
        CLONALITY_CLASSIFIER (
            CIGAR_PARSER.out.indels
            .join(CIGAR_PARSER.out.edition)
            .map { [it[0], it[1], it[4]] }
        )
        ch_versions = ch_versions.mix(CLONALITY_CLASSIFIER.out.versions.first())
    }


    //
    //
    //
    CRISPRSEQ_PLOTTER (
        CIGAR_PARSER.out.indels
        .join(ch_oriented_reference)
        .join(ch_protospacer)
    )
    ch_versions = ch_versions.mix(CRISPRSEQ_PLOTTER.out.versions.first())


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
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
