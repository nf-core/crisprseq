---
order: 1
---

# nf-core/crisprseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/crisprseq/usage](https://nf-co.re/crisprseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The **nf-core/crisprseq** pipeline allows the analysis of CRISPR edited DNA. It evaluates the quality of gene editing experiments using targeted next generation sequencing (NGS) data.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes _(see section below for an explanation of samplesheet columns)_:

```console
sample,fastq_1,fastq_2,reference,protospacer,template
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,GCT...CCT,GGGGCCACTAGGGACAGGAT,
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,GCT...CCT,GGGGCCACTAGGGACAGGAT,
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,GCT...CCT,GGGGCCACTAGGGACAGGAT,
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 6 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 3 samples, where `chr6` is single-end and has a template sequence _(this is a reduced samplesheet, please refer to the [pipeline example saplesheet](https://nf-co.re/crisprseq/1.0/assets/samplesheet.csv) to see the full version)_.

```console
sample,fastq_1,fastq_2,reference,protospacer,template
hCas9-TRAC-a,hCas9-TRAC-a_R1.fastq.gz,hCas9-TRAC-a_R2.fastq.gz,GCT...CCT,GGGGCCACTAGGGACAGGAT,
hCas9-AAVS1-a,hCas9-AAVS1-a_R1.fastq.gz,hCas9-AAVS1-a_R2.fastq.gz,GCT...CCT,GGGGCCACTAGGGACAGGAT,
chr6,chr6-61942198-61942498_R1.fastq.gz,,CAA...GGA,TTTTATGATATTTATCTTTT,TTC...CAA
```

| Column        | Description                                                                                                                                                                            |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`      | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`     | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`     | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". (Optional)                                                  |
| `reference`   | Reference sequence of the target region.                                                                                                                                               |
| `protospacer` | Sequence of the protospacer used for CRISPR editing. Must not includ the PAM.                                                                                                          |
| `template`    | Sequence of the template used in templet-based editing experiments. (Optional)                                                                                                         |

An [example samplesheet](https://nf-co.re/crisprseq/1.0/assets/samplesheet.csv) has been provided with the pipeline.

## Optional pipeline steps

### Trimming of overrepresented sequences

To trim the overrepresented sequences found with FastQC from the reads, use the parameter `--overrepresented`.
Such sequences are not trimmed by default.
When using the `--overrepresented` parameter, Cutadapt is used to trim overrepresented sequences from the input FASTQ files.

### UMI clustering

If the provided samples were sequenced using umi-molecular identifiers (UMIs), use the parameter `--umi_clustering` in order to run the clustering steps.

1. Extract UMI sequences (Python script)
2. Cluster UMI sequences ([`Vsearch`](https://github.com/torognes/vsearch))
3. Obtain the most abundant UMI sequence for each cluster ([`Vsearch`](https://github.com/torognes/vsearch))
4. Obtain a consensus for each cluster ([`minimap2`](https://github.com/lh3/minimap2))
5. Polish consensus sequence ([`racon`](https://github.com/lbcb-sci/racon))
6. Repeat a second round of consensus + polishing (`minimap2` + `racon`)
7. Obtain the final consensus of each cluster ([Medaka](https://nanoporetech.github.io/medaka/index.html))

## Other input parameters

### Reference

If you want to provide the same reference for every sample, you can select a genome with `--genome` or provide a reference FASTA file with `--reference_fasta`.
Using any of these two parameters will override any reference sequence provided through an input sample sheet.

Please refer to the [nf-core website](https://nf-co.re/usage/reference_genomes) for general usage docs and guidelines regarding reference genomes.

### Protospacer

If you want to provide the same protospacer sequence for every sample, you can provide the sequence with the parameter `--protospacer`.
Using this parameter will override any protospacer sequence provided through an input sample sheet.

Providing a protospacer, either through a sample sheet or by using the parameter `--protospacer` is required.

## Alignment options

By default, the pipeline uses `minimap2` (i.e. `--aligner minimap2`) to map the sequenced FASTQ reads to the reference.
You also have the option to select other alignment tools by using the parameter `--alignment`. Possible options are `minimap2`, `bwa` or `bowtie2`.

The default alignment with `minimap2` uses adapted parameters which were seen to improve the alignment and reduce potential sequencing or alignment errors.
The default parameters are:

- A matching score of 29
- A mismatching penalty of 17
- A gap open penalty of 25
- A gap extension penalty of 2.

Please refer to the original [CRISPR-Analytics](https://doi.org/10.1371/journal.pcbi.1011137) publication to see the benchmarking of such parameters.

In order to customise such parameters, you can override the arguments given to `minimap2` by creating a configuration file and provide it to your nextflow run with `-c`:

```groovy
// Custom config file custom.config
process {
    withName: MINIMAP2_ALIGN_ORIGINAL {
        ext.args = '-A 29 -B 17 -O 25 -E 2'
    }
}
```

Command:

```bash
nextflow run nf-core/crisprseq --input samplesheet.csv --analysis targeted --outdir <OUTDIR> -profile docker -c custom.config
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/crisprseq --input samplesheet.csv --analysis targeted --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```
