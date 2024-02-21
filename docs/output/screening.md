---
order: 2
---

# nf-core/crisprseq: Output

## Introduction

This document describes the output produced by the pooled screens analysis of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [FastQC](#fastqc) - Read Quality Control
  - [cutadapt](#cutadapt) - Trimming reads from fastq files
- [Mapping](#mapping) - bowtie2 aligned reads
- [Counting](#counting)
  - [MAGeCK count](#mageck-count) - Mapping reads to reference library
- [CNV correction](#cnv-correction))
  - [CRISPRcleanR](#crisprcleanr-normalization) - Copy Number Variation correction and read normalization in case of knock-out screens.
- [Gene essentiality](#gene-essentiality-computation)
  - [MAGeCK rra](#mageck-rra) - modified robust ranking aggregation (RRA) algorithm
  - [MAGeCK mle](#mageck-mle) - maximum-likelihood estimation (MLE) for robust identification of CRISPR-screen hits
  - [BAGEL2](#BAGEL2) - Bayes Factor to identify essential genes
  - [MAGeCKFlute](#flutemle) - graphics to visualise MAGECK MLE output
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### cutadapt

<details markdown="1">
<summary>Output files</summary>

- `cutadapt/`
  - `*.log`: log file of the command ran and the output
  - `*.trim.fastq.gz`: trimmed fastq files

</details>

[cutadapt](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads. MAGeCK count normally automatically detects adapter sequences and trims, however if trimming lengths are different, cutadapt can be used, as mentioned [here](https://sourceforge.net/p/mageck/wiki/advanced_tutorial/).
For further reading and documentation see the [cutadapt helper page](https://cutadapt.readthedocs.io/en/stable/guide.html).

## Alignment

<details markdown="1">
<summary>Output files</summary>

- `bowtie2/`
  - `*.log`: log file of the command ran and the output
  - `*.bam`: bam file
  - `*.bowtie2`: index from bowtie2 from the provided fasta file

</details>

## Counting

### MAGeCK count

<details markdown="1">
<summary>Output files</summary>

- `mageck/count`
  - `*_count.txt`: read counts per sample per sgRNA and gene, tab separated
  - `*_count_normalized.txt`: normalized read counts, tab separated
  - `*_count_summary.txt`: tab separated summary of the quality controls of the count table
  - `*_count_table.log`: log information of the run

</details>

## CNV correction

### CRISPRcleanR normalization

<details markdown="1">
<summary>Output files</summary>

- `CRISPRcleanR/normalization`
  - `*_norm_table.tsv`: read counts normalized with crisprcleanr
  - `*.RData`: RData tables containing corrected counts, fold changes and normalized counts
  </details>

## Gene essentiality computation

### MAGeCK mle

<details markdown="1">
<summary>Output files</summary>

- `mageck/mle`
  - `*_gene_summary.txt`: ranked table of the genes and their associated p-values
  - `*_sgrna_summary.txt`: sgRNA ranking results, tab separated file
  - `*.log`: log of the run

</details>

### MAGeCK rra

<details markdown="1">
<summary>Output files</summary>

- `mageck/rra`
  - `*_gene_summary.txt`: ranked table of the genes and their associated p-values
  - `*_count_sgrna_summary.txt`: sgRNA ranking results, tab separated file containing means, p-values
  - `*.report.Rmd`: markdown report recapping essential genes
  - `*_count_table.log`: log of the run
  - `*_scatterview.png`: scatter view of the targeted genes in the library and their logFC
  - `*_rank.png`: rank view of the targeted genes in the library

</details>

[MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) is a computational tool to identify important genes from CRISPR-Cas9 screens.

### BAGEL2

<details markdown="1">
<summary>Output files</summary>

- `bagel2/fold_change`
  - `*.foldchange`: foldchange between the reference and treatment contrast provided
- `bagel2/bayes_factor`
  - `*.bf`: bayes factor per gene
- `bagel2/precision_recall`
  - `*.pr`: precision recall per gene
- `bagel2/graphs`
  - `barplot*.png`: barplot of the bayes factor distribution
  - `PR*.png`: precision recall plot (Recall vs FDR)

</details>

[bagel2](https://github.com/hart-lab/bagel) is a computational tool to identify important essential genes for CRISPR-Cas9 screening experiments.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
