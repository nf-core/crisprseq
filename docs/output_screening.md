# nf-core/crisprseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.


## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [FastQC](#sequences) - Input sequence preparation (reference, protospacer, template)
  - [FastQC](#fastqc) - Read Quality Control
- [Counting](#counting)
  - [MAGeCK count](#count) - Mapping reads to reference
- [CNV correction](#counting)
  - [CRISPRcleanR](#crisprcleanr) - Copy Number Variation correction and read normalization in case of knock-out screens.
- [Gene essentiality](#gene-essentiality)
  - [MAGeCK rra](#rra) - modified robust ranking aggregation (RRA) algorithm
  - [MAGeCK mle](#mle) -  maximum-likelihood estimation (MLE) for robust identification of CRISPR-screen hits
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

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### Adapters

<details markdown="1">
<summary>Output files</summary>



</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) finds over-represented sequences in samples. It lists all of the sequence which make up more than 0.1% of the total reads. For each over-represented sequence the program will look for matches in a database of common contaminants and will report the best hit it finds. Hits must be at least 20bp in length and have no more than 1 mismatch.

### Cutadapt

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/cutadapt/`
  - `*.cutadapt.log`: Cutadapt log file
  - `*.trim.fastq.gz`: Sample reads trimmed with overrepresented sequences removed

</details>

### Mageck MLE

<details markdown="1">
<summary>Output files</summary>

- `mageck`
  - `*.seqtk-seq.fastq.gz`: Quality filtered reads.

</details>

[Seqtk](https://github.com/lh3/seqtk) masks (converts to Ns) bases with quality lower than 20 and removes sequences shorter than 80 bases.

<!-- ### UMI -->

## Mapping

### minimap2

<details markdown="1">
<summary>Output files</summary>

- `minimap2/`
  - `*.bam`: BAM file containing aligned reads
  - `*.bai`: BAI index

</details>

[Minimap2](https://github.com/lh3/minimap2) is a sequence alignment program that aligns DNA sequences against a reference database.

### BWA

<details markdown="1">
<summary>Output files</summary>

- `bwa/`
  - `*.bam`: BAM file containing aligned reads
  - `*.bai`: BAI index

</details>

[BWA-MEM](https://github.com/lh3/bwa) BWA is a software package for mapping low-divergent sequences against a reference genome.

### bowtie2

<details markdown="1">
<summary>Output files</summary>

- `bowtie2/`
  - `*.bam`: BAM file containing aligned reads
  - `*.bai`: BAI index

</details>

[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) aligns sequencing reads to reference sequences.

## Edits calling

### CIGAR

<details markdown="1">
<summary>Output files</summary>

- `cigar/`
  - `*_cutSite.json`: Contains the protospacer cut site position in the reference.
  - `*_edition.html`: Interactive pie chart with the percentage of edition types. Reads are classified between WT (without an edit) and indels. Indes are divided between deletions, insertions and delins (deletion + insertion). Deletions and insertions can be out of frame or in frame.
    ![Test sample hCas9-AAVS1-a edition plot](images/hCas9-AAVS1-a_edition.png)
  - `*_edits.csv`: Table containing the data visualized in the pie chart.
  - `*_indels.csv`: Table containing information of all reads. Edit type, edit start and length, if the edition happens above the error rate, if it's located into the common edit window, the frequency, the percentage, the pattern, surrounding nucleotides in case of insertions, the protospacer cut site, the sample id, number of aligned reads and number of reads with and without a template modification.
  - `*_QC-indels.html`: Interactive pie chart with information about aligned reads. Reads are classified between WT and containing indels. Both types are classified between passing the filtering steps or not. Indel reads passing the filtering steps are divided in reads with a modification above the error rate and located in the common edit window, above the error rate but not in the edit region, viceversa, or any of those conditions.
    ![Test sample hCas9-AAVS1-a QC indels plot](images/hCas9-AAVS1-a_QC-indels.png)
  - `*_reads.html`: Interactive pie chart with percentage of the number of raw reads, reads merged with Pear, reads passing quality filters and UMI clustered reads.
    ![Test sample hCas9-AAVS1-a reads plot](images/hCas9-AAVS1-a_reads.png)
  - `*_subs-perc.csv`: Table containing the percentage of each nucleotide found for each reference position.

</details>

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
