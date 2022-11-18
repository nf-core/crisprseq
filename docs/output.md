# nf-core/crisprseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [Sequences](#sequences) - Input sequence preparation (reference, protospacer, template)
  - [cat](#cat) - Concatenate sample fastq files if requiered
  - [Pear](#pear) - Join double-end reads if required
  - [FastQC](#fastqc) - Read Quality Control
  - [Adapters](#adapters) - Find adapters (Overrepresented sequences) in reads
  - [Cutadapt](#cutadapt) - Trim adapters
  - [Seqtk](#seqtk) - Mask low-quality bases
  <!-- -UMI(#umi) -->
- [Mapping](#mapping)
  - [minimap2](#minimap2) - Mapping reads to reference
  - [BWA](#bwa) - Mapping reads to reference
  - [bowtie2](#bowtie2) - Mapping reads to reference
- [Edits calling](#edits-calling)
  - [CIGAR](#cigar) - Parse CIGAR to call edits
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing

### Sequences

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/sequences/`
  - `*_reference.fasta`: Sequence used as a reference.
  - `*_template.fasta`: Provided template sequence.
  - `*_correctOrient.fasta`: Reference sequence in the correct orientation.
  - `NewRef.fasta`: New reference generated from adding the changes made by the template to the original reference.
  - `*_template-align.bam`: Alignment of the new reference (with template changes) to the original reference.

</details>

Contains the input sequences (reference, protospacer and template). Sequences are preprocessed as required:
- The reference is returned in the correct orientation.
  > In order to provide the reference in the correct orientation, the protospacer is searched in the reference sequence. The reverse complement is returned if the protospacer matches the reference in reverse complement.
- The template is used to obtain a new reference with the expected changed.

### cat

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/cat/`
  - `*.merged.fastq.gz`: Concatenated fastq files

</details>

If multiple libraries/runs have been provided for the same sample in the input samplesheet (e.g. to increase sequencing depth) then these will be merged at the very beginning of the pipeline in order to have consistent sample naming throughout the pipeline. Please refer to the [usage](https://nf-co.re/crisprseq/usage) documentation to see how to specify these samples in the input samplesheet.

### Pear

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/pear/`
  - `*.assembled.fastq.gz`: Assembled paired-end reads
  - `*.discarded.fastq.gz`: Discarded reads
  - `*.unassembled.forward.fastq.gz`: Unassembled paired-end reads - forward (R1)
  - `*.unassembled.reverse.fastq.gz`: Unassembled paired-end reads - reverse (R2)

</details>

[PEAR](https://cme.h-its.org/exelixis/web/software/pear/) is a pair-end read merger.

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

- `preprocessing/adapters/`
  - `*_overrepresented.fasta`: Contains overrepresented sequences found by FastQC

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) finds overrepresented sequences found in samples. It lists all of the sequence which make up more than 0.1% of the total. For each overrepresented sequence the program will look for matches in a database of common contaminants and will report the best hit it finds. Hits must be at least 20bp in length and have no more than 1 mismatch.

### Cutadapt

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/cutadapt/`
  - `*.cutadapt.log`: Cutadapt log file
  - `*.trim.fastq.gz`: Sample reads trimmed with overrepresented sequences removed

</details>

### Seqtk

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/seqtk/`
  - `*.seqtk-seq.fastq.gz`: Quality filtered reads.

</details>

[Seqtk](https://github.com/lh3/seqtk) masks (converts to Ns) bases with quality lower than 20 and removes sequences shorter than 80 bases.

<!-- ### UMI -->
## Mapping
### minimap2

<details markdown="1">
<summary>Output files</summary>

- `minimap2/`
  - `*xx.html`: ...

</details>

xxx

### BWA

<details markdown="1">
<summary>Output files</summary>

- `bwa/`
  - `*xx.html`: ...

</details>

xxx

### bowtie2

<details markdown="1">
<summary>Output files</summary>

- `bowtie2/`
  - `*xx.html`: ...

</details>

xxx

## Edits calling
### CIGAR

<details markdown="1">
<summary>Output files</summary>

- `cigar/`
  - `*xx.html`: ...

</details>

xxx

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
