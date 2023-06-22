---
order: 2
---

# nf-core/crisprseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/crisprseq/usage](https://nf-co.re/crisprseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The **nf-core/crisprseq** pipeline allows the analysis of CRISPR edited CRISPR pooled DNA. It can evaluate important genes from knock-out or activation CRISPR-Cas9 screens.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/crisprseq --analysis screening --input samplesheet.csv --library library.csv --outdir <OUTDIR> -profile docker
```

The following required parameters are here described.

### Full samplesheet

The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

```console
sample,fastq_1,fastq_2,condition
SRR8983579,SRR8983579.small.fastq.gz,,control
SRR8983580,SRR8983580.small.fastq.gz,,treatment
```

| Column      | Description                                                                                                                |
| ----------- | -------------------------------------------------------------------------------------------------------------------------- |
| `sample`    | Custom sample name. . Spaces in sample names are automatically converted to underscores (`_`).                             |
| `fastq_1`   | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `fastq_2`   | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `condition` | Condition of the sample, for instance "treatment" or "control".                                                            |

An [example samplesheet](https://github.com/nf-core/test-datasets/blob/crisprseq/testdata/samplesheet_test.csv) has been provided with the pipeline.

The pipeline currently supports 2 algorithms to detect gene essentiality, MAGeCK rra and MAGeCK mle. MAGeCK MLE (Maximum Likelihood Estimation) and MAGeCK RRA (Robust Ranking Aggregation) are two different methods provided by the MAGeCK software package to analyze CRISPR-Cas9 screens.

### MAGeCK rra

MAGeCK RRA performs robust ranking aggregation to identify genes that are consistently ranked highly across multiple replicate screens. To run MAGeCK rra, `--rra_contrasts` should be used with a `csv` separated file stating the two conditions to be compared.

### MAGeCK mle

MAGeCK MLE uses a maximum likelihood estimation approach to estimate the effects of gene knockout on cell fitness. It models the read count data of guide RNAs targeting each gene and estimates the dropout probability for each gene. MAGeCK mle requires a design matrix. The design matrix is a `txt` file indicating the effects of different conditions on different samples.
An [example design matrix](https://github.com/nf-core/test-datasets/blob/crisprseq/testdata/design_matrix.txt) has been provided with the pipeline.
If there are several designs to be run, you can input a folder containing all the design matrices. The output results will automatically take the name of the design matrix, so make sure you give a meaningful name to the file, for instance "Drug_vs_control.txt".

### Running CRISPRcleanR

CRISPRcleanR is used for gene count normalization and the removal of biases for genomic segments for which copy numbers are amplified. Currently, the pipeline only supports annotation libraries already present in the R package and which can be found [here](https://github.com/francescojm/CRISPRcleanR/blob/master/Reference_Manual.pdf). To use CRISPRcleanR normalization, use `--crisprcleanr library`, `library` being the exact name as the library in the CRISPRcleanR documentation (e.g: "AVANA_Library").

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other hidden nextflow files, eg. history of pipeline runs and old logs.
```
