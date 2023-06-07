# nf-core/crisprseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/crisprseq/usage](https://nf-co.re/crisprseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The **nf-core/crisprseq** pipeline allows the analysis of CRISPR edited DNA. It evaluates the quality of gene editing experiments using targeted next generation sequencing (NGS) data (`targeted`) as well as important genes from knock-out or activation CRISPR-Cas9 screens using CRISPR pooled DNA (`screening`).

## Type of analysis

The `--analysis` parameter specifies whether the user intends to perform editing or screening in the crisprseq pipeline.

Please refer to the respective Usage documentation:
- [Usage targeted](usage-targeted)
- [Usage screening](usage-screening)
