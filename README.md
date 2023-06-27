# ![nf-core/crisprseq](docs/images/nf-core-crisprseq_logo_light.png#gh-light-mode-only) ![nf-core/crisprseq](docs/images/nf-core-crisprseq_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/crisprseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.7598497-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.7598497)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/crisprseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23crisprseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/crisprseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/crisprseq** is a bioinformatics best-practice analysis pipeline for the analysis of CRISPR edited data. It allows the evaluation of the quality of gene editing experiments using targeted next generation sequencing (NGS) data (`targeted`) as well as the discovery of important genes from knock-out or activation CRISPR-Cas9 screens using CRISPR pooled DNA (`screening`).

nf-core/crisprseq can be used to analyse:

- CRISPR gene knockouts (KO)
- CRISPR knock-ins (KI)
- Base editing (BE) and prime editing (PE) experiments
- CRISPR screening experiments (KO, CRISPRa (activation) or CRISPRi (interference))

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/crisprseq/results).

## Pipeline summary

For crispr targeted:

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/nf-core/crisprseq/dev/docs/images/crisprseq_targeted_metro_map_dark.png">
  <img alt="Text changing depending on mode. Light: 'So light!' Dark: 'So dark!'" src="https://raw.githubusercontent.com/nf-core/crisprseq/dev/docs/images/crisprseq_targeted_metro_map.png">
</picture>

1. Merge paired-end reads ([`Pear`](https://cme.h-its.org/exelixis/web/software/pear/doc.html))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter trimming ([`Cutadapt`](http://dx.doi.org/10.14806/ej.17.1.200))
4. Quality filtering ([`Seqtk`](https://github.com/lh3/seqtk))
5. UMI clustering (optional):
   a. Extract UMI sequences (Python script)
   b. Cluster UMI sequences ([`Vsearch`](https://github.com/torognes/vsearch))
   c. Obtain the most abundant UMI sequence for each cluster ([`Vsearch`](https://github.com/torognes/vsearch))
   d. Obtain a consensus for each cluster ([`minimap2`](https://github.com/lh3/minimap2))
   e. Polish consensus sequence ([`racon`](https://github.com/lbcb-sci/racon))
   f. Repeat a second rand of consensus + polishing (`minimap2` + `racon`)
   g. Obtain the final consensus of each cluster ([Medaka](https://nanoporetech.github.io/medaka/index.html))
6. Read mapping:
   - ([`minimap2`](https://github.com/lh3/minimap2), _default_)
   - ([`bwa`](http://bio-bwa.sourceforge.net/))
   - ([`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
7. CIGAR parsing for edit calling ([`R`](https://www.r-project.org/))

For crispr screening:

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Read mapping ([`MAGeCK count`](https://sourceforge.net/p/mageck/wiki/usage/#count))
3. Optional: CNV correction and normalization with ([`CRISPRcleanR`](https://github.com/francescojm/CRISPRcleanR))
4. Rank sgRNAs and genes ;
   a. ([MAGeCK test](https://sourceforge.net/p/mageck/wiki/usage/#test))
   b. ([MAGeCK mle](https://sourceforge.net/p/mageck/wiki/Home/#mle))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run nf-core/crisprseq -profile test_screening,YOURPROFILE --outdir <OUTDIR>
   # or
   nextflow run nf-core/crisprseq -profile test_targeted,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run nf-core/crisprseq --input samplesheet.csv --analysis <targeted/screening> --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

The nf-core/crisprseq pipeline comes with documentation about the pipeline [usage](https://nf-co.re/crisprseq/usage), [parameters](https://nf-co.re/crisprseq/parameters) and [output](https://nf-co.re/crisprseq/output).

## Credits

nf-core/crisprseq targeted is based on [CRISPR-A](https://doi.org/10.1101/2022.09.02.506351) [[Sanvicente-García, et.al. (2023)](https://doi.org/10.1371/journal.pcbi.1011137)], originally written by Marta Sanvicente García at [Translational Synthetic Biology](https://synbio.upf.edu/) from [Universitat Pompeu Fabra](https://www.upf.edu/home).

The pipeline was re-written in Nextflow DSL2 and is primarily maintained by Júlia Mir Pedrol (@mirpedrol) at [Quantitative Biology Center (QBiC)](https://www.qbic.uni-tuebingen.de/) from [Universität Tübingen](https://uni-tuebingen.de/en/).

nf-core/crisprseq screening was written and is primarly maintained by Laurence Kuhlburger (@LaurenceKuhl) at [Quantitative Biology Center (QBiC)](https://www.qbic.uni-tuebingen.de/) from [Universität Tübingen](https://uni-tuebingen.de/en/).

<!--We thank the following people for their extensive assistance in the development of this pipeline:

 TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#crisprseq` channel](https://nfcore.slack.com/channels/crisprseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/crisprseq for your analysis, please cite it using the following doi: [10.5281/zenodo.7598497](https://doi.org/10.5281/zenodo.7598497)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

Crispr-Analytics:

> Sanvicente-García M, García-Valiente A, Jouide S, Jaraba-Wallace J, Bautista E, Escobosa M, et al. (2023)
> CRISPR-Analytics (CRISPR-A): A platform for precise analytics and simulations for gene editing. PLoS Comput Biol 19(5): e1011137. https://doi.org/10.1371/journal.pcbi.1011137
