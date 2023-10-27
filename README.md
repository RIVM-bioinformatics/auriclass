![install with bioconda](https://img.shields.io/conda/v/bioconda/auriclass)

![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/rivm-bioinformatics/auriclass/test.yaml?label=Tests)
![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/rivm-bioinformatics/auriclass/super_linter.yaml?label=Linting)
![GitHub deployments](https://img.shields.io/github/deployments/RIVM-bioinformatics/Auriclass/github-pages?label=Documentation%20deployment)  


[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://github.com/RIVM-bioinformatics/auriclass/blob/master/LICENSE)
![GitHub release](https://img.shields.io/github/v/release/rivm-bioinformatics/auriclass)

# AuriClass: quick estimation of *Candida auris* clade membership

AuriClass is a small tool which predicts *Candida auris* clade based on Mash distances from reference genomes. It accepts fastq or fasta files. Analysis typically takes a minute for Fastq data and a couple of seconds for Fasta data.

Documentation is available from https://rivm-bioinformatics.github.io/auriclass.

## Installation
The easiest way to install AuriClass is through mamba/conda:

```
mamba create -n env_auriclass -c bioconda -c conda-forge auriclass
conda activate env_auriclass
```

Alternatively, [biocontainers](https://quay.io/repository/biocontainers/auriclass) has built a container from the AuriClass bioconda package. This container can be used with different container software, e.g. Singularity:

```
singularity pull docker://quay.io/biocontainers/auriclass:0.5.1--pyhdfd78af_0
```

Please note that `auriclass:latest` is not defined, so make sure to specify the latest bioconda version and build.

## Minimal examples

Running AuriClass on only the forward reads gives the best results, as the reverse reads are usually more noisy:

```
auriclass Candida_auris_R1.fq.gz
```

But it can also be run on an arbitrary number of fastq files of the same organism:

```
auriclass Candida_auris_R1.fq.gz Candida_auris_R2.fq.gz Candida_auris_unpaired.fq.gz
```

or with a fasta file:

```
auriclass Candida_auris.fasta.gz
```

A standard analysis creates two files:
- `report.tsv` which contains the clade prediction, closest reference sample and QC checks. The report contains the default "isolate" as sample name.
- `report.YYYY-mm-dd_HH-MM-SS.log` which contains messages written to STDERR with additional data.
