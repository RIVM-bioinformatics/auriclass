![install with bioconda](https://img.shields.io/conda/v/bioconda/auriclass)
![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/rivm-bioinformatics/auriclass/test.yaml?label=Tests)
![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/rivm-bioinformatics/auriclass/super_linter.yaml?label=Linting)
![GitHub deployments](https://img.shields.io/github/deployments/RIVM-bioinformatics/Auriclass/github-pages?label=Documentation%20deployment)  
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://github.com/RIVM-bioinformatics/auriclass/blob/master/LICENSE)
![GitHub release](https://img.shields.io/github/v/release/rivm-bioinformatics/auriclass)

# AuriClass: quick estimation of *Candida auris* clade membership

AuriClass is a small tool which predicts *Candida auris* clade based on Mash distances from reference genomes. It accepts fastq or fasta files. Analysis typically takes a minute for Fastq data and a couple of seconds for Fasta data.

## Quickstart

Install the tool using mamba and run an analysis on only the forward reads:

```
mamba create -n env_auriclass -c bioconda -c conda-forge auriclass
conda activate env_auriclass
auriclass -o output_report.tsv Candida_auris_R1.fq.gz
```

## Motivation

We needed a tool that:

- quickly and accurately predicts *C. auris* clade from Illumina fastq and from fasta data
- could be used on a single sample in an automated workflow

If you have fasta data, also check out the [cauris_cladetyper task](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/task_cauris_cladetyper.wdl) from TheiaEuk workflow.