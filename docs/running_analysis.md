# Running an analysis

Running AuriClass on only the forward reads gives the best results, as the reverse reads are usually more noisy:

```
auriclass Candida_auris_R1.fq.gz
```

But it can also be run on an arbitrary number of fastq files of the same sample:

```
auriclass Candida_auris_R1.fq.gz Candida_auris_R2.fq.gz Candida_auris_unpaired.fq.gz
```

or with a fasta file:

```
auriclass Candida_auris.fasta.gz
```

!!! note
    All input files are treated as a single sample, so running `auriclass input_genomes/*.fasta` will try to do a single clade prediction for all fasta files in `input_genomes`. To run AuriClass on multiple samples, check out [Running multiple samples in parallel](#running-multiple-samples-in-parallel)

A standard analysis creates two files:
- `report.tsv` which contains the clade prediction, closest reference sample and QC checks. The report contains the default "isolate" as sample name.
- `report.YYYY-mm-dd_HH-MM-SS.log` which contains messages written to STDERR with additional data.

### Specifying output

This will create `report_Candida_auris.tsv` and `report_Candida_auris.log`. The report will contain the sample name "Candida_auris".

```
auriclass --name Candida_auris -o report_Candida_auris.tsv Candida_auris_R1.fq.gz
```

### Forcing full analysis

This forces AuriClass to perform the full analysis on a sample, even if the species quality check fails.

```
auriclass --non_candida_threshold 1 Candida_auris_R1.fq.gz
```

!!! warning
    This will cause some QC checks to always pass, so use with caution.

### Running multiple samples in parallel

AuriClass currently does not support multithreaded analysis. For both Fastq and Fasta input, multi-threading `mash` commands where possible makes minimal difference. Additionally, peak memory usage is typically 100-200 Mb RAM according to `/usr/bin/time -v`. Therefore, the best option if you have to analyse a lot of files is to run multiple AuriClass analyses concurrently. 

Two good options are either:

- workflow managers such as [Snakemake](https://snakemake.readthedocs.io/en/stable/), [Nextflow](https://www.nextflow.io/), etc. or
- [GNU `parallel`](https://www.gnu.org/software/parallel/).

An example how to run AuriClass using `parallel` on a set of assemblies in the directory `fasta_input`:

```
parallel auriclass --fasta -n {/.} -o report.{/.}.tsv {} ::: fasta_input/*
```