# FAQ

- I get the error `Input file {file_path} is not a fastq or fasta file`. How do I solve this?

This happens when the `pyfastx` library cannot parse the input file(s) as Fastq or Fasta. Please check whether the files are correctly formatted and not corrupted. Plain and gzipped Fastq or Fasta files should work fine, other compression algorithms might not be supported. Check [pyfastx docs](https://pyfastx.readthedocs.io/en/latest/) for more information.

- I get the error `Mash is not installed`. How do I solve this?

AuriClass relies on the availability of a Mash executable. Dependencies like mash and pyfastx can e.g. be installed through mamba/conda.

- I expect that my sample is *C. auris*, but my test sample returns a FAIL message in the output report for the species check. How can I force AuriClass to run the full analysis anyway?

You can force AuriClass to run the whole analysis always by changing the QC thresholds, e.g. setting `--non_candida_threshold 1`. Use with caution, as an important QC step is disabled this way.

- I get a warning stating "WARN: distance {distance} to closest sample is above threshold". What does this mean?

If the lowest distance observed is above the threshold controlled by `--high_distance_threshold` (default 0.003), a warning is emitted. This might mean a couple of different things, listed below. But first, it is important to realise that with assembled fasta input all samples will typically have lower distances to reference genomes than with fastq input, especially if the genomes are closely related. The relative differences remain, although absolute thresholds, like the one set by `--high_distance_threshold` could be lower for fasta input. Additionally, only supplying Illumina forward reads will show lower mash distances compared to supplying both forward and reverse reads.

With that in mind, the warning can mean the following:
1. Your sample is part of a new clade. It seems to be part of *C. auris*, but the distance from any known clade is relatively high. This would usually only be picked up for Illumina forward and reverse read inputs.
2. Your data (probably the reverse reads) are noisy.

In both cases, a comparative method with higher resolution (e.g. reference-based mapping phylogeny) would help solve what's wrong.

- I see "SKIPPED" in the output report for QC columns. What does this mean?

If the sample fails the `qc_species` check (mash distance > `--non_candida_threshold`), the sample is assumed to be a species not related to *Candida*. Therefore, no further analysis is performed and the output report is immediately returned.

"SKIPPED" will also appear for QC_genome_size, QC_multiple_hits and QC_high_distance if you disable these check using `--no_qc`.

- I want to build my own database and use this for this tool. How do I start?

You would need to select reference genomes and define their respective clades. A sketch of the reference genomes should be supplied to AuriClass using the `-r` flag, while a CSV file following the format of `data/clade_config.csv` should be supplied using the `-c` flag. Any related species which should be excluded from the analysis should be defined as "outgroup" in the clade configuration file.

!!! warning
    Before using a new databases, you have to calibrate the thresholds and settings you're planning to use. In the calibration of the default dataset, it was clear that for example the minimal kmer coverage (`-m` in both AuriClass and Mash sketch) has a big influence on the exact value of mash distance. Other factors that should be considered are kmer size (`-k`) and sketch size (`-s`).

If you're missing certain reference genomes or if new clades need to be added, you can also open an [issue on GitHub](https://github.com/RIVM-bioinformatics/auriclass/issues).