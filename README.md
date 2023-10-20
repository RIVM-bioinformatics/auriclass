
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://github.com/RIVM-bioinformatics/auriclass/blob/master/LICENSE)
![GitHub release](https://img.shields.io/github/v/release/rivm-bioinformatics/auriclass)

# AuriClass: quick estimation of *Candida auris* clade membership

AuriClass is a small tool which predicts *Candida auris* clade based on Mash distances from reference genomes. Currently, AuriClass only supports Fastq input.

## Contents

- [Description](#Description)
- [Installation](#Installation)
- [Examples](#Examples)
- [Usage](#Usage)
- [Methodology](#Methodology)
- [FAQ](#FAQ)
- [Future plans](#Future-plans) 

## Description 

### Input

AuriClass takes an arbitrary number of Fastq files as input. AuriClass assumes all files originate from a single sample.

### Reference data (included)

Two files with reference data are needed to start the analysis. Default files are included:

| File | Description |
|---|---|
| `data/clade_config.csv` | CSV table listing two columns: filename and associated clade. All included species that are not *Candida auris* are listed as "outgroup". |
| `data/Candida_auris_clade_references.msh` | Mash sketch of reference genomes of *Candida auris* and related species. Reference genomes for *C. auris* clades were selected from [Suphavilai *et al*., Discovery of the sixth *Candida auris* clade in Singapore (medRxiv)](https://doi.org/10.1101/2023.08.01.23293435) and include two (nearly) completed genomes per clade, if available from NCBI. For Clade V only one (nearly) complete genome was available from NCBI. Related species were selected from figure 3 of [Muñoz *et al*., Genomic insights into multidrug-resistance, mating and virulence in *Candida auris* and related emerging species (NatComm)](https://doi.org/10.1038/s41467-018-07779-6). Representative genomes were downloaded from NCBI. The sketch was constructed using sketch size 50,000 and k-mer size 27. |

See the FAQ for how to build your own database or change calibrated settings.

### Output

AuriClass outputs a tab-separated text file containing the main output, and a log file which save info written to STDERR.

The columns in the output report have the following information:

| Column | Description |
|---|---|
| Sample | Sample name. Default value is "isolate". |
| Clade | Predicted clade. If the closest sample is not *Candida auris*, this field will contain "not Candida auris". |
| Mash_distance_from_closest_reference | Mash distance from the closest reference genome. |
| QC_decision | "FAIL" is QC_species is "FAIL". "WARN" if any other QC step indicates "WARN". Otherwise "PASS". |
| QC_species | "FAIL" if the mash distance from the closest reference is more than `--non_candida_threshold`. This indicates that the sample is not closely related to any *Candida* spp. |
| QC_other_Candida | "WARN" if the closest reference genome is not a *Candida auris* reference genome. This indicates that the sample is closer to another species than *Candida auris*. |
| QC_genome_size | "WARN" if the estimated genome size is outside the specified range. This is probably not accurate in case you run a diploid sample. |
| QC_multiple_hits | "WARN" if at least one sample is within the 99% Mash error bounds of the closest hit (from `mash bounds -k [KMERSIZE] -p 0.99`). This should not happen if the appropriate settings and database are used. |
| QC_possible_new_clade | "WARN" if the closest reference genome is *Candida auris*, but the Mash distance is higher than the `--new_clade_threshold`. |

## Installation
Currently, AuriClass can only be cloned from GitHub, and dependencies should be installed through `mamba`/`conda`. 

```
# Checkout the repo
git clone https://github.com/RIVM-bioinformatics/auriclass.git
# Install conda environment
cd auriclass
mamba env create -f env.yaml
# Activate conda environment
conda activate env_auriclass
```

A Conda package in the bioconda channels is planned.

## Examples

### Minimal example
This will create two files:
- `report.tsv` which contains the clade prediction, closest reference sample and QC checks. The report contains the default "isolate" as sample name.
- `report.YYYY-mm-dd_HH-MM-SS.log` which contains messages written to STDERR with additional data.

```
python auriclass.py Candida_auris_R1.fq.gz Candida_auris_R2.fq.gz
```

### Specifying output

This will create `report_Candida_auris.tsv` and `report_Candida_auris.log`. The report will contain the sample name "Candida_auris".

```
python auriclass.py --name Candida_auris -o report_Candida_auris.tsv Candida_auris_R1.fq.gz Candida_auris_R2.fq.gz
```

### Forcing full analysis

This forces AuriClass to perform the full analysis on a sample, even if the species quality check fails. This will cause some QC checks to always pass, so use with caution.

```
python auriclass.py --non_candida_threshold 1 Candida_auris_R1.fq.gz Candida_auris_R2.fq.gz
```

## Usage

```
usage: auriclass.py [-h] [-n NAME] [-o OUTPUT_REPORT_PATH] [-t N_THREADS] [--log_file_path LOG_FILE_PATH] [--verbose] [--debug] [--version] [--expected_genome_size EXPECTED_GENOME_SIZE EXPECTED_GENOME_SIZE]
                    [--non_candida_threshold NON_CANDIDA_THRESHOLD] [--new_clade_threshold NEW_CLADE_THRESHOLD] [-r REFERENCE_SKETCH_PATH] [-c CLADE_CONFIG_PATH] [-k KMER_SIZE] [-s SKETCH_SIZE] [-m MINIMAL_KMER_COVERAGE]
                    read_file_paths [read_file_paths ...]

AuriClass predicts Candida auris clade from Illumina WGS data

options:
  -h, --help            show this help message and exit

REQUIRED:
  read_file_paths       Paths to read files

Main arguments:
  -n NAME, --name NAME  Name of isolate (default: isolate)
  -o OUTPUT_REPORT_PATH, --output_report_path OUTPUT_REPORT_PATH
                        Path to output report (default: report.tsv)
  -t N_THREADS, --n_threads N_THREADS
                        Number of threads. NOTE: multithreading has minimal effect on performance as mash sketch is single-threaded (default: 1)
  --log_file_path LOG_FILE_PATH
                        Path to log file (default: None)
  --verbose             Verbose output (default: False)
  --debug               Very verbose output (default: False)
  --version             show program's version number and exit

QC arguments:
  --expected_genome_size EXPECTED_GENOME_SIZE EXPECTED_GENOME_SIZE
                        Expected genome size range. Defaults 11.4-14.6 Mb are based on 150 NCBI genomes and take mash genome size overestimation into account. (default: [11400000, 14900000])
  --non_candida_threshold NON_CANDIDA_THRESHOLD
                        If the minimal distance from a reference sample is above this threshold, the sample might not be a Candida sp. (default: 0.01)
  --new_clade_threshold NEW_CLADE_THRESHOLD
                        If the minimal distance from a reference sample is above this threshold, the sample might not be a known C. auris clade. (default: 0.005)

Other arguments
NOTE: Only change these settings if you are doing something special.
NOTE: This will require rebuilding the reference sketch and recalibration of thresholds!:
  -r REFERENCE_SKETCH_PATH, --reference_sketch_path REFERENCE_SKETCH_PATH
                        Path to reference sketch (default: data/Candida_auris_clade_references.msh)
  -c CLADE_CONFIG_PATH, --clade_config_path CLADE_CONFIG_PATH
                        Path to clade config (default: data/clade_config.csv)
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size (default: 27)
  -s SKETCH_SIZE, --sketch_size SKETCH_SIZE
                        Sketch size (default: 50000)
  -m MINIMAL_KMER_COVERAGE, --minimal_kmer_coverage MINIMAL_KMER_COVERAGE
                        Minimal kmer coverage (default: 3)
```

## Methodology

Under the hood, AuriClass does the following:

1. Read in arguments using `argparse`, with some basic validation.
2. Perform extra validation on arguments.
3. Validate input files: check if sample files exist and if input sequences are Fastq format.
4. Check whether dependencies are installed (currently only Mash)
5. Create an object of the `AuriClassAnalysis` class, setting attributes parsed from the CLI.
6. Sketch the reads of the input sample and save to a temporary file.
7. Run `mash dist` on sample and reference sketches.
8. Check whether the estimated genome size (parsed from `mash sketch` STDERR) is within the expected range.
9. Based on the closest reference genome, select which clade is most likely.
10. Check whether the sample is close enough to any reference sample to continue the analysis. If not, abort the analysis without error and set clade to "not Candida auris".
11. Check whether the sample is closest to *Candida auris*, if not set clade to "other Candida sp.".
12. Check whether the closest sample is close a clade reference genome. If not, save a warning that this sample is *Candida auris* but might be a novel clade.
13. Get error bounds from `mash bounds`, parse the output and compare this to the distances observed between test sample and references. If the mash distance between test sample and closest reference genome is within the expected error bound of the mash distance between test sample and second closest reference genome, a warning is emitted.
14. Parse QC failures & warnings and save to a TSV file.

## FAQ

- I get the error `Input file {filepath} is not a fastq file. AuriClass expects fastq files`. How do I solve this?

This happens when the `pyfastx` library cannot parse the input file(s) as FASTQ. Please check whether the files are correctly formatted and not corrupted. Plain and gzipped Fastq files should work fine, other compression algorithms might not be supported. Check [pyfastx docs](https://pyfastx.readthedocs.io/en/latest/) for more information.

- I get the error `Mash is not installed`. How do I solve this?

AuriClass relies on the availability of a Mash executable. Dependencies like mash and pyfastx can e.g. be installed through mamba/conda.

- I expect that my sample is *C. auris*, but my test sample returns a FAIL message in the output report for the species check. How can I force AuriClass to run the full analysis anyway?

You can force AuriClass to run the whole analysis always by changing the QC thresholds, e.g. setting `--non_candida_threshold 1`. Use with caution, as an important QC step is disabled this way.

- I only have Fasta data, how can I still type my *C. auris* samples?

At least two options are available:

1. Use the cauris_cladetyper task from TheiaEuk workflow, available from the [Theiagen GitHub](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/task_cauris_cladetyper.wdl). Also see the associated publication https://doi.org/10.3389/fpubh.2023.1198213.
2. Simulate read data from you Fasta files using wgsim, ART or another read simulator. 

Allowing Fasta input files is planned for AuriClass. Before this can be implemented though, extra logic and possibly extra finetuning of settings is required. Therefore, this is not yet supported.

- I want to build my own database and use this for this tool. How do I start?

You would need to select reference genomes and define their respective clades. A sketch of the reference genomes should be supplied to AuriClass using the `-r` flag, while a CSV file following the format of `data/clade_config.csv` should be supplied using the `-c` flag. Any related species which should be excluded from the analysis should be defined as "outgroup" in the clade configuration file.

Before using a new databases, it would be wise to calibrate the thresholds and settings you're planning to use. In the calibration of the default dataset, it was clear that for example the minimal kmer coverage (`-m` in both AuriClass and Mash sketch) has a big influence on the exact value of mash distance. Other factors that should be considered are kmer size (`-k`) and sketch size (`-s`).

If you're missing certain reference genomes or if new clades need to be added, you can also open an [issue on GitHub](https://github.com/RIVM-bioinformatics/auriclass/issues).

## Future plans

1. Make a bioconda package
2. Support Fasta input