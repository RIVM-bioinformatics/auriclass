# Full usage

```bash
usage: auriclass [-h] [-n NAME] [-o OUTPUT_REPORT_PATH] [--fastq] [--fasta] [--no_qc]
                 [--log_file_path LOG_FILE_PATH] [--verbose] [--debug] [--version]
                 [--expected_genome_size EXPECTED_GENOME_SIZE EXPECTED_GENOME_SIZE]
                 [--non_candida_threshold NON_CANDIDA_THRESHOLD]
                 [--high_dist_threshold HIGH_DIST_THRESHOLD] [-r REFERENCE_SKETCH_PATH]
                 [-c CLADE_CONFIG_PATH] [-k KMER_SIZE] [-s SKETCH_SIZE]
                 [-m MINIMAL_KMER_COVERAGE]
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
  --fastq               Input files are fastq files (default: False)
  --fasta               Input files are fasta files (default: False)
  --no_qc               Skip extended QC (default: False)
  --log_file_path LOG_FILE_PATH
                        Path to log file (default: None)
  --verbose             Verbose output (default: False)
  --debug               Very verbose output (default: False)
  --version             show program's version number and exit

QC arguments:
  --expected_genome_size EXPECTED_GENOME_SIZE EXPECTED_GENOME_SIZE
                        Expected genome size range. Defaults 11.4-14.6 Mb are based on 150
                        NCBI genomes and take mash genome size overestimation into
                        account. (default: [11400000, 14900000])
  --non_candida_threshold NON_CANDIDA_THRESHOLD
                        If the minimal distance from a reference sample is above this
                        threshold, the sample might not be a Candida sp. (default: 0.01)
  --high_dist_threshold HIGH_DIST_THRESHOLD
                        If the minimal distance from a reference sample is above this
                        threshold, a warning is emitted. See the docs for more info.
                        (default: 0.003)

Other arguments
NOTE: Only change these settings if you are doing something special.
NOTE: This will require rebuilding the reference sketch and recalibration of thresholds!:
  -r REFERENCE_SKETCH_PATH, --reference_sketch_path REFERENCE_SKETCH_PATH
                        Path to reference sketch (default: )
  -c CLADE_CONFIG_PATH, --clade_config_path CLADE_CONFIG_PATH
                        Path to clade config (default: )
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size (default: 27)
  -s SKETCH_SIZE, --sketch_size SKETCH_SIZE
                        Sketch size (default: 50000)
  -m MINIMAL_KMER_COVERAGE, --minimal_kmer_coverage MINIMAL_KMER_COVERAGE
                        Minimal kmer coverage (default: 3)
```