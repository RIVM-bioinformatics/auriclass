import argparse
from pathlib import Path

from .general import check_number_within_range
from .version import __description__, __package_name__, __version__


def auriclass_arg_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__description__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required arguments
    required_args = parser.add_argument_group("REQUIRED")
    required_args.add_argument(
        "read_file_paths",
        help="Paths to read files",
        nargs="+",
    )

    # Main arguments
    main_args = parser.add_argument_group("Main arguments")
    main_args.add_argument(
        "-n",
        "--name",
        help="Name of isolate",
        default="isolate",
    )
    main_args.add_argument(
        "-o",
        "--output_report_path",
        help="Path to output report",
        default="report.tsv",
        type=Path,
    )
    main_args.add_argument(
        "--fastq",
        help="Input files are fastq files",
        action="store_true",
    )
    main_args.add_argument(
        "--fasta",
        help="Input files are fasta files",
        action="store_true",
    )
    main_args.add_argument(
        "--no_qc",
        help="Skip extended QC",
        action="store_true",
        dest="no_qc",
    )
    main_args.add_argument(
        "--log_file_path",
        help="Path to log file",
        type=Path,
    )
    main_args.add_argument(
        "--verbose",
        help="Verbose output",
        action="store_true",
    )
    main_args.add_argument(
        "--debug",
        help="Very verbose output",
        action="store_true",
    )
    main_args.add_argument(
        "--version",
        action="version",
        version=f"{__package_name__} {__version__}",
    )

    # QC arguments
    qc_args = parser.add_argument_group("QC arguments")

    qc_args.add_argument(
        "--expected_genome_size",
        help="Expected genome size range. Defaults 11.4-14.6 Mb are based on"
        " 150 NCBI genomes and take mash genome size overestimation into account.",
        default=[11_400_000, 14_900_000],
        nargs=2,
        type=check_number_within_range(0, 100_000_000),
    )
    qc_args.add_argument(
        "--non_candida_threshold",
        help="If the minimal distance from a reference sample is above this threshold, the sample might not be a Candida sp.",
        default=0.01,
        type=check_number_within_range(0, 1),
    )
    qc_args.add_argument(
        "--high_dist_threshold",
        help="If the minimal distance from a reference sample is above this threshold, a warning is emitted. See the docs for more info.",
        default=0.003,
        type=check_number_within_range(0, 1),
    )

    # Other arguments
    other_args = parser.add_argument_group(
        "Other arguments\nNOTE: Only change these settings if you are doing something special.\nNOTE: This will require rebuilding the reference sketch and recalibration of thresholds!"
    )
    other_args.add_argument(
        "-r",
        "--reference_sketch_path",
        help="Path to reference sketch",
        default="data/Candida_auris_clade_references.msh",
    )
    other_args.add_argument(
        "-c",
        "--clade_config_path",
        help="Path to clade config",
        default="data/clade_config.csv",
    )
    other_args.add_argument(
        "-k",
        "--kmer_size",
        help="Kmer size",
        default=27,
        type=check_number_within_range(1, 32),
    )
    other_args.add_argument(
        "-s",
        "--sketch_size",
        help="Sketch size",
        default=50_000,
        type=check_number_within_range(1000, 1_000_000),
    )
    other_args.add_argument(
        "-m",
        "--minimal_kmer_coverage",
        help="Minimal kmer coverage",
        default=3,
        type=check_number_within_range(1, 100),
    )
    args = parser.parse_args()
    return args
