#!/usr/bin/env python3

import argparse
import logging
import subprocess
import tempfile
from datetime import datetime
from io import StringIO
from pathlib import Path

import pandas as pd

from utils.general import add_tag, check_number_within_range, is_fasta, is_fastq
from version import __description__, __package_name__, __version__


# define class Sample
class AuriClassAnalysis:
    def __init__(
        self,
        name,
        output_report_path,
        read_paths,
        reference_sketch_path,
        kmer_size,
        sketch_size,
        minimal_kmer_coverage,
        n_threads,
        clade_config_path,
        genome_size_range,
        non_candida_threshold,
        new_clade_threshold,
    ):
        self.name = name
        self.output_report_path = output_report_path
        self.read_paths = read_paths
        self.reference_sketch_path = reference_sketch_path
        self.n_threads = n_threads
        self.kmer_size = kmer_size
        self.sketch_size = sketch_size
        self.minimal_kmer_coverage = minimal_kmer_coverage
        # self.probability is used for mash bounds
        # this is not exposed in the CLI, because this should typically not be changed without extensive testing.
        self.probability = 0.99
        self.clade_dict = pd.read_csv(clade_config_path, index_col=0).to_dict(
            orient="dict"
        )
        self.genome_size_range = genome_size_range
        self.non_candida_threshold = non_candida_threshold
        self.new_clade_threshold = new_clade_threshold
        self.qc_decision = None
        self.qc_genome_size = None
        self.qc_other_candida = None
        self.qc_species = None
        self.qc_multiple_hits = None
        self.qc_new_clade = None
        self.query_sketch_path = None
        self.minimal_distance = None
        self.clade = None
        self.samples_within_error_bound = None
        self.error_bound = None
        self.stdout = None
        self.stderr = None
        self.mash_output = None
        self.distances = None

    def validate_argument_logic(self):
        # Check if specified genome size range is valid
        if self.genome_size_range[0] > self.genome_size_range[1]:
            raise ValueError(
                "Expected genome size range is invalid: lower bound is higher than upper bound"
            )
        elif (self.genome_size_range[0] < 100) & (self.genome_size_range[1] < 100):
            logging.warning(
                f"Expected genome size range boundaries {self.genome_size_range} are below 100: treating these as Mbp instead of bp"
            )
            self.genome_size_range = [
                self.genome_size_range[0] * 1_000_000,
                self.genome_size_range[1] * 1_000_000,
            ]

    def validate_input_files(self):
        for filepath in self.read_paths:
            if not Path(filepath).exists():
                raise FileNotFoundError(f"Read file {filepath} does not exist")
            if is_fastq(filepath) == False:
                if is_fasta(filepath):
                    raise ValueError(
                        f"Input file {filepath} is a fasta file. AuriClass expects fastq files"
                    )
                else:
                    raise ValueError(
                        f"Input file {filepath} is not a fastq file. AuriClass expects fastq files"
                    )

        if not Path(self.reference_sketch_path).exists():
            raise FileNotFoundError(
                f"Reference sketch {self.reference_sketch_path} does not exist"
            )

    def check_dependencies(self):
        """
        Check if dependencies are installed
        """
        # check if mash is installed
        try:
            subprocess.call(
                ["mash", "-h"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except FileNotFoundError:
            raise FileNotFoundError("mash is not installed")

    def sketch_query(self):
        """
        Sketch query genome
        """
        # run mash sketch and capture stdout and stderr
        command = [
            "mash",
            "sketch",
            "-r",
            "-m",
            str(self.minimal_kmer_coverage),
            "-o",
            self.query_sketch_path,
            "-k",
            str(self.kmer_size),
            "-s",
            str(self.sketch_size),
            *self.read_paths,
        ]
        logging.info(add_tag("mash sketch", " ".join(command)))
        output = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if "ERROR: Did not find fasta records in" in output.stderr.decode("utf-8"):
            raise ValueError(
                f"Did not find sequence records in {self.read_paths}. Please check if these are valid fastq files"
            )
        # Add stdout and stderr to object
        self.stdout = output.stdout.decode("utf-8")
        stderr = output.stderr.decode("utf-8")
        if stderr:
            for line in stderr.splitlines():
                logging.info(add_tag("mash sketch", line))
        # add_tag_to_stderr("mash sketch", stderr)

        # find line with "Estimated genome size" and get value
        for line in stderr.splitlines():
            if "Estimated genome size" in line:
                self.estimated_genome_size = float(line.split()[-1])
                break
        else:
            raise ValueError(
                "Estimated genome size could not be parsed from mash sketch STDERR"
            )

    def run_mash(self):
        """
        Run mash screen and predict clade
        """
        command = [
            "mash",
            "dist",
            "-p",
            str(self.n_threads),
            self.reference_sketch_path,
            self.query_sketch_path,
        ]
        logging.info(add_tag("mash dist", " ".join(command)))
        output = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Append stdout and stderr to object
        stderr = output.stderr.decode("utf-8")
        if stderr:
            for line in stderr.splitlines():
                logging.info(add_tag("mash dist", line))
        # Read output into pandas dataframe
        df = pd.read_csv(
            StringIO(output.stdout.decode("utf-8")),
            sep="\t",
            header=None,
            names=["Reference", "Query", "Distance", "P-value", "Matching-hashes"],
        )
        df["Clade"] = df["Reference"].map(self.clade_dict["clade"])
        self.mash_output = df

    def check_genome_size(self):
        # Compare estimated genome size with expected genome size range
        if (
            self.estimated_genome_size < self.genome_size_range[0]
            or self.estimated_genome_size > self.genome_size_range[1]
        ):
            logging.warning(
                f"AuriClass estimated genome size of {self.estimated_genome_size} is outside the expected range of {self.genome_size_range}"
            )
            self.qc_genome_size = "WARN: genome size outside expected range"

    def select_clade(self):
        """
        Select closest sample and corresponding clade
        """
        # Get highest hit, by sorting on distance and selecting first hit
        self.closest_sample = self.mash_output.sort_values("Distance").iloc[0][
            "Reference"
        ]
        # Get clade of closest sample
        self.clade = self.clade_dict["clade"][self.closest_sample]
        # Get minimal distance using self.mash_output and closest sample
        self.minimal_distance = self.mash_output.loc[
            self.mash_output["Reference"] == self.closest_sample, "Distance"
        ].values[0]

    def check_non_candida(self):
        # Check if distance is higher than self.non_candida_threshold
        # If higher, this is probably not Candida
        if self.minimal_distance > self.non_candida_threshold:
            logging.warning(
                f"AuriClass found a distance of {self.minimal_distance} to the closest sample, please ensure this is Candida auris"
            )
            self.qc_species = f"WARN: distance {self.minimal_distance} to closest sample is above threshold"
            return False
        else:
            return True

    def check_for_outgroup(self):
        """
        Check if closest sample is outgroup
        """
        if self.clade == "outgroup":
            logging.warning(
                f"AuriClass found a non-Candida auris reference as closest sample, please ensure this is Candida auris"
            )
            self.qc_other_candida = (
                f"WARN: outgroup reference {self.closest_sample} as closest sample"
            )
            return False
        else:
            return True

    def check_possible_new_clade(self):
        # Check if closest sample is above 0.005 distance --> new clade?
        if self.minimal_distance > 0.005:
            logging.warning(
                f"AuriClass found a distance of {self.minimal_distance} to the closest sample, please ensure this is Candida auris"
            )
            self.qc_new_clade = f"WARN: distance {self.minimal_distance} to closest sample is above threshold"

    def get_error_bounds(self):
        """
        Get mash error bounds and return
        """
        command = [
            "mash",
            "bounds",
            "-k",
            str(self.kmer_size),
            "-p",
            str(self.probability),
        ]
        logging.info(add_tag("mash bounds", " ".join(command)))
        output = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stderr = output.stderr.decode("utf-8")
        if stderr:
            for line in stderr.splitlines():
                logging.info(add_tag("mash bounds", line))

        return output.stdout.decode("utf-8")

    def process_error_bounds(self, error_bounds_text):
        """
        Reads text content and takes lines after line which contains "Mash distance" and before "Screen distance" (irrespective of preceding whitespace) .
        Content is read into a pandas dataframe.
        """
        # Split text into lines
        lines = error_bounds_text.splitlines()
        # Find line with "Mash distance"
        for i, line in enumerate(lines):
            if "Mash distance" in line:
                start_line = i + 1
                break
        # Find line with "Screen distance"
        for i, line in enumerate(lines):
            if "Screen distance" in line:
                end_line = i
                break
        # Read lines into dataframe
        df = pd.read_csv(
            StringIO("\n".join(lines[start_line:end_line])),
            sep="\t",
        )
        # Get first column of which column name is higher than min_dist
        for col in df.columns[1:]:
            if float(col) > self.minimal_distance:
                selected_col = str(col)
                break
        self.error_bound = df.loc[
            df["Sketch"] == self.sketch_size, selected_col
        ].values[0]

    def compare_with_error_bounds(self):
        """
        Compare error_bound with distance between highest hit and other samples
        """
        # Get distance between highest hit and other samples
        self.distances = self.mash_output.loc[
            self.mash_output["Reference"] != self.closest_sample, "Distance"
        ].values
        # Check if distance is higher than error_bound
        count_above_threshold = 0
        for distance in self.distances:
            if distance < (self.minimal_distance + self.error_bound):
                logging.debug(
                    add_tag(
                        "compare_with_error_bounds",
                        f"distance {distance} is above threshold: ({self.minimal_distance} + {self.error_bound})",
                    )
                )
                count_above_threshold += 1
            else:
                logging.debug(
                    add_tag(
                        "compare_with_error_bounds",
                        f"distance {distance} is below threshold: ({self.minimal_distance} + {self.error_bound})",
                    )
                )
        self.samples_within_error_bound = count_above_threshold
        # If samples are found within error bound, warn user
        if self.samples_within_error_bound > 0:
            logging.warning(
                f"AuriClass found {self.samples_within_error_bound} sample(s) within the error bound of {self.error_bound} of the closest sample"
            )
            self.qc_multiple_hits = (
                f"WARN: {self.samples_within_error_bound} sample(s) within error bound"
            )

    def save_report(self):
        """
        Write tsv row of sample name, clade, distance and number of samples within error bound, with column names
        """
        # Check if any of the qc attributes contain "WARN"
        if any(
            [
                self.qc_species,
            ]
        ):
            self.qc_decision = "FAIL"
        elif any(
            [
                self.qc_genome_size,
                self.qc_other_candida,
                self.qc_multiple_hits,
                self.qc_new_clade,
            ]
        ):
            self.qc_decision = "WARN"
        else:
            self.qc_decision = "PASS"

        pd.DataFrame(
            [
                [
                    self.name,
                    self.clade,
                    self.minimal_distance,
                    self.qc_decision,
                    self.qc_species,
                    self.qc_other_candida,
                    self.qc_genome_size,
                    self.qc_multiple_hits,
                    self.qc_new_clade,
                ]
            ],
            columns=[
                "Sample",
                "Clade",
                "Mash_distance_from_closest_reference",
                "QC_decision",
                "QC_species",
                "QC_other_Candida",
                "QC_genome_size",
                "QC_multiple_hits",
                "QC_possible_new_clade",
            ],
        ).fillna("-").to_csv(self.output_report_path, sep="\t", index=False)


if __name__ == "__main__":
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
        "-t",
        "--n_threads",
        help="Number of threads.\nNOTE: multithreading has minimal effect on performance as mash sketch is single-threaded",
        default=1,
        type=check_number_within_range(0, 100),
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
        "--new_clade_threshold",
        help="If the minimal distance from a reference sample is above this threshold, the sample might not be a known C. auris clade.",
        default=0.005,
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

    if args.log_file_path:
        log_file_path = args.log_file_path
    else:
        log_file_path = args.output_report_path.with_suffix(
            f".{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
        )

    # Set logging format
    logging.basicConfig(
        filename=log_file_path,
        filemode="w",
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    logging.getLogger().addHandler(logging.StreamHandler())

    # If -V is specified, set logging level to DEBUG
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Create object
    sample = AuriClassAnalysis(
        name=args.name,
        output_report_path=args.output_report_path,
        read_paths=args.read_file_paths,
        reference_sketch_path=args.reference_sketch_path,
        kmer_size=int(args.kmer_size),
        sketch_size=int(args.sketch_size),
        minimal_kmer_coverage=int(args.minimal_kmer_coverage),
        n_threads=int(args.n_threads),
        clade_config_path=args.clade_config_path,
        genome_size_range=[int(size) for size in args.expected_genome_size],
        non_candida_threshold=float(args.non_candida_threshold),
        new_clade_threshold=float(args.new_clade_threshold),
    )

    # Check dependencies
    sample.validate_argument_logic()
    sample.validate_input_files()
    sample.check_dependencies()

    # Sketch query genome using tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        sample.query_sketch_path = f"{tmpdir}/tmpfile.msh"
        sample.sketch_query()

        # Run mash screen
        sample.run_mash()

    # Check if genome size is within expected range
    sample.check_genome_size()

    # Process results
    sample.select_clade()

    # Check if all samples not above certain threshold indicating other species
    probably_candida = sample.check_non_candida()

    if probably_candida:
        # Check if clade is not "outgroup"
        probably_cauris = sample.check_for_outgroup()

        if probably_cauris:
            # Check if closest sample is above 0.005 distance --> new clade?
            sample.check_possible_new_clade()

            # Check error bounds and check number of samples within error bounds
            error_bounds_text = sample.get_error_bounds()
            sample.process_error_bounds(error_bounds_text)
            sample.compare_with_error_bounds()
        else:
            sample.clade = "other Candida sp."
    else:
        sample.clade = "not Candida auris"

    # Save report
    sample.save_report()
