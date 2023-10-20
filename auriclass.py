#!/usr/bin/env python3

import argparse
import logging
import subprocess
import tempfile
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Any, Dict, Hashable, List, Union

import pandas as pd

from utils.general import (
    add_tag,
    check_dependencies,
    check_number_within_range,
    is_fastq,
    validate_argument_logic,
    validate_input_files,
)
from version import __description__, __package_name__, __version__


# define class Sample
class AuriClassAnalysis:
    def __init__(
        self,
        name: str,
        output_report_path: Path,
        read_paths: List[Path],
        reference_sketch_path: Path,
        kmer_size: int,
        sketch_size: int,
        minimal_kmer_coverage: int,
        n_threads: int,
        clade_config_path: Path,
        genome_size_range: List[int],
        non_candida_threshold: float,
        new_clade_threshold: float,
    ) -> None:
        self.name: str = name
        self.output_report_path: Path = output_report_path
        self.read_paths: List[Path] = read_paths
        self.reference_sketch_path: Path = reference_sketch_path
        self.n_threads: int = n_threads
        self.kmer_size: int = kmer_size
        self.sketch_size: int = sketch_size
        self.minimal_kmer_coverage: int = minimal_kmer_coverage
        # self.probability is used for mash bounds
        # this is not exposed in the CLI, because this should typically not be changed without extensive testing.
        self.probability: float = 0.99
        self.clade_dict: Dict[Hashable, Any] = pd.read_csv(
            clade_config_path,
            index_col=0,
            dtype=str,
        ).to_dict(orient="dict")
        self.genome_size_range: List[int] = genome_size_range
        self.non_candida_threshold: float = non_candida_threshold
        self.new_clade_threshold: float = new_clade_threshold
        self.qc_decision: str = ""
        self.qc_genome_size: str = ""
        self.qc_other_candida: str = ""
        self.qc_species: str = ""
        self.qc_multiple_hits: str = ""
        self.qc_new_clade: str = ""
        self.query_sketch_path: Path = Path()
        self.minimal_distance: float = float()
        self.clade: str = ""
        self.samples_within_error_bound: int = int()
        self.error_bound: float = float()
        self.stdout: str = ""
        self.stderr: str = ""
        self.mash_output: pd.DataFrame = pd.DataFrame()
        self.distances: List[float] = [float()]

    def sketch_query(self) -> None:
        """
        Sketches the query fastq files using mash.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If no sequence records are found in the specified fastq file(s).

            If the estimated genome size cannot be parsed from the mash sketch STDERR.

        Notes
        -----
        This function sets the following attributes of the object:
        - stdout: the STDOUT of the mash sketch command
        - estimated_genome_size: the estimated genome size of the current test sample


        The function uses the following attributes of the object:
        - minimal_kmer_coverage: the minimal kmer coverage to use for the sketch
        - query_sketch_path: the path to the query sketch
        - kmer_size: the kmer size to use for the sketch
        - sketch_size: the size of the sketch to create
        - read_paths: the paths to the fastq files to sketch
        """
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
        command_list_of_str = [str(item) for item in command]
        logging.info(add_tag("mash sketch", " ".join(command_list_of_str)))
        output = subprocess.run(
            command_list_of_str,
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

        # find line with "Estimated genome size" and get value
        for line in stderr.splitlines():
            if "Estimated genome size" in line:
                self.estimated_genome_size = float(line.split()[-1])
                break
        else:
            raise ValueError(
                "Estimated genome size could not be parsed from mash sketch STDERR"
            )

    def run_mash(self) -> pd.DataFrame:
        """
        Run Mash dist and save results.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the results of the Mash dist, including the reference and query sequences,
            the distance between them, the p-value, the number of matching hashes, and the predicted clade.

        Notes
        -----
        This function sets the following attributes of the object:
        - mash_output: a pandas DataFrame containing the output of the Mash distance calculation

        The function uses the following attributes of the object:
        - n_threads: the number of threads to use
        - reference_sketch_path: the path to the reference sketch
        - query_sketch_path: the path to the query sketch
        """
        command = [
            "mash",
            "dist",
            "-p",
            str(self.n_threads),
            self.reference_sketch_path,
            self.query_sketch_path,
        ]
        command_list_of_str = [str(item) for item in command]
        logging.info(add_tag("mash dist", " ".join(command_list_of_str)))
        output = subprocess.run(
            command_list_of_str,
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
        self.mash_output = df
        return df

    def check_genome_size(self) -> None:
        """
        Compare the estimated genome size with the expected genome size range.

        If the estimated genome size is outside the expected range, a warning message is logged and the
        `qc_genome_size` attribute is set to "WARN: genome size outside expected range".

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        None

        Notes
        -----
        This function sets the following attributes of the object:
        - qc_genome_size: a warning message if the estimated genome size is outside the expected range

        The function uses the following attributes of the object:
        - estimated_genome_size: the estimated genome size of the current sample
        """
        if (
            self.estimated_genome_size < self.genome_size_range[0]
            or self.estimated_genome_size > self.genome_size_range[1]
        ):
            logging.warning(
                f"AuriClass estimated genome size of {self.estimated_genome_size} is outside the expected range of {self.genome_size_range}"
            )
            self.qc_genome_size = "WARN: genome size outside expected range"

    def select_clade(self) -> None:
        """
        Selects the closest reference sample and corresponding clade.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        None

        Notes
        -----
        This function sets the following attributes of the object:
        - closest_sample: the closest reference sample to the current test sample
        - clade: the clade of the closest reference sample
        - minimal_distance: the minimal distance between the current test sample and the closest reference sample

        The function uses the following attributes of the object:
        - mash_output: a pandas DataFrame containing the output of the Mash distance calculation
        - clade_dict: a dictionary containing the clade information for each reference sample
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

    def check_non_candida(self) -> bool:
        """
        Check if the distance to the closest sample is higher than the non-Candida threshold.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        bool
            True if the distance to the closest sample is lower than or equal to the non-Candida threshold, False otherwise.

        Notes
        -----
        This function sets the following attributes of the object:
        - qc_species: a warning message if the distance to the closest sample is higher than the non-Candida threshold

        The function uses the following attributes of the object:
        - minimal_distance: the minimal distance between the current test sample and the closest reference sample
        - non_candida_threshold: the threshold above which the distance to the closest sample is considered to be non-Candida
        """
        if self.minimal_distance > self.non_candida_threshold:
            logging.warning(
                f"AuriClass found a distance of {self.minimal_distance} to the closest sample, please ensure this is Candida auris"
            )
            self.qc_species = f"WARN: distance {self.minimal_distance} to closest sample is above threshold"
            return False
        else:
            return True

    def check_for_outgroup(self) -> bool:
        """
        Check if the closest sample is defined as outgroup.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        bool
            True if the closest sample is not an outgroup, False otherwise.

        Notes
        -----
        This function sets the following attributes of the object:
        - qc_other_candida: a warning message if the closest sample is defined as outgroup

        The function uses the following attributes of the object:
        - clade: the clade of the closest reference sample
        - closest_sample: the closest reference sample to the current test sample
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

    def check_possible_new_clade(self) -> None:
        """
        Check if the closest sample is above 0.005 distance. If it is, it may indicate a new clade.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        None

        Notes
        -----
        This function sets the following attributes of the object:
        - qc_new_clade: a warning message if the closest sample is above 0.005 distance

        The function uses the following attributes of the object:
        - minimal_distance: the minimal distance between the current test sample and the closest reference sample
        """
        if self.minimal_distance > 0.005:
            logging.warning(
                f"AuriClass found a distance of {self.minimal_distance} to the closest sample, please ensure this is Candida auris"
            )
            self.qc_new_clade = f"WARN: distance {self.minimal_distance} to closest sample is above threshold"

    def get_error_bounds(self) -> str:
        """
        Get error bounds for the current sketch size and kmer size.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        str
            The STDOUT of the mash bounds command.

        Notes
        -----
        The function uses the following attributes of the object:
        - kmer_size: the kmer size to use for the sketch
        - probability: the probability to use for the mash bounds command
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

    def process_error_bounds(self, error_bounds_text: str) -> None:
        """
        Process the error bounds text to get the error bound for the current sketch size and kmer size.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.
        error_bounds_text : str
            The STDOUT of the mash bounds command.

        Returns
        -------
        None

        Notes
        -----
        This function sets the following attributes of the object:
        - error_bound: the error bound for the current sketch size and kmer size


        The function uses the following attributes of the object:
        - sketch_size: the size of the sketch to create
        - minimal_distance: the minimal distance between the current test sample and the closest reference sample
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

    def compare_with_error_bounds(self) -> None:
        """
        Compare the distance between the current test sample and the closest reference sample with the error bound.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        None

        Notes
        -----
        This function sets the following attributes of the object:
        - distances: the distances between the current test sample and all other reference samples
        - samples_within_error_bound: the number of samples within the error bound of the closest sample
        - qc_multiple_hits: a warning message if the number of samples within the error bound of the closest sample is higher than 0

        The function uses the following attributes of the object:
        - mash_output: a pandas DataFrame containing the output of the Mash distance calculation
        - closest_sample: the closest reference sample to the current test sample
        - minimal_distance: the minimal distance between the current test sample and the closest reference sample
        - error_bound: the error bound for the current sketch size and kmer size
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

    def save_report(self) -> None:
        """
        Save the report to a TSV file.

        Parameters
        ----------
        self : object
            The AuriClassAnalysis object.

        Returns
        -------
        None

        Notes
        -----
        This function sets the following attributes of the object:
        - qc_decision: the QC decision based on the QC attributes

        The function uses the following attributes of the object:
        - qc_species: a warning message if the distance to the closest sample is higher than the non-Candida threshold
        - qc_other_candida: a warning message if the closest sample is defined as outgroup
        - qc_genome_size: a warning message if the estimated genome size is outside the expected range
        - qc_multiple_hits: a warning message if the number of samples within the error bound of the closest sample is higher than 0
        - qc_new_clade: a warning message if the closest sample is above 0.005 distance
        - name: the name of the current test sample
        - clade: the clade of the closest reference sample
        - minimal_distance: the minimal distance between the current test sample and the closest reference sample
        - output_report_path: the path to the output report
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

    # Validate inputs before creating object
    validate_input_files(args.read_file_paths)
    validate_input_files([args.reference_sketch_path])
    validate_input_files([args.clade_config_path])

    for read_file_path in args.read_file_paths:
        is_fastq(read_file_path)

    args = validate_argument_logic(args)

    # Check dependencies
    check_dependencies()

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

    # Sketch query genome using tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        sample.query_sketch_path = Path(tmpdir).joinpath("tmpfile.msh")
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
