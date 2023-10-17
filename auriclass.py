#!/usr/bin/env python3

import argparse
import logging
import subprocess
import sys
import tempfile
from io import StringIO

import pandas as pd


# define class Sample
class AuriClassAnalysis:
    def __init__(
        self,
        name,
        output_report_path,
        fw_read_path,
        rv_read_path,
        reference_sketch_path,
        kmer_size,
        sketch_size,
        minimal_kmer_coverage,
        n_threads,
        clade_config_path,
    ):
        self.name = name
        self.output_report_path = output_report_path
        self.fw_read_path = fw_read_path
        self.rv_read_path = rv_read_path
        self.reference_sketch_path = reference_sketch_path
        self.n_threads = n_threads
        self.kmer_size = kmer_size
        self.sketch_size = sketch_size
        self.minimal_kmer_coverage = minimal_kmer_coverage
        self.probability = 0.99  # used for mash bounds
        self.clade_dict = pd.read_csv(clade_config_path, index_col=0).to_dict(
            orient="dict"
        )
        self.query_sketch_path = None
        self.minimal_distance = None
        self.clade = None
        self.samples_within_error_bound = None
        self.error_bound = None
        self.stdout = None
        self.stderr = None
        self.mash_output = None
        self.distances = None

    def check_dependencies(self):
        """
        Check if dependencies are installed
        """
        # check if mash is installed
        try:
            subprocess.run(
                ["mash", "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        except FileNotFoundError:
            raise FileNotFoundError("mash is not installed")

    def add_tag_to_stderr(self, tag, stderr):
        """
        Prepend stderr text with a tag explaining from which analysis the message came.
        """
        if stderr:
            print(
                "\n".join(
                    [f"[{tag}] {line}" for line in stderr.split("\n") if line != ""]
                ),
                file=sys.stderr,
            )

    def sketch_query(self):
        """
        Sketch query genome
        """
        # run mash sketch and capture stdout and stderr
        output = subprocess.run(
            [
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
                self.fw_read_path,
                self.rv_read_path,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Add stdout and stderr to object
        self.stdout = output.stdout.decode("utf-8")
        self.add_tag_to_stderr("mash sketch", output.stderr.decode("utf-8"))

    def run_mash(self):
        """
        Run mash screen and predict clade
        """
        output = subprocess.run(
            [
                "mash",
                "dist",
                "-p",
                str(self.n_threads),
                self.reference_sketch_path,
                self.query_sketch_path,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Append stdout and stderr to object
        self.add_tag_to_stderr("mash dist", output.stderr.decode("utf-8"))
        # Read output into pandas dataframe
        df = pd.read_csv(
            StringIO(output.stdout.decode("utf-8")),
            sep="\t",
            header=None,
            names=["Reference", "Query", "Distance", "P-value", "Matching-hashes"],
        )
        df["Clade"] = df["Reference"].map(self.clade_dict["clade"])
        self.mash_output = df

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

    def get_error_bounds(self):
        """
        Get mash error bounds and return
        """
        output = subprocess.run(
            [
                "mash",
                "bounds",
                "-k",
                str(self.kmer_size),
                "-p",
                str(self.probability),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.add_tag_to_stderr("mash bounds", output.stderr.decode("utf-8"))
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

    def check_confidence(self):
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
                    f"distance {distance} is above threshold: ({self.minimal_distance} + {self.error_bound})"
                )
                count_above_threshold += 1
            else:
                logging.debug(
                    f"distance {distance} is below threshold: ({self.minimal_distance} + {self.error_bound})"
                )
        self.samples_within_error_bound = count_above_threshold
        # self.samples_within_error_bound = sum(
        #     self.distances > (self.minimal_distance + self.error_bound)
        # )
        if self.samples_within_error_bound > 0:
            logging.warning(
                f"AuriClass found {self.samples_within_error_bound} sample(s) within the error bound of {self.error_bound}"
            )

    def save_report(self):
        """
        Write tsv row of sample name, clade, distance and number of samples within error bound, with column names
        """
        pd.DataFrame(
            [
                [
                    self.name,
                    self.clade,
                    self.mash_output.iloc[0]["Distance"],
                    self.samples_within_error_bound,
                ]
            ],
            columns=[
                "Sample",
                "Clade",
                "Mash_distance_from_ref",
                "Samples_within_error_bound",
            ],
        ).to_csv(self.output_report_path, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Predict clade of Candida auris isolate"
    )
    parser.add_argument(
        "read_file_paths",
        help="Paths to read files",
        nargs=2,
    )
    parser.add_argument(
        "-n",
        "--name",
        help="Name of isolate",
        default="isolate",
    )
    parser.add_argument(
        "-o",
        "--output_report_path",
        help="Path to output report",
        default="report.tsv",
    )
    parser.add_argument(
        "-r",
        "--reference_sketch_path",
        help="Path to reference sketch",
        default="data/Candida_auris_clade_references.msh",
    )
    parser.add_argument(
        "-k",
        "--kmer_size",
        help="Kmer size",
        default=32,
        type=int,
    )
    parser.add_argument(
        "-s", "--sketch_size", help="Sketch size", default=10_000, type=int
    )
    parser.add_argument(
        "-m",
        "--minimal_kmer_coverage",
        help="Minimal kmer coverage",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-t",
        "--n_threads",
        help="Number of threads",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-c",
        "--clade_config_path",
        help="Path to clade config",
        default="data/clade_config.csv",
    )
    args = parser.parse_args()

    # Set logging format
    logging.basicConfig(format="%(message)s")

    # Create object
    sample = AuriClassAnalysis(
        name=args.name,
        output_report_path=args.output_report_path,
        fw_read_path=args.read_file_paths[0],
        rv_read_path=args.read_file_paths[1],
        reference_sketch_path=args.reference_sketch_path,
        kmer_size=args.kmer_size,
        sketch_size=args.sketch_size,
        minimal_kmer_coverage=args.minimal_kmer_coverage,
        n_threads=args.n_threads,
        clade_config_path=args.clade_config_path,
    )

    # Check dependencies
    sample.check_dependencies()

    # Sketch query genome using tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        sample.query_sketch_path = f"{tmpdir}/tmpfile.msh"
        sample.sketch_query()

        # Run mash screen
        sample.run_mash()

    # Process results
    sample.select_clade()

    # Check error bounds and check number of samples within error bounds
    error_bounds_text = sample.get_error_bounds()
    sample.process_error_bounds(error_bounds_text)
    sample.check_confidence()

    # Save report
    sample.save_report()

    print(
        sample.name,
        sample.kmer_size,
        sample.sketch_size,
        sample.minimal_kmer_coverage,
        sample.closest_sample,
        sample.clade,
        sample.minimal_distance,
        "|".join([str(x) for x in sample.distances]),
        sample.error_bound,
        sample.samples_within_error_bound,
        sep="\t",
    )
