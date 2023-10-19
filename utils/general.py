import argparse
import logging
import subprocess
from pathlib import Path

import pyfastx


def add_tag(tag, lines):
    """
    Prepend a tag to text explaining from which analysis the message came.
    """
    full_tag = f"[{tag}]"
    if lines:
        return "\n".join(
            [f"{full_tag} {line}" for line in lines.split("\n") if line != ""]
        )
    else:
        return f"{full_tag}"


def check_number_within_range(minimum=0, maximum=1):
    """
    Creates a function to check whether a numeric value is within a range, inclusive.

    The generated function can be used by the `type` parameter in argparse.ArgumentParser.
    See https://stackoverflow.com/a/53446178.

    Args:
        value: the numeric value to check.
        minimum: minimum of allowed range, inclusive.
        maximum: maximum of allowed range, inclusive.

    Returns:
        A function which takes a single argument and checks this against the range.

    Raises:
        argparse.ArgumentTypeError: if the value is outside the range.
        ValueError: if the value cannot be converted to float.
    """

    def generated_func_check_range(value: str) -> str:
        value_f = float(value)
        if (value_f < minimum) or (value_f > maximum):
            raise argparse.ArgumentTypeError(
                f"Supplied value {value} is not within expected range {minimum} to {maximum}."
            )
        return str(value)

    return generated_func_check_range


def is_fastq(filepath):
    try:
        pyfastx.Fastq(filepath, build_index=False)
    except RuntimeError:
        raise RuntimeError(
            f"Input file {filepath} is not a fastq file. AuriClass expects fastq files"
        )


# def is_fasta(file):
#     try:
#         pyfastx.Fasta(file, build_index=False)
#         return True
#     except RuntimeError:
#         return False


def validate_input_files(list_of_files):
    for filepath in list_of_files:
        if not Path(filepath).exists():
            raise FileNotFoundError(f"Required input file {filepath} does not exist")


def validate_argument_logic(args):
    # Check if specified genome size range is valid
    if args.genome_size_range[0] > args.genome_size_range[1]:
        raise ValueError(
            "Expected genome size range is invalid: lower bound is higher than upper bound"
        )
    elif (args.genome_size_range[0] < 100) & (args.genome_size_range[1] < 100):
        logging.warning(
            f"Expected genome size range boundaries {args.genome_size_range} are below 100: treating these as Mbp instead of bp"
        )
        args.genome_size_range = [
            args.genome_size_range[0] * 1_000_000,
            args.genome_size_range[1] * 1_000_000,
        ]
    return args


def check_dependencies():
    """
    Check if dependencies are installed
    """
    try:
        subprocess.call(
            ["mash", "-h"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except FileNotFoundError:
        raise FileNotFoundError("Mash is not installed")
