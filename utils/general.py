import argparse
import logging
import subprocess
from pathlib import Path
from typing import Callable, List, Union

import pyfastx


def add_tag(tag: str, lines: str) -> str:
    """
    Add a tag to a string or list of strings.

    Parameters
    ----------
    tag : str
        The tag to add.
    lines : str or list of str
        The string or list of strings to add the tag to.

    Returns
    -------
    str
        The string with the tag added to each line.
    """
    full_tag = f"[{tag}]"
    if lines:
        return "\n".join(
            [f"{full_tag} {line}" for line in lines.split("\n") if line != ""]
        )
    else:
        return f"{full_tag}"


def check_number_within_range(
    minimum: float = 0, maximum: float = 1
) -> Union[Callable[[str], str], argparse.FileType]:
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


def is_fastq(filepath: str) -> None:
    """
    Check if a file is a fastq file.

    Parameters
    ----------
    filepath : str
        The path to the file to check.

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If the file is not a fastq file.
    """
    try:
        pyfastx.Fastq(filepath, build_index=False)
    except RuntimeError:
        raise RuntimeError(
            f"Input file {filepath} is not a fastq file. AuriClass expects fastq files"
        )


# def is_fasta(file: str) -> bool:
#     """
#     Check if a file is a fasta file.

#     Parameters
#     ----------
#     filepath : str
#         The path to the file to check.

#     Returns
#     -------
#     None

#     Raises
#     ------
#     RuntimeError
#         If the file is not a fasta file.
#     """
#     try:
#         pyfastx.Fasta(file, build_index=False)
#         return True
#     except RuntimeError:
#         return False


def validate_input_files(list_of_files: List[str]) -> None:
    """
    Check if files in a list exist.

    Parameters
    ----------
    list_of_files : list of str
        The list of filepaths to check.

    Returns
    -------
    None

    Raises
    ------
    FileNotFoundError
        If any of the files do not exist.
    """
    for filepath in list_of_files:
        if not Path(filepath).exists():
            raise FileNotFoundError(f"Required input file {filepath} does not exist")


def validate_argument_logic(args: argparse.Namespace) -> argparse.Namespace:
    """
    Check if arguments are valid, based on predefined logic.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments to check.

    Returns
    -------
    argparse.Namespace
        The arguments, if they are valid. Modified if necessary.

    Raises
    ------
    ValueError
        If the arguments are not valid.
    """
    # Convert to float
    args.expected_genome_size = [
        float(args.expected_genome_size[0]),
        float(args.expected_genome_size[1]),
    ]
    # Check if specified genome size range is valid
    if args.expected_genome_size[0] > args.expected_genome_size[1]:
        raise ValueError(
            "Expected genome size range is invalid: lower bound is higher than upper bound"
        )
    elif (args.expected_genome_size[0] < 100) & (args.expected_genome_size[1] < 100):
        logging.warning(
            f"Expected genome size range boundaries {args.expected_genome_size} are below 100: treating these as Mbp instead of bp"
        )
        args.expected_genome_size = [
            args.expected_genome_size[0] * 1_000_000,
            args.expected_genome_size[1] * 1_000_000,
        ]
    return args


def check_dependencies() -> None:
    """
    Check if dependencies are installed.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    FileNotFoundError
        If mash is not installed.
    """
    try:
        subprocess.call(
            ["mash", "-h"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except FileNotFoundError:
        raise FileNotFoundError("Mash is not installed")
