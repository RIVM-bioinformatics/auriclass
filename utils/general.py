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


def is_fastq(file: str) -> bool:
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
        pyfastx.Fastq(file, build_index=False)
        return True
    except RuntimeError:
        return False


def is_fasta(file: str) -> bool:
    """
    Check if a file is a fasta file.

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
        If the file is not a fasta file.
    """
    try:
        pyfastx.Fasta(file, build_index=False)
        return True
    except RuntimeError:
        return False


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


def guess_input_type(list_of_file_paths: List[str]) -> str:
    """
    Guess the input type of a list of files.

    Parameters
    ----------
    list_of_file_paths : list of str
        The list of filepaths to check.

    Returns
    -------
    str
        The input type, either "fastq" or "fasta".

    Raises
    ------
    ValueError
        If any of the input files can be parsed as both fastq and fasta.

        If any of the input files cannot be parsed as either fastq or fasta.

        If the input files are a mix of fastq and fasta files.

        If no input files were found.
    """
    fastq_count = 0
    fasta_count = 0
    for file_path in list_of_file_paths:
        if is_fastq(file_path):
            if is_fasta(file_path):
                raise ValueError(
                    f"Input file {file_path} can be parsed as both fastq and fasta. Please specify --fastq or --fasta"
                )
            else:
                fastq_count += 1
        elif is_fasta(file_path):
            fasta_count += 1
        else:
            raise ValueError(f"Input file {file_path} is not a fastq or fasta file")

    if fastq_count > 0 and fasta_count > 0:
        raise ValueError(f"Input files are a mix of fastq and fasta files")
    elif fastq_count > 0:
        return "fastq"
    elif fasta_count > 0:
        return "fasta"
    else:
        raise ValueError(f"No input files were found")


def confirm_input_type(list_of_file_paths, input_type):
    """
    Confirms that input files are of the specified type.

    Parameters
    ----------
    list_of_file_paths : list of str
        The list of filepaths to check.
    input_type : str
        The input type to check.

    Returns
    -------
    None

    Raises
    ------
    None

    Notes
    -----
    This function does not raise any exceptions, but it does log a warning if any of the input files are not of the specified type.
    This is intended to run if the user specifies the input type with --fastq or --fasta.
    """
    for file_path in list_of_file_paths:
        if input_type == "fastq":
            if not is_fastq(file_path):
                logging.warning(
                    f"Input file {file_path} cannot be parsed as a fastq file, please check if --fastq is appropriate"
                )
        elif input_type == "fasta":
            if not is_fasta(file_path):
                logging.warning(
                    f"Input file {file_path} cannot be parsed as a fasta file, please check if --fasta is appropriate"
                )
