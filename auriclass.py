#!/usr/bin/env python3

import logging
from datetime import datetime

from utils.args import auriclass_arg_parser
from utils.classes import FastaAuriclass, FastqAuriclass
from utils.general import (
    check_dependencies,
    confirm_input_type,
    guess_input_type,
    validate_argument_logic,
    validate_input_files,
)
from utils.version import __description__, __package_name__, __version__


def main() -> None:
    args = auriclass_arg_parser()

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

    args = validate_argument_logic(args)

    # Check dependencies
    check_dependencies()

    # Get input type from args and guess otherwise
    if args.fastq:
        input_type = "fastq"
        confirm_input_type(args.read_file_paths, input_type)
    elif args.fasta:
        input_type = "fasta"
        confirm_input_type(args.read_file_paths, input_type)
    else:
        input_type = guess_input_type(args.read_file_paths)

    # Create object
    if input_type == "fastq":
        fastq_sample = FastqAuriclass(
            name=args.name,
            output_report_path=args.output_report_path,
            read_paths=args.read_file_paths,
            reference_sketch_path=args.reference_sketch_path,
            kmer_size=int(args.kmer_size),
            sketch_size=int(args.sketch_size),
            minimal_kmer_coverage=int(args.minimal_kmer_coverage),
            clade_config_path=args.clade_config_path,
            genome_size_range=[int(size) for size in args.expected_genome_size],
            non_candida_threshold=float(args.non_candida_threshold),
            high_dist_threshold=float(args.high_dist_threshold),
        )

        # Run object
        fastq_sample.run()

    elif input_type == "fasta":
        fasta_sample = FastaAuriclass(
            name=args.name,
            output_report_path=args.output_report_path,
            read_paths=args.read_file_paths,
            reference_sketch_path=args.reference_sketch_path,
            kmer_size=int(args.kmer_size),
            sketch_size=int(args.sketch_size),
            minimal_kmer_coverage=int(args.minimal_kmer_coverage),
            clade_config_path=args.clade_config_path,
            genome_size_range=[int(size) for size in args.expected_genome_size],
            non_candida_threshold=float(args.non_candida_threshold),
            high_dist_threshold=float(args.high_dist_threshold),
        )

        # Run object
        fasta_sample.run()


if __name__ == "__main__":
    main()
