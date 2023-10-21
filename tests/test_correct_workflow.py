import argparse
import os
import tempfile

import pytest

from auriclass.classes import FastaAuriclass, FastqAuriclass
from auriclass.general import (
    check_dependencies,
    is_fastq,
    validate_argument_logic,
    validate_input_files,
)

os.makedirs("tmp_data", exist_ok=True)


mash_output_to_dict_fastq = {
    "Reference": {0: "tests/data/NC_001416.1.fasta", 1: "tests/data/NC_001604.1.fasta"},
    "Query": {0: "tests/data/NC_001416.1_1.fq.gz", 1: "tests/data/NC_001416.1_1.fq.gz"},
    "Distance": {0: 9.55405e-06, 1: 1.0},
    "P-value": {0: 0, 1: 1},
    "Matching-hashes": {0: "48451/48476", 1: "0/50000"},
    "Clade": {0: "Lambda phage", 1: "outgroup"},
}
error_bounds_text_ref = "\nParameters (run with -h for details):\n   k:   27\n   p:   0.99\n\n\tMash distance\nSketch\t0.05\t0.1\t0.15\t0.2\t0.25\t0.3\t0.35\t0.4\n100\t0.0306864\tinf\tinf\tinf\tinf\tinf\tinf\tinf\n500\t0.00994399\t0.0280701\tinf\tinf\tinf\tinf\tinf\tinf\n1000\t0.00677921\t0.0181803\t0.0545726\tinf\tinf\tinf\tinf\tinf\n5000\t0.0029432\t0.00752207\t0.0165713\t0.0384649\tinf\tinf\tinf\tinf\n10000\t0.00204623\t0.00517551\t0.0110846\t0.0266814\t0.0654553\tinf\tinf\tinf\n50000\t0.0008979\t0.0022223\t0.00466358\t0.0097208\t0.0223838\t0.0493898\tinf\tinf\n100000\t0.000633326\t0.00155922\t0.00326618\t0.00666368\t0.0141222\t0.0343733\tinf\tinf\n500000\t0.000281841\t0.000691558\t0.00143912\t0.00291027\t0.0058613\t0.0126052\t0.0289633\tinf\n1000000\t0.000199108\t0.00048843\t0.00101565\t0.00203733\t0.00412576\t0.00839607\t0.0183086\t0.0453242\n\n\tScreen distance\nSketch\t0.05\t0.1\t0.15\t0.2\t0.25\t0.3\t0.35\t0.4\n100\t0.0202309\t0.0568091\t0.85\t0.8\t0.75\t0.7\t0.65\t0.6\n500\t0.00751972\t0.0196909\t0.055602\t0.8\t0.75\t0.7\t0.65\t0.6\n1000\t0.00517769\t0.0123854\t0.0349441\t0.8\t0.75\t0.7\t0.65\t0.6\n5000\t0.00228269\t0.00515208\t0.011506\t0.0321089\t0.75\t0.7\t0.65\t0.6\n10000\t0.00160085\t0.00359078\t0.00775975\t0.018176\t0.75\t0.7\t0.65\t0.6\n50000\t0.000707304\t0.00156526\t0.00337465\t0.00742041\t0.0205406\t0.7\t0.65\t0.6\n100000\t0.00049869\t0.0011065\t0.0023513\t0.00516381\t0.0123875\t0.0471479\t0.65\t0.6\n500000\t0.000222472\t0.00049172\t0.00103321\t0.00226061\t0.0052637\t0.0140566\t0.65\t0.6\n1000000\t0.00015706\t0.000346666\t0.000728674\t0.00158227\t0.00365475\t0.00918557\t0.0349265\t0.6\n\n"

mash_output_to_dict_fasta = {
    "Reference": {0: "tests/data/NC_001416.1.fasta", 1: "tests/data/NC_001604.1.fasta"},
    "Query": {
        0: "tests/data/NC_001416.1.fasta.gz",
        1: "tests/data/NC_001416.1.fasta.gz",
    },
    "Distance": {0: 0, 1: 1},
    "P-value": {0: 0, 1: 1},
    "Matching-hashes": {0: "48476/48476", 1: "0/50000"},
    "Clade": {0: "Lambda phage", 1: "outgroup"},
}


def test_fastq():
    """
    Test if ValueError is raised when min_genome_size is larger than max_genome_size
    """
    args = argparse.Namespace(
        name="test",
        read_file_paths=[
            "tests/data/NC_001416.1_1.fq.gz",
            "tests/data/NC_001416.1_2.fq.gz",
        ],
        output_report_path="tmp_data/test_report.tsv",
        reference_sketch_path="tests/data/ref_sketch.msh",
        expected_genome_size=(40_000, 60_000),
        kmer_size=27,
        sketch_size=50_000,
        minimal_kmer_coverage=3,
        clade_config_path="tests/data/clade_config.csv",
        non_candida_threshold=0.1,
        high_dist_threshold=0.003,
        no_qc=False,
    )

    # Validate inputs before creating object
    validate_input_files(args.read_file_paths)
    validate_input_files([args.reference_sketch_path])
    validate_input_files([args.clade_config_path])

    args = validate_argument_logic(args)

    # Check dependencies
    check_dependencies()

    testsample = FastqAuriclass(
        name=args.name,
        read_paths=args.read_file_paths,
        output_report_path=args.output_report_path,
        reference_sketch_path=args.reference_sketch_path,
        genome_size_range=args.expected_genome_size,
        kmer_size=args.kmer_size,
        sketch_size=args.sketch_size,
        minimal_kmer_coverage=args.minimal_kmer_coverage,
        clade_config_path=args.clade_config_path,
        non_candida_threshold=args.non_candida_threshold,
        high_dist_threshold=args.high_dist_threshold,
        no_qc=args.no_qc,
    )

    # Sketch query genome using tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        testsample.query_sketch_path = f"{tmpdir}/tmpfile.msh"
        testsample.sketch_fastq_query()

        # Run mash screen
        testsample.run_mash_dist()
        assert testsample.mash_output.to_dict() == mash_output_to_dict_fastq

    # Check if genome size is within expected range
    testsample.check_genome_size()
    assert testsample.estimated_genome_size == 48454.7
    assert testsample.qc_genome_size == ""

    # Process results
    testsample.select_clade()
    assert testsample.clade == "Lambda phage"
    assert testsample.minimal_distance == 9.55405e-06
    assert testsample.closest_sample == "tests/data/NC_001416.1.fasta"

    # Check if all samples not above certain threshold indicating other species
    probably_candida = testsample.check_non_candida()
    assert probably_candida
    assert testsample.qc_species == ""

    if probably_candida:
        # Check if clade is not "outgroup"
        probably_cauris = testsample.check_for_outgroup()
        assert probably_cauris
        assert testsample.qc_other_candida == ""

        if probably_cauris:
            # Check if closest sample is above 0.005 distance --> new clade?
            testsample.check_high_dist()
            assert testsample.qc_high_distance == ""

            # Check error bounds and check number of samples within error bounds
            error_bounds_text = testsample.get_error_bounds()
            assert error_bounds_text == error_bounds_text_ref

            testsample.process_error_bounds(error_bounds_text)
            assert testsample.error_bound == 0.0008979
            testsample.compare_with_error_bounds()
            assert testsample.distances == [1.0]
            assert testsample.samples_within_error_bound == 0
            assert testsample.qc_multiple_hits == ""

    # Save report
    testsample.save_report()


def test_fasta():
    """
    Test if ValueError is raised when min_genome_size is larger than max_genome_size
    """
    args = argparse.Namespace(
        name="test",
        read_file_paths=[
            "tests/data/NC_001416.1.fasta.gz",
        ],
        output_report_path="tmp_data/test_report.tsv",
        reference_sketch_path="tests/data/ref_sketch.msh",
        expected_genome_size=(40_000, 60_000),
        kmer_size=27,
        sketch_size=50_000,
        minimal_kmer_coverage=3,
        clade_config_path="tests/data/clade_config.csv",
        non_candida_threshold=0.1,
        high_dist_threshold=0.005,
        no_qc=False,
    )

    # Validate inputs before creating object
    validate_input_files(args.read_file_paths)
    validate_input_files([args.reference_sketch_path])
    validate_input_files([args.clade_config_path])

    args = validate_argument_logic(args)

    # Check dependencies
    check_dependencies()

    testsample = FastaAuriclass(
        name=args.name,
        read_paths=args.read_file_paths,
        output_report_path=args.output_report_path,
        reference_sketch_path=args.reference_sketch_path,
        genome_size_range=args.expected_genome_size,
        kmer_size=args.kmer_size,
        sketch_size=args.sketch_size,
        minimal_kmer_coverage=args.minimal_kmer_coverage,
        clade_config_path=args.clade_config_path,
        non_candida_threshold=args.non_candida_threshold,
        high_dist_threshold=args.high_dist_threshold,
        no_qc=args.no_qc,
    )

    # Sketch query genome using tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        testsample.query_sketch_path = f"{tmpdir}/tmpfile.msh"
        testsample.sketch_fasta_query()

        # Run mash screen
        testsample.run_mash_dist()
        assert testsample.mash_output.to_dict() == mash_output_to_dict_fasta

    # Check if genome size is within expected range
    testsample.parse_genome_size()
    testsample.check_genome_size()
    assert testsample.estimated_genome_size == 48502
    assert testsample.qc_genome_size == ""

    # Process results
    testsample.select_clade()
    assert testsample.clade == "Lambda phage"
    assert testsample.minimal_distance == 0
    assert testsample.closest_sample == "tests/data/NC_001416.1.fasta"

    # Check if all samples not above certain threshold indicating other species
    probably_candida = testsample.check_non_candida()
    assert probably_candida
    assert testsample.qc_species == ""

    if probably_candida:
        # Check if clade is not "outgroup"
        probably_cauris = testsample.check_for_outgroup()
        assert probably_cauris
        assert testsample.qc_other_candida == ""

        if probably_cauris:
            # Check if closest sample is above 0.005 distance --> new clade?
            testsample.check_high_dist()
            assert testsample.qc_high_distance == ""

            # Check error bounds and check number of samples within error bounds
            error_bounds_text = testsample.get_error_bounds()
            assert error_bounds_text == error_bounds_text_ref

            testsample.process_error_bounds(error_bounds_text)
            assert testsample.error_bound == 0.0008979
            testsample.compare_with_error_bounds()
            assert testsample.distances == [1.0]
            assert testsample.samples_within_error_bound == 0
            assert testsample.qc_multiple_hits == ""

    # Save report
    testsample.save_report()
