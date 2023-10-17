import os
import subprocess
from pathlib import Path
from sys import path

import pytest

from auriclass import AuriClassAnalysis

os.makedirs("tmp_data", exist_ok=True)


def test_check_dependency_fail():
    """
    Test if the mash dependency can be found
    """
    testsample = AuriClassAnalysis(
        name="test",
        read_paths=["tests/data/NC_001416.1_1.fq.gz", "tests/data/NC_001416.1_2.fq.gz"],
        output_report_path="tmp_data/test_report.tsv",
        reference_sketch_path="tests/data/ref_sketch.msh",
        genome_size_range=(40_000, 60_000),
        kmer_size=27,
        sketch_size=50_000,
        minimal_kmer_coverage=3,
        n_threads=1,
        clade_config_path="tests/data/clade_config.csv",
        non_candida_threshold=0.1,
        new_clade_threshold=0.005,
    )
    testsample.check_dependencies()
