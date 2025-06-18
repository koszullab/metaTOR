# Test for CLI tools of metator, MGE related


from pathlib import Path
import metator.commands as mtc
import pandas as pd
import pytest
import os


# Use global variables for temporary tmp/out folders
@pytest.fixture(scope="session")
def tmp_dir(tmp_path_factory):
    p = tmp_path_factory.mktemp("tmp")
    os.makedirs(p, exist_ok=True)
    return str(p)


global_args = {
    "NETWORK": "tests_data/network.txt",
    "CONTIGS": "tests_data/contig_data_final.txt",
    "BINNING": "tests_data/binning.txt",
    "SUBSET_FASTA": "tests_data/subset.assembly.fa",
    "SUBSET_PAIRS": "tests_data/subset.pairs.gz",
    "SUBSET_MGES": "tests_data/subset.mges.txt",
}


def test_mge(tmp_dir):
    """Test the metator mge command with subset data from pipeline output."""
    args = (
        "-c {CONTIGS} -a {SUBSET_FASTA} -b {BINNING} -m {SUBSET_MGES} "
        "-o {OUT_TEST} --no-clean-up --threshold-bin 0.5 {SUBSET_PAIRS}"
    ).format(
        OUT_TEST=Path(tmp_dir, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator mge", *args.split(" "))
    proc = mtc.Mge(args.split(" "), {})
    proc.execute()

