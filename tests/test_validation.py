# Test for validation modules


import metator.validation as mtv
import pytest
import os, shutil

ALGORITHM = ("alg", ["louvain", "leiden"])


def test_checkm():
    ...


def test_compare_bins():
    ...


@pytest.mark.parametrize(*ALGORITHM)
def test_louvain_recursif(alg):
    """Crash test for the partition of contaminated bins."""
    os.makedirs("tests_data/out_test/", exist_ok=True)
    os.makedirs("tests_data/out_test/tmp/", exist_ok=True)
    os.makedirs("tests_data/out_test/recursive_bin/", exist_ok=True)
    mtv.louvain_recursif(
        assembly="tests_data/assembly.fa",
        iterations=5,
        overlapping_parameter=0.9,
        resolution_parameter=1.0,
        outdir="tests_data/out_test/",
        recursive_fasta_dir="tests_data/out_test/recursive_bin",
        algorithm=alg,
        tmpdir="tests_data/out_test/tmp/",
        checkm_file="tests_data/outdir2/overlapping_checkm_results.txt",
        taxonomy_file="tests_data/outdir2/overlapping_checkm_taxonomy.txt",
        contigs_data_file="tests_data/outdir2/contig_data_partition.txt",
        network_file="tests_data/outdir2/network.txt",
        size=1000000,
        threads=1,
    )
    shutil.rmtree("tests_data/out_test/")


def test_recursive_decontamination():
    ...


def update_recursif_louvain():
    ...


def write_bin_contigs():
    ...
