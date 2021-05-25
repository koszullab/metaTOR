# Test for validation module

import metator.io as mio
import metator.validation as mtv
import pytest
import os
import re
import shutil

ALGORITHM = ("alg", ["louvain", "leiden"])


def test_checkm():
    ...


def test_compare_bins():
    checkm_summary = mtv.compare_bins(
        overlapping_checkm_file="tests_data/outdir2/overlapping_checkm_results.txt",
        overlapping_taxonomy_file="tests_data/outdir2/overlapping_checkm_taxonomy.txt",
        recursive_checkm_file="tests_data/outdir2/recursif_checkm_results.txt",
        recursive_taxonomy_file="tests_data/outdir2/recursif_checkm_taxonomy.txt",
    )
    over = mio.read_results_checkm(
        "tests_data/outdir2/overlapping_checkm_results.txt",
        "tests_data/outdir2/overlapping_checkm_taxonomy.txt",
    )
    rec = mio.read_results_checkm(
        "tests_data/outdir2/recursif_checkm_results.txt",
        "tests_data/outdir2/recursif_checkm_taxonomy.txt",
    )

    # Assert there is the rigth numbers of keys for a bin:
    keys = list(checkm_summary.keys())
    assert len(checkm_summary[keys[0]]) == 10

    # Assert there are no more keys then available:
    test_over = [re.search("_0", key) is not None for key in keys]
    assert sum(test_over) <= len(over)
    assert len(checkm_summary) - sum(test_over) <= len(rec)

    # Assert there are no duplicates  (one bin from overlapping and from
    # recursive)
    over_id = []
    rec_id = []
    for i in range(len(test_over)):
        if test_over[i]:
            over_id.append(int(keys[i].split("_")[1]))
        else:
            rec_id.append(int(keys[i].split("_")[1]))
    assert not (set(over_id) & set(rec_id))


@pytest.mark.parametrize(*ALGORITHM)
def test_recursive_clustering(alg):
    """Crash test for the partition of contaminated bins."""
    os.makedirs("tests_data/out_test/", exist_ok=True)
    os.makedirs("tests_data/out_test/tmp/", exist_ok=True)
    os.makedirs("tests_data/out_test/recursive_bin/", exist_ok=True)
    mtv.recursive_clustering(
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
        cluster_matrix=True,
        size=1000000,
        threads=1,
    )
    shutil.rmtree("tests_data/out_test/")


def give_results_info():
    ...


def test_recursive_decontamination():
    ...


def update_contigs_data_recursive():
    ...


def write_bin_contigs():
    ...
