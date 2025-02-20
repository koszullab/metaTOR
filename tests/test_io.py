# Test for io module.

import metator.io as mio
import pypairix
import pytest
import os
import shutil
from . import LOUVAIN_PATH


def test_check_checkm(): ...


def test_check_fasta_index(): ...


def test_check_is_fasta(): ...


def test_check_louvain_cpp():
    test = mio.check_louvain_cpp(LOUVAIN_PATH)
    assert test


def test_check_pypairix():
    test = mio.check_pypairix()
    assert test


def test_check_pairtools():
    test = mio.check_pairtools()
    assert test


def test_generate_fasta_index(): ...


def test_generate_temp_dir(): ...


def test_get_restriction_site(): ...


def test_get_pairs(): ...


def test_process_ligation_sites(): ...


def test_read_bin_summary(): ...


def test_read_compressed(): ...


def test_read_contig_data(): ...


def test_read_results_checkm(): ...


def test_retrieve_fasta(): ...


def test_sort_pairs(): ...


def test_sort_pairs_pairtools():
    pairfile = "tests_data/outdir/alignment.pairs"
    tmp_dir = "tests_data/out_test_io"
    os.makedirs(tmp_dir, exist_ok=True)
    test_file = os.path.join(tmp_dir, "test.pairs")
    shutil.copyfile(pairfile, test_file)
    mio.sort_pairs_pairtools(test_file, 1, True, False)
    pairs_data = pypairix.open(os.path.join(tmp_dir, f"test_sorted.pairs.gz"))
    shutil.rmtree(tmp_dir)


def test_write_bin_summary(): ...
