# Test for io module.

import metator.io as mio
import numpy as np
import pypairix
import pytest
import os
import shutil
from scipy import sparse


def test_check_checkm(): ...


def test_check_fasta_index(tmp_path, monkeypatch):
    """Index basenames keep their leading "../"."""
    for suffix in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]:
        (tmp_path / f"idx.{suffix}").touch()
    for suffix in ["amb", "ann", "bwt", "pac", "sa"]:
        (tmp_path / f"ref.{suffix}").touch()
    (tmp_path / "work").mkdir()
    monkeypatch.chdir(tmp_path / "work")
    assert mio.check_fasta_index("../idx") == "../idx"
    assert mio.check_fasta_index(str(tmp_path / "idx")) == str(tmp_path / "idx")
    assert mio.check_fasta_index("../ref", mode="bwa") == "../ref"
    assert mio.check_fasta_index("../missing") is None


def test_check_is_fasta(): ...


def test_check_pypairix():
    test = mio.check_pypairix()
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
    pairfile = "tests_data/alignment.pairs"
    tmp_dir = "tests_data/out_test_io"
    os.makedirs(tmp_dir, exist_ok=True)
    test_file = os.path.join(tmp_dir, "test.pairs")
    shutil.copyfile(pairfile, test_file)
    mio.sort_pairs_pairtools(test_file, 1, True, False)
    pairs_data = pypairix.open(os.path.join(tmp_dir, f"test_sorted.pairs.gz"))
    shutil.rmtree(tmp_dir)


def test_write_bin_summary(): ...


def test_save_sparse_matrix_requires_coo(tmp_path):
    with pytest.raises(ValueError, match="coo"):
        mio.save_sparse_matrix(sparse.csr_matrix(np.eye(2)), tmp_path / "matrix.txt")
