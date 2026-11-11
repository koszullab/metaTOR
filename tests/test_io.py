# Test for io module.

import bz2
import gzip

import metator.io as mio
import numpy as np
import pypairix
import pytest
import os
import shutil
from scipy import sparse


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


def test_check_is_fasta(tmp_path):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nACGT\n")
    other = tmp_path / "contigs.txt"
    other.write_text("not a fasta\n")
    assert mio.check_is_fasta(str(fasta))
    assert not mio.check_is_fasta(str(other))
    assert not mio.check_is_fasta(str(tmp_path / "missing.fa"))


def test_check_pypairix():
    test = mio.check_pypairix()
    assert test


def test_generate_temp_dir(tmp_path):
    first = mio.generate_temp_dir(str(tmp_path / "tmp"))
    second = mio.generate_temp_dir(str(tmp_path / "tmp"))
    assert os.path.isdir(first) and os.path.isdir(second)
    assert first != second


def test_get_restriction_site():
    assert mio.get_restriction_site("DpnII") == "GATC"
    assert mio.get_restriction_site("DpnII,HinfI") == "GA.TC|GATC"


def test_read_compressed(tmp_path):
    text = "line 1\nline 2\n"
    (tmp_path / "plain.txt").write_text(text)
    with gzip.open(tmp_path / "file.gz", "wt") as handle:
        handle.write(text)
    with bz2.open(tmp_path / "file.bz2", "wt") as handle:
        handle.write(text)
    for name in ["plain.txt", "file.gz", "file.bz2"]:
        with mio.read_compressed(str(tmp_path / name)) as handle:
            assert handle.read() == text


def test_sort_pairs_pairtools():
    pairfile = "tests_data/alignment.pairs"
    tmp_dir = "out_test_io"
    os.makedirs(tmp_dir, exist_ok=True)
    test_file = os.path.join(tmp_dir, "test.pairs")
    shutil.copyfile(pairfile, test_file)
    sorted_file = mio.sort_pairs_pairtools(test_file, 1, True, False)
    assert sorted_file == os.path.join(tmp_dir, "test_sorted.pairs.gz")
    assert not os.path.exists(test_file)  # remove=True
    with open(pairfile) as pairs:
        n_input = sum(not line.startswith("#") for line in pairs)
    with gzip.open(sorted_file, "rt") as pairs:
        n_sorted = sum(not line.startswith("#") for line in pairs)
    assert n_sorted == n_input
    pypairix.open(sorted_file)  # the index exists and loads


def test_save_sparse_matrix_requires_coo(tmp_path):
    with pytest.raises(ValueError, match="coo"):
        mio.save_sparse_matrix(sparse.csr_matrix(np.eye(2)), tmp_path / "matrix.txt")
