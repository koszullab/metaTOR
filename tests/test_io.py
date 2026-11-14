# Test for io module.

import bz2
import gzip
import os
import shutil
import subprocess as sp
import types
import zipfile

import metator.io as mio
import numpy as np
import pandas as pd
import pypairix
import pytest
from scipy import sparse
from tests import synthetic


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


def test_check_checkm(monkeypatch):
    def missing(*args, **kwargs):
        raise sp.CalledProcessError(127, "checkm")

    monkeypatch.setattr(mio.sp, "check_output", missing)
    assert not mio.check_checkm()
    monkeypatch.setattr(mio.sp, "check_output", lambda *args, **kwargs: b"CheckM v1.2")
    assert mio.check_checkm()


def test_check_fasta_index_other_mode():
    assert mio.check_fasta_index("ref.fa", mode="minimap2") == "ref.fa"


def test_check_pypairix_missing(monkeypatch):
    monkeypatch.setattr(mio, "pypairix", types.SimpleNamespace())
    with pytest.raises(AttributeError):
        mio.check_pypairix()


def test_generate_temp_dir_without_permission(tmp_path, monkeypatch):
    def denied(path):
        raise PermissionError(path)

    monkeypatch.setattr(mio.os, "makedirs", denied)
    with pytest.raises(PermissionError, match="write permission"):
        mio.generate_temp_dir(str(tmp_path))


def test_read_compressed_zip(tmp_path):
    with zipfile.ZipFile(tmp_path / "one.zip", "w") as archive:
        archive.writestr("reads.fq", "@r1\nACGT\n")
    with mio.read_compressed(str(tmp_path / "one.zip")) as handle:
        assert handle.read() == "@r1\nACGT\n"
    with zipfile.ZipFile(tmp_path / "two.zip", "w") as archive:
        archive.writestr("a.fq", "a")
        archive.writestr("b.fq", "b")
    with pytest.raises(IOError, match="single"):
        mio.read_compressed(str(tmp_path / "two.zip"))


def test_retrieve_fasta(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">c1\nACGTACGTAC\n")
    assert mio.retrieve_fasta(str(fasta), "bowtie2", str(tmp_path)) == str(fasta)
    # A bowtie2 index: the fasta is rebuilt with bowtie2-inspect.
    rebuilt = mio.retrieve_fasta("tests_data/assembly", "bowtie2", str(tmp_path))
    assert rebuilt == str(tmp_path / "assembly.fa")
    assert synthetic.read_fasta(rebuilt) == synthetic.read_fasta("tests_data/assembly.fa")
    # A bwa index: the fasta must sit next to it (.fa or .fasta).
    for suffix in ["amb", "ann", "bwt", "pac", "sa"]:
        (tmp_path / f"ref.{suffix}").touch()
    assert mio.retrieve_fasta(str(tmp_path / "ref"), "bwa", str(tmp_path)) == str(fasta)
    fasta.rename(tmp_path / "ref.fasta")
    assert mio.retrieve_fasta(str(tmp_path / "ref"), "bwa", str(tmp_path)) == str(tmp_path / "ref.fasta")
    (tmp_path / "ref.fasta").unlink()
    with pytest.raises(ValueError):
        mio.retrieve_fasta(str(tmp_path / "ref"), "bwa", str(tmp_path))
    with pytest.raises(ValueError):
        mio.retrieve_fasta(str(tmp_path / "nothing"), "bowtie2", str(tmp_path))


def test_save_sparse_matrix(tmp_path):
    matrix = sparse.coo_matrix(([1.5, 2.0], ([0, 1], [1, 2])), shape=(3, 3))
    mio.save_sparse_matrix(matrix, tmp_path / "matrix.txt")
    lines = (tmp_path / "matrix.txt").read_text().splitlines()
    assert lines == ["3\t3\t2", "0\t1\t1.500", "1\t2\t2.000"]


@pytest.mark.parametrize(
    "version, parallel",
    [("sort (GNU coreutils) 9.4", True), ("sort (GNU coreutils) 8.4", False), ("2.3-Apple (154)", False)],
)
def test_sort_pairs(tmp_path, monkeypatch, version, parallel):
    """Pairs sorted by contig names whatever the sort flavour; --parallel only
    for GNU sort >= 8.23."""
    real_popen = sp.Popen
    commands = []

    def popen(cmd, *args, **kwargs):
        commands.append(cmd)
        if cmd == ["sort", "--version"]:
            return real_popen(["echo", version], *args, **kwargs)
        return real_popen(cmd, *args, **kwargs)

    monkeypatch.setattr(mio.sp, "Popen", popen)
    (tmp_path / "in.txt").write_text("#header\ncontig_10\tcontig_2\ncontig_2\tcontig_3\ncontig_1\tcontig_9\n")
    mio.sort_pairs(str(tmp_path / "in.txt"), str(tmp_path / "out.txt"), tmp_dir=str(tmp_path), threads=2)
    assert (tmp_path / "out.txt").read_text() == "contig_1\tcontig_9\ncontig_2\tcontig_3\ncontig_10\tcontig_2\n"
    sort_cmd = commands[-1]
    assert ("--parallel=2" in sort_cmd) == parallel
    assert f"--temporary-directory={tmp_path}" in sort_cmd


def test_import_mges_tables(tmp_path):
    (tmp_path / "binning.txt").write_text("c1_split_00001\tbin_1\nc2\tbin_2\n")
    (tmp_path / "mges.txt").write_text("c2 phage\nc3\n")
    pd.DataFrame({"ID": [1, 2, 3], "Name": ["c1", "c2", "c3"], "Size": [10, 20, 30]}).to_csv(
        tmp_path / "contigs.txt", sep="\t", index=False
    )
    binning = mio.import_anvio_binning(str(tmp_path / "binning.txt"))
    assert binning == {"c1": "bin_1", "c2": "bin_2"}
    mges = mio.import_mges_contigs(str(tmp_path / "mges.txt"))
    assert mges == ["c2", "c3"]
    data, mges_ids = mio.import_contig_data_mges(str(tmp_path / "contigs.txt"), binning, mges)
    assert mges_ids == [1, 2]
    assert list(data.Final_bin) == ["bin_1", "bin_2", "ND"]
    assert list(data.Binned) == [True, True, False]
    assert list(data.MGE) == [False, True, True]


@pytest.mark.parametrize("shotgun", [True, False])
def test_bin_summary_round_trip(tmp_path, shotgun):
    summary = {
        "bin_1": {"Length": 1000, "Weighted completeness": 0.912345, "HiC_abundance": 1.234567},
        "bin_2": {"Length": 2000, "Weighted completeness": 0.5, "HiC_abundance": 0.1},
    }
    if shotgun:
        for i, values in enumerate(summary.values()):
            values["SG_abundance"] = 10.0 / (i + 1)
    mio.write_bin_summary(summary, str(tmp_path / "bin_summary.txt"))
    table = mio.read_bin_summary(str(tmp_path / "bin_summary.txt"))
    assert list(table.index) == ["bin_1", "bin_2"]
    assert table.loc["bin_1", "Weighted completeness"] == 0.91
    assert table.loc["bin_1", "HiC_abundance"] == 1.2346
    assert ("SG_abundance" in table.columns) == shotgun


def test_read_contig_data(tmp_path):
    pd.DataFrame({"ID": [1, 2], "Name": ["c1", "c2"], "Size": [10, 20]}).to_csv(
        tmp_path / "contig_data.txt", sep="\t", index=False
    )
    data = mio.read_contig_data(str(tmp_path / "contig_data.txt"))
    assert data.loc["c2", "Size"] == 20
    assert list(data.Name) == ["c1", "c2"]


@pytest.mark.parametrize("with_tmp_column", [True, False])
def test_write_mge_data(tmp_path, with_tmp_column):
    data = pd.DataFrame({"Name": ["c1"], "Score": [0.123456]})
    if with_tmp_column:
        data["tmp"] = "c1_1"
    mio.write_mge_data(data, str(tmp_path / "mges.tsv"))
    assert (tmp_path / "mges.tsv").read_text() == "Name\tScore\nc1\t0.12\n"


# Pairs of the mock metagenome (tests/synthetic.py).

needs_pairtools = pytest.mark.skipif(shutil.which("pairtools") is None, reason="pairtools is not installed")


def _count_query(pairs_data, metagenome, first, second):
    """Pairs between two contigs found by pypairix, and in the simulation."""
    sizes = metagenome.truth.Length
    found = sum(1 for _ in pairs_data.query2D(first, 0, sizes[first], second, 0, sizes[second], 1))
    pairs = metagenome.pairs
    expected = (((pairs.chr1 == first) & (pairs.chr2 == second)) | ((pairs.chr1 == second) & (pairs.chr2 == first))).sum()
    return found, expected


@needs_pairtools
def test_get_pairs_data(metagenome, metagenome_files, tmp_path):
    first, second = (c.name for c in metagenome.molecule("genome_2").contigs[:2])
    pairs = tmp_path / "synthetic.pairs"
    shutil.copyfile(metagenome_files["pairs"], pairs)
    # Not indexed: sorted, compressed and indexed next to the input.
    found, expected = _count_query(mio.get_pairs_data(str(pairs)), metagenome, first, second)
    assert found == expected > 0
    sorted_pairs = tmp_path / "synthetic_sorted.pairs.gz"
    assert sorted_pairs.exists() and pairs.exists()
    # Indexed already: opened as it is.
    assert _count_query(mio.get_pairs_data(str(sorted_pairs)), metagenome, first, second) == (found, expected)
    # The sorted file next to the input is reused rather than rebuilt.
    mtime = sorted_pairs.stat().st_mtime_ns
    assert _count_query(mio.get_pairs_data(str(pairs)), metagenome, first, second) == (found, expected)
    assert sorted_pairs.stat().st_mtime_ns == mtime


@needs_pairtools
def test_sort_pairs_pairtools_synthetic(metagenome, metagenome_files, tmp_path, monkeypatch):
    pairs = tmp_path / "synthetic.pairs"
    shutil.copyfile(metagenome_files["pairs"], pairs)
    sorted_pairs = mio.sort_pairs_pairtools(str(pairs), tmp_dir=str(tmp_path))
    assert sorted_pairs == str(tmp_path / "synthetic_sorted.pairs.gz")
    with gzip.open(sorted_pairs, "rt") as handle:
        body = [line.split("\t") for line in handle if not line.startswith("#")]
    assert len(body) == len(metagenome.pairs)
    # Sorted by chr1, chr2, pos1, pos2, as pairix needs.
    keys = [(row[1], row[3], int(row[2]), int(row[4])) for row in body]
    assert keys == sorted(keys)
    # An existing output is kept unless --force.
    with pytest.raises(ValueError):
        mio.sort_pairs_pairtools(str(pairs), tmp_dir=str(tmp_path))
    # Without bgzip, pysam compresses; remove deletes the input.
    real_which = shutil.which
    monkeypatch.setattr(mio.shutil, "which", lambda name: None if name == "bgzip" else real_which(name))
    assert mio.sort_pairs_pairtools(str(pairs), remove=True, force=True, tmp_dir=str(tmp_path)) == sorted_pairs
    assert not pairs.exists()
    first, second = (c.name for c in metagenome.molecule("genome_4").contigs[:2])
    found, expected = _count_query(pypairix.open(sorted_pairs), metagenome, first, second)
    assert found == expected > 0
