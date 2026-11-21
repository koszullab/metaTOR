# Test for io module.

import bz2
import gzip
import os
import subprocess as sp
import zipfile

import metator.io as mio
import pandas as pd
import pyfastx
import pytest


def test_check_is_fasta(tmp_path):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nACGT\n")
    other = tmp_path / "contigs.txt"
    other.write_text("not a fasta\n")
    assert mio.check_is_fasta(str(fasta))
    assert not mio.check_is_fasta(str(other))
    assert not mio.check_is_fasta(str(tmp_path / "missing.fa"))
    # A gzipped fasta is a fasta.
    compressed = tmp_path / "contigs.fa.gz"
    with gzip.open(compressed, "wt") as gz:
        gz.write(">c1\nACGT\n")
    assert mio.check_is_fasta(str(compressed))


def test_generate_temp_dir(tmp_path):
    first = mio.generate_temp_dir(str(tmp_path / "tmp"))
    second = mio.generate_temp_dir(str(tmp_path / "tmp"))
    assert os.path.isdir(first) and os.path.isdir(second)
    assert first != second


def test_generate_temp_dir_leaves_the_random_stream_of_igraph(tmp_path):
    """The name of a temporary directory is not drawn from Python's random
    module, which igraph draws from: each command makes one, so a run step
    by step (network, partition, validation) drew another stream than
    pipeline, and partitioned the network otherwise under one seed."""
    import random

    random.seed(0)
    state = random.getstate()
    mio.generate_temp_dir(str(tmp_path / "tmp"))
    assert random.getstate() == state


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


def test_file_compression(tmp_path):
    """The compression of a file from its first bytes, whatever its name, as
    read_compressed finds it."""
    text = b"line 1\n"
    (tmp_path / "plain.gz").write_bytes(text)
    (tmp_path / "gzip.txt").write_bytes(gzip.compress(text))
    (tmp_path / "bzip2.gz").write_bytes(bz2.compress(text))
    with zipfile.ZipFile(tmp_path / "archive.txt", "w") as archive:
        archive.writestr("file.txt", text)
    found = {name: mio.file_compression(str(tmp_path / name)) for name in ["plain.gz", "gzip.txt", "bzip2.gz", "archive.txt"]}
    assert found == {"plain.gz": "uncompressed", "gzip.txt": "gz", "bzip2.gz": "bz2", "archive.txt": "zip"}


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


def test_remove_pyfastx_index(tmp_path):
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">c1\nACGT\n")
    index = tmp_path / "assembly.fa.fxi"
    index.touch()
    mio.remove_pyfastx_index(str(fasta))
    assert fasta.exists() and not index.exists()
    # No index: nothing to remove.
    mio.remove_pyfastx_index(str(fasta))
    assert fasta.exists()


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


def test_rename_1x_columns(logged):
    """The columns of the contig data that metaTOR 1.x named after bins get
    their 2.0 names (for host.py, which reads the tables of 1.x). A table of
    2.0 is returned as it is, and a column of a 2.0 name wins over the 1.x
    column of the same meaning, which is left as it is."""
    old = [
        "Core_bin_ID", "Core_bin_contigs", "Core_bin_size", "Overlapping_bin_ID", "Overlapping_bin_contigs",
        "Overlapping_bin_size", "Recursive_bin_ID", "Recursive_bin_contigs", "Recursive_bin_size", "Final_bin",
    ]  # fmt: skip
    new = [
        "Core_cluster_ID", "Core_cluster_contigs", "Core_cluster_size", "Overlapping_cluster_ID",
        "Overlapping_cluster_contigs", "Overlapping_cluster_size", "Recursive_cluster_ID", "Recursive_cluster_contigs",
        "Recursive_cluster_size", "Final_cluster",
    ]  # fmt: skip
    table = pd.DataFrame([range(len(old) + 1)], columns=["Name", *old])
    renamed = mio.rename_1x_columns(table)
    assert list(renamed.columns) == ["Name", *new]
    assert renamed.to_numpy().tolist() == table.to_numpy().tolist()
    assert mio.rename_1x_columns(renamed) is renamed
    both = pd.DataFrame({"Final_bin": ["a"], "Final_cluster": ["b"]})
    assert mio.rename_1x_columns(both) is both
    message = f"Contig data of metaTOR 1.x: its columns {', '.join(old)} are read as {', '.join(new)}."
    assert logged() == [("INFO", message)]


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
    assert list(data.Final_cluster) == ["bin_1", "bin_2", "ND"]
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
    mio.write_bin_summary(summary, str(tmp_path / "cluster_summary.txt"))
    table = pd.read_csv(tmp_path / "cluster_summary.txt", sep="\t", index_col=0)
    assert list(table.index) == ["bin_1", "bin_2"]
    assert table.loc["bin_1", "Weighted completeness"] == 0.91
    assert table.loc["bin_1", "HiC_abundance"] == 1.2346
    assert ("SG_abundance" in table.columns) == shotgun


@pytest.mark.parametrize("with_tmp_column", [True, False])
def test_write_mge_data(tmp_path, with_tmp_column):
    data = pd.DataFrame({"Name": ["c1"], "Score": [0.123456]})
    if with_tmp_column:
        data["tmp"] = "c1_1"
    mio.write_mge_data(data, str(tmp_path / "mges.tsv"))
    assert (tmp_path / "mges.tsv").read_text() == "Name\tScore\nc1\t0.12\n"


def test_write_contigs(tmp_path):
    """The records of the contigs, in the order given, each as it is in the
    fasta file (description and line breaks kept), as `pyfastx extract`
    writes them. A contig that is not in the file is an error, not a missing
    record."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">c1 first contig\nACGT\nAC\n>c2\nGGGG\n>c3\nTT\n")
    out = tmp_path / "out.fa"
    mio.write_contigs(pyfastx.Fasta(str(fasta)), ["c3", "c1"], str(out))
    assert out.read_text() == ">c3\nTT\n>c1 first contig\nACGT\nAC\n"
    with pytest.raises(KeyError, match="c4"):
        mio.write_contigs(pyfastx.Fasta(str(fasta)), ["c2", "c4"], str(out))


def test_write_contigs_numeric_names(tmp_path):
    """Contigs named by numbers, as pandas reads them from a contig data
    file (integers): written by name, not taken for the index of a record
    (pyfastx reads an integer key so)."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">1\nAAAA\n>2\nCCCC\n>3\nGGGG\n")
    out = tmp_path / "out.fa"
    mio.write_contigs(pyfastx.Fasta(str(fasta)), pd.Series([1, 2]), str(out))
    assert out.read_text() == ">1\nAAAA\n>2\nCCCC\n"
