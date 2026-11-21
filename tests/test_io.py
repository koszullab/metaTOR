# Test for io module.

import bz2
import gzip
import os
import zipfile

import metator.io as mio
import pandas as pd
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


def test_mge_contigs_that_the_contig_data_lacks(tmp_path, logged):
    """The contig data tables list the contigs of the map of the run only:
    the MGE contigs that a map made with a contig filter lacks are left out,
    with a warning that says how many and names the first."""
    pd.DataFrame({"ID": [1, 2], "Name": ["c1", "c2"], "Size": [10, 20]}).to_csv(tmp_path / "contigs.txt", sep="\t", index=False)
    mges = ["c2", *(f"short_{i}" for i in range(7)), "short_0"]
    data, mges_ids = mio.import_contig_data_mges(str(tmp_path / "contigs.txt"), {}, mges)
    assert mges_ids == [1] and list(data.MGE) == [False, True]
    assert logged() == [
        (
            "WARNING",
            f"7 MGE contig(s) are not in the contig data {tmp_path / 'contigs.txt'}, which lists the contigs of the map "
            "of the run only (a map made with a contig filter lacks some): they are left out of the MGE clusters: "
            "'short_0', 'short_1', 'short_2', 'short_3', 'short_4', ... (+2 more).",
        )
    ]


@pytest.mark.parametrize("shotgun", [True, False])
def test_bin_summary_round_trip(tmp_path, shotgun):
    summary = {
        "bin_1": {"Length": 1000, "Weighted completeness": 0.912345, "HiC_abundance": 1.234567},
        "bin_2": {"Length": 2000, "Weighted completeness": 0.5, "HiC_abundance": 0.1},
    }
    if shotgun:
        for i, values in enumerate(summary.values()):
            values["SG_abundance"] = 10.0 / (i + 1)
    layout = mio.write_bin_summary(summary, str(tmp_path / "cluster_summary.txt"))
    table = pd.read_csv(tmp_path / "cluster_summary.txt", sep="\t", index_col=0)
    assert list(table.index) == ["bin_1", "bin_2"]
    assert table.loc["bin_1", "Weighted completeness"] == 0.91
    assert table.loc["bin_1", "HiC_abundance"] == 1.2346
    assert ("SG_abundance" in table.columns) == shotgun
    assert layout.rows == ["bin_1", "bin_2"] and layout.columns == list(table.columns)


def test_bin_summary_layout_writes_it_again(tmp_path):
    """The layout of the summary (its rows, columns and formats) writes it
    again, byte for byte, from its values as they were (or as the genomes
    table of the map gives them back: numbers, and text): two decimals for
    the decimals, four for the abundances, the whole numbers as they are,
    the text and the booleans as written, and an empty field for a missing
    value."""
    summary = {
        "b_2": {"Length": 1500.0, "GC-content": 40.123, "recursive": True, "step": 1, "parent": "b_1",
                "HiC_abundance": 1.234567, "SG_abundance": 2.5},
        "b_1": {"Length": 3000.0, "GC-content": float("nan"), "recursive": False, "step": 0, "parent": None,
                "HiC_abundance": 0.1, "SG_abundance": 0.0},
    }  # fmt: skip
    path = tmp_path / "cluster_summary.txt"
    layout = mio.write_bin_summary(summary, str(path))
    assert layout.rows == ["b_2", "b_1"]
    assert layout.formats == {
        "Length": "%.2f", "GC-content": "%.2f", "recursive": "%s", "step": "%d", "parent": "%s",
        "HiC_abundance": "%.4f", "SG_abundance": "%.4f",
    }  # fmt: skip
    lines = ["\t".join(["", *layout.columns])]
    for row in layout.rows:
        values = [summary[row][column] for column in layout.columns]
        lines.append("\t".join([row, *("" if pd.isna(value) else layout.formats[column] % value
                                        for column, value in zip(layout.columns, values, strict=True))]))  # fmt: skip
    assert path.read_text() == "\n".join(lines) + "\n"


@pytest.mark.parametrize("with_tmp_column", [True, False])
def test_write_mge_data(tmp_path, with_tmp_column):
    data = pd.DataFrame({"Name": ["c1"], "Score": [0.123456]})
    if with_tmp_column:
        data["tmp"] = "c1_1"
    mio.write_mge_data(data, str(tmp_path / "mges.tsv"))
    assert (tmp_path / "mges.tsv").read_text() == "Name\tScore\nc1\t0.12\n"
