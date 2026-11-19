# Test for io module.

import bz2
import collections
import gzip
import os
import shutil
import subprocess as sp
import types
import warnings
import zipfile

import metator.io as mio
import numpy as np
import pandas as pd
import pyfastx
import pypairix
import pysam
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


def test_check_fasta_index_of_the_reference_only(tmp_path):
    """Only the index of the reference itself is its index: another index of
    its directory, that of another assembly, was taken for it."""
    for suffix in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]:
        (tmp_path / f"other.fa.{suffix}").touch()
        (tmp_path / f"assembly.fa2.{suffix}").touch()
    for suffix in ["amb", "ann", "bwt", "pac", "sa"]:
        (tmp_path / f"other.fa.{suffix}").touch()
    fasta = str(tmp_path / "assembly.fa")
    assert mio.check_fasta_index(fasta, mode="bwa") is None
    assert mio.check_fasta_index(fasta, mode="bowtie2") is None
    # A large bowtie2 index.
    for suffix in ["1", "2", "3", "4", "rev.1", "rev.2"]:
        (tmp_path / f"assembly.fa.{suffix}.bt2l").touch()
    assert mio.check_fasta_index(fasta, mode="bowtie2") == fasta
    # An incomplete bwa index.
    for suffix in ["amb", "ann", "bwt", "pac"]:
        (tmp_path / f"assembly.fa.{suffix}").touch()
    assert mio.check_fasta_index(fasta, mode="bwa") is None
    (tmp_path / "assembly.fa.sa").touch()
    assert mio.check_fasta_index(fasta, mode="bwa") == fasta


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


def test_check_pypairix():
    test = mio.check_pypairix()
    assert test


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


def test_read_fasta_names(tmp_path):
    """The names of the sequences of a fasta file, as Biopython names its
    records, which create_contig_data reads: the first word of each header
    line, a header without one included."""
    from Bio import SeqIO

    fasta = tmp_path / "contigs.fa"
    fasta.write_bytes(">c1 a description\nACGT\nACGT\n>  c2\tx\r\nAC\r\n>\nAC\n>NODE_522_é\nAC\n".encode())
    for path in [fasta, "tests_data/assembly.fa", "tests_data/subset.assembly.fa"]:
        assert list(mio.read_fasta_names(str(path))) == [record.id for record in SeqIO.parse(str(path), "fasta")]
    assert list(mio.read_fasta_names(str(fasta))) == ["c1", "c2", "", "NODE_522_é"]


def _pairs_key(row):
    """The order of pairtools sort, and of hicstuff: contig names byte by byte
    (NODE_10 before NODE_2), then positions as numbers."""
    return row[1].encode(), row[3].encode(), int(row[2]), int(row[4])


def _read_pairs(path):
    """Header lines and rows (split on tabs) of a pairs file, compressed or not."""
    with mio.read_compressed(str(path)) as handle:
        lines = handle.read().splitlines()
    return [line for line in lines if line.startswith("#")], [line.split("\t") for line in lines if not line.startswith("#")]


def test_sort_and_index_pairs(tmp_path):
    test_file = tmp_path / "test.pairs"
    shutil.copyfile("tests_data/alignment.pairs", test_file)
    sorted_file = mio.sort_and_index_pairs(str(test_file), 1, True, False)
    assert sorted_file == str(tmp_path / "test_sorted.pairs.gz")
    assert os.path.isfile(f"{sorted_file}.px2")
    assert not test_file.exists()  # remove=True
    input_header, input_body = _read_pairs("tests_data/alignment.pairs")
    header, body = _read_pairs(sorted_file)
    # The header is kept, with the new sort order.
    assert "#sorted: readID" in input_header
    assert header == [("#sorted: chr1-chr2-pos1-pos2" if line.startswith("#sorted") else line) for line in input_header]
    # The same pairs, by contig names in byte order (NODE_1404 before
    # NODE_522), then positions: the order of pairtools sort.
    assert sorted(body) == sorted(input_body)
    keys = [_pairs_key(row) for row in body]
    assert keys == sorted(keys)
    # pypairix finds the pairs between two contigs.
    sizes = dict(line.split()[1:] for line in header if line.startswith("#chromsize"))
    pairs = collections.Counter((row[1], row[3]) for row in input_body if row[1] != row[3])
    pairs_data = pypairix.open(sorted_file)

    def query(chr1, chr2, autoflip):
        return [row[0] for row in pairs_data.query2D(chr1, 0, int(sizes[chr1]), chr2, 0, int(sizes[chr2]), autoflip)]

    # Contigs whose pairs all have the same orientation: asked the other way
    # round, as metaTOR's readers ask, pypairix flips the query.
    first, second = max((key for key in pairs if key[::-1] not in pairs), key=pairs.get)
    found = query(first, second, 0)
    assert len(found) == pairs[first, second] > 1
    with warnings.catch_warnings(record=True) as caught:
        assert query(second, first, 0) == []
    assert "flipped order" in str(caught[0].message)
    assert sorted(query(second, first, 1)) == sorted(found)
    # Contigs with pairs in both orientations: one block for each.
    first, second = max((key for key in pairs if key[::-1] in pairs), key=lambda key: pairs[key] + pairs[key[::-1]])
    assert len(query(first, second, 0)) == pairs[first, second] > 0
    assert len(query(second, first, 0)) == pairs[second, first] > 0


# Contigs whose names sort differently in byte order (c1, c10, c2), the
# order of pairtools sort, and in version order (c1, c2, c10).
SMALL_PAIRS = [
    "r1\tc2\t30\tc10\t5\t+\t-",
    "r2\tc1\t50\tc2\t7\t+\t+",
    "r3\tc2\t10\tc10\t9\t-\t-",
    "r4\tc1\t20\tc2\t3\t+\t-",
    "r5\tc10\t4\tc10\t8\t+\t+",
]
SMALL_PAIRS_ORDER = ["r4", "r2", "r5", "r3", "r1"]
FORMAT_LINE = "## pairs format v1.0"
COLUMNS_LINE = "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2"
SORTED_LINE = "#sorted: chr1-chr2-pos1-pos2"


def _write_small_pairs(path, header=(FORMAT_LINE, "#sorted: readID", COLUMNS_LINE), compression=None):
    text = "".join(f"{line}\n" for line in [*header, *SMALL_PAIRS])
    if compression == "bgzip":
        path.with_suffix("").write_text(text)
        pysam.tabix_compress(str(path.with_suffix("")), str(path))
        path.with_suffix("").unlink()
    else:
        with {None: open, "gzip": gzip.open, "bz2": bz2.open}[compression](path, "wt") as handle:
            handle.write(text)
    return path


@pytest.mark.parametrize(
    "header, expected",
    [
        ([FORMAT_LINE, "#sorted: readID", COLUMNS_LINE], [FORMAT_LINE, SORTED_LINE, COLUMNS_LINE]),
        ([FORMAT_LINE, COLUMNS_LINE, SORTED_LINE], [FORMAT_LINE, COLUMNS_LINE, SORTED_LINE]),
        # No #sorted line: one is added, after the format line.
        ([FORMAT_LINE, COLUMNS_LINE], [FORMAT_LINE, SORTED_LINE, COLUMNS_LINE]),
        ([COLUMNS_LINE], [SORTED_LINE, COLUMNS_LINE]),
        ([], [SORTED_LINE]),
    ],
)
def test_sort_and_index_pairs_header(tmp_path, header, expected):
    pairs = _write_small_pairs(tmp_path / "sample.pairs", header)
    sorted_pairs = mio.sort_and_index_pairs(str(pairs))
    header, body = _read_pairs(sorted_pairs)
    assert header == expected
    assert [row[0] for row in body] == SMALL_PAIRS_ORDER
    assert sorted(body) == sorted(line.split("\t") for line in SMALL_PAIRS)
    assert not os.path.exists(sorted_pairs[: -len(".gz")])
    assert not os.path.exists(f"{sorted_pairs[: -len('.gz')]}.tmp")
    assert pairs.exists()  # remove=False
    assert [row[0] for row in pypairix.open(sorted_pairs).query2D("c2", 0, 100, "c10", 0, 100, 0)] == ["r3", "r1"]


@pytest.mark.parametrize(
    "compression, name, sorted_name",
    [
        (None, "sample.pairs", "sample_sorted.pairs.gz"),
        ("gzip", "sample.pairs.gz", "sample.pairs_sorted.pairs.gz"),
        ("bgzip", "sample.pairs.gz", "sample.pairs_sorted.pairs.gz"),
        ("bz2", "sample.pairs.bz2", "sample.pairs_sorted.pairs.gz"),
    ],
)
def test_sort_and_index_pairs_compressed_input(tmp_path, compression, name, sorted_name):
    """Compressed inputs are sorted alike; the output is named after the
    input without its last extension."""
    pairs = _write_small_pairs(tmp_path / name, compression=compression)
    sorted_pairs = mio.sort_and_index_pairs(str(pairs), threads=2, remove=True, tmp_dir=str(tmp_path))
    assert sorted_pairs == str(tmp_path / sorted_name)
    assert not pairs.exists()
    header, body = _read_pairs(sorted_pairs)
    assert header == [FORMAT_LINE, SORTED_LINE, COLUMNS_LINE]
    assert [row[0] for row in body] == SMALL_PAIRS_ORDER
    assert [row[0] for row in pypairix.open(sorted_pairs).query2D("c2", 0, 100, "c1", 0, 100, 1)] == ["r4", "r2"]


@pytest.mark.parametrize(
    "bgzip",
    [pytest.param(True, marks=pytest.mark.skipif(shutil.which("bgzip") is None, reason="bgzip is not installed")), False],
)
def test_sort_and_index_pairs_force(tmp_path, monkeypatch, bgzip):
    """An existing output is kept unless force; bgzip compresses, or pysam
    without it."""
    if not bgzip:
        real_which = shutil.which
        monkeypatch.setattr(mio.shutil, "which", lambda name: None if name == "bgzip" else real_which(name))
    pairs = _write_small_pairs(tmp_path / "sample.pairs")
    sorted_pairs = mio.sort_and_index_pairs(str(pairs))
    expected = _read_pairs(sorted_pairs)
    # BGZF, the gzip variant pypairix needs: "BC" extra subfield.
    with open(sorted_pairs, "rb") as handle:
        assert handle.read(14)[12:] == b"BC"
    with pytest.raises(ValueError):
        mio.sort_and_index_pairs(str(pairs))
    # The uncompressed file of an interrupted run blocks it too.
    os.rename(sorted_pairs, f"{sorted_pairs}.bak")
    (tmp_path / "sample_sorted.pairs").write_text("unfinished\n")
    with pytest.raises(ValueError):
        mio.sort_and_index_pairs(str(pairs))
    os.rename(f"{sorted_pairs}.bak", sorted_pairs)
    # force overwrites both, and the index; remove deletes the input.
    _write_small_pairs(pairs, header=[FORMAT_LINE])
    assert mio.sort_and_index_pairs(str(pairs), remove=True, force=True) == sorted_pairs
    assert not pairs.exists()
    assert not (tmp_path / "sample_sorted.pairs").exists()
    assert _read_pairs(sorted_pairs) == ([FORMAT_LINE, SORTED_LINE], expected[1])
    assert [row[0] for row in pypairix.open(sorted_pairs).query2D("c1", 0, 100, "c2", 0, 100, 0)] == ["r4", "r2"]


def test_sort_and_index_pairs_leading_zeros(tmp_path):
    """contig_01 and contig_1, one name to version order (sort -V of GNU
    coreutils >= 9.1), are two in byte order: each gets its own block,
    which pypairix indexes."""
    pairs = tmp_path / "sample.pairs"
    pairs.write_text(
        f"{FORMAT_LINE}\nr1\tcontig_01\t1\tcontig_1\t5\t+\t+\nr2\tcontig_1\t5\tcontig_1\t9\t+\t+\n"
        "r3\tcontig_01\t9\tcontig_1\t2\t+\t+\n"
    )
    sorted_pairs = mio.sort_and_index_pairs(str(pairs))
    _, body = _read_pairs(sorted_pairs)
    assert [row[0] for row in body] == ["r1", "r3", "r2"]
    pairs_data = pypairix.open(sorted_pairs)
    for (chr1, chr2), expected in {("contig_01", "contig_1"): ["r1", "r3"], ("contig_1", "contig_1"): ["r2"]}.items():
        assert [row[0] for row in pairs_data.query2D(chr1, 0, 100, chr2, 0, 100, 0)] == expected
    assert pairs.exists()


def test_save_sparse_matrix_requires_coo(tmp_path):
    with pytest.raises(ValueError, match="coo"):
        mio.save_sparse_matrix(sparse.csr_matrix(np.eye(2)), tmp_path / "matrix.txt")


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


def test_rename_1x_columns(logged):
    """The columns of the contig data that metaTOR 1.x named after bins get
    their 2.0 names. A table of 2.0 is returned as it is, and a column of a
    2.0 name wins over the 1.x column of the same meaning, which is left as
    it is."""
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


@pytest.mark.parametrize("name", ["contig_data_partition.txt", "contig_data_final.txt"])
def test_contig_data_of_1x(name, tmp_path):
    """A contig data file of metaTOR 1.x reads as the same file of 2.0, in
    read_contig_data (qc) and import_contig_data_mges (mge)."""
    table = f"tests_data/{name}"
    old = synthetic.as_1x(table, tmp_path / name)
    pd.testing.assert_frame_equal(mio.read_contig_data(old), mio.read_contig_data(table))
    binning, mges = {"Contig_0": "metator_00009_00000"}, ["Contig_57977"]
    for got, expected in zip(
        mio.import_contig_data_mges(old, binning, mges), mio.import_contig_data_mges(table, binning, mges), strict=True
    ):
        if isinstance(expected, pd.DataFrame):
            pd.testing.assert_frame_equal(got, expected)
        else:
            assert got == expected


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
    table = mio.read_bin_summary(str(tmp_path / "cluster_summary.txt"))
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


def _count_query(pairs_data, metagenome, first, second):
    """Pairs between two contigs found by pypairix, and in the simulation."""
    sizes = metagenome.truth.Length
    found = sum(1 for _ in pairs_data.query2D(first, 0, sizes[first], second, 0, sizes[second], 1))
    pairs = metagenome.pairs
    expected = (((pairs.chr1 == first) & (pairs.chr2 == second)) | ((pairs.chr1 == second) & (pairs.chr2 == first))).sum()
    return found, expected


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


def test_sort_and_index_pairs_synthetic(metagenome, metagenome_files, tmp_path):
    pairs = tmp_path / "synthetic.pairs"
    shutil.copyfile(metagenome_files["pairs"], pairs)
    sorted_pairs = mio.sort_and_index_pairs(str(pairs), tmp_dir=str(tmp_path))
    assert sorted_pairs == str(tmp_path / "synthetic_sorted.pairs.gz")
    _, body = _read_pairs(sorted_pairs)
    assert len(body) == len(metagenome.pairs)
    # Sorted by chr1, chr2, pos1, pos2, as pairix needs.
    keys = [_pairs_key(row) for row in body]
    assert keys == sorted(keys)
    first, second = (c.name for c in metagenome.molecule("genome_4").contigs[:2])
    found, expected = _count_query(pypairix.open(sorted_pairs), metagenome, first, second)
    assert found == expected > 0


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
