# Test for network module

import collections
import gzip
import metator.io as mio
import numpy as np
import metator.network as mtn
import pandas as pd
import pytest
import os
import shutil
from os.path import join
from pathlib import Path


assembly = "tests_data/assembly.fa"
depth_file = "tests_data/depth.txt"
alignment_file = "tests_data/alignment.pairs"
contigs_data = {
    "NODE_1": {
        "id": 1,
        "Name": "NODE_1",
        "length": 642311,
        "GC": 38.6876450815882,
        "hit": 30795,
        "coverage": 41.1565,
        "RS": 0,
    },
    "NODE_6": {
        "id": 6,
        "Name": "NODE_6",
        "length": 505338,
        "GC": 53.0124391991103,
        "hit": 5533054,
        "coverage": 1330.17,
        "RS": 2505,
    },
    "NODE_8": {
        "id": 8,
        "Name": "NODE_8",
        "length": 416773,
        "GC": 26.977035460550468,
        "hit": 206333,
        "coverage": 97.6355,
        "RS": 3496,
    },
    "NODE_10": {
        "id": 10,
        "Name": "NODE_10",
        "length": 396325,
        "GC": 36.0124391991103,
        "hit": 306312,
        "coverage": 0,
        "RS": 2578,
    },
}


def test_alignment_to_contacts():
    # Test alignment reader.
    tmp_dir = "tmp_network"
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_file_net = "network.txt"
    tmp_file_contig = "contig_data_network.txt"
    contig_data, hit_data = mtn.create_contig_data(assembly, nb_alignment=2)
    mtn.alignment_to_contacts(
        [alignment_file, alignment_file],
        contig_data,
        0,
        hit_data,
        tmp_dir,
        tmp_file_net,
        tmp_file_contig,
        tmp_dir,
        8,
        "empirical_hit",
        True,
    )
    shutil.rmtree(tmp_dir)


def _expected_network(pairs_file, contig_data, edge=0):
    """Inter-contig contacts per contig pair, counted straight from a pairs
    file: not a pair with a read past the end of its contig, and with edge
    0, every other pair (no edge region)."""
    expected = collections.Counter()
    with open(pairs_file) as pairs:
        for line in pairs:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            c1, pos1, c2, pos2 = fields[1], int(fields[2]), fields[3], int(fields[4])
            len1, len2 = contig_data[c1]["length"], contig_data[c2]["length"]
            id1, id2 = contig_data[c1]["id"], contig_data[c2]["id"]
            if pos1 > len1 or pos2 > len2:
                continue
            inside = not edge or (pos1 > edge and pos2 > edge and len1 - pos1 > edge and len2 - pos2 > edge)
            if id1 != id2 and inside:
                expected[(min(id1, id2), max(id1, id2))] += 1
    return expected


def _read_network(network_file):
    return pd.read_csv(network_file, sep="\t", header=None, names=["id1", "id2", "weight"])


def test_compute_network():
    # Test computing network.
    tmp_dir = "tmp_network"
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_file_pre = join(tmp_dir, "prenetwork.txt")
    tmp_file_sor = join(tmp_dir, "prenetwork_sorted.txt")
    tmp_file_net = join(tmp_dir, "network.txt")
    contig_data, hit_data = mtn.create_contig_data(assembly, nb_alignment=1)
    contig_data, _out_files_list = mtn.precompute_network(
        [alignment_file],
        contig_data,
        0,
        hit_data,
        tmp_file_pre,
        tmp_dir,
    )
    expected = _expected_network(alignment_file, contig_data)
    # Case without normalization: one line per contig pair, with its count.
    mtn.compute_network(
        tmp_file_pre,
        tmp_file_net,
        contig_data,
        tmp_dir,
        tmp_file_sor,
        1,
        "None",
    )
    alpha = pd.read_csv(tmp_file_sor, sep="\t", header=None).iloc[0, 1]
    assert alpha == "NODE_1404"
    network = _read_network(tmp_file_net)
    assert not network.duplicated(["id1", "id2"]).any()
    assert {(a, b): w for a, b, w in network.itertuples(index=False)} == expected
    # Case with normalization, reusing the same file names: counts must not
    # accumulate across calls.
    mtn.compute_network(
        tmp_file_pre,
        tmp_file_net,
        contig_data,
        tmp_dir,
        tmp_file_sor,
        8,
        "length",
    )
    network = _read_network(tmp_file_net)
    assert len(network) == len(expected)
    by_id = {data["id"]: data for data in contig_data.values()}
    for a, b, weight in network.itertuples(index=False):
        factor = np.sqrt(by_id[a]["hit"] / by_id[a]["length"] * by_id[b]["hit"] / by_id[b]["length"])
        assert weight == pytest.approx(expected[(a, b)] / factor, rel=1e-6)
    shutil.rmtree(tmp_dir)


def test_alignment_to_contacts_merges_libraries(tmp_path):
    """With two libraries, network.txt sums them: one line per contig pair."""
    contig_data, hit_data = mtn.create_contig_data(assembly, nb_alignment=2)
    network_file, _ = mtn.alignment_to_contacts(
        [alignment_file, alignment_file],
        contig_data,
        0,
        hit_data,
        str(tmp_path),
        "network.txt",
        "contig_data_network.txt",
        str(tmp_path),
        1,
        "None",
        False,
    )
    merged = _read_network(network_file).set_index(["id1", "id2"])["weight"]
    single = _read_network(tmp_path / "network_0.txt").set_index(["id1", "id2"])["weight"]
    assert merged.index.is_unique
    pd.testing.assert_series_equal(merged.sort_index(), 2 * single.sort_index())


def test_create_contig_data():
    # Test contig data builder.
    # Case without depth file.
    contig_data, hit_data = mtn.create_contig_data(assembly, nb_alignment=3)
    assert len(contig_data) == 1219
    assert contig_data["NODE_522"] == {
        "id": 1,
        "length": 22786,
        "GC": 62.0775915035548,
        "hit": 0,
        "coverage": "-",
        "RS": "-",
    }
    assert hit_data["NODE_114957"] == {"id": 643, "hit": [0, 0, 0]}
    # Case with depth file.
    contig_data, hit_data = mtn.create_contig_data(assembly, 1, depth_file, "DpnII")
    assert len(contig_data) == 1219
    assert contig_data["NODE_522"] == {
        "id": 1,
        "length": 22786,
        "GC": 62.0775915035548,
        "hit": 0,
        "coverage": 4.76595,
        "RS": 162,
    }
    assert hit_data is None


def test_normalize_pair():
    # Test for edges normalisation methods.
    length = mtn.normalize_pair(contigs_data, ("NODE_1", "NODE_6"), 52, "length")
    abundance = mtn.normalize_pair(contigs_data, ("NODE_1", "NODE_8"), 2004, "abundance")
    abundance_null = mtn.normalize_pair(contigs_data, ("NODE_1", "NODE_10"), 275, "abundance")
    RS = mtn.normalize_pair(contigs_data, ("NODE_6", "NODE_8"), 800, "RS")
    EH = mtn.normalize_pair(contigs_data, ("NODE_8", "NODE_1"), 2022, "empirical_hit")
    TH = mtn.normalize_pair(contigs_data, ("NODE_6", "NODE_8"), 3333, "theoritical_hit")
    TH_null = mtn.normalize_pair(contigs_data, ("NODE_1", "NODE_8"), 2004, "theoritical_hit")
    assert length == pytest.approx(71.77, abs=1e-2)
    assert abundance == pytest.approx(31.61, abs=1e-2)
    assert abundance_null == 0
    assert RS == pytest.approx(221.57, abs=1e-2)
    assert EH == pytest.approx(2.54e-2, abs=1e-4)
    assert TH == pytest.approx(2.51e-05, abs=1e-7)
    assert TH_null == 0


def test_precompute_network():
    # Test precompute network.
    tmp_dir = "tmp_network"
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_file = "tmp_precompute_network.tsv"
    # Case of single file with edges.
    contig_data, hit_data = mtn.create_contig_data(assembly, nb_alignment=1)
    contig_data, out_files_list = mtn.precompute_network(
        [alignment_file],
        contig_data,
        2500,
        hit_data,
        tmp_file,
        tmp_dir,
    )
    alpha = pd.read_csv(out_files_list[0], sep="\t", header=None).iloc[105, 0]
    assert contig_data["NODE_522"]["hit"] == 1288
    assert out_files_list == ["tmp_network/prenetwork0.txt"]
    assert alpha == "NODE_2210"
    # Case of multiple files with self contacts.
    contig_data, hit_data = mtn.create_contig_data(assembly, nb_alignment=2)
    contig_data, out_files_list = mtn.precompute_network(
        [alignment_file, alignment_file],
        contig_data,
        0,
        hit_data,
        tmp_file,
        tmp_dir,
        True,
    )
    beta = pd.read_csv(out_files_list[0], sep="\t", header=None).iloc[105, 0]
    assert contig_data["NODE_522"]["hit"] == 2576
    assert out_files_list == [
        "tmp_network/prenetwork0.txt",
        "tmp_network/prenetwork1.txt",
    ]
    assert beta == "NODE_11843"
    os.remove(tmp_file)
    shutil.rmtree(tmp_dir)


@pytest.mark.parametrize("edge", [0, 10])
def test_edge_region(edge, tmp_path):
    """At edge 0, there is no edge region: a pair with a read at position 0
    or on the last base of its contig is a contact, as in a .metacool map
    of the pairs. A pair with a read past the end of its contig is none, at
    any edge, and its hits do not count either (metacooler leaves it out of
    a map). Above 0, the reads within edge bp of a contig end make no
    edge, but count their hits."""
    assembly_file = tmp_path / "assembly.fa"
    assembly_file.write_text(">c1\n" + "A" * 100 + "\n>c2\n" + "C" * 200 + "\n")
    pairs = tmp_path / "test.pairs"
    pairs.write_text(
        "## pairs format v1.0\n"
        "r1\tc1\t50\tc2\t100\t+\t-\n"  # inside
        "r2\tc1\t0\tc2\t100\t+\t-\n"  # position 0
        "r3\tc1\t100\tc2\t200\t+\t-\n"  # last bases
        "r4\tc1\t5\tc2\t195\t+\t-\n"  # within 10 bp of an end
        "r5\tc1\t101\tc2\t50\t+\t-\n"  # past the end of c1
    )
    contig_data, hit_data = mtn.create_contig_data(str(assembly_file), 1)
    mtn.alignment_to_contacts([str(pairs)], contig_data, edge, hit_data, str(tmp_path), "network.txt", "contigs.txt",
                              str(tmp_path), 1, "None", False)  # fmt: skip
    assert (tmp_path / "network.txt").read_text() == f"1\t2\t{4 if edge == 0 else 1}\n"
    assert contig_data["c1"]["hit"] == contig_data["c2"]["hit"] == 4


@pytest.mark.parametrize("edge", [0, 10])
def test_edge_region_with_self_contacts(edge, tmp_path):
    """With self contacts, the pairs within a contig are kept, and --edge
    leaves out the contacts between contigs near a contig end, as without:
    it was ignored with --self-contacts."""
    assembly_file = tmp_path / "assembly.fa"
    assembly_file.write_text(">c1\n" + "A" * 100 + "\n>c2\n" + "C" * 200 + "\n")
    pairs = tmp_path / "test.pairs"
    pairs.write_text(
        "## pairs format v1.0\n"
        "r1\tc1\t50\tc2\t100\t+\t-\n"  # inside
        "r4\tc1\t5\tc2\t195\t+\t-\n"  # within 10 bp of an end
        "r6\tc1\t5\tc1\t60\t+\t-\n"  # within c1, near its start
    )
    contig_data, hit_data = mtn.create_contig_data(str(assembly_file), 1)
    mtn.alignment_to_contacts([str(pairs)], contig_data, edge, hit_data, str(tmp_path), "network.txt", "contigs.txt",
                              str(tmp_path), 1, "None", True)  # fmt: skip
    assert (tmp_path / "network.txt").read_text() == f"1\t1\t1\n1\t2\t{2 if edge == 0 else 1}\n"


def test_create_contig_data_reads_the_depth_file_by_name(tmp_path):
    """The depth file is read by contig name: in another order, it gives the
    same lengths and coverages. It was read line by line, in the order of
    the assembly, and gave each contig the coverage of another."""
    header, *lines = Path(depth_file).read_text().splitlines(keepends=True)
    reversed_depth = tmp_path / "depth.txt"
    reversed_depth.write_text(header + "".join(reversed(lines)))
    expected, _ = mtn.create_contig_data(assembly, 1, depth_file, "DpnII")
    found, _ = mtn.create_contig_data(assembly, 1, str(reversed_depth), "DpnII")
    assert found == expected
    assert found["NODE_522"]["coverage"] == 4.76595


@pytest.mark.parametrize(
    "case, message",
    [
        ("missing", "1 contigs of the assembly are not in the depth file"),
        ("other length", "1 contigs of the assembly have another length in the depth file"),
    ],
)
def test_create_contig_data_refuses_the_depth_file_of_another_assembly(case, message, tmp_path, caplog):
    header, first, *lines = Path(depth_file).read_text().splitlines(keepends=True)
    if case == "other length":
        name, length, *rest = first.split("\t")
        lines.append("\t".join([name, str(int(length) + 1), *rest]))
    other = tmp_path / "depth.txt"
    other.write_text(header + "".join(lines))
    with caplog.at_level("INFO", logger="metator_logger"), pytest.raises(ValueError, match="is not that of the assembly"):
        mtn.create_contig_data(assembly, 1, str(other))
    assert f"{message} {other}, such as NODE_522" in caplog.text


def test_create_contig_data_from_a_compressed_assembly(tmp_path):
    """A gzipped assembly gives the same contig data as the uncompressed one."""
    compressed = tmp_path / "assembly.fa.gz"
    with open(assembly, "rb") as fasta, gzip.open(compressed, "wb") as gz:
        shutil.copyfileobj(fasta, gz)
    assert mtn.create_contig_data(str(compressed), 2, depth_file, "DpnII") == mtn.create_contig_data(
        assembly, 2, depth_file, "DpnII"
    )


def test_write_contig_data():
    # Test contig data writer.
    tmp_file = "tmp_contig_data_test.tsv"
    mtn.write_contig_data(contigs_data, tmp_file)
    data = mio.read_contig_data(tmp_file)
    assert data.loc["NODE_1", "Size"] == 642311
    assert data.loc["NODE_8", "Hit"] == 206333
    assert data.loc["NODE_6", "Restriction_site"] == 2505
    os.remove(tmp_file)


def test_write_hit_data():
    # Test for hit data writer.
    tmp_file = "tmp_hit_data_test.tsv"
    hit_data = {
        "contig_1": {"id": 1, "hit": [50, 25, 88]},
        "contig_2": {"id": 2, "hit": [2, 4, 18]},
    }
    mtn.write_hit_data(hit_data, tmp_file)
    data = pd.read_csv(tmp_file, sep="\t", header=None, index_col=0)
    assert data.loc[2, 3] == 4
    assert data.loc[1, 1] == "contig_1"
    os.remove(tmp_file)
