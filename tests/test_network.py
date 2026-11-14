# Test for network module

import collections
import metator.io as mio
import numpy as np
import metator.network as mtn
import pandas as pd
import pytest
import os
import shutil
from os.path import join


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
    """Inter-contig contacts per contig pair, counted straight from a pairs file."""
    expected = collections.Counter()
    with open(pairs_file) as pairs:
        for line in pairs:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            c1, pos1, c2, pos2 = fields[1], int(fields[2]), fields[3], int(fields[4])
            len1, len2 = contig_data[c1]["length"], contig_data[c2]["length"]
            id1, id2 = contig_data[c1]["id"], contig_data[c2]["id"]
            inside = pos1 > edge and pos2 > edge and len1 - pos1 > edge and len2 - pos2 > edge
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
