# Test for partition module

import random

import metator.clusters as mtc
import metator.partition as mtp
import numpy as np
import pandas as pd
import pytest

assembly = "tests_data/assembly.fa"
# The network of a run, as metator.network.map_network gives it.
network = pd.read_csv("tests_data/outdir/network.txt", sep="\t", header=None, names=["contig1", "contig2", "weight"])
iterations = 5
resolution_parameter = 0.9
spins = 2
threads = 8

overlapping_parameter = 0.6

output_partition = {
    1: "0;13;0;10;5",
    2: "0;6;3;10;5",
    3: "7;0;7;0;0",
    4: "7;0;7;0;0",
    5: "1;1;7;10;0",
    6: "2;8;1;13;1",
    7: "8;11;6;4;1",
    8: "6;9;15;4;13",
    9: "10;7;3;2;8",
    10: "2;8;1;13;1",
}
core_bins_contigs = {
    0: [1],
    1: [2],
    2: [3, 4],
    3: [5],
    4: [6, 10],
    5: [7],
    6: [8],
    7: [9],
}
core_bins_iterations = pd.DataFrame(
    [
        [0, 13, 0, 10, 5],
        [0, 6, 3, 10, 5],
        [7, 0, 7, 0, 0],
        [1, 1, 7, 10, 0],
        [2, 8, 1, 13, 1],
        [8, 11, 6, 4, 1],
        [6, 9, 15, 4, 13],
        [10, 7, 3, 2, 8],
    ]
)
overlapping_bins = {
    1: [1, 2],
    2: [3, 4],
    3: [5],
    4: [6, 10],
    5: [7],
    6: [8],
    7: [9],
}
cols = [
    "ID",
    "Name",
    "Size",
    "GC_content",
    "Hit",
    "Shotgun_coverage",
    "Restriction_site",
    "Core_cluster_ID",
    "Core_cluster_contigs",
    "Core_cluster_size",
    "Overlapping_cluster_ID",
    "Overlapping_cluster_contigs",
    "Overlapping_cluster_size",
]


def test_algo_partition():
    # Test algo partition choice.
    for algorithm in ["louvain", "leiden", "error"]:
        try:
            partition = mtp.algo_partition(algorithm, network, iterations, resolution_parameter)
            assert len(partition) == 1058
        except ValueError:
            assert algorithm == "error"


def test_defined_overlapping_bins():
    # Test overlapping bin detection.
    hamming_distance = mtp.get_hamming_distance(core_bins_iterations, threads)
    ob = mtp.defined_overlapping_bins(
        overlapping_parameter,
        hamming_distance,
        core_bins_contigs,
    )
    assert ob == overlapping_bins


def test_detect_core_bins():
    # Test core bin detection.
    cc_contigs, cc_iterations = mtp.detect_core_bins(output_partition, iterations)
    assert cc_contigs == core_bins_contigs
    assert (cc_iterations == core_bins_iterations).all().all()


# Number of the 5 iterations in which two core bins of core_bins_iterations
# share their community: 5 times one minus their normalized Hamming distance.
shared_iterations = np.array(
    [
        [5, 3, 0, 1, 0, 0, 0, 0],
        [3, 5, 0, 1, 0, 0, 0, 1],
        [0, 0, 5, 2, 0, 0, 0, 0],
        [1, 1, 2, 5, 0, 0, 0, 0],
        [0, 0, 0, 0, 5, 1, 0, 0],
        [0, 0, 0, 0, 1, 5, 1, 0],
        [0, 0, 0, 0, 0, 1, 5, 0],
        [0, 1, 0, 0, 0, 0, 0, 5],
    ]
)


def test_get_distances_splitmat():
    # Test hamming distance computation worker.
    x = mtp.get_distances_splitmat(core_bins_iterations[0:1], core_bins_iterations)
    assert np.sum(x.data) == pytest.approx(1.8, abs=1e-5)
    assert x.shape == (8, 1)
    assert x.nnz == 3


def test_get_distances_splitmat_values():
    """One column per core bin of the slice, one row per core bin, with the
    values scikit-learn's 1 - pairwise_distances(metric="hamming") gave."""
    x = mtp.get_distances_splitmat(core_bins_iterations[2:5], core_bins_iterations)
    assert x.format == "csr"
    np.testing.assert_allclose(x.toarray(), shared_iterations[:, 2:5] / 5, rtol=0, atol=1e-12)
    # Exact zeros (no shared iteration) are not stored.
    assert x.nnz == np.count_nonzero(shared_iterations[:, 2:5])


def test_get_hamming_distance():
    # Test hamming distance computation.
    hamming_distance = mtp.get_hamming_distance(core_bins_iterations, threads)
    assert np.sum(hamming_distance.data) == pytest.approx(12, abs=1e-5)
    assert hamming_distance.shape == (8, 8)
    assert hamming_distance.nnz == 22
    np.testing.assert_allclose(hamming_distance.toarray(), shared_iterations / 5, rtol=0, atol=1e-12)


def test_leiden_iterations():
    # Test leiden partition.
    partition = mtp.leiden_iterations(network, iterations, resolution_parameter)
    _val = int(partition[1].split(";")[0])
    assert len(partition) == 1058
    assert len(partition[1].split(";")) == iterations


def test_louvain_iterations():
    # Test louvain partition.
    partition = mtp.louvain_iterations(network, iterations)
    _val = int(partition[1].split(";")[0])
    assert len(partition) == 1058
    assert len(partition[1].split(";")) == iterations


def test_network_graph():
    """The vertices are the contigs of the edges, in the order of their IDs;
    the edges are the rows, in their order."""
    graph, contigs = mtp.network_graph(network)
    assert graph.vcount() == len(contigs) == 1058
    assert graph.ecount() == 13415
    assert graph.is_weighted()
    assert list(contigs) == sorted(set(network.contig1) | set(network.contig2))
    first = graph.es[0]
    assert (contigs[first.source], contigs[first.target]) == (network.contig1[0], network.contig2[0])
    assert graph.es["weight"] == list(network.weight)


def test_partition_table():
    """The contig data table of a partition, as contig_data_partition.txt held
    it: the contig data, left as it is, with the core and the overlapping
    cluster of each contig, its number of contigs and its size."""
    contig_data = pd.read_csv("tests_data/outdir/contig_data_network.txt", sep="\t")
    before = contig_data.copy()
    partition = mtc.Partition(
        core={core_id + 1: contigs for core_id, contigs in core_bins_contigs.items()},
        overlapping=overlapping_bins,
        names={number: mtc.overlapping_name("metator", number) for number in overlapping_bins},
    )
    contigs_data = mtp.partition_table(contig_data, partition)
    # The contig data given is left as it is.
    pd.testing.assert_frame_equal(contig_data, before)
    assert list(contigs_data.columns) == cols
    assert contigs_data.loc[3, "Overlapping_cluster_ID"] == "00002"
    assert contigs_data.loc[0, "Overlapping_cluster_contigs"] == 2
    assert contigs_data.loc[1, "Overlapping_cluster_size"] == 37801
    assert list(contigs_data.Core_cluster_ID[:4]) == ["00001", "00002", "00003", "00003"]
    # A contig in no cluster (no edge in the network): "-".
    assert (contigs_data.loc[10:, ["Core_cluster_ID", "Overlapping_cluster_ID"]] == "-").all().all()


def test_partition():
    """The partition of the network: its core clusters, numbered from 1, and
    its overlapping clusters, each of the contigs of its core clusters in the
    order of their numbers, then of their IDs, named after the prefix; every
    contig of the network in one core cluster and in one overlapping
    cluster."""
    random.seed(0)
    partition = mtp.partition("louvain", network, iterations, overlapping_parameter, resolution_parameter, 1, "sample")
    nodes = sorted(set(network.contig1) | set(network.contig2))
    assert sorted(contig for contigs in partition.core.values() for contig in contigs) == nodes
    assert sorted(contig for contigs in partition.overlapping.values() for contig in contigs) == nodes
    assert sorted(partition.core) == list(range(1, len(partition.core) + 1))
    core = {contig: number for number, contigs in partition.core.items() for contig in contigs}
    for number, contigs in partition.overlapping.items():
        assert contigs == sorted(contigs, key=lambda contig: (core[contig], contig))
        assert partition.names[number] == f"sample_{number:05d}_00000"
    assert all(contigs == sorted(contigs) for contigs in partition.core.values())
