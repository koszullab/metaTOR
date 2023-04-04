# Test for partition module

import metator.io as mio
import metator.partition as mtp
import networkx as nx
import numpy as np
import os
import pandas as pd
import pytest
import shutil

assembly = "tests_data/assembly.fa"
network_file = "tests_data/outdir/network.txt"
iterations = 5
resolution_parameter = 0.9
spins = 2
threads = 8
LEIDEN_PATH = os.environ["LEIDEN_PATH"]
LOUVAIN_PATH = os.environ["LOUVAIN_PATH"]
overlapping_parameter = 0.6

tmp_dir = "tmp_partition_clustering"
os.makedirs(tmp_dir, exist_ok=True)
partition = mtp.louvain_iterations_cpp(
    network_file, iterations, tmp_dir, LOUVAIN_PATH
)
shutil.rmtree(tmp_dir)

contigs_data = pd.read_csv(
    "tests_data/outdir/contig_data_partition.txt", sep="\t"
)
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
    "Core_bin_ID",
    "Core_bin_contigs",
    "Core_bin_size",
    "Overlapping_bin_ID",
    "Overlapping_bin_contigs",
    "Overlapping_bin_size",
]


def test_algo_partition():
    # Test algo partition choice.
    tmp_dir = "tmp_partition_partition"
    os.makedirs(tmp_dir, exist_ok=True)
    network = nx.read_edgelist(
        network_file, nodetype=int, data=(("weight", float),)
    )
    subnetwork = network.subgraph(np.arange(1, 5))
    for algorithm in ["louvain", "leiden", "error"]:
        try:
            mtp.algo_partition(
                algorithm,
                network_file,
                subnetwork,
                iterations,
                resolution_parameter,
                tmp_dir,
                spins,
            )
        except ValueError:
            assert algorithm == "error"
    shutil.rmtree(tmp_dir)


def test_build_clustering_matrix():
    # Test clustering matrix building.
    hamming_distance = mtp.get_hamming_distance(core_bins_iterations, threads)
    M = mtp.build_clustering_matrix(core_bins_contigs, hamming_distance, 10)
    assert np.sum(M.data) == pytest.approx(14.6, abs=1e-5)
    assert M.shape == (11, 11)
    assert M.nnz == 21


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
    cc_contigs, cc_iterations = mtp.detect_core_bins(
        output_partition, iterations
    )
    assert cc_contigs == core_bins_contigs
    assert (cc_iterations == core_bins_iterations).all().all()


def test_generate_fasta():
    # Test fasta generation.
    tmp_dir = "tmp_partition_fasta"
    os.makedirs(tmp_dir, exist_ok=True)
    mtp.generate_fasta(
        assembly,
        overlapping_bins,
        contigs_data,
        50_000,
        tmp_dir,
        tmp_dir,
        "MetaTOR",
    )
    assert len(os.listdir(tmp_dir)) == 10
    shutil.rmtree(tmp_dir)


def test_get_distances_splitmat():
    # Test hamming distance computation worker.
    x = mtp.get_distances_splitmat(
        core_bins_iterations[0:1], core_bins_iterations
    )
    assert np.sum(x.data) == pytest.approx(1.8, abs=1e-5)
    assert x.shape == (8, 1)
    assert x.nnz == 3


def test_get_hamming_distance():
    # Test hamming distance computation.
    hamming_distance = mtp.get_hamming_distance(core_bins_iterations, threads)
    assert np.sum(hamming_distance.data) == pytest.approx(12, abs=1e-5)
    assert hamming_distance.shape == (8, 8)
    assert hamming_distance.nnz == 22


def test_leiden_iterations_java():
    # Test leiden partition.
    tmp_dir = "tmp_partition_clustering"
    os.makedirs(tmp_dir, exist_ok=True)
    partition = mtp.leiden_iterations_java(
        network_file, iterations, resolution_parameter, tmp_dir, LEIDEN_PATH
    )
    _val = int(partition[1].split(";")[0])
    assert len(partition) == 1058
    assert len(partition[1].split(";")) == iterations
    shutil.rmtree(tmp_dir)


def test_louvain_iterations_cpp():
    # Test louvain partition.
    _val = int(partition[1].split(";")[0])
    assert len(partition) == 1058
    assert len(partition[1].split(";")) == iterations


def test_partition():
    ...


def test_remove_isolates():
    # Test isolate removing from partition.
    partition1 = dict()
    for i in range(1, 1219):
        try:
            partition1[i] = partition[i]
        except KeyError:
            partition1[i] = ";".join(
                map(str, map(int, np.ones(iterations) * i))
            )
    partition2 = mtp.remove_isolates(partition1, network_file)
    assert partition2 == partition


# def test_spinglass_partition():
# Test spinglass partition.
# network = nx.read_edgelist(
#     network_file, nodetype=int, data=(("weight", float),)
# )
# subnetwork = network.subgraph(np.arange(1, 5))
# partition = mtp.spinglass_partition(subnetwork, spins)
# assert partition == {1: "0", 2: "0", 3: "0", 4: "0"}


def test_update_contigs_data():
    # Updating contigs data with partition values.
    tmp_dir = "tmp_partition_contigs"
    os.makedirs(tmp_dir, exist_ok=True)
    contigs_data, _ = mtp.update_contigs_data(
        "tests_data/outdir/contig_data_network.txt",
        core_bins_contigs,
        overlapping_bins,
        tmp_dir,
    )
    assert list(contigs_data.columns) == cols
    assert contigs_data.loc[3, "Overlapping_bin_ID"] == "00002"
    assert contigs_data.loc[0, "Overlapping_bin_contigs"] == 2
    assert contigs_data.loc[1, "Overlapping_bin_size"] == 37801
    shutil.rmtree(tmp_dir)
