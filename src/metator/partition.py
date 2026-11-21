#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Partition the metaHiC network with the Louvain or Leiden algorithm.

The partition iterates Louvain or Leiden on the network of the contigs of a
map (see metator.network.map_network), finds the core clusters (contigs
always together), then merges core clusters whose Hamming distance is below
the overlap threshold (so as not to split genomes). metator partition MAP
stores the result (metator.clusters.Partition) in the map as clustering
rounds; the validation reads it back (partition_table).
"""

import igraph as ig
import multiprocessing
import numpy as np
import pandas as pd

from functools import partial
import metator.clusters as mtc
from metator.log import logger
from scipy import sparse
from scipy.spatial.distance import cdist


def algo_partition(
    algorithm="louvain",
    network=None,
    iterations=10,
    resolution_parameter=1.0,
):
    """Function to partition the network depednding on the used algorithm.

    Parameters
    ----------
    algorithm : str
        Algorithm to use to partition network.
    network : pandas.DataFrame
        The network: one row per edge, with the ids of its two contigs and its
        weight, normalized or not (columns contig1, contig2 and weight, see
        metator.network.map_network).
    iterations : int
        Number of iterations of the algorithm of Leiden or Louvain.
    resolution_parameter : float
        Resolution parameter for Leiden clustering.

    Returns
    -------
    dict:
        Dictionnary with the id of the contig as key and the list of the
        results of each iterations separated by a semicolon as values."""
    if algorithm == "leiden":
        output_partition = leiden_iterations(network, iterations, resolution_parameter)
    elif algorithm == "louvain":
        output_partition = louvain_iterations(network, iterations)
    else:
        logger.error('algorithm should be either "louvain" or "leiden"')
        raise ValueError
    return output_partition


def defined_overlapping_bins(
    overlap,
    hamming_distance,
    core_bins_contigs,
):
    """This function extract the overlapped bins

    Two core bins overlap when the fraction of iterations that put them
    together reaches the threshold.

    Parameters
    ----------
    overlap : float
        hamming distance threshold use to consider that two bins are
        overlapping.
    hamming_distance : scipy.sparse.csr.csr_matrix
        Matrix with all the previously computed hamming distance between two
        core bins.
    core_bins_contigs : dict
        Dictionnary which has as keys the core bins id and as value the id of
        the contigs of the core bin.

    Returns
    -------
    dict:
        A dictionnary with the id of the overlapping bins as keys and the list
        of id of their contigs as values."""
    # Extract bins which are connected, i.e. bins with an hamming distance
    # superior than the threshold given. The small variation is necessary as
    # python give a float not really equal to the true value
    # (i.e. 0.1 -> 0.09999999999999998)
    connections = hamming_distance >= (overlap - 1e-10)
    overlapping_bins_id = sparse.csgraph.connected_components(connections, directed=False)[1]

    # Create a dictionnary of the overlapped bins (ID from the previous file)
    # with the ID of their contigs as value
    overlapping_bins = {}
    cc_id = 0
    # Iterate on each core bins.
    for oc_id in overlapping_bins_id:
        # Extract contig ID from the core bin.
        core_bin_contigs = core_bins_contigs[cc_id].copy()
        # Add the contig ID on the overlapping bin.
        if oc_id + 1 not in overlapping_bins:
            overlapping_bins[oc_id + 1] = core_bin_contigs
        else:
            overlapping_bins[oc_id + 1] += core_bin_contigs
        cc_id += 1

    logger.info(f"{len(overlapping_bins)} overlapping clusters were found.")

    return overlapping_bins


def detect_core_bins(output_partition, iterations):
    """Detect core bins from the output of the partition algorithm.

    A core bin holds the contigs that every iteration puts in the same
    community.

    Parameters
    ----------
    output_partition : dict
        Dictionnary with the id of the contig as key and the list of the
        results of each iterations separated by a semicolon as values.
    iterations : int
        Number of iterations made previously with the partition algorithm.

    Returns
    -------
    dict:
        Dictionnary which has as keys the core bins id and as value the id of
        the contigs of the core bin.
    pandas.core.frame.DataFrame:
        Table with the id of the core bin and their values for each iterations."""
    # finding duplicate values in the output of louvain or leiden using a
    # flipped dictionary.

    # Create dictionnary for core bins
    core_bins = {}
    core_bins_contigs = {}
    core_bins_rows = []
    core_bin_id = 0
    for key, value in output_partition.items():
        if value not in core_bins:
            # Create an entry in a dictionnary with all the contigs with
            # iterations list as a key.
            core_bins[value] = core_bin_id

            # Create an entry in a dictionnary with all the contigs with core
            # bin id as a key.
            core_bins_contigs[core_bin_id] = [key]
            core_bin_id += 1
            # Add a line to compute the array used to compute the distance
            # between two core bins
            core_bins_rows.append(list(map(int, value.split(";"))))
        # If already an entry created for this bin add a contig in the lists.
        else:
            core_bins_contigs[core_bins[value]].append(key)

    # Transform the array in a dataframe
    core_bins_iterations = pd.DataFrame(np.array(core_bins_rows, dtype=int).reshape(-1, iterations))

    logger.info(f"{len(core_bins)} core clusters were found.")

    return core_bins_contigs, core_bins_iterations


def get_distances_splitmat(bins, core_bins_iterations):
    """This function takes a segment of the full iterative clustering matrix and
    computes, for each index (i.e. contig), the hamming distance to each of the
    other indices.

    The values are one minus the normalized Hamming distance: the fraction of
    the iterations in which both core bins fell in the same community.

    Parameters
    ----------
    bins : pandas.core.frame.DataFrame
        Slice of the table with the id of the core bin and their values for
        each iterations.
    core_bins_iterations : pandas.core.frame.DataFrame
        Table with the id of the core bin and their values for each iterations.

    Returns
    -------
    scipy.sparse.csr.csr_matrix:
        matrix of the distance of the possible pairs from the slice of the
        table and the table itself."""
    x = sparse.csr_matrix(1 - cdist(core_bins_iterations, bins.values, metric="hamming"))
    return x


def get_hamming_distance(core_bins_iterations, threads):
    """Generate matrix of Hamming distances between all pairs of core bins.

    Parameters
    ----------
    core_bins_iterations : pandas.core.frame.DataFrame
        Table with the id of the core bin as index and their values for each
        iterations.
    threads : int
        Number of cores to parallelize computation.

    Returns
    -------
    scipy.sparse.csr.csr_matrix:
        Matrix with all the previously computed hamming distance between two
        core bins."""

    # Compute Hamming distances in the core-bin-level iterative clustering
    # matrix, in parallel
    step = 1000
    steps = np.arange(step, len(core_bins_iterations.index) + step, step)
    split_core_bins = [core_bins_iterations[(k - step) : k] for k in steps]
    pool = multiprocessing.Pool(processes=threads)
    res = pool.map(
        partial(
            get_distances_splitmat,
            core_bins_iterations=core_bins_iterations,
        ),
        split_core_bins,
    )
    res = sparse.hstack(res)
    pool.close()
    return res.tocsr()


def leiden_iterations(network, iterations, resolution_parameter):
    """Partition the network with the igraph implementation of Leiden.

    Parameters
    ----------
    network : pandas.DataFrame
        The network (see network_graph).
    iterations : int
        Number of iterations of the algorithm of Leiden.
    resolution_parameter : float
        Resolution parameter for Leiden clustering.

    Returns
    -------
    dict:
        Dictionnary with the id of the contig as key and the list of the
        results of each iterations separated by a semicolon as values."""
    graph, contigs = network_graph(network)
    memberships = [
        graph.community_leiden(
            objective_function="modularity",
            weights="weight",
            resolution=resolution_parameter,
            n_iterations=4,
        ).membership
        for _ in range(iterations)
    ]
    return merge_iterations(contigs, memberships)


def louvain_iterations(network, iterations):
    """Partition the network with the igraph implementation of Louvain.

    Parameters
    ----------
    network : pandas.DataFrame
        The network (see network_graph).
    iterations : int
        Number of iterations of the algorithm of Louvain.

    Returns
    -------
    dict:
        Dictionnary with the id of the contig as key and the list of the
        results of each iterations separated by a semicolon as values."""
    graph, contigs = network_graph(network)
    memberships = [graph.community_multilevel(weights="weight").membership for _ in range(iterations)]
    return merge_iterations(contigs, memberships)


def merge_iterations(contigs, memberships):
    """Merge the clusters of each iteration in one dictionnary.

    Parameters
    ----------
    contigs : numpy.ndarray
        Contig ids of the network vertices.
    memberships : list of list
        Cluster of each vertex, for each iteration.

    Returns
    -------
    dict:
        Dictionnary with the id of the contig as key and the list of the
        results of each iterations separated by a semicolon as values."""
    memberships = np.array(memberships).T
    return {int(contig): ";".join(map(str, clusters)) for contig, clusters in zip(contigs, memberships, strict=True)}


def partition(
    algorithm,
    network,
    iterations,
    overlapping_parameter,
    resolution_parameter,
    threads,
    prefix,
):
    """Partition the network of the contigs of a map: its core clusters, then its
    overlapping clusters.

    Parameters
    ----------
    algorithm : str
        Algorithm to use to partition the network.
    network : pandas.DataFrame
        The network of the contigs of the map (see
        metator.network.map_network).
    iterations : int
        Number of iterations to use for the partition.
    overlapping_parameter : int
        Hamming distance threshold to use to merge bins (percentage).
    resolution_parameter : float
        Resolution parameter to use if Leiden algorithm is chosen.
    threads : int
        Number of threads to use.
    prefix : str
        Sample prefix to use, which names the overlapping clusters (see
        metator.clusters.overlapping_name).

    Returns
    -------
    metator.clusters.Partition:
        The core and the overlapping clusters, numbered from 1, each with the
        IDs of its contigs."""

    # Perform the iterations of Louvain or Leiden to partition the network.
    logger.info("Start iterations:")
    output_partition = algo_partition(algorithm, network, iterations, resolution_parameter)

    # Detect core bins
    logger.info("Detect core clusters:")
    (
        core_bins_contigs,
        core_bins_iterations,
    ) = detect_core_bins(output_partition, iterations)

    # Compute the Hamming distance between core bins.
    logger.info("Detect overlapping clusters:")
    hamming_distance = get_hamming_distance(
        core_bins_iterations,
        threads,
    )

    # Defined overlapping bins according to the threshold
    overlapping_bins = defined_overlapping_bins(
        overlapping_parameter,
        hamming_distance,
        core_bins_contigs,
    )

    # The core bins numbered from 1, as the contig data tables number them,
    # and the overlapping bins, numbered from 1 already.
    return mtc.Partition(
        core={core_id + 1: contigs for core_id, contigs in core_bins_contigs.items()},
        overlapping=dict(overlapping_bins),
        names={number: mtc.overlapping_name(prefix, number) for number in overlapping_bins},
    )


def network_graph(network):
    """The network as an igraph graph: its vertices are the contigs of its edges,
    in the order of their ids, and its edges are in the order of its rows.

    Parameters
    ----------
    network : pandas.DataFrame
        One row per edge: the ids of its two contigs (columns contig1 and
        contig2) and its weight, normalized or not (column weight), as
        metator.network.map_network gives them.

    Returns
    -------
    igraph.Graph:
        Weighted undirected network.
    numpy.ndarray:
        Contig ids of the vertices."""
    contigs, vertices = np.unique(network[["contig1", "contig2"]].to_numpy().ravel(), return_inverse=True)
    graph = ig.Graph(
        n=len(contigs),
        edges=vertices.reshape(-1, 2).tolist(),
        edge_attrs={"weight": network["weight"].tolist()},
    )
    return graph, contigs


def partition_table(contig_data, partition):
    """The contig data table of a partition: the contig data with the core and the
    overlapping cluster of each contig, as contig_data_partition.txt held it,
    and as the validation reads it.

    Parameters
    ----------
    contig_data : pandas.DataFrame
        The contig data, which is left as it is.
    partition : metator.clusters.Partition
        The core and the overlapping clusters of the partition.

    Returns
    -------
    pandas.core.frame.DataFrame:
        Table with all the information on the contigs included their
        appartenance to the bins."""

    # A copy of the table, from 0 to n - 1: the contig of ID i is row i - 1.
    contigs_data = contig_data.reset_index(drop=True).copy()
    sizes = contigs_data["Size"].to_numpy()
    for kind, clusters in (("Core", partition.core), ("Overlapping", partition.overlapping)):
        ids, counts, lengths = _cluster_columns(clusters, sizes)
        # Object columns, as they mix "-" and numbers.
        contigs_data[f"{kind}_cluster_ID"] = pd.Series(ids, index=contigs_data.index, dtype=object)
        contigs_data[f"{kind}_cluster_contigs"] = pd.Series(counts, index=contigs_data.index, dtype=object)
        contigs_data[f"{kind}_cluster_size"] = pd.Series(lengths, index=contigs_data.index, dtype=object)

    return contigs_data


def _cluster_columns(clusters, sizes):
    """The cluster of each contig of a contig data table (see partition_table):
    its ID (the number of its cluster, with five digits), its number of contigs
    and its size in bp (sizes: the size of each contig, by ID - 1), ``-`` for a
    contig in none. clusters: the IDs of the contigs of each cluster, by
    number."""
    ids, counts, lengths = (np.full(len(sizes), "-", dtype=object) for _ in range(3))
    for number, contigs in clusters.items():
        rows = np.asarray(contigs, dtype=np.int64) - 1
        ids[rows] = f"{number:05d}"
        counts[rows] = len(rows)
        lengths[rows] = sizes[rows].sum()
    return ids, counts, lengths
