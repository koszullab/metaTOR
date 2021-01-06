#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Partition HiC network using Louvain function.
General utility function for generating communities from the network file.

Core functions to partition the network are:
    - defined_overlapping_communities
    - detect_core_communities
    - extract_contigs
    - generate_fasta
    - get_distances_splitmat
    - hamming_distance
    - louvain_partition_cpp
    - louvain_partition_py
    - update_contigs_data
"""


import multiprocessing
import community as community_louvain
import metator.io as mio
import networkx as nx
import numpy as np
import pandas as pd
import subprocess as sp
import sys
from Bio import SeqIO
from functools import partial
from metator.log import logger
from os.path import join
from scipy import sparse
from sklearn import metrics


def defined_overlapping_communities(
    overlap, hamming_distance, core_communities, core_communities_iterations
):
    """This function extract the overlapped communities

    From the hamming distances between the core communities, the function
    identifies the overlapping communities and create a dictionnary with the
    list of the contigs ID for each core community.

    Two core communities are considered overlapping if there have a percentage
    of identity superior or equal to the threshold given.

    Parameters:
    -----------
    overlap : float
        Threshold use to consider that two communities are overlapping.
    hamming_distance : scipy.sparse.csr.csr_matrix
        Matrix with all the previously computed hamming distance between two
        core communities
    core_communities : dict
        Dictionaary which has as keys the values of the iterations from Louvain
        separated by a semicolon and as value the id of the contigs of the core
        community.
    core_communities_iteration : pandas.core.frame.DataFrame
        Table with the id of the core community and their values for each
        iterations.

    Returns:
    --------
    dict:
        A dictionnary with the id of the overlapping communities as keys and
        the list of id of their contigs as values.
    """
    # Extract communities which are connected, i.e. communities with an hamming
    # distance superior than the threshold given.
    connections = hamming_distance >= overlap
    overlapping_communities_id = sparse.csgraph.connected_components(
        connections, directed=False
    )[1]

    # Create a dictionnary of the overlapped communities (ID from the previous
    # file) with the ID of their contigs as value
    overlapping_communities = {}
    cc_id = 0
    # Iterate on each core communities.
    for oc_id in overlapping_communities_id:
        # Extract contig ID from the core community.
        core_community = core_communities_iterations.iloc[cc_id]
        core_community = map(str, list(core_community))
        core_community_contigs = core_communities[";".join(core_community)]
        # Add the contig ID on the overlapping community.
        if oc_id + 1 not in overlapping_communities:
            overlapping_communities[oc_id + 1] = core_community_contigs
        else:
            overlapping_communities[oc_id + 1] += core_community_contigs
        cc_id += 1
    return overlapping_communities


def detect_core_communities(output_louvain, iterations):
    """Detect core commmunities from the output of louvain

    The function search for duplicated values in the output of Louvain algorithm
    in order to find contigs which are always in the same community. The
    communities find with this method are called the core communities.

    Parameters
    ----------
    output_louvain : dict
        Dictionnary with the id of the contig as key and the list of the results
        of each iterations separated by a semicolon as values.
    iterations : int
        Number of iterations made previously with Louvain algorithm.

    Returns
    -------
    dict:
        Dictionaary which has as keys the values of the iterations from Louvain
        separated by a semicolon and as value the id of the contigs of the core
        community.
    pandas.core.frame.DataFrame:
        Table with the id of the core community and their values for each
        iterations.
    """
    # finding duplicate values in the output of louvain using a flipped
    # dictionary.

    # Create dictionnary for core communities
    core_communities = {}
    core_communities_iterations = np.empty((0, iterations), int)
    for key, value in output_louvain.items():
        if value not in core_communities:
            core_communities[value] = [key]
            # Add a line to compute the array used to compute the distance
            # between two core communities
            core_communities_iterations = np.append(
                core_communities_iterations,
                np.array([list(map(int, value.split(";")))]),
                axis=0,
            )
        else:
            core_communities[value].append(key)

    # Transform the array in a dataframe
    core_communities_iterations = pd.DataFrame(core_communities_iterations)
    return core_communities, core_communities_iterations


def extract_contigs(assembly, list_contigs, output_file):
    """Extract the contigs of the list given from a fasta assembly file and
    write them in a new fasta file.

    Parameters:
    -----------
    assembly : str
        Path to the fasta file of the original assembly.
    list_contigs : list of str
        List with the names of the contigs used in the assembly file.
    output_file : str
        Path to the output fasta file.
    """
    # Select contigs of the list.
    records = (
        contigs
        for contigs in SeqIO.parse(assembly, "fasta")
        if contigs.id in list_contigs
    )
    # Write the new fasta file.
    SeqIO.write(records, output_file, "fasta")
    return 0


def generate_fasta(assembly, communities, contigs_data, output_dir):
    """Generate the fasta files of each communities from the assembly.

    Parameters:
    -----------
    assembly : str
        Path to the fasta file of the original assembly.
    communities : dict
        A dictionnary with the id of the overlapping communities as keys and
        the list of id of their contigs as values.
    contigs_data : str
        Path to the file where the contigs data are stored.
    output_dir : str
        Path to the output directory where the fasta of all the community will
        be written.
    """
    contigs_data = pd.read_csv(contigs_data, sep="\t", header=None)
    # For each community create a list of the contigs and extract them from the
    # assembly to create a new fasta file with only the community.
    for community in communities:
        # Extract the list of the contigs from the contigs data file.
        list_contigs_id = communities[community]
        list_contigs = list_contigs_id
        for indice, value in enumerate(list_contigs_id):
            list_contigs[indice] = contigs_data.iloc[value - 1, 1]
        # Define the output file.
        output_file = join(output_dir, "MetaTOR_{0}_0.fa".format(community))
        # Create the fasta file.
        extract_contigs(assembly, list_contigs, output_file)
    return 0


def get_distances_splitmat(comm, core_communities_iterations):
    """This function takes a segment of the full iterative clustering matrix
    and computes, for each index (i.e. contig), the hamming distance to each
    of the other indices.

    Parameters:
    -----------
    comm : pandas.core.frame.DataFrame
        Slice of the table with the id of the core community and their values
        for each iterations.
    core_communities_iterations : pandas.core.frame.DataFrame
        Table with the id of the core community and their values for each
        iterations.

    Returns:
    --------
    scipy.sparse.csr.csr_matrix:
        matrix of the distance of the possible pairs from the slice of the table
        and the table itself.
    """
    x = sparse.csr_matrix(
        1
        - metrics.pairwise_distances(
            core_communities_iterations, comm.values, metric="hamming"
        )
    )
    return x


def hamming_distance(core_communities_iterations, n_iter, threads):
    """Generate matrix of Hamming distances between all pairs of core
    communities

    Parameters
    ----------
    core_communities_iterations : pandas.core.frame.DataFrame
        Table with the id of the core community and their values for each
        iterations.
    core_communities : dict
        Dictionaary which has as keys the values of the iterations from Louvain
        separated by a semicolon and as value the id of the contigs of the core
        community.
    n_iter : int
        Number of iterations of Louvain made previously.
    threads : int
        Number of cores to parallelize computation.

    Returns
    -------
    scipy.sparse.csr.csr_matrix:
        Matrix with all the previously computed hamming distance between two
        core communities
    """

    # Compute Hamming distances in the core-community-level iterative clustering
    # matrix, in parallel
    step = 1000
    steps = np.arange(step, len(core_communities_iterations.index) + step, step)
    split_core_communities = [
        core_communities_iterations[(k - step) : k] for k in steps
    ]
    pool = multiprocessing.Pool(processes=threads)
    res = pool.map(
        partial(
            get_distances_splitmat,
            core_communities_iterations=core_communities_iterations,
        ),
        split_core_communities,
    )
    res = sparse.hstack(res)
    pool.close()
    return res


# def louvain_iterations_cpp(network_file, iterations, output_dir):
#     """Use the cpp original Louvain to partition the network."""
#     # Check if louvain cpp is available in the computer. If it's not available
#     # launch python_louvain instead.
#     if not mio.check_louvain_cpp():
#         logger.warning(
#             "Louvain cpp was not found. Louvain from python is used instead"
#         )
#         return louvain_partition_python

#     # Convert the file in binary file for Louvain partitionning.
#     convert_args = {
#         "net_txt": network_file,
#         "net_bin": network_bin,
#         "net_weights": network_weights,
#     }

#     cmd = ("convert_net -i {net_txt} -o {net_bin} -w {net_weight}").format(
#         **convert.args
#     )
#     sp.Popen(cmd, shell=True)
#     out, err = process.communicate()

#     # convert_net
#     #     -i "$projects"/network/"$library"_network_norm.txt
#     #     -o "$projects"/binning/"$library"_net.bin
#     #     -w "$projects"/binning/"$library"_net.weights

#     # Run the iterations of Louvain
#     for i in range(iterations):
#         logger.info("Iteration in progress: {0}".format(i))
#         # Partiotining with weights using louvain and compute the community
#         # tree.
#         # louvain
#         #     "$projects"/binning/"$library"_net.bin
#         #     -l
#         #     -1
#         #     -w "$projects"/binning/"$library"_net.weights
#         #     > "$projects"/binning/"$library"_net.tree
#         # hierarchy
#         #     "$projects"/binning/"$library"_net.tree
#         #     > "$projects"/binning/"$library"_level_louvain.txt
#         # level=$(tail -1 "$projects"/binning/"$library"_level_louvain.txt
#         #         | awk '{print $2}')
#         # hierarchy
#         #     "$projects"/binning/"$library"_net.tree
#         #     -l "$level"
#         #     > "$projects"/tmp/"$library"_output_louvain_"$iteration".txt
#         # # Save in a temporary file the bin id for each contig at each iterations
#         # # (ordered by the contig id).
#         # cat
#         #     "$projects"/tmp/"$library"_output_louvain_"$iteration".txt
#         #     | awk '{print $2";"}'
#         #     > "$projects"/tmp/"$library"_bin_idx_"$iteration".txt
#     return 0


def louvain_iterations_py(network_file, iterations):
    """Use python-louvain algorithm to partition the network.

    The fonction will make ietrations of the algorithm of Louvain to partition
    the given network. The iterations will allow to select the nodes which will
    always be associated together.

    Parameters:
    -----------
    network_file : str
        Path to the network computed previously. The file is 3 columns table
        separated by a tabulation with the id of the first contigs the id of the
        second one and the weights of the edge normalized or not.
    iterations : int
        Number of iterations of the algorithm of Louvain.

    Returns:
    --------
    dict:
        Dictionnary with the id of the contig as key and the list of the results
        of each iterations separated by a semicolon as values.
    """
    # Convert the file in a networkx graph for Louvain partitionning.
    network = nx.read_edgelist(
        network_file, nodetype=int, data=(("weight", float),)
    )

    # Initiation (first iteration).
    logger.info("Iteration in progress: 1")
    output_louvain = community_louvain.best_partition(network)
    modularity = community_louvain.modularity(output_louvain, network)
    logger.info("Modularity: {0}".format(modularity))

    # Run the iterations of Louvain.
    for iteration in range(2, iterations + 1):
        logger.info("Iteration in progress: {0}".format(iteration))
        partition = community_louvain.best_partition(network)
        modularity = community_louvain.modularity(partition, network)
        logger.info("Modularity: {0}".format(modularity))
        for j in partition:
            output_louvain[j] = "{0};{1}".format(
                output_louvain[j], partition[j]
            )
    return output_louvain


def update_contigs_data(contid_data_file, partition):
    """Add community information in the contigs data file."""
    return 0
