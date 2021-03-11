#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Partition HiC network using Louvain function.

General utility function for generating bins from the network file using Louvain
algorithm. If an overlapping threshold is given, it will use the Hamming 
distance to group together bins relatively closed (to avaoid to split genomes in
two or more bins).

Core functions to partition the network are:
    - defined_overlapping_bins
    - detect_core_bins
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


def defined_overlapping_bins(
    overlap, hamming_distance, core_bins, core_bins_iterations
):
    """This function extract the overlapped bins

    From the hamming distances between the core bins, the function identifies
    the overlapping bins and create a dictionnary with the list of the contigs
    ID for each core bin.

    Two core bins are considered overlapping if there have a percentage of
    identity superior or equal to the threshold given.

    Parameters:
    -----------
    overlap : float
        Threshold use to consider that two bins are overlapping.
    hamming_distance : scipy.sparse.csr.csr_matrix
        Matrix with all the previously computed hamming distance between two
        core bins
    core_bins : dict
        Dictionaary which has as keys the values of the iterations from Louvain
        separated by a semicolon and as value the id of the contigs of the core
        bin.
    core_bins_iteration : pandas.core.frame.DataFrame
        Table with the id of the core bin and their values for each iterations.

    Returns:
    --------
    dict:
        A dictionnary with the id of the overlapping bins as keys and the list
        of id of their contigs as values.
    """
    # Extract bins which are connected, i.e. bins with an hamming distance
    # superior than the threshold given.
    connections = hamming_distance >= overlap
    overlapping_bins_id = sparse.csgraph.connected_components(
        connections, directed=False
    )[1]

    # Create a dictionnary of the overlapped bins (ID from the previous file)
    # with the ID of their contigs as value
    overlapping_bins = {}
    cc_id = 0
    # Iterate on each core bins.
    for oc_id in overlapping_bins_id:
        # Extract contig ID from the core bin.
        core_bin = core_bins_iterations.iloc[cc_id]
        core_bin = map(str, list(core_bin))
        core_bin_contigs = core_bins[";".join(core_bin)].copy()
        # Add the contig ID on the overlapping bin.
        if oc_id + 1 not in overlapping_bins:
            overlapping_bins[oc_id + 1] = core_bin_contigs
        else:
            overlapping_bins[oc_id + 1] += core_bin_contigs
        cc_id += 1

    logger.info(
        "{0} overlapping bins were found.".format(len(overlapping_bins))
    )

    return overlapping_bins


def detect_core_bins(output_louvain, iterations):
    """Detect core bins from the output of louvain

    The function search for duplicated values in the output of Louvain algorithm
    in order to find contigs which are always in the same bin. The bins find
    with this method are called the core bins.

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
        Dictionnary which has as keys the values of the iterations from Louvain
        separated by a semicolon and as value the id of the contigs of the core
        bin.
    pandas.core.frame.DataFrame:
        Table with the id of the core bin and their values for each iterations.
    """
    # finding duplicate values in the output of louvain using a flipped
    # dictionary.

    # Create dictionnary for core bins
    core_bins = {}
    core_bins_iterations = np.empty((0, iterations), int)
    for key, value in output_louvain.items():
        if value not in core_bins:
            core_bins[value] = [key]
            # Add a line to compute the array used to compute the distance
            # between two core bins
            core_bins_iterations = np.append(
                core_bins_iterations,
                np.array([list(map(int, value.split(";")))]),
                axis=0,
            )
        else:
            core_bins[value].append(key)

    # Transform the array in a dataframe
    core_bins_iterations = pd.DataFrame(core_bins_iterations)

    logger.info("{0} core bins were found.\n".format(len(core_bins)))

    return core_bins, core_bins_iterations


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


def generate_fasta(assembly, bins, contigs_data, size, output_dir):
    """Generate the fasta files of each bins from the assembly.

    Parameters:
    -----------
    assembly : str
        Path to the fasta file of the original assembly.
    bins : dict
        A dictionnary with the id of the overlapping bins as keys and the list
        of id of their contigs as values.
    contigs_data : pandas.core.frame.DataFrame
        Table with all the information on the contigs included their
        appartenance to the bins.
    size :  int
        Thrshold size chosen to write the bins.
    output_dir : str
        Path to the output directory where the fasta of all the bin will be
        written.
    """

    nb_bins = 0
    length_bins = 0
    # For each bin create a list of the contigs and extract them from the
    # assembly to create a new fasta file with only the bin.
    for bin in bins:
        # Extract the list of the contigs from the contigs data file.
        list_contigs_id = bins[bin]
        list_contigs = list_contigs_id
        # Test if the bin is bigger than the size threshold given.
        length_bin = contigs_data.iloc[list_contigs[0] - 1, 11]
        if length_bin >= size:
            nb_bins += 1
            length_bins += length_bin
            for indice, value in enumerate(list_contigs_id):
                list_contigs[indice] = contigs_data.iloc[value - 1, 1]
            # Define the output file.
            output_file = join(output_dir, "MetaTOR_{0}_0.fa".format(bin))
            # Create the fasta file.
            extract_contigs(assembly, list_contigs, output_file)
    logger.info("{0} bins have been extracted".format(nb_bins))
    logger.info(
        "Total size of the extracted bins: {0}Mb".format(
            round(length_bins / 10 ** 6, 3)
        )
    )
    return 0


def get_distances_splitmat(bins, core_bins_iterations):
    """This function takes a segment of the full iterative clustering matrix and
    computes, for each index (i.e. contig), the hamming distance to each of the
    other indices.

    Parameters:
    -----------
    bins : pandas.core.frame.DataFrame
        Slice of the table with the id of the core bin and their values for each
        iterations.
    core_bins_iterations : pandas.core.frame.DataFrame
        Table with the id of the core bin and their values for each iterations.

    Returns:
    --------
    scipy.sparse.csr.csr_matrix:
        matrix of the distance of the possible pairs from the slice of the table
        and the table itself.
    """
    x = sparse.csr_matrix(
        1
        - metrics.pairwise_distances(
            core_bins_iterations, bins.values, metric="hamming"
        )
    )
    return x


def hamming_distance(core_bins_iterations, n_iter, threads):
    """Generate matrix of Hamming distances between all pairs of core bins

    Parameters
    ----------
    core_bins_iterations : pandas.core.frame.DataFrame
        Table with the id of the core bin and their values for each iterations.
    core_bins : dict
        Dictionaary which has as keys the values of the iterations from Louvain
        separated by a semicolon and as value the id of the contigs of the core
        bin.
    n_iter : int
        Number of iterations of Louvain made previously.
    threads : int
        Number of cores to parallelize computation.

    Returns
    -------
    scipy.sparse.csr.csr_matrix:
        Matrix with all the previously computed hamming distance between two
        core bins
    """

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
    return res


# TODO
def louvain_iterations_cpp(network_file, iterations, tmp_dir, louvain_path):
    """Use the cpp original Louvain to partition the network.

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

    # Check if louvain cpp is available in the computer. If it's not available
    # launch python_louvain instead.
    if not mio.check_louvain_cpp(louvain_path):
        logger.warning(
            "Louvain cpp was not found. Louvain from python is used instead"
        )
        return louvain_iterations_py(network_file, iterations)

    # Defined temporary files and args for louvain fonction calling and path to
    # the variables to call.
    network_bin = join(tmp_dir, "net_bin")
    network_weight = join(tmp_dir, "net_weight")
    network_tree = join(tmp_dir, "net_tree")
    network_labels = join(tmp_dir, "labels.txt")
    level_louvain = join(tmp_dir, "level.txt")
    output = join(tmp_dir, "output_louvain_")
    louvain = join(louvain_path, "louvain")
    convert_net = join(louvain_path, "convert_net")
    hierarchy = join(louvain_path, "hierarchy")
    output_louvain = dict()

    # Create dictionnary of all arguments
    louvain_args = {
        "net_txt": network_file,
        "net_bin": network_bin,
        "net_weight": network_weight,
        "net_tree": network_tree,
        "net_labels": network_labels,
        "level_file": level_louvain,
        "output": output,
        "level": 0,
        "iteration": 0,
        "convert_net": convert_net,
        "louvain": louvain,
        "hierarchy": hierarchy,
    }

    # Convert the file in binary file for Louvain partitionning.
    cmd = ("{convert_net} -i {net_txt} -o {net_bin} -r {net_labels} -w {net_weight}").format(
        **louvain_args
    )
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()

    # Create a dictionary of Louvain labels and original contig id.
    labels = dict()
    with open(louvain_args["net_labels"]) as label_file:
        for label in label_file:
            label = label.split()
            labels[label[1]] = int(label[0])

    # Run the iterations of Louvain
    for i in range(iterations):
        logger.info("Iteration in progress: {0}".format(i))

        louvain_args["iteration"] = i

        # Partiotining with weights using louvain and compute the bin tree.
        cmd = ("{louvain} {net_bin} -l -1 -w {net_weight} > {net_tree}").format(
            **louvain_args
        )
        process = sp.Popen(cmd, shell=True)
        out, err = process.communicate()

        cmd = ("{hierarchy} {net_tree} > {level_file}").format(**louvain_args)
        process = sp.Popen(cmd, shell=True)
        out, err = process.communicate()

        level_file = open(level_louvain, "r")
        louvain_args["level"] = level_file.readlines()[-1][6]
        level_file.close()

        cmd = (
            "{hierarchy} {net_tree} -l {level} > {output}{iteration}.txt"
        ).format(**louvain_args)
        process = sp.Popen(cmd, shell=True)
        out, err = process.communicate()

        # Save the results in a dictionnary
        if iterations == 1:
            with open(output + str(i) + ".txt", "r") as out:
                for line in out:
                    result = line.split(" ")
                    output_louvain[labels[result[0]]] = result[1][:-1]
        elif i == 0:
            with open(output + str(i) + ".txt", "r") as out:
                for line in out:
                    result = line.split(" ")
                    output_louvain[labels[result[0]]] = result[1][:-1] + ";"
        elif i == iterations - 1:
            with open(output + str(i) + ".txt", "r") as out:
                for line in out:
                    result = line.split(" ")
                    output_louvain[labels[result[0]]] += result[1][:-1]
        else:
            with open(output + str(i) + ".txt", "r") as out:
                for line in out:
                    result = line.split(" ")
                    output_louvain[labels[result[0]]] += result[1][:-1] + ";"

    return output_louvain


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


def update_contigs_data(contig_data_file, core_bins, overlapping_bins, outdir):
    """Add bin information in the contigs data file.

    This function allow to update the contigs data file which were created
    previously in the network functions with the columns: contig id, contig
    name, contig length, GC content, hit, coverage. The function will add six
    columns: core bin id, core bin number of contigs, core bin length,
    overlapping bin id, overlapping bin number of contigs, overlapping bin
    length.

    The previous file will be overwritten.

    Parameters:
    -----------
    contig_data_file : str
        Path to the contigs data file.
    core_bins : dict
        Dictionnary which has as keys the values of the iterations from Louvain
        separated by a semicolon and as values the list of the id of the
        contigs.
    overlapping_bins : dict
        A dictionnary with the id of the overlapping bins as keys and the list
        of id of their contigs as values.
    outdir : str
        Path of the output directory to write the update contigs data file

    Returns:
    --------
    pandas.core.frame.DataFrame:
        Table with all the information on the contigs included their
        appartenance to the bins.
    """

    # Read the table
    contigs_data = pd.read_csv(
        contig_data_file, sep="\t", header=None, index_col=False
    )
    contigs_data.columns = [
        "id",
        "name",
        "length",
        "GC_content",
        "hit",
        "coverage",
    ]

    # Add new empty columns
    contigs_data["cc_id"] = "-"
    contigs_data["cc_nb"] = "-"
    contigs_data["cc_length"] = "-"
    contigs_data["oc_id"] = "-"
    contigs_data["oc_nb"] = "-"
    contigs_data["oc_length"] = "-"

    # Add core bin information
    n = 1
    for i in core_bins:
        # Extract contigs of the bin
        core_bin = [id - 1 for id in core_bins[i]]
        core_bin_data = contigs_data.iloc[core_bin]
        core_bin_contigs_number = len(core_bin)
        core_bin_length = sum(core_bin_data.length)
        # Write the new information
        contigs_data.iloc[core_bin, 6] = n
        contigs_data.iloc[core_bin, 7] = core_bin_contigs_number
        contigs_data.iloc[core_bin, 8] = core_bin_length
        n += 1

    # Add overlapping information
    for i in overlapping_bins:
        # Extract contigs of the bin
        overlapping_bin = [id - 1 for id in overlapping_bins[i]]
        overlapping_bin_data = contigs_data.iloc[overlapping_bin]
        overlapping_bin_contigs_number = len(overlapping_bin)
        overlapping_bin_length = sum(overlapping_bin_data.length)
        # Write the new information
        contigs_data.iloc[overlapping_bin, 9] = i
        contigs_data.iloc[overlapping_bin, 10] = overlapping_bin_contigs_number
        contigs_data.iloc[overlapping_bin, 11] = overlapping_bin_length

    # Write the new file
    contig_data_file_2 = join(outdir, "contig_data_partition.txt")
    contigs_data.to_csv(contig_data_file_2, sep="\t", header=None, index=False)

    return contigs_data
