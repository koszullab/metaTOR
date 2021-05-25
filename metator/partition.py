#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Partition metaHiC network using Louvain or Leiden algorithm.

General utility function for generating bins from the network file using Louvain
algorithm. If an overlapping threshold is given, it will use the Hamming
distance to group together bins relatively closed (to avaoid to split genomes in
different bins).

Core functions to partition the network are:
    - build_clustering_matrix
    - defined_overlapping_bins
    - detect_core_bins
    - generate_fasta
    - get_distances_splitmat
    - get_hamming_distance
    - leiden_partition_java
    - louvain_partition_cpp
    - partition
    - remove_isolates
    - update_contigs_data
"""

import metator.io as mio
import multiprocessing
import numpy as np
import os
import pandas as pd
import subprocess as sp
from functools import partial
from metator.log import logger
from os.path import join
from scipy import sparse
from sklearn import metrics


def build_clustering_matrix(core_bins_contigs, hamming_distance, N):
    """Function to return the clustering matrix in sparse format.

    For each contigs, the value correspond to the number of iterations where the
    contigs are clusterized together divided by the number of iterations. A
    value of 1 means that the contigs are in the same core bin.

    Parameters:
    -----------
    core_bins_contigs : dict
        Dictionnary which has as keys the core bins id and as value the id of
        the contigs of the core bin.
    hamming_distance : scipy.sparse.csr.csr_matrix:
        Matrix with all the previously computed hamming distance between two
        core bins.
    N : int
        Number of contigs in the assembly.

    Returns:
    --------
    scipy.sparse.coo.coo_matrix:
        Matrix with all the previously computed hamming distance between two
        contigs.
    """
    # To do it we build a transition matrix T which look like the identity
    # matrix but not square to extend our matrix as fellow: B = T.T * A * T
    rows = []
    cols = []
    values = []
    for core_bin in core_bins_contigs:
        for contig_id in core_bins_contigs[core_bin]:
            rows.append(core_bin)
            cols.append(contig_id)
            values.append(1)
    transition_matrix = sparse.coo_matrix(
        (values, (rows, cols)),
        shape=(len(core_bins_contigs), N + 1),
        dtype=np.int32,
    )
    # Compute the clustering matrix on only the upper triangle of the hamming
    # distance matrix as it's symmetric to reduce memory usage.
    hamming_distance = sparse.triu(hamming_distance, k=1)
    clustering_matrix = transition_matrix.T.dot(hamming_distance).dot(
        transition_matrix
    )
    return ((clustering_matrix + clustering_matrix) / 2).tocoo()


def defined_overlapping_bins(
    overlap, hamming_distance, core_bins_contigs, core_bins_iterations
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
        hamming distance threshold use to consider that two bins are
        overlapping.
    hamming_distance : scipy.sparse.csr.csr_matrix
        Matrix with all the previously computed hamming distance between two
        core bins.
    core_bins_contigs : dict
        Dictionnary which has as keys the core bins id and as value the id of
        the contigs of the core bin.
    core_bins_iteration : pandas.core.frame.DataFrame
        Table with the id of the core bin and their values for each iterations.

    Returns:
    --------
    dict:
        A dictionnary with the id of the overlapping bins as keys and the list
        of id of their contigs as values.
    """
    # Extract bins which are connected, i.e. bins with an hamming distance
    # superior than the threshold given. The small variation is necessary as
    # python give a float not really equal to the true value
    # (i.e. 0.1 -> 0.09999999999999998)
    connections = hamming_distance >= (overlap - 1e-10)
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
        core_bin_contigs = core_bins_contigs[cc_id].copy()
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


def detect_core_bins(output_partition, iterations):
    """Detect core bins from the output of the partition algorithm.

    The function search for duplicated values in the output of Louvain or Leiden
    algorithm in order to find contigs which are always in the same bin. The
    bins find with this method are called the core bins.

    Parameters:
    -----------
    output_partition : dict
        Dictionnary with the id of the contig as key and the list of the results
        of each iterations separated by a semicolon as values.
    iterations : int
        Number of iterations made previously with the partition algorithm.

    Returns:
    --------
    dict:
        Dictionnary which has as keys the core bins id and as value the id of
        the contigs of the core bin.
    pandas.core.frame.DataFrame:
        Table with the id of the core bin and their values for each iterations.
    """
    # finding duplicate values in the output of louvain or leiden using a
    # flipped dictionary.

    # Create dictionnary for core bins
    core_bins = {}
    core_bins_contigs = {}
    core_bins_iterations = np.empty((0, iterations), int)
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
            core_bins_iterations = np.append(
                core_bins_iterations,
                np.array([list(map(int, value.split(";")))]),
                axis=0,
            )
        # If already an entry created for this bin add a contig in the lists.
        else:
            core_bins_contigs[core_bins[value]].append(key)

    # Transform the array in a dataframe
    core_bins_iterations = pd.DataFrame(core_bins_iterations)

    logger.info("{0} core bins were found.\n".format(len(core_bins)))

    return core_bins_contigs, core_bins_iterations


def generate_fasta(
    assembly, overlapping_bins, contigs_data, size, output_dir, tmpdir
):
    """Generate the fasta files of each bins from the assembly.

    Parameters:
    -----------
    assembly : str
        Path to the fasta file of the original assembly.
    overlapping_bins : dict
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
    tmpdir : str
        Path to the temporary directory to write the temporary contigs list
        files.
    """

    nb_bins = 0
    length_bins = 0
    # For each bin create a list of the contigs and extract them from the
    # assembly to create a new fasta file with only the bin.
    for bin_id in overlapping_bins:
        # Extract the list of the contigs from the contigs data file.
        list_contigs_id = overlapping_bins[bin_id]
        list_contigs_name = []
        # Test if the bin is bigger than the size threshold given.
        length_bin = contigs_data.loc[
            list_contigs_id[0] - 1, "Overlapping_bin_size"
        ]
        if length_bin >= size:
            nb_bins += 1
            length_bins += length_bin
            for contig_id in list_contigs_id:
                list_contigs_name.append(
                    contigs_data.loc[contig_id - 1, "Name"]
                )
            # Define the output file.
            output_file = join(output_dir, "MetaTOR_{0}_0.fa".format(bin_id))
            # Create the fasta file.
            contigs_file = join(tmpdir, "MetaTOR_{0}_0.txt".format(bin_id))
            with open(contigs_file, "w") as f:
                for contig_name in list_contigs_name:
                    f.write("%s\n" % contig_name)
            cmd = "pyfastx extract {0} -l {1} > {2}".format(
                assembly, contigs_file, output_file
            )
            process = sp.Popen(cmd, shell=True)
            process.communicate()
    logger.info("{0} bins have been extracted".format(nb_bins))
    logger.info(
        "Total size of the extracted bins: {0}Mb".format(
            round(length_bins / 10 ** 6, 3)
        )
    )


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


def get_hamming_distance(core_bins_iterations, n_iter, threads):
    """Generate matrix of Hamming distances between all pairs of core bins.

    Parameters:
    -----------
    core_bins_iterations : pandas.core.frame.DataFrame
        Table with the id of the core bin and their values for each iterations.
    core_bins : dict
        Dictionnary which has as keys the values of the iterations from Louvain
        or Leiden separated by a semicolon and as value the id of the contigs of
        the core bin.
    n_iter : int
        Number of iterations made previously.
    threads : int
        Number of cores to parallelize computation.

    Returns:
    --------
    scipy.sparse.csr.csr_matrix:
        Matrix with all the previously computed hamming distance between two
        core bins.
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
    return res.tocsr()


def leiden_iterations_java(
    network_file, iterations, resolution_parameter, tmp_dir, leiden_path
):
    """Use the java implementation of Leiden to prtition the network.

    Parameters:
    -----------
    network_file : str
        Path to the network computed previously. The file is 3 columns table
        separated by a tabulation with the id of the first contigs the id of the
        second one and the weights of the edge normalized or not.
    iterations : int
        Number of iterations of the algorithm of Leiden.
    resolution_parameter : float
        Resolution parameter for Leiden clustering.
    tmp_dir : str
        Path to the temporary directory.
    leiden_path : str
        Path to the directory with network analysis java implementation.

    Returns:
    --------
    dict:
        Dictionnary with the id of the contig as key and the list of the results
        of each iterations separated by a semicolon as values.
    """
    output_partition = dict()

    # Run the iterations of Leiden
    for i in range(iterations):
        logger.info("Iteration in progress: {0}".format(i))

        output = join(tmp_dir, "partition_{0}.txt".format(i))

        # Clusterize the network using Leiden.
        cmd = (
            " java -cp {0} nl.cwts.networkanalysis.run.RunNetworkClustering -i 4 -r {1} -w -o {2} -q Modularity -a Leiden {3}"
        ).format(leiden_path, resolution_parameter, output, network_file)
        process = sp.Popen(cmd, shell=True, stderr=sp.DEVNULL)
        process.communicate()

        # Save the results in a dictionnary
        if i == 0:
            with open(output, "r") as out:
                for line in out:
                    result = line.split("\t")
                    output_partition[int(result[0])] = result[1][:-1]

        else:
            with open(output, "r") as out:
                for line in out:
                    result = line.split("\t")
                    output_partition[int(result[0])] += ";" + result[1][:-1]

    # Remove isolates (nodes with no contacts):
    output_partition.pop(0)
    output_partition = remove_isolates(output_partition, network_file)

    return output_partition


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
    tmp_dir : str
        Path to the temporary directory.
    louvain_path : str
        Path to the directory with louvain functions.

    Returns:
    --------
    dict:
        Dictionnary with the id of the contig as key and the list of the results
        of each iterations separated by a semicolon as values.
    """

    # Check if louvain cpp is available in the computer. If it's not available
    # launch python_louvain instead.
    if not mio.check_louvain_cpp(louvain_path):
        logger.error("Louvain implementation was not found.")
        logger.error(
            "You should have a LOUVAIN_PATH variable in your environnement"
        )
        raise NameError

    # Defined temporary files and args for louvain fonction calling and path to
    # the variables to call.
    network_bin = join(tmp_dir, "net_bin")
    network_weight = join(tmp_dir, "net_weight")
    network_tree = join(tmp_dir, "net_tree")
    network_labels = join(tmp_dir, "labels.txt")
    level_louvain = join(tmp_dir, "level.txt")
    output = join(tmp_dir, "output_louvain_")
    louvain = join(louvain_path, "louvain")
    convert = join(louvain_path, "convert")
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
        "convert": convert,
        "louvain": louvain,
        "hierarchy": hierarchy,
    }

    # Convert the file in binary file for Louvain partitionning.
    cmd = (
        "{convert} -i {net_txt} -o {net_bin} -r {net_labels} -w {net_weight}"
    ).format(**louvain_args)
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
        if i == 0:
            with open(output + str(i) + ".txt", "r") as out:
                for line in out:
                    result = line.split(" ")
                    output_louvain[labels[result[0]]] = result[1][:-1]

        else:
            with open(output + str(i) + ".txt", "r") as out:
                for line in out:
                    result = line.split(" ")
                    output_louvain[labels[result[0]]] += ";" + result[1][:-1]

    return output_louvain


def partition(
    algorithm,
    assembly,
    cluster_matrix,
    contig_data_file,
    iterations,
    network_file,
    outdir,
    fasta_dir,
    overlapping_parameter,
    resolution_parameter,
    size,
    temp_directory,
    threads,
):
    """Function to call the others functions to partition the network.

    Parameters:
    -----------
    algorithm : str
        Algorithm to use to partition the network. Either leiden or louvain.
    assembly : str
        Path to the assembly file used for the partition.
    cluster_matrix : bool
        If True, build and save the clustering matrix.
    contig_data_file : str
        Path to the contig data table to update.
    iterations : int
        Number of iterations to use for the partition.
    network_file : str
        Path to the network file.
    outdir : str
        Path to the output directory where to write the output files.
    fasta_dir : str
        Path to directory where to write the fasta files.
    overlapping_parameter : int
        Hamming distance threshold to use to merge bins (percentage).
    resolution_parameter : float
        Resolution parameter to use if Leiden algorithm is chosen. It will be a
        factor of the cost function used. A resolution parameter of 1 will be
        equivalent as the modularity function used in Louvain. Higher these
        parameters, smaller the bins will be in the output.
    size : int
        Threshold size in base pair of the output bins.
    temp_directory : str
        Path to the directory used to write temporary files.
    threads : int
        Number of threads to use.

    Returns:
    --------
    scipy.sparse.coo.coo_matrix:
        Matrix with all the previously computed hamming distance between two
        contigs.
    str:
        Path to the new contig data file with the bin informations in it.
    """

    # Create partition folders in the temporary directory
    temp_directory = join(temp_directory, "partition")
    os.makedirs(temp_directory, exist_ok=True)
    temp_directory_clustering = join(temp_directory, "clustering")
    os.makedirs(temp_directory_clustering, exist_ok=True)
    temp_directory_bins = join(temp_directory, "partition_bins")
    os.makedirs(temp_directory_bins, exist_ok=True)

    # Perform the iterations of Louvain or Leiden to partition the network.
    logger.info("Start iterations:")
    if algorithm == "leiden":
        LEIDEN_PATH = os.environ["LEIDEN_PATH"]
        output_partition = leiden_iterations_java(
            network_file,
            iterations,
            resolution_parameter,
            temp_directory_clustering,
            LEIDEN_PATH,
        )
    elif algorithm == "louvain":
        LOUVAIN_PATH = os.environ["LOUVAIN_PATH"]
        output_partition = louvain_iterations_cpp(
            network_file,
            iterations,
            temp_directory_clustering,
            LOUVAIN_PATH,
        )
    else:
        logger.error('algorithm should be either "louvain" or "leiden"')
        raise ValueError

    # Detect core bins
    logger.info("Detect core bins:")
    (
        core_bins_contigs,
        core_bins_iterations,
    ) = detect_core_bins(output_partition, iterations)

    # Compute the Hamming distance between core bins.
    logger.info("Detect overlapping bins:")
    hamming_distance = get_hamming_distance(
        core_bins_iterations,
        iterations,
        threads,
    )

    # Defined overlapping bins according to the threshold
    overlapping_bins = defined_overlapping_bins(
        overlapping_parameter,
        hamming_distance,
        core_bins_contigs,
        core_bins_iterations,
    )

    # Update the contigs_data_file.
    logger.info("Extract bins:")
    contigs_data, contigs_data_file = update_contigs_data(
        contig_data_file,
        core_bins_contigs,
        overlapping_bins,
        outdir,
    )

    # Generate Fasta file
    generate_fasta(
        assembly,
        overlapping_bins,
        contigs_data,
        size,
        fasta_dir,
        temp_directory_bins,
    )

    if cluster_matrix:
        # Build clustering matrix and save it.
        logger.info("Build  clustering matrix")
        clustering_matrix = build_clustering_matrix(
            core_bins_contigs, hamming_distance, len(contigs_data.ID)
        )
        # Save the clustering matrix
        clustering_matrix_file = join(outdir, "clustering_matrix_partition")
        sparse.save_npz(clustering_matrix_file, clustering_matrix)
    else:
        clustering_matrix_file = None

    return clustering_matrix_file, contigs_data_file


def remove_isolates(output_partition, network_file):
    """Remove isolates, i.e. nodes without any contacts in the network in the
    partition. This step is necessary as it will slow the further process of the
    communities. This function is only useful while using Leiden algorithm.

    Parameters:
    -----------
    output_partition : dict
        Dictionnary with the id of the contig as key and the list of the results
        of each iterations separated by a semicolon as values.
    network_file : str
        Path to the network computed previously. The file is 3 columns table
        separated by a tabulation with the id of the first contigs the id of the
        second one and the weights of the edge normalized or not.

    Returns:
    --------
    dict:
        Dictionnary with the id of the contig as key and the list of the results
        of each iterations separated by a semicolon as values without isolates.
    """
    nodes_presents = []
    with open(network_file, "r") as network:
        for line in network:
            line = line.split("\t")
            nodes_presents.append(int(line[0]))
            nodes_presents.append(int(line[1]))
    for i in range(1, max(nodes_presents)):
        if i not in nodes_presents:
            output_partition.pop(i)
    return output_partition


def update_contigs_data(
    contig_data_file, core_bins_contigs, overlapping_bins, outdir
):
    """Add bin information in the contigs data file.

    This function allow to update the contigs data file which were created
    previously in the network functions with the columns: contig id, contig
    name, contig length, GC content, hit, coverage, restriction site. The
    function will add six columns: core bin id, core bin number of contigs, core
    bin length, overlapping bin id, overlapping bin number of contigs,
    overlapping bin length.

    Parameters:
    -----------
    contig_data_file : str
        Path to the contigs data file.
    core_bins_contigs : dict
        Dictionnary which has as keys the core bins id and as value the id of
        the contigs of the core bin.
    overlapping_bins : dict
        A dictionnary with the id of the overlapping bins as keys and the list
        of id of their contigs as values.
    outdir : str
        Path of the output directory to write the update contigs data file.

    Returns:
    --------
    pandas.core.frame.DataFrame:
        Table with all the information on the contigs included their
        appartenance to the bins.
    """

    # Read the table
    contigs_data = pd.read_csv(
        contig_data_file, sep="\t", header=0, index_col=False
    )

    # Add new empty columns
    contigs_data["Core_bin_ID"] = "-"
    contigs_data["Core_bin_contigs"] = "-"
    contigs_data["Core_bin_size"] = "-"
    contigs_data["Overlapping_bin_ID"] = "-"
    contigs_data["Overlapping_bin_contigs"] = "-"
    contigs_data["Overlapping_bin_size"] = "-"

    # Add core bin information
    for i in core_bins_contigs:
        # Extract contigs of the bin
        core_bin = [id - 1 for id in core_bins_contigs[i]]
        core_bin_data = contigs_data.iloc[core_bin]
        core_bin_contigs_number = len(core_bin)
        core_bin_length = sum(core_bin_data.Size)
        # Write the new information
        contigs_data.loc[core_bin, "Core_bin_ID"] = i + 1
        contigs_data.loc[core_bin, "Core_bin_contigs"] = core_bin_contigs_number
        contigs_data.loc[core_bin, "Core_bin_size"] = core_bin_length

    # Add overlapping information
    for i in overlapping_bins:
        # Extract contigs of the bin
        overlapping_bin = [id - 1 for id in overlapping_bins[i]]
        overlapping_bin_data = contigs_data.iloc[overlapping_bin]
        overlapping_bin_contigs_number = len(overlapping_bin)
        overlapping_bin_length = sum(overlapping_bin_data.Size)
        # Write the new information
        contigs_data.loc[overlapping_bin, "Overlapping_bin_ID"] = i
        contigs_data.loc[
            overlapping_bin, "Overlapping_bin_contigs"
        ] = overlapping_bin_contigs_number
        contigs_data.loc[
            overlapping_bin, "Overlapping_bin_size"
        ] = overlapping_bin_length

    # Write the new file
    contig_data_file_2 = join(outdir, "contig_data_partition.txt")
    contigs_data.to_csv(contig_data_file_2, sep="\t", header=True, index=False)

    return contigs_data, contig_data_file_2
