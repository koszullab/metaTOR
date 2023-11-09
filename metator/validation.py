#!/usr/bin/env python3
# coding: utf-8

"""Validates bins using miComplete and make a recursive partition to try to
decontaminate them.

General functions to validate bins completion using miComplete and make 
recursive iterations of Louvain or Leiden to try to partition contaminated bins.
Only bins with more than 50% completion and 5% contamination are subject to the
recursive step. If the recursive step gave worst results than the first 
(decrease of the completion with no decrease of the contamination), it will keep
the original bin.


Functions in this module:
    - correct_final_bin
    - get_bin_coverage
    - give_results_info
    - merge_micomplete
    - micomplete_compare_bins
    - micomplete_quality
    - recursive_clustering
    - recursive_clustering_worker
    - recursive_decontamination
    - update_contigs_data_recursive
    - write_bins_contigs

CheckM deprecated functions:
    - checkm
    - checkm_compare_bins
"""

import logging
import metator.io as mio
import metator.figures as mtf
import metator.partition as mtp
import micomplete
import multiprocessing
import networkx as nx
import numpy as np
import os
import pandas as pd
import pyfastx
import shutil
import subprocess as sp
from functools import partial
from metator.log import logger
from os.path import join
from scipy import sparse


def correct_final_bin(contigs_data, final_fasta_dir, bin_summary):
    """Function to compute the coverage of each bin.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the informations of the final bins kept by MetaTOR.
    contigs_data : pandas.core.frame.DataFrame
        Dataframe with the contigs informations.

    Returns:
    --------
    pandas.core.frame.DataFrame
        Dataframe with the contigs informations.
    """
    contigs_data = contigs_data.set_index("Name", drop=False)
    for bin_name in bin_summary:
        fasta_file = join(final_fasta_dir, f"{bin_name}.fa")
        fasta = pyfastx.Fasta(fasta_file)
        for seq in fasta:
            contigs_data.loc[seq.name, "Final_bin"] = bin_name
    contigs_data = contigs_data.reset_index(drop=True)
    return contigs_data


def get_bin_coverage(bin_summary, contigs_data):
    """Function to compute the coverage of each bin.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the informations of the final bins kept by MetaTOR.
    contigs_data : pandas.core.frame.DataFrame
        Dataframe with the contigs informations.

    Returns:
    --------
    dict:
        Dictionnary with the informations of the final bins kept by MetaTOR with
        the coverage.
    """
    # Compute HiC_coverage
    total_hic_hit = 0
    total_sg_hit = 0
    for i in contigs_data.index:
        bin_name = contigs_data.loc[i, "Final_bin"]
        if contigs_data.loc[contigs_data.index[0], "Shotgun_coverage"] != "-":
            total_sg_hit += (
                contigs_data.loc[i, "Size"]
                * contigs_data.loc[i, "Shotgun_coverage"]
            )
        total_hic_hit += contigs_data.loc[i, "Hit"]
        if bin_name != "ND":
            try:
                bin_summary[bin_name]["HiC_abundance"] += (
                    100 * contigs_data.loc[i, "Hit"]
                )
            except KeyError:
                bin_summary[bin_name]["HiC_abundance"] = (
                    100 * contigs_data.loc[i, "Hit"]
                )

            # If no depth files were given do not compute the Shotgun coverage.
            if (
                contigs_data.loc[contigs_data.index[0], "Shotgun_coverage"]
                != "-"
            ):
                try:
                    bin_summary[bin_name]["SG_abundance"] += (
                        contigs_data.loc[i, "Size"]
                        * contigs_data.loc[i, "Shotgun_coverage"]
                        * 100
                    )
                except KeyError:
                    bin_summary[bin_name]["SG_abundance"] = (
                        contigs_data.loc[i, "Size"]
                        * contigs_data.loc[i, "Shotgun_coverage"]
                        * 100
                    )
    for bin_name in bin_summary:
        # Divide the HiC abundance by two as the hit are counted twice.
        bin_summary[bin_name]["HiC_abundance"] /= 2 * total_hic_hit
        if contigs_data.loc[contigs_data.index[0], "Shotgun_coverage"] != "-":
            bin_summary[bin_name]["SG_abundance"] /= total_sg_hit
    return bin_summary


def give_results_info(bin_summary):
    """Function to return the general information about the binning results.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the summary results of the kept bins.
    """

    # Defined categories of the bins
    HQ = 0  # Completness >= 90 and Contamination <= 5
    total_size_HQ = 0
    MQ = 0  # Completness >= 70 and Contamination <= 10
    total_size_MQ = 0
    LQ = 0  # Completness >= 50 and Contamination <= 10
    total_size_LQ = 0
    conta_bins = 0  # Completness >= 50 and Contamination > 10
    total_size_conta_bins = 0
    others = 0  # Not determined bins.
    total_size_others = 0

    # Class each bin in a category
    for bin_name in bin_summary:
        completness = float(bin_summary[bin_name]["Weighted completeness"])
        contamination = float(bin_summary[bin_name]["Weighted redundancy"])
        size = int(bin_summary[bin_name]["Length"])
        if completness >= 0.5:
            if ((contamination - 1) / completness) > 0.1:
                conta_bins += 1
                total_size_conta_bins += size
            else:
                if (
                    completness >= 0.9
                    and ((contamination - 1) / completness) > 0.05
                ):
                    HQ += 1
                    total_size_HQ += size
                elif completness >= 0.7:
                    MQ += 1
                    total_size_MQ += size
                else:
                    LQ += 1
                    total_size_LQ += size
        else:
            others += 1
            total_size_others += size
    total = HQ + MQ + LQ + conta_bins + others
    total_size = (
        total_size_HQ
        + total_size_MQ
        + total_size_LQ
        + total_size_conta_bins
        + total_size_others
    )

    # Return info in the logger:
    logger.info(
        "{0} bins have been kept after the recursive iterations.".format(total)
    )
    logger.info("Total size of the extracted bins: {0}".format(total_size))
    logger.info("HQ MAGs: {0}\tTotal Size: {1}".format(HQ, total_size_HQ))
    logger.info("MQ MAGs: {0}\tTotal Size: {1}".format(MQ, total_size_MQ))
    logger.info("LQ MAGs: {0}\tTotal Size: {1}".format(LQ, total_size_LQ))
    logger.info(
        "Contaminated potential MAGs: {0}\tTotal Size: {1}".format(
            conta_bins, total_size_conta_bins
        )
    )
    logger.info(
        "Others bins: {0}\tTotal Size: {1}".format(others, total_size_others)
    )


def merge_micomplete(out_bact105, out_arch131, outfile):
    """Function to merge bacterial and archaeal output of micomplete into one
    file.

    Parameters:
    out_bact105 : str
        Output of micomplete using 105 bacterial markers.
    out_arch131 : str
        Output of micomplete using 131 arcaheal markers.
    outfile : str
        Final merged output.
    """
    # Reads both files.
    bact105 = pd.read_csv(out_bact105, sep="\t", comment="#", index_col=0).loc[
        :, :"CDs"
    ]
    arch131 = pd.read_csv(out_arch131, sep="\t", comment="#", index_col=0).loc[
        :, :"CDs"
    ]

    # Write header
    with open(outfile, "w") as out:
        out.write("## miComplete\n")
        out.write(f"## v{micomplete.__version__}\n")
        out.write(f"## Weights: Bact105 and Arch131\n")
        out.write(
            "Name\t{0}\tMarkers\n".format("\t".join(list(arch131.columns)))
        )
        for bin_id in bact105.index:
            if (
                bact105.loc[bin_id, "Weighted completeness"]
                >= arch131.loc[bin_id, "Weighted completeness"]
            ):
                out.write(
                    "{0}\t{1}\tBacteria\n".format(
                        bin_id, "\t".join(map(str, list(bact105.loc[bin_id])))
                    )
                )
            else:
                out.write(
                    "{0}\t{1}\tArchaea\n".format(
                        bin_id, "\t".join(map(str, list(arch131.loc[bin_id])))
                    )
                )


def micomplete_compare_bins(
    recursive_micomplete_file, bin_summary, parent_dict, step,
):
    """Compare the completness and contamination of the bins from the first step
    and from the recursive step. If the recursive step decrease the completion
    of one bin without decreasing its contamination, it will kept the bin before
    the recursive step. Moreover if the completion goes below 50% with the
    recursive it will not take the recursive bins. The goal is to keep complete
    bins even if there are contaminated that the users can manually curate.

    Parameters:
    -----------
    recursive_micomplete_file : str
        Path to the miComplete summary from the recursive step.
    bin_summary : dict
        Dictionnary containing iinformation about the bins.
    parent_dict : dict
        Dictionnary with recursive bin_id as key and parent bin as values.
    step : int
        Recursive step.

    Returns:
    --------
    dict:
        Dictionnary with the informations of the final bins kept by MetaTOR.
    """
    # Setup contamination
    contamination = False

    # Load the miComplete summary
    micomplete_recursive_summary = pd.read_csv(
        recursive_micomplete_file, sep="\t", comment="#", index_col=0,
    ).iloc[:, :13]

    bin_summary_tab = pd.DataFrame.from_dict(bin_summary, orient="index")

    # Create new columns for recursive values
    bin_summary_tab["max_rec_completness"] = np.nan
    bin_summary_tab["rec_id"] = np.nan

    # Retrieve maximum completness of the recursive bins.
    for recursive_bin in micomplete_recursive_summary.index:
        parent_bin = parent_dict[recursive_bin]
        bin_summary_tab.loc[parent_bin, "max_rec_completness"] = max(
            float(
                micomplete_recursive_summary.loc[
                    recursive_bin, "Weighted completeness"
                ]
            ),
            bin_summary_tab.loc[parent_bin, "max_rec_completness"],
        )
        try:
            bin_summary_tab.loc[parent_bin, "rec_id"].append(recursive_bin)
        except AttributeError:
            bin_summary_tab.loc[parent_bin, "rec_id"] = [recursive_bin]

    # If there are some recursive bins which have not loose too much completion
    # write their information otherwise keep the original bins.
    for overlapping_bin in bin_summary_tab.index:
        if np.isnan(
            bin_summary_tab.loc[overlapping_bin, "max_rec_completness"]
        ):
            bin_summary[overlapping_bin]["recursive"] = False
        else:
            max_rec = bin_summary_tab.loc[
                overlapping_bin, "max_rec_completness"
            ]
            over = float(
                bin_summary_tab.loc[overlapping_bin, "Weighted completeness"]
            )
            if (max_rec - 0.4) > ((over - 0.4) / 1.6):
                rec_ids = bin_summary_tab.loc[overlapping_bin, "rec_id"]
                # Case of only one bin.
                if isinstance(rec_ids, str):
                    rec_ids = [rec_ids]
                for rec_id in rec_ids:
                    bin_summary[rec_id] = micomplete_recursive_summary.loc[
                        rec_id, :"CDs"
                    ].to_dict()
                    bin_summary[rec_id]["recursive"] = True
                    bin_summary[rec_id]["step"] = step
                    bin_summary[rec_id]["parent"] = overlapping_bin
                    # If one bin from recursive step contamination set to true
                    # to continue the while loop.
                    contamination = True
                bin_summary.pop(overlapping_bin, "None")
            else:
                bin_summary[overlapping_bin]["recursive"] = False

    return bin_summary, contamination


def micomplete_quality(fasta_dir, outfile, threads):
    """Function to evaluate fasta bins using miComplete. Write the bins quality
    summury in the outfile.

    Parameters:
    -----------
    fasta_dir : str
        Path to the input fasta of the bins to evaluate.
    outfile : str
        Path to the file where the results of miComplete will be written.
    threads : int
        Numbers of threads to use for miComplete.
    """
    # Prepare input table for micomplete.
    list_fasta = filter(
        lambda x: ".fa" in x,
        [join(fasta_dir, path) for path in os.listdir(fasta_dir)],
    )
    tmp_seq_tab = join(fasta_dir, "micomplete_seq.tsv")
    with open(tmp_seq_tab, "w") as tab:
        for fasta in list_fasta:
            tab.write(f"{fasta}\tfna\n")

    # Create two temporary output for archaea and bacteria.
    out_bact105 = join(fasta_dir, "micomplete_bact105.tsv")
    out_arch131 = join(fasta_dir, "micomplete_arch131.tsv")

    # Launch miComplete using subprocess for bacteria.
    cmd = f"miComplete {tmp_seq_tab} --hmms Bact105 --weights Bact105 --threads {threads} --outfile {out_bact105}"
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()

    # Launch miComplete using subprocess for archaea.
    cmd = f"miComplete {tmp_seq_tab} --hmms Arch131 --weights Arch131 --threads {threads} --outfile {out_arch131}"
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()

    # Merge bacteria and archaea
    merge_micomplete(out_bact105, out_arch131, outfile)

    # Remove micomplete temporary files.
    for path in filter(lambda x: ".fa" in x, os.listdir(fasta_dir)):
        name = path.split(".")[0]
        os.remove(f"{name}_prodigal.faa")
        os.remove(f"{name}.tblout")
    os.remove(tmp_seq_tab)
    os.remove(out_bact105)
    os.remove(out_arch131)
    try:
        os.remove("miComplete.log")
    except FileNotFoundError:
        pass


def recursive_clustering(
    assembly,
    iterations,
    overlapping_parameter,
    resolution_parameter,
    outdir,
    recursive_fasta_dir,
    algorithm,
    tmpdir,
    bin_summary,
    contigs_data,
    network,
    cluster_matrix,
    size,
    threads,
    prefix,
):
    """Function to run recursive iterations on contaminated bins in order to try
    to improve the quality of the bins using Louvain or Leiden algorithm.

    Parameters:
    -----------
    assembly : str
        Path to the fasta file used as assembly.
    iterations : int
        Number of iterations to use for recursive iterations of Louvain or
        Leiden.
    overlapping_parameter : float
        Hamming distance threshold to consider two bins as the same bin.
    resolution parameter : float
        Resolution parameter of Leiden algorithm.
    outdir : str
        Path to the output directory.
    recursive_fasta_dir : str
        Path to the directory where to write the decontaminated fasta.
    algorithm : str
        Algorithm to use, either louvain, leiden or spinglass.
    tmpdir : str
        Path the temp directory.
    bin_summary : dict
        Dictionary containing information about the bins.
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs.
    network : str
        Metator full network.
    cluster_matrix : bool
        If True, build the clustering matrix and save it.
    size : int
        Size threshold in base pairs of the bins.
    threads : int
        Number of threads to use.
    prefix : str
        Sample prefix to use.

    Returns:
    --------
    boolean:
        True if at least one new recursive bin has been generated.
    pandas.DataFrame
        Updated dictionnary which has as keys the values of the iterations from
        the recursive partition separated by a semicolon and as values the list
        of the id of the contigs.
    scipy.sparse.coo.coo_matrix:
        Matrix with all the previously computed hamming distance between two
        contigs.
    dict
        Dictionnary with recursive bin_id as key and parent bin as values.
    """

    # Create temporary folders.
    tmpdir_binning = join(tmpdir, "recursive_bins")
    os.makedirs(tmpdir_binning, exist_ok=True)

    # Default no contamination
    contamination = False

    # Create an empty matrix
    N = len(contigs_data.ID)
    clustering_matrix = sparse.coo_matrix((N + 1, N + 1), dtype=np.float32)

    # Put iterations to on eif spinglass partition used.
    if algorithm == "spinglass":
        iterations = 1
        logger.error("Spinglass is no longer maintained.")
        raise ValueError

    # Check bin_id to decontaminate.
    bin_ids = []
    for bin_id in bin_summary:
        try:
            completness = float(bin_summary[bin_id]["Weighted completeness"])
            conta = float(bin_summary[bin_id]["Weighted redundancy"])
            recursive = bin_summary[bin_id]["recursive"]
            if recursive:
                if completness >= 0.33:
                    if (conta - 1) / completness >= 0.05:
                        bin_ids.append(bin_id)
        except KeyError:
            continue

    # Iterate on cmicomplete summary to find conatminated bins:
    # pool = multiprocessing.Pool(processes=threads)
    output_partitions = list(
        map(
            partial(
                recursive_clustering_worker,
                bin_summary=bin_summary,
                tmpdir=tmpdir,
                network=network,
                algorithm=algorithm,
                iterations=iterations,
                resolution_parameter=resolution_parameter,
                contigs_data=contigs_data,
            ),
            bin_ids,
        )
    )

    parent_dict = dict()

    for i, bin_id in enumerate(bin_ids):
        output_partition = output_partitions[i]

        # Detect core bins
        (
            recursive_core_bins,
            recursive_bins_iterations,
        ) = mtp.detect_core_bins(output_partition, iterations)

        # Compute the Hamming distance between core bins.
        hamming_distance = mtp.get_hamming_distance(
            recursive_bins_iterations, threads,
        )

        # Defined overlapping bins according to the threshold
        recursive_bins = mtp.defined_overlapping_bins(
            overlapping_parameter, hamming_distance, recursive_core_bins,
        )

        # update bin data and generate fasta
        (
            contamination,
            contigs_data,
            parent_dict,
        ) = update_contigs_data_recursive(
            bin_id,
            contigs_data,
            recursive_bins,
            assembly,
            recursive_fasta_dir,
            tmpdir_binning,
            size,
            contamination,
            parent_dict,
            prefix,
        )

        # Build the clustering matrix of the subnetwork and add it.
        if cluster_matrix > 0:
            clustering_matrix += mtp.build_clustering_matrix(
                recursive_core_bins, hamming_distance, cluster_matrix
            )

        logger.info("Recursive step for {0} is done.".format(bin_id))

    # Save the clustering matrix
    if cluster_matrix:
        clustering_matrix_file = join(outdir, "clustering_matrix_recursive")
        sparse.save_npz(clustering_matrix_file, clustering_matrix)
    else:
        clustering_matrix_file = None

    return contamination, contigs_data, clustering_matrix_file, parent_dict


def recursive_clustering_worker(
    bin_id,
    bin_summary,
    tmpdir,
    network,
    algorithm,
    iterations,
    resolution_parameter,
    contigs_data,
):
    """Worker to partition one bin if it's contaminated.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary of the output of micomplete as values and the bin id as
        keys.
    tmpdir : str
        Path the temp directory.
    algorithm : str
        Algorithm to use, either louvain, leiden or spinglass.
    network : networkx.classes.graph.Graph
        Network of contigs from HiC librairies.
    iterations : int
        Number of iterations to use for recursive iterations of Louvain or
        Leiden.
    resolution parameter : float
        Resolution parameter of Leiden algorithm.
    contigs_data_file : str
        Path to the contigs data file from metator partition.
    """
    # Create temporary folders.
    tmpdir_subnetwork = join(tmpdir, "recursive_bins", bin_id)
    os.makedirs(tmpdir_subnetwork, exist_ok=True)
    tmpdir_clustering = join(tmpdir, "recursive_clustering", bin_id)
    os.makedirs(tmpdir_clustering, exist_ok=True)

    logger.info("Bin in progress: {0}".format(bin_id))
    subnetwork_file = join(tmpdir_subnetwork, "subnetwork_" + bin_id + ".txt")
    over_bin_id = str(bin_id.split("_")[-2])
    rec_bin_id = str(bin_id.split("_")[-1])

    # Extract contigs
    mask = (contigs_data["Overlapping_bin_ID"] == over_bin_id) & (
        contigs_data["Recursive_bin_ID"].apply(str) == rec_bin_id
    )
    list_contigs = list(contigs_data.loc[mask, "ID"])

    # Extract subnetwork
    subnetwork = network.subgraph(list_contigs)

    # Write the new subnetwork
    nx.write_edgelist(
        subnetwork, subnetwork_file, delimiter="\t", data=["weight"]
    )

    # Compute spin prediction on the completion/contamination values.
    spin = max(2, int(1 + float(bin_summary[bin_id]["Weighted redundancy"])),)
    # Partition the subnetwork.

    output_partition = mtp.algo_partition(
        algorithm,
        subnetwork_file,
        subnetwork,
        iterations,
        resolution_parameter,
        tmpdir_clustering,
        spin,
    )
    logger.info("Output partition done for {0}.".format(bin_id))

    return output_partition


def recursive_decontamination(
    algorithm,
    assembly,
    cluster_matrix,
    contig_data_file,
    final_fasta_dir,
    input_fasta_dir,
    iterations,
    network_file,
    outdir,
    overlapping_parameter,
    recursive_fasta_dir,
    resolution_parameter,
    size,
    temp_directory,
    threads,
    prefix,
):
    """Function to validate bins do the recursive decontamination using Louvain
    or Leiden algorithm

    Parameters:
    -----------
    algorithm : str
        Algorithm to use to recursively partition the network. Either leiden,
        louvain or spinglass.
    assembly : str
        Path to the assembly file used for the partition.
    cluster_matrix : bool
        If True, build the clustering matrix and save it.
    contig_data_file : str
        Path to the contig data table to update.
    final_fasta_dir : str
        Path to write the final fasta decontaminated bins.
    input_fasta_dir : str
        Path to the directory where the fasta bin from the partition are.
    iterations : int
        Number of iterations to use for the recursive partition.
    network_file : str
        Path to the network file.
    outdir : str
        Path to the output directory where to write the output files.
    overlapping_parameter : int
        Hamming distance threshold in percentage to use to consider to bins as
        one in the recursive partition.
    recursive_fasta_dir : str
        Path to write the fasta decontaminated bins.
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
    prefix : str
        Sample prefix to use.

    Returns:
    --------
    scipy.sparse.coo.coo_matrix:
        Matrix with all the previously computed hamming distance between two
        contigs.
    """

    # Create folders in the temporary directory
    tmpdir_recursive_clustering = join(temp_directory, "recursive_clustering")
    os.makedirs(tmpdir_recursive_clustering, exist_ok=True)

    # Defined miComplete output file path
    overlapping_micomplete_file = join(
        outdir, "overlapping_micomplete_results.txt"
    )

    logger.info("Lauch miComplete quality check.")

    # Launch miComplete
    micomplete_quality(
        input_fasta_dir, overlapping_micomplete_file, threads,
    )

    # Load network:
    network = nx.read_edgelist(
        network_file, nodetype=int, data=(("weight", float),)
    )

    # Load contigs data:
    contigs_data = pd.read_csv(
        contig_data_file,
        sep="\t",
        header=0,
        converters={"Overlapping_bin_ID": str},
    )
    contigs_data["index"] = contigs_data["ID"] - 1
    contigs_data = contigs_data.set_index("index")

    # Add new coulumns for recursive information.
    contigs_data["Recursive_bin_ID"] = f"{0:05d}"
    contigs_data["Recursive_bin_contigs"] = "-"
    contigs_data["Recursive_bin_size"] = "-"
    contigs_data["Final_bin"] = "ND"

    # Load miComplete result:
    bin_summary = mio.micomplete_results_to_dict(overlapping_micomplete_file)

    # Add columns about recursive step.
    for bin_id in bin_summary:
        bin_summary[bin_id]["recursive"] = True
        bin_summary[bin_id]["step"] = 0
        bin_summary[bin_id]["parent"] = None

    # Recursively remove contamination.
    contamination = True
    step = 1

    logger.info("Starts recursive decontamition step:")

    while contamination == True:
        # Create fasta dir.
        recursive_fasta_dir_step = join(recursive_fasta_dir, f"step_{step}")
        os.makedirs(recursive_fasta_dir_step, exist_ok=True)

        # Stop to report info log
        logger.setLevel(logging.WARNING)

        # Iterates Louvain or Leiden on contaminated and complete bins.
        (
            contamination,
            contigs_data,
            clustering_matrix_file,
            parent_dict,
        ) = recursive_clustering(
            assembly,
            iterations,
            overlapping_parameter,
            resolution_parameter,
            outdir,
            recursive_fasta_dir_step,
            algorithm,
            tmpdir_recursive_clustering,
            bin_summary,
            contigs_data,
            network,
            cluster_matrix,
            size,
            threads,
            prefix,
        )

        # Put back the info log
        logger.setLevel(logging.INFO)

        # Recursive iterations of Louvain or Leiden on the contaminated bins.
        # Save bin information if the new bins have the same quality otherwise
        # keep the original bin information.
        if contamination:
            # Run miComplete on the recursive bins.
            recursive_micomplete_file = join(
                outdir, f"recursive_micomplete_step_{step}.txt"
            )

            micomplete_quality(
                recursive_fasta_dir_step, recursive_micomplete_file, threads,
            )

            # Compare
            bin_summary, contamination = micomplete_compare_bins(
                recursive_micomplete_file, bin_summary, parent_dict, step,
            )

        # Keep overlapping bin information
        else:
            logger.info(
                f"No more contaminated bin have been found after {step} steps."
            )
        # Increase count
        logger.info(f"Recursive step {step} done.")
        step += 1
        if step == 10:
            contamination = False

    # Create fasta directory and copy final bins.
    for bin_name in bin_summary:
        dst = join(final_fasta_dir, bin_name + ".fa")
        step = bin_summary[bin_name]["step"]
        if step == 0:
            src = join(input_fasta_dir, bin_name + ".fa")
        else:
            src = join(recursive_fasta_dir, f"step_{step}", f"{bin_name}.fa")
        shutil.copyfile(src, dst)

    # Return some values of efficiency of the binning.
    give_results_info(bin_summary)

    # Correct final bin value in contigs data
    contigs_data = correct_final_bin(contigs_data, final_fasta_dir, bin_summary)

    # Write relevant bins/contigs information for anvio.
    binning_file = join(outdir, "binning.txt")
    contigs_data = write_bins_contigs(
        bin_summary, contigs_data, binning_file, prefix
    )

    # Compute the abundance of the mags.
    bin_summary = get_bin_coverage(bin_summary, contigs_data)

    # Save bin information in final file
    bin_summary_file = join(outdir, "bin_summary.txt")
    mio.write_bin_summary(bin_summary, bin_summary_file)

    # Write the new file
    contig_data_file_final = join(outdir, "contig_data_final.txt")
    contigs_data.to_csv(
        contig_data_file_final, sep="\t", header=True, index=True
    )

    # Plot some figures of contigs distribution inside bins:
    mtf.plot_figures(outdir, contigs_data, bin_summary, size)

    return clustering_matrix_file


def update_contigs_data_recursive(
    bin_id,
    contigs_data,
    recursive_bins,
    assembly,
    outdir,
    tmpdir,
    size,
    contamination,
    parent_dict,
    prefix,
):
    """Update the data of the bin according to the recursive step and generated
    their fasta.

    Parameters:
    -----------
    bin_id : str
        Name of the parental bin.
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs.
    recursive_bins : dict
        Dictionnary which has as keys the values of the recursive iterations
        from Louvain or Leiden separated by a semicolon and as values the list
        of the id of the contigs.
    assembly : str
        Path to the fasta file.
    outdir : str
        Path to the output directory to write the new fasta.
    tmpdir : str
        Path to the file where to write the temporary contigs list.
    size : int
        Size threshold to generate fasta.
    contamination : boolean
        True if one bin has already been generated, false otherwise.
    parent_dict : dict
        Dictionnary with recursive bin_id as key and parent bin as values.
    prefix : str
        Sample prefix to use.

    Returns:
    --------
    boolean:
        True if one bin has already been generated, false otherwise.
    pandas.DataFrame
        Updated dictionnary which has as keys the values of the iterations from
        the recursive partition separated by a semicolon and as values the list
        of the id of the contigs.
    dict
        Dictionnary with recursive bin_id as key and parent bin as values.
    """

    # Add recursive bin information
    rec_id = 1

    # Extract last recursive ID.
    over_id = contigs_data.loc[recursive_bins[1][0] - 1, "Overlapping_bin_ID"]
    max_rec_id = max(
        map(
            int,
            contigs_data[contigs_data.Overlapping_bin_ID == over_id][
                "Recursive_bin_ID"
            ],
        )
    )

    for i in recursive_bins:
        # Extract contigs of the bin
        recursive_bin = [id - 1 for id in recursive_bins[i]]
        recursive_bin_data = contigs_data.loc[recursive_bin]
        recursive_bin_contigs_number = len(recursive_bin)
        recursive_bin_length = sum(recursive_bin_data.Size)

        if recursive_bin_length > size:
            # If one new bin is generated change the boolean value to True
            contamination = True
            rec_final_id = max_rec_id + rec_id
            # Write the new information
            contigs_data.loc[
                recursive_bin, "Recursive_bin_ID"
            ] = f"{rec_final_id:05d}"
            contigs_data.loc[
                recursive_bin, "Recursive_bin_contigs"
            ] = recursive_bin_contigs_number
            contigs_data.loc[
                recursive_bin, "Recursive_bin_size"
            ] = recursive_bin_length

            # Defined name of the recursive bin
            oc_id = int(
                contigs_data.loc[recursive_bin[0], "Overlapping_bin_ID"]
            )
            output_file = join(
                outdir, f"{prefix}_{oc_id:05d}_{rec_final_id:05d}.fa"
            )

            # Update bin_summary:
            parent_dict[f"{prefix}_{oc_id:05d}_{rec_final_id:05d}"] = bin_id

            # Retrieve names of the contigs
            list_contigs = list(contigs_data.loc[recursive_bin, "Name"])

            # Generate the fasta
            contigs_file = join(
                tmpdir, f"{prefix}_{oc_id:05d}_{rec_final_id:05d}.txt"
            )
            with open(contigs_file, "w") as f:
                for contig in list_contigs:
                    f.write("{0}\n".format(contig))
            cmd = "pyfastx extract {0} -l {1} > {2}".format(
                assembly, contigs_file, output_file
            )
            process = sp.Popen(cmd, shell=True)
            process.communicate()

            # Add one to the recursive id
            rec_id += 1

    return contamination, contigs_data, parent_dict


def write_bins_contigs(bin_summary, contigs_data, outfile, prefix):
    """Function to write a table with the nodes kept in the bins and their bin
    id. The file is adapted to be added in anvio.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the informations of the final bins kept by MetaTOR.
    contigs_data : pandas.core.frame.DataFrame
        Dataframe with the contigs informations.
    outfile : str
        Path where to write the output file.
    prefix : str
        Sample prefix to use.

    Returns:
    --------
    pandas.core.frame.DataFrame
        Dataframe with the contigs informations with the final bin information.
    """

    # Create a list with the id of the bins
    list_bin_id = dict()
    for bin_name in bin_summary:
        over_id = bin_name.split("_")[-2]
        rec_id = bin_name.split("_")[-1]
        try:
            list_bin_id[over_id].append(rec_id)
        except KeyError:
            list_bin_id[over_id] = [rec_id]

    # Write the contigs id with their bins id in table file
    with open(outfile, "w") as f:
        for i in contigs_data.index:
            over_id = str(contigs_data.loc[i, "Overlapping_bin_ID"])
            rec_id = str(contigs_data.loc[i, "Recursive_bin_ID"])
            try:
                rec_ids = list_bin_id[over_id]
                binned = False
                # Case of a recursive bin
                if rec_id in rec_ids:
                    binned = True
                # Case where the recursive bins where not kept.
                elif rec_ids == ["0"]:
                    rec_id = "0"
                    binned = True

                if binned:
                    final_bin = f"{prefix}_{over_id}_{rec_id}"
                    contigs_data.loc[i, "Final_bin"] = final_bin
                    f.write(
                        "{0}\t{1}\n".format(
                            contigs_data.loc[i, "Name"], final_bin
                        )
                    )
            except KeyError:
                pass
    return contigs_data


#############################################
# Deprecated functions for checkM validation.
#############################################


def checkm(fasta_dir, outfile, taxonomy_file, tmpdir, threads):
    """Function to evaluate fasta bins using CheckM. Write the checkM results
    summary in the outfile and the taxonomy results in the the taxonomy file.

    Parameters:
    -----------
    fasta_dir : str
        Path to the input fasta of the bins to evaluate.
    outfile : str
        Path to the file where the results of checkm will be written.
    taxonomy_file : str
        path to the file where checkm taxonomy results will be written.
    tmpdir : str
        Path to the temporary directory where CheckM intermediary files will be
        written.
    threads : int
        Numbers of threads to use for CheckM.
    """

    logger.info("Start CheckM validation.")

    # Build CheckM tree
    cmd = "checkm tree -q -t {0} -x fa {1} {2}".format(
        threads, fasta_dir, tmpdir
    )
    logger.info(cmd)
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()

    # Build taxonomy values of the bins
    cmd = "checkm tree_qa {0} -q -o 1 -f {1}".format(tmpdir, taxonomy_file)
    logger.info(cmd)
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()

    # Build lineage marker set
    markers_set = join(tmpdir, "markers.txt")
    cmd = "checkm lineage_set -q {0} {1}".format(tmpdir, markers_set)
    logger.info(cmd)
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()

    # Compute the analysis
    cmd = "checkm analyze -q -x fa -t {0} {1} {2} {3}".format(
        threads, markers_set, fasta_dir, tmpdir
    )
    logger.info(cmd)
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()

    # Write the summary file
    cmd = "checkm qa -q {0} {1} -o 2 > {2}".format(markers_set, tmpdir, outfile)
    logger.info(cmd)
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()


def checkm_compare_bins(
    overlapping_checkm_file,
    overlapping_taxonomy_file,
    recursive_checkm_file,
    recursive_taxonomy_file,
    prefix,
):
    """Compare the completness and contamination of the bins from the first step
    and from the recursive step. If the recursive step decrease the completion
    of one bin without decreasing its contamination, it will kept the bin before
    the recursive step. Moreover if the completion goes below 50% with the
    recursive it will not take the recursive bins.

    Parameters:
    -----------
    overlapping_checkm_file : str
        Path to the checkm summary from the overlapping step.
    overlapping_taxonomy_file : str
        path to the overlapping checkm taxonomy results file.
    recursive_checkm_file : str
        Path to the checkm summary from the recursive step.
    recursive_taxonomy_file : str
        path to the recursive checkm taxonomy results file.
    prefix : str
        Sample prefix to use.

    Returns:
    --------
    dict:
        Dictionnary with the informations of the final bins kept by MetaTOR.
    """

    # Load the checkm summary
    checkm_summary_overlapping = mio.read_results_checkm(
        overlapping_checkm_file, overlapping_taxonomy_file
    )
    checkm_summary_recursive = mio.read_results_checkm(
        recursive_checkm_file, recursive_taxonomy_file
    )

    # Prepare a dictionnary for a final summary.
    checkm_summary = dict()

    # Retrieve maximum completness of the recursive bins.
    for recursive_bin in checkm_summary_recursive:
        overlapping_bin = "_".join(
            [f"{prefix}", recursive_bin.split("_")[-2], "0"]
        )
        try:
            checkm_summary_overlapping[overlapping_bin][
                "max_rec_completness"
            ] = max(
                float(checkm_summary_recursive[recursive_bin]["completness"]),
                checkm_summary_overlapping[overlapping_bin][
                    "max_rec_completness"
                ],
            )
            checkm_summary_overlapping[overlapping_bin]["rec_id"].append(
                recursive_bin
            )
        except KeyError:
            checkm_summary_overlapping[overlapping_bin][
                "max_rec_completness"
            ] = float(checkm_summary_recursive[recursive_bin]["completness"])
            checkm_summary_overlapping[overlapping_bin]["rec_id"] = [
                recursive_bin
            ]

    # If there are some recursive bins which have not loose too much completion
    # write their information otherwise keep the original bins.
    for overlapping_bin in checkm_summary_overlapping:
        try:
            max_rec = checkm_summary_overlapping[overlapping_bin][
                "max_rec_completness"
            ]
            over = float(
                checkm_summary_overlapping[overlapping_bin]["completness"]
            )
            if (max_rec > (over / 2)) & (max_rec > 50):
                for rec_id in checkm_summary_overlapping[overlapping_bin][
                    "rec_id"
                ]:
                    checkm_summary[rec_id] = checkm_summary_recursive[rec_id]
            else:
                checkm_summary_overlapping[overlapping_bin].pop(
                    "max_rec_completness"
                )
                checkm_summary_overlapping[overlapping_bin].pop("rec_id")
                checkm_summary[overlapping_bin] = checkm_summary_overlapping[
                    overlapping_bin
                ]
        except KeyError:
            checkm_summary[overlapping_bin] = checkm_summary_overlapping[
                overlapping_bin
            ]

    return checkm_summary
