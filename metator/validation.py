#!/usr/bin/env python3
# coding: utf-8

"""Validates bins using CheckM and make a recursive partition to try to
decontaminate them.

General functions to validate bins completion using checkM and make recursive
iterations of Louvain or Leiden to try to partition contaminated bins. Only bins
with more than 50% completion and 5% contamination are subject to the recursive
step. If the recursive step gave worst results than the first (decrease of the
completion with no decrease of the contamination), it will keep the original
bin.


Functions in this module:
    - checkM
    - compare_bins
    - louvain_recursif
    - give_results_info
    - recursive_decontamination
    - update_recursif_louvain
    - write_bin_contigs
"""

import logging
import metator.io as mio
import metator.partition as mtp
import networkx as nx
import os
import pandas as pd
import subprocess as sp
from metator.log import logger
from os.path import join


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


def compare_bins(
    overlapping_checkm_file,
    overlapping_taxonomy_file,
    recursif_checkm_file,
    recursif_taxonomy_file,
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
    recursif_checkm_file : str
        Path to the checkm summary from the recurisf step.
    recursif_taxonomy_file : str
        path to the recursif checkm taxonomy results file.

    Returns:
    --------
    dict:
        Dictionnary with the informations of the final bins kept by MetaTOR.
    """

    # Load the checkm summary
    checkm_summary_overlapping = mio.read_results_checkm(
        overlapping_checkm_file, overlapping_taxonomy_file
    )
    checkm_summary_recursif = mio.read_results_checkm(
        recursif_checkm_file, recursif_taxonomy_file
    )

    # Prepare a dictionnary for a final summary.
    checkm_summary = dict()

    # Retrieve maximum completness of the recursif bins.
    for recursive_bin in checkm_summary_recursif:
        overlapping_bin = "_".join(
            ["MetaTOR", recursive_bin.split("_")[1], "0"]
        )
        try:
            checkm_summary_overlapping[overlapping_bin][
                "max_rec_completness"
            ] = max(
                float(checkm_summary_recursif[recursive_bin]["completness"]),
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
            ] = float(checkm_summary_recursif[recursive_bin]["completness"])
            checkm_summary_overlapping[overlapping_bin]["rec_id"] = [
                recursive_bin
            ]

    # If there are some recursif bins which have not loose too much completion
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
                    checkm_summary[rec_id] = checkm_summary_recursif[rec_id]
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


def louvain_recursif(
    assembly,
    iterations,
    overlapping_parameter,
    resolution_parameter,
    outdir,
    recursive_fasta_dir,
    algorithm,
    tmpdir,
    checkm_file,
    taxonomy_file,
    contigs_data_file,
    network_file,
    size,
    threads,
):
    """Function to run recursive iterations on contaminated bins in order to try
    to improve the quality of the bins using Louvain or Leiden algorthm.

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
        Algorithm to use, either louvain or leiden.
    tmpdir : str
        Path the temp directory.
    checkm_file : str
        Path to the output file of CheckM from checkm function.
    taxonomy_file : str
        Path to the taxonomy CheckM file.
    contigs_data_file : str
        Path to the contigs data file from metator partition.
    network_file : str
        Path to the network file from metator network.
    size : int
        Size threshodl in base pairs of the bins.
    threads : int
        Number of threads to use.

    Return:
    boolean:
        True if at least one new recursive bin has been generated.
    """

    # Load CheckM result:
    checkm_summary = mio.read_results_checkm(checkm_file, taxonomy_file)

    # Load network:
    network = nx.read_edgelist(
        network_file, nodetype=int, data=(("weight", float),)
    )

    # Load contigs data:
    contigs_data = pd.read_csv(
        contigs_data_file, sep="\t", header=0, index_col=False
    )

    # Add new coulumns for recursive information.
    contigs_data["Recursive_bin_ID"] = "0"
    contigs_data["Recursive_bin_contigs"] = "-"
    contigs_data["Recursive_bin_length"] = "-"

    # Default no contamination
    contamination = False

    # Iterate on chcekm summary to find conatminated bins:
    for bin_id in checkm_summary:
        if (float(checkm_summary[bin_id]["completness"]) >= 50) & (
            float(checkm_summary[bin_id]["contamination"]) >= 5
        ):

            logger.info("Bin in progress: {0}".format(bin_id))
            subnetwork_file = join(tmpdir, "subnetwork_" + bin_id + ".txt")
            bin_id = int(bin_id.split("_")[1])

            # Extract contigs
            list_contigs = list(
                contigs_data["ID"][contigs_data["Overlapping_bin_ID"] == bin_id]
            )

            # Extract subnetwork
            subnetwork = network.subgraph(list_contigs)

            # Write the new subnetwork
            nx.write_edgelist(
                subnetwork, subnetwork_file, delimiter="\t", data=["weight"]
            )

            # Stop to report info log
            # logger.setLevel(logging.WARNING)

            # Use Louvain algorithmon the subnetwork.
            if algorithm == "leiden":
                LEIDEN_PATH = os.environ["LEIDEN_PATH"]
                output_partition = mtp.leiden_iterations_java(
                    subnetwork_file,
                    iterations,
                    resolution_parameter,
                    tmpdir,
                    LEIDEN_PATH,
                )
            elif algorithm == "louvain":
                LOUVAIN_PATH = os.environ["LOUVAIN_PATH"]
                output_partition = mtp.louvain_iterations_cpp(
                    subnetwork_file,
                    iterations,
                    tmpdir,
                    LOUVAIN_PATH,
                )
            else:
                logger.error('algorithm should be either "louvain" or "leiden"')
                raise ValueError

            # Detect core bins
            (
                recursif_core_bins,
                recursif_bins_iterations,
            ) = mtp.detect_core_bins(output_partition, iterations)

            # Compute the Hamming distance between core bins.
            hamming_distance = mtp.get_hamming_distance(
                recursif_bins_iterations,
                iterations,
                threads,
            )

            # Defined overlapping bins according to the threshold
            recursif_bins = mtp.defined_overlapping_bins(
                overlapping_parameter,
                hamming_distance,
                recursif_core_bins,
                recursif_bins_iterations,
            )

            # update bin data and generate fasta
            contamination, contigs_data = update_contigs_data_recursif(
                contigs_data,
                recursif_bins,
                assembly,
                recursive_fasta_dir,
                tmpdir,
                size,
                contamination,
            )

            # Put back the info log
            logger.setLevel(logging.INFO)

    # Write the new file
    contig_data_file_2 = join(outdir, "contig_data_recursif.txt")
    contigs_data.to_csv(contig_data_file_2, sep="\t", header=True, index=False)

    return contamination, contigs_data


def give_results_info(bin_summary):
    """Function to return the general information about the binning results.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the summary results of the kept bins.
    """

    # Defined categories of the bins
    HQ = 0  # Completness >= 90 and Contamination <= 5
    MQ = 0  # Completness >= 70 and Contamination <= 10
    LQ = 0  # Completness >= 50 and Contamination <= 10
    conta_bins = 0  # Completness >= 50 and Contamination > 10
    others = 0
    total_size = 0  # Total size of the extracted bins.

    # Class each bin in a category
    for bin_name in bin_summary:
        completness = float(bin_summary[bin_name]["completness"])
        contamination = float(bin_summary[bin_name]["contamination"])
        size = int(bin_summary[bin_name]["size"])
        if completness >= 50:
            if contamination > 10:
                conta_bins += 1
            else:
                if completness >= 90 and contamination <= 5:
                    HQ += 1
                elif completness >= 70:
                    MQ += 1
                else:
                    LQ += 1
        else:
            others += 1
        total_size += size
    total = HQ + MQ + LQ + conta_bins + others

    # Return info in the logger:
    logger.info(
        "{0} bins have been kept after the recursive iterations.".format(total)
    )
    logger.info("Total size of the extracted bins: {0}".format(total_size))
    logger.info("HQ MAGs: {0}".format(HQ))
    logger.info("MQ MAGs: {0}".format(MQ))
    logger.info("LQ MAGs: {0}".format(LQ))
    logger.info("Contaminated potential MAGs: {0}".format(conta_bins))
    logger.info("Others bins: {0}".format(others))


def recursive_decontamination(
    algorithm,
    assembly,
    contig_data_file,
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
):
    """Function to validate bins do the recursive decontamination using Louvain
    or Leiden algorithm

    Parameters:
    -----------
    algorithm : str
        Algorithm to use to recursively partition the network. Either leiden or
        louvain.
    assembly : str
        Path to the assembly file used for the partition.
    contig_data_file : str
        Path to the contig data table to update.
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
    """

    # Defined checkm output file path
    overlapping_checkm_file = join(outdir, "overlapping_checkm_results.txt")
    overlapping_taxonomy_file = join(outdir, "overlapping_checkm_taxonomy.txt")
    recursif_checkm_file = join(outdir, "recursif_checkm_results.txt")
    recursif_taxonomy_file = join(outdir, "recursif_checkm_taxonomy.txt")

    # Launch checkM
    checkm(
        input_fasta_dir,
        overlapping_checkm_file,
        overlapping_taxonomy_file,
        temp_directory,
        threads,
    )

    # Iterates Louvain on contaminated and complete bins
    contamination, contigs_data = louvain_recursif(
        assembly,
        iterations,
        overlapping_parameter,
        resolution_parameter,
        outdir,
        recursive_fasta_dir,
        algorithm,
        temp_directory,
        overlapping_checkm_file,
        overlapping_taxonomy_file,
        contig_data_file,
        network_file,
        size,
        threads,
    )

    # Recursive iterations of Louvain on the contaminated bins. Save bin
    # information if the new bins have the same quality otherwise keep the
    # original bin information.
    if contamination:

        # Run checkm on the recursif bins.
        temp_directory = join(temp_directory, "checkm2")
        checkm(
            recursive_fasta_dir,
            recursif_checkm_file,
            recursif_taxonomy_file,
            temp_directory,
            threads,
        )

        # Compare
        bin_summary = compare_bins(
            overlapping_checkm_file,
            overlapping_taxonomy_file,
            recursif_checkm_file,
            recursif_taxonomy_file,
        )

    # Keep overlapping bin information
    else:
        logger.info("No contaminated bin have been found")
        bin_summary = mio.read_results_checkm(
            overlapping_checkm_file, overlapping_taxonomy_file
        )

    # Retrurn some values of efficiency of the binning.
    give_results_info(bin_summary)

    # Save bin information in final file
    bin_summary_file = join(outdir, "bin_summary.txt")
    mio.write_checkm_summary(bin_summary, bin_summary_file)

    # Write relevant bins/contigs information for anvio.
    binning_file = join(outdir, "binning.txt")
    write_bins_contigs(bin_summary, contigs_data, binning_file)


def update_contigs_data_recursif(
    contigs_data, recursif_bins, assembly, outdir, tmpdir, size, contamination
):
    """Update the data of the bin according to the recursive step and generated
    their fasta.

    Parameters:
    -----------
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs.
    recursif_bins : dict
        Dictionnary which has  as keys the values of the recursive iterations
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

    Returns:
    --------
    boolean:
        True if one bin has already been generated, false otherwise.
    pandas.DataFrame
        Updated dictionnary which has as keys the values of the iterations from
        the recursive partition separated by a semicolon and as values the list
        of the id of the contigs.
    """

    # Add recursif bin information
    rec_id = 1
    for i in recursif_bins:
        # Extract contigs of the bin
        recursif_bin = [id - 1 for id in recursif_bins[i]]
        recursif_bin_data = contigs_data.iloc[recursif_bin]
        recursif_bin_contigs_number = len(recursif_bin)
        recursif_bin_length = sum(recursif_bin_data.Size)

        if recursif_bin_contigs_number > 1:
            # Write the new information
            contigs_data.iloc[recursif_bin, 13] = rec_id
            contigs_data.iloc[recursif_bin, 14] = recursif_bin_contigs_number
            contigs_data.iloc[recursif_bin, 15] = recursif_bin_length

            if recursif_bin_length > size:

                # If one new bin is generated change the boolean value to True
                contamination = True

                # Defined name of the recursif bin
                oc_id = contigs_data.iloc[recursif_bin[0], 10]
                output_file = join(
                    outdir, "MetaTOR_{0}_{1}.fa".format(oc_id, rec_id)
                )

                # Retrieve names of the contigs
                list_contigs = list(contigs_data.iloc[recursif_bin, 1])

                # Generate the fasta
                contigs_file = join(
                    tmpdir, "MetaTOR_{0}_{1}.txt".format(oc_id, rec_id)
                )
                with open(contigs_file, "w") as f:
                    for contig in list_contigs:
                        f.write("{0}\n".format(contig))
                cmd = "pyfastx extract {0} -l {1} > {2}".format(
                    assembly, contigs_file, output_file
                )
                process = sp.Popen(cmd, shell=True)
                process.communicate()
                rec_id += 1

    return contamination, contigs_data


def write_bins_contigs(bin_summary, contigs_data, outfile):
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
    """

    # Create a list with the id of the bins
    list_bin_id = dict()
    for bin_name in bin_summary:
        over_id = bin_name.split("_")[1]
        rec_id = bin_name.split("_")[2]
        try:
            list_bin_id[over_id].append(rec_id)
        except KeyError:
            list_bin_id[over_id] = [rec_id]

    # Write the contigs id with their bins id in table file
    with open(outfile, "w") as f:
        for i in range(len(contigs_data)):
            over_id = str(
                contigs_data.iloc[
                    i,
                ]["Overlapping_bin_ID"]
            )
            rec_id = str(
                contigs_data.iloc[
                    i,
                ]["Recursive_bin_ID"]
            )
            try:
                rec_ids = list_bin_id[over_id]
                if rec_id in rec_ids:
                    f.write(
                        "{0}\tMetaTOR_{1}_{2}\n".format(
                            contigs_data.iloc[
                                i,
                            ]["Name"],
                            over_id,
                            rec_id,
                        )
                    )
                # Case where the recursif bins where not kept.
                elif rec_ids == ["0"]:
                    rec_id = "0"
                    f.write(
                        "{0}\tMetaTOR_{1}_{2}\n".format(
                            contigs_data.iloc[
                                i,
                            ]["Name"],
                            over_id,
                            rec_id,
                        )
                    )
            except KeyError:
                pass
