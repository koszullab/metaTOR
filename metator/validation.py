#!/usr/bin/env python3
# coding: utf-8

"""Validation of the bins and recursive used of Louvain to try to removed
contaminated bins.

General functions to validate bins completion using checkM and make recursive
iterations of Louvain to try to partition contaminated bins.

Functions in this module:
    - checkM
    - compare_bins
    - louvain_recursif
    - update_recursif_louvain
"""


import metator.io as mio
import metator.partition as mtp
import networkx as nx
import pandas as pd
import subprocess as sp
import sys
from metator.log import logger
from os.path import join


def checkm(fasta_dir, outfile, outdir, tmpdir, threads):
    """Function to evaluate fasta bins using CheckM. Write a result summary in a
    text file.

    Parameters:
    -----------
    fasta_dir : str
        Path to the input fasta of the bins to evaluate.
    outfile : str
        Path to the file where the results of checkm will be written.
    tmpdir : str
        Path to the temporary directory where CheckM intermediary files will be
        written.
    threads : int
        Numbers of threads to use for CheckM
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
    taxonomy_file = join(outdir, "taxonomy.txt")
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


def compare_bins(overlapping_checkm_file, recursif_checkm_file):
    """Compare the completness and contamination of the bins and choose which
    are the most relevant bins.

    Parameters:
    -----------
    overlapping_checkm_file : str
        Path to the checkm summary from the overlapping step.
    recursif_checkm_file : str
        Path to the checkm summary from the recurisf step.

    Returns:
    --------
    dict:
        Dictionnary with the informations of the final bins kept by MetaTOR.
    """

    # Load the checkm summary
    checkm_summary_overlapping = mio.read_results_checkm(
        overlapping_checkm_file
    )
    checkm_summary_recursif = mio.read_results_checkm(recursif_checkm_file)

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
    outdir,
    louvain,
    tmpdir,
    checkm_file,
    contigs_data_file,
    network_file,
    size,
):
    """Function to run recursive iterations on contaminated bins in order to try
    to improve the quality of the bins using Louvain algorthm.

    Parameters:
    -----------
    assembly : str
        Path to the fasta file used as assembly.
    iterations : int
        Number of iterations to use for recursive itarations of Louvain.
    outdir : str
        Path to the output directory.
    louvain : str
        Path to the directory to found louvain functions if you want to use the
        cpp version instead of the python one.
    tmpdir : str
        Path the temp directory.
    checkm_file : str
        Path to the output file of checkm from checkm function.
    contigs_data_file : str
        Path to the contigs data file from metator partition.
    network_file : str
        Path to the network file from metator network.
    size : int
        Size threshodl in base pairs of the bins.

    Return:
    boolean:
        True if at least one new bin has been generated.
    """

    threads = 1

    # Load CheckM result:
    checkm_summary = mio.read_results_checkm(checkm_file)

    # Load network:
    network = nx.read_edgelist(
        network_file, nodetype=int, data=(("weight", float),)
    )

    # Load contigs data:
    contigs_data = pd.read_csv(
        contigs_data_file, sep="\t", header=None, index_col=False
    )
    contigs_data.columns = [
        "id",
        "name",
        "length",
        "GC_content",
        "hit",
        "coverage",
        "cc_id",
        "cc_nb",
        "cc_length",
        "oc_id",
        "oc_nb",
        "oc_length",
    ]

    # Add new coulumns for recursive information.
    contigs_data["rec_id"] = "0"
    contigs_data["rec_nb"] = "-"
    contigs_data["rec_length"] = "-"

    # Default no contamination
    contamination = False

    # Iterate on chcekm summary to find conatminated bins:
    for bin_id in checkm_summary:
        if (float(checkm_summary[bin_id]["completness"]) >= 50) & (
            float(checkm_summary[bin_id]["contamination"]) >= 5
        ):

            logger.info("Bin in progress: {0}".format(bin_id))
            subnetwork_file = join(tmpdir, "subnetwork_" + bin_id + ".txt")
            bin_id = bin_id.split("_")[1]

            # Extract contigs
            list_contigs = list(
                contigs_data["id"][contigs_data["oc_id"] == bin_id]
            )

            # Extract subnetwork
            subnetwork = network.subgraph(list_contigs)

            # Write the new subnetwork
            nx.write_edgelist(
                subnetwork, subnetwork_file, delimiter="\t", data=["weight"]
            )

            # Use Louvain algorithmon the subnetwork.
            if louvain == "None":
                output_louvain = mtp.louvain_iterations_py(
                    subnetwork_file,
                    iterations,
                )
            else:
                output_louvain = mtp.louvain_iterations_cpp(
                    subnetwork_file,
                    iterations,
                    tmpdir,
                    louvain,
                )

            # Detect core bins
            logger.info("Detect recursive bins:")
            (
                recursif_core_bins,
                recursif_bins_iterations,
            ) = mtp.detect_core_bins(output_louvain, iterations)

            # Compute the Hamming distance between core bins.
            logger.info("Detect overlapping bins:")
            hamming_distance = mtp.hamming_distance(
                recursif_bins_iterations,
                iterations,
                threads,
            )

            # Defined overlapping bins according to the threshold
            recursif_bins = mtp.defined_overlapping_bins(
                0.9,
                hamming_distance,
                recursif_core_bins,
                recursif_bins_iterations,
            )

            # update bin data
            contamination, contigs_data = update_contigs_data_recursif(
                contigs_data,
                recursif_bins,
                assembly,
                outdir,
                tmpdir,
                size,
                contamination,
            )

    # Write the new file
    contig_data_file_2 = join(outdir, "contig_data_recursif.txt")
    contigs_data.to_csv(contig_data_file_2, sep="\t", header=None, index=False)

    return contamination, contigs_data


def update_contigs_data_recursif(
    contigs_data, recursif_bins, assembly, outdir, tmpdir, size, contamination
):
    """Update the data of the bin according to the recursif step and generated
    their fasta.

    Parameters:
    -----------
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs.
    recursif_bins : dict
        Dictionnary which has  as keys the values of the iterations from Louvain
        separated by a semicolon and as values the list of the id of the
        contigs.
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
        Louvain separated by a semicolon and as values the list of the id of the
        contigs.
    """

    # Add recursif bin information
    rec_id = 1
    for i in recursif_bins:
        # Extract contigs of the bin
        recursif_bin = [id - 1 for id in recursif_bins[i]]
        recursif_bin_data = contigs_data.iloc[recursif_bin]
        recursif_bin_contigs_number = len(recursif_bin)
        recursif_bin_length = sum(recursif_bin_data.length)

        if recursif_bin_contigs_number > 1:
            # Write the new information
            contigs_data.iloc[recursif_bin, 12] = rec_id
            contigs_data.iloc[recursif_bin, 13] = recursif_bin_contigs_number
            contigs_data.iloc[recursif_bin, 14] = recursif_bin_length

            if recursif_bin_length > size:

                # If one new bin is generated change the boolean value to True
                contamination = True

                # Defined name of the recursif bin
                oc_id = contigs_data.iloc[recursif_bin[0], 9]
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
    list_bin_id = []
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
                ]["oc_id"]
            )
            rec_id = str(
                contigs_data.iloc[
                    i,
                ]["rec_id"]
            )
            try:
                rec_ids = list_bin_id[over_id]
                if rec_id in rec_ids:
                    f.write(
                        "{0}\tMetaTOR_{1}_{2}\n".format(
                            contigs_data.iloc[
                                i,
                            ]["name"],
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
                            ]["name"],
                            over_id,
                            rec_id,
                        )
                    )
            except KeyError:
                pass
