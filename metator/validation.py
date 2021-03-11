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


def checkm(fasta_dir, outfile, tmpdir, threads):
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
    cmd = "checkm qa -q {0} {1} -o 1 > {2}".format(markers_set, tmpdir, outfile)
    logger.info(cmd)
    process = sp.Popen(cmd, shell=True)
    out, err = process.communicate()


# TODO
def compare_bins(overlapping_checkm_file, recursif_checkm_file):
    """Compare the completness and contamination of the bins and choose which
    are the most relevant bins.

    Parameters:
    -----------

    Recursif:
    ---------


    """
    return 0


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
    """Function to run recursive iterations on contaminated bins."""

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
            float(checkm_summary[bin_id]["contamination"]) >= 10
        ):

            # Add a boolean to say that there is at least one bin contaminated
            contamination = True
            logger.info("Bin in progress: {0}".format(bin_id))
            subnetwork_file = join(tmpdir, "subnetwork_" + bin_id + ".txt")
            bin_id = bin_id.split("_")[1]

            # Extract contigs
            list_contigs = list(contigs_data["id"][contigs_data["oc_id"] == bin_id])
            
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
                recursif_bins,
                recursif_bins_iterations,
            ) = mtp.detect_core_bins(output_louvain, iterations)

            # update bin data
            contigs_data = update_contigs_data(
                contigs_data, recursif_bins, assembly, outdir, size
            )

    # Write the new file
    contig_data_file_2 = join(outdir, "contig_data_rec.txt")
    contigs_data.to_csv(contig_data_file_2, sep="\t", header=None, index=False)

    return contamination


def update_contigs_data(contigs_data, recursif_bins, assembly, outdir, size):
    """Update the data of the bin according to the recursif step and generated
    their fasta.

    Paramaters:
    -----------

    Returns:
    --------
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
            rec_id += 1

            if recursif_bin_length > size:

                # Defined name of the recursif bin
                oc_id = contigs_data.iloc[recursif_bin[0], 9]
                output_file = join(
                    outdir, "MetaTOR_{0}_{1}.fa".format(oc_id, rec_id)
                )

                # Retrieve names of the contigs
                list_contigs = list(contigs_data.iloc[recursif_bin, 1])

                # Generate the fasta
                mtp.extract_contigs(assembly, list_contigs, output_file)

    return contigs_data
