#!/usr/bin/env python3
# coding: utf-8

"""Validates bins using miComplete and make a recursive partition to try to
decontaminate them.

General functions to validate bins completion using miComplete and make
recursive iterations of Louvain or Leiden to try to partition contaminated bins.
Only bins with a weighted completeness of at least 33% and a contamination,
(weighted redundancy - 1) / weighted completeness, of at least 5% are subject to
the recursive step. Its bins replace the original bin when the most complete of
them keeps enough of its completeness, more than 0.4 + (c - 0.4) / 1.6 for a
bin of completeness c; otherwise the original bin is kept. The step repeats on
the new bins, 9 times at most.


Functions in this module:
    - bin_quality
    - check_partition_table
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
"""

import logging
import metator.clusters as mtc
import metator.io as mio
import metator.partition as mtp
import os
import pandas as pd
import pyfastx
import shutil
import subprocess as sp
import sys
import tempfile
from functools import partial
from importlib.metadata import version
from metator.log import logger
from os.path import join

# The fasta files of the bins that miComplete evaluates.
FASTA_EXTENSIONS = (".fa", ".fasta", ".fna")

#: The columns of the contig data table of a partition that the validation
#: reads (see check_partition_table).
PARTITION_COLUMNS = ("ID", "Name", "Size", "Hit", "Shotgun_coverage", "Overlapping_cluster_ID")


def check_partition_table(contig_data_file):
    """Refuse a contig data table that lacks a column the validation reads
    (PARTITION_COLUMNS): the error, logged, names them. Only its header is
    read: a table refused here costs no miComplete run.

    Parameters:
    -----------
    contig_data_file : str
        Path to the contig data table of a partition
        (contig_data_partition.txt).

    Raises:
    -------
    ValueError
        If the table lacks one of PARTITION_COLUMNS.
    """
    header = pd.read_csv(contig_data_file, sep="\t", nrows=0).columns
    missing = [column for column in PARTITION_COLUMNS if column not in header]
    if missing:
        logger.error(
            f"{contig_data_file} has no {', '.join(missing)} column: give the contig data table of a partition "
            "(contig_data_partition.txt)."
        )
        raise ValueError(f"{contig_data_file} lacks the columns {', '.join(missing)}.")


def bin_quality(completeness, redundancy):
    """The quality class of a bin, from miComplete's weighted completeness and
    weighted redundancy, its contamination being the redundancy minus 1:

        - "HQ" (high quality MAG): >= 90% completion ; <= 5% contamination
        - "MQ" (medium quality MAG): >= 70% completion ; <= 10% contamination
        - "LQ" (low quality MAG): >= 50% completion ; <= 10% contamination
        - "Contaminated": >= 50% completion ; > 10% contamination
        - "Other": < 50% completion

    The log of the validation counts the bins in these classes (see
    give_results_info).

    Parameters:
    -----------
    completeness : float
        Weighted completeness of the bin, from 0 to 1.
    redundancy : float
        Weighted redundancy of the bin, 1 for no duplicated marker.

    Returns:
    --------
    str:
        The class of the bin: "HQ", "MQ", "LQ", "Contaminated" or "Other".
    """
    completeness, redundancy = float(completeness), float(redundancy)
    if completeness < 0.5:
        return "Other"
    if redundancy > 1.1:
        return "Contaminated"
    if completeness >= 0.9 and redundancy <= 1.05:
        return "HQ"
    if completeness >= 0.7:
        return "MQ"
    return "LQ"


def correct_final_bin(contigs_data, final_fasta_dir, bin_summary):
    """Function to set the final bin of each contig from the fasta files of
    the final bins: the contigs of <final_fasta_dir>/<bin>.fa get the bin as
    Final_cluster.

    Parameters:
    -----------
    contigs_data : pandas.core.frame.DataFrame
        Dataframe with the contigs informations.
    final_fasta_dir : str
        Path to the directory of the fasta files of the final bins.
    bin_summary : dict
        Dictionnary with the informations of the final bins kept by MetaTOR.

    Returns:
    --------
    pandas.core.frame.DataFrame
        Dataframe with the contigs informations.
    """
    contigs_data = contigs_data.set_index("Name", drop=False)
    for bin_name in bin_summary:
        fasta_file = join(final_fasta_dir, f"{bin_name}.fa")
        # Read without an index, which would be left next to the bin.
        for name, _ in pyfastx.Fasta(fasta_file, build_index=False):
            contigs_data.loc[name, "Final_cluster"] = bin_name
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
        bin_name = contigs_data.loc[i, "Final_cluster"]
        if contigs_data.loc[contigs_data.index[0], "Shotgun_coverage"] != "-":
            total_sg_hit += contigs_data.loc[i, "Size"] * contigs_data.loc[i, "Shotgun_coverage"]
        total_hic_hit += contigs_data.loc[i, "Hit"]
        if bin_name != "ND":
            try:
                bin_summary[bin_name]["HiC_abundance"] += 100 * contigs_data.loc[i, "Hit"]
            except KeyError:
                bin_summary[bin_name]["HiC_abundance"] = 100 * contigs_data.loc[i, "Hit"]

            # If no depth files were given do not compute the Shotgun coverage.
            if contigs_data.loc[contigs_data.index[0], "Shotgun_coverage"] != "-":
                try:
                    bin_summary[bin_name]["SG_abundance"] += (
                        contigs_data.loc[i, "Size"] * contigs_data.loc[i, "Shotgun_coverage"] * 100
                    )
                except KeyError:
                    bin_summary[bin_name]["SG_abundance"] = (
                        contigs_data.loc[i, "Size"] * contigs_data.loc[i, "Shotgun_coverage"] * 100
                    )
    for bin_name in bin_summary:
        # A pair gives a hit to the contig of each of its reads: the hits of
        # the bin and the total both count read ends, so their ratio is
        # already the share of the Hi-C reads.
        bin_summary[bin_name]["HiC_abundance"] /= total_hic_hit
        if contigs_data.loc[contigs_data.index[0], "Shotgun_coverage"] != "-":
            bin_summary[bin_name]["SG_abundance"] /= total_sg_hit
    return bin_summary


def give_results_info(bin_summary):
    """Function to return the general information about the binning results:
    the number and total size of the bins of each quality class (see
    bin_quality), in the log.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the summary results of the kept bins.
    """
    # Count the bins and their size in each quality class.
    counts = {quality: [0, 0] for quality in ("HQ", "MQ", "LQ", "Contaminated", "Other")}
    for info in bin_summary.values():
        quality = bin_quality(info["Weighted completeness"], info["Weighted redundancy"])
        counts[quality][0] += 1
        counts[quality][1] += int(info["Length"])
    total = sum(count for count, _ in counts.values())
    total_size = sum(size for _, size in counts.values())

    # Return info in the logger:
    logger.info("{0} clusters have been kept after the recursive iterations.".format(total))
    logger.info("Total size of the extracted clusters: {0}".format(total_size))
    for quality, label in (
        ("HQ", "HQ MAGs"),
        ("MQ", "MQ MAGs"),
        ("LQ", "LQ MAGs"),
        ("Contaminated", "Contaminated potential MAGs"),
        ("Other", "Other clusters"),
    ):
        logger.info("{0}: {1}\tTotal Size: {2}".format(label, *counts[quality]))


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
    bact105 = pd.read_csv(out_bact105, sep="\t", comment="#", index_col=0).loc[:, :"CDs"]
    arch131 = pd.read_csv(out_arch131, sep="\t", comment="#", index_col=0).loc[:, :"CDs"]

    # Write header
    with open(outfile, "w") as out:
        out.write("## miComplete\n")
        out.write(f"## v{version('micomplete')}\n")
        out.write("## Weights: Bact105 and Arch131\n")
        out.write("Name\t{0}\tMarkers\n".format("\t".join(list(arch131.columns))))
        for bin_id in bact105.index:
            if bact105.loc[bin_id, "Weighted completeness"] >= arch131.loc[bin_id, "Weighted completeness"]:
                out.write("{0}\t{1}\tBacteria\n".format(bin_id, "\t".join(map(str, list(bact105.loc[bin_id])))))
            else:
                out.write("{0}\t{1}\tArchaea\n".format(bin_id, "\t".join(map(str, list(arch131.loc[bin_id])))))


def micomplete_compare_bins(
    recursive_micomplete_file,
    bin_summary,
    parent_dict,
    step,
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
        recursive_micomplete_file,
        sep="\t",
        comment="#",
        index_col=0,
    ).iloc[:, :13]

    # Group the recursive bins by parent bin, keeping the best completeness of
    # each family. Plain dicts: a list stored in a DataFrame cell is silently
    # turned into a string under pandas >= 3, which kept only the last child.
    children = {}
    max_rec_completeness = {}
    for recursive_bin in micomplete_recursive_summary.index:
        parent_bin = parent_dict[recursive_bin]
        children.setdefault(parent_bin, []).append(recursive_bin)
        completeness = float(micomplete_recursive_summary.loc[recursive_bin, "Weighted completeness"])
        max_rec_completeness[parent_bin] = max(completeness, max_rec_completeness.get(parent_bin, completeness))

    # If there are some recursive bins which have not loose too much completion
    # write their information otherwise keep the original bins.
    for overlapping_bin in list(bin_summary):
        if overlapping_bin not in children:
            bin_summary[overlapping_bin]["recursive"] = False
            continue
        max_rec = max_rec_completeness[overlapping_bin]
        over = float(bin_summary[overlapping_bin]["Weighted completeness"])
        if (max_rec - 0.4) > ((over - 0.4) / 1.6):
            for rec_id in children[overlapping_bin]:
                bin_summary[rec_id] = micomplete_recursive_summary.loc[rec_id, :"CDs"].to_dict()
                bin_summary[rec_id]["recursive"] = True
                bin_summary[rec_id]["step"] = step
                bin_summary[rec_id]["parent"] = overlapping_bin
            # A bin was replaced: run another recursive step.
            contamination = True
            bin_summary.pop(overlapping_bin)
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
    # miComplete writes its intermediate files to the current directory, named
    # after each fasta: run it in a private directory, removed afterwards, so
    # that runs started from the same directory cannot clobber each other.
    fasta_dir = os.path.abspath(fasta_dir)
    workdir = tempfile.mkdtemp(prefix="micomplete_", dir=os.path.dirname(os.path.abspath(outfile)))
    tmp_seq_tab = join(workdir, "micomplete_seq.tsv")
    with open(tmp_seq_tab, "w") as tab:
        for path in sorted(os.listdir(fasta_dir)):
            if path.endswith(FASTA_EXTENSIONS):
                tab.write(f"{join(fasta_dir, path)}\tfna\n")

    # Run miComplete (through metator._micomplete, which patches it for recent
    # Biopython and Python) with the bacterial and the archaeal marker sets.
    out_bact105 = join(workdir, "micomplete_bact105.tsv")
    out_arch131 = join(workdir, "micomplete_arch131.tsv")
    for markers, out in (("Bact105", out_bact105), ("Arch131", out_arch131)):
        cmd = [
            sys.executable, "-m", "metator._micomplete", tmp_seq_tab,
            "--hmms", markers,
            "--weights", markers,
            "--threads", str(threads),
            "--outfile", out,
        ]  # fmt: skip
        # In the process group of the run: a Ctrl-C, or the signal that a job
        # scheduler sends to the job, reaches it too.
        sp.run(cmd, cwd=workdir, check=True)

    # Merge bacteria and archaea
    merge_micomplete(out_bact105, out_arch131, outfile)
    shutil.rmtree(workdir)


def recursive_clustering(
    assembly,
    iterations,
    overlapping_parameter,
    resolution_parameter,
    outdir,
    recursive_fasta_dir,
    algorithm,
    bin_summary,
    contigs_data,
    network,
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
        Algorithm to use, either louvain or leiden.
    bin_summary : dict
        Dictionary containing information about the bins.
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs.
    network : pandas.DataFrame
        The network of the contigs of the map: one row per edge, with the
        columns contig1, contig2 and weight (see
        metator.network.map_network).
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
    dict
        Dictionnary with recursive bin_id as key and parent bin as values.
    """

    # Default no contamination
    contamination = False

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
            recursive_bins_iterations,
            threads,
        )

        # Defined overlapping bins according to the threshold
        recursive_bins = mtp.defined_overlapping_bins(
            overlapping_parameter,
            hamming_distance,
            recursive_core_bins,
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
            size,
            contamination,
            parent_dict,
            prefix,
        )

        logger.info("Recursive step for {0} is done.".format(bin_id))

    return contamination, contigs_data, parent_dict


def _bin_contigs(bin_id, contigs_data):
    """The names of the contigs of a bin (see bin_contig_ids)."""
    return list(contigs_data.loc[[contig_id - 1 for contig_id in bin_contig_ids(bin_id, contigs_data)], "Name"])


def bin_contig_ids(bin_id, contigs_data):
    """Function to list the IDs of the contigs of a bin from the contigs data:
    those of its overlapping bin and recursive bin, given by its name
    (<prefix>_<overlapping bin ID>_<recursive bin ID>).

    Parameters:
    -----------
    bin_id : str
        Name of the bin.
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs.

    Returns:
    --------
    list of int:
        IDs of the contigs of the bin.
    """
    over_bin_id = str(bin_id.split("_")[-2])
    rec_bin_id = str(bin_id.split("_")[-1])
    mask = (contigs_data["Overlapping_cluster_ID"] == over_bin_id) & (
        contigs_data["Recursive_cluster_ID"].apply(str) == rec_bin_id
    )
    return list(contigs_data.loc[mask, "ID"])


def recursive_clustering_worker(
    bin_id,
    network,
    algorithm,
    iterations,
    resolution_parameter,
    contigs_data,
):
    """Worker to partition one bin if it's contaminated: its subnetwork, the
    edges of the network between two of its contigs, in their order.

    Parameters:
    -----------
    bin_id : str
        Name of the bin.
    network : pandas.DataFrame
        Network of contigs from HiC librairies: one row per edge, with the
        columns contig1, contig2 and weight.
    algorithm : str
        Algorithm to use, either louvain or leiden.
    iterations : int
        Number of iterations to use for recursive iterations of Louvain or
        Leiden.
    resolution parameter : float
        Resolution parameter of Leiden algorithm.
    contigs_data : pandas.DataFrame
        The contig data of the partition, with the recursive bins so far.
    """
    logger.info("Cluster in progress: {0}".format(bin_id))

    # Extract contigs
    list_contigs = bin_contig_ids(bin_id, contigs_data)

    # Extract subnetwork: the edges between two contigs of the bin.
    subnetwork = network[network["contig1"].isin(list_contigs) & network["contig2"].isin(list_contigs)]

    # Partition the subnetwork.
    output_partition = mtp.algo_partition(
        algorithm,
        subnetwork,
        iterations,
        resolution_parameter,
    )
    logger.info("Output partition done for {0}.".format(bin_id))

    return output_partition


def recursive_decontamination(
    algorithm,
    assembly,
    contig_data_file,
    final_fasta_dir,
    input_fasta_dir,
    iterations,
    network,
    outdir,
    overlapping_parameter,
    recursive_fasta_dir,
    resolution_parameter,
    size,
    threads,
    prefix,
):
    """Function to validate bins do the recursive decontamination using Louvain
    or Leiden algorithm

    Parameters:
    -----------
    algorithm : str
        Algorithm to use to recursively partition the network. Either leiden
        or louvain.
    assembly : str
        Path to the assembly file used for the partition.
    contig_data_file : str
        Path to the contig data table to update.
    final_fasta_dir : str
        Path to write the final fasta decontaminated bins.
    input_fasta_dir : str
        Path to the directory where the fasta bin from the partition are.
    iterations : int
        Number of iterations to use for the recursive partition.
    network : pandas.DataFrame
        The network of the contigs of the map, in the IDs of the contigs in
        contig_data_file (see metator.network.relabel_network). The IDs of
        metator.network.map_network are the positions of the contigs in the
        map as it is read, which a run that regroups the map (validation,
        pipeline) moves: join by contig name.
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
    threads : int
        Number of threads to use.
    prefix : str
        Sample prefix to use.
    """

    # A table without the columns of a partition: refused before miComplete
    # runs on every cluster.
    check_partition_table(contig_data_file)

    # No bin to validate: miComplete would write no result to read.
    if not any(path.endswith(FASTA_EXTENSIONS) for path in os.listdir(input_fasta_dir)):
        logger.error(
            f"{input_fasta_dir} holds no cluster to validate (no .fa, .fasta or .fna file): partition writes none when "
            "no cluster reaches its minimum size (--size of partition or pipeline)."
        )
        raise ValueError

    # Defined miComplete output file path
    overlapping_micomplete_file = join(outdir, "overlapping_micomplete_results.txt")

    logger.info("Lauch miComplete quality check.")

    # Launch miComplete
    micomplete_quality(
        input_fasta_dir,
        overlapping_micomplete_file,
        threads,
    )

    # Load contigs data:
    contigs_data = pd.read_csv(
        contig_data_file,
        sep="\t",
        header=0,
        # As written: "00008", not 8.
        converters={"Overlapping_cluster_ID": str},
    )
    contigs_data["index"] = contigs_data["ID"] - 1
    contigs_data = contigs_data.set_index("index")

    # Add new coulumns for recursive information.
    contigs_data["Recursive_cluster_ID"] = f"{0:05d}"
    contigs_data["Recursive_cluster_contigs"] = pd.Series("-", index=contigs_data.index, dtype=object)
    contigs_data["Recursive_cluster_size"] = pd.Series("-", index=contigs_data.index, dtype=object)
    contigs_data["Final_cluster"] = "ND"

    # Load miComplete result:
    bin_summary = mio.micomplete_results_to_dict(overlapping_micomplete_file)

    # The clustering rounds of the run, for its map (see metator.clusters):
    # first the overlapping clusters that miComplete checked.
    rounds = mtc.Rounds()
    for bin_id, info in bin_summary.items():
        rounds.add(mtc.OVERLAPPING, bin_id, _bin_contigs(bin_id, contigs_data), columns=info)
    rounds.close_round(mtc.OVERLAPPING)

    # Add columns about recursive step.
    for bin_id in bin_summary:
        bin_summary[bin_id]["recursive"] = True
        bin_summary[bin_id]["step"] = 0
        bin_summary[bin_id]["parent"] = None

    # Recursively remove contamination.
    contamination = True
    step = 1

    logger.info("Starts recursive decontamition step:")

    while contamination:
        logger.info(f"------------------ Recursive step {step} started. --------------------- ")
        # Create fasta dir.
        recursive_fasta_dir_step = join(recursive_fasta_dir, f"step_{step}")
        os.makedirs(recursive_fasta_dir_step, exist_ok=True)

        # Stop to report info log
        logger.setLevel(logging.WARNING)

        # Iterates Louvain or Leiden on contaminated and complete bins.
        (
            contamination,
            contigs_data,
            parent_dict,
        ) = recursive_clustering(
            assembly,
            iterations,
            overlapping_parameter,
            resolution_parameter,
            outdir,
            recursive_fasta_dir_step,
            algorithm,
            bin_summary,
            contigs_data,
            network,
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
            recursive_micomplete_file = join(outdir, f"recursive_micomplete_step_{step}.txt")

            micomplete_quality(
                recursive_fasta_dir_step,
                recursive_micomplete_file,
                threads,
            )

            # The clusters of this step, with their contigs now: a later step
            # gives the contigs of a cluster it splits other IDs.
            step_round = mtc.recursive_round(step)
            step_quality = mio.micomplete_results_to_dict(recursive_micomplete_file)
            for bin_id, parent in parent_dict.items():
                rounds.add(
                    step_round,
                    bin_id,
                    _bin_contigs(bin_id, contigs_data),
                    parent=parent,
                    columns=step_quality.get(bin_id),
                )
            rounds.close_round(step_round)

            # Compare
            bin_summary, contamination = micomplete_compare_bins(
                recursive_micomplete_file,
                bin_summary,
                parent_dict,
                step,
            )

        # Keep overlapping bin information
        else:
            logger.info(f"No more contaminated cluster has been found after {step} steps.")
        # Increase count
        logger.info(f"Recursive step {step} done.")
        step += 1
        if step == 10:
            contamination = False

    # Whether each final bin comes from the recursive partition: the flag
    # marked the bins to examine at the next step.
    for info in bin_summary.values():
        info["recursive"] = info["step"] > 0

    # Create fasta directory and copy final bins.
    if not os.path.exists(final_fasta_dir):
        os.makedirs(final_fasta_dir)

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
    contigs_data = write_bins_contigs(bin_summary, contigs_data, binning_file)

    # Compute the abundance of the mags.
    bin_summary = get_bin_coverage(bin_summary, contigs_data)

    # The final clusters, each the cluster of the same name of the round it
    # comes from, with the columns of the bin summary (but its parent, which
    # the round says).
    final_contigs = contigs_data.groupby("Final_cluster")["Name"]
    for bin_name, info in bin_summary.items():
        rounds.add(
            mtc.FINAL,
            bin_name,
            final_contigs.get_group(bin_name) if bin_name in final_contigs.groups else [],
            parent=bin_name,
            columns={key: value for key, value in info.items() if key not in ("parent", "recursive")},
        )
    rounds.write(outdir)

    # Save bin information in final file
    bin_summary_file = join(outdir, "cluster_summary.txt")
    mio.write_bin_summary(bin_summary, bin_summary_file)

    # Write the new file
    contig_data_file_final = join(outdir, "contig_data_final.txt")
    contigs_data.to_csv(contig_data_file_final, sep="\t", header=True, index=True)


def update_contigs_data_recursive(
    bin_id,
    contigs_data,
    recursive_bins,
    assembly,
    outdir,
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
    # Indexed (<assembly>.fxi) only when a bin is extracted from it.
    fasta = None

    # Extract last recursive ID.
    over_id = contigs_data.loc[recursive_bins[1][0] - 1, "Overlapping_cluster_ID"]
    max_rec_id = max(
        map(
            int,
            contigs_data[contigs_data.Overlapping_cluster_ID == over_id]["Recursive_cluster_ID"],
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
            contigs_data.loc[recursive_bin, "Recursive_cluster_ID"] = f"{rec_final_id:05d}"
            contigs_data.loc[recursive_bin, "Recursive_cluster_contigs"] = recursive_bin_contigs_number
            contigs_data.loc[recursive_bin, "Recursive_cluster_size"] = recursive_bin_length

            # Defined name of the recursive bin
            oc_id = int(contigs_data.loc[recursive_bin[0], "Overlapping_cluster_ID"])
            output_file = join(outdir, f"{prefix}_{oc_id:05d}_{rec_final_id:05d}.fa")

            # Update bin_summary:
            parent_dict[f"{prefix}_{oc_id:05d}_{rec_final_id:05d}"] = bin_id

            # Retrieve names of the contigs
            list_contigs = list(contigs_data.loc[recursive_bin, "Name"])

            # Generate the fasta
            if fasta is None:
                fasta = pyfastx.Fasta(assembly)
            mio.write_contigs(fasta, list_contigs, output_file)

            # Add one to the recursive id
            rec_id += 1

    return contamination, contigs_data, parent_dict


def write_bins_contigs(bin_summary, contigs_data, outfile):
    """Function to write a table with the contigs of the final bins and their
    bin, without header: the binning that anvio imports.

    The bin of each contig is the Final_cluster column of the contigs data, which
    correct_final_bin sets from the fasta files of the final bins. Deriving it
    from the IDs of the overlapping and recursive bins of the contigs left out
    the contigs of a bin that was partitioned again but kept whole.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the informations of the final bins kept by MetaTOR.
    contigs_data : pandas.core.frame.DataFrame
        Dataframe with the contigs informations, with their final bin.
    outfile : str
        Path where to write the output file.

    Returns:
    --------
    pandas.core.frame.DataFrame
        Dataframe with the contigs informations, as given.
    """
    binned = contigs_data[contigs_data["Final_cluster"].isin(list(bin_summary))]
    binned[["Name", "Final_cluster"]].to_csv(outfile, sep="\t", header=False, index=False)
    return contigs_data
