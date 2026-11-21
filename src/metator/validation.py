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


The clusters that miComplete checks first are the overlapping clusters of the
partition that the map stores (see metator.metacool.read_partition) that
reach --size: their FASTA files, written by metacooler from the assembly in
the temporary directory of the run (see write_overlapping_fasta), as are
those of the clusters of each recursive step (see write_recursive_fasta). The
FASTA files of the final clusters are those of the clusters they are.

The tables of a validation (miComplete's results, the contig data and the
summary of the final clusters, the clusters of every round, the binning of
the contigs) are written in the temporary directory of the run, which -N
keeps: the map of the clusters stores what they hold (see
metator.metacool.write_map), and metator dump writes them from it.

Functions in this module:
    - bin_quality
    - correct_final_bin
    - get_bin_coverage
    - give_results_info
    - merge_micomplete
    - micomplete_compare_bins
    - micomplete_quality
    - name_recursive_clusters
    - recursive_clustering
    - recursive_clustering_worker
    - recursive_decontamination
    - write_overlapping_fasta
    - write_recursive_fasta
"""

import gzip
import logging
import metator.clusters as mtc
import metator.io as mio
import metator.partition as mtp
import os
import pandas as pd
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


def correct_final_bin(contigs_data, written, bin_summary):
    """Function to set the final bin of each contig from the contigs of the
    FASTA files of the final bins: the contigs written in <bin>.fa (written,
    the contigs of each FASTA file of the run, by cluster name: see
    write_overlapping_fasta and write_recursive_fasta) get the bin as
    Final_cluster.

    Parameters:
    -----------
    contigs_data : pandas.core.frame.DataFrame
        Dataframe with the contigs informations.
    written : dict
        The names of the contigs of the FASTA file of each cluster that the
        run checked, by cluster name.
    bin_summary : dict
        Dictionnary with the informations of the final bins kept by MetaTOR.

    Returns:
    --------
    pandas.core.frame.DataFrame
        Dataframe with the contigs informations.

    Raises:
    -------
    ValueError
        If a contig written in the FASTA file of a final bin is not in
        contigs_data.
    """
    from metator.metacool import _shown

    contigs_data = contigs_data.set_index("Name", drop=False)
    # A contig written that the table lacks would be added to it, with no
    # other value (a name that the table read otherwise: 00123 as 123).
    absent = [name for bin_name in bin_summary for name in written[bin_name] if name not in contigs_data.index]
    if absent:
        raise ValueError(
            f"{len(absent)} contig(s) of the FASTA files of the final clusters are not in the contig data of the "
            f"validation: {_shown(absent)}."
        )
    for bin_name in bin_summary:
        for name in written[bin_name]:
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


def write_overlapping_fasta(metacool, assembly, partition, contig_data, size, outdir, temp_dir=None):
    """Write the FASTA file of each overlapping cluster of a partition that
    reaches size bp (--size), for miComplete to check: <outdir>/<name>.fa,
    the records of its contigs as they are in the assembly, in the order of
    the partition (that of its core clusters, then of their contigs: see
    metator.clusters.Partition), through metacooler.write_fasta.

    Parameters:
    -----------
    metacool : str
        The map of the partition (a path, or an entry of a collection).
    assembly : str
        Path to the assembly that the map was made from.
    partition : metator.clusters.Partition
        The core and the overlapping clusters.
    contig_data : pandas.DataFrame
        The contig data of the contigs of the map, in its order (see
        metator.network.contig_table): their names and sizes.
    size : int
        The size in bp that an overlapping cluster reaches to be checked.
    outdir : str
        The directory to write the files in, made if missing.
    temp_dir : str, optional
        Where metacooler copies a compressed assembly, plain: the temporary
        directory of the run.

    Returns:
    --------
    dict:
        The names of the contigs of each cluster written, in the order of
        their records, by cluster name, the clusters by number.

    Raises:
    -------
    ValueError
        If the assembly lacks a contig of a cluster, or gives it another
        length than the map (the reason is logged as an error).
    """
    names = contig_data["Name"].to_numpy()
    sizes = partition.sizes(contig_data["Size"].to_numpy())
    kept = [number for number in partition.overlapping if sizes[number] >= size]
    groups = {partition.names[number]: [str(names[contig - 1]) for contig in partition.overlapping[number]] for number in kept}
    _write_fasta(metacool, assembly, outdir, groups, temp_dir, "The clusters to check")
    logger.info(f"{len(kept)} clusters have been extracted")
    logger.info(f"Total size of the extracted clusters: {round(sum(sizes[number] for number in kept) / 10**6, 3)}Mb")
    return groups


def plain_assembly(assembly, directory):
    """The assembly as the run hands it to metacooler.write_fasta: the file
    itself, when it is plain, or its copy decompressed in directory (the
    temporary directory of the run), when it is compressed with gzip, made
    once for the run: write_fasta would copy it again for every FASTA file
    write (the overlapping clusters, each recursive step). The copy is named
    after the assembly, without its .gz: an error of write_fasta names it."""
    if mio.file_compression(assembly) != "gz":
        return assembly
    name = os.path.basename(assembly)
    plain = join(directory, name[: -len(".gz")] if name.endswith(".gz") else f"{name}.fa")
    logger.info(f"Decompressing {assembly} into {plain}, once, for the FASTA files of the clusters.")
    with gzip.open(assembly, "rb") as source, open(plain, "wb") as copy:
        shutil.copyfileobj(source, copy)
    return plain


def write_recursive_fasta(metacool, assembly, groups, outdir, temp_dir=None):
    """Write the FASTA file of each cluster of a recursive step for miComplete
    to check: <outdir>/<name>.fa, the records of its contigs (groups: their
    names, by cluster name) as they are in the assembly, in their order,
    through metacooler.write_fasta, as write_overlapping_fasta writes those
    of the overlapping clusters.

    Raises:
    -------
    ValueError
        If the assembly lacks a contig of a cluster, or gives it another
        length than the map (the reason is logged as an error).
    """
    _write_fasta(metacool, assembly, outdir, groups, temp_dir, "The clusters of the recursive step")


def _write_fasta(metacool, assembly, outdir, groups, temp_dir, what):
    """Write the FASTA file of each group of contigs (by cluster name) in
    outdir, made if missing, with metacooler.write_fasta: a ValueError
    naming what the clusters are when it cannot (logged as an error)."""
    import metacooler

    os.makedirs(outdir, exist_ok=True)
    if not groups:
        return
    try:
        metacooler.write_fasta(metacool, assembly, outdir, groups=groups, force=True, temp_dir=temp_dir)
    except (KeyError, ValueError, OSError) as error:
        message = f"{what} cannot be written from {assembly}: {error}"
        logger.error(message)
        raise ValueError(message) from error


def recursive_clustering(
    iterations,
    overlapping_parameter,
    resolution_parameter,
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
    iterations : int
        Number of iterations to use for recursive iterations of Louvain or
        Leiden.
    overlapping_parameter : float
        Hamming distance threshold to consider two bins as the same bin.
    resolution parameter : float
        Resolution parameter of Leiden algorithm.
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
    dict
        The names of the contigs of each new recursive bin, in the order of
        the records of its FASTA file (see write_recursive_fasta), by bin
        name.
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
    groups = dict()

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

        # The new bins: their IDs in the contig data, and their contigs.
        contamination = (
            name_recursive_clusters(
                bin_id,
                contigs_data,
                recursive_bins,
                size,
                parent_dict,
                prefix,
                groups,
            )
            or contamination
        )

        logger.info("Recursive step for {0} is done.".format(bin_id))

    return contamination, contigs_data, parent_dict, groups


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
    partition,
    *,
    metacool,
    written,
    temp_dir=None,
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
        Path to the contig data table of the partition (see
        metator.partition.partition_table), to update.
    final_fasta_dir : str
        Path to write the final fasta decontaminated bins.
    input_fasta_dir : str
        Path to the directory of the FASTA files of the overlapping clusters
        to check (see write_overlapping_fasta).
    iterations : int
        Number of iterations to use for the recursive partition.
    network : pandas.DataFrame
        The network of the contigs of the map, in the IDs of the contigs in
        contig_data_file: their positions in the map as the run read it (see
        metator.network.map_network), which the contig data of the partition
        numbers them by too (see metator.partition.partition_table). A run
        that regroups the map (validation, pipeline) moves its contigs: join
        the network of the map it wrote to these tables by contig name (see
        metator.network.relabel_network).
    outdir : str
        Path to the directory where to write the tables of the validation:
        the temporary directory of the run (see the module's docstring).
    overlapping_parameter : int
        Hamming distance threshold in percentage to use to consider to bins as
        one in the recursive partition.
    recursive_fasta_dir : str
        Path to write the fasta decontaminated bins, one directory per
        recursive step (step_<k>): in the temporary directory of the run.
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
    partition : metator.clusters.Partition
        The core and the overlapping clusters of the partition, recorded as
        the first rounds of the run (see metator.clusters.record_partition):
        every cluster, those that miComplete checked with its columns.
    metacool : str
        The map of the run, which metacooler writes the FASTA files of the
        recursive clusters from the assembly with (see
        write_recursive_fasta).
    written : dict
        The names of the contigs of the FASTA file of each overlapping
        cluster in input_fasta_dir, by cluster name (see
        write_overlapping_fasta).
    temp_dir : str, optional
        Where metacooler copies a compressed assembly, plain: the temporary
        directory of the run.

    Returns:
    --------
    pandas.DataFrame:
        The contig data of the validation, with the final cluster of each
        contig (Final_cluster, ND for none), as contig_data_final.txt holds
        it.
    metator.io.SummaryLayout:
        How the summary of the final clusters (cluster_summary.txt) is
        written: its rows, its columns and their format.
    """

    # No bin to validate: miComplete would write no result to read.
    if not any(path.endswith(FASTA_EXTENSIONS) for path in os.listdir(input_fasta_dir)):
        logger.error(
            "No overlapping cluster of the partition reaches the minimum size of a cluster to check (--size of "
            "validation or pipeline): there is no cluster to validate."
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
        # Names as written: a contig named 00123 is not the number 123, nor
        # one named 1e5 a number at all.
        dtype={"Name": str},
        # As written: "00008", not 8.
        converters={"Overlapping_cluster_ID": str},
        # A contig named NA (or null, N/A...) is a name: only an empty field,
        # as pandas writes a missing value, is one.
        keep_default_na=False,
        na_values=[""],
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
    # first those of the partition, every core and overlapping cluster, the
    # overlapping clusters that miComplete checked with its columns.
    rounds = mtc.record_partition(mtc.Rounds(), partition, list(contigs_data.sort_index()["Name"]), quality=bin_summary)

    # The contigs of the FASTA file of each cluster checked, by name: those of
    # the final clusters are theirs.
    written = dict(written)

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
            groups,
        ) = recursive_clustering(
            iterations,
            overlapping_parameter,
            resolution_parameter,
            algorithm,
            bin_summary,
            contigs_data,
            network,
            size,
            threads,
            prefix,
        )

        # The FASTA files of the new bins, for miComplete.
        write_recursive_fasta(metacool, assembly, groups, recursive_fasta_dir_step, temp_dir)
        written.update(groups)

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
    contigs_data = correct_final_bin(contigs_data, written, bin_summary)

    # The contigs of the final bins and their bin, without header: the
    # binning that anvio imports.
    binned = contigs_data[contigs_data["Final_cluster"].isin(list(bin_summary))]
    binned[["Name", "Final_cluster"]].to_csv(join(outdir, "binning.txt"), sep="\t", header=False, index=False)

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
    summary = mio.write_bin_summary(bin_summary, join(outdir, "cluster_summary.txt"))

    # Write the new file
    contigs_data.to_csv(join(outdir, "contig_data_final.txt"), sep="\t", header=True, index=True)
    return contigs_data, summary


def name_recursive_clusters(bin_id, contigs_data, recursive_bins, size, parent_dict, prefix, groups):
    """Record the bins of the recursive partition of a contaminated bin that
    are larger than size: their IDs in the contig data, their parent, and
    their contigs, for their FASTA files (see write_recursive_fasta).

    Parameters:
    -----------
    bin_id : str
        Name of the parental bin.
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs, updated in place: the
        Recursive_cluster_ID, Recursive_cluster_contigs and
        Recursive_cluster_size of the contigs of each new bin.
    recursive_bins : dict
        Dictionnary which has as keys the values of the recursive iterations
        from Louvain or Leiden separated by a semicolon and as values the list
        of the id of the contigs.
    size : int
        Size threshold of a bin.
    parent_dict : dict
        Dictionnary with recursive bin_id as key and parent bin as values,
        updated in place.
    prefix : str
        Sample prefix to use.
    groups : dict
        The names of the contigs of each new bin, in the order of the
        partition, by bin name (<prefix>_<overlapping ID>_<recursive ID>),
        updated in place.

    Returns:
    --------
    boolean:
        True if one bin is larger than size.
    """
    contamination = False
    # Add recursive bin information
    rec_id = 1

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
            name = f"{prefix}_{oc_id:05d}_{rec_final_id:05d}"

            # Update bin_summary:
            parent_dict[name] = bin_id

            # The contigs of its FASTA file.
            groups[name] = [str(contig) for contig in contigs_data.loc[recursive_bin, "Name"]]

            # Add one to the recursive id
            rec_id += 1

    return contamination
