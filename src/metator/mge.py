#!/usr/bin/env python3
# coding: utf-8

"""Module with the MGE clustering functions.

It clusters the contigs of a .metacool contact map that are mobile genetic
elements (MGEs) from the Hi-C contacts between them: those that the map
holds between two contigs, and those that its contact statistics count
within each contig, between reads more than 1 kb apart.

Core functions to cluster the MGE contigs:
    - check_mge_map
    - generate_bin_summary
    - map_matrix
    - mge_binning
    - mge_contigs
    - resolve_matrix
    - shuffle_mge_bins
    - update_mge_data
    - write_mges_fasta
"""

import numpy as np
import pandas as pd
from metator.log import logger
from typing import List, Tuple

import metator.metacool as mtmc

# The contacts within a contig that mge counts: those of reads more than
# 1 kb apart, closer ones saying nothing about the element. metacooler's
# contact statistics count them, per contig, in the contigs table of a map
# that it made from the pairs (qc_intra_contig_long: reads at least
# distance_threshold apart), at 1001 bp by default from metacooler 0.4.0 on.
QC_SOURCE = "pairs"
DISTANCE_THRESHOLD = 1001
INTRA_LONG = "qc_intra_contig_long"

# The columns of the contigs table of the map that mge reads, but its
# contacts: those that metator network stores.
MGE_COLUMNS = ["GC_content", "Hit"]

# N between two contigs of an MGE cluster, in its record of mges_binned.fa.
SPACER = 200

# How to make a map that mge reads, said when it refuses one.
HOW_TO_MAKE_A_MAP = (
    "Make the map from the pairs files of the run, with the contact statistics of its reads at 1001 bp: hicstuff "
    "pipeline -M metacool, or metacooler create pairs, with metacooler 0.4.0 or later (metacooler create pairs "
    "--qc-threshold 1001 with an earlier default), then run metator network on it (or pipeline)."
)

# Where the network of a run is, said when mge refuses a map without it: in
# the map that the run wrote, which --output-map put elsewhere.
WHERE_THE_NETWORK_IS = (
    "run metator network on the map first (or pipeline), or give the map that such a run wrote (the --output-map "
    "of network, pipeline or validation)."
)


def check_mge_map(metacool, mc):
    """Refuse a map that mge cannot cluster the MGE contigs of, before any
    work.

    The contacts within an MGE contig are those of its reads more than 1 kb
    apart: the contact statistics of a map made from the pairs files of the
    run tell them (its QC block, measured from the pairs, at 1001 bp, with
    its counts per contig), which its pixels do not. The GC content and hits
    of each contig are those that metator network stores in the map. A map
    without clusters is fine: mge reads none of them.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    mc : metacooler.MetaCooler
        The map, opened.

    Raises
    ------
    metator.metacool.MetacoolError
        If the map holds no contact statistics of its pairs, has them
        measured from its pixels or at another distance than 1001 bp, or
        lacks the columns of metator network: the reason, and how to make a
        map that mge reads.
    """
    with mtmc.read_errors(metacool):
        qc = mc.metadata.get("qc")
        columns = set(mc.contigs.columns)
    if not isinstance(qc, dict) or INTRA_LONG not in columns:
        raise mtmc.MetacoolError(
            f"{metacool} holds no contact statistics of the reads of its contigs (metacooler's QC block, with "
            f"{INTRA_LONG} in its contigs table): metator mge counts the contacts within each MGE contig between reads "
            f"more than 1 kb apart, which the pixels of a map do not tell. {HOW_TO_MAKE_A_MAP}"
        )
    if qc.get("source") != QC_SOURCE:
        raise mtmc.MetacoolError(
            f"{metacool} has the contact statistics of its contigs measured from its pixels (QC source "
            f"{qc.get('source')!r}), not from the pairs: metator mge counts the contacts within each MGE contig "
            f"between reads more than 1 kb apart, which the pixels of a map do not tell. {HOW_TO_MAKE_A_MAP}"
        )
    if qc.get("distance_threshold") != DISTANCE_THRESHOLD:
        raise mtmc.MetacoolError(
            f"{metacool} splits the contacts within a contig at {qc.get('distance_threshold')} bp (its QC "
            f"distance_threshold): metator mge counts those between reads more than 1 kb apart, at "
            f"{DISTANCE_THRESHOLD} bp and more. {HOW_TO_MAKE_A_MAP}"
        )
    missing = [column for column in MGE_COLUMNS if column not in columns]
    if missing:
        raise mtmc.MetacoolError(
            f"{metacool} lacks the {' and '.join(missing)} column(s) of its contigs table that metator network "
            f"stores: {WHERE_THE_NETWORK_IS}"
        )


def mge_contigs(metacool, mc, mges_list):
    """The MGE contigs of a map, the longest first (ties by name), with the
    values of each that mge reads.

    That order, which resolve_matrix breaks the ties of its scores in and
    update_mge_data numbers the clusters of single contigs in, depends
    neither on the order of the map nor on its clusters: a map of the
    clusters, grouped by cluster, gives the MGE clusters, and their names,
    of the map of the contigs it was written from, longest first, as hicstuff
    and metacooler make it (which the contig data tables of a run followed).

    An MGE contig that the map lacks (a map made with a contig filter, or
    another assembly) is left out, with a warning that says how many and
    names the first.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    mc : metacooler.MetaCooler
        The map, opened (see check_mge_map).
    mges_list : list of str
        Names of the MGE contigs (see metator.io.import_mges_contigs).

    Returns
    -------
    pandas.DataFrame
        One row per MGE contig of the map, the longest first, ties by name:
        Name, Size (its length), GC_content and Hit (as metator network
        stored them).
    numpy.ndarray
        Their positions in the map, from 0, in that order.
    numpy.ndarray
        The contacts within each, between reads more than 1 kb apart (its
        qc_intra_contig_long).

    Raises
    ------
    metator.metacool.MetacoolError
        If none of the MGE contigs is in the map: there is nothing to
        cluster.
    """
    with mtmc.read_errors(metacool):
        contigs = mc.contigs
    names = contigs["name"].astype(str)
    in_map = set(names)
    missing = [name for name in dict.fromkeys(mges_list) if name not in in_map]
    if missing:
        logger.warning(
            f"{len(missing)} MGE contig(s) are not in {metacool}: they are left out of the MGE clusters: {mtmc._shown(missing)}."
        )
    positions = np.flatnonzero(names.isin(set(mges_list)).to_numpy())
    if not len(positions):
        raise mtmc.MetacoolError(f"None of the MGE contigs is in {metacool}: there is no MGE contig to cluster.")
    # The longest first, ties by name: whatever the order of the map.
    lengths = contigs["length"].to_numpy(dtype=np.int64)[positions]
    positions = positions[np.lexsort((names.to_numpy()[positions], -lengths))]
    selected = contigs.iloc[positions]
    data = pd.DataFrame(
        {
            "Name": names.iloc[positions].to_numpy(),
            "Size": selected["length"].to_numpy(dtype=np.int64),
            "GC_content": selected["GC_content"].to_numpy(dtype=np.float64),
            "Hit": selected["Hit"].to_numpy(dtype=np.int64),
        }
    )
    return data, positions, selected[INTRA_LONG].to_numpy(dtype=np.int64)


def map_matrix(metacool, n_contigs, positions, intra):
    """The contacts of the MGE contigs of a map, as resolve_matrix reads
    them: the matrix of their raw contacts, in the order of positions.

    Between two contigs, the contacts are those of the map, summed by pair of
    contigs (metator.network.map_network, without weights: the counts),
    whatever its bin size; within a contig, those of its reads more than 1
    kb apart (intra, see mge_contigs). Only the upper triangle is filled,
    which resolve_matrix reads: the contacts of the i-th and j-th contigs
    of positions at [min(i, j), max(i, j)], whatever their order in the map.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    n_contigs : int
        Number of contigs of the map.
    positions : numpy.ndarray
        Positions in the map of the MGE contigs, from 0, in the order of the
        matrix (see mge_contigs).
    intra : numpy.ndarray
        Contacts within each, between reads more than 1 kb apart.

    Returns
    -------
    numpy.ndarray
        Raw matrix of contacts between the MGE contigs.
    """
    import metator.network as mtn

    network = mtn.map_network(metacool, None)
    index = np.full(n_contigs, -1, dtype=np.int64)
    index[positions] = np.arange(len(positions))
    # IDs from 1; each pair of contigs once, in the upper triangle of the
    # matrix whatever their order in the map.
    i = index[network["contig1"].to_numpy() - 1]
    j = index[network["contig2"].to_numpy() - 1]
    between = (i >= 0) & (j >= 0)
    i, j = i[between], j[between]
    mat = np.zeros((len(positions), len(positions)))
    mat[np.minimum(i, j), np.maximum(i, j)] = network["weight"].to_numpy()[between]
    mat[np.arange(len(positions)), np.arange(len(positions))] = intra
    return mat


def generate_bin_summary(contigs_data: "pd.DataFrame", mge_bins: dict) -> "pd.DataFrame":
    """The summary of the MGE clusters: one row per cluster, with its length,
    GC content, hits, clustering score and contigs.

    Parameters
    ----------
    contigs_data : pd.DataFrame
        The MGE contigs (see mge_contigs): Name, Size, GC_content and Hit.
    mge_bins : dict
        Dictionary with the mge bin id as key and its contigs ("Contigs")
        and score ("Score") as value.

    Return
    ------
    pd.DataFrame :
        Summary table of the MGE clusters, one row per cluster in the order
        of mge_bins: BinName (MetaTOR_MGE_<ID>), BinLength, GC, Hit,
        BinningScore, ContigsNumber and Contigs (by name, comma-separated).
    """
    # Create output empty DataFrame
    cols = {
        "BinName": pd.Series(dtype="str"),
        "BinLength": pd.Series(dtype="int"),
        "GC": pd.Series(dtype="float"),
        "Hit": pd.Series(dtype="int"),
        "BinningScore": pd.Series(dtype="float"),
        "ContigsNumber": pd.Series(dtype="int"),
        "Contigs": pd.Series(dtype="str"),
    }

    summary = pd.DataFrame(cols, index=mge_bins.keys())
    contigs_data = contigs_data.set_index("Name", drop=False)
    # Iterates on the bins to fill the table.
    for bin_id in mge_bins.keys():
        summary.loc[bin_id, "BinName"] = f"MetaTOR_MGE_{bin_id:05d}"
        contigs = mge_bins[bin_id]["Contigs"]
        summary.loc[bin_id, "ContigsNumber"] = int(f"{len(contigs):1d}")
        summary.loc[bin_id, "Contigs"] = ",".join(contigs)
        summary.loc[bin_id, "BinningScore"] = mge_bins[bin_id]["Score"]
        length, hit, gc = 0, 0, 0
        for contig in contigs:
            length += contigs_data.loc[contig, "Size"]
            hit += contigs_data.loc[contig, "Hit"]
            gc += contigs_data.loc[contig, "GC_content"] * contigs_data.loc[contig, "Size"]
        summary.loc[bin_id, "BinLength"] = int(f"{length:1d}")
        summary.loc[bin_id, "Hit"] = int(f"{hit:1d}")
        summary.loc[bin_id, "GC"] = gc / length
    return summary.reset_index(drop=True)


def mge_binning(
    metacool: str,
    mges_list: List[str],
    threshold_bin: float = 0.8,
    random: bool = False,
) -> Tuple[pd.DataFrame, dict]:
    """Main function to cluster the MGE contigs of a map.

    The map is checked first (see check_mge_map): the refusals come before
    any work. The MGE contigs that the map lacks are left out, with a
    warning (see mge_contigs).

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``: a map that metator
        network stored the network of its contigs in (that of a run of
        network, pipeline or validation).
    mges_list : list of str
        Names of the MGE contigs.
    threshold_bin : float
        Threshold of score to bin contigs. [Default: .8]
    random : bool
        If enabled, make a random shuffling of the bins.

    Returns:
    --------
    pandas.DataFrame:
        The summary of the MGE clusters (see generate_bin_summary).
    dict:
        The MGE clusters: their id as key, their contigs ("Contigs", the
        longest first, ties by name: see mge_contigs) and score ("Score")
        as value.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map cannot be read, cannot be clustered (see check_mge_map),
        or holds none of the MGE contigs.
    """
    with mtmc._metacooler_logging(), mtmc.read_errors(metacool):
        mc = mtmc.open_map(metacool)
        check_mge_map(metacool, mc)
        mges_data, positions, intra = mge_contigs(metacool, mc, mges_list)
        mat = map_matrix(metacool, len(mc.contigs), positions, intra)
    logger.info(
        f"{len(positions)} MGE contig(s) in {metacool}: {int(np.triu(mat, 1).sum())} contacts between them, "
        f"{int(intra.sum())} within them (between reads more than 1 kb apart)."
    )

    # Associates contigs
    bins = resolve_matrix(mat, threshold_bin)

    # Update mges_data and generates bins.
    mges_data, mge_bins = update_mge_data(mges_data, bins)

    if random:
        # Shuffle to simulate random bins.
        mges_data, mge_bins = shuffle_mge_bins(mges_data)

    return generate_bin_summary(mges_data, mge_bins), mge_bins


def write_mges_fasta(metacool, assembly, mge_bins, out_file, temp_dir=None):
    """Write the FASTA file with one record per MGE cluster, so that the tools
    that check the quality of an MGE genome, which expect one sequence per
    genome, read each cluster as one: named after the cluster
    (MetaTOR_MGE_<ID>), its contigs, in the order of the cluster (the
    longest first, ties by name: see mge_binning), joined by 200 N, in lines
    of 80 bases (metacooler.write_fasta, stitched). The file is replaced if
    it exists.

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: the map whose contigs the MGE clusters are.
    assembly : str
        Path to the assembly the map was made from: FASTA, plain or
        compressed with gzip.
    mge_bins : dict
        The MGE clusters (see mge_binning).
    out_file : str
        Path to the FASTA file to write.
    temp_dir : str, optional
        Where a compressed assembly is copied, plain: the temporary
        directory of the run.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the file cannot be written: the assembly lacks a contig, gives
        one another length than the map, or is not FASTA (a ValueError,
        KeyError or OSError of metacooler, message kept).
    """
    import metacooler

    groups = {f"MetaTOR_MGE_{bin_id:05d}": [str(contig) for contig in cluster["Contigs"]] for bin_id, cluster in mge_bins.items()}
    with mtmc._metacooler_logging(), mtmc._as_metacool_errors(), mtmc.read_errors(metacool):
        metacooler.write_fasta(
            metacool,
            assembly,
            out_file,
            groups=groups,
            stitch=True,
            spacer=SPACER,
            single_file=True,
            force=True,
            temp_dir=temp_dir,
        )
    logger.info(f"{len(groups)} clusters have been extracted")


def resolve_matrix(mat: "np.ndarray", threshold: float = 0.8) -> List[Tuple]:
    """Main function to bin mges contigs.

    From the marix of contacts associates the contigs with a lot of
    interactions. To normalize the contacts we divide the contacts by the
    geometric mean of the contacts in intra. To avoid the noise of the close
    range contacts we have remove the contacts below 1000bp. When two contigs
    are binned, they are fused together.

    Parameters:
    -----------
    mat : np.array
        Matrix of the raw contacts between the contigs. Upper triangle and the
        contacts in intra below 1000bp are not kept.
    threshold : float
        Threshold of score to bin contigs. [Default: .8]

    Returns:
    List of tuple:
        List of the Tuple of the associated bins.
    """

    bins: list = []
    n = len(mat)
    # Save a copy of the initial matrix to keep the intra initial values.
    mat0 = np.copy(mat)

    # Normalize values
    for i in range(n):
        for j in range(i + 1, n):
            if (mat[i, i] > 0) and (mat[j, j] > 0):
                mat[i, j] = mat[i, j] / np.sqrt(mat[i, i] * mat[j, j])
            else:
                mat[i, j] = 0

    # Remove the intra contacts.
    for i in range(n):
        mat[i, i] = 0

    # While there is an association bigger than 1 associates the contigs.
    maxi = np.max(mat)
    while maxi > threshold:
        # Handle the case where we have multiple points at the maximum value.
        try:
            i, j = map(int, np.where(mat == maxi))
        except TypeError:
            i, j = map(int, [np.where(mat == maxi)[0][0], np.where(mat == maxi)[1][0]])

        # Compute the sum of the count to create the merge vector.
        a = mat0[i, :] + mat0[:, i] + mat0[j, :] + mat0[:, j]
        # Compute a new intra count.
        intra = mat0[i, i] + mat0[j, j] + mat0[i, j]
        mat0[i, i] = intra
        # Normalize the new vector.
        for k in range(n):
            if k < i:
                mat0[k, i] = a[k]
                if mat0[k, k] > 0:
                    mat[k, i] = a[k] / np.sqrt(intra * mat0[k, k])
                else:
                    mat[k, i] = 0
            elif k > i:
                mat0[i, k] = a[k]
                if mat0[k, k] > 0:
                    mat[i, k] = a[k] / np.sqrt(intra * mat0[k, k])
                else:
                    mat[k, i] = 0
        # Remove the contig j as it has been merged with the contig i.
        mat[j, :] = np.zeros(n)
        mat[:, j] = np.zeros(n)
        mat0[j, :] = np.zeros(n)
        mat0[:, j] = np.zeros(n)

        bins.append([i, j, maxi])
        maxi = np.max(mat)
    return bins


def shuffle_mge_bins(
    mges_data: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Function to shuffle id to imitate a random binning with the same bins
    distribution as the one created by Metator MGE.

    Parameters:
    -----------
    mges_data : pandas.DataFrame
        Table with the contigs name as index and with information from
        metator partition or validation.

    Returns:
    --------
    pandas.DataFrame:
        Input table with the mge bin id column randomly shuffled.
    dict:
        Dictionary with the mge bin id as key and the list of the contigs
        name as value from the shuffle contigs bins.
    """

    # Shuffle the ids of the dataframe
    mges_bin_ids = mges_data.MetaTOR_MGE_bin
    shuffle_ids = np.random.permutation(mges_bin_ids)
    mges_data["shuffle"] = shuffle_ids

    mges_bins: dict = {}
    # Rebuild the bins from the shuffled ids, in the same shape as the real
    # binning: {bin id: {"Contigs": [...], "Score": ...}}.
    for index in mges_data.index:
        contig = mges_data.loc[index, "Name"]
        mge_bin_id = mges_data.loc[index, "shuffle"]
        mge_bin = mges_bins.setdefault(mge_bin_id, {"Contigs": [], "Score": np.nan})
        mge_bin["Contigs"].append(contig)
    return mges_data, mges_bins


def update_mge_data(mges_data: pd.DataFrame, bins: List[Tuple]) -> pd.DataFrame:
    """Function to update the mge bins data.

    Parameters
    ----------
    mges_data : pd.DataFrame
        Table wit the mge contig information.
    bins : list of tuple
        List of pairs of contigs with their score of association.

    Returns
    -------
    pandas.DataFrame :
        Table wit the mge contig information updated.
    dictionnary :
        Dictionary of the mge bins.
    """
    # Initiation
    mges_data["MetaTOR_MGE_bin"] = 0
    mges_data["MetaTOR_MGE_Score"] = 0.0
    bin_id = 0
    mge_bins: dict = {}
    mges_data = mges_data.set_index(np.arange(len(mges_data)))
    for contig_tuple in bins:
        i, j, score = contig_tuple
        # If no existing bin, creates one.
        if (mges_data.loc[i, "MetaTOR_MGE_bin"] == 0) and (mges_data.loc[j, "MetaTOR_MGE_bin"] == 0):
            bin_id += 1
            current_bin = bin_id
            current_score = score
        # If one existing bin, append it.
        elif (mges_data.loc[i, "MetaTOR_MGE_bin"] == 0) or (mges_data.loc[j, "MetaTOR_MGE_bin"] == 0):
            current_bin = max(
                mges_data.loc[i, "MetaTOR_MGE_bin"],
                mges_data.loc[j, "MetaTOR_MGE_bin"],
            )
            current_score = min(
                score,
                max(
                    mges_data.loc[i, "MetaTOR_MGE_Score"],
                    mges_data.loc[j, "MetaTOR_MGE_Score"],
                ),
            )
        # Complicated case as we have to fuse two existing bins. The bin j
        # is not reused.
        else:
            current_bin = mges_data.loc[i, "MetaTOR_MGE_bin"]
            bin_j = mges_data.loc[j, "MetaTOR_MGE_bin"]
            if current_bin > bin_j:
                current_bin -= 1
            bin_id -= 1
            for k in range(len(mges_data)):
                if mges_data.loc[k, "MetaTOR_MGE_bin"] == bin_j:
                    mges_data.loc[k, "MetaTOR_MGE_bin"] = current_bin
                elif mges_data.loc[k, "MetaTOR_MGE_bin"] > bin_j:
                    mges_data.loc[k, "MetaTOR_MGE_bin"] -= 1
            current_score = min(
                score,
                mges_data.loc[i, "MetaTOR_MGE_Score"],
                mges_data.loc[j, "MetaTOR_MGE_Score"],
            )
        mges_data.loc[i, "MetaTOR_MGE_bin"] = current_bin
        mges_data.loc[j, "MetaTOR_MGE_bin"] = current_bin
        mges_data.loc[i, "MetaTOR_MGE_Score"] = current_score
        mges_data.loc[j, "MetaTOR_MGE_Score"] = current_score

    # Update unbinned contigs and build mge bins.
    for k in mges_data.index:
        if mges_data.loc[k, "MetaTOR_MGE_bin"] == 0:
            bin_id += 1
            mges_data.loc[k, "MetaTOR_MGE_bin"] = bin_id
            mge_bins[bin_id] = {
                "Contigs": [mges_data.loc[k, "Name"]],
                "Score": np.nan,
            }
        else:
            curr_bin = mges_data.loc[k, "MetaTOR_MGE_bin"]
            try:
                mge_bins[curr_bin]["Contigs"].append(mges_data.loc[k, "Name"])
                mge_bins[curr_bin]["Score"] = min(
                    mge_bins[curr_bin]["Score"],
                    mges_data.loc[k, "MetaTOR_MGE_Score"],
                )
            except KeyError:
                mge_bins[curr_bin] = {
                    "Contigs": [mges_data.loc[k, "Name"]],
                    "Score": mges_data.loc[k, "MetaTOR_MGE_Score"],
                }
    return mges_data, mge_bins
