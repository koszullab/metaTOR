#!/usr/bin/env python3
# coding: utf-8

"""Module with the mges binning functions.

It generates mges bins from the contigs annotated as mges based on the
metaHiC pairs clusterization based on the sequence and
shotgun co-abundance information.

Core function to partition mges contigs:
    - build_matrix
    - generate_bin_summary
    - generate_mges_bins_pairs
    - generate_mges_fasta
    - mge_binning
    - resolve_matrix
    - shuffle_mge_bins
    - update_mge_data
"""


import numpy as np
import pandas as pd
import pypairix
import subprocess as sp
from metator.log import logger
from os.path import join
from typing import List, Tuple


def build_matrix(contigs: List[str], contigs_size: List[int], pairs_files: List[str]) -> "np.ndarray":
    """Function to extract the pairs from a set of contigs from pairs files. Run
    faster if the files are indexed. The contacts are stored in a raw matrix of
    contacts between the list of given contigs. The contacts beloz 1kb are
    removed.

    Parameters:
    -----------
    contigs : List of str
        List of the mge contigs names uses in the alignment.
    contigs_size : List of int
        List of the size in bp of the contigs (same order as the contigs).
    pairs_files : List of str
        List of the path of the pairs file from the alignment. If possible index
        them first using pypairix.

    Return:
    -------
    np.array:
        Raw matrix of contacts between the given contigs.
    """

    # Initiation
    npairs = 0
    n = len(contigs)
    mat = np.zeros((n, n))
    # Write one pair file for all the ones given.
    for pairs_file in pairs_files:
        # Check if the pypairix index exist
        try:
            pairs_data = pypairix.open(pairs_file)
            pypairix_index = True
        except pypairix.PairixError:
            logger.warning("No pypairix index found. Iterates on the pairs.")
            pypairix_index = False
        # Need a sorted (chr1 chr2 pos1 pos2) pair file indexed with pypairix.
        if pypairix_index:
            for i, contig in enumerate(contigs):
                # Only need to retrieve the upper triangle.
                for j in range(i, len(contigs)):
                    pairs_lines = pairs_data.query2D(
                        contig,
                        0,
                        contigs_size[i],
                        contigs[j],
                        0,
                        contigs_size[j],
                        1,
                    )
                    for p in pairs_lines:
                        npairs += 1
                        # The threshold of 1000 is to remove the close range
                        # contacts.
                        if i == j:
                            if np.abs(int(p[2]) - int(p[4])) > 1000:
                                mat[i, i] += 1
                        else:
                            mat[i, j] += 1
        # else Iterates on the input pairs file (take much longer than with
        # the index).
        else:
            with open(pairs_file, "r") as input_pairs:
                for pairs_line in input_pairs:
                    # Ignore header lines.
                    if pairs_line.startswith("#"):
                        continue
                    # Split the line on the tabulation and check if both contigs
                    # are in the bin.
                    pairs = pairs_line.split("\t")
                    if pairs[1] in contigs and pairs[3] in contigs:
                        npairs += 1
                        i = contigs.index(pairs[1])
                        j = contigs.index(pairs[3])
                        # The threshold of 1000 is to remove the close range
                        # contacts.
                        if i == j:
                            if np.abs(int(pairs[2]) - int(pairs[4])) > 1000:
                                mat[i, i] += 1
                        else:
                            mat[i, j] += 1
    logger.info(f"{npairs} pairs extracted.")
    return mat


def generate_bin_summary(contigs_data: "pd.DataFrame", mge_bins: dict, outfile: str) -> "pd.DataFrame":
    """Function to generate and write the mge binning summary.

    Parameters
    ----------
    contigs_data : pd.DataFrame
        Table with contig information from metator.
    mge_bins : dict
        Dictionary with the mge bin id as key and the list of the contigs
        name as value.
    outfile : str
        Output file where to write the summary table.

    Return
    ------
    pd.DataFrame :
        Summary table of the viral bins.
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
        "MetagenomicBin": pd.Series(dtype="str"),
        "AssociationScore": pd.Series(dtype="float"),
        "AssociatedBins": pd.Series(dtype="str"),
    }

    summary = pd.DataFrame(cols, index=mge_bins.keys())
    contigs_data.set_index("Name", drop=False, inplace=True)
    # Iterates on the bins to fill the table.
    for bin_id in mge_bins.keys():
        summary.loc[bin_id, "BinName"] = f"MetaTOR_MGE_{bin_id:05d}"
        contigs = mge_bins[bin_id]["Contigs"]
        summary.loc[bin_id, "ContigsNumber"] = int(f"{len(contigs):1d}")
        summary.loc[bin_id, "Contigs"] = ",".join(contigs)
        summary.loc[bin_id, "BinningScore"] = mge_bins[bin_id]["Score"]
        summary.loc[bin_id, "MetagenomicBin"] = "None"
        summary.loc[bin_id, "AssociationScore"] = np.nan
        summary.loc[bin_id, "AssociatedBins"] = "NA"
        length, hit, gc = 0, 0, 0
        for contig in contigs:
            length += contigs_data.loc[contig, "Size"]
            hit += contigs_data.loc[contig, "Hit"]
            gc += contigs_data.loc[contig, "GC_content"] * contigs_data.loc[contig, "Size"]
        summary.loc[bin_id, "BinLength"] = int(f"{length:1d}")
        summary.loc[bin_id, "Hit"] = int(f"{hit:1d}")
        summary.loc[bin_id, "GC"] = gc / length

    # Write the summary.
    summary.to_csv(
        outfile,
        sep="\t",
        na_rep="NA",
        float_format="%.2f",
        header=True,
        index=False,
    )
    return summary


def generate_mge_bins_pairs(
    mges_data: pd.DataFrame,
    pairs_files: List[str],
    threshold: float = 0.8,
) -> Tuple[pd.DataFrame, dict]:
    """Generates the binning of the mges contigs based on HiC pairs.

    Parameters:
    -----------
    mges_data : pandas.DataFrame
        Table with the contigs name as index and with information from
        metator partition or validation.
    pairs_file : List of str
        List of the path of the pairs file from the alignment. If possible index
        them first using pypairix.
    threshold : float
        Threshold of score to bin contigs. [Default: .8]

    Returns:
    --------
    pandas.DataFrame:
        Input table with the mge bin id column added.
    dict:
        Dictionary with the mge bin id as key and the list of the contigs
        name as value.
    """

    # Extract contigs and contigs size.
    contigs = list(mges_data.Name)
    contigs_size = list(mges_data.Size)

    # Build matrix
    mat = build_matrix(contigs, contigs_size, pairs_files)

    # Associates contigs
    bins = resolve_matrix(mat, threshold)

    # Update mges_data and generates bins.
    mges_data, mge_bins = update_mge_data(mges_data, bins)

    return mges_data, mge_bins


def generate_mges_fasta(fasta: str, mge_bins: dict, out_file: str, tmp_dir: str):
    """Generate the fasta file with one sequence entry for each bin. The
    sequences are generated like that to be accepted by QC tools as one MgeMAG. In
    the sequences 180 "N" spacers are added between two contigs.

    Parameters:
    -----------
    fasta : str
        Path to the fasta file with the original mges contigs sequences.
    mge_bins : dict
        A dictionnary with the id of the mge bins as keys and the list
        of name of their contigs as values.
    out_file : str
        Path to the output file where the fasta of all the mges bins will be
        written.
    tmp_dir : str
        Path to the temporary directory to write the temporary contigs list
        files.
    """

    nb_bins = 0
    # For each bin create a list of the contigs and extract them from the
    # fasta to create a new fasta file with only the bin.
    with open(out_file, "w") as out:
        for bin_id in mge_bins:
            # Extract the list of the contigs from the contigs data file.
            list_contigs_name = mge_bins[bin_id]["Contigs"]
            nb_bins += 1
            # Create a temporary fasta file.
            contigs_file = join(tmp_dir, f"MetaTOR_MGE_{bin_id:05d}.txt")
            temp_file = join(tmp_dir, f"MetaTOR_MGE_{bin_id:05d}.fa")
            with open(contigs_file, "w") as f:
                for contig_name in list_contigs_name:
                    f.write("%s\n" % contig_name)
            cmd = f"pyfastx extract {fasta} -l {contigs_file} > {temp_file}"
            process = sp.Popen(cmd, shell=True)
            process.communicate()

            # Concatenated the fasta in one sequence entry with 180
            # "N" spacer between contigs.
            with open(temp_file, "r") as tmp:
                start = True
                for line in tmp:
                    if line.startswith(">"):
                        if start:
                            start = False
                            out.write((f">MetaTOR_MGE_{bin_id:05d}\n"))
                        else:
                            out.write("N" * 200)
                    else:
                        out.write(line)
    logger.info(f"{nb_bins} bins have been extracted")

    return mge_bins


def mge_binning(
    fasta_mges_contigs: str,
    contigs_data: pd.DataFrame,
    mges_list_id: List[int],
    out_dir: str,
    pairs_files: List[str],
    tmp_dir: str,
    threshold_bin: float = 0.8,
    random: bool = False,
):
    """Main function to bin mges contigs.

    Generates a fasta file where each entry is a mge bin, with 180bp "N" as
    spacers between contigs.

    Parameters:
    -----------
    fasta_mges_contigs : str
        Path to the fasta containing the mges sequences. It could contain
        other sequences.
    contigs_data : pandas.DataFrame
        Table with the contig name as keys and with the values of the
        contig id the associated bin name, and either if the contig is binned
        and if it's a mge contig. The name of the keys are "id", "bin",
        "binned", "mge".
    mges_list_id : list of int
        List of the mges contigs IDs.
    out_dir : str
        Path to the directory where to write the output data.
    pairs_files : List of str
        List of the path of the pairs file from the alignment. If possible index
        them first using pypairix.
    tmp_dir : str
        Path to temporary directory for intermediate files.
    threshold_bin : float
        Threshold of score to bin contigs. [Default: .8]
    random : bool
        If enabled, make a random shuffling of the bins.
    """

    # Create output and temporary files.
    contigs_file = join(tmp_dir, "mge_contigs.txt")
    temp_fasta = join(tmp_dir, "mges.fa")
    mge_data_file = join(out_dir, "mges_bin_summary.tsv")
    fasta_mges_bins = join(out_dir, "mges_binned.fa")

    # Import contigs data from `metator partition/validation`.
    mges_data = pd.DataFrame(contigs_data.loc[mges_list_id, :])

    with open(contigs_file, "w") as f:
        for contig_name in list(mges_data.Name):
            f.write("%s\n" % contig_name)

    # Extract fasta to have sequences at the same order as the depth file.
    cmd = "pyfastx extract {0} -l {1} > {2}".format(fasta_mges_contigs, contigs_file, temp_fasta)
    process = sp.Popen(cmd, shell=True)
    process.communicate()
    mges_data, mge_bins = generate_mge_bins_pairs(mges_data, pairs_files, threshold_bin)

    if random:
        # Shuffle to simulate random bins. Uncomment to do it
        mges_data, mge_bins = shuffle_mge_bins(mges_data)

    # Generate MGE fasta (one sequence per mgeMAG, with 180bp "N" as spacers).
    mge_bins = generate_mges_fasta(fasta_mges_contigs, mge_bins, fasta_mges_bins, tmp_dir)

    summary = generate_bin_summary(mges_data, mge_bins, mge_data_file)
    return summary


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
    # Shuffle mges bins according to the shuffle dataframe.
    for index in mges_data.index:
        contig = mges_data.loc[index, "Name"]
        mge_bin_id = mges_data.loc[index, "shuffle"]
        try:
            mges_bins[mge_bin_id]["Contig"].append(contig)
        except KeyError:
            mges_bins[mge_bin_id]["Contig"] = [contig]
            mges_bins[mge_bin_id]["Score"] = np.nan
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
    mges_data.set_index(np.arange(len(mges_data)), inplace=True)
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
