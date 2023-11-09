#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Figures foor metaTOR output.

General utility functions to plot some figures to describe metaTOR output.

Core functions to plot the firgures are:
    - barplot_bins_number
    - barplot_bins_size
    - build_matrix_network
    - build_vmags_summary
    - figures_bins_distribution
    - figures_bins_size_distribution
    - figure_camembert_quality
    - figures_mags_GC_boxplots
    - figures_mags_HiC_cov_boxplots
    - figures_mags_SG_cov_boxplots
    - generates_frags_network
    - network_heatmap
    - pie_bins_size_distribution
    - plot_figures
    - reindex_df
"""


import hicstuff.hicstuff as hcs
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import seaborn as sns
from metator.log import logger
from os.path import join


def barplot_bins_number(labels, list_checkv_summary, out_file):
    """Plot the number distribution of the bins quality of different experiments
    or methods.

    Parameters:
    -----------
    labels : list os str
        List of the labels of the different experiments/methods.
    list_checkv_summary : list of pandas.DataFrame
        List of the checkV quality summary tables.
    out_file : str
        Path to the output file.
    """
    # Create an empty dataframe.
    stacked_data = pd.DataFrame(
        [[0] * 5] * len(labels),
        columns=[">90", ">70", ">50", ">30", "<30"],
        index=labels,
    )

    # Extract vMAGs summary and add it to the table.
    for i, checkv_summary in enumerate(list_checkv_summary):
        mags_summary = build_vmags_summary(checkv_summary)
        stacked_data.loc[labels[i], :] = list(mags_summary["bins"])

    # Plot the table.
    stacked_data = stacked_data.apply(lambda x: x * 100 / sum(x), axis=1)
    stacked_data.plot(
        kind="bar",
        stacked=True,
        color=[
            "#313695",
            "#4575b4",
            "#abd9e9",
            "#fdae61",
            "#a50026",
        ],
    )
    plt.title("MGEs bins quality distribution")
    plt.xlabel("Experiment")
    plt.ylabel("Percentage bins number (%)")
    plt.legend(loc="upper right", bbox_to_anchor=(0.2, 0.0, 1.0, 1.0))
    # Save the file
    plt.savefig(out_file, dpi=200, bbox_inches="tight")


def barplot_bins_size(labels, list_checkv_summary, out_file):
    """Plot the size distribution of the bins quality of different experiments
    or methods.

    Parameters:
    -----------
    labels : list os str
        List of the labels of the different experiments/methods.
    list_checkv_summary : list of pandas.DataFrame
        List of the checkV quality summary tables.
    out_file : str
        Path to the output file.
    """
    # Create an empty dataframe.
    stacked_data = pd.DataFrame(
        [[0] * 5] * len(labels),
        columns=[">90", ">70", ">50", ">30", "<30"],
        index=labels,
    )

    # Extract vMAGs summary and add it to the table.
    for i, checkv_summary in enumerate(list_checkv_summary):
        mags_summary = build_vmags_summary(checkv_summary)
        stacked_data.loc[labels[i], :] = list(mags_summary["size"])

    # Plot the table.
    stacked_data = stacked_data.apply(lambda x: x * 100 / sum(x), axis=1)
    stacked_data.plot(
        kind="bar",
        stacked=True,
        color=[
            "#313695",
            "#4575b4",
            "#abd9e9",
            "#fdae61",
            "#a50026",
        ],
    )
    plt.title("MGEs bins quality distribution")
    plt.xlabel("Experiment")
    plt.ylabel("Percentage bins size (%)")
    plt.legend(loc="upper right", bbox_to_anchor=(0.2, 0.0, 1.0, 1.0))
    # Save the file
    plt.savefig(out_file, dpi=200, bbox_inches="tight")


def build_matrix_network(pairs_files, info_contigs, N, bin_size):
    """Build a binned matrix from the pairs alignement files and the final
    binning information in the contig data file.

    Parameters:
    -----------
    pairs_files : list of str
        List of path of the alignment ".pairs" file.
    info_contigs: dict
        Dictionnary with frags id of the contigs.
    N : int
        Size of the final bin matrix.
    bin_size : int
        Size of the bin to plot the heatmap.

    Returns:
    --------
    numpy.array:
        Binned dense matrix of the network.
    """
    # Initiate the matix
    matrix = np.zeros((N, N))

    # Iterates on the input pairs file
    for pairs_file in pairs_files:
        with open(pairs_file, "r") as input_pairs:
            for pairs_line in input_pairs:
                # Ignore header lines.
                if pairs_line.startswith("#"):
                    continue
                # Split the line on the tabulation and check if both contigs
                # are in the bin.
                pairs = pairs_line.split("\t")
                contig1, contig2 = pairs[1], pairs[3]
                pos1, pos2 = pairs[2], pairs[4]
                try:
                    frag1 = (
                        info_contigs[contig1]["start"]
                        + (int(pos1) + info_contigs[contig1]["init"])
                        // bin_size
                    )
                    frag2 = (
                        info_contigs[contig2]["start"]
                        + (int(pos2) + info_contigs[contig2]["init"])
                        // bin_size
                    )
                    if frag1 < frag2:
                        matrix[frag1, frag2] += 1
                    elif frag1 > frag2:
                        matrix[frag2, frag1] += 1
                except KeyError:
                    continue
    matrix_norm = hcs.normalize_dense(
        matrix + matrix.T, norm="SCN", iterations=100
    )
    return matrix_norm


def build_vmags_summary(checkv_summary):
    """Build vMAGs quality summary table.

    Parameters:
    -----------
    checkv_summary : pandas.core.frame.DataFrame
        Table with the informations about the final MGE bins.
    """
    # Add a qualitive quality column in bin_summary with the quality of the MAGs
    # checkv_summary["MAG_quality"] = "ND"
    # Change NA value to 0
    mask = checkv_summary.completeness == "NA"
    checkv_summary.loc[mask, "completeness"] = 0
    checkv_summary = checkv_summary.loc[checkv_summary.provirus == "No", :]
    checkv_summary = checkv_summary.reset_index()
    # Build a small table with the sum for each quality category.
    mags_summary = pd.DataFrame(
        [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]],
        columns=["bins", "size"],
    )
    for i in range(len(checkv_summary)):
        completness = float(checkv_summary.loc[i, "completeness"])
        size = int(checkv_summary.loc[i, "contig_length"])
        if completness >= 90:
            mags_summary.loc[0, "bins"] += 1
            mags_summary.loc[0, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = ">90"
        if completness >= 70:
            mags_summary.loc[1, "bins"] += 1
            mags_summary.loc[1, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = ">70"
        elif completness >= 50:
            mags_summary.loc[2, "bins"] += 1
            mags_summary.loc[2, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = ">50"
        elif completness >= 30:
            mags_summary.loc[3, "bins"] += 1
            mags_summary.loc[3, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = ">30"
        else:
            mags_summary.loc[4, "bins"] += 1
            mags_summary.loc[4, "size"] += size
            # checkv_summary.loc[i, "MAG_quality"] = "<30"
    return mags_summary


def figures_bins_distribution(bin_summary, out_file):
    """Function to plot the distribution of the estimated completion and
    contamination of the final bins generated by metaTOR. The bins are ordered
    by decreasing completion and contamination values superior to 105% are set
    to 105%.

    Parameters:
    -----------
    bin_summary : pandas.core.frame.DataFrame
        Table with the informations about the final bins.
    out_file : str
        Path where to save the figure.
    """
    # Transform values as float
    bin_summary["Weighted completeness"] = bin_summary[
        "Weighted completeness"
    ].apply(float)
    bin_summary["Weighted redundancy"] = bin_summary[
        "Weighted redundancy"
    ].apply(float)
    # Sort the values by decreasing completion.
    bin_summary = bin_summary.sort_values(
        by="Weighted completeness", ascending=False
    )
    # Put the contamination values bigger than 105 to 105.
    mask = bin_summary["Weighted redundancy"] >= 105
    bin_summary.loc[mask, "Weighted redundancy"] = 105
    # Plot the contamination and the completion.
    _fig, ax = plt.subplots()
    ax.tick_params(axis="both", which="both", length=0)
    plt.box(on=None)
    plt.scatter(
        y=bin_summary["Weighted completeness"] * 100,
        x=range(len(bin_summary)),
        color="r",
        label="Completion",
        s=2.5,
    )
    plt.scatter(
        y=bin_summary["Weighted redundancy"] * 100 - 100,
        x=range(len(bin_summary)),
        color="k",
        label="Contamination",
        s=2.5,
    )
    plt.axhline(y=5, color="k", linestyle=(0, (5, 10)))
    plt.axhline(y=10, color="k", linestyle=(0, (5, 10)))
    plt.axhline(y=50, color="k", linestyle=(0, (5, 10)))
    plt.axhline(y=70, color="k", linestyle=(0, (5, 10)))
    plt.axhline(y=90, color="k", linestyle=(0, (5, 10)))
    plt.axhline(y=100, color="k", linestyle=(0, (5, 10)))
    plt.yticks([5, 10, 50, 70, 90, 100])
    plt.xticks(np.arange(0, len(bin_summary), 50))
    plt.xlabel("Rank index")
    plt.ylabel("Completness/Contamination (%)")
    plt.legend(
        bbox_to_anchor=(0.0, -0.25, 1.0, 1.0),
        loc="lower center",
        ncol=2,
        borderaxespad=0.0,
    )
    # Save the file
    plt.savefig(out_file, dpi=200, bbox_inches="tight")
    plt.close()


def figures_bins_size_distribution(
    bin_summary, total_size, threshold, out_file
):
    """Function to plot a camembert of the fraction of size corresponding
    to each quality of bins. We defined 6 categories:
        - High quality MAGs (HQ MAGs): >= 90% completion ; < 5% contamination
        - Medium quality MAGs (MQ MAGs): >= 70% completion ; < 10% contamination
        - Low quality MAGs (LQ MAGs): >= 50% completion ; < 10% contamination
        - Contaminated bins: >= 50% completion ; >= 10% contamination
        - Bins superior to size threshold but with less 50% completion.
        - Unbinned contigs

    Parameters:
    -----------
    bin_summary : pandas.core.frame.DataFrame
        Table with the informations about the final bins.
    total_size : int
        Size of the whole assembly.
    threshold : int
        Minimun size of the bins considered.
    out_file : str
        Path where to save the figure.
    """
    # Add a qualitive quality column in bin_summary with the quality of the MAGs
    bin_summary["MAG_quality"] = "ND"
    # Build a small table with the sum for each quality category.
    mags_summary = pd.DataFrame(
        [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [None, 0]],
        columns=["bins", "size"],
    )
    for i in range(len(bin_summary)):
        completness = float(bin_summary.loc[i, "Weighted completeness"])
        contamination = float(bin_summary.loc[i, "Weighted redundancy"])
        size = int(bin_summary.loc[i, "Length"])
        if completness >= 0.5:
            if contamination > 1.1:
                mags_summary.loc[3, "bins"] += 1
                mags_summary.loc[3, "size"] += size
                bin_summary.loc[i, "MAG_quality"] = "Contaminated"
            else:
                if completness >= 0.9 and contamination <= 1.05:
                    mags_summary.loc[0, "bins"] += 1
                    mags_summary.loc[0, "size"] += size
                    bin_summary.loc[i, "MAG_quality"] = "HQ"
                elif completness >= 0.7:
                    mags_summary.loc[1, "bins"] += 1
                    mags_summary.loc[1, "size"] += size
                    bin_summary.loc[i, "MAG_quality"] = "MQ"
                else:
                    mags_summary.loc[2, "bins"] += 1
                    mags_summary.loc[2, "size"] += size
                    bin_summary.loc[i, "MAG_quality"] = "LQ"
        else:
            mags_summary.loc[4, "bins"] += 1
            mags_summary.loc[4, "size"] += size
            bin_summary.loc[i, "MAG_quality"] = "Other"
    mags_summary.loc[5, "size"] = total_size - sum(mags_summary["size"])
    # Plot the camembert of the size ratio.
    labels = [
        "HQ MAGs: {0} - {1}Mb".format(
            int(mags_summary.loc[0, "bins"]),
            round(mags_summary.loc[0, "size"] / 1000000, 2),
        ),
        "MQ MAGs: {0} - {1}Mb".format(
            int(mags_summary.loc[1, "bins"]),
            round(mags_summary.loc[1, "size"] / 1000000, 2),
        ),
        "LQ MAGs: {0} - {1}Mb".format(
            int(mags_summary.loc[2, "bins"]),
            round(mags_summary.loc[2, "size"] / 1000000, 2),
        ),
        "Contaminated bins: {0} - {1}Mb".format(
            int(mags_summary.loc[3, "bins"]),
            round(mags_summary.loc[3, "size"] / 1000000, 2),
        ),
        "Others bins: {0} - {1}Mb".format(
            int(mags_summary.loc[4, "bins"]),
            round(mags_summary.loc[4, "size"] / 1000000, 2),
        ),
        "Others contigs - {0}Mb".format(
            round(mags_summary.loc[5, "size"] / 1000000, 2),
        ),
    ]
    plt.pie(
        mags_summary["size"],
        colors=[
            "#103b6f",
            "#6666ff",
            "#ccccff",
            "#ed3139",
            "#c0bcbc",
            "#eaeded",
        ],
    )
    plt.legend(labels, bbox_to_anchor=(0.9, 0.0, 1.0, 1.0), loc="upper right")
    plt.text(
        -1.5,
        -1.2,
        "Total size of the assembly: {0}Mb".format(
            round(total_size / 1000000, 2)
        ),
        fontdict=None,
    )
    plt.text(
        -1.5,
        -1.35,
        "Percentage of MAGs: {0}%".format(
            round(sum(mags_summary["size"][0:3]) / total_size * 100, 2)
        ),
        fontdict=None,
    )
    plt.text(
        -1.5,
        -1.5,
        "Bins threshold: {0}kb".format(round(threshold / 1000, 2),),
        fontdict=None,
    )
    plt.title("Size proportion of bins depending on their quality")
    # Save the file
    plt.savefig(out_file, dpi=200, bbox_inches="tight")
    plt.close()
    return bin_summary


def figure_camembert_quality(
    out_file,
    prefix,
    n_religated,
    n_loops,
    n_weirds,
    n_informative_intra,
    n_informative_inter,
    n_intra_mags,
    n_inter_mags,
    uncut_thr,
    loop_thr,
):
    """Functon to plot nice hicstuff camembert from metator high quality
    contigs.
    """
    plt.figure(2, figsize=(6, 6))
    # The slices will be ordered and plotted counter-clockwise.
    total = n_inter_mags + n_intra_mags
    fracs = [
        n_religated,
        n_loops,
        n_weirds,
        n_informative_intra,
        n_informative_inter,
        n_inter_mags,
    ]
    # Format labels to include event names and proportion
    labels = list(
        map(
            lambda x: (x[0] + ": %.2f%%") % (100 * x[1] / total),
            [
                ("Religated", n_religated),
                ("Loops", n_loops),
                ("Weirds", n_weirds),
                ("Informative intracontigs", n_informative_intra),
                ("Informative intercontigs", n_informative_inter),
                ("Noise inter", n_inter_mags),
            ],
        )
    )
    colors = ["#D55E00", "#E69F00", "#999999", "#103b6f", "#6096fd", "#cc0000"]
    patches, _ = plt.pie(fracs, colors=colors, startangle=90)
    plt.legend(
        patches, labels, loc="upper left", bbox_to_anchor=(-0.1, 1.0),
    )
    if prefix:
        plt.title(
            "Distribution of library events in {}".format(prefix),
            bbox={"facecolor": "1.0", "pad": 5},
        )
    plt.text(
        0.3, 1.15, "Threshold Uncuts = " + str(uncut_thr), fontdict=None,
    )
    plt.text(
        0.3, 1.05, "Threshold Loops = " + str(loop_thr), fontdict=None,
    )
    plt.text(
        -1.5,
        -1.2,
        f"Total number of reads in the estimation: {total / 1_000_000:.2f} millions reads",
        fontdict=None,
    )
    noise = 100 * n_inter_mags / (n_intra_mags + n_inter_mags)
    plt.text(
        -1.5, -1.3, f"Estimated noise signal = {noise:.2f}%", fontdict=None,
    )
    informative = (
        100
        * (n_informative_intra + n_informative_inter)
        / (n_intra_mags + n_inter_mags)
    )
    plt.text(
        -1.5, -1.4, f"Informative reads = {informative:.2f}%", fontdict=None,
    )
    plt.savefig(out_file)
    plt.close()


def figures_mags_GC_boxplots(contigs_data, out_file):
    """Function to plot barplots GC content of the contigs for each bins. The
    count are weighted by the size of the contigs.

    Parameters:
    -----------
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs, the ubinned contigs should
        have removed and a column size_weight should have been add to weight
        with the size.
    out_file : str
        Path where to write the figure.
    """
    # Sort by decreasing GC.
    contigs_data = contigs_data.sort_values(by="GC-content", ascending=False)
    # Build the palette
    my_pal = {
        "HQ": "#103b6f",
        "MQ": "#6666ff",
        "LQ": "#ccccff",
        "Contaminated": "#ed3139",
        "Other": "#c0bcbc",
    }
    # Plot the figure.
    _fig, ax = plt.subplots(figsize=(20, 10))
    boxplot = sns.boxplot(
        y="GC_content",
        x="Final_bin",
        data=reindex_df(contigs_data, "size_weight"),
        palette=my_pal,
        hue="MAG_quality",
        hue_order=["HQ", "MQ", "LQ", "Contaminated", "Other"],
        flierprops=dict(
            markerfacecolor="0.5",
            markersize=1.0,
            marker="o",
            markeredgewidth=0.0,
        ),
        linewidth=1.25,
        width=1.0,
        dodge=False,
    )
    ax.set(xlabel="", ylabel="GC content distribution", xticklabels=[])
    labels = [
        "HQ MAGs (Completion > 90, Contamination < 5)",
        "MQ MAGs (Completion > 70, Contamination < 10)",
        "LQ MAGs (Completion > 50, Contamination < 10)",
        "Contaminated (Completion > 50, Contamination < 10)",
        "Other (Completion < 50)",
    ]
    handles, _ = boxplot.get_legend_handles_labels()
    boxplot.legend(handles, labels)
    # Save the file.
    plt.savefig(out_file, dpi=200, bbox_inches="tight")
    plt.close()


def figures_mags_HiC_cov_boxplots(contigs_data, out_file):
    """Function to plot barplots HiC_coverage coverage of the contigs for each
    bins. The count are weighted by the size of the contigs.

    Parameters:
    -----------
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs, the ubinned contigs should
        have removed and a column size_weight should have been add to weight
        with the size.
    out_file : str
        Path where to write the figure.
    """
    # Compute the HiC coverage.
    contigs_data["HiC_cov"] = 100 * contigs_data["Hit"] / contigs_data["Size"]
    # Sort by decreasing HiC coverage.
    contigs_data = contigs_data.sort_values(by="HiC_cov", ascending=False)
    # Build the palette
    my_pal = {
        "HQ": "#103b6f",
        "MQ": "#6666ff",
        "LQ": "#ccccff",
        "Contaminated": "#ed3139",
        "Other": "#c0bcbc",
    }
    # Plot the figure.
    _fig, ax = plt.subplots(figsize=(20, 10))
    boxplot = sns.boxplot(
        y="HiC_cov",
        x="Final_bin",
        data=reindex_df(contigs_data, "size_weight"),
        palette=my_pal,
        hue="MAG_quality",
        hue_order=["HQ", "MQ", "LQ", "Contaminated", "Other"],
        flierprops=dict(
            markerfacecolor="0.5",
            markersize=1.0,
            marker="o",
            markeredgewidth=0.0,
        ),
        linewidth=1.25,
        width=1.0,
        dodge=False,
    )
    plt.yscale("log")
    ax.set(xlabel="", ylabel="HiC coverage distribution", xticklabels=[])
    labels = [
        "HQ MAGs (Completion > 90, Contamination < 5)",
        "MQ MAGs (Completion > 70, Contamination < 10)",
        "LQ MAGs (Completion > 50, Contamination < 10)",
        "Contaminated (Completion > 50, Contamination < 10)",
        "Other (Completion < 50)",
    ]
    handles, _ = boxplot.get_legend_handles_labels()
    boxplot.legend(handles, labels)
    # Save the file.
    plt.savefig(out_file, dpi=200, bbox_inches="tight")
    plt.close()


def figures_mags_SG_cov_boxplots(contigs_data, out_file):
    """Function to plot barplots assembly coverage of the contigs for each bins.
    The count are weighted by the size of the contigs.

    Parameters:
    -----------
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs, the ubinned contigs should
        have removed and a column size_weight should have been add to weight
        with the size.
    out_file : str
        Path where to write the figure.
    """
    # Sort by decreasing Shotgun coverage.
    contigs_data = contigs_data.sort_values(by="SG_Coverage", ascending=False)
    # Build the palette
    my_pal = {
        "HQ": "#103b6f",
        "MQ": "#6666ff",
        "LQ": "#ccccff",
        "Contaminated": "#ed3139",
        "Other": "#c0bcbc",
    }
    # Plot the figure.
    _fig, ax = plt.subplots(figsize=(20, 10))
    boxplot = sns.boxplot(
        y="Shotgun_coverage",
        x="Final_bin",
        data=reindex_df(contigs_data, "size_weight"),
        palette=my_pal,
        hue="MAG_quality",
        hue_order=["HQ", "MQ", "LQ", "Contaminated", "Other"],
        flierprops=dict(
            markerfacecolor="0.5",
            markersize=1.0,
            marker="o",
            markeredgewidth=0.0,
        ),
        linewidth=1.25,
        width=1.0,
        dodge=False,
    )
    plt.yscale("log")
    ax.set(xlabel="", ylabel="Assembly coverage distribution", xticklabels=[])
    labels = [
        "HQ MAGs (Completion > 90, Contamination < 5)",
        "MQ MAGs (Completion > 70, Contamination < 10)",
        "LQ MAGs (Completion > 50, Contamination < 10)",
        "Contaminated (Completion > 50, Contamination < 10)",
        "Other (Completion < 50)",
    ]
    handles, _ = boxplot.get_legend_handles_labels()
    boxplot.legend(handles, labels)
    # Save the file.
    plt.savefig(out_file, dpi=200, bbox_inches="tight")
    plt.close()


def generates_frags_network(contig_data, bin_summary, bin_size):
    """Generates info fragments for all contigs and start position all bins.

    Parameters:
    -----------
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs.
    bin_summary : dict
        Dictionnary with the informations of the final bins kept by MetaTOR.
    bin_size : int
        Size of the bin to plot the heatmap.

    Returns:
    --------
    dict:
        Dictionnary with frags id of the contigs.
    list of int:
        List of the start positions of the mags
    int:
        Size of the final binned matrix
    """

    # Define start frag id for each MAGs. The rest column is the size still
    # available from the previous contigs.
    frag_id = 0
    mag_starts = []
    for bin_name in bin_summary:
        bin_summary[bin_name]["start"] = frag_id
        mag_starts.append(frag_id)
        bin_summary[bin_name]["current"] = frag_id
        bin_summary[bin_name]["rest"] = 0
        frag_id += (int(bin_summary[bin_name]["Length"]) // bin_size) + 1

    # Define start frag id for each bins
    info_contigs = dict()
    for i in contig_data.index:
        # Only add the contig if it's in a final bin.
        if contig_data.loc[i, "Final_bin"] != "ND":
            # Call bin name dict to have the status of start frag id and rest.
            bin_name = contig_data.loc[i, "Final_bin"]
            size = contig_data.loc[i, "Size"]
            info_contigs[contig_data.loc[i, "Name"]] = {
                "start": bin_summary[bin_name]["current"],
                "init": bin_summary[bin_name]["rest"],
            }
            # Update rest and start position.
            bin_summary[bin_name]["current"] += (
                size + bin_summary[bin_name]["rest"]
            ) // bin_size
            bin_summary[bin_name]["rest"] = (
                size + bin_summary[bin_name]["rest"]
            ) % bin_size

    return info_contigs, mag_starts, frag_id


def network_heatmap(matrix, out_file=None, mag_starts=None):
    """Plot the heatmap of the network, with the contig ordering by their bin
    attribution.

    Parameters:
    -----------
    matrix : numpy.array
        Binned dense matrix of the network.
    out_file : str
        Name of the file to save the plot.
    mag_starts : list of int
        List of bin positions where to draw mags starts as dotted lines.
    """

    # Plot the heatmap.
    vmax = np.percentile(matrix, 99.5)
    im_kwargs = {
        "vmin": 0,
        "vmax": vmax,
        "cmap": "viridis",
        "interpolation": "none",
    }
    li_kwargs = {"ls": ":", "alpha": 0.8, "c": "white", "linewidth": 0.1}
    plt.figure()
    plt.imshow(matrix, **im_kwargs)
    plt.colorbar()
    plt.axis("off")

    # Add line between MAGs.
    if mag_starts is not None:
        for pos in mag_starts:
            plt.axvline(pos, **li_kwargs)
            plt.axhline(pos, **li_kwargs)

    # Save or show the file.
    if out_file:
        plt.savefig(out_file, bbox_inches="tight", pad_inches=0.0, dpi=2000)
        plt.close()
    else:
        plt.show()


def pie_bins_size_distribution(checkv_summary, out_file):
    """Function to plot a camembert of the fraction of size corresponding
    to their completness.

    Parameters:
    -----------
    checkv_summary : pandas.core.frame.DataFrame
        Table with the informations about the final MGE bins.
    out_file : str
        Path where to save the figure.
    """
    mags_summary = build_vmags_summary(checkv_summary)
    # Plot the camembert of the size ratio.
    labels = [
        ">90%: {0} - {1}Mb".format(
            int(mags_summary.loc[0, "bins"]),
            round(mags_summary.loc[0, "size"] / 1000000, 2),
        ),
        ">70%: {0} - {1}Mb".format(
            int(mags_summary.loc[1, "bins"]),
            round(mags_summary.loc[1, "size"] / 1000000, 2),
        ),
        ">50%: {0} - {1}Mb".format(
            int(mags_summary.loc[2, "bins"]),
            round(mags_summary.loc[2, "size"] / 1000000, 2),
        ),
        ">30%: {0} - {1}Mb".format(
            int(mags_summary.loc[3, "bins"]),
            round(mags_summary.loc[3, "size"] / 1000000, 2),
        ),
        "Others contigs: {0} - {1}Mb".format(
            int(mags_summary.loc[4, "bins"]),
            round(mags_summary.loc[4, "size"] / 1000000, 2),
        ),
    ]
    total_size = np.sum(list(mags_summary["size"]))
    fig, ax = plt.subplots()
    plt.pie(
        mags_summary["size"],
        colors=[
            "#313695",
            "#4575b4",
            "#abd9e9",
            "#fdae61",
            "#a50026",
        ],
    )
    plt.legend(labels, bbox_to_anchor=(0.9, 0.0, 1.0, 1.0), loc="upper right")
    plt.text(
        -1.5,
        -1.2,
        "Total size of the assembly: {0}Mb".format(
            round(total_size / 1000000, 2)
        ),
        fontdict=None,
    )
    plt.text(
        -1.5,
        -1.35,
        "Percentage of MGE 50% complete: {0}%".format(
            round(sum(mags_summary["size"][0:3]) / total_size * 100, 2)
        ),
        fontdict=None,
    )
    plt.title("Size proportion of MGE bins depending on their completness.")
    # Save the file
    plt.savefig(out_file, dpi=200, bbox_inches="tight")


def plot_figures(out_dir, contigs_data, bin_summary, threshold):
    """Function to generates all figures.

    Parameters:
    -----------
    out_dir : str
        Path of the output directory of MetaTOR.
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs.
    bin_summary : dict
        Dictionnary with the informations of the final bins kept by MetaTOR.
    threshold : int
        Minimun size of the bins considered.
    """
    # Create plot directory
    plot_dir = join(out_dir, "plot")
    os.makedirs(plot_dir, exist_ok=True)

    # Test if Shotgun coverage have been given
    if contigs_data.loc[contigs_data.index[0], "Shotgun_coverage"] == "-":
        SG_cov = False
    else:
        SG_cov = True

    # Create outfiles
    outfile_bins_distribution = join(plot_dir, "bins_distribution.png")
    outfile_bins_size_distribution = join(
        plot_dir, "bins_size_distribution.png"
    )
    outfile_MAGs_GC = join(plot_dir, "MAGs_GC_distribution.png")
    outfile_MAGs_HiC_cov = join(plot_dir, "MAGs_HiC_cov_distribution.png")
    if SG_cov:
        outfile_MAGs_SG_cov = join(plot_dir, "MAGs_SG_cov_distribution.png")

    # Look for pairs file in out_dir
    pairs_files = []
    list_files = os.listdir(out_dir)
    for file in list_files:
        if re.search("\.pairs$", file):
            pairs_files.append(join(out_dir, file))

    # Plot heatmap with a binning of 50kb.
    if len(pairs_files) > 0:
        heatmap_file = join(plot_dir, "network_heatmap.png")
        bin_size = 50000
        info_contigs, mag_starts, N = generates_frags_network(
            contigs_data, bin_summary, bin_size
        )
        matrix = build_matrix_network(pairs_files, info_contigs, N, bin_size)
        network_heatmap(matrix, heatmap_file, mag_starts)
    else:
        logger.warning(
            "No pairs alignment files found in %s, no heatmap will be generated.",
            out_dir,
        )

    # Transform dictionnary to pandas DataFrame.
    bin_summary = pd.DataFrame.from_dict(bin_summary, orient="index")
    bin_summary["Bin"] = bin_summary.index
    bin_summary.index = range(len(bin_summary))
    # Compute size of the assembly.
    total_size = sum(contigs_data["Size"])

    # Plot bin distribution
    figures_bins_distribution(bin_summary, outfile_bins_distribution)
    bin_summary = figures_bins_size_distribution(
        bin_summary, total_size, threshold, outfile_bins_size_distribution
    )

    # Remove unbinned contigs.
    contigs_data = contigs_data[contigs_data["Final_bin"] != "ND"]
    # Create a column of size divided by the minimum size of the contigs.
    min_size = min(contigs_data["Size"])
    contigs_data["size_weight"] = (contigs_data["Size"] / min_size).apply(int)
    # Merge the bin_summary and the contigs data file.
    contigs_data = pd.merge(
        contigs_data, bin_summary, left_on="Final_bin", right_on="Bin"
    )

    # Plots the distribution of GC and coverage inside the bins.
    figures_mags_GC_boxplots(contigs_data, outfile_MAGs_GC)
    figures_mags_HiC_cov_boxplots(contigs_data, outfile_MAGs_HiC_cov)
    if SG_cov:
        figures_mags_SG_cov_boxplots(contigs_data, outfile_MAGs_SG_cov)


def reindex_df(dataframe, weight_col):
    """Expand the dataframe to prepare for resampling result is 1 row per count
    per sample.

    Parameters:
    -----------
    dataframe : pandas.core.frame.DataFrame
        Table to expand.
    weight_col : str
        Name of the weighted column to use.

    Returns:
    --------
    pandas.core.frame.DataFrame:
        Table expanded with new index.
    """
    dataframe = dataframe.reindex(dataframe.index.repeat(dataframe[weight_col]))
    dataframe.reset_index(drop=True, inplace=True)
    return dataframe
