#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Figures foor metator output.

General utility functions to plot some figures to describe metator output.

Core functions to plot the firgures are:
    - build_matrix_map
    - build_matrix_network
    - figures_bins_distribution
    - figures_bins_size_distribution
    - figure_camembert_quality
    - figures_mags_GC_boxplots
    - figures_mags_HiC_cov_boxplots
    - figures_mags_SG_cov_boxplots
    - generates_frags_network
    - heatmap_windows
    - network_heatmap
    - plot_figures
    - plot_heatmap
    - reindex_df
"""

import metator.io as mio
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from matplotlib.patches import Patch
from metator.log import logger
from os.path import join


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
        # metator network leaves bgzip-compressed pairs files in its outdir.
        with mio.read_compressed(pairs_file) as input_pairs:
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
                    frag1 = info_contigs[contig1]["start"] + (int(pos1) + info_contigs[contig1]["init"]) // bin_size
                    frag2 = info_contigs[contig2]["start"] + (int(pos2) + info_contigs[contig2]["init"]) // bin_size
                    if frag1 < frag2:
                        matrix[frag1, frag2] += 1
                    elif frag1 > frag2:
                        matrix[frag2, frag1] += 1
                except KeyError:
                    continue
    return normalize_scn(matrix + matrix.T)


def heatmap_windows(bins, contigs_data, bin_summary, bin_size):
    """The windows of the heatmap that the bins of a map cover.

    The MAGs are laid end to end in the order of bin_summary, the contigs of
    each in the order of contigs_data, and cut into windows of bin_size bp,
    the windows of build_matrix_network (generates_frags_network): MAG m
    starts at the window after the last of the MAG before it, and spans
    Length // bin_size + 1 windows. A bin of the map covers the windows that
    its stretch of its contig falls in, each for the share of the bin that
    falls there: a map holds no read positions, so the contacts of a bin are
    spread over its windows as its length is.

    Parameters:
    -----------
    bins : pandas.DataFrame
        The bins of the map, in its order: chrom (the contig name), start
        and end.
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs: Name, Size and Final_cluster
        (ND for a contig in no MAG).
    bin_summary : dict
        The MAGs, in their order, by name, with their Length.
    bin_size : int
        Size of a window, in bp.

    Returns:
    --------
    scipy.sparse.csr_matrix:
        One row per bin of the map, one column per window: the share of the
        bin in the window. The bins of the contigs in no MAG have none.
    list of int:
        The first window of each MAG.
    int:
        The number of windows.
    """
    from scipy.sparse import csr_matrix

    spans = np.array([int(info["Length"]) // bin_size + 1 for info in bin_summary.values()], dtype=np.int64)
    mag_starts = np.cumsum(spans) - spans
    first_window = dict(zip(bin_summary, mag_starts.tolist(), strict=True))
    # Where each contig of a MAG starts in the MAG laid end to end.
    final = contigs_data[contigs_data["Final_cluster"] != "ND"]
    sizes = final["Size"].astype(np.int64)
    offsets = sizes.groupby(final["Final_cluster"].to_numpy()).cumsum() - sizes
    placed = pd.DataFrame(
        {"offset": offsets.to_numpy(), "window": final["Final_cluster"].map(first_window).to_numpy()},
        index=final["Name"].astype(str).to_numpy(),
    )
    on_contig = placed.reindex(bins["chrom"].astype(str).to_numpy())
    kept = on_contig["window"].notna().to_numpy() & (bins["end"] > bins["start"]).to_numpy()
    rows = np.flatnonzero(kept)
    begin = (on_contig["offset"].to_numpy()[kept] + bins["start"].to_numpy()[kept]).astype(np.int64)
    end = (on_contig["offset"].to_numpy()[kept] + bins["end"].to_numpy()[kept]).astype(np.int64)
    first, last = begin // bin_size, (end - 1) // bin_size
    # One entry per window that a bin covers: most cover one or two.
    counts = (last - first + 1).astype(np.int64)
    bin_rows = np.repeat(rows, counts)
    within = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
    windows = np.repeat(first, counts) + within
    lows = np.maximum(np.repeat(begin, counts), windows * bin_size)
    highs = np.minimum(np.repeat(end, counts), (windows + 1) * bin_size)
    shares = (highs - lows) / np.repeat(end - begin, counts)
    columns = windows + np.repeat(on_contig["window"].to_numpy()[kept].astype(np.int64), counts)
    size = int(sum(spans))
    return csr_matrix((shares, (bin_rows, columns)), shape=(len(bins), size)), mag_starts.tolist(), size


def build_matrix_map(metacool, contigs_data, bin_summary, bin_size, chunksize=10_000_000):
    """Build the matrix of the heatmap from the pixels of a .metacool map of
    the contigs, at any bin size: the matrix that build_matrix_network
    builds from the pairs the map holds, but for where the contacts fall
    within a bin of the map, which holds no read positions: the contacts of
    a bin, with another bin or within itself, are spread over its windows
    by their shares (see heatmap_windows), as reads that fall anywhere in
    it. As in build_matrix_network, the contacts within a window are left
    out. The matrix is normalized as build_matrix_network normalizes it.

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection.
    contigs_data : pandas.DataFrame
        Table with all the data from the contigs (see heatmap_windows).
    bin_summary : dict
        The MAGs, in their order, by name, with their Length.
    bin_size : int
        Size of a window of the heatmap, in bp.
    chunksize : int
        Number of pixels read at a time.

    Returns:
    --------
    numpy.array:
        Normalized dense matrix of the windows of the MAGs.
    list of int:
        The first window of each MAG.
    """
    import metator.metacool as mtmc

    with mtmc._metacooler_logging():
        mc = mtmc.open_map(metacool)
        windows, mag_starts, size = heatmap_windows(mc.bins()[:], contigs_data, bin_summary, bin_size)
        matrix = np.zeros((size, size))
        for chunk in mc.iter_pixels(chunksize=chunksize):
            counts = chunk["count"].to_numpy().astype(np.float64)
            first = windows[chunk["bin1_id"].to_numpy()].multiply(counts[:, None]).tocsr()
            matrix += (first.T @ windows[chunk["bin2_id"].to_numpy()]).toarray()
    # The map stores each contact once, in the upper triangle of its bins
    # (its diagonal included), as build_matrix_network counts each pair
    # once; the contacts within a window are left out, as there.
    matrix += matrix.T
    np.fill_diagonal(matrix, 0)
    return normalize_scn(matrix), mag_starts


def normalize_scn(matrix, iterations=100):
    """The SCN normalization of hicstuff (hicstuff.hicstuff.normalize_dense,
    norm "SCN") of a symmetric matrix, to within rounding, but in place:
    each iteration divides the rows, then the columns, by their sums, those
    of a sum of 0 left as they are, and the result is made symmetric, the
    mean of the matrix and its transpose. It makes two copies of the
    matrix, where hicstuff makes several temporary ones at each iteration
    (a mask, a matrix of the sums, the masked values)."""
    matrix = np.array(matrix, np.float64)
    for _ in range(iterations):
        sums = matrix.sum(axis=1)
        np.divide(matrix, sums[:, None], out=matrix, where=sums[:, None] != 0)
        sums = matrix.sum(axis=0)
        np.divide(matrix, sums[None, :], out=matrix, where=sums[None, :] != 0)
    symmetric = matrix + matrix.T
    symmetric /= 2
    return symmetric


def figures_bins_distribution(bin_summary, out_file):
    """Function to plot the distribution of the estimated completion and
    contamination of the final bins generated by metator. The bins are ordered
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
    bin_summary["Weighted completeness"] = bin_summary["Weighted completeness"].apply(float)
    bin_summary["Weighted redundancy"] = bin_summary["Weighted redundancy"].apply(float)
    # Sort the values by decreasing completion.
    bin_summary = bin_summary.sort_values(by="Weighted completeness", ascending=False)
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


def bin_quality(completeness, redundancy):
    """The quality class of a bin, from miComplete's weighted completeness and
    weighted redundancy, its contamination being the redundancy minus 1:

        - "HQ" (high quality MAG): >= 90% completion ; <= 5% contamination
        - "MQ" (medium quality MAG): >= 70% completion ; <= 10% contamination
        - "LQ" (low quality MAG): >= 50% completion ; <= 10% contamination
        - "Contaminated": >= 50% completion ; > 10% contamination
        - "Other": < 50% completion

    The plot of the quality of the bins and the log of the validation both
    count the bins in these classes.

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


def figures_bins_size_distribution(bin_summary, total_size, threshold, out_file):
    """Function to plot a camembert of the fraction of size corresponding
    to each quality of bins: the classes of bin_quality (HQ, MQ and LQ MAGs,
    contaminated bins, and other bins, with less than 50% completion), and the
    unbinned contigs.

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
    rows = {"HQ": 0, "MQ": 1, "LQ": 2, "Contaminated": 3, "Other": 4}
    for i in range(len(bin_summary)):
        quality = bin_quality(bin_summary.loc[i, "Weighted completeness"], bin_summary.loc[i, "Weighted redundancy"])
        mags_summary.loc[rows[quality], "bins"] += 1
        mags_summary.loc[rows[quality], "size"] += int(bin_summary.loc[i, "Length"])
        bin_summary.loc[i, "MAG_quality"] = quality
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
        "Contaminated clusters: {0} - {1}Mb".format(
            int(mags_summary.loc[3, "bins"]),
            round(mags_summary.loc[3, "size"] / 1000000, 2),
        ),
        "Other clusters: {0} - {1}Mb".format(
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
        "Total size of the assembly: {0}Mb".format(round(total_size / 1000000, 2)),
        fontdict=None,
    )
    plt.text(
        -1.5,
        -1.35,
        "Percentage of MAGs: {0}%".format(round(sum(mags_summary["size"][0:3]) / total_size * 100, 2)),
        fontdict=None,
    )
    plt.text(
        -1.5,
        -1.5,
        "Clusters threshold: {0}kb".format(
            round(threshold / 1000, 2),
        ),
        fontdict=None,
    )
    plt.title("Size proportion of clusters depending on their quality")
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
    # The slices will be ordered and plotted counter-clockwise.
    total = n_inter_mags + n_intra_mags
    if not total:
        logger.warning(f"No pair to plot: {out_file} is not drawn.")
        return
    plt.figure(2, figsize=(6, 6))
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
        patches,
        labels,
        loc="upper left",
        bbox_to_anchor=(-0.1, 1.0),
    )
    if prefix:
        plt.title(
            "Distribution of library events in {}".format(prefix),
            bbox={"facecolor": "1.0", "pad": 5},
        )
    plt.text(
        0.3,
        1.15,
        "Threshold Uncuts = " + str(uncut_thr),
        fontdict=None,
    )
    plt.text(
        0.3,
        1.05,
        "Threshold Loops = " + str(loop_thr),
        fontdict=None,
    )
    plt.text(
        -1.5,
        -1.2,
        f"Total number of reads in the estimation: {total / 1_000_000:.2f} millions reads",
        fontdict=None,
    )
    noise = 100 * n_inter_mags / (n_intra_mags + n_inter_mags)
    plt.text(
        -1.5,
        -1.3,
        f"Estimated noise signal = {noise:.2f}%",
        fontdict=None,
    )
    informative = 100 * (n_informative_intra + n_informative_inter) / (n_intra_mags + n_inter_mags)
    plt.text(
        -1.5,
        -1.4,
        f"Informative reads = {informative:.2f}%",
        fontdict=None,
    )
    plt.savefig(out_file)
    plt.close()


# Colours and legend labels of the MAG quality classes, in legend order.
MAG_QUALITY_COLOURS = {
    "HQ": "#103b6f",
    "MQ": "#6666ff",
    "LQ": "#ccccff",
    "Contaminated": "#ed3139",
    "Other": "#c0bcbc",
}
MAG_QUALITY_LABELS = {
    "HQ": "HQ MAGs (Completion > 90, Contamination < 5)",
    "MQ": "MQ MAGs (Completion > 70, Contamination < 10)",
    "LQ": "LQ MAGs (Completion > 50, Contamination < 10)",
    "Contaminated": "Contaminated (Completion > 50, Contamination < 10)",
    "Other": "Other (Completion < 50)",
}
# Outlines, whiskers, caps and medians of the boxes.
BOX_LINE_COLOUR = "0.15"


def _quality_boxplots(contigs_data, column, ylabel, out_file, log_scale=False):
    """Plot one box per bin of a contig value, coloured by the bin quality.

    Each contig counts once per unit of its size weight (see reindex_df).
    The bins are laid out from left to right in their order of first
    appearance in the table, so the caller chooses it by sorting the
    table. Contigs of a bin whose MAG_quality is not a known class are
    left out, and so are missing and infinite values.

    Parameters:
    -----------
    contigs_data : pandas.DataFrame
        Table of the binned contigs, with the columns Final_cluster,
        MAG_quality, size_weight and the plotted column.
    column : str
        Column of the value to plot.
    ylabel : str
        Label of the y axis.
    out_file : str
        Path where to write the figure.
    log_scale : bool
        Whether to plot the values on a log scale.
    """
    data = reindex_df(contigs_data, "size_weight")
    data = data[data["MAG_quality"].isin(MAG_QUALITY_COLOURS)]
    bins = data.groupby("Final_cluster", sort=False)
    values = [group[column].dropna().to_numpy(dtype=float) for _, group in bins]
    values = [box[np.isfinite(box)] for box in values]
    colours = [MAG_QUALITY_COLOURS[group["MAG_quality"].iloc[0]] for _, group in bins]

    fig, ax = plt.subplots(figsize=(20, 10))
    lines = {"color": BOX_LINE_COLOUR, "linewidth": 1.25}
    boxes = ax.boxplot(
        values,
        positions=range(len(values)),
        widths=1.0,
        patch_artist=True,
        boxprops={"edgecolor": BOX_LINE_COLOUR, "linewidth": 1.25},
        whiskerprops=lines,
        capprops=lines,
        medianprops=lines,
        flierprops={
            "marker": "o",
            "markersize": 1.0,
            "markerfacecolor": "0.5",
            "markeredgewidth": 0.0,
        },
    )
    for box, colour in zip(boxes["boxes"], colours, strict=True):
        box.set_facecolor(colour)
    if log_scale:
        ax.set_yscale("log")
    ax.set(xlabel="", ylabel=ylabel, xticklabels=[])
    # Every class is in the legend, whether or not a bin has it.
    ax.legend(
        handles=[
            Patch(
                facecolor=colour,
                edgecolor=BOX_LINE_COLOUR,
                linewidth=1.25,
                label=MAG_QUALITY_LABELS[quality],
            )
            for quality, colour in MAG_QUALITY_COLOURS.items()
        ]
    )
    fig.savefig(out_file, dpi=200, bbox_inches="tight")
    plt.close(fig)


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
    _quality_boxplots(contigs_data, "GC_content", "GC content distribution", out_file)


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
    _quality_boxplots(contigs_data, "HiC_cov", "HiC coverage distribution", out_file, log_scale=True)


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
    contigs_data = contigs_data.sort_values(by="Shotgun_coverage", ascending=False)
    _quality_boxplots(
        contigs_data,
        "Shotgun_coverage",
        "Assembly coverage distribution",
        out_file,
        log_scale=True,
    )


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
        if contig_data.loc[i, "Final_cluster"] != "ND":
            # Call bin name dict to have the status of start frag id and rest.
            bin_name = contig_data.loc[i, "Final_cluster"]
            size = contig_data.loc[i, "Size"]
            info_contigs[contig_data.loc[i, "Name"]] = {
                "start": bin_summary[bin_name]["current"],
                "init": bin_summary[bin_name]["rest"],
            }
            # Update rest and start position.
            bin_summary[bin_name]["current"] += (size + bin_summary[bin_name]["rest"]) // bin_size
            bin_summary[bin_name]["rest"] = (size + bin_summary[bin_name]["rest"]) % bin_size

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
    plt.figure()
    plt.imshow(matrix, **im_kwargs)
    plt.colorbar()
    plt.axis("off")
    # Saved at about one pixel per window, from 100 to 2000 dpi. It was
    # saved at 2000 dpi whatever its size, several pixels a window, which
    # took gigabytes of memory to draw a few thousand windows.
    axes = plt.gca().get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
    dpi = int(min(2000, max(100, np.ceil(len(matrix) / axes.height))))

    # Add line between MAGs: a pixel wide at that dpi, at least 0.1 pt (they
    # were 0.1 pt, which is no pixel below 720 dpi).
    li_kwargs = {"ls": ":", "alpha": 0.8, "c": "white", "linewidth": max(0.1, 72 / dpi)}
    if mag_starts is not None:
        for pos in mag_starts:
            plt.axvline(pos, **li_kwargs)
            plt.axhline(pos, **li_kwargs)

    # Save or show the file.
    if out_file:
        plt.savefig(out_file, bbox_inches="tight", pad_inches=0.0, dpi=dpi)
        plt.close()
    else:
        plt.show()


def plot_heatmap(out_dir, plot_dir, contigs_data, bin_summary, metacool=None):
    """Draw the heatmap of the network of the MAGs, network_heatmap.png in
    plot_dir, in windows of 50 kb: from metacool, the .metacool map of the
    contigs of the run (see build_matrix_map), or, without one, from the
    pairs files of out_dir (see build_matrix_network).

    It is a figure: one that cannot be drawn is a warning, and the run goes
    on. The warning gives the message of a MetacoolError, or the type,
    message and traceback of another error, which is unexpected. With no
    MAG, none is drawn.
    """
    from metator.metacool import MetacoolError

    if not bin_summary:
        logger.info("No MAG: the heatmap of the clusters is not drawn.")
        return
    heatmap_file = join(plot_dir, "network_heatmap.png")
    bin_size = 50000
    pairs_files = [join(out_dir, file) for file in os.listdir(out_dir) if file.endswith(".pairs.gz")]
    try:
        if metacool is not None:
            matrix, mag_starts = build_matrix_map(metacool, contigs_data, bin_summary, bin_size)
            network_heatmap(matrix, heatmap_file, mag_starts)
        elif pairs_files:
            info_contigs, mag_starts, N = generates_frags_network(contigs_data, bin_summary, bin_size)
            matrix = build_matrix_network(pairs_files, info_contigs, N, bin_size)
            network_heatmap(matrix, heatmap_file, mag_starts)
        else:
            logger.warning(
                "No pairs alignment files found in %s, no heatmap will be generated.",
                out_dir,
            )
    except Exception as error:
        plt.close("all")
        unexpected = not isinstance(error, MetacoolError)
        reason = f"{type(error).__name__}: {error}" if unexpected else str(error)
        logger.warning(f"The heatmap of the clusters is not drawn: {reason}", exc_info=unexpected)


def plot_figures(out_dir, contigs_data, bin_summary, threshold, metacool=None):
    """Function to generates all figures.

    The heatmap of the network of the MAGs is drawn from metacool, the
    .metacool map of the contigs of the run (see build_matrix_map), whatever
    the run started from; without one, from the pairs files of out_dir
    (see build_matrix_network), and none without those.

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
    metacool : str or None
        Path to the .metacool map of the contigs of the run, or URI of an
        entry of a .mmetacool collection.
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
    outfile_bins_size_distribution = join(plot_dir, "bins_size_distribution.png")
    outfile_MAGs_GC = join(plot_dir, "MAGs_GC_distribution.png")
    outfile_MAGs_HiC_cov = join(plot_dir, "MAGs_HiC_cov_distribution.png")
    if SG_cov:
        outfile_MAGs_SG_cov = join(plot_dir, "MAGs_SG_cov_distribution.png")

    plot_heatmap(out_dir, plot_dir, contigs_data, bin_summary, metacool=metacool)

    # Transform dictionnary to pandas DataFrame.
    bin_summary = pd.DataFrame.from_dict(bin_summary, orient="index")
    bin_summary["Bin"] = bin_summary.index
    bin_summary.index = range(len(bin_summary))
    # Compute size of the assembly.
    total_size = sum(contigs_data["Size"])

    # Plot bin distribution
    figures_bins_distribution(bin_summary, outfile_bins_distribution)
    bin_summary = figures_bins_size_distribution(bin_summary, total_size, threshold, outfile_bins_size_distribution)

    # Remove unbinned contigs.
    contigs_data = contigs_data[contigs_data["Final_cluster"] != "ND"]
    # Create a column of size divided by the minimum size of the contigs.
    min_size = min(contigs_data["Size"])
    contigs_data["size_weight"] = (contigs_data["Size"] / min_size).apply(int)
    # Merge the bin_summary and the contigs data file.
    contigs_data = pd.merge(contigs_data, bin_summary, left_on="Final_cluster", right_on="Bin")

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
    return dataframe.reset_index(drop=True)
