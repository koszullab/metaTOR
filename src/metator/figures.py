#!/usr/bin/env python3
# coding: utf-8

"""The plots of the final clusters of a validation.

``validation`` and ``pipeline`` draw them in ``plot/``, in their output
directory, once the clusters are final (see plot_figures):

- ``bins_distribution.png``: the completeness and the contamination of each
  cluster, the clusters by decreasing completeness (see
  figures_bins_distribution);
- ``bins_size_distribution.png``: the share of the contigs of the run in each
  quality class of the clusters (see metator.validation.bin_quality), and in
  no cluster (see figures_bins_size_distribution).

They are drawn with matplotlib's Figure, without pyplot: no window, whatever
the backend of matplotlib.
"""

import os
from os.path import join

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from metator.log import logger

# The quality classes of the clusters (see metator.validation.bin_quality), in
# the order of the plots, with their colour.
MAG_QUALITY_COLOURS = {
    "HQ": "#103b6f",
    "MQ": "#6666ff",
    "LQ": "#ccccff",
    "Contaminated": "#ed3139",
    "Other": "#c0bcbc",
}
# The name of each class in the legend of the pie chart.
MAG_QUALITY_NAMES = {
    "HQ": "HQ MAGs",
    "MQ": "MQ MAGs",
    "LQ": "LQ MAGs",
    "Contaminated": "Contaminated clusters",
    "Other": "Other clusters",
}
# The contigs in no final cluster, in the pie chart.
UNCLUSTERED_COLOUR = "#eaeded"

# The thresholds of the quality classes, in percent (contamination: 5 and 10;
# completeness: 50, 70 and 90), and 100: the dashed lines of the dot plot.
QUALITY_LINES = (5, 10, 50, 70, 90, 100)
# The contamination drawn above it is drawn at it, in percent.
MAX_CONTAMINATION = 105

# The columns of the summary of the clusters that the plots read.
PLOTTED_COLUMNS = ["Length", "Weighted completeness", "Weighted redundancy"]


def _saved(fig, out_file):
    """Write the figure fig to out_file, as a PNG of 200 dpi cut to its content."""
    fig.savefig(out_file, dpi=200, bbox_inches="tight")


def figures_bins_distribution(bin_summary, out_file):
    """Plot the completeness and the contamination of each final cluster: a dot
    each, in percent, the clusters by decreasing completeness (their rank on
    the x axis), with dashed lines at the thresholds of the quality classes
    and at 100%.

    The contamination is the weighted redundancy minus 1 (see
    metator.validation.bin_quality), drawn from 0 to 105%: a cluster in which
    miComplete finds no marker gene, of a redundancy of 0, at 0, and a
    contamination of more than 105% at 105%.

    Parameters
    ----------
    bin_summary : pandas.DataFrame
        The final clusters, one per row, with miComplete's Weighted
        completeness and Weighted redundancy, as fractions.
    out_file : str
        Path of the PNG file to write.
    """
    completeness = bin_summary["Weighted completeness"].astype(float).to_numpy() * 100
    contamination = bin_summary["Weighted redundancy"].astype(float).to_numpy() * 100 - 100
    contamination = np.clip(contamination, 0, MAX_CONTAMINATION)
    # By decreasing completeness; clusters of the same completeness in the
    # order of the summary.
    order = np.argsort(-completeness, kind="stable")
    ranks = np.arange(len(order))

    fig = Figure()
    ax = fig.subplots()
    ax.tick_params(axis="both", which="both", length=0)
    ax.set_frame_on(False)
    ax.scatter(ranks, completeness[order], color="r", label="Completion", s=2.5)
    ax.scatter(ranks, contamination[order], color="k", label="Contamination", s=2.5)
    for line in QUALITY_LINES:
        ax.axhline(y=line, color="k", linestyle=(0, (5, 10)))
    ax.set_yticks(QUALITY_LINES)
    ax.set_xticks(np.arange(0, len(order), 50))
    ax.set_xlabel("Rank index")
    ax.set_ylabel("Completeness/Contamination (%)")
    ax.legend(bbox_to_anchor=(0.0, -0.25, 1.0, 1.0), loc="lower center", ncol=2, borderaxespad=0.0)
    _saved(fig, out_file)


def figures_bins_size_distribution(bin_summary, total_size, threshold, out_file):
    """Plot a pie chart of the size of the contigs in each quality class of the
    final clusters (see metator.validation.bin_quality: HQ, MQ and LQ MAGs,
    contaminated clusters and other clusters), and in no final cluster.

    The legend counts the clusters of each class and their size, in Mb; below
    the pie, the total size of the contigs, the share of the MAGs (HQ, MQ and
    LQ) in it, and the size threshold of the clusters.

    Parameters
    ----------
    bin_summary : pandas.DataFrame
        The final clusters, one per row, with their Length and miComplete's
        Weighted completeness and Weighted redundancy.
    total_size : int
        The size of the contigs of the run, in bp: those of its map, the
        whole assembly for a map of all its contigs.
    threshold : int
        The size threshold of the clusters (--size), in bp.
    out_file : str
        Path of the PNG file to write.

    Returns
    -------
    pandas.DataFrame
        A copy of bin_summary, with the quality class of each cluster in a
        MAG_quality column.
    """
    from metator.validation import bin_quality

    bin_summary = bin_summary.copy()
    bin_summary["MAG_quality"] = [
        bin_quality(completeness, redundancy)
        for completeness, redundancy in zip(bin_summary["Weighted completeness"], bin_summary["Weighted redundancy"], strict=True)
    ]
    qualities = bin_summary["MAG_quality"]
    lengths = bin_summary["Length"].astype(float).astype(np.int64)
    counts = {quality: int((qualities == quality).sum()) for quality in MAG_QUALITY_COLOURS}
    sizes = {quality: int(lengths[qualities == quality].sum()) for quality in MAG_QUALITY_COLOURS}
    unclustered = total_size - sum(sizes.values())
    mags = sizes["HQ"] + sizes["MQ"] + sizes["LQ"]

    fig = Figure()
    ax = fig.subplots()
    ax.pie([*sizes.values(), unclustered], colors=[*MAG_QUALITY_COLOURS.values(), UNCLUSTERED_COLOUR])
    labels = [
        f"{MAG_QUALITY_NAMES[quality]}: {counts[quality]} - {round(sizes[quality] / 1_000_000, 2)}Mb"
        for quality in MAG_QUALITY_COLOURS
    ]
    labels.append(f"Other contigs - {round(unclustered / 1_000_000, 2)}Mb")
    ax.legend(labels, bbox_to_anchor=(0.9, 0.0, 1.0, 1.0), loc="upper right")
    ax.text(-1.5, -1.2, f"Total size of the assembly: {round(total_size / 1_000_000, 2)}Mb")
    ax.text(-1.5, -1.35, f"Percentage of MAGs: {round(mags / total_size * 100, 2)}%")
    ax.text(-1.5, -1.5, f"Clusters threshold: {round(threshold / 1000, 2)}kb")
    ax.set_title("Size proportion of clusters depending on their quality")
    _saved(fig, out_file)
    return bin_summary


def plot_figures(plot_dir, contigs_data, bin_summary, threshold):
    """Draw the plots of the final clusters of a validation in plot_dir, made if
    missing: bins_distribution.png (see figures_bins_distribution) and
    bins_size_distribution.png (see figures_bins_size_distribution).

    They are figures: one that cannot be drawn is a warning, with the type,
    message and traceback of the error, and the others are drawn; a run goes
    on to write the map of its clusters.

    Parameters
    ----------
    plot_dir : str
        The directory of the plots: ``plot/`` in the output directory of the
        run.
    contigs_data : pandas.DataFrame
        The contig data of the validation: the Size of each contig of the run,
        those of its map.
    bin_summary : dict
        The final clusters, by name, with their Length and miComplete's
        Weighted completeness and Weighted redundancy.
    threshold : int
        The size threshold of the clusters (--size), in bp.
    """
    table = pd.DataFrame(
        [[info[column] for column in PLOTTED_COLUMNS] for info in bin_summary.values()],
        index=list(bin_summary),
        columns=PLOTTED_COLUMNS,
    )
    total_size = int(contigs_data["Size"].sum())
    plots = {
        "bins_distribution.png": lambda path: figures_bins_distribution(table, path),
        "bins_size_distribution.png": lambda path: figures_bins_size_distribution(table, total_size, threshold, path),
    }
    try:
        os.makedirs(plot_dir, exist_ok=True)
    except OSError as error:
        logger.warning(f"The plots of the clusters are not drawn: {plot_dir} cannot be made ({error}).")
        return
    for name, draw in plots.items():
        path = join(plot_dir, name)
        try:
            draw(path)
        except Exception as error:
            logger.warning(f"{path} is not drawn: {type(error).__name__}: {error}", exc_info=True)
        else:
            logger.debug(f"{path} is drawn.")
