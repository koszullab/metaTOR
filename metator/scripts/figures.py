#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quickly draw figures from metaTOR calls.
"""

import matplotlib

matplotlib.use("Agg")

import argparse

import numpy as np

import matplotlib.backends.backend_pdf
from metator.scripts import hicstuff as hcs
from matplotlib import pyplot as plt
from scipy import sparse


# Use (prettier) seaborn if available
SEABORN = False

try:
    import seaborn as sns

    SEABORN = True
except ImportError:
    pass

DEFAULT_INTERVALS = (0, 2, 10, 50, 100, 250, 500, 1000, np.inf)
DEFAULT_SATURATION_THRESHOLD = 80
DEFAULT_MAX_SIZE_MATRIX = 2000
DEFAULT_DPI = 200


def spaceless_pdf_plot_maker(array, filename, vmax=None, dpi=DEFAULT_DPI):
    """Draw a pretty plot from an array

    A function that performs all the tedious matplotlib
    magic to draw a 2D array with as few parameters and
    as little whitespace as possible.

    Parameters
    ----------
    array : array_like
        The input array to draw.
    filename : file, str or pathlib.Path
        The output image to save the array into.
    vmax : float, optional
        The default saturation threshold for the array. If set to None, the
        80th percentile value of the array is chosen. Default is None.
    dpi : int, optional
        Dots per inch (DPI) of the output image. Default is 200.

    Returns
    -------
    None
    """

    if vmax is None:
        vmax = np.percentile(array, DEFAULT_SATURATION_THRESHOLD)
    plt.gca().set_axis_off()
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.figure()
    if SEABORN:
        sns.heatmap(array, vmax=vmax, cmap="Reds")
    else:
        plt.imshow(array, vmax=vmax, cmap="Reds", interpolation="none")
        plt.colorbar()
    plt.savefig(filename, bbox_inches="tight", pad_inches=0.0, dpi=dpi)
    plt.close()


def draw_sparse_matrix(
    array_filename,
    output_image,
    vmax=DEFAULT_SATURATION_THRESHOLD,
    max_size_matrix=DEFAULT_MAX_SIZE_MATRIX,
):
    """Draw a quick preview of a sparse matrix with automated
    binning and normalization.
    """

    matrix = np.loadtxt(array_filename, dtype=np.int32, skiprows=1)
    try:
        row, col, data = matrix.T
    except ValueError:
        row, col, data = matrix
    size = max(np.amax(row), np.amax(col)) + 1
    S = sparse.coo_matrix((data, (row, col)), shape=(size, size))
    if max_size_matrix <= 0:
        binning = 1
    else:
        binning = (size // max_size_matrix) + 1
    binned_S = hcs.bin_sparse(S, subsampling_factor=binning)
    dense_S = binned_S.todense()
    dense_S = dense_S + dense_S.T - np.diag(np.diag(dense_S))
    normed_S = hcs.normalize_dense(dense_S)
    spaceless_pdf_plot_maker(normed_S, output_image, vmax=vmax)


def make_barplots(sizes_file, output, intervals=DEFAULT_INTERVALS):

    sizes = np.loadtxt(sizes_file)

    data_for_barplot = [sum(sizes)]

    n = len(intervals) - 1
    for i in range(n):
        fits_interval = (sizes > intervals[i]) * (sizes < intervals[i + 1])
        chunks_in_interval = sizes[fits_interval].sum()
        data_for_barplot.append(chunks_in_interval)

    labels = ["Total"]

    for i in range(n):
        if intervals[i + 1] == np.inf:
            labels.append("> {}".format(intervals[i]))
        elif intervals[i] == 0:
            labels.append("1")
        else:
            labels.append("{} - {}".format(intervals[i], intervals[i + 1]))

    if SEABORN:
        sns.barplot(labels, data_for_barplot)
    else:
        y = np.arange(n + 1)
        plt.bar(y, data_for_barplot, align="center", alpha=0.5)
        plt.xticks(y, labels)

    plt.xlabel("Size of parent core (in chunks)")
    plt.ylabel("Number of chunks")
    plt.title("Distribution of chunks by size of parent core")
    plt.savefig(output, bbox_inches="tight", pad_inches=0.0)


def draw_regression(sizes_file, output):

    try:
        size_label = "{}".format(sizes_file.split(".")[0].split("_")[-1])
    except ValueError:
        size_label = ""

    sizes_array = np.loadtxt(sizes_file, usecols=(0, 1))
    indices = sizes_array[:, 0]
    values = sizes_array[:, 1]

    if SEABORN:
        sns.regplot(indices, values, fit_reg=False)
    else:
        plt.scatter(indices, values)

    plt.xlabel("Number of iterations")
    plt.ylabel("Number of bins")

    plt.title("Evolution of bins > {}".format(size_label))
    plt.savefig(output, bbox_inches="tight", pad_inches=0.0)


def draw_enrichments(output, *files):

    sizes = []
    labels = []

    n = len(files)

    all_files = list(files)

    for enrichment_file in all_files:
        data_enrichments = np.loadtxt(enrichment_file, usecols=(1, 2))
        label = enrichment_file.split("/")[-1].split("_")[0]

        my_size = data_enrichments[:, 0]
        my_hits = data_enrichments[:, 1]

        my_total_hits = np.int32(np.sum(my_hits))

        unweighted_sizes = np.zeros(my_total_hits)

        cursor = 0
        for u, v in zip(my_size, my_hits):
            slice_size = np.int32(v)
            unweighted_sizes[cursor : cursor + slice_size] = u
            cursor += slice_size

        sizes.append(np.log10(unweighted_sizes[unweighted_sizes > 1]))
        labels.append(label)

    if SEABORN:
        sizes = [size.tolist() for size in sizes]
        sns.violinplot(
            data=sizes,
            width=1,
            inner="box",
            cut=0,
            linewidth=.3,
            bw=.18,
            color="grey",
        )
    else:
        labels = [""] + labels
        plt.violinplot(
            sizes,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            widths=1,
            bw_method=.18,
        )

    y = np.arange(n)

    plt.xticks(y, labels)

    plt.xlabel("Number of hits")
    plt.ylabel("log(Sizes)")

    plt.title("Distribution of hits vs bin sizes")
    plt.savefig(output, bbox_inches="tight", pad_inches=0.0)


def draw_logplots(enrichment_files, output):

    proportional = "proportion" in enrichment_files or "proportion" in output

    enrichment_data = np.loadtxt(enrichment_files, usecols=(1, 2))
    sizes = enrichment_data[:, 0]
    hits = enrichment_data[:, 1]

    if SEABORN:
        sns.regplot(sizes, hits, fit_reg=False)
    else:
        plt.scatter(sizes, hits)

    plt.xlabel("Bin size")
    plt.ylabel("Number of hits")

    if proportional:
        plt.xscale("log")
        plt.title("Bin size vs relative hits")
    else:
        plt.loglog()
        plt.title("Bin size vs hits")
    plt.savefig(
        "{}_nhit.pdf".format(output), bbox_inches="tight", pad_inches=0.0
    )
    plt.close()


def main():

    parser = argparse.ArgumentParser(
        description="Draw and save matrices," "barplots, boxplots and so on."
    )

    parser.add_argument(
        "-p", "--plots", help="Draw scatterplots for core size evolution"
    )

    parser.add_argument("-m", "--matrices", help="Draw contact maps")

    parser.add_argument(
        "-b",
        "--barplots",
        help="Draw barplots for" "parent core size distribution",
    )

    parser.add_argument(
        "-v",
        "--violins",
        help="Draw violin plots to view" "hits vs. size distribution",
        nargs="*",
    )

    parser.add_argument(
        "-l",
        "--logplots",
        help="Draw log plots to view"
        "bin size distribution for different hits",
    )

    parser.add_argument("-s", "--sparse", help="Draw sparse matrix")

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output image file to draw in any extension",
    )

    args = parser.parse_args()

    plots = args.plots
    matrices = args.matrices
    barplots = args.barplots
    logplots = args.logplots
    violins = args.violins
    output_file = args.output
    sparse_matrix = args.sparse

    if barplots:
        input_file = barplots
        make_barplots(input_file, output_file)

    if matrices:
        input_file = matrices
        spaceless_pdf_plot_maker(input_file, output_file)

    if plots:
        input_file = plots
        draw_regression(input_file, output_file)

    if logplots:
        input_file = logplots
        draw_logplots(input_file, output_file)

    if violins:
        input_files = violins
        draw_enrichments(output_file, *input_files)

    if sparse_matrix:
        input_file = sparse_matrix
        draw_sparse_matrix(input_file, output_file)


if __name__ == "__main__":

    main()
