# Test for figures module

import metator.figures as mtf
import pandas as pd
import shutil


def test_figures_bins_distribution():
    ...


def test_figures_bins_size_distribution():
    ...


def test_figures_mags_GC_boxplots():
    ...


def test_figures_mags_HiC_cov_boxplots():
    ...


def test_figures_mags_SG_cov_boxplots():
    ...


def test_plot_figures():
    # Import all files and put them as the same format as
    bin_summary_file = "tests_data/outdir/bin_summary.txt"
    contigs_data_file = "tests_data/outdir/contig_data_final.txt"
    out_dir = "tests_data/outdir"
    bin_summary = pd.read_csv(bin_summary_file, sep="\t")
    bin_summary.index = bin_summary['Unnamed: 0']
    bin_summary = bin_summary.to_dict(orient="index")
    contigs_data = pd.read_csv(contigs_data_file, sep="\t")
    mtf.plot_figures(out_dir, contigs_data, bin_summary, 10000)


def test_reindex_df():
    ...
