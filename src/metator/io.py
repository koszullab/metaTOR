#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Core tools to build I/O for metator

This mdoule contains all core I/O functions:
    - check_is_fasta
    - file_compression
    - generate_temp_dir
    - get_restriction_site
    - import_mges_contigs
    - process_ligation_sites
    - read_compressed
    - rename_1x_columns
    - write_bin_summary
"""

import bz2
import gzip
import io
import os
import pandas as pd
import zipfile
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
import metator.clusters as mtcl
from metator.log import logger
from os.path import join, exists
from secrets import randbits
from typing import NamedTuple


def check_is_fasta(in_file):
    """
    Function from hicstuff.io (https://github.com/koszullab/hicstuff/)
    Checks whether input file is in fasta format.

    Parameters:
    -----------
    in_file : str
        Path to the input file.

    Returns:
    --------
    bool :
        True if the input is in fasta format, False otherwise
    """
    try:
        with read_compressed(in_file) as handle:
            test = any(SeqIO.parse(handle, "fasta"))
    # Biopython >= 1.87 raises ValueError on any text before the first '>'
    # line, blank or '#' lines included. Earlier versions skipped that text
    # (1.85 and 1.86 with a warning), so text with no '>' line gave no record.
    except (FileNotFoundError, ValueError):
        test = False

    return test


# The first bytes of the compressed files that read_compressed opens.
COMPRESSION_MAGIC = {
    b"\x1f\x8b\x08": "gz",
    b"\x42\x5a\x68": "bz2",
    b"\x50\x4b\x03\x04": "zip",
}


def file_compression(filename):
    """The compression of a file, found from its first bytes, whatever its
    name, as read_compressed opens it: "gz", "bz2", "zip", or
    "uncompressed" for any other file.

    Function from hicstuff.io (https://github.com/koszullab/hicstuff/)
    """
    with open(filename, "rb") as f:
        file_start = f.read(max(len(magic) for magic in COMPRESSION_MAGIC))
    for magic, filetype in COMPRESSION_MAGIC.items():
        if file_start.startswith(magic):
            return filetype
    return "uncompressed"


def generate_temp_dir(path):
    """Temporary directory generation

    Function from hicstuff.io (https://github.com/koszullab/hicstuff/)
    Generates a temporary file with a random name at the input path.

    The name is drawn from the random numbers of the system, not from
    Python's random module, which igraph draws from too: the partition of a
    run does not depend on how many of these directories it made (one per
    command), so that the steps run one by one in one process, after one
    seed, draw the stream of pipeline.

    Parameters:
    -----------
    path : str
        The path at which the temporary directory will be created.

    Returns:
    --------
    str
        The path of the newly created temporary directory.
    """
    exist = True
    while exist:
        # Keep trying random directory names if they already exist
        directory = str(hex(randbits(32)))[2:]
        full_path = join(path, directory)
        if not exists(full_path):
            exist = False
    try:
        os.makedirs(full_path)
    except PermissionError:
        raise PermissionError(
            "The temporary directory cannot be created in {}. Make sure you have write permission.".format(path)
        ) from None
    return full_path


def get_restriction_site(enzyme):
    """Function to return a regex which corresponds to all possible restriction
    sites given a set of enzyme.

    Parameters:
    -----------
    enzyme : str
        String that contains the names of the enzyme separated by a comma.

    Returns:
    --------
    str :
        Regex that corresponds to all possible restriction sites given a set of
        enzyme.

    Examples:
    ---------
    >>> get_restriction_site('DpnII')
    'GATC'
    >>> get_restriction_site('DpnII,HinfI')
    'GA.TC|GATC'
    """

    # Split the str on the comma to separate the different enzymes.
    enzyme = enzyme.split(",")

    # Check on Biopython dictionnary the enzyme.
    rb = RestrictionBatch(enzyme)

    # Initiation:
    restriction_list = []

    # Iterates on the enzymes.
    for enz in rb:
        # Extract restriction sites and look for cut sites.
        restriction_list.append(enz.site.replace("N", "."))

    # Build the regex for all retsriction sites.
    pattern = "|".join(sorted(list(set(restriction_list))))
    return pattern


def import_mges_contigs(mges_file):
    """Import list of mges contigs.

    Parameters:
    -----------
    mges_file : str
        Path to the mges file which contains the list of mges contigs. One
        contig per line.

    Returns:
    --------
    list:
        List of the mges contigs names.
    """
    mges_list = []
    with open(mges_file, "r") as mges:
        for contig in mges:
            mges_list.append(contig.split()[0])
    return mges_list


def micomplete_results_to_dict(micomplete_file):
    """Read micomplte output file and transfrom it as a dictionnary with the
    bin name as keys and bin information as values.

    Parameters
    ----------
    micomplete_file : str
        Path to the micomplete output file.

    Returns
    -------
    dict
        Dictionnary of the output of miComplete as values and the bin id as
        keys.
    """
    # Read table.
    micomplete_summary = pd.read_csv(
        micomplete_file,
        sep="\t",
        comment="#",
        index_col=0,
        **mtcl.EXACT_FLOATS,
    ).iloc[:, :13]

    # Transform to dictionnary.
    micomplete_summary = micomplete_summary.to_dict(orient="index")

    return micomplete_summary


def read_compressed(filename):
    """Read compressed file

    Function from hicstuff.io (https://github.com/koszullab/hicstuff/)
    Opens the file in read mode with appropriate decompression algorithm.

    Parameters:
    -----------
    filename : str
        The path to the input file

    Returns:
    --------
    file-like object
        The handle to access the input file's content
    """

    # Open file with appropriate function
    comp = file_compression(filename)
    if comp == "gz":
        return gzip.open(filename, "rt")
    elif comp == "bz2":
        return bz2.open(filename, "rt")
    elif comp == "zip":
        zip_arch = zipfile.ZipFile(filename, "r")
        if len(zip_arch.namelist()) > 1:
            raise IOError("Only a single fastq file must be in the zip archive.")
        else:
            # ZipFile opens as bytes by default, using io to read as text
            zip_content = zip_arch.open(zip_arch.namelist()[0], "r")
            return io.TextIOWrapper(zip_content, encoding="utf-8")
    else:
        return open(filename, "r")


#: The columns of the contig data files that metaTOR 1.x named after bins,
#: with their names since metaTOR 2.0, which calls its groups of contigs
#: clusters. No command reads the tables of 1.x anymore: host.py renames
#: their columns (see rename_1x_columns), and the writer of the map of the
#: clusters drops those that a map made from a table of 1.x carries
#: (metator.metacool._write_back_columns).
CONTIG_DATA_1X_COLUMNS = {
    "Core_bin_ID": "Core_cluster_ID",
    "Core_bin_contigs": "Core_cluster_contigs",
    "Core_bin_size": "Core_cluster_size",
    "Overlapping_bin_ID": "Overlapping_cluster_ID",
    "Overlapping_bin_contigs": "Overlapping_cluster_contigs",
    "Overlapping_bin_size": "Overlapping_cluster_size",
    "Recursive_bin_ID": "Recursive_cluster_ID",
    "Recursive_bin_contigs": "Recursive_cluster_contigs",
    "Recursive_bin_size": "Recursive_cluster_size",
    "Final_bin": "Final_cluster",
}


def rename_1x_columns(contig_data):
    """The table of contig data contig_data, its columns named after bins by
    metaTOR 1.x renamed as metaTOR 2.0 names them (CONTIG_DATA_1X_COLUMNS).
    Kept for host.py, which calls it to read the contig data files of a 1.x
    run as those of a 2.0 run: no command reads the tables of 1.x anymore. A
    column already named as in 2.0 is kept, and the 1.x column of the same
    meaning then left as it is.

    Parameters
    ----------
    contig_data : pandas.DataFrame
        Table of contig data, as read from a contig data file.

    Returns
    -------
    pandas.DataFrame
        The table, its 1.x columns renamed, or the table itself if it has
        none.
    """
    renamed = {
        old: new for old, new in CONTIG_DATA_1X_COLUMNS.items() if old in contig_data.columns and new not in contig_data.columns
    }
    if not renamed:
        return contig_data
    logger.info(
        "Contig data of metaTOR 1.x: its columns %s are read as %s.",
        ", ".join(renamed),
        ", ".join(renamed.values()),
    )
    return contig_data.rename(columns=renamed)


# The columns of the summary of the clusters that write_bin_summary writes with
# four decimals: the abundances. Its other decimal numbers have two, and its
# whole numbers none.
ABUNDANCE_COLUMNS = ("HiC_abundance", "SG_abundance")


class SummaryLayout(NamedTuple):
    """How write_bin_summary wrote the summary of the final clusters
    (cluster_summary.txt), which the map of the clusters records with its
    values (see metator.cli.validation_entry), to write it again: a header
    line, an empty field then the columns, then a line per row, its name
    then a field per column, tab-separated.

    Attributes
    ----------
    rows : list of str
        The names of the clusters, in the order of the lines.
    columns : list of str
        The columns, in order, after the name of the cluster.
    formats : dict
        The format of each column, by name, as printf writes it: "%.2f" (a
        number with two decimals), "%.4f" (four: the abundances), "%d" (a
        whole number) or "%s" (text as it is: True or False for a boolean).
        A missing value is an empty field.
    """

    rows: list
    columns: list
    formats: dict


def _summary_format(column, values):
    """The format of a column of the summary of the clusters, as
    write_bin_summary writes it (see SummaryLayout): values, the column
    before it is written."""
    if column in ABUNDANCE_COLUMNS:
        return "%.4f"
    if pd.api.types.is_bool_dtype(values):
        return "%s"
    if pd.api.types.is_integer_dtype(values):
        return "%d"
    if pd.api.types.is_float_dtype(values):
        return "%.2f"
    return "%s"


def write_bin_summary(bin_summary, bin_summary_file):
    """Function to write the bin summary from dictionnary to table text file.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the output of the miComplete of the bins.
    bin_summary_file : str
        Path to the output file to write the summary informations of the bins.

    Returns
    -------
    SummaryLayout
        How the summary is written: its rows, its columns and their format,
        which the map of the clusters records (see metator.cli.validation_entry).
    """

    # Transform dictionnary to pandas DataFrame.
    bin_summary = pd.DataFrame.from_dict(bin_summary, orient="index")
    formats = {str(column): _summary_format(column, bin_summary[column]) for column in bin_summary.columns}

    # Change float format of the coverage.
    bin_summary["HiC_abundance"] = bin_summary["HiC_abundance"].map(lambda x: "%.4f" % x)
    try:
        bin_summary["SG_abundance"] = bin_summary["SG_abundance"].map(lambda x: "%.4f" % x)
    except KeyError:
        pass

    # Write the file.
    bin_summary.to_csv(bin_summary_file, sep="\t", float_format="%.2f")
    return SummaryLayout([str(row) for row in bin_summary.index], list(formats), formats)
