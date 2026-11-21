#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Core tools to build I/O for metator

This mdoule contains all core I/O functions:
    - check_is_fasta
    - file_compression
    - generate_temp_dir
    - get_restriction_site
    - import_anvio_binning
    - import_contig_data_mges
    - import_mges_contigs
    - process_ligation_sites
    - read_compressed
    - remove_pyfastx_index
    - rename_1x_columns
    - sort_pairs
    - write_bin_summary
    - write_contigs
    - write_mge_data
"""

import bz2
import gzip
import io
import os
import pandas as pd
import subprocess as sp
import zipfile
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from metator.log import logger
from os.path import join, exists
from secrets import randbits


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


def import_anvio_binning(binning_file):
    """Import Anvio binning file.

    Parameters:
    -----------
    binning_file : str
        Path to the binning file from anvio to import.

    Returns:
    --------
    dict:
        Dictionnary with contig name as keys and bin name as value.
    """
    binning_result = {}
    with open(binning_file, "r") as binning:
        for line in binning:
            line = line.split()
            # Remove the split name add by anvio
            binning_result[line[0].split("_split")[0]] = line[1]
    return binning_result


def import_contig_data_mges(contig_data_file, binning_result, mges_list):
    """Import contigs data.

    Parameters:
    -----------
    contig_data_file : str
        Path to the contigs data file from MetaTOR.
    binning_result : dict
        Dictionnary with contig name as keys and bin name as value.
    mges_list : list
        List of the mges contigs names.

    Returns:
    --------
    pandas.DataFrame:
        Table with the contig information given with more column with the given
        anvio binning and the mge annotation.
    list
        List of the mge contigs ID.
    """
    contig_data = pd.read_csv(contig_data_file, sep="\t")
    contig_data["Binned"] = False
    contig_data["Final_cluster"] = "ND"
    contig_data["MGE"] = False
    mges_list_id = []
    for i in contig_data.index:
        if contig_data.loc[i, "Name"] in mges_list:
            contig_data.loc[i, "MGE"] = True
            mges_list_id.append(contig_data.index[i])
        try:
            contig_data.loc[i, "Final_cluster"] = binning_result[contig_data.loc[i, "Name"]]
            contig_data.loc[i, "Binned"] = True
        except KeyError:
            continue
    return contig_data, mges_list_id


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
#: clusters. Only host.py reads the tables of 1.x now (see
#: rename_1x_columns): every command reads those of 2.0.
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


def remove_pyfastx_index(fasta):
    """Remove the index pyfastx writes next to a fasta file, if there is
    one: pyfastx only indexes the assembly when a bin is extracted from it.

    Parameters:
    -----------
    fasta : str
        Path to the fasta file.
    """
    try:
        os.remove(f"{fasta}.fxi")
    except FileNotFoundError:
        pass


def write_contigs(fasta, names, out_file):
    """Write the records of some contigs of a fasta file, in the given order,
    each as it is in the file, as `pyfastx extract` writes them.

    Parameters:
    -----------
    fasta : pyfastx.Fasta
        The fasta file, read through its pyfastx index (<fasta>.fxi, which
        pyfastx writes next to it).
    names : iterable
        Names of the contigs. Read as text: pyfastx reads an integer key as
        the index of a record, and pandas reads a column of names such as
        1, 2... as integers.
    out_file : str
        Path to the fasta file to write.

    Raises:
    -------
    KeyError
        If a contig is not in the fasta file.
    """
    with open(out_file, "w") as out:
        for name in names:
            out.write(fasta[str(name)].raw)


def sort_pairs(in_file, out_file, tmp_dir=None, threads=1, buffer="2G"):
    """
    Adapted function from hicstuff.io (https://github.com/koszullab/hicstuff/)
    Sort a pairs file in batches using UNIX sort.

    Parameters:
    -----------
    in_file : str
        Path to the unsorted input file
    out_file : str
        Path to the sorted output file.
    keys : list of str
        list of columns to use as sort keys. Each column can be one of readID,
        chr1, pos1, chr2, pos2, frag1, frag2. Key priorities are according to
        the order in the list.
    tmp_dir : str
        Path to the directory where temporary files will be created. Defaults
        to current directory.
    threads : int
        Number of parallel sorting threads.
    buffer : str
        Buffer size used for sorting. Consists of a number and a unit.
    """
    # TODO: Write a pure python implementation to drop GNU coreutils depencency,
    # could be inspired from: https://stackoverflow.com/q/14465154/8440675

    # Check if UNIX sort version supports parallelism
    parallel_ok = True
    sort_ver = sp.Popen(["sort", "--version"], stdout=sp.PIPE)
    sort_ver = sort_ver.communicate()[0].decode().split("\n")[0].split(" ")[-1].split(".")
    # If so, specify threads, otherwise don't mention it in the command line
    try:
        sort_ver = list(map(int, sort_ver))
        if sort_ver[0] < 8 or (sort_ver[0] == 8 and sort_ver[1] < 23):
            logger.warning(
                "GNU sort version is {0} but >8.23 is required for parallel sort. Sorting on a single thread.".format(
                    ".".join(map(str, sort_ver))
                )
            )
            parallel_ok = False
    # BSD sort has a different format and will throw error upon parsing. It does
    # not support parallel processes anyway.
    except ValueError:
        logger.warning("Using BSD sort instead of GNU sort, sorting on a single thread.")
        parallel_ok = False

    # Sort pairs into a fresh output file.
    with open(out_file, "w") as output:
        grep_proc = sp.Popen(["grep", "-v", "^#", in_file], stdout=sp.PIPE)
        sort_cmd = ["sort", "-S %s" % buffer, "-k1,1V", "-k2,2V"]
        if tmp_dir is not None:
            sort_cmd.append("--temporary-directory={0}".format(tmp_dir))
        if parallel_ok:
            sort_cmd.append("--parallel={0}".format(threads))
        sort_proc = sp.Popen(sort_cmd, stdin=grep_proc.stdout, stdout=output)
        sort_proc.communicate()


def write_bin_summary(bin_summary, bin_summary_file):
    """Function to write the bin summary from dictionnary to table text file.

    Parameters:
    -----------
    bin_summary : dict
        Dictionnary with the output of the miComplete of the bins.
    bin_summary_file : str
        Path to the output file to write the summary informations of the bins.
    """

    # Transform dictionnary to pandas DataFrame.
    bin_summary = pd.DataFrame.from_dict(bin_summary, orient="index")

    # Change float format of the coverage.
    bin_summary["HiC_abundance"] = bin_summary["HiC_abundance"].map(lambda x: "%.4f" % x)
    try:
        bin_summary["SG_abundance"] = bin_summary["SG_abundance"].map(lambda x: "%.4f" % x)
    except KeyError:
        pass

    # Write the file.
    bin_summary.to_csv(bin_summary_file, sep="\t", float_format="%.2f")


def write_mge_data(mge_data, out_file):
    """Write mge binning information.

    Parameters:
    -----------
    mges_data : pandas.DataFrame
        Table with the mge binning information.
    out_file : str
        Path to write the mge data information.
    """

    # Drop the mge id column which is the concatenation of the two previous
    # ones.
    try:
        mge_data.drop("tmp", inplace=True, axis=1)
    except KeyError:
        pass

    # Write the data frame
    mge_data.to_csv(out_file, sep="\t", index=False, float_format="%.2f")
