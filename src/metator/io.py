#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Core tools to build I/O for metator

This mdoule contains all core I/O functions:
    - check_fasta_index
    - check_is_fasta
    - check_pypairix
    - file_compression
    - generate_fasta_index
    - generate_temp_dir
    - get_pairs_data
    - get_restriction_site
    - import_anvio_binning
    - import_contig_data_mges
    - import_mges_contigs
    - process_ligation_sites
    - read_bin_summary
    - read_compressed
    - read_contig_data
    - read_fasta_names
    - remove_pyfastx_index
    - rename_1x_columns
    - retrieve_fasta
    - sort_and_index_pairs
    - sort_pairs
    - write_bin_summary
    - write_contigs
    - write_mge_data
"""

import bz2
import gzip
import hicstuff.io as hio
import io
import numpy as np
import os
import pandas as pd
import pathlib
import pypairix
import pysam
import shutil
import subprocess as sp
import zipfile
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from metator.log import logger
from os.path import join, exists, isfile
from secrets import randbits


def check_fasta_index(ref, mode="bowtie2"):
    """
    Function from hicstuff.io (https://github.com/koszullab/hicstuff/)
    Checks for the existence of a bowtie2 or bwa index based on the reference
    file name.

    Parameters:
    -----------
    ref : str
        Path to the reference genome.
    mode : str
        The alignment software used to build the index. bowtie2 or bwa. If any
        other value is given, the function returns the reference path.

    Returns:
    --------
    str
        The bowtie2 or bwa index basename. None if no index was found
    """
    ref = pathlib.Path(ref)
    # The files of the index of ref itself: those of another index of its
    # directory (another assembly's) are not its index.
    if mode == "bowtie2":
        # Bowtie2 has 6 index files, .bt2l for a large index.
        suffixes = ["1", "2", "3", "4", "rev.1", "rev.2"]
        index = None
        for ext in ["bt2", "bt2l"]:
            bt2_idx_files = [pathlib.Path(f"{ref}.{suffix}.{ext}") for suffix in suffixes]
            if all(path.is_file() for path in bt2_idx_files):
                index = bt2_idx_files
                break
    elif mode == "bwa":
        bwa_idx_files = [pathlib.Path(f"{ref}.{ext}") for ext in ["amb", "ann", "bwt", "pac", "sa"]]
        index = bwa_idx_files if all(path.is_file() for path in bwa_idx_files) else None
    else:
        index = [ref]
    if index is not None:
        # Convert the PosixPath objects to strings and get the longest common
        # prefix to obtain index basename (without the dot)
        # A character-wise prefix is intended: "idx.1.bt2", "idx.rev.1.bt2"... -> "idx."
        index = os.path.commonprefix(list(map(str, index))).rstrip(".")  # noqa: RUF071
    return index


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


def check_pypairix():
    """
    Function to test if pypairix is available.

    Returns:
    --------
    bool:
        True if pypairix is available.
    """
    if not hasattr(pypairix, "build_index"):
        logger.error("Cannot find 'pypairix' installed.")
        raise AttributeError("pypairix has no build_index")
    return True


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


def generate_fasta_index(fasta, aligner, outdir):
    """Generate fasta index.

    Parameters:
    -----------
    fasta : str
        Path to the fasta reference to index.
    aligner : str
        Aligner to use to build the index.
    outdir : str
        Path to the directory to write the index.

    Returns:
    --------
    str:
        Path to the bowtie2 index build
    """
    logger.info("Build index from the given fasta.\n")
    index = join(outdir, "index")
    if aligner == "bowtie2":
        cmd = "bowtie2-build -q {0} {1}".format(fasta, index)
    elif aligner == "bwa":
        cmd = "bwa index -p {1} {0}".format(fasta, index)
    process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    _out, _err = process.communicate()
    return index


def generate_temp_dir(path):
    """Temporary directory generation

    Function from hicstuff.io (https://github.com/koszullab/hicstuff/)
    Generates a temporary file with a random name at the input path.

    The name is drawn from the random numbers of the system, not from
    Python's random module, which igraph draws from too: the partition of a
    run does not depend on how many of these directories it made (one per
    command), so that the steps run one by one in one process, after one
    seed, draw the stream of pipeline. The alignment of hicstuff still
    draws from it.

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


def get_pairs_data(pairfile, threads=1, remove=False, force=False):
    """Extract pairs data from pypairix indexed pairs file. If no pypairix
    index is found, sort, compress and index the pairs file with
    sort_and_index_pairs.

    Parameters
    ----------
    pairfile : str
        Path to the pairfile to extract data.
    threads : int
        Number of threads to use to sort pairs if necessary. [Default: 1]
    remove : bool
        If set to true, it will remove the unsorted pair file. [Default: False]
    force : bool
        If set overwrite existing files. [Default: False]

    Returns
    -------
    str :
        Path to the sorted and indexed pair file.
    """
    # Check if pypairix index exists, generate it otherwise.
    try:
        pairs_data = pypairix.open(pairfile)
    except pypairix.PairixError:
        try:
            pairfile_sorted = f"{os.path.splitext(pairfile)[0]}_sorted.pairs.gz"
            pairs_data = pypairix.open(pairfile_sorted)
        except pypairix.PairixError:
            logger.warning("No pypairix index found. Build the index.")
            # Sort's temporary files go to the working directory: /tmp can be
            # small, or in memory.
            pairfile = sort_and_index_pairs(pairfile, threads, remove, force, tmp_dir=".")
            pairs_data = pypairix.open(pairfile)
    return pairs_data


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
    contig_data = rename_1x_columns(pd.read_csv(contig_data_file, sep="\t"))
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


def read_bin_summary(bin_summary_file):
    """Read bin summary file from metator pipeline.

    Parameters
    ----------
    bin_summary_file : str
        Path to the bin summary file.

    Returns
    -------
    pandas.DataFrame
        Table from bin summary with bin name as index.
    """
    return pd.read_csv(bin_summary_file, sep="\t", index_col=0)


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


def read_contig_data(contig_data_file):
    """Read a contig data file from metator pipeline, of metaTOR 1.x too
    (see rename_1x_columns).

    Parameters
    ----------
    contig_data_file : str
        Path to the contig data file.

    Returns
    -------
    pandas.DataFrame
        Table from contig data with contig name as index.
    """
    data = rename_1x_columns(pd.read_csv(contig_data_file, sep="\t"))
    data = data.set_index("Name", drop=False)
    return data


#: The columns of the contig data files that metaTOR 1.x named after bins,
#: with their names since metaTOR 2.0, which calls its groups of contigs
#: clusters.
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
    metaTOR 1.x renamed as metaTOR 2.0 names them (CONTIG_DATA_1X_COLUMNS),
    so that every command reads the contig data files of a 1.x run as those
    of a 2.0 run. A column already named as in 2.0 is kept, and the 1.x
    column of the same meaning then left as it is.

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


def read_fasta_names(fasta):
    """The names of the sequences of a fasta file, in their order, as
    Biopython names its records (SeqIO.parse, which create_contig_data reads
    the assembly with): the first word of each header line. Only the header
    lines are looked at."""
    with open(fasta, "r") as handle:
        for line in handle:
            if line.startswith(">"):
                words = line[1:].split(None, 1)
                yield words[0] if words else ""


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


def retrieve_fasta(in_file, aligner, tmpdir):
    """
    Function to retrieve fasta from the given reference file. If index is given
    retrieve it using bowtie2 inspect. Thraw an error if not a fasta or bowtie2
    index.

    Parameters:
    -----------
    in_file : str
        Path to the reference file given.
    aligner : str
        Name of the aligner used. Either 'bowtie2' or 'bwa'.
    tmpdir : str
        Path to the temp directory to write the fasta if necessary.

    Returns:
    --------
    str:
        Path to the fasta file.
    """
    if check_is_fasta(in_file):
        fasta = in_file
    else:
        if check_fasta_index(in_file, aligner):
            if aligner == "bowtie2":
                logger.info("Retrieve fasta from bowtie2 index.\n")
                fasta = join(tmpdir, "assembly.fa")
                cmd = "bowtie2-inspect {0} > {1}".format(in_file, fasta)
                process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
                _out, _err = process.communicate()
            elif aligner == "bwa":
                if isfile(in_file + ".fa"):
                    if check_is_fasta(in_file + ".fa"):
                        fasta = in_file + ".fa"
                elif isfile(in_file + ".fasta"):
                    if check_is_fasta(in_file + ".fasta"):
                        fasta = in_file + ".fasta"
                else:
                    logger.error("If you give bwa index, please make sure the fasta exists with the same prefix.")
                    raise ValueError
        else:
            logger.error("Please give as a reference a bowtie2 index or a fasta.")
            raise ValueError
    return fasta


def save_sparse_matrix(s_mat, path):
    """Save a sparse matrix

    Saves a sparse matrix object into tsv format.

    Parameters
    ----------
    s_mat : scipy.sparse.coo_matrix
        The sparse matrix to save on disk
    path : str
        File path where the matrix will be stored
    """
    if s_mat.format != "coo":
        raise ValueError("Sparse matrix must be in coo format")
    _dtype = s_mat.dtype
    sparse_arr = np.vstack([s_mat.row, s_mat.col, s_mat.data]).T

    np.savetxt(
        path,
        sparse_arr,
        header="{nrows}\t{ncols}\t{nonzero}".format(nrows=s_mat.shape[0], ncols=s_mat.shape[1], nonzero=s_mat.nnz),
        comments="",
        fmt=["%i", "%i", "%1.3f"],
        delimiter="\t",
    )


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


# The order pypairix needs: the pairs of each (chr1, chr2) block together,
# sorted by pos1.
PAIRS_SORT_KEYS = ["chr1", "chr2", "pos1", "pos2"]


def sort_and_index_pairs(pairfile, threads=1, remove=False, force=False, tmp_dir=None):
    """Sort a pairs file by chr1, chr2, pos1 and pos2 with hicstuff, compress
    it with bgzip (or pysam if bgzip is not in the path) and index it with
    pypairix, which only reads bgzip compressed pairs files.

    Contig names sort byte by byte (chr1, chr10, chr2), as with pairtools
    sort, so that each one forms a single block; positions sort as numbers.

    The header of the input is kept, with "#sorted: chr1-chr2-pos1-pos2"
    in place of its "#sorted" line. A header without one gets it after its
    "## pairs format" line.

    Parameters
    ----------
    pairfile : str
        Path to the pairfile to sort, compress and index. It can be
        compressed (gzip, bgzip, bzip2 or zip).
    threads : int
        Number of threads to use. [Default: 1]
    remove : bool
        If set to true, it will remove the unsorted pair file. [Default: False]
    force : bool
        If set overwrite existing files. [Default: False]
    tmp_dir : str
        Path to the directory where sort writes its temporary files.
        [Default: None, the default of sort: $TMPDIR, else /tmp]

    Returns
    -------
    str :
        Path to the sorted and indexed pair file, named after the input
        without its last extension ("sample_sorted.pairs.gz" for
        "sample.pairs", "sample.pairs_sorted.pairs.gz" for
        "sample.pairs.gz"), with its pypairix index (the same name with
        ".px2") next to it.
    """
    # Extract basename of the file.
    basename = os.path.splitext(pairfile)[0]
    sorted_pairs = f"{basename}_sorted.pairs"

    # Test if pypairix is installed.
    _ = check_pypairix()

    # Delete existing files or raise an error accordingly.
    if force:
        if os.path.isfile(sorted_pairs):
            os.remove(sorted_pairs)
    elif os.path.isfile(sorted_pairs) or os.path.isfile(f"{sorted_pairs}.gz"):
        logger.error(f"The {sorted_pairs} exists. Do not overwrite existing, use --force to overwrite or use another location.")
        raise ValueError

    # Sort pairs using hicstuff. It rewrites the "#sorted" line of the header,
    # but adds none.
    hio.sort_pairs(pairfile, sorted_pairs, keys=PAIRS_SORT_KEYS, tmp_dir=tmp_dir, threads=threads)
    if not any(line.startswith("#sorted") for line in hio.get_pairs_header(sorted_pairs)):
        _add_sorted_line(sorted_pairs)

    # Compress pairs.
    if shutil.which("bgzip"):
        sp.run(["bgzip", "-f", "-@", str(threads), sorted_pairs], check=True)
    else:
        pysam.tabix_compress(sorted_pairs, f"{sorted_pairs}.gz", force=True)
        os.remove(sorted_pairs)

    # Index pairs.
    pypairix.build_index(f"{sorted_pairs}.gz", force=1 if force else 0)

    # Remove original pairfile if remove setup.
    if remove:
        os.remove(pairfile)
    return f"{sorted_pairs}.gz"


def _add_sorted_line(pairs_file):
    """Add "#sorted: chr1-chr2-pos1-pos2" to the header of a pairs file:
    second, after the "## pairs format" line that starts a pairs file, or
    first if there is none.

    Parameters
    ----------
    pairs_file : str
        Path to the uncompressed pairs file to rewrite.
    """
    sorted_line = "#sorted: {}\n".format("-".join(PAIRS_SORT_KEYS)).encode()
    rewritten = f"{pairs_file}.tmp"
    with open(pairs_file, "rb") as pairs, open(rewritten, "wb") as output:
        first_line = pairs.readline()
        if first_line.startswith(b"##"):
            output.write(first_line + sorted_line)
        else:
            output.write(sorted_line + first_line)
        shutil.copyfileobj(pairs, output)
    os.replace(rewritten, pairs_file)


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
