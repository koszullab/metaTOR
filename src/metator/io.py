#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Core tools to build I/O for metaTOR

This mdoule contains all core I/O functions:
    - check_checkm
    - check_fasta_index
    - check_is_fasta
    - check_louvain_cpp
    - check_pypairix
    - check_pairtools
    - generate_fasta_index
    - generate_temp_dir
    - get_pairs_data
    - get_restriction_site
    - import_anvio_binning
    - import_contig_data_mges
    - import_network
    - import_mges_contigs
    - process_ligation_sites
    - read_bin_summary
    - read_compressed
    - read_contig_data
    - read_results_checkm
    - retrieve_fasta
    - sort_pairs
    - sort_pairs_pairtools
    - write_bin_summary
    - write_mge_data
"""

import bz2
import gzip
import io
import networkx as nx
import numpy as np
import os
import pandas as pd
import pathlib
import pypairix
import pairtools
import re
import subprocess as sp
import zipfile
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from metator.log import logger
from os.path import join, exists, isfile
from random import getrandbits
from packaging import version
from packaging.version import Version


def check_checkm():
    """
    Function to test if CheckM is in the path.

    Returns:
    --------
    bool:
        True if checkM found in the path, False otherwise.
    """
    try:
        checkm = sp.check_output("checkm", stderr=sp.STDOUT, shell=True)
    except sp.CalledProcessError:
        logger.error("Cannot find 'checkm' in your path please install it or add it in your path.")
        return False
    return True


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
    if mode == "bowtie2":
        # Bowtie2 should have 6 index files
        bt2_idx_files = list(ref.parent.glob("{}*bt2*".format(ref.name)))
        index = None if len(bt2_idx_files) < 6 else bt2_idx_files
    elif mode == "bwa":
        refdir = str(ref.parent)
        refdir_files = os.listdir(refdir)
        bwa_idx_files = [join(refdir, f) for f in refdir_files if re.search(r".*\.(sa|pac|bwt|ann|amb)$", f)]
        index = None if len(bwa_idx_files) < 5 else bwa_idx_files
    else:
        index = [ref]
    if index is not None:
        # Convert the PosixPath objects to strings and get the longest common
        # prefix to obtain index basename (without the dot)
        index = os.path.commonprefix(list(map(str, index))).strip(".")
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
        with open(in_file, "r") as handle:
            test = any(SeqIO.parse(handle, "fasta"))
    except FileNotFoundError:
        test = False

    return test


def check_louvain_cpp(louvain_path):
    """Function to check is the Louvain functions are callable.

    Parameters:
    -----------
    louvain_path : str
        Path of the directory where the Louvain functions are.

    Returns:
    --------
    bool:
        Boolean value describing either the Louvain functions are callable or
        not.
    """

    # Look for the path to call the functions:
    louvain = join(louvain_path, "louvain")
    convert = join(louvain_path, "convert")
    hierarchy = join(louvain_path, "hierarchy")

    # Check convert:
    try:
        convert = sp.check_output(f"{convert} --help", stderr=sp.STDOUT, shell=True)
    except sp.CalledProcessError:
        logger.error("Cannot find the 'convert' function from Louvain path.")
        return False
        raise ImportError

    # Check louvain:
    try:
        louvain = sp.check_output(f"{louvain} --help", stderr=sp.STDOUT, shell=True)
    except sp.CalledProcessError:
        logger.error("Cannot find the 'louvain' function from Louvain path.")
        return False
        raise ImportError

    # Check hierarchy:
    try:
        hierarchy = sp.check_output(f"{hierarchy} --help", stderr=sp.STDOUT, shell=True)
    except sp.CalledProcessError:
        logger.error("Cannot find the convert_net function from Louvain path.")
        return False
        raise ImportError

    return True


def check_pypairix():
    """
    Function to test if pypairix is available.

    Returns:
    --------
    bool:
        True if pypairix is available.
    """
    try:
        v = version.parse(pypairix.__version__)
    except AttributeError:
        logger.error("Cannot find 'pypairix' installed.")
        raise AttributeError
    return True


def check_pairtools():
    """
    Function to test if pairtools is in the path.

    Returns:
    --------
    bool:
        True if pairtools found in the path, False otherwise.
    """
    try:
        pairtools = sp.check_output("pairtools", stderr=sp.STDOUT, shell=True)
    except sp.CalledProcessError:
        logger.error("Cannot find 'pairtools' in your path please install it or add it in your path.")
        raise ImportError
        return False
    return True


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
        directory = str(hex(getrandbits(32)))[2:]
        full_path = join(path, directory)
        if not exists(full_path):
            exist = False
    try:
        os.makedirs(full_path)
    except PermissionError:
        raise PermissionError(
            "The temporary directory cannot be created in {}. " "Make sure you have write permission.".format(path)
        )
    return full_path


def get_pairs_data(pairfile, threads=1, remove=False, force=False):
    """Extract pairs data from pypairix indexed pairs file. If no pypairix
    indexed found, sort pairs files using pairtools executable.

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
            pairfile = sort_pairs_pairtools(pairfile, threads, remove, force)
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
    contig_data = pd.read_csv(contig_data_file, sep="\t")
    contig_data["Binned"] = False
    contig_data["Final_bin"] = "ND"
    contig_data["MGE"] = False
    mges_list_id = []
    for i in contig_data.index:
        if contig_data.loc[i, "Name"] in mges_list:
            contig_data.loc[i, "MGE"] = True
            mges_list_id.append(contig_data.index[i])
        try:
            contig_data.loc[i, "Final_bin"] = binning_result[contig_data.loc[i, "Name"]]
            contig_data.loc[i, "Binned"] = True
        except KeyError:
            continue
    return contig_data, mges_list_id


def import_network(network_file):
    """Import MetaTOR network file.

    Parameters:
    -----------
    network_file : str
        Path to the network file to import.

    Returns:
    --------
    networkx.classes.graph.Graph:
        Network as networkx class.
    """
    network = nx.read_edgelist(network_file, nodetype=int, data=(("weight", float),))
    return network


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

    # Standard header bytes for diff compression formats
    comp_bytes = {
        b"\x1f\x8b\x08": "gz",
        b"\x42\x5a\x68": "bz2",
        b"\x50\x4b\x03\x04": "zip",
    }

    max_len = max(len(x) for x in comp_bytes)

    def file_type(filename):
        """Guess file type

        Compare header bytes with those in the file and return type.
        """
        with open(filename, "rb") as f:
            file_start = f.read(max_len)
        for magic, filetype in comp_bytes.items():
            if file_start.startswith(magic):
                return filetype
        return "uncompressed"

    # Open file with appropriate function
    comp = file_type(filename)
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
    """Read bin summary file from metator pipeline.

    Parameters
    ----------
    contig_data_file : str
        Path to the contig data file.

    Returns
    -------
    pandas.DataFrame
        Table from contig data with contig name as index.
    """
    data = pd.read_csv(contig_data_file, sep="\t")
    data = data.set_index("Name", drop=False)
    return data


def read_results_checkm(checkm_file, checkm_taxonomy_file):
    """Function to transform the output summary file of checkm into a
    dictionnary.

    Parameters:
    -----------
    checkm_file : str
        Path to the summary output file of CheckM.

    Returns:
    --------
    dict:
        Dictionnary of the output of checkm with binID as keys and with three
        values lineage, completness and contamination.
    """

    # Create an empty dictionnary
    checkm_summary = dict()

    # Read the checkm summary file.
    with open(checkm_file, "r") as checkm_lines:
        for line in checkm_lines:
            # Only keep informative lines which start with a space.
            if line[0] == " ":
                line = line.split()
                checkm_summary[line[0]] = {
                    "lineage": line[1],
                    "completness": line[6],
                    "contamination": line[7],
                    "size": line[9],
                    "contigs": line[12],
                    "N50": line[14],
                    "longest_contig": line[18],
                    "GC": line[19],
                    "coding_density": line[21],
                    "taxonomy": "-",
                }

    # Read the taxonomy file.
    with open(checkm_taxonomy_file, "r") as checkm_lines:
        for line in checkm_lines:
            # Only keep informative lines which start with a space.
            if line[0] == " ":
                line = line.split()
                checkm_summary[line[0]]["taxonomy"] = line[3]

    # remove the bin header
    checkm_summary.pop("Bin")

    return checkm_summary


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
        ValueError("Sparse matrix must be in coo format")
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
                "GNU sort version is {0} but >8.23 is required for parallel "
                "sort. Sorting on a single thread.".format(".".join(map(str, sort_ver)))
            )
            parallel_ok = False
    # BSD sort has a different format and will throw error upon parsing. It does
    # not support parallel processes anyway.
    except ValueError:
        logger.warning("Using BSD sort instead of GNU sort, sorting on a single thread.")
        parallel_ok = False

    # Sort pairs and append to file.
    with open(out_file, "a") as output:
        grep_proc = sp.Popen(["grep", "-v", "^#", in_file], stdout=sp.PIPE)
        sort_cmd = ["sort", "-S %s" % buffer] + ["-k1,1V", "-k2,2V"]
        if tmp_dir is not None:
            sort_cmd.append("--temporary-directory={0}".format(tmp_dir))
        if parallel_ok:
            sort_cmd.append("--parallel={0}".format(threads))
        sort_proc = sp.Popen(sort_cmd, stdin=grep_proc.stdout, stdout=output)
        sort_proc.communicate()


def sort_pairs_pairtools(pairfile, threads=1, remove=False, force=False):
    """Sort pairs files using pairtools executable. Pairix only works with
    compressed pair files. So we use bgzip to compress them.

    Parameters
    ----------
    pairfile : str
        Path to the pairfile to sort, compress and index.
    threads : int
        Number of threads to use. [Default: 1]
    remove : bool
        If set to true, it will remove the unsorted pair file. [Default: False]
    force : bool
        If set overwrite existing files. [Default: False]

    Returns
    -------
    str :
        Path to the sorted and indexed pair file.
    """
    # Extract basename of the file.
    basename = os.path.splitext(pairfile)[0]

    # Test if pypairix and pairtools are installed and in the path.
    _ = check_pypairix()
    _ = check_pairtools()

    # Set the force parameter and delete files or raise an error accodringly.
    if force:
        force = " -f"
        if os.path.isfile(f"{basename}_sorted.pairs"):
            os.remove(f"{basename}_sorted.pairs")
    else:
        force = ""
        if os.path.isfile(f"{basename}_sorted.pairs") or os.path.isfile(f"{basename}_sorted.pairs.gz"):
            logger.error(
                f"The {basename}_sorted.pairs exists. Do not overwrite existing, use --force to overwrite or use another location."
            )
            raise ValueError

    # Sort pairs using pairtools.
    cmd = f"set -eu ; pairtools sort {pairfile} --nproc {threads} -o {basename}_sorted.pairs"
    if Version(pairtools.__version__) >= Version("1.1.0"):
        logger.debug("pairtools version >= 1.1.0. Use new options.")
        cmd = cmd + " --c1 chr1 --c2 chr2 --p1 pos1 --p2 pos2 --pt strand1"
    process = sp.Popen(cmd, shell=True)
    _out, _err = process.communicate()
    # Compressed pairs.
    cmd = f"set -eu ; bgzip {basename}_sorted.pairs -@ {threads}{force}"
    process = sp.Popen(cmd, shell=True)
    _out, _err = process.communicate()
    # Indexed pairs.
    force_pypairix = 1 if force else 0
    pypairix.build_index(f"{basename}_sorted.pairs.gz", force=force_pypairix)

    # Remove original pairfile if remove setup.
    if remove:
        os.remove(pairfile)
    return f"{basename}_sorted.pairs.gz"


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
