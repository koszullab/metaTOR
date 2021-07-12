#!/usr/bin/env python3
# coding: utf-8

"""Abstract command classes for metaTOR

This module contains all classes related to metaTOR commands:
    - network
    - partition
    - pipeline
    - validation

Note
----
Structure based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
: https://github.com/rgreinho/docopt-subcommands-example
abignaud, 20201118

Raises
------
NotImplementedError
    Will be raised if AbstractCommand is called for some reason instead of one
    of its children.
"""

import os
import shutil
import time
import metator.align as mta
import metator.io as mio
import metator.log as mtl
import metator.network as mtn
import metator.partition as mtp
import metator.validation as mtv
import metator.contact_map as mtc
from docopt import docopt
from metator.log import logger
from os.path import exists, dirname, join
from scipy.sparse import save_npz, load_npz


class AbstractCommand:
    """Abstract base command class

    Base class for the commands from which other metaTOR commands derive.
    """

    def __init__(self, command_args, global_args):
        """Initialize the commands"""
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError

    def check_output_path(self, path, force=False):
        """Throws error if the output file exists. Create required file tree
        otherwise.
        """
        # Get complete output filename and prevent overwriting unless force is
        # enabled
        if not force and os.path.exists(path):
            raise IOError(
                "Output file already exists. Use --force to overwrite"
            )
        if dirname(path):
            os.makedirs(dirname(path), exist_ok=True)


class Network(AbstractCommand):
    """Generation of network command

    Generate metaHiC contigs network from fastq reads or bam files.Generates a
    network file (in edgelist form) from either fastq files or bam files. For
    both starting input forward and reverse reads or alignment need to be in
    seprated files with the same read names. In the output network contigs are
    the network nodes and the edges are the contact counts.

    If multiple fastq are given will align them seprately and return a hit file
    with a number of hit for each sample.

    Note: the network is in a strict barebone form so that it can be reused and
    imported quickly into other applications etc. Verbose information about
    every single node in the network is written on a 'contig data' file.

    usage:
        network --forward=STR --assembly=FILE [--reverse=STR] [--depth=FILE]
        [--enzyme=STR] [--normalization=empirical_hit] [--no-clean-up]
        [--outdir=DIR] [--min-quality=30] [--self-contacts] [--start=fastq]
        [--threads=1] [--tempdir=DIR]

    options:
        -1, --forward=STR       Fastq file or list of Fastq separated by a comma
                                containing the forward reads to be aligned or
                                their corresponding bam or pairs files.
        -2, --reverse=STR       Fastq file or list of Fastq separated by a comma
                                containing the reverse reads to be aligned or
                                their corresponding bam files. Forward and
                                reverse reads need to have the same identifier
                                (read names). If start is set to pair, no
                                argument is necessary.
        -a, --assembly=FILE     The initial assembly path acting as the
                                alignment file's reference genome or the
                                basename of the bowtie2 index.
        -d, --depth=FILE        The depth.txt file from the shotgun reads used
                                to made the assembly computed by
                                jgi_summarize_bam_contig_depths from metabat2
                                pipeline.
        -e, --enzyme=STR        The list of restriction enzyme used to digest
                                the contigs separated by a comma. Example:
                                HpaII,MluCI.
        -n, --normalization=STR If None, do not normalized the count of a
                                contact by the geometric mean of the coverage of
                                the contigs. Otherwise it's the type of
                                normalization. 6 values are possible "None",
                                "abundance", "length", "RS", "empirical_hit",
                                "theoritical_hit". [Default: empirical_hit]
        -N, --no-clean-up       Do not remove temporary files.
        -o, --outdir=DIR        The output directory to write the bam files the
                                network and contig data into. Default: current
                                directory.
        -q, --min-quality=INT   Threshold of quality necessary to considered a
                                read properly aligned. [Default: 30]
        -s, --self-contacts     If enabled, count alignments between a contig
                                and itself (intracontigs contacts).
        -S, --start=STR         Start stage of the pipeline. Either "fastq",
                                "bam", or "pair". [Default: fastq]
        -t, --threads=INT       Number of parallel threads allocated for the
                                alignement. [Default: 1]
        -T, --tempdir=DIR       Temporary directory. Default to current
                                directory. [Default: ./tmp]
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        os.makedirs(self.args["--outdir"], exist_ok=True)

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(
            self.args["--outdir"], ("metator_network_" + now + ".log")
        )
        mtl.set_file_handler(log_file)

        # Transform integer variables as integer.
        min_qual = int(self.args["--min-quality"])

        # Defined boolean variables:
        self_contacts = self.args["--self-contacts"]

        # Check if forward and reverse arguments are given:
        if self.args["--start"] != "pair":
            if not self.args["--forward"] or not self.args["--reverse"]:
                logger.error(
                    "Forward and reverse arguments are necessary for fastq or bam start."
                )
                raise ValueError

        # Check if normalization in the list of possible normalization.
        list_normalization = [
            "None",
            "abundance",
            "length",
            "RS",
            "empirical_hit",
            "theoritical_hit",
        ]
        if self.args["--normalization"] not in list_normalization:
            logger.error(
                'Normalization should be among this list: "None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"'
            )
            raise ValueError
        enzyme_required = ["RS", "theoritical_hit"]
        if (
            self.args["--normalization"] in enzyme_required
            and not self.args["--enzyme"]
        ):
            logger.error(
                'For "RS" and "theoritical_hit" normalization, enzyme is required.'
            )
            raise ValueError
        depth_required = ["abundance", "theoritical_hit"]
        if (
            self.args["--normalization"] in depth_required
            and not self.args["--depth"]
        ):
            logger.error(
                'For "abundance" and "theoritical_hit" normalization, depth is required.'
            )
            raise ValueError
        if self.args["--start"] not in ["fastq", "bam", "pair", "network"]:
            logger.error(
                "Start argument should be 'fastq', 'bam', 'pair' or 'network'."
            )
            raise ValueError
        # Extract index and genome file
        assembly = self.args["--assembly"]
        # Check what is the reference. If a fasta is given build the index. If a
        # bowtie2 index is given, retreive the fasta.
        index = mio.check_fasta_index(assembly, mode="bowtie2")
        if index is None:
            if mio.check_is_fasta(assembly):
                fasta = assembly
                # If start at bam could skip the index generation.
                if self.args["--start"] == "fastq":
                    index = mio.generate_fasta_index(fasta, temp_directory)
            else:
                logger.error(
                    "Please give as assembly argument a bowtie2 index or a fasta."
                )
                raise ValueError
        else:
            fasta = mio.retrieve_fasta(index, temp_directory)

        # Print information of teh workflow:
        logger.info("Enzyme: {0}".format(self.args["--enzyme"]))
        logger.info("Normalization: {0}".format(self.args["--normalization"]))
        logger.info(
            "Minimum mapping quality: {0}".format(self.args["--min-quality"])
        )

        # Do not align if pair start
        if self.args["--start"] == "pair":
            alignment_files = self.args["--forward"].split(",")
            nb_alignment = len(alignment_files)
            contig_data, hit_data = mtn.create_contig_data(
                fasta,
                nb_alignment,
                self.args["--depth"],
                self.args["--enzyme"],
            )

        else:
            # Align pair-end reads with bowtie2
            alignment_files, contig_data, hit_data = mta.get_contact_pairs(
                self.args["--forward"],
                self.args["--reverse"],
                index,
                fasta,
                min_qual,
                self.args["--start"],
                self.args["--depth"],
                self.args["--enzyme"],
                self.args["--outdir"],
                temp_directory,
                self.args["--threads"],
            )

        # Build the network
        mtn.alignment_to_contacts(
            alignment_files,
            contig_data,
            hit_data,
            self.args["--outdir"],
            "network.txt",
            "contig_data_network.txt",
            temp_directory,
            self.args["--threads"],
            self.args["--normalization"],
            self_contacts,
        )

        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)


class Partition(AbstractCommand):
    """Partition the network using Louvain algorithm

    Partition the network file using iteratively the Louvain or Leiden
    algorithm. Then looks for 'cores' bins that are constituted by the contigs
    which are always in the same cluster at each iterations. Then using the
    Hamming distance from these cores bins, merege the bins with a small Hamming
    distance (close to 1), bigger than the percentage (overlapping parameter)
    given.

    It will also update the file to integrate the bins information of the
    contigs.

    Furthermore, both Leiden and Louvain algorithm are available here. However,
    the benchmark made show that here the Louvain algorithm have better
    performance and is faster on seawater and gut metagenomic samples.

    Note that the Louvain or Leiden software are not, in the strictest sense,
    necessary. Any program that assigns a node to a bin, does so non
    deterministically and solely outputs a list in the form: 'node_id bin_id'
    could be plugged instead.

    usage:
        partition  --assembly=FILE --contigs=FILE --network=FILE
        [--algorithm=louvain] [--cluster-matrix] [--force] [--iterations=100]
        [--no-clean-up] [--outdir=DIR] [--overlap=80] [--res-param=1.0]
        [--size=500000] [--threads=1] [--tempdir=DIR]

    options:
        -a, --assembly=FILE     The path to the assembly fasta file used to do
                                the alignment.
        -A, --algorithm=STR     Either "louvain" or "leiden", algorithm to use
                                to partition the network. [Default: louvain]
        -c, --contigs=FILE      The path to the tsv file containing the data of
                                the contigs (ID, Name, Length, GC content, Hit,
                                Coverage, Restriction Site).
        -C, --cluster-matrix    If enabled, save the clustering matrix.
        -F, --force             If enabled, would remove directory of
                                overlapping bins in the output directory.
        -i, --iterations=INT    Number of iterations of Louvain. [Default: 100]
        -n, --network=FILE      Path to the file containing the network
                                information from the meta HiC experiment compute
                                in network function previously.
        -N, --no-clean-up       Do not remove temporary files.
        -o, --outdir=DIR        Path to the directory to write the output.
                                Default to current directory. [Default: ./]
        -O, --overlap=INT       Hamming distance threshold to use to merge bins
                                (percentage). [Default: 80]
        -r, --res-param=FLOAT   Resolution paramter to use for Leiden algorithm.
                                [Default: 1.0]
        -s, --size=INT          Threshold size to keep bins in base pair.
                                [Default: 500000]
        -t, --threads=INT       Number of parallel threads allocated for the
                                partition. [Default: 1]
        -T, --tempdir=DIR       Temporary directory. Default to current
                                directory. [Default: ./tmp]
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        os.makedirs(self.args["--outdir"], exist_ok=True)
        fasta_dir = join(self.args["--outdir"], "overlapping_bin")
        if not exists(fasta_dir):
            os.makedirs(fasta_dir)
        else:
            if self.args["--force"]:
                shutil.rmtree(fasta_dir)
                os.makedirs(fasta_dir)
            else:
                logger.error(
                    "{0} already existed. Remove directory or use -F argument to overwrite it.".format(
                        fasta_dir
                    )
                )
                raise ValueError

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(
            self.args["--outdir"], ("metator_partition_" + now + ".log")
        )
        mtl.set_file_handler(log_file)

        # Transform numeric variable as numeric
        iterations = int(self.args["--iterations"])
        overlapping_parameter = int(self.args["--overlap"]) / 100
        size = int(self.args["--size"])
        threads = int(self.args["--threads"])
        resolution_parameter = float(self.args["--res-param"])

        # Check correct algorithm value
        if self.args["--algorithm"] not in ["louvain", "leiden"]:
            logger.error('algorithm should be either "louvain" or "leiden"')
            raise ValueError

        # Partition the network
        clustering_matrix_file, contigs_data_file = mtp.partition(
            self.args["--algorithm"],
            self.args["--assembly"],
            self.args["--cluster-matrix"],
            self.args["--contigs"],
            iterations,
            self.args["--network"],
            self.args["--outdir"],
            fasta_dir,
            overlapping_parameter,
            resolution_parameter,
            size,
            temp_directory,
            threads,
        )

        # Delete pyfastx index:
        os.remove(self.args["--assembly"] + ".fxi")
        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)


class Validation(AbstractCommand):
    """Use CheckM to validate the bins.

    Recursively decontaminate the contaminated bins evaluated by checkM. The
    script returns new decontaminated fasta and summary files of the
    decontamination.

    Only bins with more than 50% completion and 5% contamination are subject to
    the recursive step. If the recursive step gave worst results than the first
    (decrease of the completion with no decrease of the contamination), it will
    keep the original bin.

    usage:
        validation --assembly=FILE --contigs=FILE --fasta=DIR --network=FILE
        [--algorithm=louvain] [--cluster-matrix] [--force] [--iterations=10]
        [--no-clean-up] [--outdir=DIR] [--overlap=90] [--res-param=1.0]
        [--size=500000] [--threads=1] [--tempdir=DIR]

    options:
        -a, --assembly=FILE     The path to the assembly fasta file used to do
                                the alignment.
        -A, --algorithm=STR     Algorithm to use. Either "louvain" or "leiden".
                                [Default: louvain]
        -c, --contigs=FILE      The path to the file containing the data of the
                                contigs from the partition step (13 columns).
        -C, --cluster-matrix    If enabled, save the clustering matrix.
        -f, --fasta=DIR         Path to the directory containing the input fasta
                                files of the bins to decontaminate.
        -F, --force             If enable, would remove directory of recursive
                                bins in the output directory.
        -i, --iterations=INT    Number of recursive iterations of Louvain.
                                [Default: 10]
        -n, --network=FILE      Path to the file containing the network
                                information from the meta HiC experiment compute
                                in network function previously.
        -N, --no-clean-up       Do not remove temporary files.
        -o, --outdir=DIR        Path to the directory to write the output.
                                Default to current directory. [Default: ./]
        -O, --overlap=INT       Hamming distance threshold to use to merge bins
                                (percentage). [Default: 90]
        -r, --res-param=FLOAT   Resolution paramter to use for Leiden
                                algorithm. [Default: 1.0]
        -s, --size=INT          Threshold size to keep bins in base pair.
                                [Default: 500000]
        -t, --threads=INT       Number of parallel threads allocated for the
                                partition. [Default: 1]
        -T, --tempdir=DIR       Temporary directory. Default to current
                                directory. [Default: ./tmp]
    """

    # Launch checkM to evaluate the completion and the contamination. If asked
    # rerun Louvain to try to reduce the contamination, rerun checkM if the
    # contamination decrease without a huge decrease of the completion keep the
    # new bins. Otherwise go back to the old state.
    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        os.makedirs(self.args["--outdir"], exist_ok=True)
        recursive_fasta_dir = join(self.args["--outdir"], "recursive_bin")
        final_fasta_dir = join(self.args["--outdir"], "final_bin")
        if not exists(recursive_fasta_dir):
            os.makedirs(recursive_fasta_dir)
        else:
            if self.args["--force"]:
                shutil.rmtree(recursive_fasta_dir)
                os.makedirs(recursive_fasta_dir)
            else:
                logger.error(
                    "{0} already existed. Remove directory or use -F argument to overwrite it.".format(
                        recursive_fasta_dir
                    )
                )
                raise ValueError
        if not exists(final_fasta_dir):
            os.makedirs(final_fasta_dir)
        else:
            if self.args["--force"]:
                shutil.rmtree(final_fasta_dir)
                os.makedirs(final_fasta_dir)
            else:
                logger.error(
                    "{0} already existed. Remove directory or use -F argument to overwrite it.".format(
                        final_fasta_dir
                    )
                )
                raise ValueError

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(
            self.args["--outdir"], ("metator_validation_" + now + ".log")
        )
        mtl.set_file_handler(log_file)

        # Transform numeric variable as numeric
        iterations = int(self.args["--iterations"])
        size = int(self.args["--size"])
        threads = int(self.args["--threads"])
        overlapping_parameter = int(self.args["--overlap"]) / 100
        resolution_parameter = float(self.args["--res-param"])

        # Check checkM availability
        if not mio.check_checkm():
            logger.error(
                "CheckM is not in the path. Could not make the iterations"
            )
            raise NameError

        # Check correct algorithm value
        if self.args["--algorithm"] not in ["louvain", "leiden"]:
            logger.error('algorithm should be either "louvain" or "leiden"')
            raise ValueError

        clustering_matrix_file = mtv.recursive_decontamination(
            self.args["--algorithm"],
            self.args["--assembly"],
            self.args["--cluster-matrix"],
            self.args["--contigs"],
            final_fasta_dir,
            self.args["--fasta"],
            iterations,
            self.args["--network"],
            self.args["--outdir"],
            overlapping_parameter,
            recursive_fasta_dir,
            resolution_parameter,
            size,
            temp_directory,
            threads,
        )

        # Delete pyfastx index:
        os.remove(self.args["--assembly"] + ".fxi")
        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)


class Pipeline(AbstractCommand):
    """Launch the full metator pipeline

    Partition the contigs from teh  assembly in bins from the metaHiC reads.

    It's possible to start from the fastq, the bam, the pair or the network
    files. It's will also possible to ask or not to run a validation step which
    will decontaminate the bins when it's necessary. However as it's the
    critical step for memory usage (~40G), it's possible to skip these step.

    usage:
        pipeline --assembly=FILE [--forward=STR] [--reverse=STR]
        [--algorithm=louvain] [--cluster-matrix] [--contigs=FILE] [--depth=FILE]
        [--enzyme=STR] [--force] [--iterations=100] [--rec-iter=10]
        [--network=FILE] [--no-clean-up] [--normalization=empirical_hit]
        [--outdir=DIR] [--overlap=80] [--rec-overlap=90]  [--min-quality=30]
        [--res-param=1.0] [--size=500000] [--start=fastq] [--threads=1]
        [--tempdir=DIR] [--skip-validation]

    options:
        -1, --forward=STR       Fastq file or list of Fastq separated by a comma
                                containing the forward reads to be aligned or
                                their corresponding bam or pairs files.
        -2, --reverse=STR       Fastq file or list of Fastq separated by a comma
                                containing the reverse reads to be aligned or
                                their corresponding bam files. Forward and
                                reverse reads need to have the same identifier
                                (read names). If pair start is used no argument
                                is necessary.
        -a, --assembly=FILE     The initial assembly path acting as the
                                alignment file's reference genome or the
                                basename of the bowtie2 index.
        -A, --algorithm=STR     Algorithm to use. Either "louvain" or "leiden".
                                [Default: louvain]
        -c, --contigs=FILE      The path to the file containing the data ofthe
                                contigs (ID, Name, Length, GC content, Hit,
                                Coverage, Restriction site).
        -C, --cluster-matrix    If enabled, save the clustering matrix.
        -d, --depth=FILE        The depth.txt file from the shotgun reads used
                                to made the assembly computed by
                                jgi_summarize_bam_contig_depths from metabat2
                                pipeline.
        -e, --enzyme=STR        The list of restriction enzyme used to digest
                                the contigs separated by a comma. Example:
                                HpaII,MluCI.
        -F, --force             If enable, would remove directory of overlapping
                                bins in the output directory.
        -i, --iterations=INT    Number of iterations of Louvain for the
                                partition step. [Default: 100]
        -j, --rec-iter=INT      Number of iterations of Louvain for the
                                recursive step. [Default: 10]
        -n, --network=FILE      Path to the file containing the network
                                information from the meta HiC experiment compute
                                in network function previously.
        -N, --no-clean-up       Do not remove temporary files.
        -m, --normalization=STR If None, do not normalized the count of a
                                contact by the geometric mean of the coverage of
                                the contigs. Otherwise it's the type of
                                normalization. 6 values are possible None,
                                abundance, length, RS, empirical_hit,
                                theoritical_hit. [Default: empirical_hit]
        -o, --outdir=DIR        The output directory to write the bam files the
                                network and contig data into. Default: current
                                directory.
        -O, --overlap=INT       Hamming distance threshold to use to merge bins
                                for the first step (percentage). [Default: 80]
        -P, --rec-overlap=INT   Hamming distance threshold to use to merge bins
                                for the recursive step (percentage).
                                [Default: 90]
        -q, --min-quality=INT   Threshold of quality necessary to considered a
                                read properly aligned. [Default: 30]
        -r, --res-param=FLOAT   Resolution paramter to use for Leiden
                                algorithm. [Default: 1.0]
        -s, --size=INT          Threshold size to keep bins in base pair.
                                [Default: 500000]
        -S, --start=STR         Start stage of the pipeline. Either fastq, bam
                                pair, or network. [Default: fastq]
        -t, --threads=INT       Number of parallel threads allocated for the
                                alignement. [Default: 1]
        -T, --tempdir=DIR       Temporary directory. Default to current
                                directory. [Default: ./tmp]
        -v, --skip-validation   If  enables do not do the validation step which
                                have an high memory usage (checkM ~ 40G)
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        os.makedirs(self.args["--outdir"], exist_ok=True)
        overlapping_fasta_dir = join(self.args["--outdir"], "overlapping_bin")
        if not exists(overlapping_fasta_dir):
            os.makedirs(overlapping_fasta_dir)
        else:
            if self.args["--force"]:
                shutil.rmtree(overlapping_fasta_dir)
                os.makedirs(overlapping_fasta_dir)
            else:
                print(self.args["--force"])
                logger.error(
                    "{0} already existed. Remove directory or use -F argument to overwrite it.".format(
                        overlapping_fasta_dir
                    )
                )
                raise ValueError

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(self.args["--outdir"], ("metator_" + now + ".log"))
        mtl.set_file_handler(log_file)

        # Define variable
        min_qual = int(self.args["--min-quality"])
        iterations = int(self.args["--iterations"])
        recursive_iterations = int(self.args["--rec-iter"])
        overlapping_parameter = int(self.args["--overlap"]) / 100
        recursive_overlapping_parameter = int(self.args["--rec-overlap"]) / 100
        size = int(self.args["--size"])
        threads = int(self.args["--threads"])
        resolution_parameter = float(self.args["--res-param"])

        # Check correct algorithm value
        if self.args["--algorithm"] not in ["louvain", "leiden"]:
            logger.error('algorithm should be either "louvain" or "leiden"')
            raise ValueError

        # Check if normalization in the list of possible normalization.
        list_normalization = [
            "None",
            "abundance",
            "length",
            "RS",
            "empirical_hit",
            "theoritical_hit",
        ]
        if self.args["--normalization"] not in list_normalization:
            logger.error(
                'Normalization should be among this list: "None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"'
            )
            raise ValueError
        enzyme_required = ["RS", "theoritical_hit"]
        if (
            self.args["--normalization"] in enzyme_required
            and not self.args["--enzyme"]
        ):
            logger.error(
                'For "RS" and "theoritical_hit" normalization, enzyme is required.'
            )
            raise ValueError
        depth_required = ["abundance", "theoritical_hit"]
        if (
            self.args["--normalization"] in depth_required
            and not self.args["--depth"]
        ):
            logger.error(
                'For "abundance" and "theoritical_hit" normalization, depth is required.'
            )
            raise ValueError

        # Sanity check for validation
        if not self.args["--skip-validation"]:
            recursive_fasta_dir = join(self.args["--outdir"], "recursive_bin")
            if not exists(recursive_fasta_dir):
                os.makedirs(recursive_fasta_dir)
            else:
                if self.args["--force"]:
                    shutil.rmtree(recursive_fasta_dir)
                    os.makedirs(recursive_fasta_dir)
                else:
                    logger.error(
                        "{0} already existed. Remove directory or use -F argument to overwrite it.".format(
                            recursive_fasta_dir
                        )
                    )
                    raise ValueError
            final_fasta_dir = join(self.args["--outdir"], "final_bin")
            if not exists(final_fasta_dir):
                os.makedirs(final_fasta_dir)
            else:
                if self.args["--force"]:
                    shutil.rmtree(final_fasta_dir)
                    os.makedirs(final_fasta_dir)
                else:
                    logger.error(
                        "{0} already existed. Remove directory or use -F argument to overwrite it.".format(
                            final_fasta_dir
                        )
                    )
                    raise ValueError

            # Check checkM availability
            if not mio.check_checkm():
                logger.error(
                    "CheckM is not in the path. Could not make the iterations"
                )
                raise NameError

        # Manage start point.
        if self.args["--start"] == "fastq":
            start = 1
        elif self.args["--start"] == "bam":
            start = 2
        elif self.args["--start"] == "pair":
            start = 3
        elif self.args["--start"] == "network":
            start = 4
        else:
            logger.error(
                "Start argument should be 'fastq', 'bam', 'pair' or 'network'."
            )
            raise ValueError

        # Check if forward and reverse reads are given for fastq and bam start.
        if start <= 2:
            if not self.args["--forward"] or not self.args["--reverse"]:
                logger.error(
                    "Forward and reverse arguments are necessary for fastq or bam start."
                )
                raise ValueError

        # Print information of the workflow:
        if start == 1:
            logger.info("Minimum mapping quality: {0}".format(min_qual))
        if start <= 2:
            logger.info("Enzyme: {0}".format(self.args["--enzyme"]))
            logger.info(
                "Normalization: {0}".format(self.args["--normalization"])
            )
        logger.info("Partition algorithm: {0}".format(self.args["--algorithm"]))
        logger.info("Partition iterations: {0}".format(iterations))
        logger.info("Overlapping parameter: {0}".format(overlapping_parameter))
        if not self.args["--skip-validation"]:
            logger.info(
                "Recursive partition iterations: {0}".format(
                    recursive_iterations
                )
            )
            logger.info(
                "Recursive overlapping parameter: {0}".format(
                    recursive_overlapping_parameter
                )
            )

        # Extract index and genome file
        assembly = self.args["--assembly"]
        # Check what is the reference. If a fasta is given build the index. If a
        # bowtie2 index is given, retreive the fasta.
        index = mio.check_fasta_index(assembly, mode="bowtie2")
        if index is None:
            if mio.check_is_fasta(assembly):
                fasta = assembly
                if start == 1:
                    index = mio.generate_fasta_index(fasta, temp_directory)
            else:
                logger.error(
                    "Please give as assembly argument a bowtie2 index or a fasta."
                )
                raise ValueError
        else:
            fasta = mio.retrieve_fasta(index, temp_directory)

        # Run the whole workflow
        if start <= 3:
            if start <= 2:
                # Align pair-end reads with bowtie2
                alignment_files, contig_data, hit_data = mta.get_contact_pairs(
                    self.args["--forward"],
                    self.args["--reverse"],
                    index,
                    fasta,
                    min_qual,
                    self.args["--start"],
                    self.args["--depth"],
                    self.args["--enzyme"],
                    self.args["--outdir"],
                    temp_directory,
                    self.args["--threads"],
                )
            else:
                alignment_files = self.args["--forward"].split(",")
                nb_alignment = len(alignment_files)
                contig_data, hit_data = mtn.create_contig_data(
                    fasta,
                    nb_alignment,
                    self.args["--depth"],
                    self.args["--enzyme"],
                )
            # Build the network
            network_file, contigs_data_file = mtn.alignment_to_contacts(
                alignment_files,
                contig_data,
                hit_data,
                self.args["--outdir"],
                "network.txt",
                "contig_data_network.txt",
                temp_directory,
                self.args["--threads"],
                self.args["--normalization"],
                False,
            )
        else:
            contigs_data_file = self.args["--contigs"]
            network_file = self.args["--network"]

        # Partition the network
        clustering_matrix_partition_file, contigs_data_file = mtp.partition(
            self.args["--algorithm"],
            fasta,
            self.args["--cluster-matrix"],
            contigs_data_file,
            iterations,
            network_file,
            self.args["--outdir"],
            overlapping_fasta_dir,
            overlapping_parameter,
            resolution_parameter,
            size,
            temp_directory,
            threads,
        )

        # remove contig_data_network if not an input
        if start <= 2:
            contig_data_network_file = join(
                self.args["--outdir"], "contig_data_network.txt"
            )
            os.remove(contig_data_network_file)

        # Launch validation if desired.
        if not self.args["--skip-validation"]:
            clustering_matrix_recursive_file = mtv.recursive_decontamination(
                self.args["--algorithm"],
                fasta,
                self.args["--cluster-matrix"],
                contigs_data_file,
                final_fasta_dir,
                overlapping_fasta_dir,
                recursive_iterations,
                network_file,
                self.args["--outdir"],
                recursive_overlapping_parameter,
                recursive_fasta_dir,
                resolution_parameter,
                size,
                temp_directory,
                threads,
            )

            if self.args["--cluster-matrix"]:
                # Make the sum with the partiton clustering matrix and save it.
                clustering_matrix = load_npz(
                    clustering_matrix_partition_file + ".npz"
                )
                clustering_matrix_recursive = load_npz(
                    clustering_matrix_recursive_file + ".npz"
                )
                clustering_matrix = (
                    (clustering_matrix + clustering_matrix_recursive) / 2
                ).tocoo()
                clustering_matrix_file = join(
                    self.args["--outdir"], "clustering_matrix"
                )
                save_npz(clustering_matrix_file, clustering_matrix)

            # Remove contig_data_partition file
            contig_data_partition_file = join(
                self.args["--outdir"], "contig_data_partition.txt"
            )
            os.remove(contig_data_partition_file)

        # Delete pyfastx index:
        os.remove(fasta + ".fxi")
        # Delete the temporary folder.
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)


class Contactmap(AbstractCommand):
    """Generates a HiC contact map of a MetaTOR object from the pairs files, the
    contig data file and a fasta file containing the contigs of interest.

    Generates the Hi-C matrix of one MetaTOR object form the pair alignment file
    of metaTOR. MetaTOR object integrated are contigs, and core, overlapping,
    recursive, final or personalized bins. For the personalized bins you should
    add a column in the contig_data_final.txt file with a header called "Other".

    Do not display the matrix, you have to run another pipeline to display it
    such hicstuff view or cooler show depending on the metrix format chosen.
    It's also possible to try to scaffold the bin using instagraal pipeline
    using the output.

    usage:
        contactmap --assembly=FILE --contig-data=FILE --enzyme=STR --name=STR --pairs=FILE
        [--filter] [--force] [--mat-fmt=graal] [--min-size=5000]
        [--no-clean-up] [--object=final_bin] [--outdir=DIR] [--pcr-dup]
        [--tmpdir=DIR] [--threads=1]

    options:
        -a, --assembly=FILE     Path to the fasta file containing the contigs of
                                interest. Could be the whole or the extracted
                                contigs of one bin.
        -c, --contig-data=FILE  Path to the contig_data_final.txt file form
                                MetaTOR output.
        -e, --enzyme=STR        The list of restriction enzyme used to digest
                                the contigs separated by a comma. Example:
                                HpaII,MluCI.
        -n, --name=STR          Name of the metator or its numerical ID for the
                                core bin, overlapping bin or recursive bin
                                objects. Example: "NODE_1", "MetaTOR_1_0" or 8.
        -p, --pairs=FILE        Path of the ".pairs" file. If more than one is
                                given, files should be separated by a comma.
        -D, --pcr-dup,          Filter out PCR duplicates based on read
                                positions.
        -f, --filter            Filter out spurious 3C events (loops and uncuts)
                                using hicstuff filter. For more informations,
                                see Cournac et al. BMC Genomics, 2012.
        -F, --force             Write files even if the output files already
                                exists.
        -m, --mat-fmt=STR       The format of the output sparse matrix. Can be
                                "bg2" for 2D Bedgraph format, "cool" for
                                Mirnylab's cooler software, or "graal" for
                                graal-compatible plain text COO format.
                                [default: graal]
        -N, --no-clean-up       If enabled intermediary files will be kept.
        -o, --outdir=DIR        Output directory. Default creates a new
                                directory contact_map in the current directory.
        -O, --object=STR        Type of metator object used to build a contact
                                map. Either "contig", "core_bin",
                                "overlapping_bin", "recursive_bin", "final_bin"
                                or "other". [Default: final_bin]
        -s, --min-size=INT      Minimum size threshold to consider contigs.
                                [Default: 5000]
        -t, --threads=INT       Number of threads to allocate. [Default: 1]
        -T, --tmpdir=DIR        Directory for storing intermediary files and
                                temporary files. Default creates a "tmp" folder
                                in the current directory.
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tmpdir"]:
            self.args["--tmpdir"] = "./tmp"
        tmp_dir = mio.generate_temp_dir(self.args["--tmpdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = join(
                ".", "contact_map_" + self.args["--name"]
            )
        os.makedirs(self.args["--outdir"], exist_ok=True)

        mtc.generate_contact_map(
            self.args["--assembly"],
            self.args["--contig-data"],
            self.args["--enzyme"],
            self.args["--name"],
            self.args["--pairs"],
            self.args["--outdir"],
            tmp_dir,
            self.args["--filter"],
            self.args["--force"],
            self.args["--mat-fmt"],
            self.args["--object"],
            int(self.args["--min-size"]),
            self.args["--pcr-dup"],
            int(self.args["--threads"]),
        )

        # Delete pyfastx index:
        os.remove(self.args["--assembly"] + ".fxi")
        # Delete the temporary folder.
        if not self.args["--no-clean-up"]:
            shutil.rmtree(tmp_dir)
