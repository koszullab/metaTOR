#!/usr/bin/env python3
# coding: utf-8

"""Abstract command classes for metaTOR

This module contains all classes related to metaTOR commands:
    - network
    - partition
    - validation
    - pipeline
    - qc
    - contactmap
    - scaffold
    - pairs    
    - host
    - binning

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

import logging
import multiprocessing as mp
import os
import shutil
import time
import metator.align as mta
import metator.contact_map as mtc
import metator.host as mth
import metator.io as mio
import metator.log as mtl
import metator.mge as mtm
import metator.quality_check as mtq
import metator.network as mtn
import metator.partition as mtp
import metator.scaffold as mts
import metator.validation as mtv
from docopt import docopt
from functools import partial
from metator.log import logger
from metator.version import __version__
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
        network --forward=STR --assembly=FILE [--reverse=STR]
        [--aligner=bowtie2] [--aligner-mode=normal] [--depth=FILE] [--edge=0]
        [--enzyme=STR] [--normalization=empirical_hit] [--no-clean-up]
        [--outdir=DIR] [--min-quality=30] [--self-contacts] [--start=fastq]
        [--threads=1] [--tmpdir=DIR]

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
        -b, --aligner=STR       Aligner algorithm to use. Either "bwa" or
                                "bowtie2". [Default: bowtie2]
        -B, --aligner-mode=STR  Mode of alignment from hicstuff. Either normal,
                                iterative or cutsite. [Default: normal]
        -d, --depth=FILE        The depth.txt file from the shotgun reads used
                                to made the assembly computed by
                                jgi_summarize_bam_contig_depths from metabat2
                                pipeline.
        -E, --edge=INT          Distance of the edge region in base pair on the
                                contigs where the mapping reads are not
                                considered as inter contigs. [Default: 0]
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
        -T, --tmpdir=DIR       Temporary directory. Default to current
                                directory. [Default: ./tmp]
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tmpdir"]:
            tmp_dir = mio.generate_temp_dir("./tmp")
        else:
            tmp_dir = self.args["--tmpdir"]
            os.makedirs(tmp_dir, exist_ok=True)

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        os.makedirs(self.args["--outdir"], exist_ok=True)

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(self.args["--outdir"], (f"metator_network_{now}.log"))
        generate_log_header(log_file, cmd="network", args=self.args)

        # Transform integer variables as integer.
        min_qual = int(self.args["--min-quality"])
        edge = int(self.args["--edge"])

        # Defined boolean variables:
        self_contacts = self.args["--self-contacts"]

        # Check if forward and reverse arguments are given:
        if (
            self.args["--start"] == "fastq"
            or (
                self.args["--start"] == "bam"
                and self.args["--aligner"] == "bowtie2"
            )
        ) and not self.args["--reverse"]:
            logger.error(
                "Forward and reverse arguments are necessary for fastq with %s start and %s aligner.",
                self.args["--start"],
                self.args["--aligner"],
            )
            raise ValueError

        # Check aligner.
        if self.args["--aligner"] not in ["bowtie2", "bwa"]:
            logger.error('Aligner should be either "bowtie2" or "bwa".')
            raise ValueError

        # Check aligner mode.
        if self.args["--aligner-mode"] not in [
            "normal",
            "iterative",
            "cutsite",
        ]:
            logger.error(
                'Aligner mode should be either "normal", "iterative" or "cutsite".'
            )
            raise ValueError
        if (
            self.args["--aligner-mode"] == "cutsite"
            and not self.args["--enzyme"]
        ):
            logger.warning(
                "cutsite mode required an enzyme. Iterative mode will be use instead."
            )
            self.args["--aligner-mode"] = "iterative"

        # Check correct algorithm value.

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
        index = mio.check_fasta_index(assembly, mode=self.args["--aligner"])
        if index is None:
            if mio.check_is_fasta(assembly):
                fasta = assembly
                # If start at bam could skip the index generation.
                if self.args["--start"] == "fastq":
                    index = mio.generate_fasta_index(
                        fasta, self.args["--aligner"], tmp_dir
                    )
            else:
                logger.error(
                    "Please give as assembly argument a %s index or a fasta.",
                    self.args["--aligner"],
                )
                raise ValueError
        else:
            fasta = mio.retrieve_fasta(index, self.args["--aligner"], tmp_dir)

        # Print information of teh workflow:
        logger.info("Aligner algorithm: %s", self.args["--aligner"])
        logger.info("Enzyme: %s", self.args["--enzyme"])
        logger.info("Normalization: %s", self.args["--normalization"])
        logger.info("Minimum mapping quality: %s", self.args["--min-quality"])

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
                self.args["--aligner"],
                self.args["--aligner-mode"],
                min_qual,
                self.args["--start"],
                self.args["--depth"],
                self.args["--enzyme"],
                self.args["--outdir"],
                tmp_dir,
                self.args["--threads"],
            )

        # Build the network
        mtn.alignment_to_contacts(
            alignment_files,
            contig_data,
            edge,
            hit_data,
            self.args["--outdir"],
            "network.txt",
            "contig_data_network.txt",
            tmp_dir,
            self.args["--threads"],
            self.args["--normalization"],
            self_contacts,
        )

        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(tmp_dir)

        generate_log_footer(log_file)


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
        [--no-clean-up] [--outdir=DIR] [--overlap=80] [--prefix=STR]
        [--res-param=1.0] [--size=500000] [--threads=1] [--tmpdir=DIR]

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
        -p, --prefix=STR        Prefix to use for fasta files. By default just
                                'metator_', otherwise 'STR_metator_'.
        -r, --res-param=FLOAT   Resolution paramter to use for Leiden algorithm.
                                [Default: 1.0]
        -s, --size=INT          Threshold size to keep bins in base pair.
                                [Default: 500000]
        -t, --threads=INT       Number of parallel threads allocated for the
                                partition. [Default: 1]
        -T, --tmpdir=DIR       Temporary directory. Default to current
                                directory. [Default: ./tmp]
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tmpdir"]:
            tmp_dir = mio.generate_temp_dir("./tmp")
        else:
            tmp_dir = self.args["--tmpdir"]
            os.makedirs(tmp_dir, exist_ok=True)

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
                    "%s already existed. Remove directory or use -F argument to overwrite it.",
                    fasta_dir,
                )
                raise ValueError

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(self.args["--outdir"], (f"metator_partition_{now}.log"))
        generate_log_header(log_file, cmd="partition", args=self.args)

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

        # Create prefix.
        if not self.args["--prefix"]:
            prefix = "metator"
        else:
            prefix = f"{self.args['--prefix']}_metator"

        # Partition the network
        _clustering_matrix_file, _contigs_data_file = mtp.partition(
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
            tmp_dir,
            threads,
            prefix,
        )

        # Delete pyfastx index:
        os.remove(self.args["--assembly"] + ".fxi")
        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(tmp_dir)

        generate_log_footer(log_file)


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
        [--no-clean-up] [--outdir=DIR] [--overlap=90] [--prefix=STR]
        [--res-param=1.0] [--size=500000] [--threads=1] [--tmpdir=DIR]

    options:
        -a, --assembly=FILE     The path to the assembly fasta file used to do
                                the alignment.
        -A, --algorithm=STR     Algorithm to use. Either "louvain", "leiden" or
                                "spinglass". [Default: louvain]
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
        -p, --prefix=STR        Prefix to use for fasta files. By default just
                                'metator_', otherwise 'STR_metator_'.
        -r, --res-param=FLOAT   Resolution paramter to use for Leiden
                                algorithm. [Default: 1.0]
        -s, --size=INT          Threshold size to keep bins in base pair.
                                [Default: 500000]
        -t, --threads=INT       Number of parallel threads allocated for the
                                partition. [Default: 1]
        -T, --tmpdir=DIR       Temporary directory. Default to current
                                directory. [Default: ./tmp]
    """

    # Launch checkM to evaluate the completion and the contamination. If asked
    # rerun Louvain to try to reduce the contamination, rerun checkM if the
    # contamination decrease without a huge decrease of the completion keep the
    # new bins. Otherwise go back to the old state.
    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tmpdir"]:
            tmp_dir = mio.generate_temp_dir("./tmp")
        else:
            tmp_dir = self.args["--tmpdir"]
            os.makedirs(tmp_dir, exist_ok=True)

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
                    "%s already existed. Remove directory or use -F argument to overwrite it.",
                    recursive_fasta_dir,
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
                    "%s already existed. Remove directory or use -F argument to overwrite it.",
                    final_fasta_dir,
                )
                raise ValueError

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(
            self.args["--outdir"], (f"metator_validation_{now}.log")
        )
        generate_log_header(log_file, cmd="validation", args=self.args)

        # Transform numeric variable as numeric
        iterations = int(self.args["--iterations"])
        size = int(self.args["--size"])
        threads = int(self.args["--threads"])
        overlapping_parameter = int(self.args["--overlap"]) / 100
        resolution_parameter = float(self.args["--res-param"])

        # Check correct algorithm value
        if self.args["--algorithm"] not in ["louvain", "leiden", "spinglass"]:
            logger.error(
                'algorithm should be either "louvain", "leiden" or "spinglass.'
            )
            raise ValueError

        # Create prefix.
        if not self.args["--prefix"]:
            prefix = "metator"
        else:
            prefix = f"{self.args['--prefix']}_metator"

        _clustering_matrix_file = mtv.recursive_decontamination(
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
            tmp_dir,
            threads,
            prefix,
        )

        # Delete pyfastx index:
        os.remove(self.args["--assembly"] + ".fxi")
        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(tmp_dir)

        generate_log_footer(log_file)


class Pipeline(AbstractCommand):
    """Launch the full metator pipeline

    Partition the contigs from teh  assembly in bins from the metaHiC reads.

    It's possible to start from the fastq, the bam, the pair or the network
    files. It's will also possible to ask or not to run a validation step which
    will decontaminate the bins when it's necessary. However as it's the
    critical step for memory usage (~40G), it's possible to skip these step.

    usage:
        pipeline --assembly=FILE [--forward=STR] [--reverse=STR]
        [--algorithm=louvain] [--aligner=bowtie2] [--aligner-mode=normal]
        [--cluster-matrix] [--depth=FILE] [--edge=0] [--enzyme=STR] [--force]
        [--iterations=100] [--rec-iter=10] [--junctions=NNNNN] [--no-clean-up]
        [--normalization=empirical_hit] [--outdir=DIR] [--overlap=80]
        [--prefix=STR] [--rec-overlap=90]  [--min-quality=30] [--res-param=1.0]
        [--size=500000] [--start=fastq] [--scaffold] [--threads=1]
        [--tmpdir=DIR]

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
        -A, --algorithm=STR     Algorithm to use. Either "louvain", "leiden" or
                                "spinglass". If spinglass is chosen, the first
                                partition will be done using louvain
                                [Default: louvain]
        -b, --aligner=STR       Aligner algorithm to use. Either "bwa" or
                                "bowtie2". [Default: bowtie2]
        -B, --aligner-mode=STR  Mode of alignment from hicstuff. Either normal,
                                iterative or cutsite. [Default: normal]
        -C, --cluster-matrix    If enabled, save the clustering matrix.
        -d, --depth=FILE        The depth.txt file from the shotgun reads used
                                to made the assembly computed by
                                jgi_summarize_bam_contig_depths from metabat2
                                pipeline.
        -E, --edge=INT          Distance of the edge region in base pair on the
                                contigs where the mapping reads are not
                                considered as inter contigs. [Default: 0]
        -e, --enzyme=STR        The list of restriction enzyme used to digest
                                the contigs separated by a comma. Example:
                                HpaII,MluCI.
        -F, --force             If enable, would remove directory of overlapping
                                bins in the output directory.
        -i, --iterations=INT    Number of iterations of Louvain for the
                                partition step. [Default: 100]
        -j, --rec-iter=INT      Number of iterations of Louvain for the
                                recursive step. [Default: 10]
        -J, --junctions=STR     Sequences to use as junction between contigs.
                                [Default: NNNNN]
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
        -p, --prefix=STR        Prefix to use for fasta files. By default just
                                'metator_', otherwise 'STR_metator_'.
        -P, --rec-overlap=INT   Hamming distance threshold to use to merge bins
                                for the recursive step (percentage).
                                [Default: 90]
        -q, --min-quality=INT   Threshold of quality necessary to considered a
                                read properly aligned. [Default: 30]
        -r, --res-param=FLOAT   Resolution parameter to use for Leiden
                                algorithm. [Default: 1.0]
        -s, --size=INT          Threshold size to keep bins in base pair.
                                [Default: 500000]
        -S, --start=STR         Start stage of the pipeline. Either fastq, bam,
                                or pair. [Default: fastq]
        -t, --threads=INT       Number of parallel threads allocated for the
                                alignement. [Default: 1]
        -T, --tmpdir=DIR        Temporary directory. Default to current
                                directory. [Default: ./tmp]
        -v, --scaffold          If enables, it will scaffold the genomes at the
                                end.
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tmpdir"]:
            tmp_dir = mio.generate_temp_dir("./tmp")
        else:
            tmp_dir = self.args["--tmpdir"]
            os.makedirs(tmp_dir, exist_ok=True)

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
                    "%s already existed. Remove directory or use -F argument to overwrite it.",
                    overlapping_fasta_dir,
                )
                raise ValueError

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(self.args["--outdir"], (f"metator_pipeline_{now}.log"))
        generate_log_header(log_file, cmd="pipeline", args=self.args)

        # Define variable
        min_qual = int(self.args["--min-quality"])
        iterations = int(self.args["--iterations"])
        edge = int(self.args["--edge"])
        recursive_iterations = int(self.args["--rec-iter"])
        overlapping_parameter = int(self.args["--overlap"]) / 100
        recursive_overlapping_parameter = int(self.args["--rec-overlap"]) / 100
        size = int(self.args["--size"])
        threads = int(self.args["--threads"])
        resolution_parameter = float(self.args["--res-param"])

        # Check aligner.
        if self.args["--aligner"] not in ["bowtie2", "bwa"]:
            logger.error('Aligner should be either "bowtie2" or "bwa".')
            raise ValueError

        # Check aligner mode.
        if self.args["--aligner-mode"] not in [
            "normal",
            "iterative",
            "cutsite",
        ]:
            logger.error(
                'Aligner mode should be either "normal", "iterative" or "cutsite".'
            )
            raise ValueError
        if (
            self.args["--aligner-mode"] == "cutsite"
            and not self.args["--enzyme"]
        ):
            logger.warning(
                "cutsite mode required an enzyme. Iterative mode will be use instead."
            )
            self.args["--aligner-mode"] = "iterative"

        # Check correct algorithm value.
        if self.args["--algorithm"] not in ["louvain", "leiden", "spinglass"]:
            logger.error(
                'algorithm should be either "louvain", "leiden" or "spinglass.'
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

        # Sanity check for validation
        recursive_fasta_dir = join(self.args["--outdir"], "recursive_bin")
        if not exists(recursive_fasta_dir):
            os.makedirs(recursive_fasta_dir)
        else:
            if self.args["--force"]:
                shutil.rmtree(recursive_fasta_dir)
                os.makedirs(recursive_fasta_dir)
            else:
                logger.error(
                    "%s already existed. Remove directory or use -F argument to overwrite it",
                    recursive_fasta_dir,
                )
                raise ValueError
        final_fasta_dir = join(self.args["--outdir"], "final_bin_unscaffold")
        if not exists(final_fasta_dir):
            os.makedirs(final_fasta_dir)
        else:
            if self.args["--force"]:
                shutil.rmtree(final_fasta_dir)
                os.makedirs(final_fasta_dir)
            else:
                logger.error(
                    "%s already existed. Remove directory or use -F argument to overwrite it.",
                    final_fasta_dir,
                )
                raise ValueError

        # Sanity check for scaffolding
        if self.args["--scaffold"]:
            scaffold_fasta_dir = join(self.args["--outdir"], "scaffold_bin")
            if not exists(scaffold_fasta_dir):
                os.makedirs(scaffold_fasta_dir)
            else:
                if self.args["--force"]:
                    shutil.rmtree(scaffold_fasta_dir)
                    os.makedirs(scaffold_fasta_dir)
                else:
                    logger.error(
                        "%s already existed. Remove directory or use -F argument to overwrite it.",
                        scaffold_fasta_dir,
                    )
                    raise ValueError

        # Create prefix.
        if not self.args["--prefix"]:
            prefix = "metator"
        else:
            prefix = f"{self.args['--prefix']}_metator"

        # Manage start point.
        if self.args["--start"] == "fastq":
            start = 1
        elif self.args["--start"] == "bam":
            start = 2
        elif self.args["--start"] == "pair":
            start = 3
        else:
            logger.error("Start argument should be 'fastq', 'bam', or 'pair'.")
            raise ValueError

        # Check if forward and reverse reads are given for fastq and bam start.
        if (
            self.args["--start"] == "fastq"
            or (
                self.args["--start"] == "bam"
                and self.args["--aligner"] == "bowtie2"
            )
        ) and not self.args["--reverse"]:
            logger.error(
                "Forward and reverse arguments are necessary for fastq with %s start and %s aligner.",
                self.args["--start"],
                self.args["--aligner"],
            )
            raise ValueError

        # Print information of the workflow:
        if start == 1:
            logger.info(f"Minimum mapping quality: {min_qual}")
        if start <= 2:
            logger.info(f"Enzyme: {self.args['--enzyme']}")
            logger.info(f"Normalization: {self.args['--normalization']}")
        logger.info(f"Aligner algorithm: {self.args['--aligner']}")
        logger.info(f"Partition algorithm: {self.args['--algorithm']}")
        logger.info(f"Partition iterations: {iterations}")
        logger.info(f"Overlapping parameter: {overlapping_parameter}")
        logger.info(f"Recursive partition iterations: {recursive_iterations}")
        logger.info(
            f"Recursive overlapping parameter: {recursive_overlapping_parameter}\n"
        )

        # Extract index and genome file
        assembly = self.args["--assembly"]
        # Check what is the reference. If a fasta is given build the index. If a
        # bowtie2 index is given, retreive the fasta.
        index = mio.check_fasta_index(assembly, mode=self.args["--aligner"])
        if index is None:
            if mio.check_is_fasta(assembly):
                fasta = assembly
                if start == 1:
                    index = mio.generate_fasta_index(
                        fasta, self.args["--aligner"], tmp_dir
                    )
            else:
                logger.error(
                    "Please give as assembly argument a bowtie2 index or a fasta."
                )
                raise ValueError
        else:
            fasta = mio.retrieve_fasta(index, self.args["--aligner"], tmp_dir)

        # Run the whole workflow
        if start <= 2:
            # Align pair-end reads with bowtie2
            alignment_files, contig_data, hit_data = mta.get_contact_pairs(
                self.args["--forward"],
                self.args["--reverse"],
                index,
                fasta,
                self.args["--aligner"],
                self.args["--aligner-mode"],
                min_qual,
                self.args["--start"],
                self.args["--depth"],
                self.args["--enzyme"],
                self.args["--outdir"],
                tmp_dir,
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
            edge,
            hit_data,
            self.args["--outdir"],
            "network.txt",
            "contig_data_network.txt",
            tmp_dir,
            self.args["--threads"],
            self.args["--normalization"],
            False,
        )

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
            tmp_dir,
            threads,
            prefix,
        )

        # remove contig_data_network if not an input
        if start <= 2:
            contig_data_network_file = join(
                self.args["--outdir"], "contig_data_network.txt"
            )
            os.remove(contig_data_network_file)

        # Launch validation.
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
            tmp_dir,
            threads,
            prefix,
        )

        if self.args["--cluster-matrix"]:
            # Make the sum with the partition clustering matrix and save it.
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

        # Launch the scaffold
        if self.args["--scaffold"]:
            bin_summary = mio.read_bin_summary(
                join(self.args["--outdir"], "bin_summary.txt")
            )
            task = partial(
                mts.parallel_scaffold,
                final_fasta_dir=final_fasta_dir,
                alignment_files=alignment_files,
                scaffold_fasta_dir=scaffold_fasta_dir,
                junctions=self.args["--junctions"],
            )
            pool = mp.Pool(processes=int(self.args["--threads"]))
            pool.map(task, bin_summary.index)

        # Delete pyfastx index:
        os.remove(fasta + ".fxi")
        # Delete the temporary folder.
        if not self.args["--no-clean-up"]:
            shutil.rmtree(tmp_dir)
            if start < 3:
                for pair_file in alignment_files:
                    os.remove(pair_file)

        generate_log_footer(log_file)


class Qc(AbstractCommand):
    """Generates some quality check on the output of metator.

    The quality check is done on the high quality MAGs and only contigs witha
    size bigger than 100kb. The goal is to assess HiC quality metrics on genomes
    we kinda know to have a estimation of the quality of the library.

    If the plot are set it will return some plot of the event distribution on
    the high quality contigs.

    usage:
        qc --assembly=FILE --enzyme=STR [--bin-summary=FILE]
        [--contig-data=FILE] [--metator-dir=DIR] [--no-clean-up] [--outdir=DIR]
        [--prefix=STR] [--plot] [--threshold=STR] [--tmpdir=DIR] <pairsfile>...

    arguments:
        pairsfile           File(s) containing pairs information.

    options:
        -a, --assembly=FILE     Path to the fasta file containing the contigs of
                                interest. Could be the whole or the extracted
                                contigs of one bin.
        -b, --bin-summary=FILE  Path to the bin_summary.txt file from MetaTOR
                                output.
        -c, --contig-data=FILE  Path to the contig_data_final.txt file from
                                MetaTOR output.
        -e, --enzyme=STR        The list of restriction enzyme used to digest
                                the contigs separated by a comma. Example:
                                HpaII,MluCI.
        -n, --no-clean-up       If enabled intermediary files will be kept.
        -O, --metator-dir=DIR   Output directory from metator pipeline. If set
                                contigs data and bin summary will be found
                                automatically.
        -o, --outdir=DIR        Directory to save output plots and log.
        -p, --prefix=STR        Name of the sample to add on plot and files.
        -P, --plot              If enable display some plots.
        -t, --threshold=STR     Hicstuff religation and loop thresholds.
                                Two integers seperated by a coma.
        -T, --tmpdir=DIR        Temporary directory to save pairs files
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tmpdir"]:
            self.args["--tmpdir"] = "./tmp"
        tmp_dir = mio.generate_temp_dir(self.args["--tmpdir"])

        # Set output directory.
        if not self.args["--outdir"]:
            if self.args["--prefix"]:
                prefix = self.args["--prefix"]
                self.args["--outdir"] = join(".", f"metator_qc_{prefix}")
            else:
                self.args["--outdir"] = join(".", "metator_qc")
        os.makedirs(self.args["--outdir"], exist_ok=True)

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(self.args["--outdir"], (f"metator_qc_{now}.log"))
        generate_log_header(log_file, cmd="qc", args=self.args)

        # Try to find input files.
        if self.args["--metator-dir"]:
            metator_dir = self.args["--metator-dir"]
        else:
            metator_dir = "."
        if self.args["--bin-summary"]:
            bin_summary_file = self.args["--bin-summary"]
        else:
            bin_summary_file = join(metator_dir, "bin_summary.txt")
            if not exists(bin_summary_file):
                logger.error(
                    "Please give a bin summary file, or the ouput directory of metator or launch this command from the output directory of metator."
                )
                logger.error(
                    "You can also check that a bin_summary.txt is present in metator output directory"
                )
                raise FileNotFoundError
        if self.args["--contig-data"]:
            contig_data_file = self.args["--contig-data"]
        else:
            contig_data_file = join(metator_dir, "contig_data_final.txt")
            if not exists(contig_data_file):
                logger.error(
                    "Please give a contig data file, or the ouput directory of metator or launch this command from the output directory of metator."
                )
                logger.error(
                    "You can also check that a contig_data_final.txt is present in metator output directory"
                )
                raise FileNotFoundError

        # Set parameters.
        if self.args["--prefix"]:
            prefix = self.args["--prefix"]
        else:
            prefix = "metator_qc"
        if self.args["--enzyme"]:
            self.args["--enzyme"] = self.args["--enzyme"].split(",")
        if self.args["--threshold"]:
            self.args["--threshold"] = self.args["--threshold"].split(",")

        # Launch quality check
        mtq.quality_check(
            contig_data_file,
            bin_summary_file,
            self.args["--assembly"],
            self.args["<pairsfile>"],
            self.args["--outdir"],
            tmp_dir,
            prefix,
            self.args["--plot"],
            enzyme=self.args["--enzyme"],
            threshold=self.args["--threshold"],
        )

        # Delete the temporary folder.
        if not self.args["--no-clean-up"]:
            shutil.rmtree(tmp_dir)

        generate_log_footer(log_file)


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
        contactmap --assembly=FILE --contig-data=FILE --enzyme=STR --name=STR
        [--filter] [--force] [--mat-fmt=graal] [--min-size=5000]
        [--no-clean-up] [--object=final_bin] [--outdir=DIR] [--pcr-dup]
        [--tmpdir=DIR] [--threads=1] <pairsfile>...

    argument:
        pairsfile               File(s) containing pairs information.

    options:
        -a, --assembly=FILE     Path to the fasta file containing the contigs of
                                interest. Could be the whole or the extracted
                                contigs of one bin.
        -c, --contig-data=FILE  Path to the contig_data_final.txt file from
                                MetaTOR output.
        -e, --enzyme=STR        The list of restriction enzyme used to digest
                                the contigs separated by a comma. Example:
                                HpaII,MluCI.
        -n, --name=STR          Name of the metator or its numerical ID for the
                                core bin, overlapping bin or recursive bin
                                objects. Example: "NODE_1", "MetaTOR_1_0" or 8.
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

        # Enable file logging
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = join(self.args["--outdir"], (f"metator_contactmap_{now}.log"))
        generate_log_header(log_file, cmd="contactmap", args=self.args)

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
            self.args["<pairsfile>"],
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

        # Delete the temporary folder.
        if not self.args["--no-clean-up"]:
            shutil.rmtree(tmp_dir)
            # Delete pyfastx index:
            os.remove(self.args["--assembly"] + ".fxi")

        generate_log_footer(log_file)

class Scaffold(AbstractCommand):
    """Scaffold a bin from metator.

    Use the pairs to build a scaffold from the binned contigs. The scaffolding
    is based on contact between the edges of the contigs. There are normalized
    based on intra contacts of the contigs except for the small contigs where
    a low score is given.

    usage:
        scaffold --bin-name=STR --input-fasta=FILE [--junctions=NNNNN]
        [--out-fasta=FILE] [--out-frags=FILE] [--threads=1] [--threshold=0.05]
        [--window-size=5000] <pairsfile>...

    arguments:
        pairsfile               File(s) containing pairs information.

    options:
        -b, --bin-name=STR      Name of the bin to scaffold.
        -i, --input-fasta=FILE  Path to the input fasta.
        -j, --junctions=STR     Sequences to use as junction between contigs.
                                [Default: NNNNN]
        -o, --out-fasta=FILE    Path to write the output fasta. Default in the
                                current directory: '"$bin_name"_scaffolded.fa'.
        -O, --out-frags=FILE    Path to write the fragments information. Default
                                in the current directory:
                                '"$bin_name"_info_frags.txt'.
        -t, --threads=INT       Numbers of threads to allocate if pairs need to
                                be sorted. [Default: 1]
        -T, --threshold=FLOAT   Threshold score to consider an association
                                between two contigs. [Default: 0.05]
        -w, --window-size=INT   Size of the window in base pair to use as edge
                                window to scaffold. [Default: 5000]
    """

    def execute(self):

        # Generate log
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = f"metator_scaffold_{now}.log"
        generate_log_header(log_file, cmd="scaffold", args=self.args)

        # Create output files:
        if self.args["--out-fasta"] is None:
            out_fasta = join(".", f'{self.args["--bin-name"]}_scaffolded.fa')
        else:
            out_fasta = self.args["--out-fasta"]
        if self.args["--out-frags"] is None:
            out_info = join(".", f'{self.args["--bin-name"]}_info_frags.txt')
        else:
            out_info = self.args["--out-frags"]

        # Get the scaffolds
        mts.get_scaffolds(
            self.args["--bin-name"],
            self.args["--input-fasta"],
            self.args["<pairsfile>"],
            out_fasta,
            out_info,
            threshold=float(self.args["--threshold"]),
            threads=int(self.args["--threads"]),
            window_size=int(self.args["--window-size"]),
            junctions=self.args["--junctions"],
        )

        generate_log_footer(log_file)


class Pairs(AbstractCommand):
    """Sort, compress and index pairs files for faster assess to the data.

    Sort the pairs file using pairtools. Compress them using bgzip. Index them
    using pypairix.

    usage:
        pairs [--force] [--remove] [--threads=1] <pairsfile>...

    arguments:
        pairsfile           File(s) containing pairs information.

    options:
        -F, --force         Write files even if the output files already exists.
        -r, --remove        Remove the input file at the end to keep only the
                            sorted, compressed and indexed pairs file.
        -t, --threads=INT   Numbers of thread to allocate. [Default: 1]
    """

    def execute(self):

        # Generate log
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = f"metator_pairs_{now}.log"
        generate_log_header(log_file, cmd="pairs", args=self.args)

        # Iterates on pairfiles given.
        pairsfiles = self.args["<pairsfile>"]
        for pairsfile in pairsfiles:
            logger.info(f"Processing {pairsfile}...")
            # Run the sort/compress/index command.
            _ = mio.sort_pairs_pairtools(
                pairsfile,
                threads=int(self.args["--threads"]),
                remove=self.args["--remove"],
                force=self.args["--force"],
            )

        generate_log_footer(log_file)



class Host(AbstractCommand):
    """Detect host of mge annotated contigs.

    It will return an output file with the mges information from MetaTOR
    binning with the added information from the anvio binning and the detected
    bacterial MAG host of the mges.

    usage:
        host --network=FILE --binning=FILE --mges=FILE --contig-data=FILE
        [--outfile=FILE] [--threshold=0.1]

    options:
        -b, --binning=FILE      Path to the anvio binning file.
        -c, --contig-data=FILE  Path to the MetaTOR contig data file.
        -m, --mges=FILE       Path to the file with mges contigs list.
        -n, --network=FILE      Path to the network file.
        -o, --outfile=FILE      Path where to write the output file.
        -t, --threshold=FLOAT   Threshold to consider an association with a MAG.
                                [Default: 0.1]
    """

    def execute(self):

        # Generate log
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = f"metator_host_{now}.log"
        generate_log_header(log_file, cmd="host", args=self.args)

        # Defined the output file if none are given
        if not self.args["--outfile"]:
            self.args["--outfile"] = "./mges_data_host.tsv"

        # Import the files
        binning_result = mio.import_anvio_binning(self.args["--binning"])
        mges_list = mio.import_mges_contigs(self.args["--mges"])
        contig_data, mges_list_id = mio.import_contig_data_mges(
            self.args["--contig-data"], binning_result, mges_list
        )
        network = mio.import_network(self.args["--network"])

        # Run the host detection
        mth.host_detection(
            network,
            contig_data,
            mges_list,
            mges_list_id,
            self.args["--outfile"],
            self.args["--threshold"],
        )

        generate_log_footer(log_file)



class Mge(AbstractCommand):
    """Bin mges contigs.

    Bin the mge contigs in mges MAGs using the bacterial host detection from
    MetaVir and the metagenomic binning base on sequences and coverage from
    metabat2.

    The results are then checked using checkV and some plots are displayed to
    viusalize them. It will return the updated mges data, the mges fasta,
    the detailed ouputs of checKV and the plots.

    The mge fasta will contain one entry by mge MAG, with 180 "N" spacers
    between contigs.

    usage:
        mge --network=FILE --binning=FILE --mges=FILE --contigs-data=FILE --fasta=FILE
        [--checkv-db=DIR] [--depth=FILE] [--method=pairs] [--no-clean-up]
        [--outdir=DIR] [--plot] [--random] [--threads=1] [--tmpdir=DIR] 
        [--threshold-bin=0.8] [--threshold-asso=0.1] <pairsfile>...

    arguments:
        pairsfile               File(s) containing pairs information.

    options:
        -b, --binning=FILE      Path to the anvio binning file.
        -c, --contigs-data=FILE  Path to the MetaTOR contig data file.
        --checkv-db=DIR         Directory where the checkV database is stored.
                                By default the CHECKVDB environment variable is
                                used.
        -d, --depth=FILE        Path to the depth file from metabat2 script:
                                jgi_summarize_bam_contig_depths.
        -f, --fasta=FILE        Path to the fasta file with tha mge contigs
                                sequences.
        -m, --mges=FILE         Path to the file with mges contigs list.
        -M, --method=STR        Method for the binning. Either 'metabat' or
                                'pairs' [Default: pairs].
        -n, --network=FILE      Path to the network file.
        -N, --no-clean-up       If enabled, remove the temporary files.
        -o, --outdir=DIR        Path to the output directory where the output
                                will be written. Default current directory.
        -p, --plot              If enable, make summary plots.
        -r, --random            If enable, make a random binning.
        -s, --threshold-bin=FLOAT       Threshold to use for binning. 
                                        [Default: 0.8]
        -S, --threshold-asso=FLOAT      Threshold to use for association. If 
                                several MAGs have value higher than this ratio 
                                of total contatcs several association are 
                                considered. [Default: 0.1]
        -t, --threads=INT       Number of threads to use for checkV.
                                [Default: 1]
        -T, --tmpdir=DIR        Path to temporary directory. [Default: ./tmp]
    """

    def execute(self):

        # Generate log
        now = time.strftime("%Y%m%d%H%M%S")
        log_file = f"metator_mge_{now}.log"
        generate_log_header(log_file, cmd="mge", args=self.args)

        # Defined the temporary directory.
        if not self.args["--tmpdir"]:
            self.args["--tmpdir"] = "./tmp"
        tmp_dir = mio.generate_temp_dir(self.args["--tmpdir"])
        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        os.makedirs(self.args["--outdir"], exist_ok=True)

        # Set remove tmp for checkV.
        if not self.args["--no-clean-up"]:
            remove_tmp = True
        else:
            remove_tmp = False

        # Set checkV database path
        if not self.args["--checkv-db"]:
            self.args["--checkv-db"] = os.getenv("CHECKVD")

        # Sanity check
        if self.args["--method"] == "metabat" and not self.args["--depth"]:
            logger.error("Depth file is necessary if method is metabat.")
            raise ValueError

        pairs_files = self.args["<pairsfile>"]

        # Import the files
        binning_result = mio.import_anvio_binning(self.args["--binning"])
        mges_list = mio.import_mges_contigs(self.args["--mges"])
        contigs_data, mges_list_id = mio.import_contig_data_mges(
            self.args["--contigs-data"], binning_result, mges_list
        )
        network = mio.import_network(self.args["--network"])

        # Run the mges binning
        mtm.mge_binning(
            checkv_db=self.args["--checkv-db"],
            depth_file=self.args["--depth"],
            fasta_mges_contigs=self.args["--fasta"],
            network=network,
            contigs_data=contigs_data,
            mges_list_id=mges_list_id,
            out_dir=self.args["--outdir"],
            pairs_files=pairs_files,
            tmp_dir=tmp_dir,
            threshold_bin=float(self.args["--threshold-bin"]),
            threshold_asso=float(self.args["--threshold-asso"]),
            association=True,
            plot=self.args["--plot"],
            remove_tmp=remove_tmp,
            threads=int(self.args["--threads"]),
            method=self.args["--method"],
            random=self.args["--random"],
        )

        # Delete the temporary folder.
        if remove_tmp:
            shutil.rmtree(tmp_dir)
            # Delete pyfastx index:
            os.remove(self.args["--fasta"] + ".fxi")

        generate_log_footer(log_file)

def generate_log_header(log_path, cmd, args):
    mtl.set_file_handler(log_path, formatter=logging.Formatter(""))
    logger.info(f"## MetaTOR: v{__version__} log file")
    logger.info(f"## date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(
        f"## Command used: metator {cmd} {' '.join(f'{key} {value}' for key, value in args.items())}"
    )
    logger.info("\n---\n")
    mtl.set_file_handler(log_path, formatter=mtl.logfile_formatter)


def generate_log_footer(log_path):
    mtl.set_file_handler(log_path, formatter=logging.Formatter(""))
    logger.info("\n---")
    logger.info(f"## date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
