#!/usr/bin/env python3
# coding: utf-8

"""Abstract command classes for metaTOR

This module contains all classes related to metaTOR commands:

    -alignment
    -binning
    -network 
    -partition
    -pipeline

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

from docopt import docopt
import os
import shutil
import metator.align as mta
import metator.io as mio
import metator.network as mtn


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


class Align(AbstractCommand):
    """Alignment command

    Align reads from froward and reverse fastq files. Multiple fastq could be
    given separated by commas. Re-aligned the unmapped reads using ligation
    sites to optimize the proportion of uniquely mapped reads.

    usage:
        align [--ligation-sites=STR] [--tempdir=DIR] [--threads=1]
        [--min-quality=30] --genome=FILE --out=FILE --forward
        reads_for.fastq[,reads_for2.fastq...] --reverse
        reads_rev.fastq[,reads_rev2.fastq...]

    options:
        -1, --forward=STR           Fastq file or list of Fastq separated by a
                                    comma containing the forward reads to be
                                    aligned.
        -2, --reverse=STR           Fastq file or list of Fastq separated by a
                                    comma containing the reverse reads to be
                                    aligned. Forward and reverse reads need to
                                    have the same identifier.
        -g, --genome=FILE           The genome on which to map the reads. Must
                                    be the path to the bowtie2/bwa index.
        -l, --ligation-sites=STR    Ligation site or list of ligation sites
                                    separated by a comma of the restriction
                                    enzyme(s). For example
                                    GATCGATC,GANTGATC,GANTANTC,GATCANTC
                                    for DpnII and HinfI. If no option given, it
                                    will align only once the reads. [Default:
                                    None]
        -o, --out=FILE              Path where the alignment will be written in
                                    bed2D format.
        -q, --min-quality=INT       Threshold of quality necessary to considered
                                    a read properly aligned. [Default: 30]
        -t, --threads=INT           Number of parallel threads allocated for the
                                    alignement. [Default: 1]
        -T, --tempdir=DIR           Temporary directory. [Default: ./tmp]
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Transform integer variables as integer.
        min_qual = int(self.args["--min-quality"])

        # Align pair-end reads with bowtie2
        mta.pairs_alignment(
            self.args["--forward"],
            self.args["--reverse"],
            min_qual,
            temp_directory,
            self.args["--genome"],
            self.args["--ligation-sites"],
            self.args["--out"],
            self.args["--threads"],
        )

        # Delete the temporary folder
        shutil.rmtree(temp_directory)


class Network(AbstractCommand):
    """Generation of network command

    Generates a network file (in edgelist form) from an alignment in bed2D
    format. Contigs are the network nodes and the edges are the contact counts.

    The network is in a strict barebone form so that it can be reused and
    imported quickly into other applications etc. Verbose information about
    every single node in the network is written on a 'contig data' file, by
    default called 'idx_contig_length_GC_hit_cov.txt'

    usage:
        network --genome=FILE --outdir=DIR --input=FILE [--normalized]
        [--output-file-contig-data=STR] [--output-file-network=STR]
        [--read-size=150] [--self-contacts] [--tempdir=DIR] [--threads=1]

    options:
        -g, --genome=FILE               The initial assembly path acting as the
                                        alignment file's reference genome.
        -i, --input=FILE                Path to the bed2D file used as input
        -n, --normalized                If enabled,  normalize contacts between
                                        contigs by their geometric mean
                                        coverage.
        -o, --outdir=DIR                The output directory to write the
                                        network and contig data into. Default:
                                        current directory.
        --output-file-contig-data=STR   The specific file name for the output
                                        chunk data file. [Default:
                                        'idx_contig_length_GC_hit_cov.txt']
        --output-file-network=STR       The specific file name for the output
                                        network file. Default is network.txt
        -r, --read-size=INT             Size of reads used for mapping.
                                        [Default: 150]
        -s, --self-contacts             If enabled, count alignments between a
                                        contig and itself.
        -t, --threads=INT               Number of parallel threads allocated for
                                        the alignement. [Default: 1]
        -T, --tempdir=DIR               Temporary directory. Default to current
                                        directory. [Default: ./tmp]
    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."

        if not self.args["--output-file-contig-data"]:
            self.args[
                "--output-file-contig-data"
            ] = "idx_contig_length_GC_hit_cov.txt"

        if not self.args["--output-file-network"]:
            self.args["--output-file-network"] = "network.txt"

        # Transform integer variables as integer
        read_size = int(self.args["--read-size"])

        # Defined boolean variables
        normalized = self.args["--normalized"]
        self_contacts = self.args["--self-contacts"]

        mtn.alignment_to_contacts(
            bed2D_file=self.args["--input"],
            genome=self.args["--genome"],
            output_dir=self.args["--outdir"],
            output_file_network=self.args["--output-file-network"],
            output_file_contig_data=self.args["--output-file-contig-data"],
            tmpdir=temp_directory,
            read_size=read_size,
            n_cpus=self.args["--threads"],
            normalized=normalized,
            self_contacts=self_contacts,
        )

        # Delete the temporary folder
        shutil.rmtree(temp_directory)


# TODO: Check if the Louvain algorithm is available in Python.
class Partition(AbstractCommand):
    """"""


class Binning(AbstractCommand):
    """"""


class Validation(AbstractCommand):
    """"""


class Pipeline(AbstractCommand):
    """"""
