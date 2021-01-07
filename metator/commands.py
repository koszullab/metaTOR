#!/usr/bin/env python3
# coding: utf-8

"""Abstract command classes for metaTOR

This module contains all classes related to metaTOR commands:

    -align
    -network 
    -partition
    -pipeline
    -validation

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
import metator.align as mta
import metator.io as mio
import metator.network as mtn
import metator.partition as mtp
from docopt import docopt
from os.path import exists


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

        # Defined the output directory and output file names.
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


# TODO: Check if the Louvain algorithm is available in Python. Else keep
# something as the tools of Lyam to execute Louvain.
class Partition(AbstractCommand):
    """Partition the network using Louvain algorithm

    Partition the network file using iteratively the Louvain algorithm. Then
    looks for 'cores' that are easily found by identifying identical lines on
    the global Louvain output. Using hamming distance from these core
    communities, group the communities with more than the percentage given as
    the overlap value.

    If the contigs data information is given, it will also update the file to
    integrate the communities information of the contigs.
    If the version of Louvain is not found, the python version of Louvain will
    be used.

    Note that the Louvain software is not, in the strictest sense, necessary.
    Any program that assigns a node to a community, does so non
    deterministically and solely outputs a list in the form: 'node_id
    community_id' could be plugged instead.

    usage:
        partition  --outdir=DIR --network-file=FILE --assembly=FILE
        [--iterations=100] [--louvain=STR] [--overlap=90] [--size=300000]
        [--threads=1] [--tempdir=DIR] [--contigs-data=STR]

    options:
        -a, --assembly=FILE         The path to the assembly fasta file used to
                                    do the alignment.
        -c, --contigs-data=FILE     The path to the file containing the data of
                                    the contigs (ID, Name, Length, GC content,
                                    Hit, Coverage).
        -i, --iterations=INT        Number of iterartion of Louvain.
                                    [Default: 100]
        -l, --louvain=STR           Informatic language used for louvain
                                    algorithm. Two values are possible: "py" or
                                    "cpp". [Default: py]
        -n, --network-file=FILE     Path to the file containing the network
                                    information from the meta HiC experiment
                                    compute in network function previously.
        -o, --outdir=DIR            Path to the directory to write the output.
                                    Default to current directory. [Default: ./]
        -O, --overlap=INT           Percentage of the identity necessary to be
                                    considered as a part of the core community.
                                    [Default: 90]
        -s, --size=INT              Threshold size to keep communities in base
                                    pair. [Default: 300000]
        -t, --threads=INT           Number of parallel threads allocated for the
                                    partition. [Default: 1]
        -T, --tempdir=DIR           Temporary directory. Default to current
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
        if not exists(self.args["--outdir"]):
            os.makedirs(self.args["--outdir"])

        # Transform numeric variable as numeric
        if self.args["--iterations"]:
            iterations = int(self.args["--iterations"])
        if self.args["--overlap"]:
            overlap = float(self.args["--overlap"])
        if self.args["--size"]:
            size = int(self.args["--size"])
        if self.args["--threads"]:
            threads = int(self.args["--threads"])

        # TODO: Test which function is necessary --> C or Python and call it

        # Perform the iterations of Louvain to partition the network.
        if self.args["--louvain"] == "cpp":
            louvain_iterations_cpp()
        else:
            output_louvain = mtp.louvain_iterations_py(
                self.args["--network-file"],
                iterations,
            )
        # Detect core communities
        (
            core_communities,
            core_communities_iterations,
        ) = mtp.detect_core_communities(output_louvain, iterations)

        # Compute the Hamming distance between core communities.
        hamming_distance = mtp.hamming_distance(
            core_communities_iterations,
            iterations,
            threads,
        )

        # Defined overlapping communities according to the threshold
        overlapping_communities = mtp.defined_overlapping_communities(
            overlap,
            hamming_distance,
            core_communities,
            core_communities_iterations,
        )

        # Update the contigs_data_file.
        contigs_data = mtp.update_contigs_data(
            self.args["--contigs-data"],
            core_communities,
            overlapping_communities,
        )

        # Generate Fasta file
        mtp.generate_fasta(
            self.args["--assembly"],
            overlapping_communities,
            contigs_data,
            size,
            self.args["--outdir"],
        )


class Validation(AbstractCommand):
    """Use CheckM to validate the communities.

    Use checkM to validate bacterial and archae communities. The script returns
    the output of CheckM is an output directory.

    It is possible to also partition again the contaminated communities to
    improve them. The new communities contamination and completion will be
    compute again. If there is a loss of the completion from the original
    communities, i.e. the new iterations may split the organism in multiple
    communities, go back to the original communities.

    usage: validation

    options:

    """

    # Launch checkM to evaluate the completion and the contamination. If asked
    # rerun Louvain to try to reduce the contamination, rerun checkM if the
    # contamination decrease without a huge decrease of the completion keep the
    # new communities. Otherwise go back to the old state.
    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."


class Pipeline(AbstractCommand):
    """Launch the full metator pipeline

    Partition the assembly in communities from the HiC reads of the
    metapopulation.

    It's possible to start from the fastq, the bam, the bed2D, or the network
    files. It's also possible to ask or not to run the validation step which is
    the critical step for memory usage.

    usage: pipeline

    options:

    """

    def execute(self):

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        # Launch alignment
        # Launch network
        # Launch binning
        # Launch validation if necessary
