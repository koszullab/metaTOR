#!/usr/bin/env python3
# coding: utf-8

"""Abstract command classes for metaTOR

This module contains all classes related to metaTOR commands:

    - align
    - cutsite
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
import sys
import metator.align as mta
import metator.cutsite as mtc
import metator.io as mio
import metator.network as mtn
import metator.partition as mtp
import metator.validation as mtv
from docopt import docopt
from metator.log import logger
from os.path import exists, dirname, join

from pyinstrument import Profiler
from pyinstrument.renderers import ConsoleRenderer


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
    given separated by commas. Option to digest the reads using ligation
    sites to optimize the number of mapped reads.

    usage:
        align --forward=STR --reverse=STR --genome=FILE [--tempdir=DIR]
        [--threads=1] [--min-quality=30] [--no-clean-up] [--outdir=DIR]

    options:
        -1, --forward=STR       Fastq file or list of Fastq separated by a comma
                                containing the forward reads to be aligned.
        -2, --reverse=STR       Fastq file or list of Fastq separated by a comma
                                containing the reverse reads to be aligned.
                                Forward and reverse reads need to have the same
                                identifier.
        -g, --genome=FILE       The genome on which to map the reads. Must be
                                the path to the bowtie2/bwa index. Mandatory if
                                not digestion only enabled.
        -N, --no-clean-up       Do not remove temporary files.
        -o, --outdir=DIR        Path of the directory where the alignment will
                                be written in bed2D format and the digested
                                fastq. Default to current directory.
                                [Default: .]
        -q, --min-quality=INT   Threshold of quality necessary to considered a
                                read properly aligned. [Default: 30]
        -t, --threads=INT       Number of parallel threads allocated for the
                                alignement. [Default: 1]
        -T, --tempdir=DIR       Temporary directory. [Default: ./tmp]
    """

    def execute(self):

        # Start Profiler
        profiler = Profiler()
        profiler.start()

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        if not exists(self.args["--outdir"]):
            os.makedirs(self.args["--outdir"])

        # Transform integer variables as integer.
        min_qual = int(self.args["--min-quality"])

        # Align pair-end reads with bowtie2
        mta.pairs_alignment(
            self.args["--forward"],
            self.args["--reverse"],
            min_qual,
            temp_directory,
            self.args["--genome"],
            self.args["--outdir"],
            self.args["--threads"],
        )

        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)

        session = profiler.stop()
        profile_renderer = ConsoleRenderer(
            unicode=True, color=True, show_all=True
        )
        print(profile_renderer.render(session))


class Cutsite(AbstractCommand):
    """Cutsite command

    Generates new gzipped fastq files from original fastq. The function will cut
    the reads at their ligation sites and creates new pairs of reads with the
    different fragments obtained after cutting at the digestion sites.

    There are three choices to how combine the fragments. 1. "for_vs_rev": All
    the combinations are made between one forward fragment and one reverse
    fragment. 2. "all": All 2-combinations are made. 3. "pile": Only
    combinations between adjacent fragments in the initial reads are made.

    usage:
        cutsite --forward=STR --reverse=STR --enzyme=STR [--threads=1]
        [--outdir=DIR] [--mode=for_vs_rev]

    options:
        -1, --forward=STR       Fastq file or list of Fastq separated by a comma
                                containing the forward reads to be aligned.
        -2, --reverse=STR       Fastq file or list of Fastq separated by a comma
                                containing the reverse reads to be aligned.
                                Forward and reverse reads need to have the same
                                identifier.
        -e, --enzyme=STR        The list of restriction enzyme used to digest
                                the genome separated by a comma. Example:
                                DpnII,HinfI. If no option given, it will align
                                only once the reads. [Default: None]
        -m, --mode=STR          Digestion mode. There are three possibilities:
                                "for_vs_rev", "all" and "pile".
                                The first one "for_vs_rev" makes all possible
                                contact between fragments from forward read
                                versus the fragments of the reverse reads.
                                The second one "all" consist two make all pairs
                                of fragments possible.
                                The third one "pile" will make the contacts only
                                with the adjacent fragments.
                                [Default: for_vs_rev]
        -o, --outdir=DIR        Path of the directory where the alignment will
                                be written in bed2D format and the digested
                                fastq. Default to current directory.
                                [Default: .]
        -t, --threads=INT       Number of parallel threads allocated for the
                                alignement. [Default: 1]
    """

    def execute(self):

        # Start Profiler
        profiler = Profiler()
        profiler.start()

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        if not exists(self.args["--outdir"]):
            os.makedirs(self.args["--outdir"])

        # Digestion of the reads.
        logger.info("Digestion of the reads:")
        logger.info("Enzyme used: {0}".format(self.args["--enzyme"]))
        logger.info(
            "Mode used to cut the reads: {0}".format(self.args["--mode"])
        )

        mtc.cut_ligation_sites(
            self.args["--forward"],
            self.args["--reverse"],
            self.args["--enzyme"],
            self.args["--mode"],
            self.args["--outdir"],
            int(self.args["--threads"]),
        )

        session = profiler.stop()
        profile_renderer = ConsoleRenderer(
            unicode=True, color=True, show_all=True
        )
        print(profile_renderer.render(session))


class Network(AbstractCommand):
    """Generation of network command

    Generates a network file (in edgelist form) from an alignment in bed2D
    format. Contigs are the network nodes and the edges are the contact counts.

    The network is in a strict barebone form so that it can be reused and
    imported quickly into other applications etc. Verbose information about
    every single node in the network is written on a 'contig data' file, by
    default called 'idx_contig_length_GC_hit_cov.txt'

    usage:
        network --genome=FILE --outdir=DIR --input=FILE --depth=FILE
        [--output-file-contig-data=STR] [--self-contacts] [--no-clean-up]
        [--output-file-network=STR] [--tempdir=DIR] [--normalization=STR]
        [--threads=1] [--enzyme=STR]

    options:
        -d, --depth=FILE                The depth.txt file from
                                        jgi_summarize_bam_contig_depths from
                                        metabat2 pipeline.
        -e, --enzyme=STR                The list of restriction enzyme used to
                                        digest the contigs separated by a comma.
                                        Example: DpnII,HinfI.
        -g, --genome=FILE               The initial assembly path acting as the
                                        alignment file's reference genome.
        -i, --input=FILE                Path to the bed2D file used as input
        -n, --normalization=STR             If None, do not normalized the count of
                                        a contact by the geometric mean of the
                                        coverage of the contigs. Otherwise it's
                                        the type of normalization. 7 values are
                                        possible None, abundance, length, RS,
                                        RS_length, empirical_hit,
                                        theoritical_hit. [Default: abundance]
        -N, --no-clean-up               Do not remove temporary files.
        -o, --outdir=DIR                The output directory to write the
                                        network and contig data into. Default:
                                        current directory.
        --output-file-contig-data=STR   The specific file name for the output
                                        chunk data file. [Default:
                                        contig_data_network.txt]
        --output-file-network=STR       The specific file name for the output
                                        network file. [Default: network.txt]
        -s, --self-contacts             If enabled, count alignments between a
                                        contig and itself.
        -t, --threads=INT               Number of parallel threads allocated for
                                        the alignement. [Default: 1]
        -T, --tempdir=DIR               Temporary directory. Default to current
                                        directory. [Default: ./tmp]
    """

    def execute(self):

        # Start Profiler
        profiler = Profiler()
        profiler.start()

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        if not exists(self.args["--outdir"]):
            os.makedirs(self.args["--outdir"])

        if not self.args["--output-file-contig-data"]:
            self.args["--output-file-contig-data"] = "contig_data_network.txt"

        if not self.args["--output-file-network"]:
            self.args["--output-file-network"] = "network.txt"

        # Defined boolean variables
        self_contacts = self.args["--self-contacts"]

        # Transform str input files in list splitting on comma
        alignment_files = self.args["--input"].split(",")

        mtn.alignment_to_contacts(
            alignment_files,
            self.args["--genome"],
            self.args["--depth"],
            self.args["--outdir"],
            self.args["--output-file-network"],
            self.args["--output-file-contig-data"],
            temp_directory,
            self.args["--threads"],
            self.args["--normalization"],
            self.args["--enzyme"],
            self_contacts,
        )

        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)

        session = profiler.stop()
        profile_renderer = ConsoleRenderer(
            unicode=True, color=True, show_all=True
        )
        print(profile_renderer.render(session))


class Partition(AbstractCommand):
    """Partition the network using Louvain algorithm

    Partition the network file using iteratively the Louvain or Leiden
    algorithm. Then looks for 'cores' that are easily found by identifying
    identical lines on the global output. Using hamming distance from these core
    bins, group the bins with more than the percentage (overlapping parameter)
    given.

    It will also update the file to integrate the bins information of the
    contigs. If the version of Louvain is not found, the python version of
    Louvain will be used.

    Furthermore, both Leiden and Louvain algorithm are available here. However,
    the benchmark made show taht with this pipeline the Louvain algorithm gives
    better results and is faster on seawater and gut metagenomic samples.

    Note that the Louvain software is not, in the strictest sense, necessary.
    Any program that assigns a node to a bin, does so non deterministically and
    solely outputs a list in the form: 'node_id bin_id' could be plugged
    instead.

    usage:
        partition  --outdir=DIR --network-file=FILE --assembly=FILE --contigs-data=FILE
        [--iterations=100] [--algorithm=louvain] [--overlap=80] [--size=500000]
        [--threads=1] [--tempdir=DIR] [--no-clean-up] [--res-parameter=1.0]

    options:
        -a, --assembly=FILE         The path to the assembly fasta file used to
                                    do the alignment.
        -A, --algorithm=STR         louvain|leiden, algorithm to use to
                                    partition the network. [Default: louvain]
        -c, --contigs-data=FILE     The path to the tsv file containing the data
                                    of the contigs (ID, Name, Length, GC
                                    content, Hit, Coverage).
        -i, --iterations=INT        Number of iterations of Louvain.
                                    [Default: 100]
        -n, --network-file=FILE     Path to the file containing the network
                                    information from the meta HiC experiment
                                    compute in network function previously.
        -N, --no-clean-up           Do not remove temporary files.
        -o, --outdir=DIR            Path to the directory to write the output.
                                    Default to current directory. [Default: ./]
        -O, --overlap=INT           Percentage of the identity necessary to be
                                    considered as a part of the core bin.
                                    [Default: 80]
        -r, --res-parameter=FLOAT   Resolution paramter to use for Leiden
                                    algorithm. [Default: 1.0]
        -s, --size=INT              Threshold size to keep bins in base pair.
                                    [Default: 500000]
        -t, --threads=INT           Number of parallel threads allocated for the
                                    partition. [Default: 1]
        -T, --tempdir=DIR           Temporary directory. Default to current
                                    directory. [Default: ./tmp]
    """

    def execute(self):

        # Start Profiler
        profiler = Profiler()
        profiler.start()

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
            overlapping_parameter = int(self.args["--overlap"]) / 100
        if self.args["--size"]:
            size = int(self.args["--size"])
        if self.args["--threads"]:
            threads = int(self.args["--threads"])
        if self.args["--res-parameter"]:
            resolution_parameter = float(self.args["--res-parameter"])

        # Partition the network
        mtp.partition(
            self.args["--algorithm"],
            self.args["--assembly"],
            self.args["--contigs-data"],
            iterations,
            self.args["--network-file"],
            self.args["--outdir"],
            overlapping_parameter,
            resolution_parameter,
            size,
            temp_directory,
            threads,
        )

        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)

        session = profiler.stop()
        profile_renderer = ConsoleRenderer(
            unicode=True, color=True, show_all=True
        )
        print(profile_renderer.render(session))


class Validation(AbstractCommand):
    """Use CheckM to validate the bins.

    Use checkM to validate bacterial and archae bins. The script returns the
    output of CheckM is an output directory.

    It is possible to also partition again the contaminated bins to improve
    them. The new bins contamination and completion will be compute again. If
    there is a loss of the completion from the original bins, i.e. the new
    iterations may split the organism in multiple bins, go back to the original
    bins.

    usage:
        validation --outdir=DIR --network=FILE --assembly=FILE --fasta=DIR --contigs=STR
        [--iterations=10] [--algorithm=louvain] [--size=500000] [--res-param=1.0]
        [--threads=1] [--tempdir=DIR]  [--no-clean-up] [--overlap=90]

    options:
        -a, --assembly=FILE     The path to the assembly fasta file used to do
                                the alignment.
        -A, --algorithm=STR     Algorithm to use. Either "louvain" or "leiden".
                                [Default: louvain]
        -c, --contigs=FILE      The path to the file containing the data ofthe
                                contigs (ID, Name, Length, GC content, Hit,
                                Coverage).
        -f, --fasta=DIR         Path to the directory containing the input fasta
                                files of the bins.
        -i, --iterations=INT    Number of recursive iterations of Louvain.
                                [Default: 10]
        -n, --network=FILE      Path to the file containing the network
                                information from the meta HiC experiment compute
                                in network function previously.
        -N, --no-clean-up       Do not remove temporary files.
        -o, --outdir=DIR        Path to the directory to write the output.
                                Default to current directory. [Default: ./]
        -O, --overlap=INT       Percentage of the identity necessary to be
                                considered as a part of the core bin.
                                [Default: 90]
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

        # Start Profiler
        profiler = Profiler()
        profiler.start()

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        if not exists(self.args["--outdir"]):
            os.makedirs(self.args["--outdir"])
        if not exists(self.args["--outdir"]):
            os.makedirs(os.join(self.args["--outdir"], "fasta"))

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
            sys.exit(1)

        # Check correct algorithm value
        if self.args["--algorithm"] not in ["louvain", "leiden"]:
            logger.error('algorithm should be either "louvain" or "leiden"')
            raise ValueError

        mtv.recursive_decontamination(
            self.args["--algorithm"],
            self.args["--assembly"],
            self.args["--contigs"],
            self.args["--fasta"],
            iterations,
            self.args["--network"],
            self.args["--outdir"],
            overlapping_parameter,
            resolution_parameter,
            size,
            temp_directory,
            threads,
        )

        session = profiler.stop()
        profile_renderer = ConsoleRenderer(
            unicode=True, color=True, show_all=True
        )
        print(profile_renderer.render(session))


class Pipeline(AbstractCommand):
    """Launch the full metator pipeline

    Partition the assembly in bins from the HiC reads of the metapopulation.

    It's possible to start from the fastq, the bam, the bed2D, or the network
    files. It's also possible to ask or not to run the validation step which is
    the critical step for memory usage.

    usage: pipeline  --fasta=FILE --forward
        reads_for.fastq[,reads_for2.fastq...] --reverse
        reads_rev.fastq[,reads_rev2.fastq...] [--assembly=FILE] [--tempdir=DIR]
        [--threads=1] [--normalized] [--no-clean-up] [--overlap=90]
        [--iterations=100] [--size=100] [--self-contacts] [--min-quality=30]
        [--algorithm=STR] [--outdir=DIR]

    options:
        -1, --forward=STR           Fastq file or list of Fastq separated by a
                                    comma containing the forward reads to be
                                    aligned.
        -2, --reverse=STR           Fastq file or list of Fastq separated by a
                                    comma containing the reverse reads to be
                                    aligned. Forward and reverse reads need to
                                    have the same identifier.
        -A, --algorithm=STR         Path to louvain cpp (faster than python
                                    implementation). If None given, use python
                                    implementation instead. [Default: None]
        -f, --fasta=FILE            The genome on which to map the reads. Must
                                    be the path to the bowtie2/bwa index or the
                                    fasta.
        -i, --iterations=INT        Number of iterartion of Louvain.
                                    [Default: 100]
        -n, --normalized            If enabled,  normalize contacts between
                                    contigs by their geometric mean coverage.
        -N, --no-clean-up           Do not remove temporary files.
        -o, --outdir=DIR            Path where the alignment will be written in
                                    bed2D format.
        -O, --overlap=INT           Percentage of the identity necessary to be
                                    considered as a part of the core bin.
                                    [Default: 90]
        -q, --min-quality=INT       Threshold of quality necessary to considered
                                    a read properly aligned. [Default: 30]
        -s, --size=INT              Threshold size to keep bins in base pair.
                                    [Default: 300000]
        -S, --self-contacts         If enabled, count alignments between a
                                    contig and itself.
        -t, --threads=INT           Number of parallel threads allocated for the
                                    alignement. [Default: 1]
        -T, --tempdir=DIR           Temporary directory. [Default: ./tmp]
    """

    def execute(self):

        # Start Profiler
        profiler = Profiler()
        profiler.start()

        # Defined the temporary directory.
        if not self.args["--tempdir"]:
            self.args["--tempdir"] = "./tmp"
        temp_directory = mio.generate_temp_dir(self.args["--tempdir"])

        # Defined the output directory and output file names.
        if not self.args["--outdir"]:
            self.args["--outdir"] = "."
        if not exists(self.args["--outdir"]):
            os.makedirs(self.args["--outdir"])

        # Transform integer variables as integer.
        if self.args["--min-quality"]:
            min_qual = int(self.args["--min-quality"])
        if self.args["--iterations"]:
            iterations = int(self.args["--iterations"])
        if self.args["--overlap"]:
            overlap = float(self.args["--overlap"])
        if self.args["--size"]:
            size = int(self.args["--size"])
        if self.args["--threads"]:
            threads = int(self.args["--threads"])

        # Defined boolean variables.
        normalized = self.args["--normalized"]
        self_contacts = self.args["--self-contacts"]

        # Create two path for the fasta index or the fasta assembly from the
        # given file.
        index = mio.check_fasta_index(self.args["--fasta"])
        if index is None:
            if mio.check_is_fasta(self.args["--fasta"]):
                fasta = self.args["--fasta"]
                index = self.args["--fasta"]
        else:
            fasta = mio.retrieve_fasta(temp_directory, temp_directory)

        # Align pair-end reads with bowtie2.
        pairs = mta.pairs_alignment(
            self.args["--forward"],
            self.args["--reverse"],
            min_qual,
            temp_directory,
            index,
            self.args["--outdir"],
            self.args["--threads"],
        )

        # Generate the network.
        mtn.alignment_to_contacts(
            pairs,
            fasta,
            self.args["--outdir"],
            "network.txt",
            "contigs_data_network.txt",
            temp_directory,
            self.args["--threads"],
            normalized,
            self_contacts,
        )

        network_file = join(self.args["--outdir"], "network.txt")
        contigs_data = join(self.args["--outdir"], "contigs_data_network.txt")

        # Perform iterations of Louvain.
        if self.args["--algorithm"] == "None":
            output_partition = mtp.louvain_iterations_py(
                network_file,
                iterations,
            )
        else:
            output_partition = mtp.louvain_iterations_cpp(
                network_file,
                iterations,
                temp_directory,
                self.args["--algorithm"],
            )

        # Detect core bins
        (
            core_bins,
            core_bins_iterations,
        ) = mtp.detect_core_bins(output_partition, iterations)

        # Compute the Hamming distance between core bins.
        hamming_distance = mtp.hamming_distance(
            core_bins_iterations,
            iterations,
            threads,
        )

        # Defined overlapping bins according to the threshold
        overlapping_bins = mtp.defined_overlapping_bins(
            overlap,
            hamming_distance,
            core_bins,
            core_bins_iterations,
        )

        # Update the contigs_data_file.
        contigs_data = mtp.update_contigs_data(
            contigs_data,
            core_bins,
            overlapping_bins,
            self.args["--outdir"],
        )

        # Generate Fasta file
        mtp.generate_fasta(
            fasta,
            overlapping_bins,
            contigs_data,
            size,
            self.args["--outdir"],
            temp_directory,
        )

        # TODO: Launch validation if necessary.

        # Delete the temporary folder.
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)

        session = profiler.stop()
        profile_renderer = ConsoleRenderer(
            unicode=True, color=True, show_all=True
        )
        print(profile_renderer.render(session))
