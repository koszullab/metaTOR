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
import metator.align as mta
import metator.cutsite as mtc
import metator.io as mio
import metator.network as mtn
import metator.partition as mtp
from docopt import docopt
from metator.log import logger
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
        [--self-contacts] [--tempdir=DIR] [--threads=1] [--no-clean-up]

    options:
        -g, --genome=FILE               The initial assembly path acting as the
                                        alignment file's reference genome.
        -i, --input=FILE                Path to the bed2D file used as input
        -n, --normalized                If enabled,  normalize contacts between
                                        contigs by their geometric mean
                                        coverage.
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
        normalized = self.args["--normalized"]
        self_contacts = self.args["--self-contacts"]

        mtn.alignment_to_contacts(
            bed2D_file=self.args["--input"],
            genome=self.args["--genome"],
            output_dir=self.args["--outdir"],
            output_file_network=self.args["--output-file-network"],
            output_file_contig_data=self.args["--output-file-contig-data"],
            tmpdir=temp_directory,
            n_cpus=self.args["--threads"],
            normalized=normalized,
            self_contacts=self_contacts,
        )

        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)


# TODO: Check if the Louvain algorithm is available in Python. Else keep
# something as the tools of Lyam to execute Louvain.
class Partition(AbstractCommand):
    """Partition the network using Louvain algorithm

    Partition the network file using iteratively the Louvain algorithm. Then
    looks for 'cores' that are easily found by identifying identical lines on
    the global Louvain output. Using hamming distance from these core bins,
    group the bins with more than the percentage given as the overlap value.

    If the contigs data information is given, it will also update the file to
    integrate the bins information of the contigs. If the version of Louvain is
    not found, the python version of Louvain will be used.

    Note that the Louvain software is not, in the strictest sense, necessary.
    Any program that assigns a node to a bin, does so non deterministically and
    solely outputs a list in the form: 'node_id bin_id' could be plugged
    instead.

    usage:
        partition  --outdir=DIR --network-file=FILE --assembly=FILE
        [--iterations=100] [--louvain=STR] [--overlap=INT] [--size=300000]
        [--threads=1] [--tempdir=DIR] [--contigs-data=STR] [--no-clean-up]

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
        -N, --no-clean-up           Do not remove temporary files.
        -o, --outdir=DIR            Path to the directory to write the output.
                                    Default to current directory. [Default: ./]
        -O, --overlap=INT           Percentage of the identity necessary to be
                                    considered as a part of the core bin.
                                    [Default: 90]
        -s, --size=INT              Threshold size to keep bins in base pair.
                                    [Default: 300000]
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
            overlap = int(self.args["--overlap"]) / 100
        if self.args["--size"]:
            size = int(self.args["--size"])
        if self.args["--threads"]:
            threads = int(self.args["--threads"])

        # TODO: Test which function is necessary --> C or Python and call it

        # Perform the iterations of Louvain to partition the network.
        logger.info("Start iterations of Louvain:")
        if self.args["--louvain"] == "cpp":
            louvain_iterations_cpp()
        else:
            output_louvain = mtp.louvain_iterations_py(
                self.args["--network-file"],
                iterations,
            )
        # Detect core bins
        logger.info("Detect core bins:")
        (
            core_bins,
            core_bins_iterations,
        ) = mtp.detect_core_bins(output_louvain, iterations)

        # Compute the Hamming distance between core bins.
        logger.info("Detect overlapping bins:")
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
        logger.info("Extract bins:")
        contigs_data = mtp.update_contigs_data(
            self.args["--contigs-data"],
            core_bins,
            overlapping_bins,
        )

        # Generate Fasta file
        mtp.generate_fasta(
            self.args["--assembly"],
            overlapping_bins,
            contigs_data,
            size,
            self.args["--outdir"],
        )

        # Delete the temporary folder
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)


class Validation(AbstractCommand):
    """Use CheckM to validate the bins.

    Use checkM to validate bacterial and archae bins. The script returns the
    output of CheckM is an output directory.

    It is possible to also partition again the contaminated bins to improve
    them. The new bins contamination and completion will be compute again. If
    there is a loss of the completion from the original bins, i.e. the new
    iterations may split the organism in multiple bins, go back to the original
    bins.

    usage: validation

    options:

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


class Pipeline(AbstractCommand):
    """Launch the full metator pipeline

    Partition the assembly in bins from the HiC reads of the metapopulation.

    It's possible to start from the fastq, the bam, the bed2D, or the network
    files. It's also possible to ask or not to run the validation step which is
    the critical step for memory usage.

    usage: pipeline [--tempdir=DIR] [--threads=1] [--normalized] [--no-clean-up]
        [--overlap=90] [--iterations=100] [--size=100] [--self-contacts]
        [--min-quality=30] --genome=FILE --out=DIR --forward
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
        -i, --iterations=INT        Number of iterartion of Louvain.
                                    [Default: 100]
        -n, --normalized            If enabled,  normalize contacts between
                                    contigs by their geometric mean coverage.
        -N, --no-clean-up           Do not remove temporary files.
        -o, --out=DIR               Path where the alignment will be written in
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

        # Create two path for the fasta index or the fasta assembly.
        assembly, assembly_index = mio.check_fasta_index(
            self.args["--assembly"]
        )

        # Align pair-end reads with bowtie2.
        pairs = mta.pairs_alignment(
            self.args["--forward"],
            self.args["--reverse"],
            min_qual,
            temp_directory,
            self.args["--genome"],
            None,
            self.args["--out"],
            self.args["--threads"],
        )

        # Generate the network.
        network, contigs_data = mtn.alignment_to_contacts(
            bed2D_file=pairs,
            genome=self.args["--genome"],
            output_dir=self.args["--outdir"],
            output_file_network="network.txt",
            output_file_contig_data="idx_contig_length_GC_hit_cov.txt",
            tmpdir=temp_directory,
            n_cpus=self.args["--threads"],
            normalized=normalized,
            self_contacts=self_contacts,
        )

        # Perform iterations of Louvain.
        output_louvain = mtp.louvain_iterations_py(
            network,
            iterations,
        )

        # Detect core bins
        (
            core_bins,
            core_bins_iterations,
        ) = mtp.detect_core_bins(output_louvain, iterations)

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
        )

        # Generate Fasta file
        mtp.generate_fasta(
            assembly,
            overlapping_bins,
            contigs_data,
            size,
            self.args["--outdir"],
        )

        # TODO: Launch validation if necessary.

        # Delete the temporary folder.
        if not self.args["--no-clean-up"]:
            shutil.rmtree(temp_directory)