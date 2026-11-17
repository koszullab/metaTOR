#!/usr/bin/env python3
# coding: utf-8

"""metaTOR command line interface.

One function per command: network, partition, validation, pipeline, qc,
contactmap, metacool, scaffold, pairs and mge. The ``metator`` entry point
is :func:`cli`.
"""

import logging
import multiprocessing as mp
import os
import re
import shlex
import shutil
import sys
import time
from contextlib import contextmanager
from functools import partial
from os.path import dirname, exists, join

import rich_click as click
from scipy.sparse import load_npz, save_npz

import metator.align as mta
import metator.contact_map as mtc
import metator.io as mio
import metator.log as mtl
import metator.metacool as mtmc
import metator.mge as mtm
import metator.network as mtn
import metator.partition as mtp
import metator.quality_check as mtq
import metator.scaffold as mts
import metator.validation as mtv
from metator.log import logger
from metator.version import __version__

ALIGNERS = ["bowtie2", "bwa"]
ALIGNER_MODES = ["normal", "iterative", "cutsite"]
ALGORITHMS = ["louvain", "leiden"]
NORMALIZATIONS = ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"]
STARTS = ["fastq", "bam", "pair", "metacool"]
MATRIX_FORMATS = ["graal", "bg2", "cool"]

# Panel holding the options common to most commands (output, temporary
# files, threads), with --help.
GENERAL = "Options"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def check_output_path(path, force=False):
    """Throws error if the output file exists. Create required file tree
    otherwise.
    """
    # Get complete output filename and prevent overwriting unless force is
    # enabled
    if not force and os.path.exists(path):
        raise IOError("Output file already exists. Use --force to overwrite")
    if dirname(path):
        os.makedirs(dirname(path), exist_ok=True)


def bin_prefix(user_prefix):
    """Prefix of the bin names: "metator" or "<--prefix>_metator".

    Bin names become fasta file names that miComplete renames (it drops
    every dot), so only letters, digits, "_" and "-" are accepted.
    """
    if not user_prefix:
        return "metator"
    if not re.fullmatch(r"[A-Za-z0-9_-]+", user_prefix):
        logger.error(
            "--prefix may only contain letters, digits, '_' and '-': %r",
            user_prefix,
        )
        raise ValueError(f"Invalid --prefix: {user_prefix!r}")
    return f"{user_prefix}_metator"


@contextmanager
def private_tmp_dir(tmpdir, clean_up=True):
    """This run's private temporary directory, inside --tmpdir.

    --tmpdir may be a directory the user already keeps files in (for
    example /scratch/$USER), so the run works in a fresh random
    subdirectory, and the clean-up at the end of the ``with`` block
    removes only that. --tmpdir itself is removed only if this run created
    it and nothing else is left in it.

    Nothing is removed without clean_up (--no-clean-up), nor when the run
    fails: its files are left for inspection.
    """
    base = tmpdir or "./tmp"
    base_created = not exists(base)
    tmp_dir = mio.generate_temp_dir(base)
    logger.debug("Temporary directory: %s", tmp_dir)
    yield tmp_dir
    if clean_up:
        shutil.rmtree(tmp_dir)
        if base_created and not os.listdir(base):
            os.rmdir(base)


def refuse_output_dir(path, force):
    """Refuse an existing output directory of bins, unless --force replaces
    it."""
    if exists(path) and not force:
        logger.error(
            "%s already existed. Remove directory or use -F argument to overwrite it.",
            path,
        )
        raise ValueError


def make_output_dir(path, force):
    """Create an output directory of bins; an existing one is replaced with
    --force, or refused.
    """
    refuse_output_dir(path, force)
    if exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def check_normalization(normalization, enzyme, depth):
    """The normalizations that need the restriction sites or the shotgun
    depth of the contigs."""
    if normalization in ["RS", "theoritical_hit"] and not enzyme:
        logger.error('For "RS" and "theoritical_hit" normalization, enzyme is required.')
        raise ValueError
    if normalization in ["abundance", "theoritical_hit"] and not depth:
        logger.error('For "abundance" and "theoritical_hit" normalization, depth is required.')
        raise ValueError


def check_reverse_reads(start, aligner, reverse):
    """Reverse reads (or bam files) are needed, except from pairs or from a
    bwa alignment."""
    if (start == "fastq" or (start == "bam" and aligner == "bowtie2")) and not reverse:
        logger.error(
            "Forward and reverse arguments are necessary for fastq with %s start and %s aligner.",
            start,
            aligner,
        )
        raise ValueError


def check_metacool_start(start, forward, edge, scaffold=False):
    """A .metacool start: one contact map, without the edge region of the
    contigs, nor the scaffolding (pipeline), which read the pairs files."""
    if start != "metacool":
        return
    if not forward:
        logger.error("--start metacool needs the .metacool file, given with -1/--forward.")
        raise ValueError
    if "," in forward:
        logger.error(
            "--start metacool takes one .metacool file (its libraries are already summed), given with "
            "-1/--forward: got %d files.",
            len(forward.split(",")),
        )
        raise ValueError
    if edge > 0:
        logger.error(
            "--edge %d cannot be applied with --start metacool: a contact map holds the contacts of bins, not the "
            "positions of the reads. Leave --edge at 0, or start from the pairs files.",
            edge,
        )
        raise ValueError
    if scaffold:
        logger.error("--scaffold reads the pairs files: it cannot run with --start metacool.")
        raise ValueError


def contig_network(start, forward, alignment_files, contig_data, hit_data, edge, outdir, tmp_dir, threads, *args):
    """network.txt and contig_data_network.txt of a run (args: normalization
    and self contacts): from its pairs files, or, with --start metacool,
    from the .metacool map given to --forward (see
    metator.network.metacool_to_contacts). A map that cannot give them
    stops the run, as metaTOR stops on an input it cannot use: the reason
    is logged as an error, in the log file of the run, and its
    MetacoolError is raised, as `metator metacool` does."""
    files = ("network.txt", "contig_data_network.txt")
    if start != "metacool":
        return mtn.alignment_to_contacts(alignment_files, contig_data, edge, hit_data, outdir, *files, tmp_dir, threads, *args)
    try:
        return mtn.metacool_to_contacts(forward, contig_data, outdir, *files, tmp_dir, threads, *args)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise


def aligner_mode_for(aligner_mode, enzyme):
    """The cutsite mode digests the reads: without an enzyme, align
    iteratively instead."""
    if aligner_mode == "cutsite" and not enzyme:
        logger.warning("cutsite mode required an enzyme. Iterative mode will be use instead.")
        return "iterative"
    return aligner_mode


def reference_files(assembly, aligner, start, tmp_dir):
    """Index and fasta of the assembly.

    The assembly is a fasta, whose index is built when the reads are
    aligned (fastq start), or the basename of a bowtie2 or bwa index, whose
    fasta is retrieved.
    """
    index = mio.check_fasta_index(assembly, mode=aligner)
    if index is None:
        if mio.check_is_fasta(assembly):
            fasta = assembly
            # Starting from bam or pairs files skips the index generation.
            if start == "fastq":
                index = mio.generate_fasta_index(fasta, aligner, tmp_dir)
        else:
            logger.error("Please give as assembly argument a %s index or a fasta.", aligner)
            raise ValueError
    else:
        fasta = mio.retrieve_fasta(index, aligner, tmp_dir)
    return index, fasta


def generate_log_header(log_path):
    """Start the log file with the version, the date and the command line."""
    mtl.set_file_handler(log_path, formatter=logging.Formatter(""))
    logger.info(f"## MetaTOR: v{__version__} log file")
    logger.info(f"## date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"## Command used: {shlex.join(['metator', *sys.argv[1:]])}")
    logger.info("\n---\n")
    mtl.set_file_handler(log_path, formatter=mtl.logfile_formatter)


def generate_log_footer(log_path):
    mtl.set_file_handler(log_path, formatter=logging.Formatter(""))
    logger.info("\n---")
    logger.info(f"## date: {time.strftime('%Y-%m-%d %H:%M:%S')}")


def log_file_name(directory, command):
    """Log file of a run: <directory>/metator_<command>_<date>.log."""
    now = time.strftime("%Y%m%d%H%M%S")
    return join(directory, f"metator_{command}_{now}.log")


def run_metadata():
    """The command line and the parameters of the running command, as stored
    in a .metacool file."""
    params = click.get_current_context().params
    return {
        "command": shlex.join(["metator", *sys.argv[1:]]),
        "parameters": {name: list(value) if isinstance(value, tuple) else value for name, value in params.items()},
    }


def assembly_basename(assembly):
    """Name of an assembly: its fasta file name without the extensions, or
    the basename of its index."""
    return re.sub(r"(\.(fa|fas|fasta|fna))?(\.gz)?$", "", os.path.basename(os.path.normpath(assembly)))


def pipeline_metacool(pairs_files, outdir, prefix, assembly, threads, metacool=None):
    """Write the contact map of a pipeline run: <outdir>/<prefix>.metacool,
    or <outdir>/metator.metacool without --prefix.

    From the pairs files of the run or, for a run that started from a
    .metacool map (metacool), from that map: its bins and contacts, grouped
    by final bin (see metator.metacool.reassign_metacool).

    The map is an extra output, written last: when it cannot be written, a
    warning says why, and the run keeps its other outputs and its exit
    status. A map of the same name left by an earlier run is removed first:
    this run has replaced the tables it was made from, so it must not stay
    behind, passing for this run's, when the new map cannot be written. An
    outfile that is the input map of the run (the same path, or a symbolic
    link to it or from it) is not: it is only replaced once the new map has
    been written and validated, and stays as it was when it cannot be. A
    symbolic link at the outfile is replaced, not its target.
    """
    outfile = join(outdir, f"{prefix or 'metator'}.metacool")
    # What the warning adds when the new map is not written.
    kept = ""
    try:
        if metacool and mtmc.is_input_map(metacool, outfile):
            kept = " It is still the input map of the run."
            logger.info(f"{outfile} is the input map {metacool}: it is replaced once the map of the bins is validated.")
        elif os.path.lexists(outfile):
            logger.info(f"Removing {outfile}, the contact map of an earlier run.")
            os.remove(outfile)
        contig_data = join(outdir, "contig_data_final.txt")
        shared = {"bin_summary": join(outdir, "bin_summary.txt"), "force": True, "metadata": run_metadata()}
        if metacool:
            mtmc.reassign_metacool(metacool, contig_data, outfile, **shared)
        else:
            mtmc.write_metacool(
                pairs_files, contig_data, outfile, assembly_name=assembly_basename(assembly), threads=threads, **shared
            )
    except Exception as error:
        # A MetacoolError says why; another error is unexpected: its type
        # and traceback go to the log too.
        unexpected = not isinstance(error, mtmc.MetacoolError)
        reason = f"{type(error).__name__}: {error}" if unexpected else str(error)
        logger.warning(
            f"{outfile} is not written: {reason.rstrip('.')}.{kept} The other outputs of the run are kept.",
            exc_info=unexpected,
        )


# ---------------------------------------------------------------------------
# Options shared by several commands
# ---------------------------------------------------------------------------


def _option(*names, panel=GENERAL, **attrs):
    """A command option, shown in the given help panel."""
    return click.option(*names, panel=panel, **attrs)


def _options(*decorators):
    """Several options at once, listed in the order of the help."""

    def decorate(f):
        for decorator in reversed(decorators):
            f = decorator(f)
        return f

    return decorate


def pairsfile_argument():
    return click.argument("pairsfile", nargs=-1, required=True, help="File(s) containing pairs information.", panel="Input")


ASSEMBLY_INDEX_HELP = (
    "The initial assembly path acting as the alignment file's reference genome or the basename of the bowtie2 or bwa index."
)
ASSEMBLY_PARTITION_HELP = "The path to the assembly fasta file used to do the alignment."


def assembly_option(help):
    return _option("-a", "--assembly", required=True, metavar="FILE", panel="Input", help=help)


def reads_options(forward_required):
    """Inputs of the network step: reads, alignments or pairs."""
    return _options(
        _option(
            "-1",
            "--forward",
            required=forward_required,
            metavar="STR",
            panel="Input",
            help="Fastq file or list of Fastq separated by a comma containing the forward reads to be aligned or their "
            "corresponding bam or pairs files. With --start metacool, one .metacool contact map (or an entry of a "
            ".mmetacool collection: file.mmetacool::entries/NAME) of the contigs of the assembly.",
        ),
        _option(
            "-2",
            "--reverse",
            metavar="STR",
            panel="Input",
            help="Fastq file or list of Fastq separated by a comma containing the reverse reads to be aligned or their "
            "corresponding bam files. Forward and reverse reads need to have the same identifier (read names). Not "
            "needed from pairs files.",
        ),
        assembly_option(ASSEMBLY_INDEX_HELP),
        _option(
            "-S",
            "--start",
            type=click.Choice(STARTS),
            default="fastq",
            show_default=True,
            panel="Input",
            help="Start stage of the pipeline. [cyan]metacool[/] builds the network from a .metacool contact map "
            "(see --forward) instead of reads, alignments or pairs.",
        ),
        _option(
            "-d",
            "--depth",
            metavar="FILE",
            panel="Input",
            help="The depth.txt file from the shotgun reads used to made the assembly computed by "
            "[cyan]jgi_summarize_bam_contig_depths[/] from the metabat2 pipeline.",
        ),
    )


def alignment_options():
    return _options(
        _option(
            "-b",
            "--aligner",
            type=click.Choice(ALIGNERS),
            default="bowtie2",
            show_default=True,
            panel="Alignment",
            help="Aligner algorithm to use.",
        ),
        _option(
            "-B",
            "--aligner-mode",
            type=click.Choice(ALIGNER_MODES),
            default="normal",
            show_default=True,
            panel="Alignment",
            help="Mode of alignment from hicstuff. [cyan]cutsite[/] needs --enzyme: without it, the reads are "
            "aligned iteratively.",
        ),
        _option(
            "-e",
            "--enzyme",
            metavar="STR",
            panel="Alignment",
            help="The list of restriction enzyme used to digest the contigs separated by a comma. Example: HpaII,MluCI.",
        ),
        _option(
            "-q",
            "--min-quality",
            type=click.INT,
            default=30,
            show_default=True,
            metavar="INT",
            panel="Alignment",
            help="Threshold of quality necessary to considered a read properly aligned.",
        ),
    )


def network_options(normalization_short):
    """Options of the network step; the short alias of --normalization is -n
    in network, -m in pipeline."""
    return _options(
        _option(
            normalization_short,
            "--normalization",
            type=click.Choice(NORMALIZATIONS),
            default="empirical_hit",
            show_default=True,
            # The choices are listed in the help: they are too long for the
            # value column.
            metavar="STR",
            panel="Network",
            help=f"Normalization of the contact counts, one of {', '.join(NORMALIZATIONS)}. [cyan]None[/] keeps the raw "
            "counts; [cyan]RS[/] and [cyan]theoritical_hit[/] need --enzyme, [cyan]abundance[/] and "
            "[cyan]theoritical_hit[/] need --depth.",
        ),
        _option(
            "-E",
            "--edge",
            type=click.INT,
            default=0,
            show_default=True,
            metavar="INT",
            panel="Network",
            help="Distance of the edge region in base pair on the contigs where the mapping reads are not considered as "
            "inter contigs. Not available with --start metacool.",
        ),
    )


def partition_options(iterations, overlap, iterations_help, overlap_help):
    """Options of the Louvain or Leiden partition of the network."""
    return _options(
        _option(
            "-A",
            "--algorithm",
            type=click.Choice(ALGORITHMS),
            default="louvain",
            show_default=True,
            panel="Partition",
            help="Algorithm to use to partition the network.",
        ),
        _option(
            "-i",
            "--iterations",
            type=click.INT,
            default=iterations,
            show_default=True,
            metavar="INT",
            panel="Partition",
            help=iterations_help,
        ),
        _option(
            "-O",
            "--overlap",
            type=click.INT,
            default=overlap,
            show_default=True,
            metavar="INT",
            panel="Partition",
            help=overlap_help,
        ),
        _option(
            "-r",
            "--res-param",
            type=click.FLOAT,
            default=1.0,
            show_default=True,
            metavar="FLOAT",
            panel="Partition",
            help="Resolution parameter to use for the Leiden algorithm.",
        ),
        _option(
            "-s",
            "--size",
            type=click.INT,
            default=500000,
            show_default=True,
            metavar="INT",
            panel="Partition",
            help="Threshold size to keep bins in base pair.",
        ),
        _option(
            "-C",
            "--cluster-matrix",
            is_flag=True,
            panel="Partition",
            help="If enabled, save the clustering matrix.",
        ),
    )


def outdir_option(default=None, help="The output directory. Default: current directory."):
    return _option("-o", "--outdir", default=default, show_default=default is not None, metavar="DIR", help=help)


def prefix_option():
    return _option(
        "-p",
        "--prefix",
        metavar="STR",
        help="Prefix to use for fasta files. By default just 'metator_', otherwise 'STR_metator_'. Only letters, "
        "digits, '_' and '-'.",
    )


def force_option(help):
    return _option("-F", "--force", is_flag=True, help=help)


def tmpdir_option(default="./tmp", help="Temporary directory."):
    return _option("-T", "--tmpdir", default=default, show_default=default is not None, metavar="DIR", help=help)


def threads_option(help):
    return _option("-t", "--threads", type=click.INT, default=1, show_default=True, metavar="INT", help=help)


def no_clean_up_option(short="-N", help="Do not remove temporary files."):
    return _option(short, "--no-clean-up", is_flag=True, help=help)


# ---------------------------------------------------------------------------
# CLI group
# ---------------------------------------------------------------------------


@click.group("metator", context_settings={"help_option_names": ["-h", "--help"]})
@click.rich_config(help_config=click.RichHelpConfiguration(text_markup="rich", show_arguments=True, max_width=100))
@click.command_panel("Binning", commands=["network", "partition", "validation", "pipeline"])
@click.command_panel("Other commands", commands=["qc", "contactmap", "metacool", "scaffold", "pairs", "mge"])
@click.version_option(__version__, "-v", "--version", prog_name="metator", message="%(version)s")
def cli():
    """MetaHiC pipeline for generating and manipulating MAGs.

    The binning runs in three steps, in this order: [cyan]network[/],
    [cyan]partition[/] and [cyan]validation[/]; [cyan]pipeline[/] runs them
    all in one command line.
    """


# ---------------------------------------------------------------------------
# network
# ---------------------------------------------------------------------------


@cli.command(short_help="Generate the metaHiC contig network from fastq, bam, pairs or .metacool files and normalize it.")
@reads_options(forward_required=True)
@alignment_options()
@network_options("-n")
@_option(
    "-s",
    "--self-contacts",
    is_flag=True,
    panel="Network",
    help="If enabled, count alignments between a contig and itself (intracontigs contacts).",
)
@outdir_option(help="The output directory to write the bam files the network and contig data into. Default: current directory.")
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated for the alignment.")
@no_clean_up_option()
def network(
    forward,
    reverse,
    assembly,
    start,
    depth,
    aligner,
    aligner_mode,
    enzyme,
    min_quality,
    normalization,
    edge,
    self_contacts,
    outdir,
    tmpdir,
    threads,
    no_clean_up,
):
    """Generate the metaHiC contig network.

    Generates a network file (in edgelist form) from either fastq, bam or
    pairs files. For both starting inputs, forward and reverse reads or
    alignments need to be in separated files with the same read names. In
    the output network, contigs are the network nodes and the edges are the
    contact counts.

    If multiple fastq are given, they are aligned separately and the hit
    file has a number of hits for each sample.

    With [cyan]--start metacool[/], the network is built from a .metacool
    contact map of the contigs, whose libraries are summed (see the
    README): that of the pairs it was made of, but for the pairs with a read
    on the last base of a contig or at position 0, which it keeps. --edge
    must then be 0.

    Note: the network is in a strict barebone form so that it can be reused
    and imported quickly into other applications etc. Verbose information
    about every single node in the network is written on a 'contig data'
    file.
    """
    # Private temporary directory inside --tmpdir.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir:
        # Defined the output directory and output file names.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)

        # Enable file logging
        log_file = log_file_name(outdir, "network")
        generate_log_header(log_file)

        # Check if forward and reverse arguments are given:
        check_reverse_reads(start, aligner, reverse)
        check_metacool_start(start, forward, edge)
        aligner_mode = aligner_mode_for(aligner_mode, enzyme)
        check_normalization(normalization, enzyme, depth)

        # Check what is the reference. If a fasta is given build the index. If a
        # bowtie2 index is given, retreive the fasta.
        index, fasta = reference_files(assembly, aligner, start, tmp_dir)

        # Print information of the workflow:
        logger.info("Aligner algorithm: %s", aligner)
        logger.info("Enzyme: %s", enzyme)
        logger.info("Normalization: %s", normalization)
        logger.info("Minimum mapping quality: %s", min_quality)

        # Do not align if pair or metacool start
        if start in ["pair", "metacool"]:
            # A map has no pairs files: its libraries are summed.
            alignment_files = forward.split(",") if start == "pair" else []
            contig_data, hit_data = mtn.create_contig_data(fasta, len(alignment_files) or 1, depth, enzyme)
        else:
            # Align pair-end reads with bowtie2
            alignment_files, contig_data, hit_data = mta.get_contact_pairs(
                forward,
                reverse,
                index,
                fasta,
                aligner,
                aligner_mode,
                min_quality,
                start,
                depth,
                enzyme,
                outdir,
                tmp_dir,
                threads,
            )

        # Build the network
        contig_network(
            start, forward, alignment_files, contig_data, hit_data, edge, outdir, tmp_dir, threads, normalization,
            self_contacts,
        )  # fmt: skip

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# partition
# ---------------------------------------------------------------------------


@cli.command(short_help="Partition the metaHiC network with the Louvain or Leiden algorithm.")
@assembly_option(ASSEMBLY_PARTITION_HELP)
@_option(
    "-c",
    "--contigs",
    required=True,
    metavar="FILE",
    panel="Input",
    help="The path to the tsv file containing the data of the contigs (ID, Name, Length, GC content, Hit, Coverage, "
    "Restriction Site).",
)
@_option(
    "-n",
    "--network",
    required=True,
    metavar="FILE",
    panel="Input",
    help="Path to the file containing the network information from the meta HiC experiment computed by the network command.",
)
@partition_options(
    iterations=100,
    overlap=80,
    iterations_help="Number of iterations of Louvain or Leiden.",
    overlap_help="Hamming distance threshold to use to merge bins (percentage).",
)
@outdir_option(default="./", help="Path to the directory to write the output.")
@prefix_option()
@force_option("If enabled, would remove directory of overlapping bins in the output directory.")
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated for the partition.")
@no_clean_up_option()
def partition(
    assembly,
    contigs,
    network,
    algorithm,
    iterations,
    overlap,
    res_param,
    size,
    cluster_matrix,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
):
    """Partition the network using the Louvain or Leiden algorithm.

    Partition the network file using iteratively the Louvain or Leiden
    algorithm. Then looks for 'cores' bins that are constituted by the
    contigs which are always in the same cluster at each iterations. Then
    using the Hamming distance from these cores bins, merge the bins with a
    small Hamming distance (close to 1), bigger than the percentage
    (overlapping parameter) given.

    It will also update the file to integrate the bins information of the
    contigs.

    Furthermore, both Leiden and Louvain algorithm are available here.
    However, the benchmark made show that here the Louvain algorithm have
    better performance and is faster on seawater and gut metagenomic
    samples.
    """
    # Private temporary directory inside --tmpdir.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir:
        # Defined the output directory.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)
        fasta_dir = join(outdir, "overlapping_bin")
        make_output_dir(fasta_dir, force)

        # Enable file logging
        log_file = log_file_name(outdir, "partition")
        generate_log_header(log_file)

        # Create prefix.
        bins_prefix = bin_prefix(prefix)

        # Partition the network
        mtp.partition(
            algorithm,
            assembly,
            cluster_matrix,
            contigs,
            iterations,
            network,
            outdir,
            fasta_dir,
            overlap / 100,
            res_param,
            size,
            tmp_dir,
            threads,
            bins_prefix,
        )

        # Delete pyfastx index (only made if a bin passed the size threshold):
        mio.remove_pyfastx_index(assembly)

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


@cli.command(short_help="Validate the bins with miComplete and decontaminate them with a recursive partition.")
@assembly_option(ASSEMBLY_PARTITION_HELP)
@_option(
    "-c",
    "--contigs",
    required=True,
    metavar="FILE",
    panel="Input",
    help="The path to the file containing the data of the contigs from the partition step (13 columns).",
)
@_option(
    "-f",
    "--fasta",
    required=True,
    metavar="DIR",
    panel="Input",
    help="Path to the directory containing the input fasta files of the bins to decontaminate.",
)
@_option(
    "-n",
    "--network",
    required=True,
    metavar="FILE",
    panel="Input",
    help="Path to the file containing the network information from the meta HiC experiment computed by the network command.",
)
@partition_options(
    iterations=10,
    overlap=90,
    iterations_help="Number of recursive iterations of Louvain or Leiden.",
    overlap_help="Hamming distance threshold to use to merge bins (percentage).",
)
@outdir_option(default="./", help="Path to the directory to write the output.")
@prefix_option()
@force_option("If enabled, would remove directory of recursive bins in the output directory.")
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated for the partition.")
@no_clean_up_option()
def validation(
    assembly,
    contigs,
    fasta,
    network,
    algorithm,
    iterations,
    overlap,
    res_param,
    size,
    cluster_matrix,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
):
    """Use miComplete to validate the bins.

    Recursively decontaminate the contaminated bins evaluated by
    miComplete. The command returns new decontaminated fasta and summary
    files of the decontamination.

    Only bins with more than 50% completion and 5% contamination are subject
    to the recursive step. If the recursive step gave worst results than the
    first (decrease of the completion with no decrease of the
    contamination), it will keep the original bin.
    """
    # Launch miComplete to evaluate the completion and the contamination. If
    # asked rerun Louvain to try to reduce the contamination, rerun miComplete
    # if the contamination decrease without a huge decrease of the completion
    # keep the new bins. Otherwise go back to the old state.

    # Private temporary directory inside --tmpdir.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir:
        # Defined the output directory and output file names.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)
        recursive_fasta_dir = join(outdir, "recursive_bin")
        final_fasta_dir = join(outdir, "final_bin")
        make_output_dir(recursive_fasta_dir, force)
        make_output_dir(final_fasta_dir, force)

        # Enable file logging
        log_file = log_file_name(outdir, "validation")
        generate_log_header(log_file)

        # Create prefix.
        bins_prefix = bin_prefix(prefix)

        mtv.recursive_decontamination(
            algorithm,
            assembly,
            cluster_matrix,
            contigs,
            final_fasta_dir,
            fasta,
            iterations,
            network,
            outdir,
            overlap / 100,
            recursive_fasta_dir,
            res_param,
            size,
            tmp_dir,
            threads,
            bins_prefix,
        )

        # Delete pyfastx index (only made if a bin was re-partitioned):
        mio.remove_pyfastx_index(assembly)

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------


@cli.command(short_help="Run network, partition and validation (and scaffold the bins) in one command line.")
@reads_options(forward_required=False)
@alignment_options()
@network_options("-m")
@partition_options(
    iterations=100,
    overlap=80,
    iterations_help="Number of iterations of Louvain or Leiden for the partition step.",
    overlap_help="Hamming distance threshold to use to merge bins for the first step (percentage).",
)
@_option(
    "-j",
    "--rec-iter",
    type=click.INT,
    default=10,
    show_default=True,
    metavar="INT",
    panel="Validation",
    help="Number of iterations of Louvain or Leiden for the recursive step.",
)
@_option(
    "-P",
    "--rec-overlap",
    type=click.INT,
    default=90,
    show_default=True,
    metavar="INT",
    panel="Validation",
    help="Hamming distance threshold to use to merge bins for the recursive step (percentage).",
)
@_option(
    "-v",
    "--scaffold",
    is_flag=True,
    panel="Scaffolding",
    help="If enabled, it will scaffold the genomes at the end.",
)
@_option(
    "-J",
    "--junctions",
    default="NNNNN",
    show_default=True,
    metavar="STR",
    panel="Scaffolding",
    help="Sequences to use as junction between contigs.",
)
@outdir_option(help="The output directory to write the bam files the network and contig data into. Default: current directory.")
@prefix_option()
@force_option("If enabled, would remove the directories of bins in the output directory.")
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated for the alignment.")
@no_clean_up_option()
@_option(
    "--no-metacool",
    is_flag=True,
    help="Do not write the contact map of the contigs, grouped by bin, as <prefix>.metacool (metator.metacool "
    "without --prefix) at the end.",
)
def pipeline(
    forward,
    reverse,
    assembly,
    start,
    depth,
    aligner,
    aligner_mode,
    enzyme,
    min_quality,
    normalization,
    edge,
    algorithm,
    iterations,
    overlap,
    res_param,
    size,
    cluster_matrix,
    rec_iter,
    rec_overlap,
    scaffold,
    junctions,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
    no_metacool,
):
    """Launch the full metaTOR pipeline.

    Partition the contigs from the assembly in bins from the metaHiC reads.

    It's possible to start from the fastq, the bam or the pairs files, or
    from a .metacool contact map of the contigs ([cyan]--start metacool[/],
    see [cyan]metator network[/]). A validation step decontaminates the bins
    when it's necessary; it is the critical step for memory usage (~40G).

    At the end, the contact map of the contigs, grouped by bin, is written
    as a .metacool file (see [cyan]metator metacool[/] and the README),
    unless --no-metacool is given: from the pairs files, or from the input
    map of a run that starts from a .metacool file, which it may replace. A
    run from a .metacool file cannot scaffold.
    """
    # Private temporary directory inside --tmpdir.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir:
        # Defined the output directory and output file names.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)

        # Enable file logging
        log_file = log_file_name(outdir, "pipeline")
        generate_log_header(log_file)

        # Define variable
        overlapping_parameter = overlap / 100
        recursive_overlapping_parameter = rec_overlap / 100

        # Check the arguments, and the assembly below, before replacing the
        # bins of an earlier run (--force): a refused run leaves them as they
        # were.
        aligner_mode = aligner_mode_for(aligner_mode, enzyme)
        check_normalization(normalization, enzyme, depth)
        check_metacool_start(start, forward, edge, scaffold)
        check_reverse_reads(start, aligner, reverse)
        bins_prefix = bin_prefix(prefix)

        # Manage start point.
        start_step = {"fastq": 1, "bam": 2, "pair": 3, "metacool": 3}[start]

        # Print information of the workflow:
        if start_step == 1:
            logger.info(f"Minimum mapping quality: {min_quality}")
        if start_step <= 2:
            logger.info(f"Enzyme: {enzyme}")
            logger.info(f"Normalization: {normalization}")
        logger.info(f"Aligner algorithm: {aligner}")
        logger.info(f"Partition algorithm: {algorithm}")
        logger.info(f"Partition iterations: {iterations}")
        logger.info(f"Overlapping parameter: {overlapping_parameter}")
        logger.info(f"Recursive partition iterations: {rec_iter}")
        logger.info(f"Recursive overlapping parameter: {recursive_overlapping_parameter}\n")

        # The bins of an earlier run: refused now without --force, before
        # the index of the assembly is built, and replaced only once the
        # assembly is checked.
        overlapping_fasta_dir = join(outdir, "overlapping_bin")
        recursive_fasta_dir = join(outdir, "recursive_bin")
        final_fasta_dir = join(outdir, "final_bin_unscaffold")
        scaffold_fasta_dir = join(outdir, "scaffold_bin")
        bin_dirs = [overlapping_fasta_dir, recursive_fasta_dir, final_fasta_dir] + ([scaffold_fasta_dir] if scaffold else [])
        for path in bin_dirs:
            refuse_output_dir(path, force)

        # Check what is the reference. If a fasta is given build the index. If a
        # bowtie2 index is given, retreive the fasta.
        index, fasta = reference_files(assembly, aligner, start, tmp_dir)
        for path in bin_dirs:
            make_output_dir(path, force)

        # Run the whole workflow
        if start_step <= 2:
            # Align pair-end reads with bowtie2
            alignment_files, contig_data, hit_data = mta.get_contact_pairs(
                forward,
                reverse,
                index,
                fasta,
                aligner,
                aligner_mode,
                min_quality,
                start,
                depth,
                enzyme,
                outdir,
                tmp_dir,
                threads,
            )
        else:
            # A map has no pairs files: its libraries are summed.
            alignment_files = forward.split(",") if start == "pair" else []
            contig_data, hit_data = mtn.create_contig_data(fasta, len(alignment_files) or 1, depth, enzyme)
        # Build the network
        network_file, contigs_data_file = contig_network(
            start, forward, alignment_files, contig_data, hit_data, edge, outdir, tmp_dir, threads, normalization, False
        )

        # Partition the network
        clustering_matrix_partition_file, contigs_data_file = mtp.partition(
            algorithm,
            fasta,
            cluster_matrix,
            contigs_data_file,
            iterations,
            network_file,
            outdir,
            overlapping_fasta_dir,
            overlapping_parameter,
            res_param,
            size,
            tmp_dir,
            threads,
            bins_prefix,
        )

        # remove contig_data_network if not an input
        if start_step <= 2:
            os.remove(join(outdir, "contig_data_network.txt"))

        # Launch validation.
        clustering_matrix_recursive_file = mtv.recursive_decontamination(
            algorithm,
            fasta,
            cluster_matrix,
            contigs_data_file,
            final_fasta_dir,
            overlapping_fasta_dir,
            rec_iter,
            network_file,
            outdir,
            recursive_overlapping_parameter,
            recursive_fasta_dir,
            res_param,
            size,
            tmp_dir,
            threads,
            bins_prefix,
        )

        if cluster_matrix:
            # Make the sum with the partition clustering matrix and save it.
            clustering_matrix = load_npz(clustering_matrix_partition_file + ".npz")
            clustering_matrix_recursive = load_npz(clustering_matrix_recursive_file + ".npz")
            clustering_matrix = ((clustering_matrix + clustering_matrix_recursive) / 2).tocoo()
            save_npz(join(outdir, "clustering_matrix"), clustering_matrix)

        # Remove contig_data_partition file
        os.remove(join(outdir, "contig_data_partition.txt"))

        # Launch the scaffold
        if scaffold:
            bin_summary = mio.read_bin_summary(join(outdir, "bin_summary.txt"))
            task = partial(
                mts.parallel_scaffold,
                final_fasta_dir=final_fasta_dir,
                alignment_files=alignment_files,
                scaffold_fasta_dir=scaffold_fasta_dir,
                junctions=junctions,
            )
            # Workers started by forkserver (Linux's default from Python 3.14)
            # or spawn (macOS) do not inherit the log file handler.
            with mp.Pool(processes=threads, initializer=mtl.set_file_handler, initargs=(log_file,)) as pool:
                pool.map(task, bin_summary.index)

        # Delete pyfastx index:
        os.remove(fasta + ".fxi")

        # Write the contact map last, while the pairs files are still there:
        # from the input map for a run that started from one.
        if not no_metacool:
            pipeline_metacool(
                alignment_files, outdir, prefix, assembly, threads, metacool=forward if start == "metacool" else None
            )

    # Delete the pairs files made from the reads or the alignments, and
    # their pypairix index.
    if not no_clean_up and start_step < 3:
        for pair_file in alignment_files:
            os.remove(pair_file)
            if os.path.exists(f"{pair_file}.px2"):
                os.remove(f"{pair_file}.px2")

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# qc
# ---------------------------------------------------------------------------


@cli.command(short_help="Generate some quality checks on the output of metaTOR.")
@pairsfile_argument()
@_option(
    "-a",
    "--assembly",
    required=True,
    metavar="FILE",
    panel="Input",
    help="Path to the fasta file containing the contigs of interest. Could be the whole or the extracted contigs of one bin.",
)
@_option(
    "-e",
    "--enzyme",
    required=True,
    metavar="STR",
    panel="Input",
    help="The list of restriction enzyme used to digest the contigs separated by a comma. Example: HpaII,MluCI.",
)
@_option(
    "-O",
    "--metator-dir",
    metavar="DIR",
    panel="Input",
    help="Output directory from metator pipeline. If set, contigs data and bin summary will be found automatically.",
)
@_option(
    "-b",
    "--bin-summary",
    metavar="FILE",
    panel="Input",
    help="Path to the bin_summary.txt file from MetaTOR output.",
)
@_option(
    "-c",
    "--contig-data",
    metavar="FILE",
    panel="Input",
    help="Path to the contig_data_final.txt file from MetaTOR output.",
)
@_option(
    "-t",
    "--threshold",
    metavar="STR",
    panel="Quality check",
    help="Hicstuff religation and loop thresholds. Two integers separated by a comma.",
)
@_option("-P", "--plot", is_flag=True, panel="Quality check", help="If enabled, display some plots.")
@_option("-p", "--prefix", metavar="STR", panel="Quality check", help="Name of the sample to add on plot and files.")
@outdir_option(help="Directory to save output plots and log. Default: ./metator_qc or ./metator_qc_<prefix>.")
@tmpdir_option(default=None, help="Temporary directory to save pairs files. Default to ./tmp in the current directory.")
@no_clean_up_option(short="-n", help="If enabled, intermediary files will be kept.")
def qc(
    pairsfile,
    assembly,
    enzyme,
    metator_dir,
    bin_summary,
    contig_data,
    threshold,
    plot,
    prefix,
    outdir,
    tmpdir,
    no_clean_up,
):
    """Generate some quality checks on the output of metaTOR.

    The quality check is done on the high quality MAGs and only contigs with
    a size bigger than 100kb. The goal is to assess HiC quality metrics on
    genomes we kinda know to have a estimation of the quality of the
    library.

    If the plot are set it will return some plot of the event distribution
    on the high quality contigs.
    """
    # Private temporary directory inside --tmpdir.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir:
        # Set output directory.
        if not outdir:
            outdir = join(".", f"metator_qc_{prefix}" if prefix else "metator_qc")
        os.makedirs(outdir, exist_ok=True)

        # Enable file logging
        log_file = log_file_name(outdir, "qc")
        generate_log_header(log_file)

        # Try to find input files.
        metator_dir = metator_dir or "."
        if bin_summary:
            bin_summary_file = bin_summary
        else:
            bin_summary_file = join(metator_dir, "bin_summary.txt")
            if not exists(bin_summary_file):
                logger.error(
                    "Please give a bin summary file, or the ouput directory of metator or launch this command from the output directory of metator."
                )
                logger.error("You can also check that a bin_summary.txt is present in metator output directory")
                raise FileNotFoundError
        if contig_data:
            contig_data_file = contig_data
        else:
            contig_data_file = join(metator_dir, "contig_data_final.txt")
            if not exists(contig_data_file):
                logger.error(
                    "Please give a contig data file, or the ouput directory of metator or launch this command from the output directory of metator."
                )
                logger.error("You can also check that a contig_data_final.txt is present in metator output directory")
                raise FileNotFoundError

        # Launch quality check
        mtq.quality_check(
            contig_data_file,
            bin_summary_file,
            assembly,
            list(pairsfile),
            outdir,
            tmp_dir,
            prefix or "metator_qc",
            plot,
            enzyme=enzyme.split(",") if enzyme else enzyme,
            threshold=threshold.split(",") if threshold else threshold,
        )

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# contactmap
# ---------------------------------------------------------------------------


@cli.command(short_help="Generate the HiC contact map of one metaTOR object (contig or bin).")
@pairsfile_argument()
@_option(
    "-a",
    "--assembly",
    required=True,
    metavar="FILE",
    panel="Input",
    help="Path to the fasta file containing the contigs of interest. Could be the whole or the extracted contigs of one bin.",
)
@_option(
    "-c",
    "--contig-data",
    required=True,
    metavar="FILE",
    panel="Input",
    help="Path to the contig_data_final.txt file from MetaTOR output.",
)
@_option(
    "-e",
    "--enzyme",
    required=True,
    metavar="STR",
    panel="Input",
    help="The list of restriction enzyme used to digest the contigs separated by a comma. Example: HpaII,MluCI.",
)
@_option(
    "-n",
    "--name",
    required=True,
    metavar="STR",
    panel="Input",
    help="Name of the metator object or its numerical ID for the core bin, overlapping bin or recursive bin objects. "
    'Example: "NODE_1", "MetaTOR_1_0" or 8.',
)
@_option(
    "-O",
    "--object",
    "metator_object",
    default="final_bin",
    show_default=True,
    metavar="STR",
    panel="Input",
    help='Type of metator object used to build a contact map. Either "contig", "core_bin", "overlapping_bin", '
    '"recursive_bin", "final_bin" or the name of another column of the contig data file, such as "Other".',
)
@_option(
    "-m",
    "--mat-fmt",
    type=click.Choice(MATRIX_FORMATS),
    default="graal",
    show_default=True,
    panel="Contact map",
    help='The format of the output sparse matrix: "bg2" for 2D Bedgraph format, "cool" for Mirnylab\'s cooler '
    'software, or "graal" for graal-compatible plain text COO format.',
)
@_option(
    "-f",
    "--filter",
    "filter_events",
    is_flag=True,
    panel="Contact map",
    help="Filter out spurious 3C events (loops and uncuts) using hicstuff filter. For more informations, see Cournac "
    "et al. BMC Genomics, 2012.",
)
@_option(
    "-D",
    "--pcr-dup",
    is_flag=True,
    panel="Contact map",
    help="Filter out PCR duplicates based on read positions.",
)
@_option(
    "-s",
    "--min-size",
    type=click.INT,
    default=5000,
    show_default=True,
    metavar="INT",
    panel="Contact map",
    help="Minimum size threshold to consider contigs.",
)
@outdir_option(help="Output directory. Default creates a new directory contact_map_<name> in the current directory.")
@force_option("Write files even if the output files already exists.")
@tmpdir_option(
    default=None,
    help='Directory for storing intermediary files and temporary files. Default creates a "tmp" folder in the current directory.',
)
@threads_option(help="Number of threads to allocate.")
@no_clean_up_option(help="If enabled, intermediary files will be kept.")
def contactmap(
    pairsfile,
    assembly,
    contig_data,
    enzyme,
    name,
    metator_object,
    mat_fmt,
    filter_events,
    pcr_dup,
    min_size,
    outdir,
    force,
    tmpdir,
    threads,
    no_clean_up,
):
    """Generate a HiC contact map of a MetaTOR object from the pairs files,
    the contig data file and a fasta file containing the contigs of
    interest.

    Generates the Hi-C matrix of one MetaTOR object from the pair alignment
    file of metator. MetaTOR object integrated are contigs, and core,
    overlapping, recursive, final or personalized bins. For the personalized
    bins you should add a column in the contig_data_final.txt file with a
    header called "Other".

    Do not display the matrix, you have to run another pipeline to display it
    such hicstuff view or cooler show depending on the matrix format chosen.
    It's also possible to try to scaffold the bin using instagraal pipeline
    using the output.
    """
    # Defined the output directory first: the log file is written there.
    outdir = outdir or join(".", "contact_map_" + name)
    os.makedirs(outdir, exist_ok=True)

    # Enable file logging
    log_file = log_file_name(outdir, "contactmap")
    generate_log_header(log_file)

    # Private temporary directory inside --tmpdir.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir:
        mtc.generate_contact_map(
            assembly,
            contig_data,
            enzyme,
            name,
            list(pairsfile),
            outdir,
            tmp_dir,
            filter_events,
            force,
            mat_fmt,
            metator_object,
            min_size,
            pcr_dup,
            threads,
        )

    if not no_clean_up:
        # Delete pyfastx index:
        os.remove(assembly + ".fxi")

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# metacool
# ---------------------------------------------------------------------------


@cli.command(short_help="Write the contact map of the contigs, grouped by bin, as a .metacool file.")
@pairsfile_argument()
@_option(
    "-c",
    "--contig-data",
    required=True,
    metavar="FILE",
    panel="Input",
    help="Path to the contig_data_final.txt file from MetaTOR output: the contigs, their final bins and what "
    "metaTOR measured on them.",
)
@_option(
    "-b",
    "--bin-summary",
    metavar="FILE",
    panel="Input",
    help="Path to the bin_summary.txt file from MetaTOR output: its columns are added to the genomes table.",
)
@_option("-o", "--outfile", required=True, metavar="FILE", help="Path to the .metacool file to write.")
@_option("-n", "--assembly-name", metavar="STR", help="Name of the assembly, stored in the file.")
@threads_option(help="Number of threads binning the pairs.")
@force_option("Overwrite the output file if it exists.")
def metacool(pairsfile, contig_data, bin_summary, outfile, assembly_name, threads, force):
    """Write the contact map of a metaTOR run as a .metacool file.

    A .metacool file is the contact map of metacooler: a cooler whose
    chromosomes are the contigs, one bin per contig, grouped by genome: one
    genome per final bin of the contig data file, and [cyan]__unbinned__[/]
    for the other contigs, with all their contacts. The README says what
    else it holds. The pairs files of several libraries are summed.

    The contig and bin names must be ASCII, and no bin may be named as
    metacooler marks a contig in no bin (NA, NaN, N/A, N/D, -1 or
    __unbinned__ in any case, and nd in any case other than ND, which marks
    a contig in no bin). [cyan]pipeline[/] writes this file at the end of a
    run.
    """
    # Generate log
    log_file = log_file_name("", "metacool")
    generate_log_header(log_file)

    try:
        mtmc.write_metacool(
            list(pairsfile),
            contig_data,
            outfile,
            bin_summary=bin_summary,
            assembly_name=assembly_name,
            threads=threads,
            force=force,
            metadata=run_metadata(),
        )
    except mtmc.MetacoolError as error:
        # The reason in the log file of the run, then the error.
        logger.error(str(error))
        raise

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# scaffold
# ---------------------------------------------------------------------------


@cli.command(short_help="Scaffold a metaTOR bin from pairs files.")
@pairsfile_argument()
@_option("-b", "--bin-name", required=True, metavar="STR", panel="Input", help="Name of the bin to scaffold.")
@_option("-i", "--input-fasta", required=True, metavar="FILE", panel="Input", help="Path to the input fasta.")
@_option(
    "-T",
    "--threshold",
    type=click.FLOAT,
    default=0.05,
    show_default=True,
    metavar="FLOAT",
    panel="Scaffolding",
    help="Threshold score to consider an association between two contigs.",
)
@_option(
    "-w",
    "--window-size",
    type=click.INT,
    default=5000,
    show_default=True,
    metavar="INT",
    panel="Scaffolding",
    help="Size of the window in base pair to use as edge window to scaffold.",
)
@_option(
    "-j",
    "--junctions",
    default="NNNNN",
    show_default=True,
    metavar="STR",
    panel="Scaffolding",
    help="Sequences to use as junction between contigs.",
)
@_option(
    "-o",
    "--out-fasta",
    metavar="FILE",
    help="Path to write the output fasta. Default in the current directory: <bin_name>_scaffolded.fa.",
)
@_option(
    "-O",
    "--out-frags",
    metavar="FILE",
    help="Path to write the fragments information. Default in the current directory: <bin_name>_info_frags.txt.",
)
@threads_option(help="Numbers of threads to allocate if pairs need to be sorted.")
def scaffold(pairsfile, bin_name, input_fasta, threshold, window_size, junctions, out_fasta, out_frags, threads):
    """Scaffold a bin from metaTOR.

    Use the pairs to build a scaffold from the binned contigs. The
    scaffolding is based on contact between the edges of the contigs. There
    are normalized based on intra contacts of the contigs except for the
    small contigs where a low score is given.
    """
    # Generate log
    log_file = log_file_name("", "scaffold")
    generate_log_header(log_file)

    # Get the scaffolds
    mts.get_scaffolds(
        bin_name,
        input_fasta,
        list(pairsfile),
        join(".", f"{bin_name}_scaffolded.fa") if out_fasta is None else out_fasta,
        join(".", f"{bin_name}_info_frags.txt") if out_frags is None else out_frags,
        threshold=threshold,
        threads=threads,
        window_size=window_size,
        junctions=junctions,
    )

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# pairs
# ---------------------------------------------------------------------------


@cli.command(short_help="Sort, compress and index pairs files.")
@pairsfile_argument()
@_option("-F", "--force", is_flag=True, help="Write files even if the output files already exists.")
@_option(
    "-r",
    "--remove",
    is_flag=True,
    help="Remove the input file at the end to keep only the sorted, compressed and indexed pairs file.",
)
@threads_option(help="Numbers of thread to allocate.")
def pairs(pairsfile, force, remove, threads):
    """Sort, compress and index pairs files for faster access to the data.

    Sort the pairs files (compressed or not) by chr1, chr2, pos1 and pos2
    using hicstuff. Compress them using bgzip. Index them using pypairix.
    Each output goes next to its input, named after it without its last
    extension: sample.pairs gives sample_sorted.pairs.gz and its index
    sample_sorted.pairs.gz.px2, sample.pairs.gz gives
    sample.pairs_sorted.pairs.gz.
    """
    # Generate log
    log_file = log_file_name("", "pairs")
    generate_log_header(log_file)

    # Iterates on pairfiles given.
    for pairs_file in pairsfile:
        logger.info(f"Processing {pairs_file}...")
        # Run the sort/compress/index command. Sort's temporary files go to
        # the working directory: /tmp can be small, or in memory.
        mio.sort_and_index_pairs(pairs_file, threads=threads, remove=remove, force=force, tmp_dir=".")

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# mge
# ---------------------------------------------------------------------------


@cli.command(short_help="Bin the MGE contigs in MGE MAGs from the metaHiC pairs.")
@pairsfile_argument()
@_option("-a", "--assembly", required=True, metavar="FILE", panel="Input", help="Path to the assembly fasta file.")
@_option("-c", "--contigs", required=True, metavar="FILE", panel="Input", help="Path to the MetaTOR contig data file.")
@_option("-b", "--binning", required=True, metavar="FILE", panel="Input", help="Path to the anvio binning file.")
@_option(
    "-m",
    "--mges",
    required=True,
    metavar="FILE",
    panel="Input",
    help="Path to the file listing MGE-containing contigs.",
)
@_option(
    "-s",
    "--threshold-bin",
    type=click.FLOAT,
    default=0.8,
    show_default=True,
    metavar="FLOAT",
    panel="Binning",
    help="Threshold to use for binning.",
)
@_option("-r", "--random", is_flag=True, panel="Binning", help="If enabled, make a random binning.")
@outdir_option(help="Path to the output directory where the output will be written. Default: current directory.")
@tmpdir_option(help="Path to temporary directory.")
@no_clean_up_option(help="If enabled, intermediary files will be kept.")
def mge(pairsfile, assembly, contigs, binning, mges, threshold_bin, random, outdir, tmpdir, no_clean_up):
    """Bin MGE contigs.

    Bin the MGE contigs in MGE MAGs using pairs, host MAG sequence and
    information.

    Command will return the updated MGE data and the MGE fasta.

    The MGE fasta will contain one entry by MGE MAG, with 180 "N" spacers
    between contigs.
    """
    # Generate log
    log_file = log_file_name("", "mge")
    generate_log_header(log_file)

    # Private temporary directory inside --tmpdir.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir:
        # Defined the output directory and output file names.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)

        # Import the files
        binning_result = mio.import_anvio_binning(binning)
        mges_list = mio.import_mges_contigs(mges)
        contigs_data, mges_list_id = mio.import_contig_data_mges(contigs, binning_result, mges_list)

        # Run the mges binning
        mtm.mge_binning(
            fasta_mges_contigs=assembly,
            contigs_data=contigs_data,
            mges_list_id=mges_list_id,
            out_dir=outdir,
            pairs_files=list(pairsfile),
            tmp_dir=tmp_dir,
            threshold_bin=threshold_bin,
            random=random,
        )

    if not no_clean_up:
        # Delete pyfastx index:
        os.remove(assembly + ".fxi")

    generate_log_footer(log_file)
