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
import signal
import sys
import threading
import time
from contextlib import ExitStack, contextmanager
from functools import partial
from os.path import dirname, exists, join

import rich_click as click
from click.core import ParameterSource
from scipy.sparse import load_npz, save_npz

import metator.align as mta
import metator.clusters as mtcl
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
STARTS = ["fastq", "bam", "pair"]
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


def check_normalization(normalization, enzyme, depth, metacool=None, read=None):
    """The normalizations that need the restriction sites or the shotgun
    coverage of the contigs: measured on the assembly with --enzyme and
    --depth or, in a run from a .metacool map (metacool), read from the
    Restriction_site and Shotgun_coverage columns of its contigs table (see
    metator.network.contig_data_from_map). Before the run reads the
    assembly, the map is checked for the column (read None); once it has,
    read lists the columns that it could read."""
    needs = [
        (("RS", "theoritical_hit"), enzyme, "enzyme", "Restriction_site"),
        (("abundance", "theoritical_hit"), depth, "depth", "Shotgun_coverage"),
    ]
    for normalizations, given, option, column in needs:
        if normalization not in normalizations or given:
            continue
        message = f'For "{normalizations[0]}" and "{normalizations[1]}" normalization, {option} is required'
        if metacool is None:
            message += "."
        elif read is None:
            try:
                columns = mtn.map_columns(metacool)
            except mtmc.MetacoolError as error:
                logger.error(str(error))
                raise
            if column in columns:
                continue
            what = mtn.MAP_COLUMNS[column][2]
            message += f", as {metacool} gives no {what} (no {column} column with a valid value)."
        elif column in read:
            continue
        else:
            what = mtn.MAP_COLUMNS[column][2]
            message += f", as the {column} column of {metacool} does not give every contig of the assembly its {what}."
        logger.error(message)
        raise ValueError


def contig_data_from_map(metacool, contig_data, normalization, enzyme, depth):
    """Read the restriction sites and the shotgun coverage of the contigs from
    the .metacool map that a run starts from, where --enzyme and --depth
    measure none (see metator.network.contig_data_from_map), then check the
    normalization against what the run has. A map that cannot be read, or
    that is not a map of the contigs of the assembly, stops the run, its
    reason logged. Returns the columns read from the map."""
    try:
        read = mtn.contig_data_from_map(metacool, contig_data, enzyme, depth)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    check_normalization(normalization, enzyme, depth, metacool, read)
    return read


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


def run_input(input_map, forward, reverse, start):
    """What a network or pipeline run starts from, as (start, forward): the
    .metacool map given as the first argument (input_map), as the metacool
    start with the map for forward, or the reads, alignments or pairs of
    -1/--forward, as --start says. A run takes one or the other; the
    command line refuses a first argument that is no map given with
    -1/--forward, as one argument too many (see Command).

    Refused before any work, the reason logged as an error: no input; a
    first argument that is not there or is no map (reads given without -1,
    for example); a map with -1/--forward, -2/--reverse or --start, even at
    its default; and a map that metacooler cannot open, with its reason (a
    .mmetacool collection rather than one of its entries, for example), as
    MetacoolError. A pipeline run checks its input before it replaces the
    clusters of an earlier run (-F)."""
    if input_map is None:
        if not forward:
            logger.error(
                "Give the run a .metacool contact map, as the first argument, or its reads, alignments or pairs, "
                "with -1/--forward (and --start)."
            )
            raise ValueError
        return start, forward
    # A file whose name holds :: is that file, not the URI of a map inside
    # one (see metator.metacool.is_map).
    path = input_map if os.path.isfile(input_map) else mtmc.map_file(input_map)
    if not os.path.isfile(path):
        reason = "no such file" if not os.path.exists(path) else f"{path} is not a file"
        message = f"Cannot read the contact map {input_map}: {reason}."
        logger.error(message)
        raise mtmc.MetacoolError(message)
    if not mtmc.is_map(input_map):
        logger.error(
            f"{input_map} is not a .metacool contact map: give reads, alignments or pairs with -1/--forward (and --start)."
        )
        raise ValueError
    if path == input_map and "::" in input_map:
        message = (
            f"Cannot read the contact map {input_map}: metacooler reads a path that holds :: as file::group. "
            "Rename the map without ::."
        )
        logger.error(message)
        raise mtmc.MetacoolError(message)
    ctx = click.get_current_context(silent=True)
    source = ctx.get_parameter_source("start") if ctx is not None else None
    explicit_start = source not in (None, ParameterSource.DEFAULT)
    given = [
        option for option, value in [("-1/--forward", forward), ("-2/--reverse", reverse), ("--start", explicit_start)] if value
    ]
    if given:
        logger.error(
            f"A run from the contact map {input_map} takes no {' or '.join(given)}: these give the reads, alignments "
            "or pairs of a run that does not start from a map."
        )
        raise ValueError
    try:
        mtmc.open_map(input_map)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    return "metacool", input_map


def check_metacool_start(start, scaffold=False):
    """A run from a map cannot scaffold: the scaffolding reads the pairs
    files."""
    if start == "metacool" and scaffold:
        logger.error("--scaffold reads the pairs files: it cannot run from a contact map.")
        raise ValueError


# The compressions that metaTOR finds from the first bytes of a file
# (metator.io.file_compression), whatever its name: it reads no zstd or lz4.
SNIFFED_CODECS = {"gz": "gzip", "bz2": "bzip2", "zip": "zip"}
# The first bytes of a zstd file and of an LZ4 frame, which metacooler reads
# (it finds the compression of a pairs file from its first bytes too,
# whatever its name), but not metaTOR.
UNREAD_BY_METATOR = {b"\x28\xb5\x2f\xfd": "zstd", b"\x04\x22\x4d\x18": "lz4"}


def check_pairs_compression(start, forward, scaffold=False):
    """Refuse, before any work, a pairs file given to a run from pairs (-1)
    that metacooler, which builds the contact maps of the contigs from it,
    or metaTOR, which reads the pairs files of --scaffold, could not read.
    Both find the compression of a pairs file from its first bytes,
    whatever its name; metaTOR reads no zstd or lz4 (UNREAD_BY_METATOR).
    Refused: an empty file; a zip archive, which metacooler cannot read;
    and, with --scaffold, a zstd or lz4 file (the network comes from the
    maps, whatever the number of libraries). The reason, with what to do, is
    logged as an error, then ValueError is raised. A file that is not there
    is left to the map, which stops on it."""
    if start != "pair":
        return
    for path in forward.split(","):
        if not os.path.isfile(path):
            continue
        with open(path, "rb") as handle:
            unread = UNREAD_BY_METATOR.get(handle.read(4))
        if os.path.getsize(path) == 0:
            reason = "is empty: give the pairs file of its library, or leave it out of -1."
        elif SNIFFED_CODECS.get(mio.file_compression(path)) == "zip":
            reason = (
                "is a zip archive, which metacooler, which builds the contact map of the contigs, cannot read: "
                "decompress it, and give the pairs file itself, plain or compressed with gzip."
            )
        elif unread and scaffold:
            reason = (
                f"is compressed with {unread}, which metacooler, which builds the contact map of the contigs, reads, "
                "but not metaTOR, which reads the pairs files for --scaffold: decompress it, or compress it with gzip "
                "(zstd and lz4 work without --scaffold)."
            )
        else:
            continue
        logger.error(f"{path} {reason}")
        raise ValueError


def check_inputs_not_the_map(start, forward, reverse, metacool):
    """The files given to a pipeline run from reads, alignments or pairs
    (-1, -2): none may be its contact map (metacool: the same path, or a
    symbolic link to it or from it), which the run removes before it builds
    the map of its contigs there. And a pairs file may not be a contact
    map, which a run takes as its first argument."""
    if start == "metacool":
        return
    for path in [path for paths in (forward, reverse) if paths for path in paths.split(",")]:
        if mtmc.is_input_map(path, metacool):
            logger.error(
                "%s is the contact map of this run (%s), which the run removes before it builds the map of its "
                "contigs there: give it as the first argument to start from it, or write the run to another output "
                "directory or under another --prefix.",
                path,
                metacool,
            )
            raise ValueError
    if start == "pair" and forward:
        for path in forward.split(","):
            if mtmc.is_map(path):
                logger.error("%s is a contact map, not a pairs file: give it as the first argument to start from it.", path)
                raise ValueError


def contig_network(metacool, contig_data, outdir, tmp_dir, threads, *args, libraries=None):
    """network.txt and contig_data_network.txt of a run (args: normalization
    and self contacts), from the .metacool map metacool (see
    metator.network.metacool_to_contacts), and, from the maps of its
    libraries (libraries), the network of each library (network_{i}.txt)
    and its hits (hit_data_alignment.txt). The map is the one a run starts
    from (its first argument), or the map of its contigs that a pipeline
    run built from its pairs files (see network_map). A map that cannot
    give them stops the run, as metaTOR stops on an input it cannot use:
    the reason is logged as an error, in the log file of the run, and its
    MetacoolError is raised, as `metator metacool` does."""
    files = ("network.txt", "contig_data_network.txt")
    try:
        return mtn.metacool_to_contacts(metacool, contig_data, outdir, *files, tmp_dir, threads, *args, libraries=libraries)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise


def network_map(start, forward, alignment_files, metacool):
    """The .metacool map that a pipeline run reads its network from, and the
    maps of its libraries: the map that the run starts from (its first
    argument; its libraries summed: none), or metacool, the map of its
    contigs that a run from reads, alignments or pairs built from its pairs
    files (see pipeline_contig_map): the network of those pairs, byte for
    byte, as both count every pair. With several libraries, the network of
    each library (network_{i}.txt) and its hits (hit_data_alignment.txt)
    come from the map of that library, an entry of the collection next to
    metacool (see library_maps_path); with one, the run writes
    network_0.txt as a copy of network.txt, whose bytes the pairs path
    writes there too."""
    if start == "metacool":
        return forward, None
    if len(alignment_files) == 1:
        return metacool, None
    collection = library_maps_path(metacool)
    return metacool, [mtmc.library_entry(collection, i) for i in range(len(alignment_files))]


def library_maps_path(metacool):
    """The .mmetacool collection of the maps of the libraries of a run with
    several libraries, next to the map of their sum: <prefix>.mmetacool."""
    return f"{os.path.splitext(metacool)[0]}.mmetacool"


def aligner_mode_for(aligner_mode, enzyme):
    """The cutsite mode digests the reads: without an enzyme, align
    iteratively instead."""
    if aligner_mode == "cutsite" and not enzyme:
        logger.warning("cutsite mode required an enzyme. Iterative mode will be use instead.")
        return "iterative"
    return aligner_mode


def check_ascii_contigs(fasta, index=None, aligner=None):
    """Refuse an assembly whose contig names are not ASCII, before any
    alignment: the contact map of the contigs, which a pipeline run builds,
    needs them in ASCII (metator.metacool.check_ascii_names, whose message
    is logged), and could only be built once the reads are aligned. The
    names are those that the network step reads
    (metator.io.read_fasta_names). The reason, with what to do (rename the
    contigs), is logged as an error, then ValueError is raised. The
    assembly is named as it was given: the fasta, or the index of aligner
    whose fasta was retrieved (index), which is to be built again once its
    contigs are renamed."""
    try:
        mtmc.check_ascii_names(mio.read_fasta_names(fasta), [])
    except mtmc.MetacoolError as error:
        reason = str(error)
    else:
        return
    if index is None:
        assembly, rebuild = fasta, ""
    else:
        assembly, rebuild = f"the {aligner} index {index}", ", build the index again"
    logger.error(
        f"The contig names of {assembly} cannot be used: {reason} Rename these contigs to ASCII names, in the "
        f"assembly and in the alignments or pairs files made from it{rebuild}, then run pipeline again."
    )
    raise ValueError


def reference_files(assembly, aligner, start, tmp_dir, ascii_names=False):
    """Index and fasta of the assembly.

    The assembly is a fasta, whose index is built when the reads are
    aligned (fastq start), or the basename of a bowtie2 or bwa index, whose
    fasta is retrieved. With ascii_names (a pipeline run that builds the
    map of its contigs), its contig names are checked as soon as its fasta
    is known, before the index is built (see check_ascii_contigs).
    """
    index = mio.check_fasta_index(assembly, mode=aligner)
    if index is None:
        if mio.check_is_fasta(assembly):
            fasta = assembly
            if ascii_names:
                check_ascii_contigs(fasta)
            # Starting from bam or pairs files skips the index generation.
            if start == "fastq":
                index = mio.generate_fasta_index(fasta, aligner, tmp_dir)
        else:
            logger.error("Please give as assembly argument a %s index or a fasta.", aligner)
            raise ValueError
    else:
        fasta = mio.retrieve_fasta(index, aligner, tmp_dir)
        if ascii_names:
            check_ascii_contigs(fasta, index=assembly, aligner=aligner)
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


def run_metadata(columns_from_map=()):
    """The command line and the parameters of the running command, as stored
    in a .metacool file, and the columns that a run from a map read from it
    (columns_from_map, see metator.network.contig_data_from_map): where the
    restriction sites and the coverage that the map then stores came from,
    which --enzyme and --depth do not say. A run from a map (its first
    argument, map) took no --start: its start is none, not the default of
    the option."""
    params = dict(click.get_current_context().params)
    if params.get("map") is not None:
        params["start"] = None
    metadata = {
        "command": shlex.join(["metator", *sys.argv[1:]]),
        "parameters": {name: list(value) if isinstance(value, tuple) else value for name, value in params.items()},
    }
    if columns_from_map:
        metadata["columns_from_map"] = list(columns_from_map)
    return metadata


def assembly_basename(assembly):
    """Name of an assembly: its fasta file name without the extensions, or
    the basename of its index."""
    return re.sub(r"(\.(fa|fas|fasta|fna))?(\.gz)?$", "", os.path.basename(os.path.normpath(assembly)))


@contextmanager
def terminated_as_interrupted():
    """Stop a run on SIGTERM as on Ctrl-C, so that it cleans up after itself:
    the lock of its map, the map it was writing next to it, its temporary
    files. A job scheduler sends SIGTERM at its time limit, or when a job is
    cancelled; by default, Python dies on it without a word, as on SIGKILL.
    The exit status is that of a process killed by SIGTERM, 143. Only the
    main thread can set a signal handler: elsewhere, the run is as it is.
    A process that the run forks dies of SIGTERM, as it would without the
    handler: a multiprocessing pool ends its workers so. Its SIGTERM is set
    back to the default before it runs any Python (_block_sigterm_to_fork).
    Where the handler runs in code whose exceptions CPython ignores, the run
    dies of SIGTERM at once, as without the handler (_die_of_an_ignored_stop).
    """
    if threading.current_thread() is not threading.main_thread() or signal.getsignal(signal.SIGTERM) is None:
        # Off the main thread, or under a handler that a host embedding
        # Python set from C, which could not be restored.
        yield
        return
    previous_hook = sys.unraisablehook
    sys.unraisablehook = partial(_die_of_an_ignored_stop, previous_hook)
    previous = signal.signal(signal.SIGTERM, _stop_run)
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, previous)
        sys.unraisablehook = previous_hook


class _Terminated(SystemExit):
    """The SystemExit that the SIGTERM handler of a run raises (_stop_run)."""


def _stop_run(signum, frame):
    """The SIGTERM handler of a run (terminated_as_interrupted)."""
    raise _Terminated(128 + signal.SIGTERM)


def _die_of_an_ignored_stop(previous_hook, unraisable):
    """sys.unraisablehook under the SIGTERM handler of a run: the SystemExit of
    the handler, raised where CPython ignores what is raised, kills the run
    with SIGTERM at once, as without the handler; it let the run go on. That
    is the handler run in a finalizer, a weakref callback or a fork hook: a
    multiprocessing pool that is dropped ends its workers in its finalizer
    (metator.partition.get_hamming_distance), and a SIGTERM that comes as the
    run forks is handled in a fork hook (_unblock_sigterm_in_parent). The
    other exceptions go to the previous hook."""
    if isinstance(unraisable.exc_value, _Terminated):
        if threading.current_thread() is threading.main_thread():
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
        os.kill(os.getpid(), signal.SIGTERM)
    previous_hook(unraisable)


_forking = threading.local()


def _block_sigterm_to_fork():
    """Before a fork, under the SIGTERM handler of a run
    (terminated_as_interrupted): block SIGTERM in the thread that forks, so
    that the child is born with it blocked, and unblocks it only once it is
    set back to its default (_default_sigterm_in_child). A SIGTERM that
    reached the child before its after-fork code ran was lost otherwise:
    the handler of the run caught it, then CPython dropped it. A
    multiprocessing pool that is closed and dropped at once
    (metator.partition.get_hamming_distance) ends its workers with SIGTERM
    as soon as they are forked: it then waited forever for the worker that
    missed it, and the run with it."""
    _forking.ours = signal.getsignal(signal.SIGTERM) is _stop_run
    _forking.unblock = _forking.ours and signal.SIGTERM not in signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM})


def _unblock_sigterm_in_parent():
    """After a fork, in the parent: SIGTERM is unblocked, if it was not
    blocked before the fork. One that the run received meanwhile is handled
    then, in this hook, whose exceptions CPython ignores: it kills the run at
    once (_die_of_an_ignored_stop)."""
    if _forking.unblock:
        signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGTERM})


def _default_sigterm_in_child():
    """After a fork, in the child of a run: SIGTERM is set back to its
    default, then unblocked if it was not blocked before the fork; a SIGTERM
    that the child received meanwhile kills it then."""
    if _forking.ours:
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        if _forking.unblock:
            signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGTERM})


os.register_at_fork(
    before=_block_sigterm_to_fork,
    after_in_parent=_unblock_sigterm_in_parent,
    after_in_child=_default_sigterm_in_child,
)


def run_map_path(start, forward, outdir, prefix, output_map=None, force=False, no_metacool=False):
    """The map a pipeline run writes: output_map when given, else the file of
    the input map of a run from a map (updated in place: a map given as
    file::/, the root of the file, is the file), else the map of a run from
    reads, alignments or pairs (metacool_path). The directory of output_map
    is made, as the map will be written there.

    Raises
    ------
    MetacoolError
        If output_map is given with no_metacool (no map is written), names
        an entry of a collection or a directory, is the collection of the
        input map (an entry), or exists and is not the input map, without
        force; if the input map is not there; or if the input map is an
        entry of a .mmetacool collection, without output_map: updating it in
        place would rewrite the whole collection; if the directory of
        output_map cannot be made; or if the map of a run from reads,
        alignments or pairs would be a URI (:: in its output directory or
        prefix).
    """
    reason = None
    if start == "metacool" and not os.path.exists(mtmc.map_file(forward)):
        reason = f"Cannot read the contact map {forward}: no such file."
    elif output_map and no_metacool:
        reason = f"--no-metacool writes no map: drop --output-map {output_map}, or --no-metacool."
    elif output_map and "::" in str(output_map):
        reason = f"--output-map takes the path of a .metacool file, not an entry of a collection: {output_map}."
    elif output_map and os.path.isdir(output_map):
        reason = f"{output_map} is a directory: --output-map takes the path of a .metacool file."
    elif output_map and start == "metacool" and mtmc.is_entry(forward) and mtmc.is_input_map(forward, output_map):
        reason = (
            f"{output_map} is the collection of the input map {forward}: a run does not replace a collection with "
            "the map of its clusters; give --output-map another file."
        )
    elif output_map:
        in_place = start == "metacool" and mtmc.is_input_map(forward, output_map)
        if os.path.lexists(output_map) and not force and not in_place:
            reason = f"{output_map} exists: give -F to replace it with the map of this run, or another --output-map."
    if reason:
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    if output_map:
        make_map_directory(output_map, "--output-map")
        return output_map
    if start == "metacool":
        if mtmc.is_entry(forward):
            message = (
                f"{forward} is an entry of a .mmetacool collection, which a run does not update in place: write the "
                "map of its clusters to another file with --output-map."
            )
            logger.error(message)
            raise mtmc.MetacoolError(message)
        # The file itself: file::/ (the root of the file) names it too.
        return mtmc.map_file(forward)
    path = metacool_path(outdir, prefix)
    if "::" in path:
        # Refused now, not once the map of the contigs is built (after the
        # alignment of a run from reads).
        reason = f"{path}, the map of this run, is not the path of a file: give -o/--outdir and --prefix without ::."
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    return path


def make_map_directory(path, option):
    """Make the directory of the map at path, where the run will write it
    (given with option): a path that is no file's (a URI, file::group, see
    metator.metacool._check_paths), before any directory is made, and a
    directory that cannot be made are refused, the reason logged."""
    if "::" in str(path):
        reason = f"{path} is not the path of a file: a map is written to a file of its own, not to a group inside one."
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    directory = os.path.dirname(os.path.abspath(path))
    try:
        os.makedirs(directory, exist_ok=True)
    except OSError as error:
        reason = f"The directory of {option} {path} cannot be made ({error.strerror}): write the map elsewhere."
        logger.error(reason)
        raise mtmc.MetacoolError(reason) from error


def metacool_path(outdir, prefix):
    """The contact map of a pipeline run: <outdir>/<prefix>.metacool, or
    <outdir>/metator.metacool without --prefix."""
    return join(outdir, f"{prefix or 'metator'}.metacool")


def _failure(error):
    """What the failure of a map says: the message of a MetacoolError, or,
    for another error, which is unexpected, its type and message. Returns
    it without its final period, and whether the error is unexpected (its
    traceback then goes to the log too)."""
    unexpected = not isinstance(error, mtmc.MetacoolError)
    reason = f"{type(error).__name__}: {error}" if unexpected else str(error)
    return reason.rstrip("."), unexpected


# What metacooler needs of the pairs files given to a run from pairs, which
# the run says when it stops on pairs files that metacooler cannot read.
PAIRS_FILES_NEEDS = (
    "metacooler, which builds it, needs the same number of fields on every row of each pairs file given to -1 (its "
    "header lines, which start with #, may come anywhere, as in pairs files of one assembly concatenated into one "
    "library), plain or compressed with gzip, bzip2 or, without --scaffold, zstd or lz4. Give the pairs files of "
    "several libraries comma-separated, one file per library."
)


def _pairs_unreadable(error):
    """Whether the map of the contigs failed on the content of its pairs
    files: metacooler reads them with pyarrow, whose errors (a parse error,
    on a header line amid the pairs or a row with one more field, for
    example) write_contig_map and write_library_maps turn into a
    MetacoolError, the latter through another, which names the library: the
    pyarrow error is somewhere in the chain of its causes. pyarrow is no
    dependency of metaTOR: it is looked up once metacooler has loaded it."""
    pyarrow = sys.modules.get("pyarrow")
    if pyarrow is None or not isinstance(error, mtmc.MetacoolError):
        return False
    cause = error.__cause__
    while cause is not None:
        if isinstance(cause, pyarrow.ArrowException):
            return True
        cause = cause.__cause__
    return False


def pipeline_contig_map(pairs_files, contig_data, outfile, assembly, threads, tmp_dir, given_pairs=False):
    """Build the contact map of the contigs of a pipeline run from its pairs
    files, as soon as it has them, at outfile (see metacool_path): every
    contig in no bin yet (see metator.metacool.write_contig_map). At the end
    of the run, pipeline_metacool writes the map of the bins over it or,
    with --no-metacool, the run removes it. given_pairs: the pairs files
    are those given to a run from pairs (--start pair), not those that the
    run made from reads or alignments; when metacooler cannot read them,
    the error also says what it needs of them (PAIRS_FILES_NEEDS).

    A map of the same name left by an earlier run is removed first: this
    run replaces the tables it was made from, so it must not stay behind,
    passing for this run's, when this one cannot be built. A map that
    cannot be built stops the run, --no-metacool or not, whatever the number
    of libraries: the reason is logged as an error, in the log file of the
    run (with the type and the traceback of an unexpected error, which is
    not a MetacoolError), then the error is raised as it is. Nothing is left
    at outfile: metacooler removes its partial file, and write_contig_map a
    map whose columns it cannot store. Once built, the map stays at outfile,
    in no bin, until the map of the bins replaces it: when that map cannot
    be written, and when the run stops with an error before its end
    (--no-metacool or not). A symbolic link at outfile is removed, not its
    target. The inputs of the run are not at outfile:
    check_inputs_not_the_map refuses them first.
    """
    collection = library_maps_path(outfile)
    try:
        for path in (outfile, collection):
            if os.path.lexists(path):
                logger.info(f"Removing {path}, the contact map of an earlier run.")
                os.remove(path)
        if len(pairs_files) > 1:
            mtmc.write_library_maps(
                pairs_files,
                contig_data,
                outfile,
                collection,
                assembly_name=assembly_basename(assembly),
                threads=threads,
                force=True,
                metadata=run_metadata(),
                tmp_dir=tmp_dir,
            )
        else:
            mtmc.write_contig_map(
                pairs_files,
                contig_data,
                outfile,
                assembly_name=assembly_basename(assembly),
                threads=threads,
                force=True,
                metadata=run_metadata(),
                tmp_dir=tmp_dir,
            )
    except Exception as error:
        reason, unexpected = _failure(error)
        needs = f" {PAIRS_FILES_NEEDS}" if given_pairs and _pairs_unreadable(error) else ""
        logger.error(f"The contact map of the contigs cannot be built at {outfile}: {reason}.{needs}", exc_info=unexpected)
        raise


def pipeline_metacool(metacool, outdir, outfile, contig_map=False, columns_from_map=()):
    """Write the contact map of the clusters of a pipeline run at outfile
    (see run_map_path), from a .metacool map of its contigs and the final
    tables of the run (see metator.metacool.reassign_metacool), with its
    clustering rounds (see metator.clusters): the map that the run built
    from its pairs files (contig_map, at outfile), or the input map of a
    run that started from one, which outfile is by default, whose columns
    it read (columns_from_map: stored in the metadata of the new map).

    outfile is only replaced once the new map has been written and
    validated: the map of the contigs that the run built stays there, in no
    cluster, and so does an input map at outfile (the same path, or a
    symbolic link to it or from it), when the new map cannot be written.
    Another map of the same name, left by an earlier run, is removed first:
    this run has replaced the tables it was made from, so it must not stay
    behind, passing for this run's, when the new map cannot be written. A
    symbolic link at the outfile is replaced, not its target.

    Raises
    ------
    MetacoolError
        If the map cannot be written: the reason, and what is kept, is
        logged as an error first. The other outputs of the run are kept.
    """
    # What the error adds when the new map is not written.
    kept = ""
    try:
        if contig_map:
            kept = "It is still the contact map of the contigs of the run, in no cluster. "
        elif mtmc.is_input_map(metacool, outfile):
            kept = "It is still the input map of the run, as it was. "
            logger.info(f"{outfile} is the input map {metacool}: it is replaced once the map of the clusters is validated.")
        elif os.path.lexists(outfile):
            logger.info(f"Removing {outfile}, the contact map of an earlier run.")
            os.remove(outfile)
        mtmc.reassign_metacool(
            metacool,
            join(outdir, "contig_data_final.txt"),
            outfile,
            bin_summary=join(outdir, "bin_summary.txt"),
            force=True,
            metadata=run_metadata(columns_from_map),
            rounds=mtcl.read_rounds(outdir),
        )
    except Exception as error:
        reason, unexpected = _failure(error)
        message = f"{outfile} is not written: {reason}. {kept}The other outputs of the run are kept."
        logger.error(message, exc_info=unexpected)
        raise mtmc.MetacoolError(message) from error


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


class StartChoice(click.Choice):
    """The values of --start: what -1/--forward holds. A run from a map takes
    the map as its first argument instead: --start metacool, which named it
    before metaTOR 2.0, is refused with what replaced it."""

    def convert(self, value, param, ctx):
        if value == "metacool":
            command = ctx.info_name if ctx is not None else "pipeline"
            self.fail(
                f"--start metacool is gone in metaTOR 2.0: give the map as the first argument, metator {command} "
                "MAP -a ASSEMBLY ...",
                param,
                ctx,
            )
        return super().convert(value, param, ctx)


# Options of metaTOR 1 that 2.0 removed, by name, with what replaced them: a
# command line that still gives one is refused with it (see Command).
REMOVED_OPTIONS = dict.fromkeys(
    ["-E", "--edge"],
    "--edge is gone in metaTOR 2.0: the network counts every pair, as its contact map does, which holds no read positions.",
)


class Command(click.RichCommand):
    """A command that refuses an option that metaTOR 2.0 removed with what
    replaced it (REMOVED_OPTIONS), rather than as an option it never had.
    In a command that takes a map as its first argument (network,
    pipeline), a first argument that is no map, given with -1/--forward,
    is one argument too many (the reverse reads given without -2, for
    example): refused as click refused it before the map was an argument.
    """

    def parse_args(self, ctx, args):
        try:
            rest = super().parse_args(ctx, args)
        except click.NoSuchOption as error:
            if error.option_name in REMOVED_OPTIONS:
                raise click.UsageError(REMOVED_OPTIONS[error.option_name], ctx) from error
            raise
        extra = ctx.params.get("map")
        if extra is not None and ctx.params.get("forward") and not mtmc.is_map(extra):
            raise click.UsageError(f"Got unexpected extra argument ({extra})", ctx)
        return rest


def map_argument():
    """The .metacool map of the contigs that a network or pipeline run starts
    from: the first argument, instead of the reads, alignments or pairs of
    -1/--forward."""
    return click.argument(
        # Shown as MAP: the parameter of network and pipeline is named so.
        "map",
        required=False,
        # Shown as FILE; a directory is refused as a usage error. An entry of a
        # collection (file::group) is no path: run_input checks it.
        type=click.Path(dir_okay=False),
        panel="Input",
        help="A .metacool contact map of the contigs of the assembly, or an entry of a .mmetacool collection "
        "(file.mmetacool::entries/NAME): the run starts from it, instead of -1/--forward.",
    )


def reads_options():
    """Inputs of the network step: reads, alignments or pairs."""
    return _options(
        _option(
            "-1",
            "--forward",
            metavar="STR",
            panel="Input",
            help="Fastq file or list of Fastq separated by a comma containing the forward reads to be aligned or their "
            "corresponding bam or pairs files (see --start). A run from a map takes the map as its first argument "
            "instead.",
        ),
        _option(
            "-2",
            "--reverse",
            metavar="STR",
            panel="Input",
            help="Fastq file or list of Fastq separated by a comma containing the reverse reads to be aligned or their "
            "corresponding bam files. Forward and reverse reads need to have the same identifier (read names). Not "
            "needed from pairs files, or the bam files of bwa, which hold both mates.",
        ),
        assembly_option(ASSEMBLY_INDEX_HELP),
        _option(
            "-S",
            "--start",
            type=StartChoice(STARTS),
            default="fastq",
            show_default=True,
            panel="Input",
            help="What -1/--forward holds: the reads ([cyan]fastq[/]), their alignments ([cyan]bam[/]) or their "
            "[cyan]pair[/]s. Not with a map.",
        ),
        _option(
            "-d",
            "--depth",
            metavar="FILE",
            panel="Input",
            help="The depth.txt file from the shotgun reads used to made the assembly computed by "
            "[cyan]jgi_summarize_bam_contig_depths[/] from the metabat2 pipeline. With --start metacool, the "
            "Shotgun_coverage column of the map gives the coverage without it, and --depth replaces it.",
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
            help="The list of restriction enzyme used to digest the contigs separated by a comma. Example: HpaII,MluCI. "
            "With --start metacool, the Restriction_site column of the map gives the restriction sites without it, "
            "and --enzyme replaces it.",
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
            "[cyan]theoritical_hit[/] need --depth. With --start metacool, the Restriction_site and "
            "Shotgun_coverage columns of the map stand for them (see --enzyme and --depth).",
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


@cli.command(
    cls=Command, short_help="Generate the metaHiC contig network from fastq, bam, pairs or a .metacool map and normalize it."
)
@map_argument()
@reads_options()
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
    map,  # the argument MAP
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
    self_contacts,
    outdir,
    tmpdir,
    threads,
    no_clean_up,
):
    """Generate the metaHiC contig network.

    [cyan]metator network MAP -a ASSEMBLY[/] builds the network from a
    .metacool contact map of the contigs, whose libraries are summed (see
    the README): that of the pairs it was made of.

    [cyan]metator network -1 FORWARD -2 REVERSE -a ASSEMBLY[/] builds it
    from fastq, bam ([cyan]-S bam[/]) or pairs ([cyan]-S pair[/], -1 only)
    files. Forward and reverse reads or alignments need to be in separated
    files with the same read names. If multiple fastq are given, they are
    aligned separately and the hit file has a number of hits for each
    sample.

    From a map, its Restriction_site and Shotgun_coverage columns, when
    they give every contig a value, give the restriction sites and the
    coverage of the contigs that --enzyme and --depth would measure; given,
    these options replace them.

    In the output network, contigs are the network nodes and the edges are
    the contact counts: every pair counts, as in a contact map.

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

        # The input: a map, or reads, alignments or pairs, with the reverse
        # reads or alignments that they need.
        start, forward = run_input(map, forward, reverse, start)
        check_reverse_reads(start, aligner, reverse)
        aligner_mode = aligner_mode_for(aligner_mode, enzyme)
        check_normalization(normalization, enzyme, depth, forward if start == "metacool" else None)

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
            if start == "metacool":
                contig_data_from_map(forward, contig_data, normalization, enzyme, depth)
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

        # Build the network: from the pairs files, or from the input map.
        if start == "metacool":
            contig_network(forward, contig_data, outdir, tmp_dir, threads, normalization, self_contacts)
        else:
            # Every pair counts, as in a contact map: --edge is gone.
            mtn.alignment_to_contacts(
                alignment_files, contig_data, 0, hit_data, outdir, "network.txt", "contig_data_network.txt",
                tmp_dir, threads, normalization, self_contacts,
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
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up):
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

    Only bins with a weighted completeness of at least 33% and a
    contamination, (weighted redundancy - 1) / weighted completeness, of at
    least 5% are subject to the recursive step. Its bins replace the
    original bin when the most complete of them keeps enough of its
    completeness; otherwise the original bin is kept.
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


@cli.command(cls=Command, short_help="Run network, partition and validation (and scaffold the clusters) in one command line.")
@map_argument()
@reads_options()
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
@force_option("Replace the directories of clusters of an earlier run in the output directory, and an existing --output-map.")
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated for the alignment.")
@no_clean_up_option()
@_option(
    "--no-metacool",
    is_flag=True,
    help="Do not keep the contact map of the contigs, grouped by cluster: a run from reads, alignments or pairs "
    "still builds the map of its contigs from the pairs files (a map that cannot be built stops the run), reads its "
    "networks from it, and removes it at the end, with the maps of its libraries (<prefix>.mmetacool) when it has "
    "several; a run from a map leaves that map as it was.",
)
@_option(
    "--output-map",
    metavar="FILE",
    help="Write the contact map of the clusters to FILE (replaced only with -F). By default, a run from a .metacool "
    "map updates that map in place, and any other run writes <outdir>/<prefix>.metacool (metator.metacool without "
    "--prefix).",
)
def pipeline(
    map,  # the argument MAP
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
    output_map,
):
    """Launch the full metaTOR pipeline.

    Cluster the contigs of the assembly from their metaHiC contacts.

    [cyan]metator pipeline MAP -a ASSEMBLY[/] starts from a .metacool contact
    map of the contigs: the run reads its network from that map, and the
    restriction sites and the coverage of the contigs from its
    Restriction_site and Shotgun_coverage columns, unless --enzyme and
    --depth are given; it updates the map in place with the map of its
    clusters (or writes --output-map), and cannot scaffold.

    [cyan]metator pipeline -1 FORWARD -2 REVERSE -a ASSEMBLY[/] starts from
    the fastq, bam ([cyan]-S bam[/]) or pairs ([cyan]-S pair[/], -1 only)
    files. As soon as it has the pairs files, the run builds the contact map
    of the contigs, in no cluster yet, as a .metacool file (see
    [cyan]metator metacool[/] and the README): a map that cannot be built
    stops the run. The run reads its network from it; with several
    libraries, the map of each library (<prefix>.mmetacool) gives the
    network of each. At the end, the map of the clusters (the contigs
    grouped by cluster) is written over it or, with --no-metacool, it is
    removed.

    Whatever the input, the network counts every pair, as a contact map
    does. A validation step decontaminates the clusters when it's
    necessary; it is the critical step for memory usage (~40G).
    """
    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with terminated_as_interrupted(), private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir, ExitStack() as held:
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
        start, forward = run_input(map, forward, reverse, start)
        aligner_mode = aligner_mode_for(aligner_mode, enzyme)
        # The normalization of a run from a map is checked against its map,
        # once its path is resolved (run_map_path, below).
        if start != "metacool":
            check_normalization(normalization, enzyme, depth)
        check_metacool_start(start, scaffold)
        metacool_file = run_map_path(start, forward, outdir, prefix, output_map, force, no_metacool)
        if start == "metacool":
            check_normalization(normalization, enzyme, depth, forward)
        check_inputs_not_the_map(start, forward, reverse, metacool_file)
        check_pairs_compression(start, forward, scaffold)
        check_reverse_reads(start, aligner, reverse)
        bins_prefix = bin_prefix(prefix)
        # A run that writes its map holds the map's lock from now on: a second
        # run on the same map is refused before any work.
        if not (start == "metacool" and no_metacool):
            held.enter_context(mtmc.map_lock(metacool_file))

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
        # bowtie2 index is given, retreive the fasta. Its contig names are
        # checked first, but for a run from a map: the map of the contigs
        # needs them in ASCII.
        index, fasta = reference_files(assembly, aligner, start, tmp_dir, ascii_names=start != "metacool")
        # The contig data of a run from pairs or from a map, before the bins
        # of an earlier run are replaced: a map whose columns cannot give
        # what the normalization needs is refused once the assembly is read.
        if start_step > 2:
            # A map has no pairs files: its libraries are summed.
            alignment_files = forward.split(",") if start == "pair" else []
            contig_data, _ = mtn.create_contig_data(fasta, len(alignment_files) or 1, depth, enzyme)
            columns_from_map = (
                contig_data_from_map(forward, contig_data, normalization, enzyme, depth) if start == "metacool" else []
            )
        for path in bin_dirs:
            make_output_dir(path, force)

        # Run the whole workflow
        if start_step <= 2:
            # Align pair-end reads with bowtie2
            # The hits of each library come from its map (see network_map).
            alignment_files, contig_data, _ = mta.get_contact_pairs(
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
            columns_from_map = []
        # The contact map of the contigs, built once from the pairs files as
        # soon as the run has them (with several libraries, the map of each
        # library too, and their sum): the map of the clusters is written over
        # it at the end. A map that cannot be built stops the run. A run from
        # a map writes the map of its clusters from that map.
        if start != "metacool":
            pipeline_contig_map(
                alignment_files, contig_data, metacool_file, assembly, threads, tmp_dir, given_pairs=start == "pair"
            )
        # Build the network: from the input map or the map of the contigs,
        # and the network of each library from its map (see network_map).
        network_source, libraries = network_map(start, forward, alignment_files, metacool_file)
        network_file, contigs_data_file = contig_network(
            network_source, contig_data, outdir, tmp_dir, threads, normalization, False, libraries=libraries
        )
        # network_0.txt, the network of the one library of a run that read
        # its network from the map of its contigs: a copy of network.txt,
        # whose bytes the pairs path writes there too. A run from a map
        # writes none: its libraries are summed.
        if start != "metacool" and len(alignment_files) == 1:
            shutil.copyfile(network_file, join(outdir, "network_0.txt"))

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
            # The heatmap of the clusters is drawn from the map of the run:
            # the map it started from, or the map of its contigs.
            metacool=forward if start == "metacool" else metacool_file,
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
        mio.remove_pyfastx_index(fasta)

        # The map of the bins, written last: over the map of the contigs that
        # the run built, or from the input map of a run that started from one.
        if start == "metacool":
            if not no_metacool:
                pipeline_metacool(forward, outdir, metacool_file, columns_from_map=columns_from_map)
        elif no_metacool:
            logger.info(f"Removing {metacool_file}, the contact map of the contigs (--no-metacool).")
            os.remove(metacool_file)
            if len(alignment_files) > 1:
                collection = library_maps_path(metacool_file)
                logger.info(f"Removing {collection}, the contact maps of the libraries (--no-metacool).")
                os.remove(collection)
        else:
            pipeline_metacool(metacool_file, outdir, metacool_file, contig_map=True)

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
    help="Name of the metator object, or its numerical ID for the core bin and overlapping bin objects. A recursive "
    "bin is numbered within its overlapping bin: give its name. "
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
        mio.remove_pyfastx_index(assembly)

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# metacool
# ---------------------------------------------------------------------------


@cli.command(short_help="Write the contact map of the contigs, grouped by bin, as a .metacool file.")
@click.argument(
    "pairsfile",
    nargs=-1,
    required=True,
    help="Pairs file(s) of the run, or one .metacool map of its contigs (or an entry of a .mmetacool collection: "
    "file.mmetacool::entries/NAME), the map that a run from a map started from.",
    panel="Input",
)
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
@_option(
    "-n",
    "--assembly-name",
    metavar="STR",
    help="Name of the assembly, stored in the file. A map written from a .metacool map keeps its name.",
)
@threads_option(help="Number of threads binning the pairs.")
@force_option("Overwrite the output file if it exists.")
def metacool(pairsfile, contig_data, bin_summary, outfile, assembly_name, threads, force):
    """Write the contact map of a metaTOR run as a .metacool file.

    A .metacool file is the contact map of metacooler: a cooler whose
    chromosomes are the contigs, one bin per contig, grouped by genome: one
    genome per final bin of the contig data file, and [cyan]__unbinned__[/]
    for the other contigs, with all their contacts. The README says what
    else it holds. The pairs files of several libraries are summed.

    It is written from the pairs files of the run or, for a run that
    started from a .metacool map of the contigs (run step by step:
    [cyan]metator network MAP[/], partition, validation), from that map:
    the map that [cyan]pipeline MAP[/] writes at the end of such a run.

    The contig and bin names must be ASCII, and no bin may be named as
    metacooler marks a contig in no bin (NA, NaN, N/A, N/D, -1 or
    __unbinned__ in any case, and nd in any case other than ND, which marks
    a contig in no bin). [cyan]pipeline[/] writes this file at the end of a
    run: from the map of the contigs that it builds from the pairs files,
    or from the map that it starts from.
    """
    # Generate log
    log_file = log_file_name("", "metacool")
    generate_log_header(log_file)

    # Its map is locked as a pipeline run locks the map it writes: the two
    # cannot write one map at once (map_lock logs its refusal). The map's
    # directory is made first, as the writers make it.
    make_map_directory(outfile, "-o/--outfile")
    with terminated_as_interrupted(), mtmc.map_lock(outfile, elsewhere="another -o/--outfile"):
        _write_metacool_command(pairsfile, contig_data, bin_summary, outfile, assembly_name, threads, force)

    generate_log_footer(log_file)


def _write_metacool_command(pairsfile, contig_data, bin_summary, outfile, assembly_name, threads, force):
    """The map of metator metacool, written once its lock is held."""
    try:
        maps = [path for path in pairsfile if mtmc.is_map(path)]
        if maps:
            # The map of the bins written from the map of the contigs, as
            # a pipeline run from a map writes it.
            if len(pairsfile) > 1:
                raise mtmc.MetacoolError(
                    f"metator metacool writes the map from pairs files, or from one .metacool map: got "
                    f"{len(pairsfile)} files, {len(maps)} of them .metacool maps."
                )
            if assembly_name:
                raise mtmc.MetacoolError(
                    f"--assembly-name names the assembly of a map written from pairs files: the map written from "
                    f"{maps[0]} keeps the assembly name of that map."
                )
            mtmc.reassign_metacool(maps[0], contig_data, outfile, bin_summary=bin_summary, force=force, metadata=run_metadata())
        else:
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

    The MGE fasta will contain one entry by MGE MAG, with 200 "N" spacers
    between contigs.
    """
    # Generate log
    log_file = log_file_name("", "mge")
    generate_log_header(log_file)

    # Private temporary directory inside --tmpdir.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up):
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
            threshold_bin=threshold_bin,
            random=random,
        )

    if not no_clean_up:
        # Delete pyfastx index:
        mio.remove_pyfastx_index(assembly)

    generate_log_footer(log_file)
