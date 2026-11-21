#!/usr/bin/env python3
# coding: utf-8

"""metaTOR command line interface.

One function per command: network, partition, validation, pipeline and
mge. The ``metator`` entry point is :func:`cli`.
"""

import logging
import os
import re
import shlex
import shutil
import stat
import sys
import time
from contextlib import ExitStack, contextmanager
from os.path import exists, join

import rich_click as click

import metator.clusters as mtcl
import metator.io as mio
import metator.log as mtl
import metator.metacool as mtmc
import metator.mge as mtm
import metator.network as mtn
import metator.partition as mtp
import metator.validation as mtv
from metator.log import logger
from metator.version import __version__

ALGORITHMS = ["louvain", "leiden"]
NORMALIZATIONS = ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"]

# Panel holding the options common to most commands (output, temporary
# files, threads), with --help.
GENERAL = "Options"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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
    fails: its files are left for inspection. An empty --tmpdir (-T '') is
    the default, ./tmp.
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
    it: a symbolic link too, even one to nothing."""
    if os.path.lexists(path) and not force:
        logger.error(
            "%s already existed. Remove directory or use -F argument to overwrite it.",
            path,
        )
        raise ValueError


def make_output_dir(path, force):
    """Create an output directory of bins; an existing one is replaced with
    --force, or refused. A symbolic link is removed, not its target.
    """
    refuse_output_dir(path, force)
    if os.path.islink(path):
        os.remove(path)
    elif exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def check_normalization(normalization, enzyme, depth, metacool, read=None):
    """The normalizations that need the restriction sites or the shotgun
    coverage of the contigs: measured on the assembly with --enzyme and
    --depth or read from the Restriction_site and Shotgun_coverage columns
    of the contigs table of the map of the run, metacool (see
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
        if read is None:
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


def run_input(input_map):
    """Check the .metacool map that a network or pipeline run starts from,
    its first argument (input_map), before any work. Refused, the reason
    logged as an error: a map that is not there or is no file (see
    refuse_missing_map), a file that is no map, a map named as a URI (see
    refuse_a_map_named_as_a_uri), and a map that metacooler cannot open,
    with its reason (a .mmetacool collection rather than one of its entries,
    for example), as MetacoolError. A pipeline run checks its map before it
    replaces the clusters of an earlier run (-F)."""
    path = refuse_missing_map(input_map)
    if not mtmc.is_map(input_map):
        logger.error(
            f"{input_map} is not a .metacool contact map: make the map of the contigs with hicstuff pipeline -M "
            "metacool, from the reads of the run, or metacooler create pairs, from its pairs files."
        )
        raise ValueError
    refuse_a_map_named_as_a_uri(input_map, path)
    open_input_map(input_map)


def refuse_missing_map(input_map):
    """Refuse a map whose file is not there, or is no file (a directory given
    as dir::/), before any work: the reason is logged as an error, then
    MetacoolError is raised. Returns the path of the file of the map:
    input_map itself when it is a file (a file whose name holds :: is that
    file, not the URI of a map inside one: see metator.metacool.is_map),
    else its file part."""
    path = input_map if os.path.isfile(input_map) else mtmc.map_file(input_map)
    if not os.path.isfile(path):
        reason = "no such file" if not os.path.exists(path) else f"{path} is not a file"
        message = f"Cannot read the contact map {input_map}: {reason}."
        logger.error(message)
        raise mtmc.MetacoolError(message)
    return path


def refuse_a_map_named_as_a_uri(input_map, path):
    """Refuse, before any work, a map whose file name holds :: (path, its
    file, is input_map itself), which metacooler reads as file::group: the
    reason is logged as an error, then MetacoolError is raised."""
    if path == input_map and "::" in input_map:
        message = (
            f"Cannot read the contact map {input_map}: metacooler reads a path that holds :: as file::group. "
            "Rename the map without ::."
        )
        logger.error(message)
        raise mtmc.MetacoolError(message)


def open_input_map(input_map):
    """The map, opened (a metacooler.MetaCooler). A map that metacooler
    cannot open is refused, before any work: its reason (a .mmetacool
    collection rather than one of its entries, which it names, for example)
    is logged as an error, then MetacoolError is raised."""
    try:
        return mtmc.open_map(input_map)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise


def refuse_an_input_as_the_map(output_map, inputs, assembly=None):
    """Refuse, before any work, an --output-map that is one of the files of
    the inputs of the run (the assembly, a table): the run would replace it
    with its map; or the pyfastx index of the assembly (<assembly>.fxi),
    there or not: the run removes it at its end (see
    metator.io.remove_pyfastx_index), with the map. The reason is logged as
    an error, then MetacoolError is raised."""
    if output_map is None:
        return
    if assembly is not None and os.path.realpath(output_map) == os.path.realpath(f"{assembly}.fxi"):
        reason = (
            f"--output-map {output_map} is the pyfastx index of the assembly {assembly}, which the run removes at "
            "its end: write the map elsewhere."
        )
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    if not os.path.exists(output_map):
        return
    for path in inputs:
        if path is not None and os.path.isfile(path) and os.path.samefile(path, output_map):
            reason = f"--output-map {output_map} is {path}, an input of the run: write the map elsewhere."
            logger.error(reason)
            raise mtmc.MetacoolError(reason)


# What a run that updates its map in place replaces of what the map holds
# already (see warn_of_what_is_replaced): its clusters (whether it does), and
# its clustering rounds made by these steps.
REPLACED_IN_PLACE = {
    "pipeline": (True, ("partition", "validation")),
    "partition": (False, ("partition", "validation")),
    "validation": (True, ("validation",)),
    "mge": (False, ("mge",)),
}


def warn_of_what_is_replaced(command, metacool):
    """Warn, before any work, that a run (command) that updates the map
    metacool in place replaces what the map holds already: its clusters (the
    genomes other than __unbinned__, of an earlier run, or of another
    binner), and its clustering rounds made by an earlier run (see
    REPLACED_IN_PLACE). Neither can be found again in the updated map, whose
    history says only that its contigs were reassigned. A run from a map
    replaced them without -F, and without a warning. --output-map writes
    the map of the run to another file instead, which leaves metacool as it
    is."""
    clusters, steps = REPLACED_IN_PLACE[command]
    try:
        with mtmc.read_errors(metacool):
            mc = mtmc.open_map(metacool)
            n_clusters = int((mc.genomes["name"].astype(str) != "__unbinned__").sum()) if clusters else 0
            rounds = [str(entry["name"]) for entry in mc.cluster_rounds if entry.get("step") in steps]
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    held = []
    if n_clusters:
        held.append(f"{n_clusters} cluster(s)")
    if rounds:
        held.append(f"the clustering round(s) {', '.join(rounds)}")
    if held:
        logger.warning(
            f"{metacool} holds {' and '.join(held)}, which this run replaces, in place: give --output-map to write the "
            f"map of this run to another file, and leave {metacool} as it is."
        )


def step_map(command, input_map, output_map, force, held, contigs=None, inputs=(), option=None, assembly=None):
    """The map that a partition, validation or mge run (command) writes, or
    None for a run without one. A run step by step from a map (network MAP,
    then partition MAP and validation MAP) updates the map, its first
    argument (input_map), in place, or writes output_map (--output-map)
    instead (see run_map_path), as pipeline MAP does. mge updates the map
    given with an option (option: --map) so.

    The map is checked as run_input checks the map of a network or pipeline
    run, then the run holds the lock of the map that it writes until its
    end, and opens its input map under it (in held, an ExitStack: see
    hold_map). Refused before any work, the reason logged as an error, as
    MetacoolError: output_map without input_map; an input map that is not
    there or is no file, is no map, is named as a URI, or that metacooler
    cannot open (a .mmetacool collection rather than one of its entries,
    for example); a map that is no map of the contigs of contigs, the
    contig data file of the run (other contigs or lengths: see
    metator.metacool.check_map_of_contigs); an output_map that is one of the
    files of contigs and inputs, or the pyfastx index of the assembly (see
    refuse_an_input_as_the_map); what run_map_path refuses (an entry of a
    collection without output_map, an output_map that is an entry, a
    directory, the collection of the input entry, a hard link of the input
    map, or that exists and is not the input map, without force); a map
    that another run holds; and, in place, a read-only map. A run that
    updates the map in place warns of the clusters and the rounds that it
    replaces there (see warn_of_what_is_replaced). Returns the map to write
    (None for a run without one) and the input map as the run read it once
    it held the lock (a metacooler.MetaCooler, or None), which its map is
    refused to be written from if it changed since (see
    metator.metacool.reassign_metacool)."""
    if option is None:
        updates, given, map_of = (
            f"{command} MAP",
            "as the first argument",
            f"the map of the contigs of the run as the first argument of {command}",
        )
    else:
        updates, given, map_of = f"{command} {option} MAP", f"with {option}", f"the map of the run with {option}"
    if input_map is None:
        if output_map is not None:
            reason = f"--output-map writes the map that {updates} updates: give the map {given}, or drop --output-map."
            logger.error(reason)
            raise mtmc.MetacoolError(reason)
        return None, None
    path = refuse_missing_map(input_map)
    if not mtmc.is_map(input_map):
        reason = f"{input_map} is not a .metacool contact map: give {map_of}, or leave it out."
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    refuse_a_map_named_as_a_uri(input_map, path)
    open_input_map(input_map)
    if contigs is not None:
        try:
            mtmc.check_map_of_contigs(input_map, contigs)
        except mtmc.MetacoolError as error:
            logger.error(str(error))
            raise
    refuse_an_input_as_the_map(output_map, [contigs, *inputs], assembly)
    outfile, in_place = run_map_path(input_map, output_map, force)
    return outfile, hold_map(held, command, input_map, outfile, in_place)


def _is_in(path, directory):
    """Whether path is directory, or a path inside it."""
    return path == directory or path.startswith(directory.rstrip(os.sep) + os.sep)


def refuse_an_input_in_the_outputs(path, outputs, force, option):
    """Refuse an input (validation -f) that is, or is in, an output that the
    run replaces with --force (see start_step): the run would remove it
    before reading it, or replace the output it is in, through a symbolic
    link too. Without --force, those outputs are refused anyway."""
    if not force:
        return
    for output in outputs:
        if not os.path.lexists(output):
            continue
        inside = _is_in(os.path.abspath(path), os.path.abspath(output))
        inside = inside or _is_in(os.path.realpath(path), os.path.realpath(output))
        if inside:
            logger.error(
                "%s, given as %s, is in %s, which -F removes: write the run to another directory (-o).", path, option, output
            )
            raise ValueError


def start_step(
    command,
    outdir,
    cluster_dirs,
    input_map,
    output_map,
    force,
    held,
    prefix=None,
    contigs=None,
    inputs=(),
    assembly=None,
):
    """Start a partition or validation run (command) in outdir: its log file,
    the map that it writes (see step_map) and its directories of clusters
    (cluster_dirs, see make_output_dir), refused without force when an
    earlier run left them. A run with a map (input_map, or output_map, which
    needs it) checks it, in its log file, before it replaces the
    directories of an earlier run (-F): against contigs, the contig data
    file of the run, whose contigs it must be a map of, and against
    inputs, the other input files, which output_map must not be (see
    step_map). A run without one replaces them first, then starts its log
    file, as before a run took a map. Either checks the prefix of the
    clusters (--prefix, see bin_prefix) before it replaces them. Returns the
    log file, the map to write (None without one), the prefix, and the
    input map as the run read it at its start (None without one: see
    step_map)."""
    with_a_map = input_map is not None or output_map is not None

    def replace_outputs():
        for path in cluster_dirs:
            make_output_dir(path, force)

    if not with_a_map:
        # The prefix of the clusters first: a refused run keeps them.
        bins_prefix = bin_prefix(prefix)
        replace_outputs()
    log_file = log_file_name(outdir, command)
    generate_log_header(log_file)
    metacool_file, planned = step_map(command, input_map, output_map, force, held, contigs, inputs, assembly=assembly)
    if with_a_map:
        bins_prefix = bin_prefix(prefix)
        replace_outputs()
    return log_file, metacool_file, bins_prefix, planned


def contig_network(metacool, contig_data, outdir, tmp_dir, threads, normalization):
    """network.txt and contig_data_network.txt of a run, normalized
    (normalization), from the .metacool map metacool, its first argument
    (see metator.network.metacool_to_contacts). A map that cannot give them
    stops the run, as metaTOR stops on an input it cannot use: the reason is
    logged as an error, in the log file of the run, and its MetacoolError is
    raised."""
    files = ("network.txt", "contig_data_network.txt")
    try:
        return mtn.metacool_to_contacts(metacool, contig_data, outdir, *files, tmp_dir, threads, normalization)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise


def check_assembly(assembly):
    """Refuse, before any work, an assembly (-a) that is not a FASTA file,
    plain or compressed with gzip, as pyfastx reads it in the partition and
    the validation (metator.io.check_is_fasta reads bzip2 and zip files
    too): the reason is logged as an error, then ValueError is raised."""
    compression = mio.file_compression(assembly) if os.path.isfile(assembly) else None
    if compression in ("bz2", "zip"):
        logger.error(
            f"{assembly} is compressed with {'bzip2' if compression == 'bz2' else 'zip'}: give the assembly as a "
            "FASTA file, plain or compressed with gzip."
        )
        raise ValueError
    if not mio.check_is_fasta(assembly):
        logger.error(f"{assembly} is not a FASTA file: give the assembly that the map was made from.")
        raise ValueError


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
    in a .metacool file, and the columns that a run read from its map
    (columns_from_map, see metator.network.contig_data_from_map): where the
    restriction sites and the coverage that the map then stores came from,
    which --enzyme and --depth do not say."""
    params = dict(click.get_current_context().params)
    metadata = {
        "command": shlex.join(["metator", *sys.argv[1:]]),
        "parameters": {name: list(value) if isinstance(value, tuple) else value for name, value in params.items()},
    }
    if columns_from_map:
        metadata["columns_from_map"] = list(columns_from_map)
    return metadata


def hold_lock(held, path, in_place):
    """Hold the lock of the map at path until the end of the run (in held, an
    ExitStack): metacooler's (metacooler.lock), which every writer of the map
    takes too, and which the edit session that writes the map at the end of
    the run re-enters (see metator.metacool._staged_edit). The run keeps the
    other writers out from its checks to its end. path is the map that the
    run writes: its input map, updated in place (in_place: the lock of the
    file that a symbolic link points to, as metacooler edits it), or
    --output-map's file, a map of its own, whose own name is locked (a
    symbolic link there is replaced, not its target). A run with
    --output-map does not lock its input map, so that runs from one map,
    each with its own --output-map, run side by side: their input map is
    refused at the end of the run if it changed since the run read it (see
    metator.metacool.reassign_metacool's planned).

    A map that another process holds, or whose lock cannot be created (a
    directory that is missing or not writable), is refused before any
    work: metacooler's reason (who holds the lock, or the files to delete)
    is logged as an error, and, for a run in place, how to write its map
    elsewhere, then a MetacoolError is raised. Taking the lock, metacooler
    takes over a lock of a process of this host that is gone, and removes
    what other writers left beside the map, with a warning each, in the log
    of the run. Neither metaTOR nor metacooler handles signals: a run killed
    outright (SIGKILL, or SIGTERM, as a job scheduler sends at its time
    limit) leaves its lock and its staged copy, which the next writer of
    the map on this host takes over and removes; a lock of another host is
    never taken over, and the error names the files to delete."""
    import metacooler

    try:
        # metacooler's warnings (a stale lock taken over, what a killed run
        # left beside the map) go to the log of the run.
        with mtmc._metacooler_logging():
            command = shlex.join(["metator", *sys.argv[1:]])
            held.enter_context(metacooler.lock(path, command=command, follow_symlinks=in_place))
    except OSError as error:
        # MapLockError, and MapLockedError: OSErrors of metacooler.
        elsewhere = " Or write this run's map elsewhere with --output-map." if in_place else ""
        message = f"{str(error).rstrip('.')}.{elsewhere}"
        logger.error(message)
        raise mtmc.MetacoolError(message) from error


def hold_map(held, command, input_map, outfile, in_place):
    """Hold the lock of the map that a run (command) writes, outfile, until
    the end of the run (see hold_lock), then open its input map: the map as
    the run reads it, under that lock when the run updates it in place
    (in_place), which the map of the run is refused to be written from if
    it changed since (see metator.metacool.reassign_metacool's planned). A
    run in place then refuses a read-only map, and warns of what it
    replaces there (see refuse_a_read_only_map and warn_of_what_is_replaced).
    Returns the input map, opened (a metacooler.MetaCooler)."""
    hold_lock(held, outfile, in_place)
    planned = open_input_map(input_map)
    if in_place:
        refuse_a_read_only_map(outfile)
        warn_of_what_is_replaced(command, input_map)
    return planned


def refuse_a_read_only_map(path):
    """Refuse, before any work, the map that a run updates in place (path,
    its input map) when it has no write permission (judged by its mode
    too, as a run as root may write any file): metacooler refuses to
    replace a read-only map, which it would do at the end of the run. The
    reason is logged as an error, then MetacoolError is raised."""
    mode = os.stat(path).st_mode
    if mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH) and os.access(path, os.W_OK):
        return
    message = (
        f"{path} is read-only (mode {stat.filemode(mode)}): the run would replace it with the map of its clusters. "
        "Make it writable (chmod u+w), or write the map of the run to another file with --output-map."
    )
    logger.error(message)
    raise mtmc.MetacoolError(message)


def run_map_path(input_map, output_map=None, force=False, no_metacool=False):
    """The map a run writes, and whether it updates its input map in place:
    output_map when given, else the file of its input map, the file that it
    points to when it is a symbolic link, as metacooler edits that file (a
    map given as file::/, the root of the file, is the file). output_map
    updates the input map in place when it names that file itself (see
    metator.metacool.is_input_map); any other output_map is a map of its
    own, a symbolic link there replaced, not its target. The directory of
    output_map is made, as the map will be written there.

    Raises
    ------
    MetacoolError
        If output_map is given with no_metacool (no map is written), names
        an entry of a collection or a directory, is the collection of the
        input map (an entry) or a hard link of the input map, or exists and
        is not the input map, without force; if the input map is an entry
        of a .mmetacool collection, without output_map: updating it in place
        would rewrite the whole collection; or if the directory of
        output_map cannot be made. The input map is there: run_input and
        step_map refuse a missing one first.
    """
    reason = None
    if output_map and no_metacool:
        reason = f"--no-metacool writes no map: drop --output-map {output_map}, or --no-metacool."
    elif output_map and "::" in str(output_map):
        reason = f"--output-map takes the path of a .metacool file, not an entry of a collection: {output_map}."
    elif output_map and os.path.isdir(output_map):
        reason = f"{output_map} is a directory: --output-map takes the path of a .metacool file."
    elif output_map and mtmc.is_entry(input_map) and mtmc.is_input_map(input_map, output_map):
        reason = (
            f"{output_map} is the collection of the input map {input_map}: a run does not replace a collection with "
            "the map of its clusters; give --output-map another file."
        )
    elif output_map and mtmc.is_hard_link_of(input_map, output_map):
        reason = (
            f"--output-map {output_map} is a hard link of the input map {input_map}: the same file under another "
            "name. Give --output-map another file, or leave it out to update the map in place."
        )
    elif output_map and os.path.lexists(output_map) and not force and not mtmc.is_input_map(input_map, output_map):
        reason = f"{output_map} exists: give -F to replace it with the map of this run, or another --output-map."
    if reason:
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    if output_map:
        try:
            mtmc.make_map_directory(output_map, f"--output-map {output_map}")
        except mtmc.MetacoolError as error:
            logger.error(str(error))
            raise
        return output_map, mtmc.is_input_map(input_map, output_map)
    if mtmc.is_entry(input_map):
        message = (
            f"{input_map} is an entry of a .mmetacool collection, which a run does not update in place: write the "
            "map of its clusters to another file with --output-map."
        )
        logger.error(message)
        raise mtmc.MetacoolError(message)
    # The file itself (file::/, the root of the file, names it too), or the
    # file that a symbolic link points to, which metacooler edits.
    path = mtmc.map_file(input_map)
    return (os.path.realpath(path) if os.path.islink(path) else path), True


def _failure(error):
    """What the failure of a map says: the message of a MetacoolError, or,
    for another error, which is unexpected, its type and message. Returns
    it without its final period, and whether the error is unexpected (its
    traceback then goes to the log too)."""
    unexpected = not isinstance(error, mtmc.MetacoolError)
    reason = f"{type(error).__name__}: {error}" if unexpected else str(error)
    return reason.rstrip("."), unexpected


@contextmanager
def writing_map(metacool, outfile, what="the map of the clusters"):
    """Write a map of a run at outfile in the with block (see run_map_path
    and step_map), from metacool, the input map of the run, which outfile is
    by default.

    outfile is only replaced once the new map has been written and
    validated: an input map at outfile (the file of the input map, see
    metator.metacool.is_input_map) stays as it was when the new map cannot
    be written. Another map of the same name, left by an earlier run, is
    removed first: this run has replaced the tables it was made from, so it
    must not stay behind, passing for this run's, when the new map cannot
    be written. A symbolic link at --output-map is removed, not its target.
    The log says what (the map of the clusters, by default) replaces an
    input map.

    Raises
    ------
    MetacoolError
        If the map cannot be written: the reason, and what is kept, is
        logged as an error first. The other outputs of the run are kept.
    """
    # What the error adds when the new map is not written.
    kept = ""
    try:
        if mtmc.is_input_map(metacool, outfile):
            kept = "It is still the input map of the run, as it was. "
            logger.info(f"{outfile} is the input map {metacool}: it is replaced once {what} is validated.")
        elif os.path.lexists(outfile):
            logger.info(f"Removing {outfile}, the contact map of an earlier run.")
            os.remove(outfile)
        yield
    except Exception as error:
        reason, unexpected = _failure(error)
        message = f"{outfile} is not written: {reason}. {kept}The other outputs of the run are kept."
        logger.error(message, exc_info=unexpected)
        raise mtmc.MetacoolError(message) from error


def pipeline_metacool(metacool, outdir, outfile, columns_from_map=(), normalization=None, planned=None):
    """Write the contact map of the clusters of a pipeline or validation run
    at outfile (see writing_map), from the .metacool map of its contigs, its
    input map, which outfile is by default, and the final tables of the run
    (see metator.metacool.reassign_metacool), with its clustering rounds
    (see metator.clusters), and the columns that it read from the map
    (columns_from_map: stored in the metadata of the new map). The
    normalization of the network of a pipeline run is stored in the map, as
    the weights metator_<normalization> of its bins (see
    metator.metacool.store_normalization). The map is refused if it changed
    since planned, the input map as the run read it at its start, was
    opened.

    Raises
    ------
    MetacoolError
        If the map cannot be written: the reason, and what is kept, is
        logged as an error first. The other outputs of the run are kept.
    """
    with writing_map(metacool, outfile):
        mtmc.reassign_metacool(
            metacool,
            join(outdir, "contig_data_final.txt"),
            outfile,
            bin_summary=join(outdir, "cluster_summary.txt"),
            force=True,
            metadata=run_metadata(columns_from_map),
            rounds=mtcl.read_rounds(outdir),
            normalization=normalization,
            planned=planned,
        )


def partition_metacool(metacool, outfile, contig_data_file, size, prefix, planned=None):
    """Store the overlapping clusters of a partition run from a map in the
    map at outfile (see writing_map): the round overlapping of its
    clustering rounds (see metator.clusters.overlapping_round), the
    overlapping clusters that reach size (--size), those that the run
    writes a FASTA file for, with their contigs, in place of the rounds that
    the map held. Its genomes and contigs tables are left as they were (see
    metator.metacool.write_rounds): the validation writes the map of the
    final clusters. The map is refused if it changed since planned, the
    input map as the run read it at its start, was opened.

    Raises
    ------
    MetacoolError
        If the map cannot be written: the reason, and what is kept, is
        logged as an error first. The other outputs of the run are kept.
    """
    with writing_map(metacool, outfile, what="the map with its overlapping clusters"):
        rounds = mtcl.overlapping_round(contig_data_file, size, prefix)
        mtmc.write_rounds(metacool, outfile, rounds, names=[mtcl.OVERLAPPING], force=True, planned=planned)


def mge_metacool(metacool, outfile, summary_file, planned=None):
    """Store the MGE clusters of an mge run in the map at outfile (see
    writing_map): the round mge of its clustering rounds (see
    metator.clusters.mge_round), from the summary that the run wrote
    (summary_file), in place of the round mge that the map held, after the
    rounds of the binning run, which are kept as they were. Its genomes and
    contigs tables, the final clusters, are left as they were, and so is
    everything else (see metator.metacool.write_round). The round records
    the version of metaTOR, and the command line and the parameters of the
    run. The map is refused if it changed since planned, the input map as
    the run read it at its start, was opened.

    Raises
    ------
    MetacoolError
        If the map cannot be written: the reason, and what is kept, is
        logged as an error first. The other outputs of the run are kept.
    """
    with writing_map(metacool, outfile, what="the map with its MGE clusters"):
        entry = {"name": mtcl.MGE, "tool": "metaTOR", "version": __version__, "step": "mge", **run_metadata()}
        mtmc.write_round(metacool, outfile, mtcl.mge_round(summary_file), entry, force=True, planned=planned)


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


def pairsfile_argument(help):
    return click.argument("pairsfile", nargs=-1, required=True, help=help, panel="Input")


ASSEMBLY_HELP = "The assembly that the map was made from: a FASTA file, plain or compressed with gzip."
ASSEMBLY_PARTITION_HELP = "The assembly that the contigs come from: a FASTA file."


def assembly_option(help):
    return _option("-a", "--assembly", required=True, metavar="FILE", panel="Input", help=help)


MAP_HELP = (
    "A .metacool contact map of the contigs of the assembly, made from the pairs of the run (by hicstuff pipeline -M "
    "metacool or metacooler create pairs), or an entry of a .mmetacool collection (file.mmetacool::entries/NAME): "
    "the run reads its network from it."
)


def map_argument(help=MAP_HELP, required=False):
    """The .metacool map of the contigs that a network or pipeline run starts
    from, its first argument (required); or that a partition or validation
    run updates (see step_map)."""
    return click.argument(
        # Shown as MAP: the parameter of these commands is named so.
        "map",
        required=required,
        # Shown as FILE; a directory is refused as a usage error. An entry of a
        # collection (file::group) is no path: run_input and step_map check it.
        type=click.Path(dir_okay=False),
        panel="Input",
        help=help,
    )


def output_map_option(help):
    """--output-map: the file that a run writes the map to, instead of the
    map of the run (see run_map_path)."""
    return _option("--output-map", metavar="FILE", help=help)


def network_options(normalization_short):
    """Options of the network step: the assembly, what is measured on its
    contigs, and the normalization; the short alias of --normalization is
    -n in network, -m in pipeline."""
    return _options(
        assembly_option(ASSEMBLY_HELP),
        _option(
            "-d",
            "--depth",
            metavar="FILE",
            panel="Input",
            help="The depth.txt file from the shotgun reads used to made the assembly computed by "
            "[cyan]jgi_summarize_bam_contig_depths[/] from the metabat2 pipeline. The Shotgun_coverage column of the map "
            "gives the coverage without it, and --depth replaces it.",
        ),
        _option(
            "-e",
            "--enzyme",
            metavar="STR",
            panel="Network",
            help="The list of restriction enzyme used to digest the contigs separated by a comma. Example: HpaII,MluCI. "
            "The Restriction_site column of the map gives the restriction sites without it, and --enzyme replaces it.",
        ),
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
            "[cyan]theoritical_hit[/] need --depth. The Restriction_site and Shotgun_coverage columns of the map stand "
            "for them (see --enzyme and --depth).",
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
            help="Threshold size to keep clusters in base pair.",
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


def tmpdir_option(help="Temporary directory."):
    return _option("-T", "--tmpdir", default="./tmp", show_default=True, metavar="DIR", help=help)


def threads_option(help):
    return _option("-t", "--threads", type=click.INT, default=1, show_default=True, metavar="INT", help=help)


def no_clean_up_option(help="Do not remove temporary files."):
    return _option("-N", "--no-clean-up", is_flag=True, help=help)


# ---------------------------------------------------------------------------
# CLI group
# ---------------------------------------------------------------------------


@click.group("metator", context_settings={"help_option_names": ["-h", "--help"]})
@click.rich_config(help_config=click.RichHelpConfiguration(text_markup="rich", show_arguments=True, max_width=100))
@click.command_panel("Clustering", commands=["network", "partition", "validation", "pipeline"])
@click.command_panel("Other commands", commands=["mge"])
@click.version_option(__version__, "-v", "--version", prog_name="metator", message="%(version)s")
def cli():
    """MetaHiC pipeline for generating and manipulating MAGs.

    The clustering runs in three steps, in this order: [cyan]network[/],
    [cyan]partition[/] and [cyan]validation[/]; [cyan]pipeline[/] runs them
    all in one command line.
    """


# ---------------------------------------------------------------------------
# network
# ---------------------------------------------------------------------------


@cli.command(short_help="Generate the metaHiC contig network from a .metacool map and normalize it.")
@map_argument(required=True)
@network_options("-n")
@outdir_option(help="The output directory to write the network and contig data into. Default: current directory.")
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated to sort the network.")
@no_clean_up_option()
def network(
    map,  # the argument MAP
    assembly,
    depth,
    enzyme,
    normalization,
    outdir,
    tmpdir,
    threads,
    no_clean_up,
):
    """Generate the metaHiC contig network.

    [cyan]metator network MAP -a ASSEMBLY[/] builds the network from a
    .metacool contact map of the contigs, whose libraries are summed (see
    the README): that of the pairs it was made of. hicstuff
    ([cyan]hicstuff pipeline -M metacool[/]) or metacooler
    ([cyan]metacooler create pairs[/]) makes the map from the reads or the
    pairs.

    The Restriction_site and Shotgun_coverage columns of the map, when they
    give every contig a value, give the restriction sites and the coverage
    of the contigs that --enzyme and --depth would measure; given, these
    options replace them.

    In the output network, contigs are the network nodes and the edges are
    the contact counts between two contigs: every pair counts, as in a
    contact map, but the pairs within a contig.

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

        # The map, the normalization against its columns, and the assembly.
        run_input(map)
        check_normalization(normalization, enzyme, depth, map)
        check_assembly(assembly)

        # Print information of the workflow:
        logger.info("Enzyme: %s", enzyme)
        logger.info("Normalization: %s", normalization)

        # The contig data, measured on the assembly and read from the map,
        # then the network of the map.
        contig_data = mtn.create_contig_data(assembly, depth, enzyme)
        contig_data_from_map(map, contig_data, normalization, enzyme, depth)
        contig_network(map, contig_data, outdir, tmp_dir, threads, normalization)

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# partition
# ---------------------------------------------------------------------------


@cli.command(short_help="Partition the metaHiC network with the Louvain or Leiden algorithm.")
@map_argument(
    help="The .metacool contact map of the contigs of the run (that of network MAP), or an entry of a .mmetacool "
    "collection (file.mmetacool::entries/NAME): the overlapping clusters of the run are stored in its clustering "
    "rounds, in place, or in --output-map, which an entry needs.",
)
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
    overlap_help="Hamming distance threshold to use to merge clusters (percentage).",
)
@outdir_option(default="./", help="Path to the directory to write the output.")
@prefix_option()
@force_option(
    "Replace the directory of overlapping clusters of an earlier run in the output directory (overlapping_clusters), "
    "and an existing --output-map."
)
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated for the partition.")
@no_clean_up_option()
@output_map_option(
    help="Store the overlapping clusters of the run in a copy of MAP written to FILE (replaced only with -F), and "
    "leave MAP as it was. Only with MAP, which is updated in place by default.",
)
def partition(
    map,  # the argument MAP
    assembly,
    contigs,
    network,
    algorithm,
    iterations,
    overlap,
    res_param,
    size,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
    output_map,
):
    """Partition the network using the Louvain or Leiden algorithm.

    Partition the network file using iteratively the Louvain or Leiden
    algorithm. Then looks for 'core' clusters that are constituted by the
    contigs which are always in the same community at each iterations. Then
    using the Hamming distance from these core clusters, merge the clusters
    with a small Hamming distance (close to 1), bigger than the percentage
    (overlapping parameter) given.

    It will also update the file to integrate the clusters of the contigs.

    Furthermore, both Leiden and Louvain algorithm are available here.
    However, the benchmark made show that here the Louvain algorithm have
    better performance and is faster on seawater and gut metagenomic
    samples.

    [cyan]metator partition MAP ...[/] runs on the map of the contigs of a
    run step by step ([cyan]metator network MAP[/], partition, validation):
    it stores the overlapping clusters that reach --size, with their
    contigs, in the clustering rounds of the map (its /clusters group), in
    place of those it held, and updates the map in place (or writes
    --output-map). Its genomes and contigs tables are left as they were:
    [cyan]metator validation MAP[/] writes the final clusters there.
    """
    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up), ExitStack() as held:
        # Defined the output directory.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)
        fasta_dir = join(outdir, "overlapping_clusters")

        # Enable file logging, the map of a run from a map checked and
        # locked, and the directory of the clusters made (see start_step).
        log_file, metacool_file, bins_prefix, planned = start_step(
            "partition",
            outdir,
            [fasta_dir],
            map,
            output_map,
            force,
            held,
            prefix,
            contigs,
            [assembly, network],
            assembly,
        )

        # Partition the network
        contig_data_file = mtp.partition(
            algorithm,
            assembly,
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

        # The overlapping clusters, stored in the map of a run from a map.
        if metacool_file is not None:
            partition_metacool(map, metacool_file, contig_data_file, size, bins_prefix, planned)

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


@cli.command(short_help="Validate the clusters with miComplete and decontaminate them with a recursive partition.")
@map_argument(
    help="The .metacool contact map of the contigs of the run (that of network MAP and partition MAP), or an entry "
    "of a .mmetacool collection (file.mmetacool::entries/NAME): it becomes the map of the clusters of the run, as "
    "pipeline MAP writes it, in place, or in --output-map, which an entry needs.",
)
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
    help="Path to the directory containing the input fasta files of the clusters to decontaminate.",
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
    overlap_help="Hamming distance threshold to use to merge clusters (percentage).",
)
@outdir_option(default="./", help="Path to the directory to write the output.")
@prefix_option()
@force_option(
    "Replace the directories of recursive and final clusters of an earlier run in the output directory "
    "(recursive_clusters, final_clusters), and an existing --output-map."
)
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated for the partition.")
@no_clean_up_option()
@output_map_option(
    help="Write the contact map of the clusters to FILE (replaced only with -F), and leave MAP as it was. Only with "
    "MAP, which is updated in place by default.",
)
def validation(
    map,  # the argument MAP
    assembly,
    contigs,
    fasta,
    network,
    algorithm,
    iterations,
    overlap,
    res_param,
    size,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
    output_map,
):
    """Use miComplete to validate the clusters.

    Recursively decontaminate the contaminated clusters evaluated by
    miComplete. The command returns new decontaminated fasta and summary
    files of the decontamination.

    Only clusters with a weighted completeness of at least 33% and a
    contamination, (weighted redundancy - 1) / weighted completeness, of at
    least 5% are subject to the recursive step. Its clusters replace the
    original cluster when the most complete of them keeps enough of its
    completeness; otherwise the original cluster is kept.

    [cyan]metator validation MAP ...[/] runs on the map of the contigs of a
    run step by step ([cyan]metator network MAP[/], partition MAP,
    validation): at its end, it writes the map of the clusters of the run,
    with all its clustering rounds, as [cyan]metator pipeline MAP[/] does,
    in place of the map (or to --output-map).
    """
    # Launch miComplete to evaluate the completion and the contamination. If
    # asked rerun Louvain to try to reduce the contamination, rerun miComplete
    # if the contamination decrease without a huge decrease of the completion
    # keep the new bins. Otherwise go back to the old state.

    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir, ExitStack() as held:
        # Defined the output directory and output file names.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)
        recursive_fasta_dir = join(outdir, "recursive_clusters")
        final_fasta_dir = join(outdir, "final_clusters")
        # The clusters to validate in a directory that -F removes: refused
        # before anything is removed.
        refuse_an_input_in_the_outputs(fasta, [recursive_fasta_dir, final_fasta_dir], force, "-f")

        # Enable file logging, the map of a run from a map checked and
        # locked, and the directories of the clusters made (see start_step).
        log_file, metacool_file, bins_prefix, planned = start_step(
            "validation",
            outdir,
            [recursive_fasta_dir, final_fasta_dir],
            map,
            output_map,
            force,
            held,
            prefix,
            contigs,
            [assembly, network],
            assembly,
        )

        mtv.recursive_decontamination(
            algorithm,
            assembly,
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

        # The map of the clusters, written last, as pipeline MAP writes it.
        if metacool_file is not None:
            pipeline_metacool(map, outdir, metacool_file, planned=planned)

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------


@cli.command(short_help="Run network, partition and validation in one command line.")
@map_argument(required=True)
@network_options("-m")
@partition_options(
    iterations=100,
    overlap=80,
    iterations_help="Number of iterations of Louvain or Leiden for the partition step.",
    overlap_help="Hamming distance threshold to use to merge clusters for the first step (percentage).",
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
    help="Hamming distance threshold to use to merge clusters for the recursive step (percentage).",
)
@outdir_option(help="The output directory. Default: current directory.")
@prefix_option()
@force_option("Replace the directories of clusters of an earlier run in the output directory, and an existing --output-map.")
@tmpdir_option()
@threads_option(help="Number of parallel threads allocated to sort the network, and for the partition and the validation.")
@no_clean_up_option()
@_option(
    "--no-metacool",
    is_flag=True,
    help="Do not write the contact map of the clusters: the map of the run is left as it was.",
)
@output_map_option(
    help="Write the contact map of the clusters to FILE (replaced only with -F), and leave MAP as it was. By default, "
    "the run updates MAP in place.",
)
def pipeline(
    map,  # the argument MAP
    assembly,
    depth,
    enzyme,
    normalization,
    algorithm,
    iterations,
    overlap,
    res_param,
    size,
    rec_iter,
    rec_overlap,
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
    map of the contigs, made from the pairs of the run (by
    [cyan]hicstuff pipeline -M metacool[/] or
    [cyan]metacooler create pairs[/]): the run reads its network from that
    map, and the restriction sites and the coverage of the contigs from its
    Restriction_site and Shotgun_coverage columns, unless --enzyme and
    --depth are given; it updates the map in place with the map of its
    clusters (or writes --output-map).

    The network counts every pair, as a contact map does. A validation step
    decontaminates the clusters when it's necessary; it is the critical step
    for memory usage (~40G).
    """
    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir, ExitStack() as held:
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
        # clusters of an earlier run (--force): a refused run leaves them as
        # they were.
        run_input(map)
        metacool_file, in_place = run_map_path(map, output_map, force, no_metacool)
        refuse_an_input_as_the_map(output_map, [assembly, depth], assembly)
        check_normalization(normalization, enzyme, depth, map)
        bins_prefix = bin_prefix(prefix)
        # A run that writes its map holds the map's lock from now on: a second
        # run on the same map is refused before any work. It reads its map
        # under that lock (planned).
        planned = None if no_metacool else hold_map(held, "pipeline", map, metacool_file, in_place)

        # Print information of the workflow:
        logger.info(f"Partition algorithm: {algorithm}")
        logger.info(f"Partition iterations: {iterations}")
        logger.info(f"Overlapping parameter: {overlapping_parameter}")
        logger.info(f"Recursive partition iterations: {rec_iter}")
        logger.info(f"Recursive overlapping parameter: {recursive_overlapping_parameter}\n")

        # The clusters of an earlier run: refused now without --force, before
        # the assembly is read, and replaced only once it is checked.
        overlapping_fasta_dir = join(outdir, "overlapping_clusters")
        recursive_fasta_dir = join(outdir, "recursive_clusters")
        final_fasta_dir = join(outdir, "final_clusters")
        bin_dirs = [overlapping_fasta_dir, recursive_fasta_dir, final_fasta_dir]
        for path in bin_dirs:
            refuse_output_dir(path, force)

        # The assembly, and the contig data of the run, before the clusters
        # of an earlier run are replaced: a map whose columns cannot give
        # what the normalization needs is refused once the assembly is read.
        check_assembly(assembly)
        contig_data = mtn.create_contig_data(assembly, depth, enzyme)
        columns_from_map = contig_data_from_map(map, contig_data, normalization, enzyme, depth)
        for path in bin_dirs:
            make_output_dir(path, force)

        # Build the network, from the map.
        network_file, contigs_data_file = contig_network(map, contig_data, outdir, tmp_dir, threads, normalization)

        # Partition the network
        contigs_data_file = mtp.partition(
            algorithm,
            assembly,
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

        # Launch validation.
        mtv.recursive_decontamination(
            algorithm,
            assembly,
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

        # Remove contig_data_partition file
        os.remove(join(outdir, "contig_data_partition.txt"))

        # Delete pyfastx index:
        mio.remove_pyfastx_index(assembly)

        # The map of the clusters, written last, from the input map.
        if not no_metacool:
            pipeline_metacool(
                map,
                outdir,
                metacool_file,
                columns_from_map=columns_from_map,
                normalization=normalization,
                planned=planned,
            )

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# mge
# ---------------------------------------------------------------------------


@cli.command(short_help="Cluster the MGE contigs into MGE MAGs from the metaHiC pairs.")
@pairsfile_argument(
    help="Pairs file(s) of the run: the contacts within each MGE contig are counted from the positions of their reads, "
    "which a .metacool map does not hold (give the map to store the MGE clusters in with --map)."
)
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
    "--map",
    # A file, or an entry of a collection (file::group), which is no path:
    # step_map checks it.
    type=click.Path(dir_okay=False),
    metavar="MAP",
    panel="Input",
    help="A .metacool contact map of the run (that of pipeline, for example), or an entry of a .mmetacool collection "
    "(file.mmetacool::entries/NAME): the MGE clusters are stored in its clustering rounds, as the round mge, in place "
    "of an earlier one, in place (or in --output-map). Its genomes and contigs tables are left as they were.",
)
@_option(
    "-s",
    "--threshold-bin",
    type=click.FLOAT,
    default=0.8,
    show_default=True,
    metavar="FLOAT",
    panel="Clustering",
    help="Threshold to use for clustering.",
)
@_option("-r", "--random", is_flag=True, panel="Clustering", help="If enabled, make a random clustering. Refused with --map.")
@outdir_option(help="Path to the output directory where the output will be written. Default: current directory.")
@force_option("Replace an existing --output-map.")
@tmpdir_option(help="Path to temporary directory.")
@no_clean_up_option(help="If enabled, intermediary files will be kept.")
@output_map_option(
    help="Store the MGE clusters in a copy of MAP written to FILE (replaced only with -F), and leave MAP as it was. Only "
    "with --map, whose map is updated in place by default.",
)
def mge(
    pairsfile,
    assembly,
    contigs,
    binning,
    mges,
    map,  # --map
    threshold_bin,
    random,
    outdir,
    force,
    tmpdir,
    no_clean_up,
    output_map,
):
    """Cluster MGE contigs.

    Cluster the MGE contigs into MGE MAGs using pairs, host MAG sequence and
    information.

    Command will return the updated MGE data and the MGE fasta.

    The MGE fasta will contain one entry by MGE MAG, with 200 "N" spacers
    between contigs.

    [cyan]metator mge ... --map MAP[/] also stores the MGE clusters in the
    .metacool contact map of the run, MAP: as the round mge of its
    clustering rounds (its /clusters group), with their contigs and the
    columns of mges_cluster_summary.tsv, in place of an earlier round mge, after
    the rounds of the clustering run, which are kept. Its genomes and contigs
    tables, the final clusters, are left as they were. MAP is updated in
    place (or written to --output-map), and locked while the run writes it.
    The contacts of the MGE contigs still come from the pairs files: those
    within a contig are counted only between reads more than 1 kb apart,
    which a map does not tell.
    """
    # Generate log
    log_file = log_file_name("", "mge")
    generate_log_header(log_file)

    # The contacts within each MGE contig need the positions of the reads.
    maps = [path for path in pairsfile if mtmc.is_map(path)]
    if maps:
        reason = (
            f"{maps[0]} is a contact map, not a pairs file: metator mge counts the contacts within each MGE contig "
            "between reads more than 1 kb apart, which a map does not tell. Give the pairs files of the run, and the "
            "map to store the MGE clusters in with --map."
        )
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    if random and map is not None:
        reason = "--random makes a random clustering, to compare the MGE clusters with: it is not stored in a map. Drop --map."
        logger.error(reason)
        raise mtmc.MetacoolError(reason)

    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up), ExitStack() as held:
        # The map that the run stores its MGE clusters in, checked and
        # locked.
        metacool_file, planned = step_map(
            "mge",
            map,
            output_map,
            force,
            held,
            contigs,
            # The pairs files and their pypairix index too.
            [assembly, binning, mges, *pairsfile, *(f"{path}.px2" for path in pairsfile)],
            option="--map",
            assembly=assembly,
        )

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

        # The MGE clusters, stored in the map.
        if metacool_file is not None:
            mge_metacool(map, metacool_file, join(outdir, "mges_cluster_summary.tsv"), planned)

    if not no_clean_up:
        # Delete pyfastx index:
        mio.remove_pyfastx_index(assembly)

    generate_log_footer(log_file)
