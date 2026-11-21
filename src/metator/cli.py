#!/usr/bin/env python3
# coding: utf-8

"""metaTOR command line interface.

One function per command: network, partition, validation, pipeline, mge
and dump. The ``metator`` entry point is :func:`cli`.
"""

import glob
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
import metator.dump as mtd
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

    miComplete drops every dot of a FASTA file name, so only letters, digits, _
    and - are accepted."""
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

    --tmpdir may hold the user's own files, so the run works in a fresh random
    subdirectory and removes only that, and --tmpdir only if it made it and it
    is empty."""
    base = tmpdir or "./tmp"
    base_created = not exists(base)
    tmp_dir = mio.generate_temp_dir(base)
    logger.debug("Temporary directory: %s", tmp_dir)
    yield tmp_dir
    if clean_up:
        logger.info(f"The temporary directory of the run, {tmp_dir}, is removed with its files: -N keeps them.")
        shutil.rmtree(tmp_dir)
        if base_created and not os.listdir(base):
            os.rmdir(base)
    else:
        logger.info(f"The temporary directory of the run, {tmp_dir}, and its files are kept (-N).")


# What a validation or pipeline run writes in its output directory besides its
# log: final_clusters, which -F replaces.
CLUSTER_DIRS = ("final_clusters",)

# An earlier run's results in the output directory of a validation or pipeline
# run (2.0 development builds' tables and recursive clusters, 1.x's outputs),
# which would pass for this run's: refused with or without -F, never removed.
EARLIER_OUTPUTS = (
    "contig_data_*.txt",
    "cluster_summary.txt",
    "binning.txt",
    "clusters.txt",
    "cluster_members.txt",
    "overlapping_micomplete_results.txt",
    "recursive_micomplete_step_*.txt",
    "recursive_clusters",
    "bin_summary.txt",
    "network.txt",
    "clustering_matrix*.npz",
    "recursive_bin",
    "final_bin",
    "final_bin_unscaffold",
    "overlapping_bin",
    "scaffold_bin",
    "plot",
)

# An earlier mge run's summary (1.x's, a development build's), which the round
# mge now holds: refused with or without -F, never removed.
MGE_EARLIER_OUTPUTS = ("mges_bin_summary.tsv", "mges_cluster_summary.tsv")


def _within(path, directory):
    """Whether path is directory, or in it (both absolute)."""
    return os.path.commonpath([path, directory]) == directory


def _earlier_directories():
    """The directories among EARLIER_OUTPUTS."""
    return tuple(name for name in EARLIER_OUTPUTS if "." not in name)


def refuse_tmpdir_in_the_clusters(tmpdir, outdir):
    """Refuse, as a usage error, before the run creates or removes anything,
    a --tmpdir in a directory of the output directory, outdir, that a
    validation or pipeline run makes anew (final_clusters, which -F replaces,
    with the temporary directory of the run in it) or refuses as the results
    of an earlier run (see EARLIER_OUTPUTS: the run would make it, then
    refuse it). --tmpdir is in it as written, or once the symbolic links of
    both are resolved, but for that directory itself: -F replaces a symbolic
    link, not its target."""
    base = tmpdir or "./tmp"
    for name in (*CLUSTER_DIRS, *_earlier_directories()):
        directory = join(outdir or ".", name)
        as_written = _within(os.path.abspath(base), os.path.abspath(directory))
        resolved = _within(os.path.realpath(base), join(os.path.realpath(outdir or "."), name))
        if as_written or resolved:
            where = (
                "where the run writes its clusters, which it makes anew (-F replaces an earlier one)"
                if name in CLUSTER_DIRS
                else "where an earlier version wrote its results, which the run refuses"
            )
            raise click.BadParameter(
                f"{base} is in {directory}, {where}: give a temporary directory outside it.",
                param_hint="'-T' / '--tmpdir'",
            )


def _is_or_holds(path, kept):
    """Whether path is one of the files of kept (the maps of the run), or a directory holding one."""
    real = os.path.realpath(path)
    return any(real == os.path.realpath(file) or _within(os.path.realpath(file), real) for file in kept if file)


def refuse_earlier_outputs(outdir, outputs=EARLIER_OUTPUTS, keep=()):
    """Refuse, with or without -F, before any work, an earlier version's results
    in outdir (outputs: EARLIER_OUTPUTS, or MGE_EARLIER_OUTPUTS for mge), which
    this run neither writes nor replaces: the reason, which names them, is logged
    as an error, then ValueError is raised, and nothing is removed. A symbolic
    link counts, even one to nothing; the maps of the run (keep) never do, nor a
    directory that holds one."""
    kept = [mtmc.map_file(path) if path else None for path in keep]
    found = sorted(
        path
        for path in {path for pattern in outputs for path in glob.glob(join(glob.escape(outdir), pattern))}
        if not _is_or_holds(path, kept)
    )
    if found:
        logger.error(
            f"{outdir} holds results of an earlier run, which this run does not write and does not replace: "
            f"{mtmc._shown([os.path.basename(path) for path in found])}. Move them away, or give the run another "
            "output directory with -o."
        )
        raise ValueError


def refuse_output_dir(path, force):
    """Refuse an existing output directory of bins, unless --force replaces it: a
    symbolic link too, even one to nothing."""
    if os.path.lexists(path) and not force:
        logger.error(
            "%s already existed. Remove directory or use -F argument to overwrite it.",
            path,
        )
        raise ValueError


def make_output_dir(path, force):
    """Create an output directory of bins; an existing one is replaced with
    --force, or refused."""
    refuse_output_dir(path, force)
    if os.path.islink(path):
        os.remove(path)
    elif exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def check_normalization(normalization, enzyme, depth, metacool, read=None):
    """The normalizations that need the restriction sites or the shotgun coverage
    of the contigs: measured on the assembly with --enzyme and --depth or read
    from the Restriction_site and Shotgun_coverage columns of the contigs table
    of the map of the run, metacool (see metator.network.contig_data_from_map)."""
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
    """Read the restriction sites and the shotgun coverage of the contigs from the
    .metacool map that a run starts from, where --enzyme and --depth measure
    none (see metator.network.contig_data_from_map), then check the
    normalization against what the run has."""
    try:
        read = mtn.contig_data_from_map(metacool, contig_data, enzyme, depth)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    check_normalization(normalization, enzyme, depth, metacool, read)
    return read


def run_input(input_map):
    """Check the .metacool map that a network or pipeline run starts from, its
    first argument (input_map), before any work."""
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
    """Refuse a map whose file is not there, or is no file (a directory given as
    dir::/), before any work: the reason is logged as an error, then
    MetacoolError is raised."""
    path = input_map if os.path.isfile(input_map) else mtmc.map_file(input_map)
    if not os.path.isfile(path):
        reason = "no such file" if not os.path.exists(path) else f"{path} is not a file"
        message = f"Cannot read the contact map {input_map}: {reason}."
        logger.error(message)
        raise mtmc.MetacoolError(message)
    return path


def refuse_a_map_named_as_a_uri(input_map, path):
    """Refuse, before any work, a map whose file name holds :: (path, its file, is
    input_map itself), which metacooler reads as file::group: the reason is
    logged as an error, then MetacoolError is raised."""
    if path == input_map and "::" in input_map:
        message = (
            f"Cannot read the contact map {input_map}: metacooler reads a path that holds :: as file::group. "
            "Rename the map without ::."
        )
        logger.error(message)
        raise mtmc.MetacoolError(message)


def open_input_map(input_map):
    """The map, opened (a metacooler.MetaCooler)."""
    try:
        return mtmc.open_map(input_map)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise


def refuse_an_input_as_the_map(output_map, inputs):
    """Refuse, before any work, an --output-map that is one of the files of the
    inputs of the run (the assembly, a table): the run would replace it with
    its map."""
    if output_map is None or not os.path.exists(output_map):
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
    "network": (False, ()),
    "pipeline": (True, ("partition", "validation")),
    "partition": (False, ("partition", "validation")),
    "validation": (True, ("validation",)),
    "mge": (False, ("mge",)),
}


def warn_of_what_is_replaced(command, metacool, rounds=None):
    """Warn, before any work, that a run (command) that updates the map metacool
    in place replaces what the map holds already: its clusters (the genomes
    other than __unbinned__, of an earlier run, or of another binner), and its
    clustering rounds made by an earlier run (see REPLACED_IN_PLACE).

    Neither can be found again in the updated map."""
    clusters, steps = REPLACED_IN_PLACE[command]
    replaced = (lambda entry: entry.get("step") in steps) if rounds is None else (lambda entry: entry.get("name") in rounds)
    try:
        with mtmc.read_errors(metacool):
            mc = mtmc.open_map(metacool)
            n_clusters = int((mc.genomes["name"].astype(str) != "__unbinned__").sum()) if clusters else 0
            rounds = [str(entry["name"]) for entry in mc.cluster_rounds if replaced(entry)]
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


def step_map(command, input_map, output_map, force, held, inputs=(), rounds=None):
    """The map that a partition, validation or mge run (command) writes."""
    map_of = f"the map of the {'run' if command == 'mge' else 'contigs of the run'} as the first argument of {command}"
    path = refuse_missing_map(input_map)
    if not mtmc.is_map(input_map):
        reason = f"{input_map} is not a .metacool contact map: give {map_of}."
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    refuse_a_map_named_as_a_uri(input_map, path)
    open_input_map(input_map)
    refuse_an_input_as_the_map(output_map, inputs)
    outfile, in_place = run_map_path(input_map, output_map, force)
    return outfile, hold_map(held, command, input_map, outfile, in_place, rounds)


def start_step(
    command,
    outdir,
    cluster_dirs,
    input_map,
    output_map,
    force,
    held,
    prefix=None,
    inputs=(),
    assembly=None,
    partition=False,
):
    """Start a partition or validation run (command) in outdir: its log file, the
    map that it writes (see step_map) and its directories of clusters
    (cluster_dirs, see make_output_dir), refused without force when an earlier
    run left them, as, with partition (a validation run), what a run of an
    earlier version wrote in outdir (see refuse_earlier_outputs), which force
    removes."""
    log_file = log_file_name(outdir, command)
    generate_log_header(log_file)
    metacool_file, planned = step_map(command, input_map, output_map, force, held, inputs=inputs)
    stored = stored_network(input_map)
    clusters = map_partition(input_map, planned, stored.digest) if partition else None
    if assembly is not None:
        check_assembly(assembly)
    bins_prefix = bin_prefix(prefix)
    # Refused, all of them, before anything is replaced.
    for path in cluster_dirs:
        refuse_output_dir(path, force)
    if partition:
        refuse_earlier_outputs(outdir, keep=(input_map, output_map))
    for path in cluster_dirs:
        make_output_dir(path, force)
    return log_file, metacool_file, bins_prefix, planned, stored, clusters


def stored_network(metacool):
    """The network that metator network stored in the map of a partition or
    validation run, metacool, read in full before the run replaces anything
    (see metator.network.stored_network): the contig data of its contigs, the
    network and its normalization, by the positions of the contigs in the map,
    which the partition that a validation run reads from the map numbers them
    by too (see metator.metacool.read_partition)."""
    try:
        return mtn.stored_network(metacool)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise


def map_contacts(metacool):
    """The contacts of the map of a network or pipeline run, metacool, summed by
    pair of contigs (see metator.network.contig_contacts), and counted in its
    log."""
    try:
        contacts = mtn.contig_contacts(metacool)
        mtn.log_contacts(metacool, contacts)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    return contacts


def map_partition(metacool, planned, network):
    """The partition that metator partition stored in the map of a validation run,
    metacool, opened as planned (see metator.metacool.read_partition): its core
    and overlapping clusters, which must be those of the network stored in the
    map (network, its digest: see metator.network.stored_network)."""
    try:
        return mtmc.read_partition(metacool, planned, network)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise


def partition_entry(normalization, digest, table):
    """What a map records of the partition of its network, under the partition key
    of its metator metadata (see metator.metacool.write_map): the version of
    metaTOR, the command line and the parameters of the run, what it read of
    the map (read: the network, see map_read; table: the contig data of the
    network), and the normalization and the digest of the network that it
    partitioned (see metator.network.network_digest), which a validation checks
    against the network stored in the map (see
    metator.metacool.read_partition)."""
    return {
        "version": __version__,
        **run_metadata(),
        "normalization": normalization,
        "read": map_read(table, normalization),
        mtmc.NETWORK_DIGEST: digest,
    }


def write_partition_table(contig_data, clusters, directory):
    """Write the contig data table of a partition (see
    metator.partition.partition_table) in directory, the temporary directory of
    a validation or pipeline run, as contig_data_partition.txt, which the
    validation reads (see metator.validation.recursive_decontamination)."""
    path = join(directory, "contig_data_partition.txt")
    mtp.partition_table(contig_data, clusters).to_csv(path, sep="\t", header=True, index=False)
    return path


def network_entry(columns_from_map=()):
    """What a map records of the run that stores the network of its contigs in it,
    under the network key of its metator metadata (see
    metator.metacool.Network): the version of metaTOR, the command line and the
    parameters of the run, and what it read of the map (read): the columns of
    its contigs table that it read in place of --enzyme and --depth (contigs:
    columns_from_map, see metator.network.contig_data_from_map), where the
    restriction sites and the coverage that the map then stores came from."""
    return {"version": __version__, **run_metadata(), "read": {"contigs": list(columns_from_map)}}


def check_assembly(assembly):
    """Refuse, before any work, an assembly (-a) that is not a FASTA file, plain
    or compressed with gzip, as the network, the validation, the pipeline and
    mge read it (Biopython, and metacooler, which writes the FASTA files of the
    clusters; metator.io.check_is_fasta reads bzip2 and zip files too): the
    reason is logged as an error, then ValueError is raised."""
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


def run_metadata():
    """The command line and the parameters of the running command, as stored in a
    .metacool file."""
    params = dict(click.get_current_context().params)
    return {
        "command": shlex.join(["metator", *sys.argv[1:]]),
        "parameters": {name: list(value) if isinstance(value, tuple) else value for name, value in params.items()},
    }


def hold_lock(held, path, in_place):
    """Hold the lock of the map at path until the end of the run (in held, an
    ExitStack): metacooler's (metacooler.lock), which every writer of the map
    takes too, and which the edit session that writes the map at the end of the
    run re-enters (see metator.metacool._staged_edit).

    In place, the lock is the file a symbolic link points to, as metacooler
    edits it; with --output-map, FILE's own name, so runs from one map, each
    with its own --output-map, run side by side."""
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


def hold_map(held, command, input_map, outfile, in_place, rounds=None):
    """Hold the lock of the map that a run (command) writes, outfile, until the
    end of the run (see hold_lock), then open its input map: the map as the run
    reads it, under that lock when the run updates it in place (in_place),
    which the map of the run is refused to be written from if it changed since
    (see metator.metacool.write_map's planned)."""
    hold_lock(held, outfile, in_place)
    planned = open_input_map(input_map)
    if in_place:
        refuse_a_read_only_map(outfile)
        warn_of_what_is_replaced(command, input_map, rounds)
    return planned


def refuse_a_read_only_map(path):
    """Refuse, before any work, the map that a run updates in place (path, its
    input map) when it has no write permission (judged by its mode too, as a
    run as root may write any file): metacooler refuses to replace a read-only
    map, which it would do at the end of the run."""
    mode = os.stat(path).st_mode
    if mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH) and os.access(path, os.W_OK):
        return
    message = (
        f"{path} is read-only (mode {stat.filemode(mode)}): the run would replace it with the map of its clusters. "
        "Make it writable (chmod u+w), or write the map of the run to another file with --output-map."
    )
    logger.error(message)
    raise mtmc.MetacoolError(message)


def run_map_path(input_map, output_map=None, force=False):
    """The map a run writes, and whether it updates its input map in place:
    output_map when given, else the file of its input map, the file that it
    points to when it is a symbolic link, as metacooler edits that file (a map
    given as file::/, the root of the file, is the file). output_map updates
    the input map in place when it names that file itself (see
    metator.metacool.is_input_map); any other output_map is a map of its own, a
    symbolic link there replaced, not its target.

    Raises
    ------
    MetacoolError
        If output_map names an entry of a collection or a directory, is the
        collection of the input map (an entry) or a hard link of the input map,
        or exists and is not the input map, without force; if the input map is
        an entry of a .mmetacool collection, without output_map; or if the
        directory of output_map cannot be made. It is raised before the run does
        any work, and its reason is logged as an error.
    """
    reason = None
    if output_map and "::" in str(output_map):
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
    """What the failure of a map says: the message of a MetacoolError, or, for
    another error, which is unexpected, its type and message."""
    unexpected = not isinstance(error, mtmc.MetacoolError)
    reason = f"{type(error).__name__}: {error}" if unexpected else str(error)
    return reason.rstrip("."), unexpected


def write_map(metacool, outfile, planned, what="the map of the clusters", **stores):
    """Write the map of a run at outfile (see run_map_path and step_map), from
    metacool, the input map of the run, which outfile is by default: what the
    run stores in it (stores: its network, its partition, its final clusters or
    its MGE clusters, see metator.metacool.write_map), in one edit session of
    metacooler.

    The map is refused if it changed since planned was opened, and outfile is
    only replaced once the new map is written and validated.

    Raises
    ------
    MetacoolError
        If the map cannot be written. The reason, and what is kept, is logged
        as an error first: the input map stays as it was, and so do the other
        outputs of the run.
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
        mtmc.write_map(metacool, outfile, force=True, planned=planned, **stores)
    except Exception as error:
        reason, unexpected = _failure(error)
        message = f"{outfile} is not written: {reason}. {kept}The other outputs of the run are kept."
        logger.error(message, exc_info=unexpected)
        raise mtmc.MetacoolError(message) from error


def run_network(metacool, planned, assembly, depth, enzyme, normalization, tmp_dir):
    """The network of the contigs of a network or pipeline run, from its map
    (metacool, as the run read it under its lock: planned): the contig data,
    measured on the assembly and read from the map (see contig_data_from_map),
    and the hits of the contigs of the map, counted on its contacts (see
    map_contacts)."""
    contig_data = mtn.create_contig_data(assembly, depth, enzyme)
    columns_from_map = contig_data_from_map(metacool, contig_data, normalization, enzyme, depth)
    contacts = map_contacts(metacool)
    names = planned.contigs["name"]
    table = mtn.contig_data_table(contig_data, names, mtn.contig_hits(contacts, len(names)), metacool)
    return mtmc.Network(table, normalization, network_entry(columns_from_map), tmp_dir), contacts


def stale_partition(planned, contacts, table, normalization):
    """Whether the input map of a network run, planned, holds a partition (the
    clustering rounds core and overlapping, which the run keeps, as it keeps
    every round) that was made on another network than the one the run stores
    (see metator.metacool.partition_network): its digest (see
    metator.network.network_digest), from the contacts of the map and the
    contig data of the run (table), with its normalization.

    Read before the run writes the map; metator validation refuses such a
    partition until metator partition runs again."""
    if not mtmc.holds_partition(planned):
        return False
    digest = mtn.network_digest(contacts, mtn.network_weights(table, normalization), table["Name"], normalization)
    return mtmc.partition_network(planned) != digest


def run_validation(
    metacool,
    assembly,
    table,
    network,
    partition,
    size,
    algorithm,
    iterations,
    overlap,
    res_param,
    threads,
    prefix,
    final_fasta_dir,
    tmp_dir,
    normalization,
):
    """The validation of the partition of a validation or pipeline run, from its
    map (metacool): the contig data of the contigs of the map (table, in its
    order), their network (network, by their positions in the map) and its
    partition (partition, a metator.clusters.Partition: the overlapping
    clusters that reach size are checked)."""
    partition_file = write_partition_table(table, partition, tmp_dir)
    # Decompressed once for the run, when it is compressed: the FASTA files
    # of the clusters are written from it at each step.
    assembly = mtv.plain_assembly(assembly, tmp_dir)
    overlapping_fasta_dir = join(tmp_dir, "overlapping_clusters")
    written = mtv.write_overlapping_fasta(metacool, assembly, partition, table, size, overlapping_fasta_dir, tmp_dir)
    contigs_data, summary = mtv.recursive_decontamination(
        algorithm,
        assembly,
        partition_file,
        final_fasta_dir,
        overlapping_fasta_dir,
        iterations,
        # The network of the map, by the positions of its contigs, as the
        # partition numbers them.
        network,
        tmp_dir,
        overlap,
        join(tmp_dir, "recursive_clusters"),
        res_param,
        size,
        threads,
        prefix,
        partition,
        metacool=metacool,
        written=written,
        temp_dir=tmp_dir,
    )
    final = contigs_data["Final_cluster"].astype(str).str.strip()
    binned = ~final.isin(mtmc.UNBINNED)
    assignment = dict(zip(contigs_data["Name"].astype(str)[binned], final[binned], strict=True))
    check_final_clusters(assignment, final_fasta_dir)
    return mtmc.FinalClusters(
        assignment=assignment,
        summary=join(tmp_dir, "cluster_summary.txt"),
        rounds=mtcl.read_rounds(tmp_dir),
        entry=validation_entry(normalization, table, summary),
    )


def check_final_clusters(assignment, final_fasta_dir):
    """Refuse the final clusters of a validation (assignment: the final cluster of
    each binned contig, by name) that are not those of the FASTA files it wrote
    in final_fasta_dir: a file of a cluster without a contig, or a cluster
    without a file.

    A name read as another (00123 as the number 123) would leave its cluster
    without contigs."""
    written = {name[: -len(".fa")] for name in os.listdir(final_fasta_dir) if name.endswith(".fa")}
    clustered = set(assignment.values())
    if written == clustered:
        return
    reasons = []
    if written - clustered:
        reasons.append(f"no contig in the final cluster(s) {mtmc._shown(sorted(written - clustered))}")
    if clustered - written:
        reasons.append(f"no FASTA file for the final cluster(s) {mtmc._shown(sorted(clustered - written))}")
    reason = (
        f"The final clusters of the run are not those of their FASTA files in {final_fasta_dir}: "
        f"{', and '.join(reasons)}. The map of the clusters is not written."
    )
    logger.error(reason)
    raise ValueError(reason)


def validation_entry(normalization, table, summary):
    """What a map records of the validation of its clusters, under the validation
    key of its metator metadata (see metator.metacool.write_map): the version
    of metaTOR, the command line and the parameters of the run, what it read of
    the map (read: the columns of the network of its contigs, from the contig
    data table, the weights of its normalization, and the rounds of its
    partition), and how the summary of its clusters, cluster_summary.txt, is
    written (summary, a metator.io.SummaryLayout), which its values in the
    genomes table of the map write again: its rows (summary_rows, the names of
    the clusters in order), its columns (summary_columns, in order) and their
    formats (summary_formats)."""
    return {
        "version": __version__,
        **run_metadata(),
        "normalization": normalization,
        "read": map_read(table, normalization, mtmc.PARTITION_ROUNDS),
        "summary_columns": list(summary.columns),
        "summary_rows": list(summary.rows),
        "summary_formats": dict(summary.formats),
    }


def map_read(table, normalization, rounds=()):
    """What a partition or validation run read of its map, as its entry of the
    metator metadata records it (read, see partition_entry and
    validation_entry): the columns of the contigs table that hold the network
    (contigs: those of table, the contig data of the network, that the map
    stores, see metator.metacool.network_columns), the weights of its
    normalization in the bins table (bins: none for "None"), and the clustering
    rounds (rounds: the partition, for a validation). pipeline records what its
    steps read, which it holds already."""
    weights = [] if normalization == "None" else [mtmc.normalization_weight_name(normalization)]
    read = {"contigs": mtmc.network_columns(table), "bins": weights}
    if rounds:
        read["rounds"] = list(rounds)
    return read


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


ASSEMBLY_HELP = "The assembly that the map was made from: a FASTA file, plain or compressed with gzip."
# A validation or pipeline run writes the clusters that miComplete checks
# from the assembly through metacooler (see
# metator.validation.write_overlapping_fasta), which decompresses a
# gzip-compressed one into the temporary directory of the run first, once.
CHECKED_ASSEMBLY_HELP = (
    ASSEMBLY_HELP + " A compressed one is decompressed once, into the temporary directory of the run, for the FASTA "
    "files of the clusters."
)
TMPDIR_TABLES_HELP = (
    "Temporary directory, where the run writes, in a directory of its own, the clusters that miComplete checks (at "
    "most the size of the assembly), miComplete's results and the tables of the validation, and, from a compressed "
    "assembly, its decompressed copy: give it room for both. Not in final_clusters, which the run makes anew, nor in "
    "a directory of the results of an earlier run (recursive_clusters, final_bin...), which the run refuses."
)
NO_CLEAN_UP_TABLES_HELP = (
    "Do not remove the temporary directory of the run, which the log names: it keeps miComplete's results, the FASTA "
    "files of the clusters it checked and the tables of the validation (contig_data_final.txt, cluster_summary.txt, "
    "binning.txt, clusters.txt, cluster_members.txt)."
)
SIZE_HELP = (
    "Threshold size to keep clusters in base pair: miComplete checks the overlapping clusters of the partition that "
    "reach it, and the clusters of the recursive partition larger than it."
)


def assembly_option(help):
    return _option("-a", "--assembly", required=True, metavar="FILE", panel="Input", help=help)


MAP_HELP = (
    "A .metacool contact map of the contigs of the assembly, made from the pairs of the run (by hicstuff pipeline -M "
    "metacool or metacooler create pairs), or an entry of a .mmetacool collection (file.mmetacool::entries/NAME): "
    "the run reads its network from it."
)


def map_argument(help=MAP_HELP, required=False):
    """The .metacool map of the contigs that a network or pipeline run starts
    from, its first argument (required); or that a partition or validation run
    updates (see step_map)."""
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
    """--output-map: the file that a run writes the map to, instead of the map of
    the run (see run_map_path)."""
    return _option("--output-map", metavar="FILE", help=help)


def network_options(normalization_short, assembly_help=ASSEMBLY_HELP):
    """Options of the network step: the assembly, what is measured on its contigs,
    and the normalization; the short alias of --normalization is -n in network,
    -m in pipeline."""
    return _options(
        assembly_option(assembly_help),
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


def partition_options(iterations, overlap, iterations_help, overlap_help, size_help=None):
    """Options of the Louvain or Leiden partition of the network, and, with
    size_help, its --size: the size of a cluster that the validation keeps
    (none for metator partition, which keeps every cluster)."""
    size = _option(
        "-s",
        "--size",
        type=click.INT,
        default=500000,
        show_default=True,
        metavar="INT",
        panel="Partition",
        help=size_help,
    )
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
        *([size] if size_help else []),
    )


def outdir_option(default=None, help="The output directory. Default: current directory."):
    return _option("-o", "--outdir", default=default, show_default=default is not None, metavar="DIR", help=help)


def prefix_option(what):
    return _option(
        "-p",
        "--prefix",
        metavar="STR",
        help=f"Name {what} STR_metator_... (by default, metator_...). Only letters, digits, '_' and '-'.",
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
@click.command_panel("Other commands", commands=["mge", "dump"])
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


@cli.command(short_help="Build the metaHiC network of the contigs of a .metacool map, and store it in the map.")
@map_argument(
    required=True,
    help=MAP_HELP + " The network of the run is stored in it, in place, or in --output-map, which an entry needs.",
)
@network_options("-n")
@outdir_option(help="The output directory, of the log of the run. Default: current directory.")
@force_option("Replace an existing --output-map.")
@tmpdir_option(
    help="Temporary directory, where the run writes, in a directory of its own, the contig data that it stores in the map (contig_data_network.txt)."
)
@no_clean_up_option(help="Keep the temporary directory of the run, which the log names, and its contig_data_network.txt.")
@output_map_option(
    help="Store the network in a copy of MAP written to FILE (replaced only with -F), and leave MAP as it was. By "
    "default, the run updates MAP in place.",
)
def network(
    map,  # the argument MAP
    assembly,
    depth,
    enzyme,
    normalization,
    outdir,
    force,
    tmpdir,
    no_clean_up,
    output_map,
):
    """Build the metaHiC network of the contigs, and store it in the map.

    [cyan]metator network MAP -a ASSEMBLY[/] builds the network from a
    .metacool contact map of the contigs (metacooler merge sums the maps of
    several libraries): that of the pairs it was made of. hicstuff
    ([cyan]hicstuff pipeline -M metacool[/]) or metacooler
    ([cyan]metacooler create pairs[/]) makes the map from the reads or the
    pairs.

    The Restriction_site and Shotgun_coverage columns of the map, when they
    give every contig a value, give the restriction sites and the coverage
    of the contigs that --enzyme and --depth would measure; given, these
    options replace them.

    In the network, contigs are the nodes and the edges are the contact
    counts between two contigs: every pair counts, as in a contact map, but
    the pairs within a contig. The run stores it in the map, for
    [cyan]metator partition MAP[/] and [cyan]metator validation MAP[/] to
    read: what it measured on each contig (GC_content, Hit, and
    Shotgun_coverage and Restriction_site when it has them) in its contigs
    table, and the normalization of the network as the weights
    metator_<normalization> of its bins, which replace those of an earlier
    network. The clusters of the map are kept: the partition of an earlier
    network too, with a warning, as [cyan]metator validation MAP[/] refuses
    it until [cyan]metator partition MAP[/] partitions the new network. The
    map is updated in place (or written to --output-map), and locked while
    the run writes it. [cyan]metator dump[/] writes the network, as
    network.txt, from the map of the clusters of a run; from a map that only
    metator network has run on, which dump refuses, [cyan]metacooler network
    MAP --level contig --balance metator_<normalization>[/] writes it, with a
    row per contig for the contacts within it, which are no edge.
    """
    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir, ExitStack() as held:
        # Defined the output directory.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)

        # Enable file logging
        log_file = log_file_name(outdir, "network")
        generate_log_header(log_file)

        # The map, the map that the run writes (the map itself, by default),
        # the normalization against its columns, and the assembly, checked
        # before any work. The run holds the lock of the map it writes from
        # now on: a second run on the same map is refused. It reads its map
        # under that lock (planned): the map its network is written from,
        # which is refused if it changed since.
        run_input(map)
        metacool_file, in_place = run_map_path(map, output_map, force)
        refuse_an_input_as_the_map(output_map, [assembly, depth])
        check_normalization(normalization, enzyme, depth, map)
        planned = hold_map(held, "network", map, metacool_file, in_place)
        check_assembly(assembly)

        # Print information of the workflow:
        logger.info("Enzyme: %s", enzyme)
        logger.info("Normalization: %s", normalization)

        # The contig data, measured on the assembly and read from the map,
        # the hits of the contigs of the map, counted on its contacts, then
        # the network, stored in the map.
        network, contacts = run_network(map, planned, assembly, depth, enzyme, normalization, tmp_dir)
        # A partition of another network, which the map keeps: said once the
        # map is written.
        stale = stale_partition(planned, contacts, network.table, normalization)
        write_map(map, metacool_file, planned, "the map with the network of its contigs", network=network)
        if stale:
            logger.warning(
                f"{metacool_file} holds a partition of another network than the one stored in it now, which the run "
                f"kept: metator validation refuses it until metator partition {metacool_file} runs again, which "
                "partitions this network."
            )

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# partition
# ---------------------------------------------------------------------------


@cli.command(short_help="Partition the metaHiC network of a .metacool map with the Louvain or Leiden algorithm.")
@map_argument(
    required=True,
    help="The .metacool contact map of the contigs of the run, which metator network stored the network in, or an "
    "entry of a .mmetacool collection (file.mmetacool::entries/NAME): the core and overlapping clusters of the run are "
    "stored in its clustering rounds, in place, or in --output-map, which an entry needs.",
)
@partition_options(
    iterations=100,
    overlap=80,
    iterations_help="Number of iterations of Louvain or Leiden.",
    overlap_help="Hamming distance threshold to use to merge clusters (percentage).",
)
@outdir_option(default="./", help="The output directory, of the log of the run.")
@prefix_option("the overlapping clusters")
@force_option("Replace an existing --output-map.")
@tmpdir_option(help="Temporary directory: partition writes no temporary file in it.")
@threads_option(help="Number of parallel threads allocated for the partition.")
@no_clean_up_option(help="Keep the temporary directory of the run: partition leaves nothing in it.")
@output_map_option(
    help="Store the partition of the run in a copy of MAP written to FILE (replaced only with -F), and leave MAP as it "
    "was. By default, the run updates MAP in place.",
)
def partition(
    map,  # the argument MAP
    algorithm,
    iterations,
    overlap,
    res_param,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
    output_map,
):
    """Partition the network using the Louvain or Leiden algorithm.

    Partition the network of the contigs that [cyan]metator network MAP[/]
    stored in the map iteratively with the Louvain or Leiden algorithm. Then
    looks for 'core' clusters that are constituted by the contigs which are
    always in the same community at each iterations. Then using the Hamming
    distance from these core clusters, merge the clusters with a small
    Hamming distance (close to 1), bigger than the percentage (overlapping
    parameter) given.

    Furthermore, both Leiden and Louvain algorithm are available here.
    However, the benchmark made show that here the Louvain algorithm have
    better performance and is faster on seawater and gut metagenomic
    samples.

    [cyan]metator partition MAP ...[/] runs on the map of the contigs of a
    run step by step ([cyan]metator network MAP[/], partition, validation):
    it stores every core and overlapping cluster, with its contigs, as the
    rounds core and overlapping of the clustering rounds of the map (its
    /clusters group), in place of those it held, and updates the map in
    place (or writes --output-map). Its genomes and contigs tables are left
    as they were: [cyan]metator validation MAP[/] checks the overlapping
    clusters that reach its --size, and writes the final clusters there.
    """
    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up), ExitStack() as held:
        # Defined the output directory.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)

        # Enable file logging, the map of the run checked and locked, and
        # the network that metator network stored in it read, with the data
        # of its contigs (see start_step).
        log_file, metacool_file, bins_prefix, planned, stored, _ = start_step(
            "partition", outdir, [], map, output_map, force, held, prefix
        )

        # Partition the network
        clusters = mtp.partition(algorithm, stored.network, iterations, overlap / 100, res_param, threads, bins_prefix)

        # The core and overlapping clusters, stored in the map, their
        # contigs by name: the IDs of the network are their positions in it.
        write_map(
            map,
            metacool_file,
            planned,
            "the map with its partition",
            partition=partition_entry(stored.normalization, stored.digest, stored.table),
            rounds=mtcl.partition_rounds(clusters, list(stored.table["Name"])),
            round_names=list(mtmc.PARTITION_ROUNDS),
        )

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


@cli.command(short_help="Validate the clusters with miComplete and decontaminate them with a recursive partition.")
@map_argument(
    required=True,
    help="The .metacool contact map of the contigs of the run (that of network MAP and partition MAP), which metator "
    "network stored the network in, and metator partition its partition, or an entry of a .mmetacool collection "
    "(file.mmetacool::entries/NAME): it becomes the map of the clusters of the run, as pipeline MAP writes it, in place, "
    "or in --output-map, which an entry needs.",
)
@assembly_option(CHECKED_ASSEMBLY_HELP)
@partition_options(
    iterations=10,
    overlap=90,
    iterations_help="Number of recursive iterations of Louvain or Leiden.",
    overlap_help="Hamming distance threshold to use to merge clusters (percentage).",
    size_help=SIZE_HELP,
)
@outdir_option(default="./", help="Path to the directory to write the output.")
@prefix_option("the recursive clusters")
@force_option(
    "Replace the final clusters of an earlier run in the output directory (final_clusters), and an existing "
    "--output-map. The other results of an earlier run in the output directory (its tables, recursive_clusters, "
    "1.x's outputs) are refused, with or without -F: they are never removed."
)
@tmpdir_option(help=TMPDIR_TABLES_HELP)
@threads_option(help="Number of parallel threads, for the partition and miComplete.")
@no_clean_up_option(help=NO_CLEAN_UP_TABLES_HELP)
@output_map_option(
    help="Write the contact map of the clusters to FILE (replaced only with -F), and leave MAP as it was. By default, "
    "the run updates MAP in place.",
)
def validation(
    map,  # the argument MAP
    assembly,
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

    The clusters to validate are the overlapping clusters of the partition
    that [cyan]metator partition MAP[/] stored in the map that reach
    --size: the run writes their FASTA files from the assembly, in its
    temporary directory, for miComplete to check.

    Recursively decontaminate the contaminated clusters evaluated by
    miComplete. The command writes the FASTA files of the final clusters in
    the output directory (final_clusters), and stores the clusters in the
    map: the tables of the validation stay in its temporary directory, which
    -N keeps, and [cyan]metator dump MAP[/] writes them from the map.

    Only clusters with a weighted completeness of at least 33% and a
    contamination, (weighted redundancy - 1) / weighted completeness, of at
    least 5% are subject to the recursive step. Its clusters replace the
    original cluster when the most complete of them keeps enough of its
    completeness; otherwise the original cluster is kept.

    [cyan]metator validation MAP ...[/] runs on the map of the contigs of a
    run step by step ([cyan]metator network MAP[/], partition MAP,
    validation), whose network it partitions again: at its end, it writes
    the map of the clusters of the run, its contigs grouped by final cluster
    with the summary of the clusters and all its clustering rounds, as
    [cyan]metator pipeline MAP[/] does, in place of the map (or to
    --output-map).
    """
    # Launch miComplete to evaluate the completion and the contamination. If
    # asked rerun Louvain to try to reduce the contamination, rerun miComplete
    # if the contamination decrease without a huge decrease of the completion
    # keep the new bins. Otherwise go back to the old state.

    # Private temporary directory inside --tmpdir, out of the directories of
    # the clusters, and the lock of the map the run writes, held until its
    # end.
    refuse_tmpdir_in_the_clusters(tmpdir, outdir)
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir, ExitStack() as held:
        # Defined the output directory and output file names.
        outdir = outdir or "."
        os.makedirs(outdir, exist_ok=True)
        final_fasta_dir = join(outdir, "final_clusters")

        # Enable file logging, the map of the run checked and locked, the
        # network that metator network stored in it and the partition that
        # metator partition stored in it read, the assembly checked, and the
        # directory of the final clusters made (see start_step).
        log_file, metacool_file, bins_prefix, planned, stored, clusters = start_step(
            "validation",
            outdir,
            [final_fasta_dir],
            map,
            output_map,
            force,
            held,
            prefix,
            [assembly],
            assembly,
            partition=True,
        )

        # The validation of the partition, in the temporary directory of the
        # run, from the network that metator network stored in the map, by
        # the positions of its contigs, as the partition numbers them.
        final = run_validation(
            map,
            assembly,
            stored.table,
            stored.network,
            clusters,
            size,
            algorithm,
            iterations,
            overlap / 100,
            res_param,
            threads,
            bins_prefix,
            final_fasta_dir,
            tmp_dir,
            stored.normalization,
        )

        # The map of the clusters, written last, as pipeline MAP writes it:
        # the network that the map stores, which it keeps, comes after the
        # columns and the weights that the regroup writes (see
        # metator.metacool.write_map).
        write_map(
            map,
            metacool_file,
            planned,
            clusters=final,
            network=mtmc.Network(stored.table, stored.normalization, entry=None, tmp_dir=tmp_dir),
        )

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------


@cli.command(short_help="Run network, partition and validation in one command line.")
@map_argument(required=True)
@network_options("-m", CHECKED_ASSEMBLY_HELP)
@partition_options(
    iterations=100,
    overlap=80,
    iterations_help="Number of iterations of Louvain or Leiden for the partition step.",
    overlap_help="Hamming distance threshold to use to merge clusters for the first step (percentage).",
    size_help=SIZE_HELP,
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
@prefix_option("the clusters")
@force_option(
    "Replace the final clusters of an earlier run in the output directory (final_clusters), and an existing "
    "--output-map. The other results of an earlier run in the output directory (its tables, recursive_clusters, "
    "1.x's outputs) are refused, with or without -F: they are never removed."
)
@tmpdir_option(help=TMPDIR_TABLES_HELP)
@threads_option(help="Number of parallel threads allocated for the partition and the validation.")
@no_clean_up_option(help=NO_CLEAN_UP_TABLES_HELP)
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

    The run does what [cyan]metator network MAP[/], partition MAP and
    validation MAP do one by one, and writes the map that they write, once,
    at its end: the FASTA files of the final clusters go to the output
    directory (final_clusters), and the tables of the validation stay in
    its temporary directory, which -N keeps.
    """
    # Private temporary directory inside --tmpdir, out of the directories of
    # the clusters, and the lock of the map the run writes, held until its
    # end.
    refuse_tmpdir_in_the_clusters(tmpdir, outdir)
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
        metacool_file, in_place = run_map_path(map, output_map, force)
        refuse_an_input_as_the_map(output_map, [assembly, depth])
        check_normalization(normalization, enzyme, depth, map)
        bins_prefix = bin_prefix(prefix)
        # The run holds the map's lock from now on: a second run on the same
        # map is refused before any work. It reads its map under that lock
        # (planned).
        planned = hold_map(held, "pipeline", map, metacool_file, in_place)

        # Print information of the workflow:
        logger.info(f"Partition algorithm: {algorithm}")
        logger.info(f"Partition iterations: {iterations}")
        logger.info(f"Overlapping parameter: {overlapping_parameter}")
        logger.info(f"Recursive partition iterations: {rec_iter}")
        logger.info(f"Recursive overlapping parameter: {recursive_overlapping_parameter}\n")

        # The final clusters of an earlier run, and what a run of an earlier
        # version wrote in the output directory: refused now without
        # --force, before the assembly is read, and replaced or removed only
        # once it is checked.
        final_fasta_dir = join(outdir, "final_clusters")
        refuse_output_dir(final_fasta_dir, force)
        refuse_earlier_outputs(outdir, keep=(map, output_map))

        # The network, as metator network makes it: the contig data,
        # measured on the assembly and read from the map, before the clusters
        # of an earlier run are replaced (a map whose columns cannot give what
        # the normalization needs is refused once the assembly is read), and
        # the contacts of its contigs, weighed with the normalization.
        check_assembly(assembly)
        network, contacts = run_network(map, planned, assembly, depth, enzyme, normalization, tmp_dir)
        make_output_dir(final_fasta_dir, force)
        edges = mtn.map_network(map, mtn.network_weights(network.table, normalization), contacts)

        # Partition the network, as metator partition does.
        clusters = mtp.partition(algorithm, edges, iterations, overlapping_parameter, res_param, threads, bins_prefix)

        # Validate the partition, as metator validation does.
        final = run_validation(
            map,
            assembly,
            network.table,
            edges,
            clusters,
            size,
            algorithm,
            rec_iter,
            recursive_overlapping_parameter,
            res_param,
            threads,
            bins_prefix,
            final_fasta_dir,
            tmp_dir,
            normalization,
        )

        # The map of the clusters, written last, once: what the three steps
        # store, as they store it one by one (see metator.metacool.write_map).
        write_map(
            map,
            metacool_file,
            planned,
            network=network,
            partition=partition_entry(
                normalization,
                mtn.network_digest(
                    contacts, mtn.network_weights(network.table, normalization), network.table["Name"], normalization
                ),
                network.table,
            ),
            clusters=final,
        )

    generate_log_footer(log_file)


# ---------------------------------------------------------------------------
# mge
# ---------------------------------------------------------------------------


@cli.command(short_help="Cluster the MGE contigs of a .metacool map into MGE MAGs, and store them in the map.")
@map_argument(
    required=True,
    help="The .metacool contact map of the run, which metator network stored the network of its contigs in (that of "
    "network MAP, pipeline MAP or validation MAP), or an entry of a .mmetacool collection "
    "(file.mmetacool::entries/NAME): the contacts of the MGE contigs come from it, and the MGE clusters are stored in "
    "its clustering rounds, in place, or in --output-map, which an entry needs. It must hold the contact statistics of "
    "the reads of its contigs at 1001 bp (made from the pairs by hicstuff pipeline -M metacool or metacooler create "
    "pairs, with metacooler 0.4.0 or later).",
)
@assembly_option(ASSEMBLY_HELP)
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
    panel="Clustering",
    help="Threshold to use for clustering.",
)
@_option(
    "-r",
    "--random",
    is_flag=True,
    panel="Clustering",
    help="Make a random clustering of the same shape, to compare the MGE clusters with: stored as the round "
    "mge_random, which leaves the round mge as it was, and written to mges_binned_random.fa.",
)
@outdir_option(help="Path to the output directory where the output will be written. Default: current directory.")
@force_option(
    "Replace an existing --output-map. A summary of the MGE clusters that an earlier version wrote in the output "
    "directory (mges_bin_summary.tsv or mges_cluster_summary.tsv) is refused, with or without -F: it is never removed."
)
@tmpdir_option(
    help="Path to temporary directory: where metacooler copies a compressed assembly, plain, while it reads it. mge "
    "writes no temporary file of its own."
)
@no_clean_up_option(help="Keep the temporary directory of the run: mge leaves nothing in it.")
@output_map_option(
    help="Store the MGE clusters in a copy of MAP written to FILE (replaced only with -F), and leave MAP as it was. By "
    "default, the run updates MAP in place.",
)
def mge(
    map,  # the argument MAP
    assembly,
    mges,
    threshold_bin,
    random,
    outdir,
    force,
    tmpdir,
    no_clean_up,
    output_map,
):
    """Cluster MGE contigs.

    Cluster the MGE contigs of a map into MGE MAGs, from the Hi-C contacts
    between them: those that the map holds between two contigs, and those
    within each contig between reads more than 1 kb apart, which the
    contact statistics of the map count (its QC block, made from the pairs
    at 1001 bp). The GC content and hits of each contig are those that
    metator network stored in the map.

    Writes the MGE fasta, mges_binned.fa: one record per MGE MAG, its contigs
    joined by 200 "N".

    The MGE clusters are stored in the map, MAP: as the round mge of its
    clustering rounds (its /clusters group), with their contigs and their
    length, GC content, hits and clustering score, in place of an earlier
    round mge, after the other rounds, which are kept (--random stores the
    round mge_random, and leaves the round mge as it was). Its genomes and
    contigs tables are left as they were. MAP is updated in place (or
    written to --output-map), and locked while the run writes it. A run
    that replaces the clusters of the map (partition, validation, pipeline)
    drops the rounds of mge, with a warning: run mge again then.
    """
    # Generate log
    log_file = log_file_name("", "mge")
    generate_log_header(log_file)

    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir, ExitStack() as held:
        # The map that the run stores its MGE clusters in, checked and
        # locked.
        name = mtcl.MGE_RANDOM if random else mtcl.MGE
        metacool_file, planned = step_map("mge", map, output_map, force, held, inputs=[assembly, mges], rounds=[name])
        # The assembly that metacooler writes mges_binned.fa from, refused
        # as the network, the validation and the pipeline refuse it.
        check_assembly(assembly)
        # The summary of an earlier version, which would pass for this run's.
        outdir = outdir or "."
        refuse_earlier_outputs(outdir, MGE_EARLIER_OUTPUTS, keep=(map, output_map))

        # The MGE contigs, and their clusters from the contacts of the map:
        # what the map cannot give is refused before anything is written.
        mges_list = mio.import_mges_contigs(mges)
        try:
            summary, mge_bins = mtm.mge_binning(map, mges_list, threshold_bin=threshold_bin, random=random)
            # The output directory, and one record per MGE cluster.
            os.makedirs(outdir, exist_ok=True)
            fasta = join(outdir, "mges_binned_random.fa" if random else "mges_binned.fa")
            mtm.write_mges_fasta(map, assembly, mge_bins, fasta, temp_dir=tmp_dir)
        except mtmc.MetacoolError as error:
            logger.error(str(error))
            raise

        # The MGE clusters, stored in the map: the round mge (or mge_random),
        # in place of the round of that name that the map held, after the
        # other rounds, which are kept as they were, and so is everything
        # else; the run is recorded with the round, and under its name in the
        # metator metadata of the map (see metator.metacool.write_map).
        write_map(
            map, metacool_file, planned, "the map with its MGE clusters", mge=(mtcl.mge_round(summary, name), mge_entry(name))
        )

    generate_log_footer(log_file)


def mge_entry(name):
    """What a map records of the mge run that stores its round name in it (mge, or
    mge_random), with the round and under the name of the round in its metator
    metadata (see metator.metacool.write_map): the tool, the version of
    metaTOR, the command line and the parameters of the run, and what it read
    of the map (read): the columns of its contigs table that metator network
    stored (metator.mge.MGE_COLUMNS), and the contact statistics of its reads,
    the counts of each contig (its qc_intra_contig_long column) and the QC
    block of its metadata."""
    read = {"contigs": [*mtm.MGE_COLUMNS, mtm.INTRA_LONG], "metadata": ["qc"]}
    return {"name": name, "tool": "metaTOR", "version": __version__, "step": "mge", **run_metadata(), "read": read}


# ---------------------------------------------------------------------------
# dump
# ---------------------------------------------------------------------------


@cli.command(short_help="Write the tables of a run from the map of its clusters.")
@map_argument(
    required=True,
    help="The map of the clusters of a run, which metator validation or metator pipeline wrote: a .metacool file, or "
    "an entry of a .mmetacool collection (file.mmetacool::entries/NAME). It is read only.",
)
@outdir_option(help="The output directory of the tables, made if it is not there. Default: current directory.")
@force_option("Replace the tables of an earlier dump that the output directory holds.")
def dump(
    map,  # the argument MAP
    outdir,
    force,
):
    """Write the tables of a run from the map of its clusters.

    [cyan]metator dump MAP -o DIR[/] writes, from the map that a run wrote
    ([cyan]metator validation MAP[/], or [cyan]metator pipeline MAP[/]), the
    tables of its validation, which a run keeps in its temporary directory
    only (with -N), and its network:

    - contig_data_final.txt: one row per contig of the map, in its order,
      with its contig data, its core, overlapping and recursive clusters, and
      its final cluster (no ID column);
    - cluster_summary.txt: the summary of the final clusters;
    - binning.txt: each contig of a final cluster, and its cluster;
    - clusters.txt and cluster_members.txt: every clustering round of the
      map, the MGE clusters of metator mge included, and their contigs;
    - network.txt: the network of the contigs, by contig name, weighed with
      the normalization that metator network stored.

    The map is read only, and the tables are of one state of it: a map that
    another run replaces as it is read is refused (run the dump again). A
    map without the clusters of a validation (the validation entry of its
    metator metadata, and the clustering round final) is refused, and so is
    one that a development build of metaTOR 2.0 validated, one whose
    clusters were edited after the run (by metacooler reassign), one whose
    network metator network replaced after the run (run it again with the
    normalization of the run), and an output directory that holds tables of
    an earlier dump, without -F. The tables are replaced all at once: a dump
    that fails leaves them as they were.
    """
    outdir = outdir or "."
    os.makedirs(outdir, exist_ok=True)

    # Enable file logging
    log_file = log_file_name(outdir, "dump")
    generate_log_header(log_file)

    # The map, checked before anything is read.
    run_input(map)
    try:
        mtd.dump(map, outdir, force)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise

    generate_log_footer(log_file)
