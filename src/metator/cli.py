#!/usr/bin/env python3
# coding: utf-8

"""metaTOR command line interface.

One function per command: network, partition, validation, pipeline, mge,
host and dump. The ``metator`` entry point is :func:`cli`.
"""

import glob
import logging
import os
import random
import re
import shlex
import shutil
import stat
import sys
import time
from contextlib import ExitStack, contextmanager
from os.path import exists, join
from typing import cast

import numpy as np
import rich_click as click
from click.core import ParameterSource

import metator.clusters as mtcl
import metator.dump as mtd
import metator.figures as mtf
import metator.host as mth
import metator.inputs as mti
import metator.io as mio
import metator.log as mtl
import metator.metacool as mtmc
import metator.mge as mtm
import metator.network as mtn
import metator.outputs as mto
import metator.partition as mtp
import metator.report as mtr
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
        logger.debug(f"The temporary directory of the run, {tmp_dir}, is removed with its files: -N keeps them.")
        shutil.rmtree(tmp_dir)
        if base_created and not os.listdir(base):
            os.rmdir(base)
    else:
        logger.info(f"The temporary directory of the run, {tmp_dir}, and its files are kept (-N).")


# The directories that a validation, mge or pipeline run makes anew in its
# output directory (see metator.outputs), which -F replaces: the FASTA files of
# its final clusters (validation), of its MGE clusters (mge: fasta/mges, or
# fasta/mges_random with --random), and those of its MAGs, their cool files
# and its MGE clusters (pipeline, the last with --mges). Its plots, in plots/,
# -F replaces too, but none of the other files there (see
# refuse_earlier_plots).
VALIDATION_DIRS = (mto.MAGS_WITH_MGES_DIR,)
PIPELINE_DIRS = (mto.MAGS_WITH_MGES_DIR, mto.MAGS_DIR, mto.COOL_DIR, mto.MGES_DIR)

# An earlier run's results in the output directory of a validation or pipeline
# run (2.0 development builds' tables, recursive clusters, final clusters and
# plots, 1.x's outputs, the plots of 1.x that 2.0 no longer draws among them),
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
    "final_clusters",
    "bin_summary.txt",
    "network.txt",
    "clustering_matrix*.npz",
    "recursive_bin",
    "final_bin",
    "final_bin_unscaffold",
    "overlapping_bin",
    "scaffold_bin",
    "plot/bins_distribution.png",
    "plot/bins_size_distribution.png",
    "plot/MAGs_*_distribution.png",
    "plot/network_heatmap.png",
)

# An earlier mge run's summary (1.x's, a development build's), which the round
# mge now holds, and its FASTA file of the MGE clusters (1.x's, a development
# build's), which fasta/mges replaces: refused with or without -F, never
# removed.
MGE_EARLIER_OUTPUTS = ("mges_bin_summary.tsv", "mges_cluster_summary.tsv", "mges_binned.fa", "mges_binned_random.fa")

# An earlier host run's table and plot of the hosts (a development build's),
# which tables/mge_hosts.txt and plots/mge_hosts.pdf replace: refused with or
# without -F, never removed.
HOST_EARLIER_OUTPUTS = ("mge_hosts.txt", "mge_hosts.png")
# The directories that metator export can make anew in its output directory
# (see metator.outputs): the FASTA files (with --assembly; those of the MGE
# clusters when the map holds them), and the cool files.
EXPORT_DIRS = (mto.MAGS_WITH_MGES_DIR, mto.MAGS_DIR, mto.MGES_DIR, mto.MGES_RANDOM_DIR, mto.COOL_DIR)


def _within(path, directory):
    """Whether path is directory, or in it (both absolute)."""
    return os.path.commonpath([path, directory]) == directory


def _earlier_directories(outputs=EARLIER_OUTPUTS):
    """The directories among outputs (EARLIER_OUTPUTS by default)."""
    return tuple(name for name in outputs if "." not in name)


def refuse_tmpdir_in_the_clusters(tmpdir, outdir, dirs=VALIDATION_DIRS, earlier=EARLIER_OUTPUTS):
    """Refuse, as a usage error, before the run creates or removes anything,
    a --tmpdir in a directory of the output directory, outdir, that the run
    makes anew (dirs, which -F replaces, with the temporary directory of the
    run in it: see VALIDATION_DIRS and PIPELINE_DIRS) or refuses as the
    results of an earlier run (earlier, EARLIER_OUTPUTS by default: the run
    would make it, then refuse it). --tmpdir is in it as written, or once the symbolic links
    of both are resolved, but for that directory itself: -F replaces a
    symbolic link, not its target."""
    base = tmpdir or "./tmp"
    for name in (*dirs, *_earlier_directories(earlier)):
        directory = join(outdir or ".", name)
        as_written = _within(os.path.abspath(base), os.path.abspath(directory))
        resolved = _within(os.path.realpath(base), join(os.path.realpath(outdir or "."), name))
        if as_written or resolved:
            where = (
                "which the run makes anew (-F replaces an earlier one)"
                if name in dirs
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
    in outdir (outputs: EARLIER_OUTPUTS, MGE_EARLIER_OUTPUTS for mge,
    HOST_EARLIER_OUTPUTS for host), which
    this run neither writes nor replaces: the reason, which names them, is logged
    as an error, then ValueError is raised, and nothing is removed. A symbolic
    link counts, even one to nothing; the maps of the run (keep) never do, nor a
    directory that holds one. The reason names them by their path in outdir
    (plot/network_heatmap.png)."""
    kept = [mtmc.map_file(path) if path else None for path in keep]
    found = sorted(
        path
        for path in {path for pattern in outputs for path in glob.glob(join(glob.escape(outdir), pattern))}
        if not _is_or_holds(path, kept)
    )
    if found:
        logger.error(
            f"{outdir} holds results of an earlier run, which this run does not write and does not replace: "
            f"{mtmc._shown([os.path.relpath(path, outdir) for path in found])}. Move them away, or give the run "
            "another output directory with -o."
        )
        raise ValueError


def refuse_earlier_plots(outdir, force):
    """Refuse, before any work, the plots of an earlier run in the output
    directory of a validation or pipeline run, outdir (see
    metator.figures.plot_paths), without force, which replaces them: a symbolic
    link too, even one to nothing. A plots/ that is no directory, and a
    directory where a plot goes, are refused with or without force: the run
    would draw no plot there. The reason is logged as an error, then
    ValueError is raised. The other files of plots/ are left as they are (the
    plots of an earlier version in plot/, see EARLIER_OUTPUTS, are refused).

    Returns
    -------
    list of str
        The plots of an earlier run, which the run removes once it is checked
        (see remove_earlier_plots)."""
    plot_dir = join(outdir, mtf.PLOT_DIR)
    if os.path.lexists(plot_dir) and not os.path.isdir(plot_dir):
        logger.error(
            f"{plot_dir} is not a directory, where the run draws its plots: move it away, or give the run another "
            "output directory with -o."
        )
        raise ValueError
    found = [path for path in mtf.plot_paths(outdir) if os.path.lexists(path)]
    for path in found:
        if os.path.isdir(path) and not os.path.islink(path):
            logger.error(
                f"{path} is a directory, where the run draws a plot: move it away, or give the run another output "
                "directory with -o."
            )
            raise ValueError
    if found and not force:
        logger.error(
            f"{plot_dir} holds the plots of an earlier run ({mtmc._shown([os.path.basename(path) for path in found])}), "
            "which this run draws again: give -F to replace them, or give the run another output directory with -o."
        )
        raise ValueError
    return found


def remove_earlier_plots(plots):
    """Remove the plots of an earlier run (see refuse_earlier_plots), which -F
    replaces: a symbolic link, not its target. The other files of plots/ are
    left as they are."""
    for path in plots:
        logger.debug(f"Removing {path}, a plot of an earlier run.")
        os.remove(path)


def refuse_no_directory(path, what):
    """Refuse, before any work, a path of the output directory of a run that is
    no directory (a file, a symbolic link to nothing), where the run makes one
    to write what (for the message: the FASTA files of its clusters, its
    tables...), with or without -F: the run could write nothing there. The
    reason is logged as an error, then ValueError is raised."""
    if os.path.lexists(path) and not os.path.isdir(path):
        logger.error(
            f"{path} is not a directory, where the run writes {what}: move it away, or give the run another output "
            "directory with -o."
        )
        raise ValueError


def refuse_earlier_tables(tables_dir, force):
    """Refuse, before any work, the tables of an earlier run in tables_dir, the
    tables/ of the output directory of a pipeline run, which it writes again
    (those of metator dump, see metator.dump.TABLES), without force, which
    replaces them: a symbolic link too, even one to nothing; and a tables_dir
    that is no directory, with or without force (see refuse_no_directory). The
    reason is logged as an error, then ValueError is raised. The other files of
    tables/ are left as they are.

    Returns
    -------
    list of str
        The tables of an earlier run, which the run removes once it is checked
        (see remove_earlier_tables): they would describe another map than the
        one that it writes, should it stop before it writes its own."""
    refuse_no_directory(tables_dir, "its tables")
    found = [join(tables_dir, name) for name in mtd.existing_tables(tables_dir)]
    if found and not force:
        logger.error(
            f"{tables_dir} holds the tables of an earlier run ({mtmc._shown([os.path.basename(path) for path in found])}), "
            "which this run writes again: give -F to replace them, or give the run another output directory with -o."
        )
        raise ValueError
    return found


def remove_earlier_tables(tables):
    """Remove the tables of an earlier run (see refuse_earlier_tables), which -F
    replaces: a symbolic link, not its target."""
    for path in tables:
        logger.debug(f"Removing {path}, a table of an earlier run.")
        os.remove(path)


def refuse_earlier_report(outdir, name, force):
    """Refuse, before any work, the report of an earlier run in the output
    directory, outdir, of a pipeline, export or report run named name, which it
    writes again (see metator.report.report_files), without force, which
    replaces it: a symbolic link too, even one to nothing. A directory where a
    file of the report goes is refused with or without force: the run could
    write nothing there. The reason is logged as an error, then ValueError is
    raised.

    Returns
    -------
    list of str
        The files of the report of an earlier run, which the run removes once
        it is checked (see remove_earlier_report): they would describe another
        map than the one that it writes, should it stop before it writes its
        own."""
    paths = [join(outdir, file) for file in mtr.report_files(name)]
    for path in paths:
        if os.path.isdir(path) and not os.path.islink(path):
            logger.error(
                f"{path} is a directory, where the run writes its report: move it away, or give the run another output "
                "directory with -o."
            )
            raise ValueError
    found = [path for path in paths if os.path.lexists(path)]
    if found and not force:
        logger.error(
            f"{outdir} holds the report of an earlier run ({mtmc._shown([os.path.basename(path) for path in found])}), "
            "which this run writes again: give -F to replace it, or give the run another output directory with -o."
        )
        raise ValueError
    return found


def remove_earlier_report(reports):
    """Remove the report of an earlier run (see refuse_earlier_report), which -F
    replaces: a symbolic link, not its target."""
    for path in reports:
        logger.debug(f"Removing {path}, the report of an earlier run.")
        os.remove(path)


def refuse_output_dir(path, force):
    """Refuse an existing output directory of the run, unless --force replaces it: a
    symbolic link too, even one to nothing."""
    if os.path.lexists(path) and not force:
        logger.error(
            "%s already existed. Remove directory or use -F argument to overwrite it.",
            path,
        )
        raise ValueError


def make_outdir(outdir):
    """Make the output directory of a run, outdir, if it is not there: one that
    cannot be made (a file there, a directory that is not writable) is
    refused before any work, and the reason is logged as an error."""
    try:
        os.makedirs(outdir, exist_ok=True)
    except OSError as error:
        message = f"The output directory {outdir} cannot be made ({error.strerror}): give the run another one with -o."
        logger.error(message)
        raise mtmc.MetacoolError(message) from error


def make_output_dir(path, force):
    """Create an output directory of the run (the FASTA files of its clusters,
    the cool files of its MAGs): an existing one is replaced with --force, or
    refused; a symbolic link or a file there is replaced, not its target."""
    refuse_output_dir(path, force)
    if os.path.islink(path) or os.path.isfile(path):
        os.remove(path)
    elif exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def check_normalization(normalization, enzyme, depth, metacool, read=None):
    """The normalizations that need the restriction sites or the shotgun coverage
    of the contigs: measured on the assembly with --enzyme and --depth or read
    from the Restriction_site and Shotgun_coverage columns of the contigs table
    of the map of the run, metacool (see metator.network.contig_data_from_map);
    None for the map that a run makes from pairs or BAM files, which has no
    such column."""
    needs = [
        (("RS", "theoritical_hit"), enzyme, "enzyme", "Restriction_site"),
        (("abundance", "theoritical_hit"), depth, "depth", "Shotgun_coverage"),
    ]
    for normalizations, given, option, column in needs:
        if normalization not in normalizations or given:
            continue
        message = f'For "{normalizations[0]}" and "{normalizations[1]}" normalization, {option} is required'
        if metacool is None:
            what = mtn.MAP_COLUMNS[column][2]
            message += f", as the map that the run makes from its pairs or BAM files gives no {what}."
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
            f"{input_map} is not a .metacool contact map: give the map of the contigs (made by hicstuff pipeline -M "
            "metacool or metacooler create pairs), or the pairs files or the BAM files of the run, which the run makes it "
            "from."
        )
        raise ValueError
    refuse_a_map_named_as_a_uri(input_map, path)
    open_input_map(input_map)


def run_inputs(arguments, min_quality_given=False):
    """The inputs of a network or pipeline run, its arguments (see
    metator.inputs.resolve): the .metacool map of the contigs, or the pairs
    files or the BAM files that the run makes it from. Checked before any work:
    the reason of a refusal is logged as an error, then raised (a
    MetacoolError). --min-quality (min_quality_given: given on the command
    line) filters the alignments of BAM files only, and is refused with other
    inputs."""
    try:
        inputs = mti.resolve(arguments)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    if inputs.kind != mti.BAM and min_quality_given:
        message = (
            f"--min-quality filters the alignments of BAM files, and the inputs of the run are "
            f"{'a map' if inputs.kind == mti.MAP else inputs.describe()}: leave out --min-quality."
        )
        logger.error(message)
        raise mtmc.MetacoolError(message)
    return inputs


def given_on_the_command_line(name):
    """Whether the parameter name of the running command was given on its
    command line, rather than left to its default."""
    ctx = cast(click.RichContext, click.get_current_context())
    return ctx.get_parameter_source(name) == ParameterSource.COMMANDLINE


def made_map_path(inputs, name, outdir=".", inplace=False, force=False):
    """The map that a run makes from its pairs or BAM files (inputs, see
    metator.inputs.make_map): name (see metator.inputs.map_name), in its output
    directory, outdir. Checked before any work, and its directory made: the
    reason of a refusal is logged as an error, then raised (a MetacoolError).

    It is refused with --inplace (inplace), which updates the map that a run
    starts from; when it holds :: (a URI, file::group, to metacooler), or is a
    directory; and when it exists, without force. With force, the map of an
    earlier run there is replaced once the new map is written (a symbolic link
    there is replaced, not its target)."""
    outfile = name if os.path.normpath(outdir or ".") == "." else join(outdir, name)
    reason = None
    if inplace:
        reason = (
            f"--inplace updates the map that a run starts from, and this run starts from {inputs.describe()}, which it "
            f"makes its map from, in its output directory ({outfile}): leave out --inplace."
        )
    elif "::" in outfile:
        reason = (
            f"{outfile}, the map that the run makes in its output directory, holds ::, which metacooler reads as "
            "file::group: give the run another output directory with -o."
        )
    elif os.path.isdir(outfile):
        reason = (
            f"{outfile} is a directory, where the run writes the map that it makes: move it away, or give the run another "
            "output directory with -o."
        )
    elif os.path.lexists(outfile) and not force:
        reason = (
            f"{outfile} exists, where the run writes the map that it makes from {inputs.describe()}: give -F to replace "
            "it, or another output directory with -o; or start the run from that map."
        )
    if reason:
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    make_outdir(os.path.dirname(outfile) or ".")
    logger.info(f"The run makes its map, {outfile}, from {inputs.describe()}, and updates it in place.")
    return outfile


def run_map(command, inputs, held, name, assembly, depth, enzyme, normalization, outdir, inplace, force, others=(), rounds=None):
    """The map of a network or pipeline run (command), checked before any work,
    as (map, outfile, planned): the map that the run starts from, the map that
    it writes (see run_map_path), and its map as the run read it, under the
    lock of the map that it writes, which the run holds from now on (see
    hold_map).

    From pairs or BAM files (inputs), the run makes its map first, at outfile
    (see made_map_path; name: its name), which it then updates in place: it
    holds the lock of outfile, and planned is None until the map is made (see
    make_run_map). The parameters that the map stores give the map that the run
    started from (map), whatever its inputs.

    others: the other files that the run reads (the list of the MGE contigs of
    pipeline --mges), which the map that it writes is refused to be, as its
    inputs; rounds: the clustering rounds that the run replaces in the map
    that it updates in place, which it warns of (see
    warn_of_what_is_replaced)."""
    if inputs.kind == mti.MAP:
        input_map = inputs.map
        run_input(input_map)
        outfile, in_place = run_map_path(input_map, outdir, inplace, force)
        refuse_an_input_as_the_map(outfile, in_place, [assembly, depth, *others])
        check_normalization(normalization, enzyme, depth, input_map)
        resolved_parameters(map=input_map)
        return input_map, outfile, hold_map(held, command, input_map, outfile, in_place, rounds=rounds, inplace=inplace)
    outfile = made_map_path(inputs, name, outdir, inplace, force)
    refuse_an_input_as_the_map(outfile, False, [assembly, depth, *inputs.files, *others])
    check_normalization(normalization, enzyme, depth, None)
    hold_lock(held, outfile, in_place=False)
    resolved_parameters(map=outfile)
    return outfile, outfile, None


def make_run_map(command, inputs, held, outfile, contig_data, assembly, min_quality, threads, tmp_dir):
    """Make the map of the contigs of a run (command) from its pairs or BAM files
    at outfile (see metator.inputs.make_map), the pairs of its BAM files in its
    temporary directory, tmp_dir, then open it, as the run reads it from now
    on, under its lock (see hold_map): the map that the run updates in place.
    A failure is logged as an error, then raised (a MetacoolError): an earlier
    map at outfile is then left as it was."""
    try:
        mti.make_map(
            inputs,
            contig_data,
            outfile,
            assembly=assembly,
            min_quality=min_quality,
            threads=threads,
            tmp_dir=tmp_dir,
            command=shlex.join(["metator", *sys.argv[1:]]),
        )
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    return hold_map(held, command, outfile, outfile, True)


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


def refuse_an_input_as_the_map(outfile, in_place, inputs):
    """Refuse, before any work, a map that the run writes in its output directory
    (outfile, unless it updates its input map in place) that is one of the
    files of the inputs of the run (the assembly, a table): the run would
    replace it with its map."""
    if in_place or not os.path.exists(outfile):
        return
    for path in inputs:
        if path is not None and os.path.isfile(path) and os.path.samefile(path, outfile):
            reason = (
                f"{outfile}, where the run writes its map, is {path}, an input of the run: give the run another output "
                "directory with -o."
            )
            logger.error(reason)
            raise mtmc.MetacoolError(reason)


def write_elsewhere(inplace):
    """How a run that updates its input map in place would leave it as it is, for
    the messages that refuse or warn of the update: without --inplace
    (inplace), the run writes its map in its output directory; a run whose
    output directory holds its input map already (see run_map_path) writes it
    in another one."""
    return "leave out --inplace" if inplace else "give the run another output directory with -o"


# What a run that updates its map in place replaces of what the map holds
# already (see warn_of_what_is_replaced): its clusters (whether it does), its
# clustering rounds made by these steps, and its rounds of these names, by
# default (mge replaces its own round, mge or mge_random; pipeline --mges, the
# round mge too).
REPLACED_IN_PLACE = {
    "network": (False, (), ()),
    "pipeline": (True, ("partition", "validation"), ()),
    "partition": (False, ("partition", "validation"), ()),
    "validation": (True, ("validation",), ()),
    "mge": (False, (), (mtcl.MGE,)),
    "host": (False, (), ()),
}

# The runs that replace the hosts of the MGE clusters that a map holds (see
# metator.metacool.write_map): host itself, and the runs that replace what
# they describe, the round mge (mge, but with --random), the rounds
# (partition) or the genomes (validation, pipeline). network replaces them
# only when it stores another network, which it finds once it is built.
REPLACES_HOSTS = ("host", "partition", "validation", "pipeline")


def warn_of_what_is_replaced(command, metacool, rounds=None, inplace=False):
    """Warn, before any work, that a run (command) that updates the map metacool
    in place (with --inplace: inplace, or as the map of its output directory)
    replaces what the map holds already: its clusters (the genomes other than
    __unbinned__, of an earlier run, or of another binner), and its clustering
    rounds made by an earlier run (see REPLACED_IN_PLACE): those that the steps
    of the binning made, and those that the run stores, by name (rounds: by
    default, those of REPLACED_IN_PLACE; mge gives its own round, mge or
    mge_random, and pipeline --mges the round mge).

    Neither can be found again in the updated map, and neither can the hosts
    of its MGE clusters that metator host stored in it, which the run replaces
    too (see REPLACES_HOSTS)."""
    clusters, steps, names = REPLACED_IN_PLACE[command]
    replaces_hosts = command in REPLACES_HOSTS or (command == "mge" and mtcl.MGE in (rounds or ()))
    names = names if rounds is None else tuple(rounds)
    try:
        with mtmc.read_errors(metacool):
            mc = mtmc.open_map(metacool)
            n_clusters = int((mc.genomes["name"].astype(str) != "__unbinned__").sum()) if clusters else 0
            rounds = [
                str(entry["name"]) for entry in mc.cluster_rounds if entry.get("step") in steps or entry.get("name") in names
            ]
            hosts = replaces_hosts and mtmc.holds_hosts(mc)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    held = []
    if n_clusters:
        held.append(f"{n_clusters} cluster(s)")
    if rounds:
        held.append(f"the clustering round(s) {', '.join(rounds)}")
    if hosts:
        held.append("the hosts of its MGE clusters (metator host)")
    if held:
        logger.warning(
            f"{metacool} holds {' and '.join(held)}, which this run replaces, in place: {write_elsewhere(inplace)} to "
            f"write the map of this run to another file, and leave {metacool} as it is."
        )


def step_map(command, input_map, outdir, inplace, force, held, inputs=(), rounds=None):
    """The map that a partition, validation, mge or host run (command) writes
    (see run_map_path), and the map as it reads it (see hold_map)."""
    run = command in ("mge", "host")
    map_of = f"the map of the {'run' if run else 'contigs of the run'} as the first argument of {command}"
    path = refuse_missing_map(input_map)
    if not mtmc.is_map(input_map):
        reason = f"{input_map} is not a .metacool contact map: give {map_of}."
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    refuse_a_map_named_as_a_uri(input_map, path)
    open_input_map(input_map)
    outfile, in_place = run_map_path(input_map, outdir, inplace, force)
    refuse_an_input_as_the_map(outfile, in_place, inputs)
    return outfile, hold_map(held, command, input_map, outfile, in_place, rounds, inplace)


def start_step(
    command,
    outdir,
    cluster_dirs,
    input_map,
    inplace,
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
    run left them, as, with partition (a validation run), the plots of an
    earlier run (see refuse_earlier_plots), which force replaces, and what a
    run of an earlier version wrote in outdir, with or without force (see
    refuse_earlier_outputs). Returns the time the run started too (see
    generate_log_header)."""
    log_file = log_file_name(outdir, mtr.sample_name(input_map))
    started = generate_log_header(log_file, outdir=outdir)
    metacool_file, planned = step_map(command, input_map, outdir, inplace, force, held, inputs=inputs)
    stored = stored_network(input_map)
    clusters = map_partition(input_map, planned, stored.digest) if partition else None
    if assembly is not None:
        check_assembly(assembly)
    bins_prefix = bin_prefix(prefix)
    # Refused, all of them, before anything is replaced: and the directory
    # that holds them (fasta) when it is no directory.
    for path in cluster_dirs:
        refuse_no_directory(os.path.dirname(path), "the FASTA files of its clusters")
        refuse_output_dir(path, force)
    plots = []
    if partition:
        refuse_earlier_outputs(outdir, keep=(input_map, metacool_file))
        plots = refuse_earlier_plots(outdir, force)
    for path in cluster_dirs:
        make_output_dir(path, force)
    remove_earlier_plots(plots)
    return started, metacool_file, bins_prefix, planned, stored, clusters


def stored_network(metacool, mc=None):
    """The network that metator network stored in the map of a partition,
    validation or host run, metacool, read in full before the run replaces
    anything (see metator.network.stored_network), through mc when it is given
    (the map as the run read it): the contig data of its contigs, the network
    and its normalization, by the positions of the contigs in the map, which
    the partition that a validation run reads from the map numbers them by too
    (see metator.metacool.read_partition)."""
    try:
        return mtn.stored_network(metacool, mc)
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


def generate_log_header(log_path, **values):
    """Start the log file of a run, as hicstuff starts its own: the version of
    metaTOR, the date, the command line (cmd) and the value of each parameter
    of the command (see command_parameters; values: the values that the run
    resolved, by parameter name, such as its output directory), then a line of
    dashes.

    Returns the time the run started, which the last line of its log counts
    from (see log_end)."""
    started = time.monotonic()
    mtl.set_file_handler(log_path, formatter=logging.Formatter(""))
    logger.info(f"## MetaTOR: v{__version__} log file")
    logger.info(f"## date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"## cmd: {shlex.join(['metator', *sys.argv[1:]])}")
    for name, value in command_parameters(values):
        logger.info(f"## {name}: {value}")
    logger.info("---")
    mtl.set_file_handler(log_path, formatter=mtl.logfile_formatter)
    return started


def command_parameters(values=None):
    """The parameters of the running command, as the header of its log lists
    them: in the order of its help, the name of each (its long option without
    the dashes, or MAP as map) and its value, as given or by default, or as the
    run resolved it (values, by parameter name). None outside a command."""
    ctx = click.get_current_context(silent=True)
    if ctx is None:
        return []
    given = {**ctx.params, **(values or {})}
    return [
        (next((opt[2:] for opt in param.opts if opt.startswith("--")), param.name), given[param.name])
        for param in ctx.command.params
        if param.name in given
    ]


def log_end(started, what, where=None):
    """End the log of a run, as hicstuff ends its own: what the run wrote (where:
    the map, or the directory of the tables), and how long it took since it
    started (see generate_log_header)."""
    logger.info(f"{what} after {mtl.duration(time.monotonic() - started)}" + (f": {where}" if where else "."))


def log_stage(number, name, what, stages=3):
    """The first line of the log of a step of a pipeline run, one of its stages:
    network, partition, validation and export; one more first from pairs or
    BAM files (the map), and two more after the validation with --mges (mge,
    host): which step, its name, and what it does, so that the log of the run
    reads step by step."""
    logger.info(f"[{number}/{stages}] {name}: {what}")


def log_file_name(directory, name):
    """Log file of a run named name (its map without .metacool, see
    metator.report.sample_name and run_name): <directory>/<name>.metator_<date>.log,
    as hicstuff names its own (<prefix>.hicstuff_<date>.log)."""
    now = time.strftime("%Y%m%d%H%M%S")
    return join(directory, f"{name}.metator_{now}.log")


def run_name(arguments, prefix, assembly):
    """The name of a network or pipeline run, which names its log and its report:
    that of its map without .metacool (see metator.report.sample_name), the map
    given as its arguments, or the map that it makes from pairs or BAM files,
    <prefix>.metacool or <assembly>.metacool (see metator.inputs.map_name).
    Arguments that the run refuses are refused once its log is started."""
    try:
        inputs = mti.resolve(arguments)
    except Exception:
        inputs = None
    if inputs is not None and inputs.kind == mti.MAP:
        return mtr.sample_name(inputs.map)
    return mtr.sample_name(mti.map_name(prefix, assembly))


#: The key of the context of a command (its meta) that holds the values that
#: the run resolved, which its metadata store with its parameters (see
#: resolved_parameters).
RESOLVED = "metator.resolved_parameters"


def resolved_parameters(**values):
    """Store, with the parameters of the running command, values that the run
    resolved (by name): the map that a network or pipeline run starts from
    (map), which its inputs give, or which it makes from them. run_metadata
    stores them, after its parameters."""
    cast(click.RichContext, click.get_current_context()).meta.setdefault(RESOLVED, {}).update(values)


def run_metadata():
    """The command line and the parameters of the running command, as stored in a
    .metacool file, with the values that the run resolved (see
    resolved_parameters)."""
    # Outside a command, get_current_context raises: it never returns None here.
    ctx = cast(click.RichContext, click.get_current_context())
    params = {**ctx.params, **ctx.meta.get(RESOLVED, {})}
    return {
        "command": shlex.join(["metator", *sys.argv[1:]]),
        "parameters": {name: list(value) if isinstance(value, tuple) else value for name, value in params.items()},
    }


def hold_lock(held, path, in_place, inplace=False):
    """Hold the lock of the map at path until the end of the run (in held, an
    ExitStack): metacooler's (metacooler.lock), which every writer of the map
    takes too, and which the edit session that writes the map at the end of the
    run re-enters (see metator.metacool._staged_edit).

    In place, the lock is the file a symbolic link points to, as metacooler
    edits it; otherwise, the map of the output directory, its own name, so
    runs from one map, each in its own output directory, run side by side.
    inplace: whether the run updates its map in place by --inplace, for the
    message that refuses it."""
    import metacooler

    try:
        # metacooler's warnings (a stale lock taken over, what a killed run
        # left beside the map) go to the log of the run.
        with mtmc._metacooler_logging():
            command = shlex.join(["metator", *sys.argv[1:]])
            held.enter_context(metacooler.lock(path, command=command, follow_symlinks=in_place))
    except OSError as error:
        # MapLockError, and MapLockedError: OSErrors of metacooler.
        elsewhere = f" Or write this run's map elsewhere: {write_elsewhere(inplace)}." if in_place else ""
        message = f"{str(error).rstrip('.')}.{elsewhere}"
        logger.error(message)
        raise mtmc.MetacoolError(message) from error


def hold_map(held, command, input_map, outfile, in_place, rounds=None, inplace=False):
    """Hold the lock of the map that a run (command) writes, outfile, until the
    end of the run (see hold_lock), then open its input map: the map as the run
    reads it, under that lock when the run updates it in place (in_place; by
    --inplace: inplace), which the map of the run is refused to be written
    from if it changed since (see metator.metacool.write_map's planned)."""
    hold_lock(held, outfile, in_place, inplace)
    planned = open_input_map(input_map)
    if in_place:
        refuse_a_read_only_map(outfile, inplace)
        warn_of_what_is_replaced(command, input_map, rounds, inplace)
    return planned


def refuse_a_read_only_map(path, inplace=False):
    """Refuse, before any work, the map that a run updates in place (path, its
    input map; by --inplace: inplace) when it has no write permission (judged
    by its mode too, as a run as root may write any file): metacooler refuses
    to replace a read-only map, which it would do at the end of the run."""
    mode = os.stat(path).st_mode
    if mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH) and os.access(path, os.W_OK):
        return
    message = (
        f"{path} is read-only (mode {stat.filemode(mode)}): the run would replace it with the map of its clusters. "
        f"Make it writable (chmod u+w), or write the map of the run to another file: {write_elsewhere(inplace)}."
    )
    logger.error(message)
    raise mtmc.MetacoolError(message)


def default_map(input_map, outdir="."):
    """The map that a run writes by default (see run_map_path): in its output
    directory, outdir, under the name of the file of its input map, or, for
    an entry of a .mmetacool collection (file.mmetacool::entries/NAME),
    NAME.metacool: a map of its own."""
    if mtmc.is_entry(input_map):
        name = str(input_map).split("::", 1)[1].strip().strip("/").rsplit("/", 1)[-1]
        name = name if name.endswith(".metacool") else f"{name}.metacool"
    else:
        name = os.path.basename(mtmc.map_file(input_map))
    # The current directory, as ., ./ or none: the name alone.
    return name if os.path.normpath(outdir or ".") == "." else join(outdir, name)


def run_map_path(input_map, outdir=".", inplace=False, force=False):
    """The map a run writes, and whether it updates its input map in place.

    By default, the map of the output directory of the run, outdir, of the name
    of its input map (see default_map): a new map, which the run writes from a
    copy of its input map, left as it was (metacooler stages the copy beside
    it, and renames it there once it is written and validated). A symbolic
    link there is a file of its own, which the map replaces, not its target.
    When it is the input map's file itself (its symbolic links followed, see
    metator.metacool.is_input_map), the run updates it in place: a run step by
    step, on the map that the step before wrote in the same output directory.

    With inplace (--inplace), the input map itself: the file that a symbolic
    link points to, as metacooler edits that file (a map given as file::/,
    the root of the file, is the file).

    The log says which map the run writes, and whether it updates its input map.

    Raises
    ------
    MetacoolError
        With inplace, if the input map is an entry of a .mmetacool collection;
        otherwise, if the map of the output directory holds ::, is a
        directory, the collection of the input map (an entry) or a hard link of
        the input map, or exists and is not the input map, without force, or if
        the output directory cannot be made. It is raised before the run does
        any work, and its reason is logged as an error.
    """
    if inplace:
        if mtmc.is_entry(input_map):
            message = (
                f"{input_map} is an entry of a .mmetacool collection, which a run does not update in place: leave out "
                "--inplace, and the run writes its map, a map of its own, in its output directory."
            )
            logger.error(message)
            raise mtmc.MetacoolError(message)
        logger.info(f"The run updates {input_map} in place (--inplace).")
        # The file itself (file::/, the root of the file, names it too), or the
        # file that a symbolic link points to, which metacooler edits.
        path = mtmc.map_file(input_map)
        return (os.path.realpath(path) if os.path.islink(path) else path), True
    outfile = default_map(input_map, outdir)
    reason = None
    if "::" in outfile:
        reason = (
            f"{outfile}, the map of the run in its output directory, holds ::, which metacooler reads as file::group: "
            "give the run another output directory with -o."
        )
    elif os.path.isdir(outfile):
        reason = (
            f"{outfile} is a directory, where the run writes its map: move it away, or give the run another output "
            "directory with -o."
        )
    elif mtmc.is_entry(input_map) and mtmc.is_input_map(input_map, outfile):
        reason = (
            f"{outfile} is the collection of the input map {input_map}: a run does not replace a collection with its "
            "map. Give the run another output directory with -o."
        )
    elif mtmc.is_hard_link_of(input_map, outfile):
        reason = (
            f"{outfile} is a hard link of the input map {input_map}: the same file under another name, where the run "
            "writes its map. Give the run another output directory with -o, or --inplace to update the map in place."
        )
    elif os.path.lexists(outfile) and not force and not mtmc.is_input_map(input_map, outfile):
        link = mtmc.is_input_map(input_map, os.path.realpath(outfile)) and os.path.islink(outfile)
        what = " (a symbolic link to the input map, which -F replaces, not the map)" if link else ""
        reason = f"{outfile} exists{what}: give -F to replace it with the map of this run, or another output directory with -o."
    if reason:
        logger.error(reason)
        raise mtmc.MetacoolError(reason)
    # Where metacooler locks and stages the map, before it renames it there.
    make_outdir(os.path.dirname(outfile) or ".")
    in_place = mtmc.is_input_map(input_map, outfile)
    if in_place:
        logger.info(
            f"The run updates {input_map} in place: it is the map of its output directory (give another one with -o to "
            "leave it as it is)."
        )
    else:
        logger.info(f"The run writes its map to {outfile}, and leaves {input_map} as it is (--inplace updates it).")
    return outfile, in_place


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
            logger.debug(f"{outfile} is the input map {metacool}: it is replaced once {what} is validated.")
        elif os.path.lexists(outfile):
            logger.info(f"Removing {outfile}, the contact map of an earlier run.")
            os.remove(outfile)
        mtmc.write_map(metacool, outfile, force=True, planned=planned, **stores)
    except Exception as error:
        reason, unexpected = _failure(error)
        message = f"{outfile} is not written: {reason}. {kept}The other outputs of the run are kept."
        logger.error(message, exc_info=unexpected)
        raise mtmc.MetacoolError(message) from error


def run_network(metacool, planned, assembly, depth, enzyme, normalization, tmp_dir, contig_data=None):
    """The network of the contigs of a network or pipeline run, from its map
    (metacool, as the run read it under its lock: planned): the contig data,
    measured on the assembly (contig_data, when the run measured it already,
    to make its map) and read from the map (see contig_data_from_map), and the
    hits of the contigs of the map, counted on its contacts (see
    map_contacts)."""
    if contig_data is None:
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
    return mtmc.partition_network(planned) != _network_digest(contacts, table, normalization)


def stale_hosts(planned, contacts, table, normalization):
    """Whether the input map of a network run, planned, holds hosts of its MGE
    clusters (see metator.metacool.holds_hosts) that metator host found on
    another network than the one the run stores (see
    metator.metacool.hosts_network and stale_partition): the run drops them
    (see metator.metacool.write_map)."""
    if not mtmc.holds_hosts(planned):
        return False
    return mtmc.hosts_network(planned) != _network_digest(contacts, table, normalization)


def _network_digest(contacts, table, normalization):
    """The digest of the network that a network run stores (see
    metator.network.network_digest): from the contacts of its map and its
    contig data (table), with its normalization."""
    return mtn.network_digest(contacts, mtn.network_weights(table, normalization), table["Name"], normalization)


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
    outdir,
    tmp_dir,
    normalization,
):
    """The validation of the partition of a validation or pipeline run, from its
    map (metacool): the contig data of the contigs of the map (table, in its
    order), their network (network, by their positions in the map) and its
    partition (partition, a metator.clusters.Partition: the overlapping
    clusters that reach size are checked). The FASTA files of the final
    clusters go to final_fasta_dir (fasta/mags_with-mges of outdir), and their
    plots to plots/ of outdir, the output directory of the run (see
    metator.figures.plot_figures)."""
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
        plot_outdir=outdir,
    )
    final = contigs_data["Final_cluster"].astype(str).str.strip()
    binned = ~final.isin(mtmc.UNBINNED)
    assignment = dict(zip(contigs_data["Name"].astype(str)[binned], final[binned], strict=True))
    check_final_clusters(assignment, final_fasta_dir)
    return mtmc.FinalClusters(
        assignment=assignment,
        summary=join(tmp_dir, mio.SUMMARY_FULL_FILE),
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
    "assembly, its decompressed copy: give it room for both. Not in a directory that the run makes anew "
    "(fasta/mags_with-mges...), nor in a directory of the results of an earlier run (recursive_clusters, "
    "final_clusters...), which the run refuses."
)
# -F of validation and pipeline (see refuse_earlier_plots and EARLIER_OUTPUTS).
FORCE_RUN_HELP = (
    "Replace the FASTA files of the final clusters (fasta/mags_with-mges), their plots (plots/bins_distribution.pdf "
    "and plots/bins_size_distribution.pdf: the other files of plots are left as they are) and the map of an earlier "
    "run in the output directory. The other results of an earlier run in the output directory (the final_clusters "
    "and plot of an earlier version, its tables, recursive_clusters, 1.x's outputs) are refused, with or without -F: "
    "they are never removed."
)
NO_CLEAN_UP_TABLES_HELP = (
    "Do not remove the temporary directory of the run, which the log names: it keeps miComplete's results, the FASTA "
    "files of the clusters it checked and the tables of the validation (contig_data_final.txt, cluster_summary.txt "
    "and cluster_summary_full.txt, its values at full precision, which the map stores, binning.txt, clusters.txt, "
    "cluster_members.txt)."
)
SIZE_HELP = (
    "Threshold size to keep clusters in base pair: miComplete checks the overlapping clusters of the partition that "
    "reach it, and the clusters of the recursive partition larger than it."
)
# -F of pipeline: that of validation, the outputs of the run (see
# metator.outputs), and what --mges writes, as mge and host write it (see
# MGE_EARLIER_OUTPUTS and HOST_EARLIER_OUTPUTS).
FORCE_PIPELINE_HELP = (
    FORCE_RUN_HELP + " It replaces fasta/mags, cool and the tables of an earlier run in tables too, and, with --mges, "
    f"fasta/mges; {mto.TABLES_DIR}/{mth.HOSTS_FILE} and {mtf.PLOT_DIR}/{mth.PLOT_FILE} are replaced without it, as host "
    "replaces them. The MGE files of an earlier version (mges_binned.fa, mge_hosts.txt...) are refused, with or "
    "without it."
)
# The thresholds of host, which pipeline --mges takes as --mge-min-share and
# --mge-min-contigs.
MIN_SHARE_HELP = (
    "The share of the signal of an MGE cluster with the MAGs, in percent, from which a MAG is a host. The contacts of "
    "the MGE cluster with the contigs in no MAG and with the other MGE clusters do not count."
)
MIN_CONTIGS_HELP = "The number of contigs of a MAG in contact with an MGE cluster from which the MAG is a host."


def assembly_option(help):
    return _option("-a", "--assembly", required=True, metavar="FILE", panel="Input", help=help)


MAP_HELP = (
    "A .metacool contact map of the contigs of the assembly, made from the pairs of the run (by hicstuff pipeline -M "
    "metacool or metacooler create pairs), or an entry of a .mmetacool collection (file.mmetacool::entries/NAME): "
    "the run reads its network from it."
)


def written_map_help(what):
    """What the help of MAP says of the map that a run writes (see run_map_path):
    what it stores in it."""
    return (
        f"The run writes it, with {what}, to the output directory, under its name (an entry as NAME.metacool), and "
        "leaves MAP as it was; but MAP is updated in place when it is that map already (the steps of a run, one by "
        "one, in one output directory), or with --inplace."
    )


INPLACE_HELP = (
    "Update MAP in place, rather than write the map of the run to the output directory: refused for an entry of a "
    ".mmetacool collection."
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


def input_help(made, what):
    """What the help of INPUT says (network and pipeline): the map that a run
    starts from, or the pairs or BAM files that it makes it from, named as made
    says; and what the run stores in its map (what)."""
    return (
        "The .metacool contact map of the contigs of the assembly (made by hicstuff pipeline -M metacool or metacooler "
        "create pairs), or an entry of a .mmetacool collection (file.mmetacool::entries/NAME): the run reads its network "
        f"from it, and writes it, with {what}, to the output directory, under its name (an entry as NAME.metacool), "
        "leaving it as it was, but with --inplace or when it is that map already. Or the pairs files, or the BAM files, "
        "of the Hi-C libraries of the run, summed: the run makes the map of the contigs from them first, one bin per "
        f"contig, in the output directory, {made}, then updates it in place. Each argument is a path, a "
        "comma-separated list or a glob pattern (quoted). A BAM file holds both mates of its reads (paired-end "
        "alignments, as bwa mem -5SP writes them), or one mate each (as bowtie2 aligns them): give the R1 files, then "
        "the R2 files."
    )


def input_argument(help):
    """The inputs of a network or pipeline run, its arguments: the .metacool map
    of the contigs, or the pairs files or the BAM files that the run makes it
    from (see metator.inputs.resolve and run_inputs)."""
    return click.argument(
        # Shown as INPUT...: one or several.
        "input",
        nargs=-1,
        required=True,
        # A directory is refused as a usage error. An entry of a collection
        # (file::group) is no path, nor a comma-separated list or a glob
        # pattern: run_inputs checks them.
        type=click.Path(dir_okay=False),
        panel="Input",
        help=help,
    )


def min_quality_option():
    """-q, --min-quality: the mapping quality that both mates of a read of a BAM
    file must reach for its pair to count (see metator.inputs.bam_pairs)."""
    return _option(
        "-q",
        "--min-quality",
        type=click.INT,
        default=30,
        show_default=True,
        metavar="INT",
        panel="Input",
        help="BAM files: the mapping quality (MAPQ) that both mates of a read must reach for its pair to count. Refused "
        "with a map or pairs files, whose contacts are read as they are.",
    )


#: --inplace of network and pipeline, whose map may be made from pairs or BAM
#: files rather than given.
RUN_INPLACE_HELP = (
    "Update the map given as INPUT in place, rather than write the map of the run to the output directory: refused "
    "for an entry of a .mmetacool collection, and for pairs or BAM files, whose map the run makes in its output "
    "directory."
)


def inplace_option(help):
    """--inplace: the run updates its map in place, rather than write the map of
    its output directory (see run_map_path)."""
    return _option("--inplace", is_flag=True, help=help)


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


def prefix_option(what, more=""):
    return _option(
        "-p",
        "--prefix",
        metavar="STR",
        help=f"Name {what} STR_metator_... (by default, metator_...).{more} Only letters, digits, '_' and '-'.",
    )


def force_option(help):
    return _option("-F", "--force", is_flag=True, help=help)


def tmpdir_option(help="Temporary directory."):
    return _option("-T", "--tmpdir", default="./tmp", show_default=True, metavar="DIR", help=help)


def threads_option(help):
    return _option("-t", "--threads", type=click.INT, default=1, show_default=True, metavar="INT", help=help)


def seed_option(help, panel="Partition"):
    """--seed: the seed of the random numbers of a run that clusters at random
    (see seed_random)."""
    return _option("--seed", type=click.IntRange(0, 2**32 - 1), metavar="INT", panel=panel, help=help)


def seed_random(seed):
    """Seed the random numbers that a run draws, with --seed (seed): those of
    Python's random module, which igraph draws from in Louvain and Leiden, and
    numpy's, which the random clustering of mge --random draws from. None
    leaves them as they are: the run draws other numbers each time.

    A step is seeded once, before its iterations, which draw from the seeded
    numbers one after the other: each iteration of Louvain or Leiden gets other
    numbers, and so another membership, and the step as a whole, all its
    iterations in their order, is the same from the same seed."""
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)


def no_clean_up_option(help="Do not remove temporary files."):
    return _option("-N", "--no-clean-up", is_flag=True, help=help)


# ---------------------------------------------------------------------------
# CLI group
# ---------------------------------------------------------------------------


@click.group("metator", context_settings={"help_option_names": ["-h", "--help"]})
@click.rich_config(help_config=click.RichHelpConfiguration(text_markup="rich", show_arguments=True, max_width=100))
@click.command_panel("Pipeline", commands=["pipeline"])
@click.command_panel("Steps", commands=["network", "partition", "validation", "mge", "host"])
@click.command_panel("Outputs", commands=["dump", "export", "report"])
@click.version_option(__version__, "-v", "--version", prog_name="metator", message="%(version)s")
def cli():
    """Bin metagenomes from their Hi-C contacts.

    [cyan]pipeline[/] runs the steps in one go: [cyan]network[/],
    [cyan]partition[/], [cyan]validation[/], then [cyan]mge[/] and
    [cyan]host[/] with --mges. [cyan]dump[/], [cyan]export[/] and
    [cyan]report[/] write the outputs of a run from its map.
    """


# ---------------------------------------------------------------------------
# network
# ---------------------------------------------------------------------------


@cli.command(
    short_help="Build the metaHiC network of the contigs of a .metacool map, or of pairs or BAM files, and store it in the map."
)
@input_argument(help=input_help("named after the FASTA file of the assembly (<assembly>.metacool)", "the network of its contigs"))
@network_options("-n")
@min_quality_option()
@outdir_option(help="The output directory, of the map of the run and of its log. Default: current directory.")
@force_option(
    "Replace the map of an earlier run in the output directory: that of the map given, or the map that the run makes "
    "from pairs or BAM files."
)
@tmpdir_option(
    help="Temporary directory, where the run writes, in a directory of its own, the contig data that it stores in the "
    "map (contig_data_network.txt), and the pairs of the BAM files that it makes its map from (about 25 bytes per "
    "pair, plus the names of its two contigs)."
)
@threads_option(help="Threads reading the BAM files and binning the pairs, when the run makes its map from them.")
@no_clean_up_option(help="Keep the temporary directory of the run, which the log names, and its contig_data_network.txt.")
@inplace_option(help=RUN_INPLACE_HELP)
def network(
    input,  # the arguments INPUT...
    assembly,
    depth,
    enzyme,
    normalization,
    min_quality,
    outdir,
    force,
    tmpdir,
    threads,
    no_clean_up,
    inplace,
):
    """Build the metaHiC network of the contigs, and store it in the map.

    [cyan]metator network INPUT... -a ASSEMBLY[/] builds the network from a
    .metacool contact map of the contigs (metacooler merge sums the maps of
    several libraries): that of the pairs it was made of. hicstuff
    ([cyan]hicstuff pipeline -M metacool[/]) or metacooler
    ([cyan]metacooler create pairs[/]) makes the map from the reads or the
    pairs; or the run makes it from the pairs files or the BAM files of its
    libraries, one bin per contig, as metacooler create pairs does, in the
    output directory, named after the assembly.

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
    it until [cyan]metator partition MAP[/] partitions the new network; but
    the hosts of the MGE clusters that [cyan]metator host[/] found on another
    network are dropped, with a warning. The run writes the map to the output
    directory, under the name of the map given, which it leaves as it was
    (--inplace updates it in place), or updates in place the map that it made
    from pairs or BAM files, and locks it while it writes it. [cyan]metator
    dump[/] writes the network, as
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
        make_outdir(outdir)

        # Enable file logging
        log_file = log_file_name(outdir, run_name(input, None, assembly))
        started = generate_log_header(log_file, outdir=outdir, input=shlex.join(input))

        # The inputs, the map that the run starts from (or makes from its
        # pairs or BAM files), the map that it writes (that of the output
        # directory, by default), the normalization against its columns, and
        # the assembly, checked before any work. The run holds the lock of the
        # map it writes from now on: a second run on the same map is refused.
        # It reads its map under that lock (planned): the map its network is
        # written from, which is refused if it changed since.
        inputs = run_inputs(input, given_on_the_command_line("min_quality"))
        map, metacool_file, planned = run_map(
            "network",
            inputs,
            held,
            mti.map_name(None, assembly),
            assembly,
            depth,
            enzyme,
            normalization,
            outdir,
            inplace,
            force,
        )
        check_assembly(assembly)

        # From pairs or BAM files, the map of the contigs, made in the output
        # directory, which the run then updates in place, from the contig data
        # measured on the assembly.
        contig_data = None
        if planned is None:
            contig_data = mtn.create_contig_data(assembly, depth, enzyme)
            planned = make_run_map("network", inputs, held, metacool_file, contig_data, assembly, min_quality, threads, tmp_dir)

        # The contig data, measured on the assembly and read from the map,
        # the hits of the contigs of the map, counted on its contacts, then
        # the network, stored in the map.
        network, contacts = run_network(map, planned, assembly, depth, enzyme, normalization, tmp_dir, contig_data)
        # A partition of another network, which the map keeps: said once the
        # map is written. The hosts of the MGE clusters found on another
        # network, which it drops, with a warning.
        stale = stale_partition(planned, contacts, network.table, normalization)
        drop_hosts = stale_hosts(planned, contacts, network.table, normalization)
        write_map(map, metacool_file, planned, "the map with the network of its contigs", network=network, drop_hosts=drop_hosts)
        if stale:
            logger.warning(
                f"{metacool_file} holds a partition of another network than the one stored in it now, which the run "
                f"kept: metator validation refuses it until metator partition {metacool_file} runs again, which "
                "partitions this network."
            )

    log_end(started, "Contact map with the network written", metacool_file)


# ---------------------------------------------------------------------------
# partition
# ---------------------------------------------------------------------------


@cli.command(short_help="Partition the metaHiC network of a .metacool map with the Louvain or Leiden algorithm.")
@map_argument(
    required=True,
    help="The .metacool contact map of the contigs of the run, which metator network stored the network in, or an "
    "entry of a .mmetacool collection (file.mmetacool::entries/NAME). "
    + written_map_help("the core and overlapping clusters of the run in its clustering rounds"),
)
@partition_options(
    iterations=100,
    overlap=80,
    iterations_help="Number of iterations of Louvain or Leiden.",
    overlap_help="Hamming distance threshold to use to merge clusters (percentage).",
)
@outdir_option(default="./", help="The output directory, of the map of the run and of its log.")
@prefix_option("the overlapping clusters")
@force_option("Replace the map of an earlier run in the output directory.")
@tmpdir_option(help="Temporary directory: partition writes no temporary file in it.")
@threads_option(help="Number of parallel threads allocated for the partition.")
@no_clean_up_option(help="Keep the temporary directory of the run: partition leaves nothing in it.")
@inplace_option(help=INPLACE_HELP)
@seed_option(
    "Seed of the random numbers of Louvain or Leiden, set once, before the iterations, which draw other numbers one "
    "after the other: from the same map, the same seed gives the same partition. By default, each run draws other "
    "random numbers."
)
def partition(
    map,  # the argument MAP
    algorithm,
    iterations,
    overlap,
    res_param,
    seed,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
    inplace,
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
    /clusters group), in place of those it held, in the map that it writes
    to the output directory (MAP itself when it is there already, as after
    [cyan]metator network MAP -o DIR[/] with the same DIR, or with
    --inplace). Its genomes and contigs tables are left as they were:
    [cyan]metator validation MAP[/] checks the overlapping clusters that
    reach its --size, and writes the final clusters there.
    """
    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up), ExitStack() as held:
        # Defined the output directory.
        outdir = outdir or "."
        make_outdir(outdir)

        # Enable file logging, the map of the run checked and locked, and
        # the network that metator network stored in it read, with the data
        # of its contigs (see start_step).
        started, metacool_file, bins_prefix, planned, stored, _ = start_step(
            "partition", outdir, [], map, inplace, force, held, prefix
        )

        # Partition the network
        seed_random(seed)
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

    log_end(started, "Contact map with the partition written", metacool_file)


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


@cli.command(short_help="Validate the clusters with miComplete and decontaminate them with a recursive partition.")
@map_argument(
    required=True,
    help="The .metacool contact map of the contigs of the run (that of network MAP and partition MAP), which metator "
    "network stored the network in, and metator partition its partition, or an entry of a .mmetacool collection "
    "(file.mmetacool::entries/NAME). "
    + written_map_help("the clusters of the run, as pipeline MAP writes it (the map of the clusters)"),
)
@assembly_option(CHECKED_ASSEMBLY_HELP)
@partition_options(
    iterations=10,
    overlap=90,
    iterations_help="Number of recursive iterations of Louvain or Leiden.",
    overlap_help="Hamming distance threshold to use to merge clusters (percentage).",
    size_help=SIZE_HELP,
)
@outdir_option(
    default="./",
    help="The output directory, of the map of the clusters, fasta/mags_with-mges (the FASTA files of the final clusters), "
    "plots and the log.",
)
@prefix_option("the recursive clusters")
@force_option(FORCE_RUN_HELP)
@tmpdir_option(help=TMPDIR_TABLES_HELP)
@threads_option(help="Number of parallel threads, for the partition and miComplete.")
@no_clean_up_option(help=NO_CLEAN_UP_TABLES_HELP)
@inplace_option(help=INPLACE_HELP)
@seed_option(
    "Seed of the random numbers of the recursive partition, set once, before it, its iterations drawing other "
    "numbers one after the other: from the same map, the same seed gives the same clusters, and partition --seed S "
    "then validation --seed S give those of pipeline --seed S. By default, each run draws other random numbers."
)
def validation(
    map,  # the argument MAP
    assembly,
    algorithm,
    iterations,
    overlap,
    res_param,
    size,
    seed,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
    inplace,
):
    """Use miComplete to validate the clusters.

    The clusters to validate are the overlapping clusters of the partition
    that [cyan]metator partition MAP[/] stored in the map that reach
    --size: the run writes their FASTA files from the assembly, in its
    temporary directory, for miComplete to check.

    Recursively decontaminate the contaminated clusters evaluated by
    miComplete. The command writes the FASTA files of the final clusters in
    the output directory (fasta/mags_with-mges: the contigs of each final
    cluster, the MGE contigs that it holds included), and their plots (plots:
    the completeness and contamination of each, bins_distribution.pdf, and
    the size of each quality class, bins_size_distribution.pdf), and stores
    the clusters in the map: the tables of the validation stay in its temporary
    directory, which -N keeps, and [cyan]metator dump MAP[/] writes them
    from the map.

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
    [cyan]metator pipeline MAP[/] does, to the output directory (in place of
    MAP when it is there already, or with --inplace).
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
        make_outdir(outdir)
        final_fasta_dir = join(outdir, mto.MAGS_WITH_MGES_DIR)

        # Enable file logging, the map of the run checked and locked, the
        # network that metator network stored in it and the partition that
        # metator partition stored in it read, the assembly checked, the
        # directory of the final clusters made, and the plots of an earlier
        # run removed (see start_step).
        started, metacool_file, bins_prefix, planned, stored, clusters = start_step(
            "validation",
            outdir,
            [final_fasta_dir],
            map,
            inplace,
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
        seed_random(seed)
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
            outdir,
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

    log_end(started, "Final contact map written", metacool_file)


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------


@cli.command(short_help="Run network, partition and validation, then mge and host with --mges, in one command line.")
@input_argument(
    help=input_help(
        "named <prefix>.metacool after --prefix, or after the FASTA file of the assembly (<assembly>.metacool)",
        "the clusters of the run (the map of the clusters)",
    )
)
@network_options("-m", CHECKED_ASSEMBLY_HELP)
@min_quality_option()
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
@_option(
    "--mges",
    metavar="FILE",
    panel="MGE",
    help="The MGE contigs, one name per line, as metator mge reads them: after the validation, the run goes on with two "
    "more steps, mge, which clusters them as metator mge --reassign does, each MGE cluster a genome of the map, and "
    "host, which finds their hosts as metator host does. Without it, the run has no steps mge and host.",
)
@_option(
    "--mge-threshold-bin",
    type=click.FLOAT,
    default=0.8,
    show_default=True,
    metavar="FLOAT",
    panel="MGE",
    help="The score above which two MGE contigs, or two MGE clusters, are merged. The step mge, as --threshold-bin of "
    "metator mge.",
)
@_option(
    "--mge-min-share",
    type=click.FloatRange(0, 100),
    default=mth.MIN_SHARE,
    show_default=True,
    metavar="FLOAT",
    panel="Hosts",
    help=f"{MIN_SHARE_HELP} The step host, as --min-share of metator host.",
)
@_option(
    "--mge-min-contigs",
    type=click.IntRange(min=1),
    default=mth.MIN_CONTIGS,
    show_default=True,
    metavar="INT",
    panel="Hosts",
    help=f"{MIN_CONTIGS_HELP} The step host, as --min-contigs of metator host.",
)
@outdir_option(
    help="The output directory, of the map of the clusters, the log, and the outputs of the run: fasta (mags_with-mges, "
    "mags, and mges with --mges), cool, tables and plots. Default: current directory."
)
@prefix_option(
    "the clusters",
    " The map that the run makes from pairs or BAM files is STR.metacool (by default, named after the assembly).",
)
@force_option(FORCE_PIPELINE_HELP)
@tmpdir_option(
    help=TMPDIR_TABLES_HELP + " From BAM files, their pairs too (about 25 bytes per pair, plus the names of its two "
    "contigs), which the run makes its map from."
)
@threads_option(
    help="Number of parallel threads allocated for the partition and the validation, and for reading the BAM files and "
    "binning the pairs that the run makes its map from."
)
@no_clean_up_option(help=NO_CLEAN_UP_TABLES_HELP)
@inplace_option(help=RUN_INPLACE_HELP)
@seed_option(
    "Seed of the random numbers of Louvain or Leiden, set once before the partition and again before the recursive "
    "partition of the validation, their iterations drawing other numbers one after the other: from the same map, the "
    "same seed gives the same clusters, those of partition --seed S then validation --seed S. By default, each run "
    "draws other random numbers."
)
def pipeline(
    input,  # the arguments INPUT...
    assembly,
    depth,
    enzyme,
    normalization,
    min_quality,
    algorithm,
    iterations,
    overlap,
    res_param,
    size,
    rec_iter,
    rec_overlap,
    seed,
    mges,
    mge_threshold_bin,
    mge_min_share,
    mge_min_contigs,
    outdir,
    prefix,
    force,
    tmpdir,
    threads,
    no_clean_up,
    inplace,
):
    """Launch the full metaTOR pipeline.

    Cluster the contigs of the assembly from their metaHiC contacts, and,
    with --mges, the contigs of its mobile genetic elements, and find the
    MAGs that host them.

    [cyan]metator pipeline INPUT... -a ASSEMBLY[/] starts from a .metacool
    contact map of the contigs, made from the pairs of the run (by
    [cyan]hicstuff pipeline -M metacool[/] or
    [cyan]metacooler create pairs[/]), or makes it first from the pairs
    files or the BAM files of the run: one bin per contig, as metacooler
    create pairs does, in the output directory, named after --prefix or
    after the assembly. The run reads its network from that map, and the
    restriction sites and the coverage of the contigs from its
    Restriction_site and Shotgun_coverage columns, unless --enzyme and
    --depth are given. It writes the map of its clusters to the output
    directory, under the name of the map, which it leaves as it was
    (--inplace updates a map given in place); the map that it makes, it
    updates in place.

    The network counts every pair, as a contact map does. A validation step
    decontaminates the clusters when it's necessary; it is the critical step
    for memory usage (~40G).

    The run does what [cyan]metator network INPUT...[/], partition MAP and
    validation MAP do one by one, and writes the map that they write, once,
    at the end of the validation: the FASTA files of the final clusters, with
    the MGE contigs that they hold, go to fasta/mags_with-mges, and their
    plots to plots (as PDF files); the tables of the validation stay in its
    temporary directory, which -N keeps. Without --mges, the run has no
    steps mge and host, and says so.

    With --mges, the list of the MGE contigs, the run goes on with two more
    steps on the map of its clusters, which it updates in place, as
    [cyan]metator mge --reassign[/] and [cyan]metator host[/] do one by one:
    mge clusters the MGE contigs, makes each MGE cluster a genome of the map,
    its contigs out of the final clusters that held them, and writes them to
    fasta/mges; host finds the MAGs that host the MGE clusters, and writes
    tables/mge_hosts.txt and plots/mge_hosts.pdf. A run that stops in one of
    them leaves the map of the clusters of the validation, which metator mge
    --reassign and metator host complete. What they would refuse is refused
    before any work: a map without the contact statistics of the reads of
    its contigs at 1001 bp (the map that the run makes from pairs or BAM
    files has them), a list none of whose contigs is in the map (from pairs
    or BAM files, in the assembly, before the run makes its map), and the
    MGE files of an earlier version in the output directory.

    Its last step, export, writes the outputs of the run from its map, as
    [cyan]metator export MAP -a ASSEMBLY[/] does, in place of those that
    its steps wrote: fasta/mags_with-mges, fasta/mags (the FASTA files of its
    MAGs, its genomes but the MGE genomes, each as its contigs, MAG.fa, and
    stitched, MAG.stitched.fa, its contigs joined by 200 N), fasta/mges with
    --mges; cool, a cool file per MAG, with the MGE clusters that it hosts,
    as [cyan]metacooler extract[/] writes it, balanced as
    [cyan]cooler balance[/] balances it; tables, the tables of
    [cyan]metator dump[/] (and mge_hosts.txt with --mges); plots; and, last,
    the report of the run, <name>.report.txt, <name>.report.html and
    <name>.report_mqc.yaml (<name>: the map without .metacool), as
    [cyan]metator report[/] writes it.
    """
    # The directories that the run makes anew in its output directory (see
    # metator.outputs): fasta/mges with --mges only.
    out_names = [name for name in PIPELINE_DIRS if name != mto.MGES_DIR or mges is not None]
    # Private temporary directory inside --tmpdir, out of those directories,
    # and the lock of the map the run writes, held until its end.
    refuse_tmpdir_in_the_clusters(tmpdir, outdir, out_names)
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir, ExitStack() as held:
        # Defined the output directory and output file names.
        outdir = outdir or "."
        make_outdir(outdir)

        # Enable file logging
        log_file = log_file_name(outdir, run_name(input, prefix, assembly))
        started = generate_log_header(log_file, outdir=outdir, input=shlex.join(input))

        # Define variable
        overlapping_parameter = overlap / 100
        recursive_overlapping_parameter = rec_overlap / 100

        # Check the arguments, and the assembly below, before replacing the
        # clusters of an earlier run (--force): a refused run leaves them as
        # they were. The inputs: the map that the run starts from, or the
        # pairs or BAM files that it makes it from, in its output directory,
        # named after --prefix (checked first) or the assembly.
        inputs = run_inputs(input, given_on_the_command_line("min_quality"))
        bins_prefix = bin_prefix(prefix)
        # The run holds the map's lock from now on: a second run on the same
        # map is refused before any work. It reads its map under that lock
        # (planned; once it has made it, from pairs or BAM files). With
        # --mges, the run replaces the round mge of the map too, in its step
        # mge, and the map that it writes is refused to be the list.
        map, metacool_file, planned = run_map(
            "pipeline",
            inputs,
            held,
            mti.map_name(prefix, assembly),
            assembly,
            depth,
            enzyme,
            normalization,
            outdir,
            inplace,
            force,
            others=[mges],
            rounds=[mtcl.MGE] if mges is not None else None,
        )

        # The outputs of an earlier run (the directories that the run makes
        # anew, its tables, its plots and its report): refused now without
        # --force, before the assembly is read, and replaced only once it is
        # checked. What a run of an earlier version wrote in the output
        # directory is refused, with or without it.
        final_fasta_dir = join(outdir, mto.MAGS_WITH_MGES_DIR)
        out_dirs = [join(outdir, name) for name in out_names]
        refuse_no_directory(join(outdir, mto.FASTA_DIR), "the FASTA files of its clusters")
        for path in out_dirs:
            refuse_output_dir(path, force)
        tables = refuse_earlier_tables(join(outdir, mto.TABLES_DIR), force)
        refuse_earlier_outputs(outdir, keep=(map, metacool_file))
        plots = refuse_earlier_plots(outdir, force)
        reports = refuse_earlier_report(outdir, mtr.sample_name(metacool_file), force)
        # With --mges, what the steps mge and host would refuse, refused now
        # too, before hours of work (see check_mges); the MGE contigs.
        mges_list = check_mges(map, planned, mges, outdir, keep=(map, metacool_file)) if mges is not None else None

        # The network, as metator network makes it: the contig data,
        # measured on the assembly and read from the map, before the clusters
        # of an earlier run are replaced (a map whose columns cannot give what
        # the normalization needs is refused once the assembly is read), and
        # the contacts of its contigs, weighed with the normalization.
        check_assembly(assembly)
        # From pairs or BAM files, a step more, first: the map of the contigs,
        # made in the output directory from the contig data measured on the
        # assembly, which the run then updates in place; with --mges, two
        # more after the validation: mge and host; and export, last. With
        # --mges, the MGE contigs are checked against the contigs of the
        # assembly, which the map has, before the run makes it.
        step = 0 if inputs.kind == mti.MAP else 1
        stages = step + (5 if mges_list is not None else 3) + 1
        contig_data = None
        if planned is None:
            log_stage(1, "map", f"the contact map of the contigs, from {inputs.describe()}", stages)
            contig_data = mtn.create_contig_data(assembly, depth, enzyme)
            if mges_list is not None:
                check_mges_of_the_assembly(assembly, contig_data, mges_list)
            planned = make_run_map("pipeline", inputs, held, metacool_file, contig_data, assembly, min_quality, threads, tmp_dir)
        log_stage(step + 1, "network", f"the contigs of {map} and their contacts", stages)
        network, contacts = run_network(map, planned, assembly, depth, enzyme, normalization, tmp_dir, contig_data)
        for path in out_dirs:
            make_output_dir(path, force)
        remove_earlier_plots(plots)
        remove_earlier_tables(tables)
        remove_earlier_report(reports)
        edges = mtn.map_network(map, mtn.network_weights(network.table, normalization), contacts)

        # Partition the network, as metator partition does.
        log_stage(step + 2, "partition", "the clusters of the network", stages)
        seed_random(seed)
        clusters = mtp.partition(algorithm, edges, iterations, overlapping_parameter, res_param, threads, bins_prefix)

        # Validate the partition, as metator validation does, from the
        # assembly decompressed once for the run when it is compressed: the
        # FASTA files of the clusters are written from it, and those of the
        # MGE clusters in the step mge.
        log_stage(
            step + 3,
            "validation",
            "miComplete checks the clusters, and the contaminated ones are partitioned again",
            stages,
        )
        plain = mtv.plain_assembly(assembly, tmp_dir)
        # Seeded again, as validation --seed seeds its recursive partition.
        seed_random(seed)
        final = run_validation(
            map,
            plain,
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
            outdir,
            tmp_dir,
            normalization,
        )

        # The map of the clusters, written once: what the three steps store,
        # as they store it one by one (see metator.metacool.write_map). With
        # --mges, the steps mge and host store the round mge and the hosts of
        # the MGE clusters in it again.
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
            mge_follows=mges_list is not None,
        )

        if mges_list is None:
            logger.info(
                "No list of MGE contigs (--mges): the run has no steps mge and host, which cluster the MGE contigs and "
                "find their hosts: give --mges to run them."
            )
        else:
            run_mge_steps(
                metacool_file,
                plain,
                mges,
                mges_list,
                mge_threshold_bin,
                mge_min_share,
                mge_min_contigs,
                outdir,
                tmp_dir,
                step + 3,
                stages,
            )

        # The outputs of the run, from its map as its steps left it, as
        # metator export writes them: in place of those that its steps wrote,
        # in the directories that it made anew (force).
        log_stage(stages, "export", "the outputs of the run, from its map", stages)
        run_export(metacool_file, plain, outdir, tmp_dir, force=True)

    log_end(started, "Final contact map written", metacool_file)


def check_mges(input_map, planned, mges, outdir, keep):
    """Refuse, before any work, what the steps mge and host of a pipeline run
    with --mges (metator mge --reassign and metator host, see run_mge_steps)
    would refuse once its clusters are made: the MGE files that an earlier
    version wrote in the output directory, outdir (see refuse_earlier_outputs,
    with MGE_EARLIER_OUTPUTS and HOST_EARLIER_OUTPUTS; keep: the maps of the
    run), as mge and host refuse them; a list of MGE contigs, mges, that
    cannot be read (see read_mges); and a map of the contigs (input_map, as
    the run read it: planned) without the contact statistics that mge reads,
    or that holds none of the MGE contigs (see metator.mge.check_contig_map).
    The reason is logged as an error. The map that a run makes from pairs or
    BAM files (planned is None until then) has those statistics, and the
    contigs of the assembly, which check_mges_of_the_assembly checks the MGE
    contigs against, before the run makes it.

    Returns
    -------
    list of str
        The MGE contigs (see metator.io.import_mges_contigs)."""
    refuse_earlier_outputs(outdir, MGE_EARLIER_OUTPUTS + HOST_EARLIER_OUTPUTS, keep=keep)
    mges_list = read_mges(mges)
    if planned is not None:
        try:
            mtm.check_contig_map(input_map, planned, mges_list)
        except mtmc.MetacoolError as error:
            logger.error(str(error))
            raise
    return mges_list


def check_mges_of_the_assembly(assembly, contig_data, mges_list):
    """Refuse, before a pipeline run with --mges makes its map from pairs or BAM
    files, a list of MGE contigs (mges_list) none of which the map would hold:
    its contigs are those of the assembly (contig_data, as measured on it),
    but those of length 0 (see metator.inputs.make_map), and its step mge
    would have no MGE contig to cluster. The reason is logged as an error."""
    if any(contig_data[name]["length"] > 0 for name in mges_list if name in contig_data):
        return
    message = f"None of the MGE contigs is a contig of {assembly}: there is no MGE contig to cluster."
    logger.error(message)
    raise mtmc.MetacoolError(message)


def run_mge_steps(metacool, assembly, mges, mges_list, threshold_bin, min_share, min_contigs, outdir, tmp_dir, step=3, stages=5):
    """The steps mge and host of a pipeline run with --mges, the last two of its
    stages, after its step step (the validation), on the map of its clusters,
    metacool, which they update in place, each as the run reads it once the
    step before wrote it: mge clusters the MGE contigs (mges_list, the
    contigs of the list mges) of the map, as metator mge --reassign does (see
    run_mge: each MGE cluster a genome of the map, its contigs out of the
    final clusters that held them), and host finds the MAGs that host the
    MGE clusters, as metator host does (see run_host): their outputs go to
    outdir, as theirs (fasta/mges; tables/mge_hosts.txt and
    plots/mge_hosts.pdf). A run that stops in one of them says what the map
    holds then, which metator mge and metator host complete, then metator
    export, the last step of the run (see run_export)."""
    name = "mge"
    try:
        step += 1
        log_stage(step, name, f"the MGE contigs of {mges}, clustered from their contacts, each cluster a genome", stages)
        run_mge(metacool, metacool, open_input_map(metacool), assembly, mges_list, threshold_bin, outdir, tmp_dir, reassign=True)
        step, name = step + 1, "host"
        log_stage(step, name, "the MAGs that host the MGE clusters, from the network of the map", stages)
        run_host(metacool, metacool, open_input_map(metacool), min_share, min_contigs, outdir)
    except BaseException:
        held, finish = (
            ("the clusters of the run", "metator mge --reassign then metator host")
            if name == "mge"
            else ("the clusters of the run and their MGE clusters", "metator host")
        )
        logger.info(
            f"The run stopped in its step {step}, {name}: {metacool} holds {held}; run {finish} on it to complete the "
            "map, then metator export to write the outputs of the run from it."
        )
        raise


def write_map_outputs(metacool, assembly, outdir, tmp_dir):
    """Write in the output directory of a run, outdir, from the map of its
    clusters, metacool (see metator.outputs), part of what metator export
    writes (see run_export), and the pipeline in its step export: fasta/mags, the FASTA files of
    its MAGs, each as its contigs and stitched (see
    metator.outputs.write_fasta_files, from the assembly of the run, plain: see
    metator.validation.plain_assembly; none without it); cool, the cool file of
    each, with the MGE clusters that it hosts, balanced (see
    metator.outputs.write_cool_files); and tables, the tables of metator dump
    (see metator.dump.dump). The directories are those that the run made anew
    (fasta/mags and cool), and tables holds no table of an earlier run (see
    refuse_earlier_tables). tmp_dir: the temporary directory of the run.

    Raises
    ------
    metator.metacool.MetacoolError
        If an output cannot be written."""
    with mtmc.read_errors(metacool):
        mc = mtmc.open_map(metacool)
        mags = mto.map_mags(mc)
        hosted = mto.hosted_mges(mc)
    if assembly is not None:
        mags_dir = join(outdir, mto.MAGS_DIR)
        mto.write_fasta_files(metacool, assembly, mags_dir, genomes=mags, temp_dir=tmp_dir)
        logger.info(f"The FASTA files of the {len(mags)} MAG(s) of the run, each as its contigs and stitched: {mags_dir}")
    cool_dir = join(outdir, mto.COOL_DIR)
    mto.write_cool_files(metacool, cool_dir, mags, hosted)
    with_mges = sum(1 for mag in mags if hosted.get(mag))
    logger.info(
        f"The cool files of the {len(mags)} MAG(s) of the run, {with_mges} of them with the MGE clusters that they "
        f"host, balanced: {cool_dir}"
    )
    tables_dir = join(outdir, mto.TABLES_DIR)
    os.makedirs(tables_dir, exist_ok=True)
    mtd.dump(metacool, tables_dir, force=True)


# ---------------------------------------------------------------------------
# mge
# ---------------------------------------------------------------------------


@cli.command(short_help="Cluster the MGE contigs of a .metacool map into MGE MAGs, and store them in the map.")
@map_argument(
    required=True,
    help="The .metacool contact map of the run, which metator network stored the network of its contigs in (that of "
    "network MAP, pipeline MAP or validation MAP), or an entry of a .mmetacool collection "
    "(file.mmetacool::entries/NAME): the contacts of the MGE contigs come from it. It must hold the contact statistics "
    "of the reads of its contigs at 1001 bp (made from the pairs by hicstuff pipeline -M metacool or metacooler create "
    "pairs, with metacooler 0.4.0 or later). " + written_map_help("the MGE clusters in its clustering rounds"),
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
    "mge_random, which leaves the round mge and the genomes of the map as they were, and written to "
    "fasta/mges_random.",
)
@_option(
    "--reassign",
    is_flag=True,
    panel="Clustering",
    help="Make each MGE cluster a genome of the map, named after it, with its GC content, hits and clustering score: "
    "its contigs leave the genomes that held them (the final clusters, MAGs, or none), which keep the values of their "
    "validation. The column mge of the genomes table says which genomes are MGE genomes (True). Refused with "
    "--random.",
)
@seed_option(
    "Seed of the random clustering of --random: from the same map, the same seed gives the same random clustering. "
    "By default, each run draws another one.",
    panel="Clustering",
)
@outdir_option(
    help="The output directory, of the map of the run and of fasta/mges (fasta/mges_random with --random), the FASTA "
    "files of the MGE clusters. Default: current directory."
)
@force_option(
    "Replace the map and fasta/mges (fasta/mges_random with --random) of an earlier run in the output directory. A "
    "summary or a FASTA file of the MGE clusters that an earlier version wrote in the output directory "
    "(mges_bin_summary.tsv, mges_cluster_summary.tsv, mges_binned.fa, mges_binned_random.fa) is refused, with or "
    "without -F: it is never removed."
)
@tmpdir_option(
    help="Path to temporary directory: where metacooler copies a compressed assembly, plain, while it reads it. mge "
    "writes no temporary file of its own."
)
@no_clean_up_option(help="Keep the temporary directory of the run: mge leaves nothing in it.")
@inplace_option(help=INPLACE_HELP)
def mge(
    map,  # the argument MAP
    assembly,
    mges,
    threshold_bin,
    random,
    reassign,
    seed,
    outdir,
    force,
    tmpdir,
    no_clean_up,
    inplace,
):
    """Cluster MGE contigs.

    Cluster the MGE contigs of a map into MGE MAGs, from the Hi-C contacts
    between them: those that the map holds between two contigs, and those
    within each contig between reads more than 1 kb apart, which the
    contact statistics of the map count (its QC block, made from the pairs
    at 1001 bp). The GC content and hits of each contig are those that
    metator network stored in the map.

    Writes the FASTA files of the MGE clusters in fasta/mges (fasta/mges_random
    with --random), two for each MGE cluster, named after it: MGE.fa, its
    contigs, and MGE.stitched.fa, one record, its contigs joined by 200 "N".

    The MGE clusters are stored in the map, MAP: as the round mge of its
    clustering rounds (its /clusters group), with their contigs and their
    length, GC content, hits and clustering score, in place of an earlier
    round mge, after the other rounds, which are kept (--random stores the
    round mge_random, and leaves the round mge as it was). Its genomes and
    contigs tables are left as they were, but with --reassign: each MGE
    cluster is a genome of the map then, named after it, with its GC content,
    hits and clustering score, and its contigs leave the genomes that held
    them (MAGs, or none), which keep the values of their validation; the
    column mge of the genomes table says which genomes are MGE genomes. A run
    that replaces the round mge undoes the MGE genomes of an earlier run
    first: their contigs go back to their final clusters, or to none. The
    run writes the map to the output directory, under the name of MAP, which
    it leaves as it was (--inplace updates MAP in place, and so does a run in
    the directory of MAP), and locks it while it writes it. A run that
    replaces the clusters of the map (partition, validation, pipeline) drops
    the rounds of mge, with a warning: run mge again then (validation and
    pipeline replace the MGE genomes with their final clusters). A run that
    replaces the round mge (not --random) drops the hosts that metator host
    stored in the map, with a warning.
    """
    # MGE genomes are made of the clusters of the round mge, not of a random
    # clustering.
    if random and reassign:
        raise click.BadParameter(
            "--reassign makes genomes of the MGE clusters of the round mge, and --random stores a random clustering "
            "as the round mge_random, which leaves them as they were: give one or the other.",
            param_hint="'--reassign'",
        )
    # The FASTA files of the MGE clusters, which the run makes anew: no
    # --tmpdir in them.
    fasta_dir = mto.MGES_RANDOM_DIR if random else mto.MGES_DIR
    refuse_tmpdir_in_the_clusters(tmpdir, outdir, [fasta_dir], earlier=MGE_EARLIER_OUTPUTS)
    # Generate log
    log_file = log_file_name("", mtr.sample_name(map))
    outdir = outdir or "."
    started = generate_log_header(log_file, outdir=outdir)

    # Private temporary directory inside --tmpdir, and the lock of the map
    # the run writes, held until its end.
    with private_tmp_dir(tmpdir, clean_up=not no_clean_up) as tmp_dir, ExitStack() as held:
        # The map that the run stores its MGE clusters in, checked and
        # locked.
        name = mtcl.MGE_RANDOM if random else mtcl.MGE
        metacool_file, planned = step_map("mge", map, outdir, inplace, force, held, inputs=[assembly, mges], rounds=[name])
        # The assembly that metacooler writes the FASTA files of the MGE
        # clusters from, refused as the network, the validation and the
        # pipeline refuse it.
        check_assembly(assembly)
        # The summary and the FASTA file of an earlier version, which would
        # pass for this run's, and the FASTA files of an earlier run, which -F
        # replaces.
        refuse_earlier_outputs(outdir, MGE_EARLIER_OUTPUTS, keep=(map, metacool_file))
        refuse_no_directory(join(outdir, mto.FASTA_DIR), "the FASTA files of its MGE clusters")
        refuse_output_dir(join(outdir, fasta_dir), force)

        # The MGE contigs, their clusters, their FASTA files and the map (see
        # run_mge), the random clustering of --random seeded with --seed.
        seed_random(seed)
        run_mge(map, metacool_file, planned, assembly, read_mges(mges), threshold_bin, outdir, tmp_dir, random, reassign)

    log_end(started, "Contact map with the MGE clusters written", metacool_file)


def read_mges(mges):
    """The MGE contigs of a mge or pipeline run, from their list, mges (see
    metator.io.import_mges_contigs), read before any work: a list that cannot
    be read (not there, a directory, not text) is refused, the reason logged
    as an error, then ValueError is raised."""
    try:
        return mio.import_mges_contigs(mges)
    except (OSError, UnicodeDecodeError) as error:
        reason = error.strerror if isinstance(error, OSError) and error.strerror else str(error)
        message = f"Cannot read the list of MGE contigs {mges} ({reason}): give a text file, one contig name per line."
        logger.error(message)
        raise ValueError(message) from error


def run_mge(metacool, outfile, planned, assembly, mges_list, threshold_bin, outdir, tmp_dir, random=False, reassign=False):
    """Cluster the MGE contigs of a map, as metator mge does, and the step mge
    of pipeline --mges (with reassign) on the map of its clusters.

    The MGE contigs (mges_list) of the map metacool are clustered from its
    contacts (see metator.mge.mge_binning), what the map cannot give refused
    before anything is written; each MGE cluster is written to fasta/mges
    (fasta/mges_random with random) in outdir, which it makes anew, from the
    assembly (see metator.mge.write_mges_fasta: a compressed one is copied
    plain in tmp_dir as it is read); and the MGE
    clusters are stored in the map of the run, outfile, written from
    metacool as the run read it (planned): the round mge (or mge_random), in
    place of the round of that name that the map held, after the other
    rounds, which are kept as they were, and so is everything else but the
    MGE genomes (made of the round mge with reassign, undone without); the
    run is recorded with the round, and under its name in the metator
    metadata of the map (see metator.metacool.write_map)."""
    name = mtcl.MGE_RANDOM if random else mtcl.MGE
    try:
        summary, mge_bins = mtm.mge_binning(metacool, mges_list, threshold_bin=threshold_bin, random=random)
        # The FASTA files of the MGE clusters, in a directory made anew: the
        # files of an earlier run (which -F replaces) are gone.
        fasta_dir = join(outdir, mto.MGES_RANDOM_DIR if random else mto.MGES_DIR)
        make_output_dir(fasta_dir, force=True)
        mtm.write_mges_fasta(metacool, assembly, mge_bins, fasta_dir, temp_dir=tmp_dir)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    write_map(
        metacool,
        outfile,
        planned,
        "the map with its MGE clusters",
        mge=(mtcl.mge_round(summary, name), mge_entry(name)),
        reassign=reassign,
    )


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
# host
# ---------------------------------------------------------------------------


@cli.command(short_help="Find the MAGs that host the MGE clusters of a .metacool map, and store them in the map.")
@map_argument(
    required=True,
    help="The .metacool contact map of the run, or an entry of a .mmetacool collection "
    "(file.mmetacool::entries/NAME): the map that metator mge stored its MGE clusters in (the round mge), whose "
    "genomes are the MAGs (the final clusters of metator validation or pipeline, or the bins of another binner, made "
    "its genomes by metacooler reassign), and which holds the network of its contigs (metator network). "
    + written_map_help("the hosts of its MGE clusters"),
)
@_option(
    "-s",
    "--min-share",
    type=click.FloatRange(0, 100),
    default=mth.MIN_SHARE,
    show_default=True,
    metavar="FLOAT",
    panel="Hosts",
    help=MIN_SHARE_HELP,
)
@_option(
    "-c",
    "--min-contigs",
    type=click.IntRange(min=1),
    default=mth.MIN_CONTIGS,
    show_default=True,
    metavar="INT",
    panel="Hosts",
    help=MIN_CONTIGS_HELP,
)
@outdir_option(
    help=f"The output directory, of the map of the run, {mto.TABLES_DIR}/{mth.HOSTS_FILE}, {mtf.PLOT_DIR}/{mth.PLOT_FILE} "
    "and the log. Default: current directory."
)
@force_option(
    f"Replace the map of an earlier run in the output directory ({mto.TABLES_DIR}/{mth.HOSTS_FILE} and "
    f"{mtf.PLOT_DIR}/{mth.PLOT_FILE} are replaced without it). The table and the plot of the hosts that an earlier "
    "version wrote in the output directory (mge_hosts.txt, mge_hosts.png) are refused, with or without -F: they are "
    "never removed."
)
@inplace_option(help=INPLACE_HELP)
def host(
    map,  # the argument MAP
    min_share,
    min_contigs,
    outdir,
    force,
    inplace,
):
    """Find the MAGs that host the MGE clusters.

    [cyan]metator host MAP[/] finds the hosts of the MGE clusters that
    [cyan]metator mge MAP[/] stored in the map (its round mge) among its
    MAGs: its genomes, the final clusters of [cyan]metator validation[/] or
    [cyan]pipeline[/], or the bins of another binner made its genomes by
    [cyan]metacooler reassign[/]; not __unbinned__, nor a genome named after
    an MGE cluster. The contigs of an MGE cluster count for it only, and the
    MAGs that hold some of them are taken without them.

    The signal of an MGE cluster with a MAG is the sum of the weights of the
    edges of the network between their contigs, as [cyan]metator network[/]
    stored it in the map, with its normalization. A MAG hosts the MGE cluster
    when its signal is at least --min-share percent of the signal of the MGE
    cluster with all the MAGs, and the MGE cluster is in contact with at
    least --min-contigs of its contigs: an MGE cluster can have one host,
    several, or none.

    The run writes every MGE cluster and MAG in contact to
    tables/mge_hosts.txt, with their signal, the number of contigs of the MAG
    in contact, the share of the signal of the MGE cluster with the MAGs
    (share) and with every other contig (share_of_all), and whether the MAG
    is a host; and the distribution of the signal and share of the hosts to
    plots/mge_hosts.pdf. It
    stores the hosts in the map: in the round mge, the number of hosts of
    each MGE cluster, their names and their shares (n_hosts, hosts,
    host_shares); in the genomes table, the number of MGE clusters that each
    MAG hosts and their names (n_mges, mges); and the run under the host
    entry of the metator metadata of the map. They replace those of an
    earlier run, and a run that replaces the round mge (mge), the clusters
    (partition, validation, pipeline) or the network (network, with another
    one) drops them, with a warning: run host again then. The run writes the
    map to the output directory, under the name of MAP, which it leaves as it
    was (--inplace updates MAP in place, and so does a run in the directory
    of MAP), and locks it while it writes it.
    """
    outdir = outdir or "."
    make_outdir(outdir)

    # Enable file logging
    log_file = log_file_name(outdir, mtr.sample_name(map))
    started = generate_log_header(log_file, outdir=outdir)

    # The lock of the map the run writes, held until its end.
    with ExitStack() as held:
        # The map that the run stores the hosts in, checked and locked; the
        # table and the plot of an earlier version, which would pass for this
        # run's; then the hosts, their table, their plot and the map (see
        # run_host).
        metacool_file, planned = step_map("host", map, outdir, inplace, force, held)
        refuse_earlier_outputs(outdir, HOST_EARLIER_OUTPUTS, keep=(map, metacool_file))
        refuse_no_directory(join(outdir, mto.TABLES_DIR), "the table of the hosts")
        refuse_no_directory(join(outdir, mtf.PLOT_DIR), "the plot of the hosts")
        run_host(map, metacool_file, planned, min_share, min_contigs, outdir)

    log_end(started, "Contact map with the hosts of the MGE clusters written", metacool_file)


def run_host(metacool, outfile, planned, min_share, min_contigs, outdir):
    """Find the MAGs that host the MGE clusters of a map, as metator host does,
    and the step host of pipeline --mges on the map of its clusters.

    The network of the map metacool is read through the map as the run read
    it (planned), and so are its MAGs and MGE clusters: what the map cannot
    give is refused before anything is written (see metator.host.find_hosts).
    Every MGE cluster and MAG in contact is written to tables/mge_hosts.txt in
    outdir, and the plot of the hosts to plots/mge_hosts.pdf (both made if
    missing); the hosts are stored
    in the map of the run, outfile, written from metacool: in the round mge
    and the genomes table, in place of those of an earlier run, and the run
    under the host entry of the metator metadata of the map (see
    metator.metacool.write_map)."""
    stored = stored_network(metacool, planned)
    try:
        hosts = mth.find_hosts(metacool, planned, stored.network, min_share, min_contigs)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    write_hosts_files(hosts, outdir)
    write_map(metacool, outfile, planned, "the map with the hosts of its MGE clusters", hosts=(hosts, host_entry(stored)))


def write_hosts_files(hosts, outdir):
    """Write the files of the hosts of the MGE clusters of a map (hosts, see
    metator.host.find_hosts) in the output directory of a run, outdir: every
    MGE cluster and MAG in contact in tables/mge_hosts.txt, and the plot of the
    hosts in plots/mge_hosts.pdf (both directories made if missing). Returns
    their paths."""
    tables_dir, plot_dir = join(outdir, mto.TABLES_DIR), join(outdir, mtf.PLOT_DIR)
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)
    table, plot = join(tables_dir, mth.HOSTS_FILE), join(plot_dir, mth.PLOT_FILE)
    mth.write_hosts(hosts.pairs, table)
    mtf._drawn(plot, lambda: mth.plot_hosts(hosts.pairs, plot))
    return table, plot


def host_entry(stored):
    """What a map records of the host run that stores the hosts of its MGE
    clusters in it, under the host entry of its metator metadata (see
    metator.metacool.write_map): the version of metaTOR, the command line and
    the parameters of the run, the normalization of the network that it read
    (stored, see metator.network.stored_network) and its digest (see
    metator.network.network_digest), which a network run checks the network
    that it stores against, and what the run read of the map (read): the
    columns of its network in its contigs table and the weights of its
    normalization (see map_read), its round mge, and its genomes, by name."""
    read = {**map_read(stored.table, stored.normalization, (mtcl.MGE,)), "genomes": ["name"]}
    return {
        "version": __version__,
        **run_metadata(),
        "normalization": stored.normalization,
        "read": read,
        mtmc.NETWORK_DIGEST: stored.digest,
    }


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
    network metator network replaced after the run (run it again with the
    normalization of the run), and an output directory that holds tables of
    an earlier dump, without -F. The tables follow the genomes of the map,
    clusters edited after the run (metacooler reassign) included. The tables are replaced all at once: a dump
    that fails leaves them as they were.
    """
    outdir = outdir or "."
    make_outdir(outdir)

    # Enable file logging
    log_file = log_file_name(outdir, mtr.sample_name(map))
    started = generate_log_header(log_file, outdir=outdir)

    # The map, checked before anything is read.
    run_input(map)
    try:
        mtd.dump(map, outdir, force)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise

    log_end(started, "Tables written", outdir)


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------

#: The results of an earlier version that metator export refuses in its output
#: directory, with or without -F: those of validation, pipeline, mge and host.
EXPORT_EARLIER_OUTPUTS = EARLIER_OUTPUTS + MGE_EARLIER_OUTPUTS + HOST_EARLIER_OUTPUTS

#: The files of the report of a run, as the help names them (see
#: metator.report.report_files): <name> is its map without .metacool.
REPORT_NAMES = ", ".join(mtr.report_files("<name>"))

#: What to do about a map that another run writes while metator export reads it.
EXPORT_AGAIN = "Run metator export again, with -F for what it wrote, once the run that writes the map has ended."

FORCE_EXPORT_HELP = (
    "Replace the folders that the run makes anew (fasta/mags_with-mges, fasta/mags, fasta/mges, fasta/mges_random, "
    "cool), the tables of tables/, the plots of the clusters of plots/ (whose other files are left as they are) and "
    f"the report ({REPORT_NAMES}) that an earlier run wrote in the output directory; {mto.TABLES_DIR}/{mth.HOSTS_FILE} and "
    f"{mtf.PLOT_DIR}/{mth.PLOT_FILE} are replaced without it, as metator host replaces them. The outputs of an earlier "
    "version (final_clusters, plot, mges_binned.fa, mge_hosts.txt...) are refused, with or without it."
)

#: The FASTA files of the clusters of the rounds of a map that metator export
#: writes, by directory (see write_round_fasta).
ROUND_FASTA = {
    mto.MAGS_WITH_MGES_DIR: "final cluster(s), with the MGE contigs that they hold",
    mto.MGES_DIR: "MGE cluster(s), each as its contigs and stitched",
    mto.MGES_RANDOM_DIR: "random MGE cluster(s) of metator mge --random, each as its contigs and stitched",
}


def map_state(metacool):
    """What changes when the file of the map metacool (that of its collection, for
    an entry) is replaced or written: its device and inode, its size and its
    time of modification. metator export reads the map by path, again and
    again, and refuses at its end what it wrote from a map that changed since
    it started."""
    status = os.stat(mtmc.map_file(metacool))
    return status.st_dev, status.st_ino, status.st_size, status.st_mtime_ns


def refuse_an_input_in_the_outputs(inputs, outdir, dirs):
    """Refuse, before any work, an input of the run (inputs: the file of its map,
    its assembly, None for none) in a directory of the output directory,
    outdir, that the run makes anew (dirs), which -F replaces with all that it
    holds: the reason is logged as an error, then ValueError is raised."""
    for path in inputs:
        for name in dirs:
            directory = join(outdir, name)
            if path is not None and _within(os.path.realpath(path), os.path.realpath(directory)):
                logger.error(
                    f"{path}, an input of the run, is in {directory}, which the run makes anew: move it away, or give "
                    "the run another output directory with -o."
                )
                raise ValueError


def write_round_fasta(metacool, assembly, outdir, groups, tmp_dir):
    """Write in the output directory of a run, outdir, the FASTA files of the
    clusters of the rounds of its map (metacool), from its assembly (plain: see
    metator.validation.plain_assembly), by directory (groups, see ROUND_FASTA;
    None for a round that the map lacks): the final clusters in
    fasta/mags_with-mges, each as its contigs, as the validation writes them;
    the MGE clusters in fasta/mges (fasta/mges_random), each as its contigs and
    stitched, as metator mge writes them (see metator.outputs.round_groups and
    mge_groups). tmp_dir: the temporary directory of the run."""
    for name, clusters in groups.items():
        if clusters is None:
            continue
        directory = join(outdir, name)
        stitched = name != mto.MAGS_WITH_MGES_DIR
        mto.write_fasta_files(metacool, assembly, directory, groups=clusters, temp_dir=tmp_dir, stitched=stitched)
        logger.info(f"The FASTA files of the {len(clusters)} {ROUND_FASTA[name]}: {directory}")


def run_export(metacool, assembly, outdir, tmp_dir, force):
    """Write the outputs of a run in its output directory, outdir, from the map
    of its clusters, metacool, and its assembly (None: without the FASTA
    files), as metator export does (see export), and the step export of
    pipeline, at its end; tmp_dir: the temporary directory of the run. force
    (-F) replaces the outputs of an earlier run (those that the pipeline made
    anew, in its step export). The reason of a refusal or of a failure is
    logged as an error, then raised."""
    # The map and the assembly, checked before any work.
    run_input(metacool)
    if assembly is not None:
        check_assembly(assembly)
    state = map_state(metacool)
    rounds = [entry.get("name") for entry in open_input_map(metacool).cluster_rounds]

    # The folders that the run makes anew (the FASTA files with the
    # assembly only, those of the MGE clusters when the map holds them),
    # the tables and the plots of an earlier run: refused now without -F,
    # which replaces them once the map is read. What an earlier version
    # wrote in the output directory is refused, with or without it, and so
    # are the inputs of the run in a folder that it makes anew.
    fasta = [mto.MAGS_WITH_MGES_DIR, mto.MAGS_DIR]
    fasta += [
        name for name, round_name in [(mto.MGES_DIR, mtcl.MGE), (mto.MGES_RANDOM_DIR, mtcl.MGE_RANDOM)] if round_name in rounds
    ]
    made = (fasta if assembly is not None else []) + [mto.COOL_DIR]
    refuse_an_input_in_the_outputs([mtmc.map_file(metacool), assembly], outdir, made)
    refuse_earlier_outputs(outdir, EXPORT_EARLIER_OUTPUTS, keep=(metacool,))
    if assembly is not None:
        refuse_no_directory(join(outdir, mto.FASTA_DIR), "the FASTA files of the run")
    for name in made:
        refuse_output_dir(join(outdir, name), force)
    tables = refuse_earlier_tables(join(outdir, mto.TABLES_DIR), force)
    plots = refuse_earlier_plots(outdir, force)
    reports = refuse_earlier_report(outdir, mtr.sample_name(metacool), force)

    # The map, read once, as metator dump reads it (and refused as dump
    # refuses it), and what the outputs need of it: the clusters of its
    # rounds final, mge and mge_random, the parameters of the validation
    # (the size threshold of the plots, the assembly), the hosts of its MGE
    # clusters, found again from its network with the thresholds of the host
    # run, and the report of the run.
    try:
        mc, layout, stored = mtd.read_run(metacool, EXPORT_AGAIN)
        with mtmc._metacooler_logging(), mtmc.read_errors(metacool, EXPORT_AGAIN):
            groups = {
                mto.MAGS_WITH_MGES_DIR: mto.round_groups(mc, mtcl.FINAL),
                mto.MGES_DIR: mto.mge_groups(mc, mtcl.MGE),
                mto.MGES_RANDOM_DIR: mto.mge_groups(mc, mtcl.MGE_RANDOM),
            }
            parameters = mc.metadata["metator"]["validation"].get("parameters", {})
            contigs, summary = mto.plot_inputs(mc, layout.rows)
            thresholds = mth.recorded_thresholds(mc)
            hosts = mth.find_hosts(metacool, mc, stored.network, *thresholds) if mtmc.holds_hosts(mc) else None
            report = mtr.collect(metacool, mc, layout, stored)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    for name in made:
        make_output_dir(join(outdir, name), force)
    remove_earlier_tables(tables)
    remove_earlier_plots(plots)
    remove_earlier_report(reports)

    # The outputs, from the map, which metacooler reads again by path; the
    # FASTA files from the assembly, decompressed once for the run when it
    # is compressed.
    try:
        with mtmc.read_errors(metacool, EXPORT_AGAIN):
            plain = None if assembly is None else mtv.plain_assembly(assembly, tmp_dir)
            if plain is None:
                logger.info(
                    "No assembly (--assembly): the FASTA files of the run (fasta) are not written. Give the "
                    f"assembly that the map was made from to write them: {parameters.get('assembly')}, as the run "
                    "recorded it."
                )
            else:
                write_round_fasta(metacool, plain, outdir, groups, tmp_dir)
            write_map_outputs(metacool, plain, outdir, tmp_dir)
        mtf.plot_figures(outdir, contigs, summary, parameters.get("size"))
        if hosts is not None:
            table, plot = write_hosts_files(hosts, outdir)
            logger.info(
                f"The hosts of the MGE clusters, found again from the network of the map with --min-share "
                f"{thresholds[0]:g} and --min-contigs {thresholds[1]}, as the host run recorded them: {table}, {plot}"
            )
        # The report last, which links to the other outputs.
        write_report(report, outdir)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        logger.info(
            f"The run stopped as it wrote the outputs of {metacool}, which is complete: metator export writes them "
            "from it (with -F, in place of what this run wrote)."
        )
        raise

    # What was read of the map, by path, again and again: of one state of
    # it.
    if map_state(metacool) != state:
        message = (
            f"{metacool} was replaced or written while metator export read it, by a run that updated it: what it wrote "
            f"may mix two maps. {EXPORT_AGAIN}"
        )
        logger.error(message)
        raise mtmc.MapReplacedError(message)


@cli.command(short_help="Write the outputs of a run again from the map of its clusters: fasta, cool, tables, plots and report.")
@map_argument(
    required=True,
    help="The map of the clusters of a run, which metator pipeline or metator validation wrote, and metator mge and "
    "metator host updated: a .metacool file, or an entry of a .mmetacool collection (file.mmetacool::entries/NAME). "
    "It is read only.",
)
@_option(
    "-a",
    "--assembly",
    metavar="FILE",
    panel="Input",
    help="The assembly that the map was made from: a FASTA file, plain or compressed with gzip. With it, the run writes "
    "the FASTA files too (fasta); without it, it leaves them out.",
)
@outdir_option(
    help="The output directory of the outputs, made if it is not there: fasta (with --assembly), cool, tables and "
    "plots, and the log. Default: current directory."
)
@force_option(FORCE_EXPORT_HELP)
@tmpdir_option(
    help="Temporary directory, where the run writes, in a directory of its own, the decompressed copy of a compressed "
    "assembly. Not in a folder that the run makes anew (fasta/mags, cool...), nor in a directory of the results of "
    "an earlier version (final_clusters...), which the run refuses."
)
def export(
    map,  # the argument MAP
    assembly,
    outdir,
    force,
    tmpdir,
):
    """Write the outputs of a run again from the map of its clusters.

    [cyan]metator export MAP -a ASSEMBLY -o DIR[/] writes, from the map of the
    clusters of a run ([cyan]metator pipeline[/], or [cyan]metator
    validation[/], then [cyan]metator mge[/] and [cyan]metator host[/]), the
    outputs that a pipeline run writes in its output directory:

    - fasta/mags_with-mges: the FASTA file of each final cluster (the round
      final of the map), with the MGE contigs that it holds: its records in
      the order of the map, where the validation wrote them in the order of
      its partition;
    - fasta/mags: each MAG (a genome of the map but the MGE genomes), as its
      contigs and stitched;
    - fasta/mges, fasta/mges_random: each MGE cluster of the round mge
      (mge_random), the same two ways;
    - cool: a cool file per MAG, with the MGE clusters that it hosts,
      balanced;
    - tables: the tables of [cyan]metator dump[/], and mge_hosts.txt, the
      hosts of the MGE clusters, found again from the network of the map with
      the thresholds that the host run recorded;
    - plots: the plots of the final clusters, from the summary that the map
      holds (the values of miComplete, as the validation drew them), and
      mge_hosts.pdf;
    - <name>.report.txt, <name>.report.html, <name>.report_mqc.yaml (<name>:
      MAP without .metacool): the report of the run, as
      [cyan]metator report[/] writes it, last.

    The FASTA files need the assembly that the map was made from
    (--assembly): without it, they are left out. The MGE clusters and their
    hosts are written when the map holds them. The map is read only: what
    [cyan]metator dump[/] refuses is refused before anything is written, and
    a map that another run writes while it is read is refused (run the export
    again). The folders that the run makes anew, the tables, the plots and
    the report of an earlier run in the output directory are refused without
    -F, and replaced with it; the map and the assembly may not be in a folder
    that the run makes anew.
    """
    outdir = outdir or "."
    # A --tmpdir in a folder that the run makes anew, or refuses as the results
    # of an earlier version, refused as a usage error before anything is made.
    refuse_tmpdir_in_the_clusters(tmpdir, outdir, EXPORT_DIRS, EXPORT_EARLIER_OUTPUTS)
    with private_tmp_dir(tmpdir) as tmp_dir:
        make_outdir(outdir)

        # Enable file logging
        log_file = log_file_name(outdir, mtr.sample_name(map))
        started = generate_log_header(log_file, outdir=outdir)

        run_export(map, assembly, outdir, tmp_dir, force)

    log_end(started, "Outputs written", outdir)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

#: What to do about a map that another run writes while metator report reads
#: it.
REPORT_AGAIN = "Run metator report again once the run that writes the map has ended."


def write_report(report, outdir):
    """Write the report of a run (see metator.report.collect) in its output
    directory, outdir (see metator.report.write_report), and log it. Returns
    the paths of the report."""
    paths = mtr.write_report(report, outdir)
    logger.info(f"The report of the run, from its map: {', '.join(paths)}")
    return paths


@cli.command(short_help="Write the report of a run from the map of its clusters: <name>.report.html, .txt and _mqc.yaml.")
@map_argument(
    required=True,
    help="The map of the clusters of a run, which metator pipeline or metator validation wrote, and metator mge and "
    "metator host updated: a .metacool file, or an entry of a .mmetacool collection (file.mmetacool::entries/NAME). "
    "It is read only.",
)
@outdir_option(
    help=f"The output directory of the report, made if it is not there: {REPORT_NAMES} and the log. Default: current directory."
)
@force_option(f"Replace the report of an earlier run ({REPORT_NAMES}) in the output directory.")
def report(
    map,  # the argument MAP
    outdir,
    force,
):
    """Write the report of a run from the map of its clusters.

    [cyan]metator report MAP -o DIR[/] digests what the map of the clusters
    of a run holds ([cyan]metator pipeline[/], or [cyan]metator
    validation[/], then [cyan]metator mge[/] and [cyan]metator host[/]): the
    run and its parameters, the contigs and where they went, the contacts of
    the reads, the clustering rounds, the final clusters by quality class
    (from the values of miComplete, as the log of the run counts them), the
    MAGs, the MGE clusters and their hosts. It writes, named after the map
    (<name>: MAP without .metacool):

    - <name>.report.txt: a text file, in sections, a key and its value per line,
      then tables, a row per line, in a format that MultiQC, or a script,
      parses (its version, format_version, in its section report);
    - <name>.report.html: a page to read the run with, its numbers, its tables,
      which sort and filter, and interactive plots, which plotly.js draws
      once the page has loaded it from cdn.plot.ly (it needs an internet
      connection for them), with links to the outputs of the run beside it;
    - <name>.report_mqc.yaml: the custom content of MultiQC, a row of its General
      Statistics table for the run, named after its map, which
      [cyan]multiqc[/] aggregates over the runs whose reports it finds.

    [cyan]metator export[/] writes it with the other outputs of the run, and
    so does [cyan]metator pipeline[/], in its step export. The map is read
    only: what [cyan]metator dump[/] refuses is refused, before anything is
    written, and so is a map that another run writes while it is read (run
    the report again). The report of an earlier run in the output directory
    is refused without -F, and replaced with it.
    """
    outdir = outdir or "."
    make_outdir(outdir)

    # Enable file logging
    log_file = log_file_name(outdir, mtr.sample_name(map))
    started = generate_log_header(log_file, outdir=outdir)

    # The map, checked before anything is read, and the report of an earlier
    # run, refused without -F; then the map, read as metator dump reads it
    # (and refused as dump refuses it), before the report of an earlier run
    # is replaced.
    run_input(map)
    reports = refuse_earlier_report(outdir, mtr.sample_name(map), force)
    try:
        mc, layout, stored = mtd.read_run(map, REPORT_AGAIN)
        with mtmc._metacooler_logging(), mtmc.read_errors(map, REPORT_AGAIN):
            collected = mtr.collect(map, mc, layout, stored)
    except mtmc.MetacoolError as error:
        logger.error(str(error))
        raise
    remove_earlier_report(reports)
    write_report(collected, outdir)

    log_end(started, "Report written", outdir)
