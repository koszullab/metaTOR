#!/usr/bin/env python3
# coding: utf-8

"""Contact map of a metaTOR run as a .metacool file.

A .metacool file is the HDF5 contact map of metacooler: a cooler whose
chromosomes are the contigs, grouped by genome. It stays readable by
``cooler.Cooler``, and ``metacooler.MetaCooler`` adds the genome layer.

metaTOR writes one bin per contig, and one genome per final bin: the
``Final_cluster`` of the contig data file. The unbinned contigs, whose
``Final_cluster`` is ``ND``, are kept with all their contacts, in metacooler's
``__unbinned__`` genome: a run with no final bin gives a map whose contigs
are all there. The columns of the bin summary (completeness,
contamination, abundances...) are stored in the genomes table, and what
metaTOR measured on each contig (GC content and Hi-C hits, plus shotgun
coverage and restriction sites when it measured them) in the contigs table.

A run starts from a .metacool map of the contigs, made before any binning
(by hicstuff pipeline -M metacool, or metacooler create pairs, from the
pairs of the Hi-C libraries): it reads its network from that map, and
writes the map of its clusters from it at the end, with its bins and
contacts, grouped by final cluster. A run step by step from a map stores
the overlapping clusters of its partition in the map first, as a
clustering round, then its validation writes the map of the clusters.

The map of the clusters of a pipeline run also stores the normalization of
its network, as a weight column of its bins table (see
store_normalization): the raw contacts stay as they are, and the
normalized contact of two bins is their count times both weights, the
way cooler stores its balancing weights.

Every map is written through metacooler, in an edit session
(metacooler.staged_edit): all the edits of a step go into one fresh copy of
the map, which metacooler validates in full once, then renames over the map,
or onto the map that the run writes instead (--output-map). The map is
whole whatever stops the run, and does not grow from run to run. A run
holds the lock of the map it writes from its checks to its end
(metacooler.lock, see metator.cli.hold_lock), which the session re-enters.

Core functions to write the file:
    - reassign_metacool, from a .metacool map
    - write_rounds, a .metacool map with the clustering rounds of a run
    - store_normalization, the normalization of the network in a map
"""

import gzip
import logging
import os
import tempfile
from contextlib import ExitStack, contextmanager

import numpy as np
import pandas as pd

import metator.io as mio
from metator.log import logger
from metator.version import __version__

# Final_cluster of a contig that is in no final bin: "ND" as written by metaTOR,
# "-" or empty when a table was edited by hand.
UNBINNED = {"ND", "-", ""}

# The genome names that metacooler reads as no genome, stripped and whatever
# their case (_UNBINNED_SENTINELS in metacooler/io/create.py; -1 from
# metacooler 0.2.0): it puts the contigs of a bin so named in __unbinned__,
# without an error. metaTOR's own bin names (metator_*, MetaTOR_MGE_*) are
# never one of them.
METACOOLER_UNBINNED = frozenset({"__unbinned__", "na", "nan", "n/a", "nd", "n/d", "-1", ""})

# The columns of contig_data_final.txt stored in the contigs table: what
# metaTOR measured on each contig. Never a bin column (Final_cluster, the bin IDs
# of the steps): the genome of a contig is its genome_id, and metacooler
# reassign carries the contigs columns into the map of another binning, where
# a bin column would be stale.
CONTIG_COLUMNS = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]

# The columns of the contigs table of a map that each normalization of the
# network reads (see metator.network.normalize_pair), besides the length of
# the contigs: the Hi-C hits that the network step counts on the contacts,
# and the shotgun coverage and restriction sites that it measures on the
# assembly (CONTIG_COLUMNS).
NORMALIZATION_COLUMNS = {
    "abundance": ("Shotgun_coverage",),
    "length": ("Hit",),
    "RS": ("Hit", "Restriction_site"),
    "empirical_hit": ("Hit",),
    "theoritical_hit": ("Shotgun_coverage", "Restriction_site"),
}

# The compressions that pandas infers from the end of a file name, in any
# case. Of these, metacooler reads only gzip, from a name that ends in .gz.
PANDAS_COMPRESSIONS = (".tar", ".tar.gz", ".tar.bz2", ".tar.xz", ".gz", ".bz2", ".zip", ".xz", ".zst")

# The first bytes of an HDF5 file, as h5py writes it (cooler, metacooler).
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"

# How many offending names an error message shows.
MAX_SHOWN = 5


class MetacoolError(ValueError):
    """The inputs cannot give a .metacool file; the message says why."""


@contextmanager
def _as_metacool_errors():
    """Raise the ValueError, OSError (FileExistsError among them) and
    KeyError of writing the map, metacooler's above all, as a MetacoolError,
    its message kept: the command logs the reason, then raises the
    error."""
    try:
        yield
    except MetacoolError:
        raise
    except KeyError as error:
        # The str() of a KeyError is the repr() of its message.
        raise MetacoolError(str(error.args[0]) if len(error.args) == 1 else str(error)) from error
    except (ValueError, OSError) as error:
        raise MetacoolError(str(error)) from error


@contextmanager
def _metacooler_logging():
    """Send the messages of metacooler to the handlers of metaTOR (terminal
    and log file) while the file is written."""
    mc_logger = logging.getLogger("metacooler")
    saved = (mc_logger.level, mc_logger.propagate)
    handlers = [handler for handler in logger.handlers if handler not in mc_logger.handlers]
    for handler in handlers:
        mc_logger.addHandler(handler)
    mc_logger.setLevel(logger.level)
    mc_logger.propagate = False
    try:
        yield
    finally:
        for handler in handlers:
            mc_logger.removeHandler(handler)
        mc_logger.level, mc_logger.propagate = saved


def _shown(names, show=repr):
    """A few names, for an error message: the first MAX_SHOWN, each written
    by show (quoted by default), then how many more there are (not a bare
    "...", which the period ending the message would make four dots)."""
    shown = ", ".join(show(name) for name in names[:MAX_SHOWN])
    more = len(names) - MAX_SHOWN
    return shown + (f", ... (+{more} more)" if more > 0 else "")


def read_contig_data(contig_data):
    """The contigs of a contig data file, read once: their bins (one per
    contig), their final bins and what metaTOR measured on them.

    Parameters
    ----------
    contig_data : str
        Path to the contig_data_final.txt file of a validation or pipeline
        run: its Name, Size and Final_cluster columns are read, and those of
        CONTIG_COLUMNS that it has.

    Returns
    -------
    pandas.DataFrame
        One bin per contig, from 0 to its size: columns chrom, start and
        end.
    dict
        Final bin of each binned contig, by contig name. The unbinned
        contigs (Final_cluster ND, - or empty) are left out.
    list of str
        The columns of CONTIG_COLUMNS that the file has and fills, in that
        order: a column that is ``-`` on every contig is left out, as
        metaTOR had nothing there (Shotgun_coverage without a depth file,
        Restriction_site without an enzyme, nor such a column of the map
        that a run started from, see
        metator.network.contig_data_from_map).
    list of str
        All the columns of the file, as pandas names them.

    Raises
    ------
    MetacoolError
        If the file is compressed otherwise than with gzip as .gz (pandas
        would read it, but not metacooler, which reads it again to store the
        contig columns), or pandas cannot read it as a table (empty, not
        UTF-8...), or
        would read it shifted or padded (a row with more or fewer fields
        than the header: the surplus leading fields of the first rows as
        their index, every value then under the name of the column to its
        left, or the fields a row lacks as empty, its contig then in no
        bin), or it has no Name, Size or Final_cluster column, or a Size that is
        not a whole number.
    """
    name = str(contig_data)
    if name.lower().endswith(PANDAS_COMPRESSIONS) and (not name.endswith(".gz") or name.lower().endswith(".tar.gz")):
        raise MetacoolError(
            f"{contig_data} is compressed in a way that metacooler, which reads it again to store the contig columns, "
            "cannot read: give it plain, or compressed with gzip as .gz."
        )
    try:
        _check_field_counts(contig_data)
        table = pd.read_csv(contig_data, sep="\t", dtype=str, keep_default_na=False)
    except (ValueError, OSError) as error:
        raise MetacoolError(f"{contig_data} cannot be read as a table: {str(error).strip().rstrip('.')}.") from error
    missing = [column for column in ["Name", "Size", "Final_cluster"] if column not in table.columns]
    if missing:
        raise MetacoolError(
            f"{contig_data} has no {' or '.join(missing)} column: give the contig_data_final.txt file of a "
            "validation or pipeline run."
        )
    try:
        sizes = table["Size"].astype("int64")
    except (ValueError, OverflowError) as error:
        raise MetacoolError(
            f"The Size column of {contig_data} must give the length of each contig in bp, as a whole number: {error}."
        ) from error
    bins = pd.DataFrame({"chrom": table["Name"], "start": 0, "end": sizes})
    final_bin = table["Final_cluster"].str.strip()
    binned = ~final_bin.isin(UNBINNED)
    measured = [column for column in CONTIG_COLUMNS if column in table and not table[column].str.strip().eq("-").all()]
    return bins, dict(zip(table["Name"][binned], final_bin[binned], strict=True)), measured, list(table.columns)


def _check_field_counts(path):
    """Raise a ValueError naming the rows of a tab-separated table whose
    number of fields is not that of its header (its first line that is not
    blank: pandas skips the empty lines and those of spaces).

    A header that ends in empty names (a trailing tab) may have more fields
    than the rows: pandas reads the fields they lack as empty, in columns
    that have no name, as metacooler reads such a table. The file is plain,
    or compressed with gzip if its name ends in .gz."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, ragged = None, []
    with opener(path, "rb") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip(b" \r\n"):
                continue
            count = line.count(b"\t") + 1
            if header is None:
                header = (len(line.rstrip(b"\r\n").rstrip(b"\t").split(b"\t")), count)
            elif not header[0] <= count <= header[1]:
                ragged.append(f"line {number} ({count})")
    if ragged:
        fields = str(header[1]) if header[0] == header[1] else f"{header[0]} to {header[1]}"
        raise ValueError(
            f"{len(ragged)} row(s) do not have the {fields} tab-separated fields of the header: {_shown(ragged, show=str)}"
        )


def open_map(metacool):
    """The metacooler.MetaCooler of a .metacool map, or of an entry of a
    .mmetacool collection: a map that cannot be read (a missing file, a file
    that is not a map, a collection without the entry...) is a
    MetacoolError naming it."""
    import metacooler

    try:
        return metacooler.MetaCooler(metacool)
    except (OSError, KeyError, ValueError) as error:
        reason = error.args[0] if isinstance(error, KeyError) and error.args else error
        raise MetacoolError(f"Cannot read the contact map {metacool}: {reason}") from error


@contextmanager
def read_errors(metacool):
    """Raise, as a MetacoolError, metacooler's refusal of a read of the map
    metacool that spans a replacement of its file (MapChangedError): a run
    that updates a map in place replaces it with the map of its clusters,
    whose bins are in another order, and so are the bin ids of its pixels.
    metacooler reads a map one chunk of pixels at a time, and checks its file
    before and after each: a read that spanned the replacement would mix the
    bins of one map with the pixels of the other."""
    import metacooler

    try:
        yield
    except metacooler.MapChangedError as error:
        raise MetacoolError(
            f"{metacool} was replaced or written while it was read, by a run that updated it: what was read of it "
            "would mix two maps. Read it again once that run has ended."
        ) from error


def check_map_lengths(metacool, contigs, lengths, of, where):
    """Refuse a map that gives contigs another length than lengths (by
    contig name; the contigs it lacks are not checked): it is not a map of
    the contigs of `of` (the assembly, or a contig data file), where (`where`)
    they have these lengths."""
    wrong = [
        f"{str(name)!r} ({int(length)} bp in the map, {lengths[str(name)]} bp in {where})"
        for name, length in zip(contigs["name"], contigs["length"], strict=True)
        if str(name) in lengths and int(length) != lengths[str(name)]
    ]
    if wrong:
        raise MetacoolError(
            f"{metacool} is not a contact map of the contigs of {of}: {len(wrong)} contig(s) have another length in "
            f"it: {_shown(wrong, show=str)}."
        )


def check_map_of_contigs(metacool, contig_data):
    """Refuse a map that is not a map of the contigs of a run: contigs that
    contig_data, a contig data file of the run (its Name and Size columns),
    has not, or that have another length in it. The map may lack some of
    them, as a map made with a contig filter does.

    Raises
    ------
    MetacoolError
        If the map holds other contigs or lengths, cannot be read, or the
        contig data has no Name or Size column.
    """
    import pandas as pd

    try:
        table = pd.read_csv(contig_data, sep="\t", usecols=["Name", "Size"], dtype={"Name": str}, keep_default_na=False)
        lengths = dict(zip(table["Name"], table["Size"].astype(int), strict=True))
    except (OSError, ValueError) as error:
        raise MetacoolError(f"Cannot read the contigs of {contig_data}: {error}") from error
    with _metacooler_logging():
        contigs = open_map(metacool).contigs
    extra = [str(name) for name in contigs["name"] if str(name) not in lengths]
    if extra:
        raise MetacoolError(
            f"{metacool} is not a contact map of the contigs of {contig_data}: {len(extra)} of its contig(s) are not "
            f"in it: {_shown(extra)}."
        )
    check_map_lengths(metacool, contigs, lengths, contig_data, contig_data)


def map_file(metacool):
    """The file of a .metacool map: the map, or the .mmetacool collection of
    an entry given as file.mmetacool::entries/NAME."""
    return str(metacool).split("::", 1)[0]


def _is_hdf5(path):
    """Whether path is a file that starts as an HDF5 file does."""
    try:
        with open(path, "rb") as handle:
            return handle.read(len(HDF5_SIGNATURE)) == HDF5_SIGNATURE
    except OSError:
        return False


def is_map(path):
    """Whether path is a contact map rather than a pairs file: an HDF5
    file, which a pairs file never is, or the URI of a map inside one
    (file::group). A path that is a file is a map only when it is an HDF5
    file, whatever its name: a pairs file may have :: in its name. A URI is
    one only when the path is no file and its file part is an HDF5 file. A
    .metacool map is a map; so are a .mmetacool collection and a .cool
    file, which open_map then refuses with metacooler's reason."""
    path = str(path)
    if os.path.isfile(path):
        return _is_hdf5(path)
    return "::" in path and _is_hdf5(map_file(path))


def is_entry(metacool):
    """Whether metacool names an entry inside a file (a .mmetacool
    collection): file::group with a group other than the root of the file
    (file::/ is the file itself, as metacooler reads it)."""
    return "::" in str(metacool) and str(metacool).split("::", 1)[1].strip().strip("/") != ""


def is_input_map(metacool, outfile):
    """Whether outfile is the file of the map metacool: the same path, or a
    symbolic link to it or from it."""
    source = map_file(metacool)
    return os.path.exists(outfile) and os.path.exists(source) and os.path.samefile(source, outfile)


def non_ascii_names(names):
    """The names that are not ASCII, each once, in their order."""
    return list(dict.fromkeys(name for name in names if not name.isascii()))


def check_ascii_names(contigs, bins):
    """Refuse contig or cluster names that are not ASCII.

    cooler reads the contig names (its chromosome names) as ASCII text, so
    metacooler refuses any other name.
    """
    for kind, names in [("contig", contigs), ("cluster", bins)]:
        bad = non_ascii_names(names)
        if bad:
            raise MetacoolError(
                f"A .metacool file needs ASCII {kind} names (cooler reads them as ASCII text), and "
                f"{len(bad)} {kind} name(s) are not: {_shown(bad)}."
            )


def check_bin_names(bins):
    """Refuse bin names that metacooler reads as no bin.

    metacooler would put the contigs of these bins in ``__unbinned__``
    without a word (ND, metaTOR's own mark of a contig in no bin, is read
    before: it is no bin name).
    """
    bad = list(dict.fromkeys(name for name in bins if name.strip().lower() in METACOOLER_UNBINNED))
    if bad:
        raise MetacoolError(
            "A .metacool file cannot have clusters named as metacooler marks a contig in no cluster (NA, NaN, N/A, N/D, "
            "-1 or __unbinned__ in any case, and nd in any case other than ND, which marks a contig in no cluster): it "
            f"would put their contigs in __unbinned__. Rename the {len(bad)} cluster(s): {_shown(bad)}."
        )


def _check_paths(outfile, force, inputs):
    """What every writer checks first, in this order: the outfile (a file:
    not a URI, file::group, which metacooler would read as a group inside
    the file and write there, truncating it; replaced only with force), then
    the input files (None for one not given)."""
    if "::" in str(outfile):
        raise MetacoolError(
            f"{outfile} is not the path of a file: a map is written to a file of its own, not to a group inside one."
        )
    if os.path.exists(outfile) and not force:
        raise MetacoolError(f"{outfile} already exists. Use --force to overwrite it.")
    missing = [path for path in inputs if path and not os.path.isfile(path)]
    if missing:
        raise MetacoolError(f"Input file(s) not found: {', '.join(map(str, missing))}.")


def _checked_inputs(contig_data, outfile, force, inputs):
    """What the writer of the map of the bins (reassign_metacool) checks
    before any work, in this order: the outfile and the input files (see
    _check_paths), then the contig data and the names of its clusters.
    Returns what read_contig_data reads. The contig names need no check:
    those of the map are ASCII, as metacooler refuses others, and the
    contigs of the contig data that the map lacks, a contig of the assembly
    with another name among them, are left out of the map of the bins."""
    _check_paths(outfile, force, inputs)
    bins, contig_bin, contig_columns, header = read_contig_data(contig_data)
    check_ascii_names((), contig_bin.values())
    check_bin_names(contig_bin.values())
    return bins, contig_bin, contig_columns, header


def reassign_metacool(
    metacool,
    contig_data,
    outfile,
    bin_summary=None,
    force=False,
    metadata=None,
    rounds=None,
    normalization=None,
    planned=None,
):
    """Write the contact map of the final bins from a .metacool map of the
    contigs, without pairs files.

    The bins and contacts are those of the input map, at its bin size (one
    bin per contig, or a fixed size), made before any binning (by
    metacooler create pairs without --contigs, for example). Each contig
    goes to its final bin, the unbinned ones and those that the contig data
    lacks to ``__unbinned__``, and the contigs are grouped by genome, then
    from the longest (metacooler reassign, layout genome-then-length). The
    genomes table holds the columns of the bin summary, and the contigs
    table the columns of the contig data that metaTOR measured (see
    read_contig_data), last: they replace those that the input map carried,
    metaTOR's own as those of the bins of an earlier run. Its other columns
    of the contigs table are kept. The QC block of the input map, when it has one, is carried over
    to the bins: metacooler keeps the counts within each contig (measured
    on the reads, for a map made from pairs files) and counts the contacts
    between contigs again. A map written by metacooler 0.2.0, with no count
    per contig, has its block computed again from the pixels (at one bin
    per contig, every contact within a contig then counts as short-range),
    which metacooler warns about for a map made from pairs files. The
    assembly name and the history of the commands are those of the input
    map.

    The map is written in an edit session of metacooler (see
    _staged_edit): metacooler reassign writes it from the input map into
    the staged copy, where the columns, the normalization and the rounds of
    the run are then stored, and metacooler validates it in full once,
    then renames it over the input map, when outfile is the input map
    itself (the same file, through a symbolic link to it or from it, whose
    target is then the map edited), or onto outfile, a new file (a symbolic
    link there is replaced, not followed). On any failure, outfile is left
    as it was, and the staged copy is removed.

    Parameters
    ----------
    metacool : str
        Path to the input .metacool file, or URI of an entry of a
        .mmetacool collection: ``file.mmetacool::entries/<name>``.
    contig_data : str
        Path to the contig_data_final.txt file: the final bin of the
        contigs (Name, Final_cluster), their size (Size, which must be their
        length in the map), and what metaTOR measured on them.
    outfile : str
        Path to the .metacool file to write.
    bin_summary : str, optional
        Path to the cluster_summary.txt file: its columns are added to the
        genomes table, one row per bin, named by its first column. It is not
        read when no contig of the map is in a final bin (the map then has
        no genome but ``__unbinned__``). Without it, the genomes table has
        none of the columns that the input map carried, only those that
        metacooler derives.
    force : bool
        Whether to replace an existing outfile.
    metadata : dict, optional
        Parameters of the run, stored with the version of metaTOR under the
        "metator" key of the file metadata.
    rounds : tuple of pandas.DataFrame, optional
        The clustering rounds of the run, as metator.clusters.read_rounds
        reads them: stored in the /clusters group of the map (see
        store_rounds), which replaces the rounds that the input map held.
        Without them, the map keeps none.
    normalization : str, optional
        The normalization of the network of the run (see
        metator.network.normalize_pair), stored in the map as a weight
        column of its bins (see store_normalization). None or "None":
        none is stored.
    planned : metacooler.MetaCooler, optional
        The input map as the run read it, opened when the run started: the
        map is refused if it changed since (another run replaced it, or a
        program wrote it), rather than written from another map than the
        one the clusters of the run come from. By default, the map as this
        function reads it.

    Raises
    ------
    MetacoolError
        If the outfile exists (without force) or is the .mmetacool
        collection of the input map, an input file is missing, the contig
        data cannot be read as a table or has no Name, Size or Final_cluster
        column or a Size that is not a whole number, a contig or bin name is
        not ASCII, a bin name is one that metacooler reads as no bin (NA,
        nan, -1...), the map cannot be read or gives a contig another length
        than the contig data, or metacooler cannot write the map (a
        ValueError, OSError or KeyError, message kept). An error of another
        type is raised as it is.
    """
    metacool, outfile = str(metacool), str(outfile)
    if is_entry(metacool) and is_input_map(metacool, outfile):
        raise MetacoolError(
            f"{outfile} is the .mmetacool collection of the input map {metacool}: the map of the clusters would replace "
            "the whole collection. Write it to another file."
        )
    bins, contig_bin, contig_columns, header = _checked_inputs(
        contig_data, outfile, force, [map_file(metacool), contig_data, bin_summary]
    )

    import metacooler

    with _metacooler_logging():
        mc = open_map(metacool)
        contigs = mc.contigs
        # Said once the map is written: a refused run drops nothing.
        later = _later_rounds(mc)
    in_place = is_input_map(metacool, outfile)
    names = [str(name) for name in contigs["name"]]
    # The binning of another version of the assembly (other contigs under
    # the same names) is no binning of this map.
    check_map_lengths(metacool, contigs, dict(zip(bins["chrom"], bins["end"], strict=True)), contig_data, "the contig data")
    in_map = set(names)
    absent = [name for name in contig_bin if name not in in_map]
    if absent:
        logger.warning(
            f"{len(absent)} clustered contig(s) of {contig_data} are not in {metacool}, and not in the map of the "
            f"clusters either: {_shown(absent)}."
        )
    # Every contig of the map, the unbinned ones as ND (no genome for
    # metacooler): never an empty table, which metacooler refuses.
    assignment = {name: contig_bin.get(name, "ND") for name in names}
    final_bins = set(assignment.values()) - {"ND"}
    logger.info(
        f"Writing {outfile} from the contact map {metacool}: {len(names)} contigs, "
        f"{sum(value != 'ND' for value in assignment.values())} of them in {len(final_bins)} clusters."
    )
    with _metacooler_logging(), _as_metacool_errors(), ExitStack() as session:
        _make_directory_of(outfile)
        # Every edit in one staged copy, which replaces the map once it is
        # valid (see _staged_edit): reassign writes it from the input map,
        # then the columns, the normalization and the rounds of the run are
        # stored in it.
        with read_errors(metacool):
            staged = session.enter_context(_staged_edit(metacool, outfile, planned=_read(mc, planned), force=force))
        scratch = session.enter_context(tempfile.TemporaryDirectory())
        if bin_summary and final_bins:
            genomes_table = bin_summary
        else:
            # Without a bin summary or without a bin, the bin names alone:
            # the genome columns that the input map carried (those of an
            # earlier run's bins) are then replaced by none, where metacooler
            # would keep them and warn that they no longer describe their
            # bins.
            genomes_table = os.path.join(scratch, "bins.txt")
            with open(genomes_table, "w") as handle:
                handle.write("".join(f"{name}\n" for name in ["bin", *sorted(final_bins)]))
        with read_errors(metacool):
            metacooler.reassign(
                metacool,
                assignment,
                # The staged map, written anew from the input map: its
                # history names the input map as the source.
                staged,
                # Grouped by genome, then from the longest, whatever the
                # order of the input map: the pixels are remapped.
                layout="genome-then-length",
                # The bin names are the first column of the table, unnamed
                # in cluster_summary.txt. It replaces the columns of the
                # genomes table that the input map carried.
                genomes_table=genomes_table,
                genomes_table_format="1",
                metadata={"metator": {"version": __version__, **(metadata or {})}},
                # Over the staged copy of the map.
                force=True,
                # On a failure, metacooler removes the partial file it
                # writes next to the staged map.
                keep_partial=False,
            )
            # The input map as the run read it (mc): one replaced or written
            # since, by a program that takes no lock of metacooler, is
            # refused (a read through mc checks its file).
            with mc.open():
                pass
        # The normalization of an earlier run, which metacooler reassign
        # carries over: it was computed from the columns that
        # _write_back_columns replaces.
        _drop_earlier_normalizations(staged)
        _write_back_columns(staged, contig_data, contig_columns, header)
        if normalization in (None, "None"):
            pass
        elif not stores_weights():
            logger.info(
                f"{outfile} stores no {normalization} normalization of the network: metacooler "
                f"{metacooler.__version__} cannot store weights (metacooler.store_weights comes with its next "
                "release)."
            )
        else:
            try:
                weights = _normalization_weights(staged, normalization)
            except MetacoolError as error:
                # Refused before anything is written: the map is written
                # without it, an extra output of the run, as the map itself.
                logger.warning(f"{outfile} stores no {normalization} normalization of the network. {error}")
            else:
                # A write that fails fails the session: the map is not
                # written, with the reason.
                _store_normalization(staged, normalization, weights)
        store_rounds(staged, rounds, name=metacool)
    # Once the session has validated the map and renamed it.
    logger.info(f"{outfile} is written: the contact map of the clusters.")
    with _metacooler_logging():
        _compact(outfile)
    _warn_dropped_rounds(later, metacool, outfile, in_place)


def _make_directory_of(outfile):
    """Make the directory of outfile, where metacooler locks and stages the
    map before it renames it there."""
    os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)


def _read(mc, planned):
    """The input map as it was read, which an edit session refuses to write
    from when it changed since (see metacooler.staged_edit): mc, opened to
    write the map, and planned, when the run gives the map as it read it at
    its start."""
    return [mc] if planned is None else [planned, mc]


def _staged_edit(metacool, outfile, planned, force=False):
    """The edit session of metacooler that a map of a run is written in
    (metacooler.staged_edit): one staged copy of the map metacool, where
    every edit of a step goes, validated in full once, when the session
    ends. When outfile is the input map itself (the same file, through a
    symbolic link to it or from it), the staged map then replaces the input
    map (its target, through a symbolic link, which keeps pointing at it),
    under the lock of the map, which the run holds, or which the session
    takes; otherwise it is renamed onto outfile, a new file, which it
    replaces with force (a symbolic link there is replaced, not followed):
    the copy is of the map alone, so an entry of a .mmetacool collection
    becomes a map of its own. The session refuses a map that changed since
    planned, the input map as the run read it, was opened."""
    import metacooler

    if is_input_map(metacool, outfile):
        return metacooler.staged_edit(metacool, planned=planned)
    return metacooler.staged_edit(metacool, output=outfile, force=force, planned=planned)


def _compact(outfile):
    """STOPGAP, until metacooler's edit sessions copy the map they install
    into a fresh file: the writers of a session write its staged map in
    place, and the session installs it with the space of what they replaced,
    which metacooler reassign then copies as it is into the next run's map.
    A session with no edit copies the map fresh. Without it, a map updated
    in place grew by 125 kB a run (the map of the synthetic test data, whose
    content is 121 kB). The map is written already: a failure is a
    warning."""
    import metacooler

    try:
        with metacooler.staged_edit(outfile):
            pass
    except (OSError, ValueError, KeyError) as error:
        logger.warning(
            f"{outfile} is written, but not copied into a fresh file ({error}): it holds the space that the edits "
            "making it left behind."
        )


def stores_weights():
    """Whether metacooler can store a normalization in a map:
    metacooler.store_weights, which its next release brings (metacooler
    0.3.0 has none)."""
    import metacooler

    return hasattr(metacooler, "store_weights")


def normalization_weight_name(normalization):
    """The weight column of the bins table of a map that stores a
    normalization of the network (see store_normalization)."""
    return f"metator_{normalization}"


def store_normalization(metacool, normalization):
    """Store a normalization of the network in a .metacool map, as a weight
    column of its bins table: metator_<normalization>.

    The weights are those of metator.network.normalization_weights, computed
    from the columns of the contigs table of the map that metaTOR stores
    (NORMALIZATION_COLUMNS: Hit, Shotgun_coverage, Restriction_site) and
    the length of the contigs, every bin of a contig taking the weight of
    its contig: the normalized contact of two bins is their count times
    both weights, and the contacts of two contigs normalized and summed are
    the edge of the network between them. metacooler show, dump and network
    --balance metator_<normalization> read them, as cooler, cooltools and
    HiGlass read balancing weights. The column records the normalization and
    the version of metaTOR, and what its weights depend on (the contacts,
    for a normalization by the hits, and the columns it reads): metacooler
    drops it when an edit changes one of them, as annotate does. A column of
    that name is replaced.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``. It is written in
        place.
    normalization : str
        Mode of normalization (see metator.network.normalize_pair), but
        "None".

    Raises
    ------
    MetacoolError
        If the normalization is unknown, the contigs table of the map lacks
        a column it reads, metacooler cannot store weights (metacooler
        0.3.0), or it cannot store the column (a ValueError, OSError or
        KeyError, message kept).
    """
    weights = _normalization_weights(metacool, normalization)
    with _metacooler_logging(), _as_metacool_errors():
        _store_normalization(metacool, normalization, weights)


def _normalization_weights(metacool, normalization):
    """What store_normalization checks and computes before it writes
    anything: the map, opened (a metacooler.MetaCooler), the weights of its
    contigs, and the columns of its contigs table that they were computed
    from. Raises MetacoolError when the normalization cannot be stored (see
    store_normalization); a run then writes its map without it."""
    import metacooler

    import metator.network as mtn

    if not stores_weights():
        raise MetacoolError(
            f"metacooler {metacooler.__version__} cannot store weights: metacooler.store_weights comes with its next release."
        )
    if normalization not in NORMALIZATION_COLUMNS:
        raise MetacoolError(
            f"No weights for the normalization {normalization!r}: it is none of {_shown(list(NORMALIZATION_COLUMNS))}."
        )
    columns = NORMALIZATION_COLUMNS[normalization]
    with _metacooler_logging():
        mc = open_map(metacool)
        contigs = mc.contigs
    missing = [column for column in columns if column not in contigs.columns]
    if missing:
        raise MetacoolError(
            f"The {normalization} normalization of the network reads the {' and '.join(missing)} column(s) of the "
            "contigs table, which the map lacks."
        )
    values = {column: pd.to_numeric(contigs[column], errors="coerce").to_numpy(dtype=float) for column in columns}
    weights = mtn.normalization_weights(
        contigs["length"].to_numpy(dtype=float),
        values.get("Hit"),
        values.get("Shotgun_coverage"),
        values.get("Restriction_site"),
        normalization,
    )
    return mc, weights, columns


def _store_normalization(metacool, normalization, weights):
    """Store the weights that _normalization_weights computed (weights: the
    map, the weights of its contigs and the columns they were computed
    from) with metacooler.store_weights, every bin of a contig taking the
    weight of its contig. Its errors are raised as they are: in an edit
    session of metacooler, a write that fails fails the session, and the
    map is not written."""
    import metacooler

    mc, contig_weights, columns = weights
    with mc.open() as group:
        offsets = group["indexes/chrom_offset"][()]
    metacooler.store_weights(
        metacool,
        normalization_weight_name(normalization),
        np.repeat(contig_weights, np.diff(offsets)),
        method="metator",
        # The hits are counted on the contacts that the pixels hold.
        depends_on=(["pixels"] if "Hit" in columns else []) + [f"contigs/{column}" for column in columns],
        parameters={"mode": normalization, "metator_version": __version__},
        force=True,
        # Refused if the map changed since the weights were computed from it.
        source=mc,
    )
    masked = int(np.count_nonzero(np.isnan(contig_weights)))
    logger.info(
        f"Stored the {normalization} normalization of the network in {metacool}, as the weights "
        f"{normalization_weight_name(normalization)} of its bins" + (f" ({masked} contig(s) without one)." if masked else ".")
    )


def _drop_earlier_normalizations(path):
    """Remove from a map the normalizations that an earlier run of metaTOR
    stored in it (see store_normalization): the weight columns of its bins
    whose normalization attribute (see metacooler's format) is metator.
    They describe the network of that run, from the columns that this run
    replaces. Left in the map, metacooler annotate would drop them all the
    same, with a warning to compute them again, which this run then does."""
    with open_map(path).open("r+") as group:
        bins = group["bins"]
        names = [name for name in bins if _attribute(bins[name], "normalization") == "metator"]
        for name in names:
            del bins[name]
    if names:
        logger.info(f"Removed the normalization(s) of the network of an earlier run from the map: {_shown(names)}.")


def _attribute(node, name):
    """An attribute of an HDF5 object, as text when it is stored as bytes
    (None when the object has none of that name)."""
    value = node.attrs.get(name)
    return value.decode() if isinstance(value, bytes) else value


def _write_back_columns(path, contig_data, contig_columns, header):
    """Give a map the metadata columns of this run, as metacooler annotate
    edits them: the map written back from an input map.

    The columns of the contig data that metaTOR measured (contig_columns)
    are stored last, in their order: the columns of a contig data file
    (header: all of them) that the map carried are dropped first, then
    those of this run are added. metacooler reassign stores the columns
    that it carries from its input map before genome_id or before the QC
    columns of the contigs: replaced where they were, they would come in
    another order in each map, as the input map had them. The dropped
    columns are also those that this run does not store: the bins of each step, the contig ids, Size..., which a map made
    with the contig data of an earlier run as its contigs table holds
    (metacooler create pairs --contigs-format Name,Final_cluster), and one of
    CONTIG_COLUMNS that this run had no values for (Shotgun_coverage
    without a depth file, nor a column of the input map that gives every
    contig its coverage...).
    """
    import metacooler

    stored = metacooler.MetaCooler(path)
    # The columns of a contig data file but Name (the contigs' own name
    # column) and those of the map that metacooler derives, its cluster
    # columns under their names of metaTOR 1.x too, which a map made with
    # the contig data of a 1.x run carries (see
    # metator.io.CONTIG_DATA_1X_COLUMNS): they would contradict the clusters
    # of this run.
    of_1x = {old for old, new in mio.CONTIG_DATA_1X_COLUMNS.items() if new in header}
    of_metator = (set(header) | of_1x | set(CONTIG_COLUMNS)) - {"Name"}
    carried = [f"contigs.{column}" for column in stored.contigs.columns if column in of_metator - {"name", "length", "genome_id"}]
    if carried:
        metacooler.annotate(path, drop=carried, inplace=True)
    if not contig_columns:
        return
    # Columns that the map no longer has: annotate adds them, with the values
    # of this run on every contig of the map (none on a contig that the
    # contig data lacks).
    metacooler.annotate(
        path,
        contigs_table=contig_data,
        contigs_table_format="Name",
        # Named one by one: the other columns of the table are bins, or the
        # unnamed pandas index.
        columns=[f"contigs.{column}" for column in contig_columns],
        # The contig data of a run lists every contig of the assembly, and a
        # map may lack some (a map made with a contig filter): no warning
        # about the rows that name them.
        unmatched="ignore",
        inplace=True,
    )


# The steps of the binning, which make the clustering rounds that a run
# replaces (see store_rounds): a round that another step made (that of
# metator mge --map), or another tool, is dropped with them.
BINNING_STEPS = ("partition", "validation")


def _later_rounds(mc):
    """The entries of the clustering rounds of the map mc that the binning
    did not make: that of metator mge --map, or of another tool. A run that
    replaces the clusters of the map stores its own rounds only, so its map
    does not keep them (see _warn_dropped_rounds)."""
    return [entry for entry in mc.cluster_rounds if entry.get("step") not in BINNING_STEPS]


def _made_by(entry):
    """What made a clustering round, from its entry: its tool and step
    (metaTOR mge)."""
    return " ".join(str(entry[key]) for key in ("tool", "step") if entry.get(key)) or "a tool it does not name"


def _warn_dropped_rounds(entries, metacool, outfile, in_place):
    """Warn, once outfile is written from the map metacool (in place, or to
    another file), that it does not keep the clustering rounds entries of
    the map, which the run did not make (see _later_rounds): each named
    with what made it, to be run again on the new clusters."""
    if not entries:
        return
    rounds = ", ".join(f"{entry['name']} (by {_made_by(entry)})" for entry in entries)
    where = f"{metacool} no longer holds" if in_place else f"{outfile} does not keep, of {metacool},"
    logger.warning(
        f"{where} the clustering round(s) {rounds}, which this run does not make: a run that replaces the clusters of "
        "a map stores its own rounds only. Run what made them again, on the new clusters, to store them."
    )


def store_rounds(path, rounds, names=None, name=None):
    """Store the clustering rounds of a run in the /clusters group of a map,
    in place of those it held: the overlapping, recursive and final clusters
    of metator.clusters, with their contigs, their parents and their columns
    (those of miComplete, and of the bin summary for the final clusters).

    Each round is stored with the tool and the version that made it, and
    the step that made it: partition for the round overlapping, validation
    for the others. A text column holds an empty text where a round has no
    value, as HDF5 stores no missing text. The contigs of a cluster that
    the map lacks (a map made with a contig filter) are left out, with a
    warning.

    Parameters
    ----------
    path : str
        Path to the .metacool file, written in place.
    rounds : tuple of pandas.DataFrame, or None
        The clusters and their contigs (see metator.clusters.read_rounds).
        None removes the rounds that the map held.
    names : list of str, optional
        The names of the rounds, in the order they were made: by default,
        those of the clusters, in their order. A round named here may have
        no cluster: the overlapping clusters of a partition when none of
        them reaches its --size.
    name : str, optional
        The map that path is staged for, which the warnings name (path by
        default).
    """
    import metacooler

    if rounds is None:
        metacooler.drop_clusters(path)
        return
    clusters, members = rounds
    order = list(dict.fromkeys(clusters["round"])) if names is None else list(names)
    steps = {name: "partition" if name == "overlapping" else "validation" for name in order}
    _store_clusters(
        path,
        [{"name": name, "tool": "metaTOR", "version": __version__, "step": steps[name]} for name in order],
        clusters,
        members,
        name=name,
    )


def _store_clusters(path, entries, clusters, members, name=None):
    """Write the /clusters group of the map at path, in place of the one it
    held: the rounds entries (each a dict with its name and what is stored
    with it), and their clusters and members (see store_rounds). A text
    column holds an empty text where a round has no value, and the contigs
    of a cluster that the map lacks are left out, with a warning, which
    names the map (name: the map that path is staged for; path by
    default)."""
    import metacooler

    contigs = set(metacooler.MetaCooler(path).contigs["name"].astype(str))
    absent = members["contig"].astype(str).map(lambda name: name not in contigs)
    if absent.any():
        logger.warning(
            f"{int(absent.sum())} contig(s) of the clusters of the run are not in {name or path}, and not in its clusters "
            f"either: {_shown(sorted(set(members.loc[absent, 'contig'])))}."
        )
        members = members[~absent]
    stored = clusters.copy()
    for column in stored.columns:
        if column not in ("parent_round", "parent") and not pd.api.types.is_numeric_dtype(stored[column]):
            stored[column] = stored[column].astype(object).where(stored[column].notna(), "").astype(str)
    metacooler.write_clusters(path, entries, stored, members)


def write_rounds(metacool, outfile, rounds, names=None, force=False, planned=None):
    """Write the map metacool at outfile with the clustering rounds of a run
    in its /clusters group (see store_rounds), in place of those it held,
    and everything else as it was: its bins, contacts, genomes and contigs
    tables, QC block and metadata. A partition run from a map stores its
    overlapping clusters so, before its validation writes the map of the
    final clusters (see reassign_metacool).

    The rounds are stored in an edit session of metacooler, as
    reassign_metacool writes a map (see _staged_edit): in a fresh copy of
    the map, which metacooler validates in full, then renames over the
    input map (when outfile is the input map itself: the same file, through
    a symbolic link to it or from it), which keeps its mode, or onto
    outfile, a new file. On any failure, outfile is left as it was, and the
    staged copy is removed. An entry of a .mmetacool collection is written
    as a map of its own.

    Parameters
    ----------
    metacool : str
        Path to the input .metacool file, or URI of an entry of a
        .mmetacool collection: ``file.mmetacool::entries/<name>``.
    outfile : str
        Path to the .metacool file to write.
    rounds : tuple of pandas.DataFrame
        The clusters and their contigs (see store_rounds).
    names : list of str, optional
        The names of the rounds, in their order (see store_rounds).
    force : bool
        Whether to replace an existing outfile: the input map itself, which
        a run updates in place, or another file.
    planned : metacooler.MetaCooler, optional
        The input map as the run read it (see reassign_metacool).

    Raises
    ------
    MetacoolError
        If the outfile is no file's path (a URI), exists (without force),
        or is the .mmetacool collection of the input map, the input map is
        missing, cannot be read or changed since planned was opened, or the
        rounds cannot be stored (a ValueError, OSError or KeyError of
        metacooler, message kept).
    """
    # The rounds of the map that the run does not make: said once the map
    # is written.
    later = []
    in_place = is_input_map(metacool, outfile)

    def store(staged, mc):
        later.extend(_later_rounds(mc))
        store_rounds(staged, rounds, names, name=metacool)

    _write_with_rounds(
        metacool,
        outfile,
        force,
        store,
        f"with the {len(rounds[0])} clusters of the run",
        "with the clustering rounds of the run",
        planned,
    )
    _warn_dropped_rounds(later, metacool, outfile, in_place)


def write_round(metacool, outfile, rounds, entry, force=False, planned=None):
    """Write the map metacool at outfile with one clustering round of its
    /clusters group replaced, and everything else as it was: the other
    rounds, its bins, contacts, genomes and contigs tables, QC block and
    metadata. metator mge stores its MGE clusters so, as the round mge,
    next to the rounds of the binning run, whose final clusters the genomes
    and contigs tables keep.

    The round replaces the round of its name that the map held, if any, and
    comes after the others, as the last made. The other rounds are kept as
    the map held them (see metacooler's MetaCooler.clusters): the entry of
    each, its clusters, their parents, contigs and the values of their
    columns, but for those that metacooler computes (n_contigs, length and
    genome), which it computes again. In the one table of the clusters of every round, a column that
    only one of them has is missing in the other: an empty text, or no
    number. So a column of whole numbers, or of booleans, that the round
    lacks is stored as decimals, or as text (True, False), as any column
    with a missing value is. The map is
    written as write_rounds writes it: in an edit session of metacooler,
    whose copy of the map stores the rounds, validated in full, then
    renamed over the input map, or onto outfile (see _staged_edit).

    Parameters
    ----------
    metacool : str
        Path to the input .metacool file, or URI of an entry of a
        .mmetacool collection: ``file.mmetacool::entries/<name>``.
    outfile : str
        Path to the .metacool file to write.
    rounds : tuple of pandas.DataFrame
        The clusters and the contigs of the round, as
        metator.clusters.Rounds.tables gives them: every cluster of the
        round named by entry.
    entry : dict
        The round as metacooler stores it: its name, and what is stored with
        it (tool, version, step, parameters...).
    force : bool
        Whether to replace an existing outfile: the input map itself, which
        a run updates in place, or another file.
    planned : metacooler.MetaCooler, optional
        The input map as the run read it (see reassign_metacool).

    Raises
    ------
    MetacoolError
        If the outfile is no file's path (a URI), exists (without force),
        or is the .mmetacool collection of the input map, the input map is
        missing, cannot be read or changed since planned was opened, or the
        rounds cannot be stored (a ValueError, OSError or KeyError of
        metacooler, message kept).
    """
    clusters, members = rounds
    name = entry["name"]

    def store(staged, mc):
        held = mc.clusters
        # The columns that metacooler computes, which it refuses as given
        # (the /clusters group of metacooler's format).
        computed = ["n_contigs", "length", "genome"]
        other = held[(held["round"] != name).to_numpy()].drop(columns=computed)
        held_members = mc.cluster_members()
        other_members = held_members[(held_members["round"] != name).to_numpy()]
        entries = [stored for stored in mc.cluster_rounds if stored["name"] != name]
        _store_clusters(staged, [*entries, entry], _stacked(other, clusters), _stacked(other_members, members), name=metacool)

    _write_with_rounds(
        metacool,
        outfile,
        force,
        store,
        f"with the {len(clusters)} clusters of the round {name} of its clustering rounds",
        f"with the round {name} of its clustering rounds",
        planned,
    )


def _stacked(first, second):
    """The rows of two tables, one after the other, with the columns of
    both: those of first, then the others of second."""
    return second.reset_index(drop=True) if first.empty else pd.concat([first, second], ignore_index=True)


def _write_with_rounds(metacool, outfile, force, store, writing, written, planned=None):
    """What write_rounds and write_round share: the checks of outfile, then
    an edit session on the map metacool (see _staged_edit), whose copy
    store(staged, mc) stores the rounds in (mc: the input map, a
    metacooler.MetaCooler), validated in full and renamed over the input
    map, or onto outfile. writing and written say what the map holds, in
    the messages of its writing; planned is the input map as the run read
    it (see reassign_metacool)."""
    metacool, outfile = str(metacool), str(outfile)
    if is_entry(metacool) and is_input_map(metacool, outfile):
        raise MetacoolError(
            f"{outfile} is the .mmetacool collection of the input map {metacool}: the map would replace the whole "
            "collection. Write it to another file."
        )
    _check_paths(outfile, force, [map_file(metacool)])
    with _metacooler_logging():
        mc = open_map(metacool)
        logger.info(f"Writing {outfile} from the contact map {metacool}, {writing}.")
        with _as_metacool_errors(), ExitStack() as session:
            _make_directory_of(outfile)
            with read_errors(metacool):
                staged = session.enter_context(_staged_edit(metacool, outfile, planned=_read(mc, planned), force=force))
            # The rounds of the copy are replaced by the new ones.
            store(staged, mc)
        logger.info(f"{outfile} is written: the contact map, {written}.")
        _compact(outfile)
