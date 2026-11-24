#!/usr/bin/env python3
# coding: utf-8

"""Contact map of a metaTOR run as a .metacool file.

A .metacool file is metacooler's HDF5 contact map: a cooler whose
chromosomes are the contigs, grouped by genome, readable by
``cooler.Cooler``. The map is the interface between the steps of a run:
``network`` stores the contig columns and the weights of its normalization
(the raw contacts stay as they are; a normalized contact is the count times
both weights, as cooler stores balancing weights), ``partition`` the rounds
core and overlapping (see read_partition), and ``validation`` the map of the
clusters: one genome per final cluster, ``__unbinned__`` for the others, the
summary in the genomes table and every round. ``pipeline`` stores what the
three steps store, in one map (see write_map), then what ``markers`` stores:
the marker genes of each contig, in its contigs table (see Markers).

Every map is written in one metacooler edit session, validated in full
once, then renamed onto the map of the output directory of the run (or over
the input map, with ``--inplace``): the map is whole whatever stops the run
and does not grow from run to run. A run holds the map's lock from its
checks to its end (see metator.cli.hold_lock).
"""

import logging
import os
import tempfile
from collections import Counter
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

import metator.clusters as mtcl
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

# The genome of metacooler that holds the contigs of a map that are in no
# genome.
UNBINNED_GENOME = "__unbinned__"

# The column of the genomes table of a map that metator mge --reassign writes
# (see write_map): True for the genomes that it made of the clusters of its
# round mge (its MGE genomes), False for the others. As text, as metacooler
# stores a column of a genomes table that holds no numbers (as it stores the
# column recursive of the summary of a validation).
MGE_GENOME_COLUMN = "mge"

# The columns of the network of the contigs that a map stores in its contigs
# table (see Network): what metaTOR measures on each contig, the GC content,
# the Hi-C hits counted on its contacts, and the shotgun coverage and the
# restriction sites, measured on the assembly or read from the map. Never a
# cluster column (Final_cluster, the cluster IDs of the steps): the final
# cluster of a contig is its genome, and its clusters are the clustering
# rounds of the map.
CONTIG_COLUMNS = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]

# The cluster columns of the contig data tables of metaTOR 2.0
# (contig_data_partition.txt and contig_data_final.txt), which metaTOR 1.x
# named after bins (see metator.io.CONTIG_DATA_1X_COLUMNS).
CLUSTER_COLUMNS = [
    "Core_cluster_ID",
    "Core_cluster_contigs",
    "Core_cluster_size",
    "Overlapping_cluster_ID",
    "Overlapping_cluster_contigs",
    "Overlapping_cluster_size",
    "Recursive_cluster_ID",
    "Recursive_cluster_contigs",
    "Recursive_cluster_size",
    "Final_cluster",
]

# metaTOR's own columns in the contigs table of a map, which a run that
# stores its network drops (see _write_back_columns): those of an earlier
# network, and the columns of a contig data table of metaTOR, 2.0 or 1.x, its
# unnamed index column too (Unnamed: 0, as pandas and metacooler read
# contig_data_final.txt), which a map made with that table as its contigs
# table carries (metacooler create pairs --contigs-format Name,Final_cluster
# keeps all its columns): they would contradict the network and the
# clusters of the run. Name is the name of the contigs, the map's own.
METATOR_TABLE_COLUMNS = frozenset({"Unnamed: 0", "ID", "Size", *CONTIG_COLUMNS, *CLUSTER_COLUMNS, *mio.CONTIG_DATA_1X_COLUMNS})

# The columns of the contigs table of a map that each normalization of the
# network reads (see metator.network.normalization_weights), besides the
# length of the contigs: the Hi-C hits that the network step counts on the
# contacts, and the shotgun coverage and restriction sites that it measures
# on the assembly (CONTIG_COLUMNS).
NORMALIZATION_COLUMNS = {
    "abundance": ("Shotgun_coverage",),
    "length": ("Hit",),
    "RS": ("Hit", "Restriction_site"),
    "empirical_hit": ("Hit",),
    "theoritical_hit": ("Shotgun_coverage", "Restriction_site"),
}

# The first bytes of an HDF5 file, as h5py writes it (cooler, metacooler).
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"

# How many offending names an error message shows.
MAX_SHOWN = 5


class MetacoolError(ValueError):
    """The inputs cannot give a .metacool file; the message says why."""


class MapReplacedError(MetacoolError):
    """A map was replaced or written while it was read (see read_errors): what was
    read of it would mix two maps."""


@contextmanager
def _as_metacool_errors():
    """Raise the ValueError, OSError (FileExistsError among them) and KeyError of
    writing the map, metacooler's above all, as a MetacoolError, its message
    kept: the command logs the reason, then raises the error."""
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
def _metacooler_logging(level=None):
    """Send the messages of metacooler to the handlers of metaTOR (terminal and
    log file) while the file is written: its warnings and errors; its
    accounts of what it reads and writes (INFO: the staged copy of the map,
    the contigs it reassigns...), which are details of a run of metaTOR, only
    when metaTOR's logger shows its own (DEBUG). level, when higher, raises
    the threshold (logging.WARNING: never its accounts)."""
    mc_logger = logging.getLogger("metacooler")
    saved = (mc_logger.level, mc_logger.propagate)
    handlers = [handler for handler in logger.handlers if handler not in mc_logger.handlers]
    for handler in handlers:
        mc_logger.addHandler(handler)
    threshold = logger.level if logger.level <= logging.DEBUG else logging.WARNING
    mc_logger.setLevel(threshold if level is None else max(level, threshold))
    mc_logger.propagate = False
    try:
        yield
    finally:
        for handler in handlers:
            mc_logger.removeHandler(handler)
        mc_logger.level, mc_logger.propagate = saved


def _shown(names, show=repr):
    """A few names, for an error message: the first MAX_SHOWN, each written by
    show (quoted by default), then how many more there are (not a bare "...",
    which the period ending the message would make four dots)."""
    shown = ", ".join(show(name) for name in names[:MAX_SHOWN])
    more = len(names) - MAX_SHOWN
    return shown + (f", ... (+{more} more)" if more > 0 else "")


def open_map(metacool):
    """The metacooler.MetaCooler of a .metacool map, or of an entry of a
    .mmetacool collection: a map that cannot be read (a missing file, a file
    that is not a map, a collection without the entry...) is a MetacoolError
    naming it."""
    import metacooler

    try:
        return metacooler.MetaCooler(metacool)
    except (OSError, KeyError, ValueError) as error:
        reason = error.args[0] if isinstance(error, KeyError) and error.args else error
        raise MetacoolError(f"Cannot read the contact map {metacool}: {reason}") from error


#: What to do about a map replaced while it was read (see read_errors).
READ_AGAIN = "Read it again once that run has ended."


@contextmanager
def read_errors(metacool, remedy=READ_AGAIN):
    """Raise, as a MapReplacedError (a MetacoolError), metacooler's refusal of a
    read of the map metacool that spans a replacement of its file
    (MapChangedError): a run that updates a map in place replaces it with the
    map of its clusters, whose bins are in another order, and so are the bin
    ids of its pixels. metacooler reads a map one chunk of pixels at a time,
    and checks its file before and after each: a read that spanned the
    replacement would mix the bins of one map with the pixels of the other. The
    message ends with remedy: what the reader does about it."""
    import metacooler

    try:
        yield
    except metacooler.MapChangedError as error:
        raise MapReplacedError(
            f"{metacool} was replaced or written while it was read, by a run that updated it: what was read of it "
            f"would mix two maps. {remedy}"
        ) from error


def check_map_lengths(metacool, contigs, lengths, of, where):
    """Refuse a map that gives contigs another length than lengths (by contig
    name; the contigs it lacks are not checked): it is not a map of the contigs
    of `of` (the assembly, or a contig data file), where (`where`) they have
    these lengths."""
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


def map_file(metacool):
    """The file of a .metacool map: the map, or the .mmetacool collection of an
    entry given as file.mmetacool::entries/NAME."""
    return str(metacool).split("::", 1)[0]


def _is_hdf5(path):
    """Whether path is a file that starts as an HDF5 file does."""
    try:
        with open(path, "rb") as handle:
            return handle.read(len(HDF5_SIGNATURE)) == HDF5_SIGNATURE
    except OSError:
        return False


def is_map(path):
    """Whether path is a contact map rather than a pairs file: an HDF5 file, which
    a pairs file never is, or the URI of a map inside one (file::group). A path
    that is a file is a map only when it is an HDF5 file, whatever its name: a
    pairs file may have :: in its name.

    A URI is one when its file part is an HDF5 file; so are a .mmetacool
    collection and a .cool file, which open_map then refuses with metacooler's
    reason."""
    path = str(path)
    if os.path.isfile(path):
        return _is_hdf5(path)
    return "::" in path and _is_hdf5(map_file(path))


def is_entry(metacool):
    """Whether metacool names an entry inside a file (a .mmetacool collection):
    file::group with a group other than the root of the file (file::/ is the
    file itself, as metacooler reads it)."""
    return "::" in str(metacool) and str(metacool).split("::", 1)[1].strip().strip("/") != ""


def _named(path):
    """path, absolute, the symbolic links of its directory followed but not its
    own name's: the file that a write to path replaces (a symbolic link there
    is replaced, not written through)."""
    path = os.path.abspath(path)
    return os.path.join(os.path.realpath(os.path.dirname(path)), os.path.basename(path))


def is_input_map(metacool, outfile):
    """Whether writing the map at outfile updates the map metacool in place:
    outfile, its own name not followed (see _named), names the file of the map,
    its symbolic links followed, as metacooler edits the file that a link
    points to. A symbolic link at outfile is a file of its own, which a write
    replaces, whatever it points to; and a map given as a symbolic link is
    updated in place through its target, not through outfile = the link."""
    source = map_file(metacool)
    return os.path.exists(source) and _named(outfile) == os.path.realpath(source)


def is_hard_link_of(metacool, outfile):
    """Whether outfile is the file of the map metacool under another name (a hard
    link: the same file, which is not the map's own name, nor a symbolic link).
    Written as a map of its own, it would replace that name only; metacooler
    refuses an output that is the same file as its input."""
    source = map_file(metacool)
    return (
        os.path.lexists(outfile)
        and not os.path.islink(outfile)
        and os.path.exists(source)
        and os.path.samefile(source, outfile)
        and not is_input_map(metacool, outfile)
    )


def _refuse_a_hard_link_of_the_map(metacool, outfile):
    """Refuse an outfile that is a hard link of the input map (see
    is_hard_link_of)."""
    if is_hard_link_of(metacool, outfile):
        raise MetacoolError(
            f"{outfile} is a hard link of the input map {metacool}: the same file under another name. Write the map to "
            f"another file, or to {os.path.realpath(map_file(metacool))} to update the map in place."
        )


def non_ascii_names(names):
    """The names that are not ASCII, each once, in their order."""
    return list(dict.fromkeys(name for name in names if not name.isascii()))


def check_ascii_names(contigs, bins):
    """Refuse contig or cluster names that are not ASCII."""
    for kind, names in [("contig", contigs), ("cluster", bins)]:
        bad = non_ascii_names(names)
        if bad:
            raise MetacoolError(
                f"A .metacool file needs ASCII {kind} names (cooler reads them as ASCII text), and "
                f"{len(bad)} {kind} name(s) are not: {_shown(bad)}."
            )


def check_bin_names(bins):
    """Refuse bin names that metacooler reads as no bin.

    metacooler would put the contigs of such a bin in __unbinned__ silently."""
    bad = list(dict.fromkeys(name for name in bins if name.strip().lower() in METACOOLER_UNBINNED))
    if bad:
        raise MetacoolError(
            "A .metacool file cannot have clusters named as metacooler marks a contig in no cluster (NA, NaN, N/A, N/D, "
            "-1 or __unbinned__ in any case, and nd in any case other than ND, which marks a contig in no cluster): it "
            f"would put their contigs in __unbinned__. Rename the {len(bad)} cluster(s): {_shown(bad)}."
        )


def _check_paths(outfile, force, inputs):
    """What every writer checks first, in this order: the outfile (a file: not a
    URI, file::group, which metacooler would read as a group inside the file
    and write there, truncating it; replaced only with force), then the input
    files (None for one not given)."""
    if "::" in str(outfile):
        raise MetacoolError(
            f"{outfile} is not the path of a file: a map is written to a file of its own, not to a group inside one."
        )
    if os.path.lexists(outfile) and not force:
        raise MetacoolError(f"{outfile} already exists. Use --force to overwrite it.")
    missing = [path for path in inputs if path and not os.path.isfile(path)]
    if missing:
        raise MetacoolError(f"Input file(s) not found: {', '.join(map(str, missing))}.")


def make_map_directory(outfile, named=None):
    """Make the directory of the map at outfile, where metacooler locks and stages
    the map before it renames it there. A directory that cannot be made is a
    MetacoolError, which names outfile as named says (outfile itself by
    default)."""
    try:
        os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
    except OSError as error:
        raise MetacoolError(
            f"The directory of {named or outfile} cannot be made ({error.strerror}): write the map elsewhere."
        ) from error


def _enter_session(session, metacool, outfile, planned, force):
    """Enter, in session (an ExitStack), the edit session of metacooler that
    writes the map of the run at outfile (see _staged_edit), and return its
    staged map."""
    if any(_changed(mc) for mc in planned):
        raise MetacoolError(
            f"{metacool} was replaced or written since it was read for this run (by another run, or a program that "
            "takes no lock of metacooler): the clusters of the run come from the map as it was. Run it again on the "
            "map as it is now."
        )
    return session.enter_context(_staged_edit(metacool, outfile, planned, force))


def _changed(mc):
    """Whether the file of the map mc (a metacooler.MetaCooler) was replaced or
    written since mc was opened, as metacooler checks it on every read."""
    import metacooler

    try:
        with mc.open():
            pass
    except metacooler.MapChangedError:
        return True
    return False


def _read(mc, planned):
    """The input map as it was read, which an edit session refuses to write from
    when it changed since (see metacooler.staged_edit): mc, opened to write the
    map, and planned, when the run gives the map as it read it at its start."""
    return [mc] if planned is None else [planned, mc]


def _staged_edit(metacool, outfile, planned, force=False):
    """The edit session of metacooler that a map of a run is written in
    (metacooler.staged_edit): one staged copy of the map metacool, where every
    edit of a step goes, validated in full once, when the session ends.

    In place when outfile is the input map's file (its links followed);
    otherwise the map alone is renamed onto outfile, so an entry of a
    collection becomes a map of its own."""
    import metacooler

    if is_input_map(metacool, outfile):
        return metacooler.staged_edit(metacool, planned=planned)
    return metacooler.staged_edit(metacool, output=outfile, force=force, planned=planned)


def normalization_weight_name(normalization):
    """The weight column of the bins table of a map that stores a normalization of
    the network (see store_normalization)."""
    return f"metator_{normalization}"


def store_normalization(metacool, normalization):
    """Store a normalization of the network in a .metacool map, as a weight column
    of its bins table: ``metator_<normalization>``.

    Every bin of a contig takes its contig's weight, so the normalized contact
    of two contigs is their count times both weights, as cooler stores
    balancing weights. depends_on lets metacooler drop the column when an edit
    changes what it was computed from.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    normalization : str
        Mode of normalization (see metator.network.normalization_weights), but
        "None".

    Raises
    ------
    MetacoolError
        If the normalization is unknown, the contigs table of the map lacks a
        column it reads, the map is replaced or written as it is read, or
        metacooler cannot store the column (a ValueError, OSError or KeyError,
        message kept)."""
    with _metacooler_logging(), _as_metacool_errors(), read_errors(metacool):
        _store_normalization(metacool, normalization, _normalization_weights(metacool, normalization))


def _normalization_weights(metacool, normalization):
    """What store_normalization checks and computes before it writes anything: the
    map, opened (a metacooler.MetaCooler), the weights of its contigs, and the
    columns of its contigs table that they were computed from."""
    import metator.network as mtn

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
    """Store the weights that _normalization_weights computed (weights: the map,
    the weights of its contigs and the columns they were computed from) with
    metacooler.store_weights, every bin of a contig taking the weight of its
    contig."""
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
    logger.debug(
        f"Stored the {normalization} normalization of the network in {metacool}, as the weights "
        f"{normalization_weight_name(normalization)} of its bins" + (f" ({masked} contig(s) without one)." if masked else ".")
    )


def _drop_earlier_normalizations(path, normalization=None):
    """Remove from a map the normalizations that metaTOR stored in it (see
    store_normalization): the weight columns of its bins whose normalization
    attribute, which metacooler's format documents for a weight column, is
    metator.

    They were computed from the columns the run replaces."""
    mc = open_map(path)
    with mc.open() as group:
        bins = group["bins"]
        names = [name for name in bins if _attribute(bins[name], "normalization") == "metator"]
    if not names:
        return
    with mc.open("r+") as group:
        for name in names:
            del group["bins"][name]
    kept = normalization_weight_name(normalization) if normalization not in (None, "None") else None
    removed = [name for name in names if name != kept]
    if removed:
        logger.info(f"Removed the normalization(s) of the network of an earlier run from the map: {_shown(removed)}.")


def _attribute(node, name):
    """An attribute of an HDF5 object, as text when it is stored as bytes (None
    when the object has none of that name)."""
    value = node.attrs.get(name)
    return value.decode() if isinstance(value, bytes) else value


def _write_back_columns(path, contig_data, contig_columns, metadata=None):
    """Give the map at path, the staged map of an edit session, the columns of the
    network of a run, as metacooler annotate edits them (see _store_network):
    those of the table contig_data that the run measured (contig_columns), in
    their order, after every other column of the contigs table (metacooler's QC
    columns included), whatever the map carried. metadata: keys of the metadata
    of the map that metacooler annotate sets as it adds them.

    metaTOR's own columns that the map carried (METATOR_TABLE_COLUMNS) are
    dropped first: an earlier network's, or a contig data table's, which would
    contradict this run."""
    import metacooler

    stored = metacooler.MetaCooler(path)
    carried = [f"contigs.{column}" for column in stored.contigs.columns if column in METATOR_TABLE_COLUMNS]
    if carried:
        metacooler.annotate(path, drop=carried, inplace=True)
    if not contig_columns:
        if metadata:
            _store_metadata(path, metadata)
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
        metadata=metadata,
        inplace=True,
    )


@dataclass
class Network:
    """The network of the contigs of a run, as a map stores it (see write_map).

    Attributes
    ----------
    table : pandas.DataFrame
        The contig data of the contigs of the map, in its order (see
        metator.network.contig_table): its columns that the run measured
        (GC_content and Hit, and Shotgun_coverage and Restriction_site when it
        has them) go to the contigs table of the map.
    normalization : str
        The normalization of the network, stored as the weights
        ``metator_<normalization>`` of the bins of the map (see
        store_normalization): none for "None", the raw counts.
    entry : dict or None
        What the map records of the run that stored the network, under the
        network key of its metator metadata: the version of metaTOR, the
        command line and the parameters of the run, and the columns of the map
        that it read (see metator.cli.network_entry).
    tmp_dir : str, optional
        Where the contig data is written for metacooler annotate
        (contig_data_network.txt): the temporary directory of the run."""

    table: pd.DataFrame
    normalization: str
    entry: dict | None = field(default_factory=dict)
    tmp_dir: str | None = None


def network_columns(table):
    """The columns of the contig data of a network (table, see Network) that the
    map stores (CONTIG_COLUMNS, in that order): those that the run measured.
    The shotgun coverage and the restriction sites are text (-) where it has
    none."""
    return [column for column in CONTIG_COLUMNS if column in table.columns and pd.api.types.is_numeric_dtype(table[column])]


def _store_network(path, network):
    """Store the network of the contigs of a run (network, a Network) in the map
    at path, the staged map of an edit session (see write_map): its measured
    columns in the contigs table, last, which replace those of an earlier run
    as _write_back_columns replaces them, its entry under the network key of
    the metator metadata of the map (unless it keeps the entry of the map: see
    Network), the other keys of which are kept, and the weights of its
    normalization, computed from the columns just written, which replace those
    of an earlier run (see store_normalization), last of the bins table."""
    table = network.table
    measured = network_columns(table)
    metator = dict(_metator_metadata(open_map(path)))
    if network.entry is not None:
        metator["network"] = {**network.entry, "normalization": network.normalization}
    # Computed from the columns that this run replaces (see
    # _drop_earlier_normalizations).
    _drop_earlier_normalizations(path, network.normalization)
    with ExitStack() as scratch:
        directory = network.tmp_dir or scratch.enter_context(tempfile.TemporaryDirectory())
        columns = os.path.join(directory, "contig_data_network.txt")
        table[["Name", *measured]].to_csv(columns, sep="\t", index=False)
        _write_back_columns(path, columns, measured, metadata={"metator": metator})
    if network.normalization != "None":
        _store_normalization(path, network.normalization, _normalization_weights(path, network.normalization))


def _metator_metadata(mc):
    """The metator key of the metadata of the map mc (an empty dict when it has
    none, or none that is a dict)."""
    metator = mc.metadata.get("metator")
    return metator if isinstance(metator, dict) else {}


def _store_metator_entries(path, entries):
    """Set entries (a dict: partition...) in the metator key of the metadata of
    the map at path, the staged map of an edit session, whose other keys, and
    the other keys of the metadata, are kept (see _store_metadata)."""
    _store_metadata(path, {"metator": {**_metator_metadata(open_map(path)), **entries}})


def _drop_metator_entries(path, names):
    """Remove the entries names of the metator key of the metadata of the map at
    path, the staged map of an edit session: what made the clustering rounds of
    those names, which the map no longer holds (see write_map)."""
    if not names:
        return
    metator = _metator_metadata(open_map(path))
    if any(name in metator for name in names):
        _store_metadata(path, {"metator": {key: value for key, value in metator.items() if key not in names}})


def _store_metadata(path, keys):
    """Set keys (a dict) in the metadata of the map at path, the staged map of an
    edit session, whose other keys are kept: a write of metacooler's
    MetaCooler.open("r+"), which the session validates with the rest."""
    import json

    import metacooler

    mc = metacooler.MetaCooler(path)
    metadata = {**mc.metadata, **keys}
    with mc.open("r+") as group:
        group.attrs["metadata"] = json.dumps(metadata)


# The steps of the binning, which make the clustering rounds that a run
# replaces (see store_rounds): a round that another step made (that of
# metator mge), or another tool, is dropped with them.
BINNING_STEPS = ("partition", "validation")

# The clustering rounds that metator partition stores (see write_map),
# which the validation reads (see read_partition).
PARTITION_ROUNDS = (mtcl.CORE, mtcl.OVERLAPPING)

# The key of the partition entry of the metator metadata of a map that
# records the network that the partition was made on: its digest (see
# metator.network.network_digest). The host entry records it too.
NETWORK_DIGEST = "network_digest"

# What metator host stores in a map (see write_map and metator.host): the run,
# under the host entry of its metator metadata, with the digest of the
# network that it read (NETWORK_DIGEST); the hosts of each MGE cluster, as
# columns of the clusters of the round mge; and the MGE clusters of each MAG,
# as columns of the genomes table. They describe the round mge, the genomes
# and the network of the map as the run read them: a run that replaces one of
# them drops them all, with a warning.
HOST_ENTRY = "host"
HOST_CLUSTER_COLUMNS = ("n_hosts", "hosts", "host_shares")
HOST_GENOME_COLUMNS = ("n_mges", "mges")

# What metator markers stores in a map (see write_map and metator.markers):
# the marker genes of each contig, in each marker set of miComplete, as
# columns of its contigs table, and the run, under the markers entry of its
# metator metadata, with the markers and the weights of each set. They are
# facts about the contigs: metator markers, and the step markers of pipeline,
# replace them (metator.cli.REPLACES_MARKERS); every other run keeps them,
# those that regroup the contigs too (validation, mge --reassign), as a
# curation by hand does. Their genes were found in the genomes of the map as
# the run read them, which the entry records (ASSIGNMENT_DIGEST): a run that
# regroups the contigs into other genomes says that the markers predate them
# (see write_map), and so does metator report (see markers_predate); metator
# markers run again finds them in the genomes as they are.
MARKERS_ENTRY = "markers"
MARKER_COLUMNS = ("markers_Bact105", "markers_Arch131")

# The key of the markers entry of the metator metadata of a map that records
# the binning that the marker genes of its contigs were found in: the digest
# of the genome of each contig then (see assignment_digest), as the partition
# and host entries record the network that they read (NETWORK_DIGEST). The
# genes key of the entry names the genomes that pyrodigal was trained on
# (trained_genomes).
ASSIGNMENT_DIGEST = "assignment_digest"


def holds_partition(mc):
    """Whether the map mc holds a partition of its network: the clustering round
    core or overlapping, which metator partition stores (see read_partition)."""
    return any(entry.get("name") in PARTITION_ROUNDS for entry in mc.cluster_rounds)


def partition_network(mc):
    """The digest of the network that the partition of the map mc was made on (see
    metator.network.network_digest), which metator partition records with it,
    under the partition key of the metator metadata of the map (as
    NETWORK_DIGEST): None when the map records none."""
    entry = _metator_metadata(mc).get("partition")
    digest = entry.get(NETWORK_DIGEST) if isinstance(entry, dict) else None
    return digest if isinstance(digest, str) else None


def holds_hosts(mc):
    """Whether the map mc holds the hosts of its MGE clusters that metator host
    found: the host entry of its metator metadata (see HOST_ENTRY)."""
    return isinstance(_metator_metadata(mc).get(HOST_ENTRY), dict)


def host_entry(mc):
    """The host entry of the metator metadata of the map mc (see HOST_ENTRY): what
    the run that found the hosts of its MGE clusters recorded, its command line
    and its parameters among them; None when the map holds no hosts."""
    entry = _metator_metadata(mc).get(HOST_ENTRY)
    return entry if isinstance(entry, dict) else None


def hosts_network(mc):
    """The digest of the network that metator host read to find the hosts that
    the map mc holds (see holds_hosts), as its host entry records it: None when
    the map records none."""
    entry = _metator_metadata(mc).get(HOST_ENTRY)
    digest = entry.get(NETWORK_DIGEST) if isinstance(entry, dict) else None
    return digest if isinstance(digest, str) else None


def holds_markers(mc):
    """Whether the map mc holds the marker genes of its contigs that metator
    markers found: the markers entry of its metator metadata, or a column of
    MARKER_COLUMNS in its contigs table."""
    held = isinstance(_metator_metadata(mc).get(MARKERS_ENTRY), dict)
    return held or any(column in mc.contigs.columns for column in MARKER_COLUMNS)


def markers_entry(mc):
    """The markers entry of the metator metadata of the map mc (see
    MARKERS_ENTRY): what the run that found the marker genes of its contigs
    recorded, its marker sets among them (see metator.markers.map_sets); None
    when the map holds no markers."""
    entry = _metator_metadata(mc).get(MARKERS_ENTRY)
    return entry if isinstance(entry, dict) else None


def map_assignment(mc):
    """The binning of the map mc (opened): the names of its contigs, in the order
    of its contigs table, and the name of the genome of each (UNBINNED_GENOME
    for none), as two arrays (see assignment_digest)."""
    contigs = mc.contigs
    genomes = mc.genomes["name"].astype(str).to_numpy()
    return contigs["name"].astype(str).to_numpy(), genomes[contigs["genome_id"].to_numpy()]


def assignment_digest(contigs, genomes):
    """The identity of a binning of the contigs of a map: the SHA-256 of the
    genome of each contig, by the names of the contigs, as
    metator.network.network_digest is the identity of a network. It is that
    of the assignment, not of a map: the same whatever the order of the
    contigs, another one once a contig is in another genome, or a genome is
    renamed. metator markers records the binning that it found the marker
    genes of the contigs in (see ASSIGNMENT_DIGEST), which markers_predate
    checks the binning of the map against.

    Parameters
    ----------
    contigs : sequence of str
        The names of the contigs of the map.
    genomes : sequence of str
        The name of the genome of each, in the same order: UNBINNED_GENOME for
        none (see map_assignment).

    Returns
    -------
    str
        The SHA-256, in hexadecimal, of the name of each contig and of its
        genome, the contigs by name, each field ended by a NUL, which no name
        holds."""
    import hashlib

    names = np.asarray(contigs, dtype=str)
    of = np.asarray(genomes, dtype=str)
    order = np.argsort(names, kind="stable")
    fields = zip(names[order].tolist(), of[order].tolist(), strict=True)
    return hashlib.sha256("".join(f"{name}\0{genome}\0" for name, genome in fields).encode()).hexdigest()


def markers_binning(mc):
    """The digest of the binning that metator markers found the marker genes of
    the contigs of the map mc in (see assignment_digest), as its markers entry
    records it (ASSIGNMENT_DIGEST): None when the map records none (it holds no
    markers, or markers that an earlier metaTOR found)."""
    entry = markers_entry(mc)
    digest = entry.get(ASSIGNMENT_DIGEST) if entry is not None else None
    return digest if isinstance(digest, str) else None


def markers_predate(mc):
    """Whether the marker genes of the contigs of the map mc predate its genomes,
    as metator report says: found by metator markers in another binning of its
    contigs than the one the map has (the digest that its markers entry records
    is not that of its genomes, see assignment_digest). Each contig has its own
    markers still, but summed over a genome whose contigs changed since, they
    give an estimate of miComplete's numbers, not those of the genome as
    metator markers would find its genes. None when the map cannot tell: it
    records no such digest (see markers_binning)."""
    found_in = markers_binning(mc)
    if found_in is None:
        return None
    return found_in != assignment_digest(*map_assignment(mc))


def _markers_found_in(mc):
    """The digest of the binning that the marker genes that the map mc holds were
    found in (see markers_binning) or, when its markers entry records none,
    that of its binning as it is, the last one that they can have been found
    in; None when the map holds no markers (see holds_markers)."""
    if not holds_markers(mc):
        return None
    return markers_binning(mc) or assignment_digest(*map_assignment(mc))


def _regrouped(contigs, clusters, genomes):
    """The digest of the binning of a map (see assignment_digest) once a run
    regroups its contigs (contigs: its contigs table): into the final clusters
    of a validation (clusters, a FinalClusters), or else into the genomes that
    its round mge makes (genomes, an MgeGenomes); a contig of no genome, or of
    ND, in UNBINNED_GENOME, where metacooler reassign puts it."""
    names = contigs["name"].astype(str).tolist()
    assignment = clusters.assignment if clusters is not None else genomes.assignment
    of = [assignment.get(name, "ND") for name in names]
    return assignment_digest(names, [UNBINNED_GENOME if genome == "ND" else genome for genome in of])


def _log_markers_predate(metacool, outfile, in_place):
    """Say, once outfile is written from the map metacool (in place, or to
    another file), that it keeps the marker genes of the contigs that metator
    markers stored in the map, which the run regrouped into other genomes than
    those that the genes were found in (see markers_predate)."""
    where = f"{metacool} keeps" if in_place else f"{outfile} keeps, of {metacool},"
    logger.info(
        f"{where} the marker genes of its contigs that metator markers found (the columns "
        f"{' and '.join(MARKER_COLUMNS)} of its contigs table), which predate its new genomes: each contig keeps its "
        "own, found in the genome that it was in then, but summed over a genome whose contigs changed, they give an "
        "estimate of miComplete's numbers, which miComplete's own run on that genome can differ from. Run metator "
        f"markers again on {outfile} to find them in its genomes as they are."
    )


@dataclass
class Markers:
    """The marker genes of the contigs of a map, as metator markers stores them
    (see write_map and metator.markers.contig_markers).

    Attributes
    ----------
    table : pandas.DataFrame
        One row per contig of the map: its name (Name) and its markers in each
        set (MARKER_COLUMNS, see metator.markers.encode_markers), an empty text
        for none.
    entry : dict
        What the map records of the run, under the markers key of its metator
        metadata (see metator.cli.markers_entry).
    tmp_dir : str, optional
        Where the table is written for metacooler annotate
        (contig_markers.txt): the temporary directory of the run."""

    table: pd.DataFrame
    entry: dict
    tmp_dir: str | None = None


#: The table of the markers that metator markers stores (see _store_markers),
#: in the temporary directory of the run: a row per contig, and one more,
#: NO_CONTIG, when no contig of the map has a marker of a set.
MARKERS_FILE = "contig_markers.txt"


#: The row of the table of the markers that names no contig (see
#: _store_markers), named so in the log of metacooler annotate, which leaves it
#: out, and in contig_markers.txt; each of its markers is this text.
NO_CONTIG = "(no contig)"
NO_CONTIG_TEXT = "(text)"


def _no_contig(names):
    """A name of no contig of a map (names: its contigs), nor of one by its name
    in the input of a map that metacooler collect made, which annotate matches
    too (the name after <genome>__): NO_CONTIG, or that after more
    underscores."""
    taken = set(names)
    name = NO_CONTIG
    while name in taken or any(other.endswith(f"__{name}") for other in taken):
        name = f"_{name}"
    return name


def _store_markers(path, markers):
    """Store the marker genes of the contigs of a map (markers, a Markers) in the
    map at path, the staged map of an edit session (see write_map): the columns
    MARKER_COLUMNS of its contigs table, last, in place of those of an earlier
    run, with metacooler annotate, and the run under the markers entry of its
    metator metadata, whose other entries are kept.

    Both columns hold text, an empty one for a contig without a marker, even
    when no contig of the map has a marker of a set: metacooler stores a column
    of a table as numbers when each of its values parses as one, an empty text
    as missing, as it would such a column. One row more then, which names no
    contig of the map (see _no_contig) and which annotate leaves out, keeps the
    columns text: every row of a table gives its columns their type, those that
    match no contig too."""
    import metacooler

    mc = open_map(path)
    held = [f"contigs.{column}" for column in MARKER_COLUMNS if column in mc.contigs.columns]
    metator = {**_metator_metadata(mc), MARKERS_ENTRY: markers.entry}
    if held:
        # Dropped, then added again, as a first run adds them: last in the
        # table, whatever an earlier run (or another tool) stored under their
        # names, as text or as numbers.
        metacooler.annotate(path, drop=held, inplace=True)
    table = markers.table[["Name", *MARKER_COLUMNS]]
    no_text = not all((table[column] != "").any() for column in MARKER_COLUMNS)
    if no_text:
        row = {"Name": _no_contig(table["Name"].astype(str)), **dict.fromkeys(MARKER_COLUMNS, NO_CONTIG_TEXT)}
        table = pd.concat([table, pd.DataFrame([row])], ignore_index=True)
    with ExitStack() as scratch:
        directory = markers.tmp_dir or scratch.enter_context(tempfile.TemporaryDirectory())
        path_of_table = os.path.join(directory, MARKERS_FILE)
        table.to_csv(path_of_table, sep="\t", index=False)
        metacooler.annotate(
            path,
            contigs_table=path_of_table,
            contigs_table_format="Name",
            columns=[f"contigs.{column}" for column in MARKER_COLUMNS],
            # The row that names no contig, if any: every other names one.
            unmatched="ignore" if no_text else "warn",
            metadata={"metator": metator},
            inplace=True,
        )


def _later_rounds(mc):
    """The entries of the clustering rounds of the map mc that the binning did not
    make: those of metator mge, or of another tool."""
    return [entry for entry in mc.cluster_rounds if entry.get("step") not in BINNING_STEPS]


def _made_by(entry):
    """What made a clustering round, from its entry: its tool and step (metaTOR
    mge)."""
    return " ".join(str(entry[key]) for key in ("tool", "step") if entry.get(key)) or "a tool it does not name"


def _warn_dropped_rounds(entries, metacool, outfile, in_place):
    """Warn, once outfile is written from the map metacool (in place, or to
    another file), that it does not keep the clustering rounds entries of the
    map, which the run did not make (see _later_rounds): each named with what
    made it, to be run again on the new clusters."""
    if not entries:
        return
    rounds = ", ".join(f"{entry['name']} (by {_made_by(entry)})" for entry in entries)
    where = f"{metacool} no longer holds" if in_place else f"{outfile} does not keep, of {metacool},"
    logger.warning(
        f"{where} the clustering round(s) {rounds}, which this run does not make: a run that replaces the clusters of "
        "a map stores its own rounds only. Run what made them again, on the new clusters, to store them."
    )


def _warn_dropped_hosts(metacool, outfile, in_place, why):
    """Warn, once outfile is written from the map metacool (in place, or to
    another file), that it does not keep the hosts of the MGE clusters that
    metator host stored in the map, which describe what the run replaced
    (why)."""
    where = f"{metacool} no longer holds" if in_place else f"{outfile} does not keep, of {metacool},"
    logger.warning(
        f"{where} the hosts of the MGE clusters that metator host found (the {HOST_ENTRY} entry of its metator "
        f"metadata, the columns {', '.join(HOST_CLUSTER_COLUMNS)} of its round {mtcl.MGE} and "
        f"{', '.join(HOST_GENOME_COLUMNS)} of its genomes table), which described {why}, which this run replaced. "
        f"Run metator host again on {outfile} to find them."
    )


def store_rounds(path, rounds, names=None, name=None):
    """Store the clustering rounds of a run in the /clusters group of a map, in
    place of those it held: the overlapping, recursive and final clusters of
    metator.clusters, with their contigs, their parents and their columns
    (those of miComplete, and of the bin summary for the final clusters).

    A text column holds an empty text where a round has no value (HDF5 stores
    no missing text), and the contigs of a cluster that the map lacks are left
    out, with a warning.

    Parameters
    ----------
    path : str
        Path to the .metacool file, written in place.
    rounds : tuple of pandas.DataFrame, or None
        The clusters and their contigs (see metator.clusters.read_rounds).
    names : list of str, optional
        The names of the rounds, in the order they were made: by default, those
        of the clusters, in their order.
    name : str, optional
        The map that path is staged for, which the warnings name (path by
        default)."""
    import metacooler

    if rounds is None:
        metacooler.drop_clusters(path)
        return
    clusters, members = rounds
    order = list(dict.fromkeys(clusters["round"])) if names is None else list(names)
    steps = {name: "partition" if name in PARTITION_ROUNDS else "validation" for name in order}
    _store_clusters(
        path,
        [{"name": name, "tool": "metaTOR", "version": __version__, "step": steps[name]} for name in order],
        clusters,
        members,
        name=name,
    )


def _store_clusters(path, entries, clusters, members, name=None):
    """Write the /clusters group of the map at path, in place of the one it held:
    the rounds entries (each a dict with its name and what is stored with it),
    and their clusters and members (see store_rounds)."""
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


def _computed_cluster_columns():
    """The columns of the clustering rounds that metacooler computes from their
    members (see metacooler.MetaCooler.clusters), which no tool stores, and
    which metacooler refuses as given."""
    from metacooler.schema import CLUSTER_GENOME_COLUMN

    return ["n_contigs", "length", CLUSTER_GENOME_COLUMN]


def map_rounds(mc):
    """The clustering rounds of a map (mc, opened), as metator.clusters.Rounds
    gives them (see Rounds.tables), and as _store_clusters stored them: every
    cluster of every round, in the order of the rounds and of their clusters,
    with its parent and the columns its round stored, not those that metacooler
    computes; and every contig of a cluster, in the order of the map. An empty
    text of a text column, which _store_clusters stores for a missing value
    (HDF5 holds no missing text), is a missing value again (None): a genuinely
    empty text too.

    Returns
    -------
    tuple of pandas.DataFrame
        The clusters and their contigs (columns round, cluster, contig)."""
    clusters = mc.clusters.drop(columns=_computed_cluster_columns(), errors="ignore").reset_index(drop=True)
    for column in clusters.columns:
        if not pd.api.types.is_numeric_dtype(clusters[column]):
            values = clusters[column].astype(object)
            clusters[column] = values.where(values.notna() & (values.astype(str) != ""), None)
    members = mc.cluster_members().reset_index(drop=True)
    return clusters, members[["round", "cluster", "contig"]]


def read_partition(metacool, mc=None, network=None):
    """The partition that a map stores, as metator partition MAP stores it (see
    write_map and metator.clusters.partition_rounds): its clustering rounds
    core and overlapping, read back as metator.clusters.Partition. With
    network, the digest of the network stored in the map (see
    metator.network.stored_network), the partition must have been made on that
    network, as the map records it (see partition_network): metator network can
    store another network in the map, and leaves the partition there, which is
    then not that of the network.

    Clusters are numbered as their names say (see
    metator.clusters.cluster_number) and contigs get their positions in the
    map, from 1. A core cluster's contigs come in ID order, an overlapping
    cluster's in the order of its core clusters, then of IDs: the order
    metator.partition.partition gives, which the map's sorted rounds lose.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection, for the messages.
    mc : metacooler.MetaCooler, optional
        The map, opened: by default, metacool is opened.
    network : str, optional
        The digest of the network stored in the map, which the partition must
        have been made on.

    Returns
    -------
    metator.clusters.Partition
        The core and the overlapping clusters of the map.

    Raises
    ------
    MetacoolError
        If the map lacks the round core or overlapping (metator partition has
        not run on it), or, with network, does not record that its partition
        was made on that network, or a cluster of them is not named as metator
        partition names it, or holds a contig in no core cluster (an
        overlapping cluster), or the map cannot be read."""
    with _metacooler_logging(), read_errors(metacool):
        if mc is None:
            mc = open_map(metacool)
        held = [entry["name"] for entry in mc.cluster_rounds]
        missing = [name for name in PARTITION_ROUNDS if name not in held]
        if missing:
            raise MetacoolError(
                f"{metacool} holds no partition of its network (the clustering round(s) {' and '.join(missing)}): run "
                "metator partition on it first, which stores its core and overlapping clusters in the map."
            )
        made_on = partition_network(mc)
        if network is not None and made_on != network:
            reason = (
                "does not record the network that its partition was made on, which metator partition records"
                if made_on is None
                else "holds a partition of another network than the one stored in it now: metator network stored "
                "another network in the map since metator partition ran on it"
            )
            raise MetacoolError(
                f"{metacool} {reason}. Run metator partition {metacool} again, which partitions the network stored in the map."
            )
        ids = pd.Index(mc.contigs["name"].astype(str))
        members = {name: mc.cluster_members(name) for name in PARTITION_ROUNDS}
        # Every cluster of the round, those with no contig too.
        clusters = mc.clusters
    found = {}
    for round_name in PARTITION_ROUNDS:
        named = list(clusters.loc[(clusters["round"] == round_name).to_numpy(), "name"].astype(str))
        numbers = {name: mtcl.cluster_number(round_name, name) for name in named}
        odd = [name for name, number in numbers.items() if number is None]
        # Two names of one number (core_00001 and core_000001).
        seen = {}
        for name, number in numbers.items():
            if number is not None and seen.setdefault(number, name) != name:
                odd.append(name)
        if odd:
            form = "core_<number>" if round_name == mtcl.CORE else "<prefix>_<number>_00000"
            raise MetacoolError(
                f"The clustering round {round_name} of {metacool} holds cluster(s) that metator partition does not "
                f"name so ({form}, one number each): {_shown(odd)}. Run metator partition on the map again."
            )
        held_members = members[round_name]
        # The ID of each contig of a cluster: its position in the map, from 1.
        rows = pd.Series(ids.get_indexer(held_members["contig"].astype(str)) + 1)
        by_cluster = rows.groupby(held_members["cluster"].astype(str).to_numpy()).agg(sorted).to_dict()
        found[round_name] = (
            {numbers[name]: [int(contig) for contig in by_cluster.get(name, [])] for name in named},
            {numbers[name]: name for name in named},
        )
    core, _ = found[mtcl.CORE]
    overlapping, names = found[mtcl.OVERLAPPING]
    # The core cluster of each contig, for the order of the overlapping ones.
    core_of = {contig: number for number, contigs in core.items() for contig in contigs}
    lacking = sorted({contig for contigs in overlapping.values() for contig in contigs} - core_of.keys())
    if lacking:
        raise MetacoolError(
            f"The overlapping clusters of {metacool} hold contig(s) in no core cluster: {_shown([ids[i - 1] for i in lacking])}. "
            "Run metator partition on the map again."
        )
    overlapping = {
        number: sorted(contigs, key=lambda contig: (core_of[contig], contig)) for number, contigs in overlapping.items()
    }
    return mtcl.Partition(core=dict(sorted(core.items())), overlapping=dict(sorted(overlapping.items())), names=names)


def _stacked(first, second):
    """The rows of two tables, one after the other, with the columns of both:
    those of first, then the others of second."""
    return second.reset_index(drop=True) if first.empty else pd.concat([first, second], ignore_index=True)


# The keys of the metator metadata that the development builds of metaTOR
# 2.0 wrote at its top, for the run that wrote the map of the clusters, which
# each step now records in its own entry (see write_map): dropped when a
# validation writes its entry.
LEGACY_RUN_KEYS = ("version", "command", "parameters", "columns_from_map")


@dataclass
class FinalClusters:
    """The final clusters of a validation, as the map of the clusters stores them
    (see write_map).

    Attributes
    ----------
    assignment : dict
        The final cluster of each binned contig, by name: the other contigs of
        the map go to metacooler's __unbinned__ genome.
    summary : str
        The summary of the final clusters at full precision
        (cluster_summary_full.txt, in the temporary directory of the run; see
        metator.io.SUMMARY_FULL_FILE): one row per final cluster, named by its
        first column, whose other columns go to the genomes table, with the
        values of miComplete, where cluster_summary.txt rounds them.
    rounds : tuple of pandas.DataFrame
        The clustering rounds of the run, every cluster and its contigs (see
        metator.clusters.read_rounds), which replace those of the map.
    entry : dict
        The entry validation of the metator metadata of the map: what the run
        was, what it read of the map, and the columns of its summary, in order."""

    assignment: dict
    summary: str | None
    rounds: tuple | None
    entry: dict


def write_map(
    metacool,
    outfile,
    *,
    network=None,
    partition=None,
    rounds=None,
    round_names=None,
    clusters=None,
    mge=None,
    reassign=False,
    hosts=None,
    drop_hosts=False,
    markers=None,
    force=False,
    planned=None,
    mge_follows=False,
    markers_follow=False,
):
    """Write the map of a run at outfile from its input map, metacool: what the
    run stores in it, everything else as it was. Each step of a run writes its
    map so, and pipeline writes what its three steps store, in one map, as the
    steps one by one write it (then its step markers stores the marker genes
    of the contigs in it, and, with --mges, its steps mge and host the round
    mge and the hosts, as markers, mge and host do).

    Everything goes into one edit session (see _staged_edit), validated once
    and renamed once; on any failure outfile is left as it was. In order: the
    final clusters (a map written anew from the input map, contigs grouped by
    cluster then longest first), or the genomes that the round mge makes (see
    mge_reassignment: a map written anew from the input map too), the network
    (columns and weights, last in both tables), the partition's rounds and
    entry, then the round mge, or the hosts of the MGE clusters. A run that
    replaces the rounds warns of those the binning did not make, and drops
    them.

    The hosts that metator host stored in the map (see HOST_ENTRY) describe
    its round mge, its genomes and its network as they were: a run that
    replaces the round mge (metator mge, but --random), the rounds
    (partition) or the genomes (validation, pipeline) drops them, and so does
    a network run with drop_hosts (another network than the one they were
    found on), each with a warning.

    The marker genes that metator markers stored in the map (see
    MARKERS_ENTRY) are facts about its contigs, which every write but that of
    markers keeps: a run that regroups the contigs into other genomes than
    those that the genes were found in (validation, mge --reassign) says so
    once the map is written, as they predate its new genomes (see
    markers_predate), unless markers_follow.

    Parameters
    ----------
    metacool : str
        Path to the input .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    outfile : str
        Path to the .metacool file to write.
    network : Network, optional
    partition : dict, optional
        The entry partition of the metator metadata.
    rounds : tuple of pandas.DataFrame, optional
        Clustering rounds, which replace those of the map: the rounds core and
        overlapping of a partition (see metator.clusters.partition_rounds).
    round_names : list of str, optional
        The names of the rounds, in the order they were made (see
        store_rounds).
    clusters : FinalClusters, optional
    mge : tuple, optional
        The tables of the round mge or mge_random (see
        metator.clusters.mge_round) and its entry (its name, and what made it:
        see metator.cli.mge_entry). The round mge undoes the MGE genomes of
        the map (see mge_reassignment).
    reassign : bool
        With the round mge (mge), whether its clusters are genomes of the map
        then (metator mge --reassign).
    hosts : tuple, optional
        The hosts of the MGE clusters of the map (a metator.host.Hosts) and
        the host entry of the run that found them (see metator.cli.host_entry),
        stored alone (see _store_hosts).
    drop_hosts : bool
        Whether to drop the hosts that the map holds, which another network
        than the one stored now gave (a network run).
    markers : Markers, optional
        The marker genes of the contigs of the map (metator markers), stored
        alone (see _store_markers). Every other write keeps those that the
        map holds (see MARKERS_ENTRY).
    force : bool
        Whether to replace an existing outfile: the input map itself, which a
        run updates in place, or another file.
    planned : metacooler.MetaCooler, optional
        The input map as the run read it, opened when the run started: the map
        is refused if it changed since (another run replaced it, or a program
        wrote it), rather than written from another map than the one that the
        run worked on.
    mge_follows : bool
        Whether the run stores the round mge and the hosts of its MGE clusters
        in the map once it is written (pipeline --mges, in its steps mge and
        host): the warnings of what the map no longer holds leave them out,
        which the run makes again.
    markers_follow : bool
        Whether the run says later what becomes of the marker genes of the
        contigs that the map holds, which this write does not say then: it
        finds them again once the map is written (pipeline, in its step
        markers).

    Raises
    ------
    ValueError
        If clusters come with rounds or round_names: the rounds of a validation
        are those of its clusters; or with mge: the genomes of a validation are
        its final clusters; or hosts, or markers, with anything else to store.
        If reassign comes without the round mge.
    MetacoolError
        If the outfile is no file's path (a URI), exists (without force), is
        the .mmetacool collection of the input map or a hard link of the input
        map, the input map is missing, cannot be read or changed since planned
        was opened, a contig of the clusters is not a contig of the input map,
        a cluster name is not ASCII or one that metacooler reads as no cluster
        (NA, nan, -1...), an MGE cluster is named after a genome that the map
        keeps (see mge_reassignment), or metacooler cannot write the map (a
        ValueError, OSError or KeyError, message kept)."""
    if clusters is not None and (rounds is not None or round_names is not None or mge is not None):
        raise ValueError(
            "write_map takes the clustering rounds and the genomes of a validation from its clusters "
            "(FinalClusters): give clusters without rounds, round_names or mge."
        )
    if hosts is not None and any(store is not None for store in (network, partition, rounds, clusters, mge, markers)):
        raise ValueError("write_map stores the hosts of the MGE clusters of metator host alone: give hosts alone.")
    if markers is not None and any(store is not None for store in (network, partition, rounds, clusters, mge)):
        raise ValueError("write_map stores the marker genes of the contigs of metator markers alone: give markers alone.")
    # What the hosts that the map holds describe, which the run replaces (see
    # HOST_ENTRY): none when it keeps them.
    replaced = (
        "its genomes"
        if clusters is not None
        else "its clustering rounds"
        if rounds is not None
        else f"its round {mtcl.MGE}"
        if mge is not None and mge[1]["name"] == mtcl.MGE
        else "its network"
        if drop_hosts
        else None
    )
    if reassign and (mge is None or mge[1]["name"] != mtcl.MGE):
        raise ValueError(f"reassign makes genomes of the clusters of the round {mtcl.MGE}: give them as mge.")
    metacool, outfile = str(metacool), str(outfile)
    if is_entry(metacool) and is_input_map(metacool, outfile):
        raise MetacoolError(
            f"{outfile} is the .mmetacool collection of the input map {metacool}: the map would replace the whole "
            "collection. Write it to another file."
        )
    _refuse_a_hard_link_of_the_map(metacool, outfile)
    _check_paths(outfile, force, [map_file(metacool)])
    if clusters is not None:
        check_ascii_names((), clusters.assignment.values())
        check_bin_names(clusters.assignment.values())
    in_place = is_input_map(metacool, outfile)
    with _metacooler_logging():
        mc = open_map(metacool)
        with read_errors(metacool):
            contigs = mc.contigs
            # The rounds of the map that the run does not make: said once the
            # map is written.
            later = _later_rounds(mc) if clusters is not None or rounds is not None else []
            # The hosts that the run drops: said once the map is written.
            dropped_hosts = replaced is not None and holds_hosts(mc)
            # The genomes of the map once its round mge is stored: None when
            # they stay as they are.
            genomes = mge_reassignment(mc, mge[0], reassign) if mge is not None and mge[1]["name"] == mtcl.MGE else None
            # The binning that the marker genes that the map holds were found
            # in, which a run that regroups the contigs keeps: said once the
            # map is written, if it regroups them into other genomes.
            regroups = clusters is not None or genomes is not None
            found_in = _markers_found_in(mc) if regroups and not markers_follow else None
        if clusters is not None:
            _check_clustered_contigs(metacool, contigs, clusters)
        writing, written = _what_is_written(contigs, network, rounds, clusters, mge, genomes, hosts, markers)
        source = "in place" if in_place else f"from {metacool}"
        logger.info(f"Writing {outfile} {source}{writing}.")
        with _as_metacool_errors(), ExitStack() as session:
            make_map_directory(outfile)
            # One staged copy, which replaces the map once it is valid (see
            # _staged_edit).
            read = _read(mc, planned)
            staged = _enter_session(session, metacool, outfile, read, force)
            if clusters is not None:
                _store_final_clusters(staged, metacool, read, contigs, clusters, partition)
            elif genomes is not None:
                _store_mge_genomes(staged, metacool, read, genomes)
            if network is not None:
                _store_network(staged, network)
            if clusters is not None:
                store_rounds(staged, clusters.rounds, name=metacool)
            else:
                if rounds is not None:
                    store_rounds(staged, rounds, round_names, name=metacool)
                if partition is not None:
                    _store_metator_entries(staged, {"partition": partition})
            if mge is not None:
                _store_round(staged, metacool, *mge)
                entry = mge[1]
                run = {key: value for key, value in entry.items() if key not in ("name", "tool", "step")}
                _store_metator_entries(staged, {entry["name"]: run})
            # What made the rounds of metator mge that the map no longer
            # holds, which the clusters of a validation carried over.
            _drop_metator_entries(staged, [dropped["name"] for dropped in later if _made_by(dropped) == "metaTOR mge"])
            if markers is not None:
                _store_markers(staged, markers)
            if hosts is not None:
                _store_hosts(staged, metacool, *hosts)
            elif replaced is not None:
                # The hosts that described what the run replaced: those that
                # the map held, and what was carried over of them (the host
                # entry that the clusters of a validation carry, the columns of
                # an earlier round mge that the other rounds kept).
                _drop_hosts(staged, metacool)
        logger.debug(f"{outfile} is written: the contact map{written}.")
    _warn_dropped_rounds(
        [entry for entry in later if not (mge_follows and entry["name"] == mtcl.MGE)], metacool, outfile, in_place
    )
    if dropped_hosts and not mge_follows:
        _warn_dropped_hosts(metacool, outfile, in_place, replaced)
    if found_in is not None and _regrouped(contigs, clusters, genomes) != found_in:
        _log_markers_predate(metacool, outfile, in_place)


def _what_is_written(contigs, network, rounds, clusters, mge, genomes=None, hosts=None, markers=None):
    """What the messages of write_map say of the map it writes: while it is
    written, and once it is (each a text that follows the name of the map, and
    the words "the contact map"). genomes: those that the round mge makes (see
    mge_reassignment)."""
    if markers is not None:
        return f", with the marker genes of its {len(contigs)} contigs", ", with the marker genes of its contigs"
    if hosts is not None:
        n_mges = len(hosts[0].clusters)
        return (
            f", with the hosts of the {n_mges} MGE clusters of its round {mtcl.MGE}",
            f", with the hosts of the MGE clusters of its round {mtcl.MGE}",
        )
    if clusters is not None:
        binned = [name for name in contigs["name"].astype(str) if name in clusters.assignment]
        final = set(clusters.assignment[name] for name in binned)
        return (
            f": {len(contigs)} contigs, {len(binned)} of them in {len(final)} clusters",
            " of the clusters",
        )
    if mge is not None:
        tables, entry = mge
        name = entry["name"]
        made = ""
        if genomes is not None:
            mges = (genomes.table[MGE_GENOME_COLUMN] == "True").sum() if MGE_GENOME_COLUMN in genomes.table else 0
            made = f", and {mges} MGE genomes" if mges else ", without its MGE genomes"
        return (
            f", with the {len(tables[0])} clusters of the round {name} of its clustering rounds{made}",
            f", with the round {name} of its clustering rounds{made}",
        )
    if rounds is not None:
        return f", with the {len(rounds[0])} clusters of the run", ", with the clustering rounds of the run"
    if network is not None:
        return ", with the network of its contigs", ", with the network of its contigs"
    return "", ""


def _check_clustered_contigs(metacool, contigs, clusters):
    """Refuse final clusters (clusters, a FinalClusters) with a contig that the
    map metacool lacks (contigs: its contigs table): the clusters of a run are
    clusters of the contigs of its map, by the names that it gives them, and a
    table that read a name as another (00123 as the number 123) gives another
    contig."""
    names = set(contigs["name"].astype(str))
    absent = [name for name in clusters.assignment if name not in names]
    if absent:
        raise MetacoolError(
            f"{len(absent)} contig(s) of the final clusters are not contigs of {metacool}: {_shown(absent)}. The "
            "clusters of a run are clusters of the contigs of its map, by their names in it."
        )


def _store_final_clusters(path, metacool, read, contigs, clusters, partition=None):
    """Write the map of the final clusters of a validation (clusters, a
    FinalClusters) at path, the staged map of an edit session (see write_map),
    anew from the input map metacool with metacooler reassign, which warns of
    what it leaves out. read: the maps opened on the input map (see _read), the
    last one to write the map, whose contigs table is contigs; one replaced or
    written while metacooler read the map is refused."""
    names = [str(name) for name in contigs["name"]]
    metator = {key: value for key, value in _metator_metadata(read[-1]).items() if key not in LEGACY_RUN_KEYS}
    if partition is not None:
        metator["partition"] = partition
    metator["validation"] = clusters.entry
    # Every contig of the map, the unbinned ones as ND (no genome for
    # metacooler): never an empty table, which metacooler refuses.
    assignment = {name: clusters.assignment.get(name, "ND") for name in names}
    final = set(assignment.values()) - {"ND"}
    with ExitStack() as scratch:
        if final and clusters.summary:
            genomes_table = clusters.summary
        else:
            # Without a summary or without a cluster, the cluster names
            # alone: the genome columns that the input map carried (those of
            # an earlier run's clusters) are then replaced by none, where
            # metacooler would keep them and warn that they no longer
            # describe their clusters.
            genomes_table = os.path.join(scratch.enter_context(tempfile.TemporaryDirectory()), "clusters.txt")
            with open(genomes_table, "w") as handle:
                handle.write("".join(f"{name}\n" for name in ["cluster", *sorted(final)]))
        # Grouped by cluster, then from the longest, whatever the order of the
        # input map: the pixels are remapped, the columns of the contigs table
        # and the weights of the bins follow them.
        _reassign(path, metacool, read, assignment, genomes_table, "genome-then-length", {"metator": metator})


def _reassign(path, metacool, read, assignment, genomes_table, layout, metadata=None):
    """Write the map at path, the staged map of an edit session (see write_map),
    anew from the input map metacool with metacooler reassign, which warns of
    what it leaves out: each contig in its genome of assignment (every contig
    of the map, by name: ND for none), the genomes with the values of
    genomes_table (a table that names them in its first column, which replaces
    the values that the input map carried), the contigs laid out as layout
    says (see metacooler.reassign), with metadata (keys of the metadata of the
    map). read: the maps opened on the input map (see _read), the last one to
    write the map; one replaced or written while metacooler read the map is
    refused."""
    import metacooler

    with read_errors(metacool):
        metacooler.reassign(
            metacool,
            assignment,
            # The staged map, written anew from the input map, which
            # metacooler only reads: its history names the input map as the
            # source.
            path,
            layout=layout,
            # The genome names are the first column of the table, unnamed in
            # the summary of the validation (cluster_summary_full.txt).
            genomes_table=genomes_table,
            genomes_table_format="1",
            metadata=metadata,
            # Over the staged copy of the map.
            force=True,
            # On a failure, metacooler removes the partial file it writes next
            # to the staged map.
            keep_partial=False,
        )
        # The input map as the run read it: one replaced or written while
        # metacooler read it, by a program that takes no lock of metacooler,
        # is refused (a read through a map opened on it checks its file).
        for mc in read:
            with mc.open():
                pass


def _store_round(path, metacool, rounds, entry):
    """Replace one clustering round of the map at path, the staged map of an edit
    session (see write_map), and keep everything else as it was: the other
    rounds, its bins, contacts, genomes and contigs tables, QC block and
    metadata. metator mge stores its MGE clusters so, as the round mge, next to
    the rounds of the binning run, whose final clusters the genomes and contigs
    tables keep.

    The round replaces the one of its name and comes last; the other rounds are
    kept, but for the columns metacooler computes again. A column only one
    round has is missing in the others, so whole numbers or booleans there
    become decimals or text."""
    clusters, members = rounds
    name = entry["name"]
    mc = open_map(path)
    held = mc.clusters
    held_members = mc.cluster_members()
    held_rounds = mc.cluster_rounds
    # The columns that metacooler computes, which it refuses as given.
    other = held[(held["round"] != name).to_numpy()].drop(columns=_computed_cluster_columns())
    other_members = held_members[(held_members["round"] != name).to_numpy()]
    entries = [stored for stored in held_rounds if stored["name"] != name]
    _store_clusters(path, [*entries, entry], _stacked(other, clusters), _stacked(other_members, members), name=metacool)


def _store_hosts(path, metacool, hosts, entry):
    """Store the hosts of the MGE clusters that metator host found (hosts, a
    metator.host.Hosts) in the map at path, the staged map of an edit session
    (see write_map), in place of those of an earlier run (see _drop_hosts):
    HOST_CLUSTER_COLUMNS in the clusters of the round mge (missing in those of
    the other rounds), HOST_GENOME_COLUMNS in the genomes table (missing for
    the genomes that are no MAG), and the run under the host entry of the
    metator metadata of the map (entry). The rounds keep their order, their
    entries, the columns of their clusters and their contigs, and the rest of
    the map is left as it was."""
    import metacooler

    _drop_hosts(path, metacool)
    mc = open_map(path)
    clusters = mc.clusters.drop(columns=_computed_cluster_columns())
    rows = (clusters["round"] == mtcl.MGE).to_numpy()
    values = hosts.clusters.set_index("name").reindex(clusters.loc[rows, "name"].astype(str))
    for column in HOST_CLUSTER_COLUMNS:
        found = values[column].to_numpy()
        if pd.api.types.is_numeric_dtype(values[column]):
            stored = np.full(len(clusters), np.nan)
        else:
            stored = np.full(len(clusters), None, dtype=object)
        stored[rows] = found
        clusters[column] = stored
    members = mc.cluster_members()[["round", "cluster", "contig"]]
    _store_clusters(path, mc.cluster_rounds, clusters, members, name=metacool)
    with tempfile.TemporaryDirectory() as directory:
        table = os.path.join(directory, "genomes.txt")
        hosts.genomes.rename(columns={"name": "genome"}).to_csv(table, sep="\t", index=False)
        metacooler.annotate(path, genomes_table=table, genomes_table_format="genome", inplace=True)
    _store_metator_entries(path, {HOST_ENTRY: entry})


def _drop_hosts(path, metacool=None):
    """Remove what metator host stored in the map at path, the staged map of an
    edit session (see write_map): the host entry of its metator metadata, the
    columns HOST_CLUSTER_COLUMNS of its clusters (those of the round mge, kept
    by the other rounds when a run replaced it) and HOST_GENOME_COLUMNS of its
    genomes table. Everything else is kept as it was. metacool: the map that
    path is staged for, which the warnings name."""
    import metacooler

    # Everything read before anything is written: a write of the map replaces
    # the file that mc opened.
    mc = open_map(path)
    clusters = mc.clusters
    stored = [column for column in HOST_CLUSTER_COLUMNS if column in clusters.columns]
    genomes = [f"genomes.{column}" for column in HOST_GENOME_COLUMNS if column in mc.genomes.columns]
    if stored:
        kept = clusters.drop(columns=[*_computed_cluster_columns(), *stored])
        members = mc.cluster_members()[["round", "cluster", "contig"]]
        _store_clusters(path, mc.cluster_rounds, kept, members, name=metacool)
    if genomes:
        metacooler.annotate(path, drop=genomes, inplace=True)
    _drop_metator_entries(path, [HOST_ENTRY])


def mge_genomes(mc):
    """The MGE genomes of the map mc (opened), by name, in the order of its genomes
    table: those that its column mge says are (True), which metator mge
    --reassign made of the clusters of its round mge (see write_map)."""
    genomes = mc.genomes
    if MGE_GENOME_COLUMN not in genomes.columns:
        return []
    flagged = genomes[MGE_GENOME_COLUMN].astype(object).map(lambda value: str(value) == "True").to_numpy(bool)
    return list(genomes["name"].astype(str).to_numpy()[flagged])


@dataclass
class MgeGenomes:
    """The genomes of a map once metator mge stores its round mge in it (see
    mge_reassignment and write_map).

    Attributes
    ----------
    assignment : dict
        The genome of each contig of the map, by name: ND for none.
    table : pandas.DataFrame
        The values of each genome, indexed by its name: the genomes table that
        metacooler reassign stores (see _write_genomes_table)."""

    assignment: dict
    table: pd.DataFrame


# The columns of the clustering rounds that name a cluster and its parent: the
# others are the values of the cluster.
CLUSTER_KEY_COLUMNS = ("round", "name", "parent_round", "parent")


def mge_reassignment(mc, rounds, reassign):
    """The genomes of the map mc (opened) once metator mge stores the round mge of
    its run in it (rounds: its clusters and their contigs, see
    metator.clusters.mge_round), or None when the run leaves them as they are.

    The MGE genomes that the map holds are undone first: those that its column
    mge says are (see mge_genomes), and those named after a cluster of its
    round mge (which metacooler reassign --round mge makes, say). Their
    contigs go back to their final cluster (the round final), or to none.
    With reassign (--reassign), each MGE cluster of the run is a genome then,
    named after it, with its values (GC, Hit, BinningScore): its contigs leave
    the genomes that held them, MAGs or none. The other genomes are kept, with
    their values, those of the validation, which describe them with the
    contigs that they lose, as the log says; a genome whose contigs all leave
    it is no genome any more, and its values go with it. The column mge of
    the genomes table says which genomes are MGE genomes (True): it is not
    there without reassign.

    Returns
    -------
    MgeGenomes, or None
        None without reassign, on a map that holds no MGE genome.

    Raises
    ------
    MetacoolError
        If a genome that the run keeps is named after an MGE cluster of the
        run: their contigs would make one genome."""
    from metacooler.qc import strip_qc_columns

    clusters, members = rounds
    contigs = mc.contigs
    genomes = strip_qc_columns(mc.genomes)
    names = genomes["name"].astype(str).to_numpy()
    held = dict(zip(contigs["name"].astype(str), names[contigs["genome_id"].to_numpy()], strict=True))
    stored = [str(entry.get("name")) for entry in mc.cluster_rounds]
    earlier = set(mge_genomes(mc))
    if mtcl.MGE in stored:
        earlier |= set(mc.cluster_members(mtcl.MGE)["cluster"].astype(str)) & set(names)
    if not earlier and not reassign:
        return None
    final = {}
    if mtcl.FINAL in stored:
        members_final = mc.cluster_members(mtcl.FINAL)
        final = dict(zip(members_final["contig"].astype(str), members_final["cluster"].astype(str), strict=True))
    base = {contig: final.get(contig, UNBINNED_GENOME) if genome in earlier else genome for contig, genome in held.items()}
    moved = {}
    if reassign:
        pairs = zip(members["contig"].astype(str), members["cluster"].astype(str), strict=True)
        moved = {contig: cluster for contig, cluster in pairs if contig in held}
    mges = list(dict.fromkeys(moved.values()))
    clash = sorted({genome for contig, genome in base.items() if contig not in moved} & set(mges))
    if clash:
        raise MetacoolError(
            f"The MGE cluster(s) {_shown(clash)} of the run are named after genome(s) of the map that the run keeps: "
            "their contigs would make one genome. Rename those genomes (metacooler reassign) first."
        )
    assignment = {**base, **moved}
    kept = [genome for genome in dict.fromkeys(assignment.values()) if genome != UNBINNED_GENOME and genome not in mges]
    _log_mge_genomes(held, base, moved, earlier, kept, set(names))
    # The values of the genomes that the run keeps, but those that the MGE
    # genomes of the run get and those that metacooler computes (their
    # length, contigs and bins); none for a genome that gets contigs back
    # after it lost them all.
    columns = [column for column in clusters.columns if column not in CLUSTER_KEY_COLUMNS]
    values = genomes.set_index(pd.Index(names, dtype=object))
    computed = {"name", "length", "n_contigs", "n_bins", MGE_GENOME_COLUMN, *columns}
    table = values.drop(columns=[column for column in values.columns if column in computed]).reindex(
        pd.Index([*kept, *mges], dtype=object)
    )
    if reassign:
        cluster_values = clusters.set_index(pd.Index(clusters["name"].astype(str), dtype=object))
        for column in columns:
            table[column] = cluster_values[column].reindex(table.index)
        table[MGE_GENOME_COLUMN] = ["False"] * len(kept) + ["True"] * len(mges)
    return MgeGenomes({contig: "ND" if genome == UNBINNED_GENOME else genome for contig, genome in assignment.items()}, table)


def _log_mge_genomes(held, base, moved, earlier, kept, genomes):
    """Log what mge_reassignment makes of the genomes of a map: the MGE genomes
    that it undoes (earlier), the genomes that lose contigs to the MGE
    genomes of the run (moved), those that lose them all, and the genomes
    that get contigs back without values (none in the genomes of the map)."""
    undone = sorted(earlier & set(held.values()))
    if undone:
        back = [contig for contig, genome in held.items() if genome in earlier]
        restored = sum(base[contig] != UNBINNED_GENOME for contig in back)
        logger.info(
            f"The {len(undone)} MGE genome(s) of the map are undone ({_shown(undone)}): of their {len(back)} contig(s), "
            f"{restored} go back to their final cluster, and {len(back) - restored} to none."
        )
    if moved:
        lost = Counter(base[contig] for contig in moved if base[contig] != UNBINNED_GENOME)
        taken = sum(lost.values())
        logger.info(
            f"The {len(set(moved.values()))} MGE clusters of the run are genomes of the map, of {len(moved)} contig(s): "
            f"{taken} of them leave {len(lost)} genome(s), {len(moved) - taken} were in none."
        )
        left = [genome for genome in lost if genome in kept]
        if left:
            listed = ", ".join(f"{genome} ({lost[genome]})" for genome in left[:MAX_SHOWN])
            more = f", ... (+{len(left) - MAX_SHOWN} more)" if len(left) > MAX_SHOWN else ""
            logger.info(
                f"{len(left)} genome(s) lose contigs to MGE genomes, and keep the values of their validation, which "
                f"describe them with those contigs: {listed}{more}."
            )
        gone = sorted(genome for genome in lost if genome not in kept)
        if gone:
            logger.warning(
                f"{len(gone)} genome(s) of the map are no genome any more, as all their contigs are in MGE genomes: "
                f"{_shown(gone)}. Their values leave the genomes table with them."
            )
    bare = [genome for genome in kept if genome not in genomes]
    if bare:
        logger.warning(
            f"{len(bare)} final cluster(s) are genomes of the map again, without values: {_shown(bare)}. Their values "
            "left the genomes table with their last contigs, which an MGE genome took: run metator validation again "
            "to measure them."
        )


def _write_genomes_table(table, path):
    """Write the values of the genomes of a map (table, indexed by their names) at
    path, as metacooler reads a genomes table (see _reassign): tab-separated,
    unquoted, the names in the first column, unnamed, the values of a column
    as they are (a decimal as the shortest text that reads it back), a
    missing value as an empty field. A table without columns is the names
    alone, under the header genome. metacooler reads a decimal back with
    pandas: to the bit with pandas 3, to about one part in 10^13 with pandas
    2, whose to_numeric rounds some long decimals off.

    Raises
    ------
    MetacoolError
        If a value holds a tab or a line break, which no field of such a table
        holds."""
    import csv

    if table.columns.empty:
        with open(path, "w") as handle:
            handle.write("".join(f"{name}\n" for name in ["genome", *table.index]))
        return
    try:
        table.rename_axis(None).to_csv(path, sep="\t", quoting=csv.QUOTE_NONE)
    except csv.Error as error:
        raise MetacoolError(
            f"A value of the genomes table of the map holds a tab or a line break, which a table of genomes cannot "
            f"hold ({error}): remove it from the map first (metacooler annotate)."
        ) from error


def _store_mge_genomes(path, metacool, read, genomes):
    """Write the map at path, the staged map of an edit session (see write_map),
    anew from the input map metacool with metacooler reassign: its contigs in
    the genomes of genomes (an MgeGenomes, see mge_reassignment), with their
    values, laid out as the input map is, grouped by genome again when it
    grouped them (metacooler reassign's layout auto), as a validation does.
    read: the maps opened on the input map (see _read)."""
    with tempfile.TemporaryDirectory() as directory:
        table = os.path.join(directory, "genomes.txt")
        _write_genomes_table(genomes.table, table)
        _reassign(path, metacool, read, genomes.assignment, table, "auto")
