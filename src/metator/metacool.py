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
pairs of the Hi-C libraries): the map is the interface between its steps,
and its only output but the FASTA files of its clusters. A run step by step
stores the network of its contigs in the map first, which the next steps
read, then its partition, as the clustering rounds core and overlapping,
which its validation reads (see read_partition), then its validation
writes the map of the clusters: its contigs grouped by final cluster, the
summary of the clusters, and every round of the run. pipeline stores what
the three steps store, in one map (see write_map).

The network is stored as what metaTOR measured on each contig, in the
contigs table, and as the normalization of the network, a weight column of
the bins table (see store_normalization): the raw contacts stay as they
are, and the normalized contact of two bins is their count times both
weights, the way cooler stores its balancing weights.

Every map is written through metacooler, in an edit session
(metacooler.staged_edit): all the edits of a run go into one copy of the
map, which metacooler validates in full once, then renames over the map, or
onto the map that the run writes instead (--output-map). The map is whole
whatever stops the run, and does not grow from run to run. A run
holds the lock of the map it writes from its checks to its end
(metacooler.lock, see metator.cli.hold_lock), which the session re-enters.
A map is written in place when the file to write is the file of the input
map, its symbolic links followed (see is_input_map); any other file is a
map of its own, a symbolic link there replaced, not its target, and a hard
link of the input map is refused.

Core functions to write the file:
    - write_map, the map of a run: its network, its partition, its final
      clusters, or the MGE clusters of metator mge
    - store_normalization, the normalization of the network in a map

And to read it: read_partition, the partition that a map stores.
"""

import logging
import os
import tempfile
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

# The columns of contig_data_final.txt stored in the contigs table: what
# metaTOR measured on each contig. Never a bin column (Final_cluster, the bin IDs
# of the steps): the genome of a contig is its genome_id, and metacooler
# reassign carries the contigs columns into the map of another binning, where
# a bin column would be stale.
CONTIG_COLUMNS = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]

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
def _metacooler_logging(level=None):
    """Send the messages of metacooler to the handlers of metaTOR (terminal
    and log file) while the file is written: those of at least the level of
    metaTOR's logger, or of level when it is higher (logging.WARNING: its
    warnings, but not its accounts of what it reads)."""
    mc_logger = logging.getLogger("metacooler")
    saved = (mc_logger.level, mc_logger.propagate)
    handlers = [handler for handler in logger.handlers if handler not in mc_logger.handlers]
    for handler in handlers:
        mc_logger.addHandler(handler)
    mc_logger.setLevel(logger.level if level is None else max(level, logger.level))
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
    with _metacooler_logging(), read_errors(metacool):
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


def _named(path):
    """path, absolute, the symbolic links of its directory followed but not
    its own name's: the file that a write to path replaces (a symbolic link
    there is replaced, not written through)."""
    path = os.path.abspath(path)
    return os.path.join(os.path.realpath(os.path.dirname(path)), os.path.basename(path))


def is_input_map(metacool, outfile):
    """Whether writing the map at outfile updates the map metacool in place:
    outfile, its own name not followed (see _named), names the file of the
    map, its symbolic links followed, as metacooler edits the file that a
    link points to. A symbolic link at outfile is a file of its own, which a
    write replaces, whatever it points to; and a map given as a symbolic
    link is updated in place through its target, not through outfile =
    the link."""
    source = map_file(metacool)
    return os.path.exists(source) and _named(outfile) == os.path.realpath(source)


def is_hard_link_of(metacool, outfile):
    """Whether outfile is the file of the map metacool under another name (a
    hard link: the same file, which is not the map's own name, nor a
    symbolic link). Written as a map of its own, it would replace that name
    only; metacooler refuses an output that is the same file as its input."""
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
    if os.path.lexists(outfile) and not force:
        raise MetacoolError(f"{outfile} already exists. Use --force to overwrite it.")
    missing = [path for path in inputs if path and not os.path.isfile(path)]
    if missing:
        raise MetacoolError(f"Input file(s) not found: {', '.join(map(str, missing))}.")


def make_map_directory(outfile, named=None):
    """Make the directory of the map at outfile, where metacooler locks and
    stages the map before it renames it there. A directory that cannot be
    made is a MetacoolError, which names outfile as named says (outfile
    itself by default)."""
    try:
        os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
    except OSError as error:
        raise MetacoolError(
            f"The directory of {named or outfile} cannot be made ({error.strerror}): write the map elsewhere."
        ) from error


def _enter_session(session, metacool, outfile, planned, force):
    """Enter, in session (an ExitStack), the edit session of metacooler that
    writes the map of the run at outfile (see _staged_edit), and return its
    staged map. A map that changed since it was read for this run (planned:
    the maps opened on it) is refused first, with metaTOR's reason: the
    clusters of the run come from the map as it was. metacooler's own
    refusals keep its reason (a map that a program wrote while metacooler
    copied it, for example)."""
    if any(_changed(mc) for mc in planned):
        raise MetacoolError(
            f"{metacool} was replaced or written since it was read for this run (by another run, or a program that "
            "takes no lock of metacooler): the clusters of the run come from the map as it was. Run it again on the "
            "map as it is now."
        )
    return session.enter_context(_staged_edit(metacool, outfile, planned, force))


def _changed(mc):
    """Whether the file of the map mc (a metacooler.MetaCooler) was replaced
    or written since mc was opened, as metacooler checks it on every read."""
    import metacooler

    try:
        with mc.open():
            pass
    except metacooler.MapChangedError:
        return True
    return False


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
    ends. When outfile names the file of the input map (see is_input_map),
    the staged map then replaces the input map (the file that a symbolic
    link points to, which keeps pointing at it), under the lock of the map,
    which the run holds, or which the session takes; otherwise it is
    renamed onto outfile, a map of its own, which it replaces with force (a
    symbolic link there is replaced, not followed): the copy is of the map
    alone, so an entry of a .mmetacool collection becomes a map of its own.
    The session refuses a map that changed since planned, the input map as
    the run read it, was opened."""
    import metacooler

    if is_input_map(metacool, outfile):
        return metacooler.staged_edit(metacool, planned=planned)
    return metacooler.staged_edit(metacool, output=outfile, force=force, planned=planned)


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
        Mode of normalization (see metator.network.normalization_weights),
        but "None".

    Raises
    ------
    MetacoolError
        If the normalization is unknown, the contigs table of the map lacks
        a column it reads, the map is replaced or written as it is read,
        or metacooler cannot store the column (a ValueError, OSError or
        KeyError, message kept).
    """
    with _metacooler_logging(), _as_metacool_errors(), read_errors(metacool):
        _store_normalization(metacool, normalization, _normalization_weights(metacool, normalization))


def _normalization_weights(metacool, normalization):
    """What store_normalization checks and computes before it writes
    anything: the map, opened (a metacooler.MetaCooler), the weights of its
    contigs, and the columns of its contigs table that they were computed
    from. Raises MetacoolError when the normalization cannot be stored (see
    store_normalization); a run then writes its map without it."""
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


def _drop_earlier_normalizations(path, normalization=None):
    """Remove from a map the normalizations that metaTOR stored in it (see
    store_normalization): the weight columns of its bins whose normalization
    attribute, which metacooler's format documents for a weight column, is
    metator. They were computed from the columns that the run replaces. Left
    in the map, metacooler annotate would drop them all the same, with a
    warning to compute them again, which the run then does for
    normalization, the normalization that it stores: the others are said to
    be removed. The map is opened for writing only when it holds one."""
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
    """An attribute of an HDF5 object, as text when it is stored as bytes
    (None when the object has none of that name)."""
    value = node.attrs.get(name)
    return value.decode() if isinstance(value, bytes) else value


def _write_back_columns(path, contig_data, contig_columns, header, metadata=None):
    """Give a map the columns of the network of this run, as metacooler
    annotate edits them: the map that metator network, or pipeline, stores
    the network of its contigs in (see _store_network). metadata: keys of
    the metadata of the map that metacooler annotate sets as it adds the
    columns.

    The columns of the contig data that metaTOR measured (contig_columns,
    from the table contig_data) are stored last, in their order: the columns
    of a contig data file (header: all of them) and those of CONTIG_COLUMNS
    that the map carried are dropped first, then those of this run are
    added. The dropped columns are also those that this run does not store:
    one of CONTIG_COLUMNS that it had no values for (Shotgun_coverage
    without a depth file, nor a column of the input map that gives every
    contig its coverage...), which an earlier network left.
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
        metadata=metadata,
        inplace=True,
    )


@dataclass
class Network:
    """The network of the contigs of a run, as a map stores it (see
    write_map).

    Attributes
    ----------
    table : pandas.DataFrame
        The contig data of the contigs of the map, in its order (see
        metator.network.contig_table): its columns that the run measured
        (GC_content and Hit, and Shotgun_coverage and Restriction_site when
        it has them) go to the contigs table of the map.
    normalization : str
        The normalization of the network, stored as the weights
        metator_<normalization> of the bins of the map (see
        store_normalization): none for "None", the raw counts.
    entry : dict
        What the map records of the run that stored the network, under the
        network key of its metator metadata: the version of metaTOR, the
        command line and the parameters of the run, and the columns of the
        map that it read. The normalization is recorded with it.
    tmp_dir : str, optional
        Where the contig data is written for metacooler annotate
        (contig_data_network.txt): the temporary directory of the run. By
        default, a temporary directory of the system, removed at once.
    """

    table: pd.DataFrame
    normalization: str
    entry: dict = field(default_factory=dict)
    tmp_dir: str | None = None


def network_columns(table):
    """The columns of the contig data of a network (table, see Network) that
    the map stores (CONTIG_COLUMNS, in that order): those that the run
    measured. The shotgun coverage and the restriction sites are text (-)
    where it has none."""
    return [column for column in CONTIG_COLUMNS if column in table.columns and pd.api.types.is_numeric_dtype(table[column])]


def _store_network(path, network):
    """Store the network of the contigs of a run (network, a Network) in the
    map at path, the staged map of an edit session (see write_map): its
    measured columns in the contigs table, which replace those of an earlier
    run as _write_back_columns replaces them, its entry under the network
    key of the metator metadata of the map, the other keys of which are
    kept, and the weights of its normalization, computed from the columns
    just written, which replace those of an earlier run (see
    store_normalization). A normalization that cannot be stored fails the
    session: the network of the map is its weights."""
    table = network.table
    measured = network_columns(table)
    metator = dict(_metator_metadata(open_map(path)))
    metator["network"] = {**network.entry, "normalization": network.normalization}
    # Computed from the columns that this run replaces (see
    # _drop_earlier_normalizations).
    _drop_earlier_normalizations(path, network.normalization)
    with ExitStack() as scratch:
        directory = network.tmp_dir or scratch.enter_context(tempfile.TemporaryDirectory())
        columns = os.path.join(directory, "contig_data_network.txt")
        table[["Name", *measured]].to_csv(columns, sep="\t", index=False)
        _write_back_columns(path, columns, measured, ["Name", *measured], metadata={"metator": metator})
    if network.normalization != "None":
        _store_normalization(path, network.normalization, _normalization_weights(path, network.normalization))


def _metator_metadata(mc):
    """The metator key of the metadata of the map mc (an empty dict when it
    has none, or none that is a dict)."""
    metator = mc.metadata.get("metator")
    return metator if isinstance(metator, dict) else {}


def _store_metator_entries(path, entries):
    """Set entries (a dict: partition...) in the metator key of the metadata
    of the map at path, the staged map of an edit session, whose other keys,
    and the other keys of the metadata, are kept: a write of metacooler's
    MetaCooler.open("r+"), which the session validates with the rest."""
    import json

    import metacooler

    mc = metacooler.MetaCooler(path)
    metadata = dict(mc.metadata)
    metadata["metator"] = {**_metator_metadata(mc), **entries}
    with mc.open("r+") as group:
        group.attrs["metadata"] = json.dumps(metadata)


# The steps of the binning, which make the clustering rounds that a run
# replaces (see store_rounds): a round that another step made (that of
# metator mge --map), or another tool, is dropped with them.
BINNING_STEPS = ("partition", "validation")

# The clustering rounds that metator partition stores (see write_map),
# which the validation reads (see read_partition).
PARTITION_ROUNDS = (mtcl.CORE, mtcl.OVERLAPPING)

# The key of the partition entry of the metator metadata of a map that
# records the network that the partition was made on: its digest (see
# metator.network.network_digest).
NETWORK_DIGEST = "network_digest"


def holds_partition(mc):
    """Whether the map mc holds a partition of its network: the clustering
    round core or overlapping, which metator partition stores (see
    read_partition)."""
    return any(entry.get("name") in PARTITION_ROUNDS for entry in mc.cluster_rounds)


def partition_network(mc):
    """The digest of the network that the partition of the map mc was made
    on (see metator.network.network_digest), which metator partition records
    with it, under the partition key of the metator metadata of the map (as
    NETWORK_DIGEST): None when the map records none."""
    entry = _metator_metadata(mc).get("partition")
    digest = entry.get(NETWORK_DIGEST) if isinstance(entry, dict) else None
    return digest if isinstance(digest, str) else None


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
    the step that made it: partition for the rounds core and overlapping,
    validation for the others. A text column holds an empty text where a round has no
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
        no cluster: those of a partition of a network without an edge.
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
    steps = {name: "partition" if name in PARTITION_ROUNDS else "validation" for name in order}
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


def read_partition(metacool, mc=None, network=None):
    """The partition that a map stores, as metator partition MAP stores it
    (see write_map and metator.clusters.partition_rounds): its clustering
    rounds core and overlapping, read back as metator.clusters.Partition.
    With network, the digest of the network stored in the map (see
    metator.network.stored_network), the partition must have been made on
    that network, as the map records it (see partition_network): metator
    network can store another network in the map, and leaves the partition
    there, which is then not that of the network.

    The clusters are numbered as their names say (core_<number>, and
    <prefix>_<number>_00000: see metator.clusters.cluster_number), and their
    contigs are given their IDs, their positions in the map, from 1, as the
    network of the map gives them (see metator.network.map_network). The
    contigs of a core cluster come in the order of their IDs, and those of an
    overlapping cluster in the order of its core clusters, by number, then
    of their IDs: as metator.partition.partition lists them, which the
    rounds of the map, sorted by metacooler, do not keep.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection, for the messages.
    mc : metacooler.MetaCooler, optional
        The map, opened: by default, metacool is opened.
    network : str, optional
        The digest of the network stored in the map, which the partition must
        have been made on. By default, the partition is read whatever network
        it was made on.

    Returns
    -------
    metator.clusters.Partition
        The core and the overlapping clusters of the map.

    Raises
    ------
    MetacoolError
        If the map lacks the round core or overlapping (metator partition has
        not run on it), or, with network, does not record that its partition
        was made on that network, or a cluster of them is not named as
        metator partition names it, or holds a contig in no core cluster (an
        overlapping cluster), or the map cannot be read.
    """
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
    """The rows of two tables, one after the other, with the columns of
    both: those of first, then the others of second."""
    return second.reset_index(drop=True) if first.empty else pd.concat([first, second], ignore_index=True)


# The keys of the metator metadata that the development builds of metaTOR
# 2.0 wrote at its top, for the run that wrote the map of the clusters, which
# each step now records in its own entry (see write_map): dropped when a
# validation writes its entry.
LEGACY_RUN_KEYS = ("version", "command", "parameters", "columns_from_map")


@dataclass
class FinalClusters:
    """The final clusters of a validation, as the map of the clusters stores
    them (see write_map).

    Attributes
    ----------
    assignment : dict
        The final cluster of each binned contig, by name: the other contigs
        of the map go to metacooler's __unbinned__ genome.
    summary : str
        The summary of the final clusters (cluster_summary.txt, in the
        temporary directory of the run): one row per final cluster, named by
        its first column, whose other columns go to the genomes table. None:
        the genomes table holds the names of the clusters alone.
    rounds : tuple of pandas.DataFrame
        The clustering rounds of the run, every cluster and its contigs (see
        metator.clusters.read_rounds), which replace those of the map. None:
        the map keeps no round.
    entry : dict
        The entry validation of the metator metadata of the map: what the
        run was, what it read of the map, and the columns of its summary, in
        order.
    """

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
    force=False,
    planned=None,
):
    """Write the map of a run at outfile from its input map, metacool: what
    the run stores in it, everything else as it was. Each step of a run
    writes its map so, and pipeline writes what its three steps store, in
    one map, as the steps one by one write it.

    The map is written in one edit session of metacooler (see _staged_edit):
    one copy of the map, where every edit of the run goes, which metacooler
    validates in full once, then renames over the input map, when outfile
    names the file of the input map (see is_input_map: given as a symbolic
    link, the input map is updated where it points), or onto outfile, a map
    of its own (a symbolic link there is replaced, not followed; a hard link
    of the input map is refused). An entry of a .mmetacool collection is
    written as a map of its own. On any failure, outfile is left as it was,
    and the staged copy is removed.

    What is stored, in this order:

    - network (metator network, pipeline): the network of the contigs, what
      metaTOR measured on each of them in the contigs table, and the weights
      of its normalization in the bins table, with the entry network of the
      metator metadata (see _store_network). The columns and the weights of
      an earlier network are replaced.
    - partition and rounds (metator partition): the entry partition of the
      metator metadata, what the partition was, and the core and
      overlapping clusters, every one with its contigs, as the rounds core
      and overlapping (round_names: the names of the rounds, in order) of
      the clustering rounds of the map, in place of those it held (see
      store_rounds). The validation reads them back (see read_partition).
    - clusters (metator validation, pipeline): the final clusters. Each
      contig goes to its final cluster, the unbinned ones to __unbinned__,
      and the contigs are grouped by cluster, then from the longest
      (metacooler reassign of the copy in place, layout genome-then-length):
      the contigs table keeps the columns of the network, and the bins table
      its weights, which follow their contigs and bins. The genomes table
      holds the columns of the summary of the clusters. The clustering rounds
      of the run replace those of the map, with the entry validation (and
      partition, for a pipeline run) of the metator metadata. The QC block of
      the input map, when it has one, is carried over to the clusters:
      metacooler keeps the counts within each contig and counts the contacts
      between contigs again.
    - mge (metator mge): the round mge, (its tables, its entry), which
      replaces the round of that name of the map, after the others (see
      _store_round).

    A run that replaces the clustering rounds of the map (partition,
    validation) does not keep the rounds that the binning did not make
    (that of metator mge, or of another tool): it warns of them, once the
    map is written.

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
        Clustering rounds, which replace those of the map: the rounds core
        and overlapping of a partition (see metator.clusters.partition_rounds).
    round_names : list of str, optional
        The names of the rounds, in the order they were made (see
        store_rounds).
    clusters : FinalClusters, optional
    mge : tuple, optional
        The tables of the round mge (see metator.clusters.mge_round) and its
        entry (its name, and what made it).
    force : bool
        Whether to replace an existing outfile: the input map itself, which a
        run updates in place, or another file.
    planned : metacooler.MetaCooler, optional
        The input map as the run read it, opened when the run started: the
        map is refused if it changed since (another run replaced it, or a
        program wrote it), rather than written from another map than the one
        that the run worked on. By default, the map as this function reads
        it.

    Raises
    ------
    MetacoolError
        If the outfile is no file's path (a URI), exists (without force), is
        the .mmetacool collection of the input map or a hard link of the
        input map, the input map is missing, cannot be read or changed since
        planned was opened, a cluster name is not ASCII or one that
        metacooler reads as no cluster (NA, nan, -1...), or metacooler cannot
        write the map (a ValueError, OSError or KeyError, message kept). An
        error of another type is raised as it is.
    """
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
        writing, written = _what_is_written(contigs, network, rounds, clusters, mge)
        logger.info(f"Writing {outfile} from the contact map {metacool}{writing}.")
        with _as_metacool_errors(), ExitStack() as session:
            make_map_directory(outfile)
            # One staged copy, which replaces the map once it is valid (see
            # _staged_edit).
            staged = _enter_session(session, metacool, outfile, _read(mc, planned), force)
            if network is not None:
                _store_network(staged, network)
            if clusters is not None:
                _store_final_clusters(staged, metacool, clusters, partition)
            else:
                if rounds is not None:
                    store_rounds(staged, rounds, round_names, name=metacool)
                if partition is not None:
                    _store_metator_entries(staged, {"partition": partition})
            if mge is not None:
                _store_round(staged, metacool, *mge)
        logger.info(f"{outfile} is written: the contact map{written}.")
    _warn_dropped_rounds(later, metacool, outfile, in_place)


def _what_is_written(contigs, network, rounds, clusters, mge):
    """What the messages of write_map say of the map it writes: while it is
    written, and once it is (each a text that follows the name of the map,
    and the words "the contact map")."""
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
        return (
            f", with the {len(tables[0])} clusters of the round {name} of its clustering rounds",
            f", with the round {name} of its clustering rounds",
        )
    if rounds is not None:
        return f", with the {len(rounds[0])} clusters of the run", ", with the clustering rounds of the run"
    if network is not None:
        return ", with the network of its contigs", ", with the network of its contigs"
    return "", ""


def _store_final_clusters(path, metacool, clusters, partition=None):
    """Store the final clusters of a validation (clusters, a FinalClusters)
    in the map at path, the staged map of an edit session (see write_map),
    whose metator metadata keeps the entries of the steps before it
    (network, partition: that of partition, for a pipeline run), and gets
    the entry validation. metacool names the map in the warnings."""
    import metacooler

    stored = open_map(path)
    names = [str(name) for name in stored.contigs["name"]]
    metator = {key: value for key, value in _metator_metadata(stored).items() if key not in LEGACY_RUN_KEYS}
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
        metacooler.reassign(
            path,
            assignment,
            # Grouped by cluster, then from the longest, whatever the order of
            # the input map: the pixels are remapped, the columns of the
            # contigs table and the weights of the bins follow them.
            layout="genome-then-length",
            # The cluster names are the first column of the table, unnamed in
            # cluster_summary.txt. It replaces the columns of the genomes
            # table that the input map carried.
            genomes_table=genomes_table,
            genomes_table_format="1",
            metadata={"metator": metator},
            # The staged copy of the map, which metacooler writes anew when it
            # regroups its contigs: its history names the map of the run.
            inplace=True,
            # On a failure, metacooler removes the partial file it writes
            # next to the staged map.
            keep_partial=False,
        )
    store_rounds(path, clusters.rounds, name=metacool)


def _store_round(path, metacool, rounds, entry):
    """Replace one clustering round of the map at path, the staged map of an
    edit session (see write_map), and keep everything else as it was: the
    other rounds, its bins, contacts, genomes and contigs tables, QC block
    and metadata. metator mge stores its MGE clusters so, as the round mge,
    next to the rounds of the binning run, whose final clusters the genomes
    and contigs tables keep.

    The round replaces the round of its name that the map held, if any, and
    comes after the others, as the last made. The other rounds are kept as
    the map held them (see metacooler's MetaCooler.clusters): the entry of
    each, its clusters, their parents, contigs and the values of their
    columns, but for those that metacooler computes (n_contigs, length and
    genome), which it computes again. In the one table of the clusters of
    every round, a column that only one of them has is missing in the other:
    an empty text, or no number. So a column of whole numbers, or of
    booleans, that the round lacks is stored as decimals, or as text (True,
    False), as any column with a missing value is.

    Parameters
    ----------
    path : str
        The staged map.
    metacool : str
        The map that path is staged for, which the warnings name.
    rounds : tuple of pandas.DataFrame
        The clusters and the contigs of the round, as
        metator.clusters.Rounds.tables gives them: every cluster of the
        round named by entry.
    entry : dict
        The round as metacooler stores it: its name, and what is stored with
        it (tool, version, step, parameters...).
    """
    from metacooler.schema import CLUSTER_GENOME_COLUMN

    clusters, members = rounds
    name = entry["name"]
    mc = open_map(path)
    held = mc.clusters
    held_members = mc.cluster_members()
    held_rounds = mc.cluster_rounds
    # The columns that metacooler computes, which it refuses as given.
    computed = ["n_contigs", "length", CLUSTER_GENOME_COLUMN]
    other = held[(held["round"] != name).to_numpy()].drop(columns=computed)
    other_members = held_members[(held_members["round"] != name).to_numpy()]
    entries = [stored for stored in held_rounds if stored["name"] != name]
    _store_clusters(path, [*entries, entry], _stacked(other, clusters), _stacked(other_members, members), name=metacool)
