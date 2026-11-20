#!/usr/bin/env python3
# coding: utf-8

"""Contact map of a metaTOR run as a .metacool file.

A .metacool file is the HDF5 contact map of metacooler: a cooler whose
chromosomes are the contigs, grouped by genome. It stays readable by
``cooler.Cooler``, and ``metacooler.MetaCooler`` adds the genome layer.

metaTOR writes one bin per contig, and one genome per final bin: the
``Final_bin`` of the contig data file. The unbinned contigs, whose
``Final_bin`` is ``ND``, are kept with all their contacts, in metacooler's
``__unbinned__`` genome: a run with no final bin gives a map whose contigs
are all there. The columns of the bin summary (completeness,
contamination, abundances...) are stored in the genomes table, and what
metaTOR measured on each contig (GC content and Hi-C hits, plus shotgun
coverage and restriction sites when it measured them) in the contigs table.

Each pairs file row whose two ends lie on contigs of the contig data file
(within their length) is one contact between the two contigs, an
intra-contig contact when both ends are on one contig: the map holds the
raw contacts, without the filters and normalization of the network. The
pairs files of several libraries are summed into one map.

A run that started from a .metacool map of the contigs (made before any
binning) has no pairs files: its map is written from that map instead, with
its bins and contacts, grouped by final bin in the same way. A pipeline run
from reads or pairs builds such a map itself, from its pairs files, as soon
as it has them, reads its network from it (with several libraries, from one
map per library and the map of their sum, see write_library_maps), and
writes the map of its bins from it at the end. A run step by step from a
map stores the overlapping clusters of its partition in the map first, as a
clustering round, then its validation writes the map of the clusters.

Core functions to write the file:
    - write_metacool, from the pairs files
    - write_contig_map, the map of the contigs in no bin yet, from the pairs
      files
    - reassign_metacool, from a .metacool map
    - write_rounds, a .metacool map with the clustering rounds of a run
"""

import gzip
import json
import logging
import os
import shutil
import socket
import sys
import tempfile
import time
import uuid
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

from metator.log import logger
from metator.version import __version__

# Final_bin of a contig that is in no final bin: "ND" as written by metaTOR,
# "-" or empty when a table was edited by hand.
UNBINNED = {"ND", "-", ""}

# The genome names that metacooler reads as no genome, stripped and whatever
# their case (_UNBINNED_SENTINELS in metacooler/io/create.py; -1 from
# metacooler 0.2.0): it puts the contigs of a bin so named in __unbinned__,
# without an error. metaTOR's own bin names (metator_*, MetaTOR_MGE_*) are
# never one of them.
METACOOLER_UNBINNED = frozenset({"__unbinned__", "na", "nan", "n/a", "nd", "n/d", "-1", ""})

# The columns of contig_data_final.txt stored in the contigs table: what
# metaTOR measured on each contig. Never a bin column (Final_bin, the bin IDs
# of the steps): the genome of a contig is its genome_id, and metacooler
# reassign carries the contigs columns into the map of another binning, where
# a bin column would be stale.
CONTIG_COLUMNS = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]

# Those that metaTOR measures on the assembly, by their key in the contig data
# of metator.network.create_contig_data: all but Hit, which the network step
# counts on the pairs.
ASSEMBLY_COLUMNS = {"GC_content": "GC", "Shotgun_coverage": "coverage", "Restriction_site": "RS"}

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
    its message kept: metator metacool logs the reason, then raises the
    error."""
    try:
        yield
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
        run: its Name, Size and Final_bin columns are read, and those of
        CONTIG_COLUMNS that it has.

    Returns
    -------
    pandas.DataFrame
        One bin per contig, from 0 to its size: columns chrom, start and
        end.
    dict
        Final bin of each binned contig, by contig name. The unbinned
        contigs (Final_bin ND, - or empty) are left out.
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
        bin), or it has no Name, Size or Final_bin column, or a Size that is
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
    missing = [column for column in ["Name", "Size", "Final_bin"] if column not in table.columns]
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
    final_bin = table["Final_bin"].str.strip()
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


def _contigs_in_map(mc, metacool, contigs, of, what):
    """The contigs of contigs (their lengths, by name, as the contig data
    file of gives them) that the map mc (that of metacool, a
    metacooler.MetaCooler) holds, in the order of the map, once the map is
    checked for reading their contacts.

    Refused: a map that stores both triangles of the matrix (storage-mode
    square), where each contact between two bins would be read twice, and
    a map that gives one of the contigs another length (see
    check_map_lengths), which is no map of them. The contigs that the map
    lacks (a map made with a contig filter, for example) are left out, with
    a warning (what: what the contigs are, for it): they have no contacts
    there, as contigs that no read maps to."""
    storage_mode = mc.info.get("storage-mode", "symmetric-upper")
    if storage_mode != "symmetric-upper":
        raise MetacoolError(
            f"{metacool} stores its contacts with storage-mode {storage_mode!r}: the contacts of contigs are read "
            "from a map that stores each contact once, in the upper triangle (storage-mode 'symmetric-upper'), as "
            "metacooler writes them."
        )
    check_map_lengths(metacool, mc.contigs, contigs, of, of)
    names = [str(name) for name in mc.contigs["name"]]
    in_map = set(names)
    absent = [name for name in contigs if name not in in_map]
    if absent:
        logger.warning(
            f"{len(absent)} of {what} are not in {metacool}: they have no contacts there, as contigs that no read maps "
            f"to: {_shown(absent)}."
        )
    return [name for name in names if name in contigs]


def contig_contacts(metacool, contigs, of, what, chunksize=1_000_000):
    """The contacts between contigs in a .metacool contact map, summed by
    pair of contigs: those between their bins, whatever the bin size of the
    map. The pixels are read at most chunksize at a time, and only those of
    the contigs are kept: the memory used grows with the number of their
    pairs in contact. metator qc reads the contacts of the large contigs of
    its high-quality MAGs so.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    contigs : dict
        The lengths of the contigs, by name, as the contig data file of
        gives them.
    of : str
        The contig data file, for the messages.
    what : str
        What the contigs are, for the warning about those that the map
        lacks (see _contigs_in_map), which have no contacts.
    chunksize : int
        Number of pixels read at a time.

    Returns
    -------
    pandas.DataFrame
        One row per pair of contigs in contact: contig1 and contig2 (their
        names, the first one first in the map, the same one for the contacts
        within a contig) and count, their contacts.

    Raises
    ------
    MetacoolError
        If the map cannot be read, stores both triangles of the matrix, or
        gives a contig another length than contigs.
    """
    import numpy as np

    with _metacooler_logging():
        mc = open_map(metacool)
        held = _contigs_in_map(mc, metacool, contigs, of, what)
        names = mc.contigs["name"].astype(str).to_numpy()
        selected = np.isin(names, held)
        sums = [pd.DataFrame({"id1": [], "id2": [], "count": []}, dtype="int64")]
        # No pixel to read when the map holds none of the contigs.
        if held:
            for chunk in mc.iter_pixels(chunksize=chunksize, contig_ids=True):
                pixels = pd.DataFrame({"id1": chunk["contig1_id"], "id2": chunk["contig2_id"], "count": chunk["count"]})
                pixels = pixels[selected[pixels["id1"].to_numpy()] & selected[pixels["id2"].to_numpy()]]
                sums.append(pixels.groupby(["id1", "id2"], as_index=False)["count"].sum())
    summed = pd.concat(sums).groupby(["id1", "id2"], as_index=False)["count"].sum()
    return pd.DataFrame(
        {
            "contig1": names[summed["id1"].to_numpy()],
            "contig2": names[summed["id2"].to_numpy()],
            "count": summed["count"].to_numpy(),
        }
    )


def held_contigs(metacool, contigs, of, what):
    """The contigs of contigs that the .metacool contact map metacool holds,
    in the order of the map, once the map is checked for reading their
    contacts (see _contigs_in_map): those that extract_contigs writes the
    contact map of. metator contactmap checks the map so before it writes
    any output.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    contigs : dict
        The lengths of the contigs, by name, as the contig data file of
        gives them.
    of : str
        The contig data file, for the messages.
    what : str
        What the contigs are, for the warning about those that the map
        lacks, which are left out.

    Returns
    -------
    list of str
        The contigs of contigs that the map holds, in its order.

    Raises
    ------
    MetacoolError
        If the map cannot be read, stores both triangles of the matrix, or
        gives a contig another length than contigs.
    """
    with _metacooler_logging():
        return _contigs_in_map(open_map(metacool), metacool, contigs, of, what)


def extract_contigs(metacool, contigs, outfile):
    """Write the contact map of contigs, from the .metacool contact map
    metacool, as metacooler extract writes it: a .cool file at outfile (a
    cooler, not a .metacool map) of their bins in the map, one per contig
    in a map of metaTOR, in the order of the map, with their contacts.
    metator contactmap writes the contact map of a bin so, from a map.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    contigs : list of str
        The contigs, which the map holds (see held_contigs).
    outfile : str
        Path to the .cool file to write, replaced if it exists.

    Raises
    ------
    MetacoolError
        If metacooler cannot write the file (a ValueError, OSError or
        KeyError, message kept).
    """
    import metacooler

    with _metacooler_logging(), _as_metacool_errors():
        metacooler.extract(metacool, contigs=list(contigs), output=outfile)


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
    """Refuse contig or bin names that are not ASCII.

    cooler reads the contig names (its chromosome names) as ASCII text, so
    metacooler refuses any other name.
    """
    for kind, names in [("contig", contigs), ("bin", bins)]:
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
            "A .metacool file cannot have bins named as metacooler marks a contig in no bin (NA, NaN, N/A, N/D, -1 or "
            "__unbinned__ in any case, and nd in any case other than ND, which marks a contig in no bin): it would put "
            f"their contigs in __unbinned__. Rename the {len(bad)} bin(s): {_shown(bad)}."
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
    """What both writers of the map of the bins check before any work, in
    this order: the outfile and the input files (see _check_paths), then
    the contig data and its names. Returns what read_contig_data reads."""
    _check_paths(outfile, force, inputs)
    bins, contig_bin, contig_columns, header = read_contig_data(contig_data)
    check_ascii_names(bins["chrom"], contig_bin.values())
    check_bin_names(contig_bin.values())
    return bins, contig_bin, contig_columns, header


def write_metacool(
    pairs_files,
    contig_data,
    outfile,
    bin_summary=None,
    assembly_name=None,
    threads=1,
    force=False,
    metadata=None,
):
    """Write the contact map of the contigs as a .metacool file.

    One bin per contig, one genome per final bin, plus metacooler's
    ``__unbinned__`` genome for the unbinned contigs, which keep their
    contacts. The contigs are ordered by genome, then from the longest
    (metacooler's genome-then-length ordering).

    Parameters
    ----------
    pairs_files : str or list of str
        Pairs file(s) of the Hi-C libraries, plain or compressed: their
        contacts are summed. Each row whose two ends lie on contigs of the
        contig data file (within their length) counts once; the other rows
        are left out, with a warning.
    contig_data : str
        Path to the contig_data_final.txt file: the contigs (Name, Size)
        and their final bins (Final_bin). Its GC_content, Hit,
        Shotgun_coverage and Restriction_site columns (CONTIG_COLUMNS) are
        stored in the contigs table, those it has and fills (see
        read_contig_data).
    outfile : str
        Path to the .metacool file to write.
    bin_summary : str, optional
        Path to the bin_summary.txt file: its columns are added to the
        genomes table, one row per bin, named by its first column. Each row
        must have as many fields as the header. Not read when no contig is
        in a final bin: the map has no genome but ``__unbinned__``.
    assembly_name : str, optional
        Name of the assembly, stored in the file.
    threads : int
        Number of threads binning the pairs.
    force : bool
        Whether to replace an existing outfile.
    metadata : dict, optional
        Parameters of the run, stored with the version of metaTOR under the
        "metator" key of the file metadata.

    Raises
    ------
    MetacoolError
        If threads is below 1, the outfile exists (without force), an input
        file is missing, the contig data cannot be read as a table or has
        no Name, Size or Final_bin column or a Size that is not a whole
        number, a contig or bin name is not ASCII, a bin name is one that
        metacooler reads as no bin (NA, nan, -1...), or metacooler cannot
        write the map or store the contig columns (a ValueError, OSError or
        KeyError, message kept: for example, a bin summary whose rows and
        header have different numbers of fields, or an outfile in a
        directory that cannot be made). No map is left when the contig
        columns cannot be stored, and no partial file of metacooler when it
        cannot write the map.
    """
    pairs_files = [str(pairs_files)] if isinstance(pairs_files, (str, os.PathLike)) else [str(p) for p in pairs_files]
    if threads < 1:
        raise MetacoolError(f"The number of threads must be at least 1, got {threads}.")
    bins, contig_bin, contig_columns, header = _checked_inputs(
        contig_data, outfile, force, [*pairs_files, contig_data, bin_summary]
    )

    import metacooler

    logger.info(
        f"Writing {outfile}: contact map of {len(bins)} contigs, {len(contig_bin)} of them in "
        f"{len(set(contig_bin.values()))} bins."
    )
    with _metacooler_logging(), _as_metacool_errors():
        if os.path.dirname(outfile):
            os.makedirs(os.path.dirname(outfile), exist_ok=True)
        metacooler.from_pairs(
            outfile,
            bins=bins,
            pairs=pairs_files,
            # Every contig, the unbinned ones as ND (no genome for
            # metacooler): a run with no final bin gives a table all ND,
            # not an empty table, which metacooler 0.2.0 refuses.
            contig_genome_map={name: contig_bin.get(name, "ND") for name in bins["chrom"]},
            assembly_name=assembly_name,
            # The bin names are the first column of bin_summary.txt, unnamed.
            # No final bin, no bin summary row: its columns would only be
            # empty columns of __unbinned__.
            genomes_table=bin_summary if bin_summary and contig_bin else None,
            genomes_table_format="1",
            metadata={"metator": {"version": __version__, **(metadata or {})}},
            # No P(s) table next to the file: metaTOR writes the map only.
            ps=False,
            nproc=threads,
            # On a failure, metacooler removes the partial file it writes
            # next to the outfile.
            keep_partial=False,
        )
        try:
            _write_back_columns(outfile, contig_data, contig_columns, header)
            _rewrite_in_a_fresh_file(outfile)
        except BaseException:
            # No map without the columns of its contigs.
            os.remove(outfile)
            raise


def write_contig_map(
    pairs_files,
    contig_data,
    outfile,
    assembly_name=None,
    threads=1,
    force=False,
    metadata=None,
    tmp_dir=None,
):
    """Write the contact map of the contigs, in no bin yet, as a .metacool
    file: the map that a pipeline run builds from its pairs files as soon as
    it has them. A run with one library reads its network from it
    (metator.network.metacool_to_contacts: the network of its pairs files,
    byte for byte). At the end of the run, reassign_metacool
    writes the map of its bins from it, which is the map that
    write_metacool writes from the same pairs files and the final tables,
    but for its history (and, when the final bins leave the contigs in the
    order of this map, for example with no final bin or one bin holding
    every contig, the place of the genome_id column: metacooler then
    relabels the contigs where they are, and writes genome_id again, after
    the QC columns).

    One bin per contig, every contig in metacooler's ``__unbinned__``
    genome, from the longest (metacooler's genome-then-length ordering,
    with one genome), and the contacts of the pairs files as write_metacool
    counts them, with the contact statistics that metacooler measures on
    the reads (the QC block, and the ``qc_*`` columns of the contigs, which
    reassign_metacool carries over to the bins). The contigs table holds
    what metaTOR measured on the assembly (ASSEMBLY_COLUMNS): GC_content,
    and Shotgun_coverage and Restriction_site when they were measured (with
    a depth file, with an enzyme). Not Hit, which the network step counts:
    the map of the bins takes it from the contig data file, with the other
    columns.

    Parameters
    ----------
    pairs_files : str or list of str
        Pairs file(s) of the Hi-C libraries, plain or compressed: their
        contacts are summed. Each row whose two ends lie on contigs of the
        contig data (within their length) counts once; the other rows are
        left out, with a warning.
    contig_data : dict
        The contigs of the assembly, by name, as
        metator.network.create_contig_data returns them: their "length",
        "GC", "coverage" and "RS" ("-" where nothing was measured). The bins
        are given to metacooler in their order, that of the assembly;
        metacooler orders the contigs from the longest. No fasta is given
        to metacooler, which would write the index of the assembly (.fai)
        next to it, and no P(s) table is written next to the map.
    outfile : str
        Path to the .metacool file to write.
    assembly_name : str, optional
        Name of the assembly, stored in the file: the map of the bins keeps
        it.
    threads : int
        Number of threads binning the pairs.
    force : bool
        Whether to replace an existing outfile.
    metadata : dict, optional
        Parameters of the run, stored with the version of metaTOR under the
        "metator" key of the file metadata.
    tmp_dir : str, optional
        Directory of the table of the contig columns, which metacooler
        annotate reads: the temporary directory of the system by default.

    Raises
    ------
    MetacoolError
        If threads is below 1, the outfile exists (without force), a pairs
        file is missing, a contig name is not ASCII, or metacooler cannot
        write the map or store the contig columns (a ValueError, OSError or
        KeyError, message kept). No map is left when the contig columns
        cannot be stored, and no partial file of metacooler when it cannot
        write the map.
    """
    pairs_files = [str(pairs_files)] if isinstance(pairs_files, (str, os.PathLike)) else [str(p) for p in pairs_files]
    if threads < 1:
        raise MetacoolError(f"The number of threads must be at least 1, got {threads}.")
    _check_paths(outfile, force, pairs_files)
    names = list(contig_data)
    check_ascii_names(names, [])

    import metacooler

    logger.info(f"Writing {outfile}: contact map of {len(names)} contigs, in no bin yet.")
    with _metacooler_logging(), _as_metacool_errors():
        if os.path.dirname(outfile):
            os.makedirs(os.path.dirname(outfile), exist_ok=True)
        metacooler.from_pairs(
            outfile,
            bins=pd.DataFrame({"chrom": names, "start": 0, "end": [contig_data[name]["length"] for name in names]}),
            pairs=pairs_files,
            # No binning yet: every contig in __unbinned__, the longest first.
            # The order of the assembly (contig_ordering="fasta") would need
            # the fasta, or #chromsize lines that pairs files may lack; the
            # map of the bins is grouped by bin, whatever the order of this
            # map.
            contig_genome_map=None,
            assembly_name=assembly_name,
            metadata={"metator": {"version": __version__, **(metadata or {})}},
            # The QC block, measured on the reads, which only the pairs give.
            qc=True,
            # No P(s) table next to the file: metaTOR writes the map only.
            ps=False,
            nproc=threads,
            # On a failure, metacooler removes the partial file it writes
            # next to the outfile.
            keep_partial=False,
        )
        try:
            _store_assembly_columns(outfile, contig_data, tmp_dir)
        except BaseException:
            # No map without the columns of its contigs.
            os.remove(outfile)
            raise
    logger.info(f"{outfile} is written: the contact map of the contigs, in no bin yet.")


def library_entry(collection, i):
    """The URI of the map of library i (from 0) in the collection of the
    maps of the libraries of a run (see write_library_maps)."""
    return f"{collection}::entries/library_{i}"


def write_library_maps(
    pairs_files,
    contig_data,
    outfile,
    collection,
    assembly_name=None,
    threads=1,
    force=False,
    metadata=None,
    tmp_dir=None,
):
    """Write the contact map of each library of a run, as the entries of a
    .mmetacool collection, and the map of their sum, merged from them: the
    maps that a pipeline run with several libraries builds from its pairs
    files as soon as it has them.

    The map of library i (library_entry(collection, i), from 0, in the
    order of pairs_files) is the map that write_contig_map writes from its
    pairs file alone, with the contact statistics measured on its reads,
    but for the columns measured on the assembly, which the map of the sum
    holds. The map of the sum (outfile) is merged from them with metacooler
    merge, which sums their contacts and their contact statistics: the map
    that write_contig_map writes from all the pairs files at once, but for
    its history. A run reads its network from the sum, and the network and
    the hits of each library from its map (see
    metator.network.metacool_to_contacts). At the end of the run, the map of
    the clusters is written over the sum.

    Parameters
    ----------
    pairs_files : list of str
        Pairs file of each Hi-C library.
    contig_data, assembly_name, threads, metadata, tmp_dir
        As for write_contig_map.
    outfile : str
        Path to the .metacool file of the sum.
    collection : str
        Path to the .mmetacool file of the libraries.
    force : bool
        Whether to replace an existing outfile or collection.

    Raises
    ------
    MetacoolError
        As write_contig_map, for either file. When a map cannot be written,
        neither file is left.
    """
    pairs_files = [str(path) for path in pairs_files]
    if threads < 1:
        raise MetacoolError(f"The number of threads must be at least 1, got {threads}.")
    _check_paths(outfile, force, pairs_files)
    _check_paths(collection, force, [])
    names = list(contig_data)
    check_ascii_names(names, [])

    import metacooler

    bins = pd.DataFrame({"chrom": names, "start": 0, "end": [contig_data[name]["length"] for name in names]})
    os.makedirs(os.path.dirname(os.path.abspath(collection)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
    logger.info(f"Writing {collection}: contact map of each of the {len(pairs_files)} libraries.")
    run = {"metator": {"version": __version__, **(metadata or {})}}

    def contig_map(path, pairs):
        # As write_contig_map: no binning yet, the longest first.
        metacooler.from_pairs(path, bins=bins, pairs=pairs, contig_genome_map=None, assembly_name=assembly_name,
                              metadata=run, qc=True, ps=False, nproc=threads, keep_partial=False)  # fmt: skip

    try:
        # The map of each library is written in the run's temporary
        # directory, then copied into the collection.
        with _metacooler_logging(), _as_metacool_errors(), tempfile.TemporaryDirectory(dir=tmp_dir) as scratch:
            maps, contacts = {}, []
            for i, pairs_file in enumerate(pairs_files):
                logger.info(f"Library {i} ({library_entry(collection, i)}): {pairs_file}.")
                path = os.path.join(scratch, f"library_{i}.metacool")
                with _named_library(i, pairs_file):
                    contig_map(path, [pairs_file])
                maps[f"library_{i}"] = Path(path)
                contacts.append(int(metacooler.MetaCooler(path).info["nnz"]))
            metacooler.create_mmetacool(Path(collection), maps, force=True)
            logger.info(f"Writing {outfile}: contact map of the {len(pairs_files)} libraries, summed.")
            if os.path.lexists(outfile):
                os.remove(outfile)
            if any(contacts):
                metacooler.merge([library_entry(collection, i) for i in range(len(pairs_files))], outfile)
            else:
                # No library holds a contact, and cooler merges no map without
                # one: the map of the sum is the map of all the pairs files at
                # once, as write_contig_map writes it, with no contact either.
                contig_map(outfile, pairs_files)
            # The run's metadata, which metacooler merge does not carry over.
            _store_assembly_columns(outfile, contig_data, tmp_dir, metadata=run)
    except BaseException:
        # Both or neither.
        for path in (outfile, collection):
            if os.path.lexists(path):
                os.remove(path)
        raise
    logger.info(f"{outfile} is written: the contact map of the contigs, in no cluster yet.")


def _store_assembly_columns(path, contig_data, tmp_dir=None, metadata=None):
    """Give the map of the contigs the columns that metaTOR measured on the
    assembly (ASSEMBLY_COLUMNS), as metacooler annotate adds them: those
    measured on a contig at least (read_contig_data leaves out a column that
    is ``-`` on every contig), GC_content always, with the values that the
    contig data files of the run hold (write_contig_data writes their
    str). metadata: keys to store in the metadata of the map, as metacooler
    annotate merges them."""
    import metacooler

    columns = [column for column, key in ASSEMBLY_COLUMNS.items() if any(str(data[key]) != "-" for data in contig_data.values())]
    with tempfile.TemporaryDirectory(dir=tmp_dir) as scratch:
        table = os.path.join(scratch, "contig_columns.txt")
        with open(table, "w") as handle:
            handle.write("\t".join(["Name", *columns]) + "\n")
            for name, data in contig_data.items():
                handle.write("\t".join([name, *(str(data[ASSEMBLY_COLUMNS[column]]) for column in columns)]) + "\n")
        metacooler.annotate(
            path,
            contigs_table=table,
            contigs_table_format="Name",
            # metacooler leaves out a contig of no length, which the
            # assembly may hold: no warning about its row.
            unmatched="ignore",
            inplace=True,
            metadata=metadata,
        )


@contextmanager
def _named_library(i, pairs_file):
    """Name the library and its pairs file in an error that its map raises:
    each library is its own map, written from a temporary file."""
    try:
        yield
    except (ValueError, OSError, KeyError) as error:
        reason = error.args[0] if isinstance(error, KeyError) and len(error.args) == 1 else error
        raise MetacoolError(f"Library {i}, {pairs_file}: {reason}") from error


def reassign_metacool(metacool, contig_data, outfile, bin_summary=None, force=False, metadata=None, rounds=None):
    """Write the contact map of the final bins from a .metacool map of the
    contigs, without pairs files.

    The bins and contacts are those of the input map, at its bin size (one
    bin per contig, or a fixed size), made before any binning (by
    metacooler create pairs without --contigs, for example). Each contig
    goes to its final bin, the unbinned ones and those that the contig data
    lacks to ``__unbinned__``, and the contigs are grouped by genome, then
    from the longest, as write_metacool groups them (metacooler reassign,
    layout genome-then-length): a map of the same genomes, in the same
    order. The genomes table holds the columns of the bin summary, and the
    contigs table the columns of the contig data that write_metacool stores
    (see read_contig_data), last, as write_metacool stores them: they
    replace those that the input map carried, metaTOR's own as those of the
    bins of an earlier run. Its other columns of the contigs table are
    kept. The QC block of the input map, when it has one, is carried over
    to the bins: metacooler keeps the counts within each contig (measured
    on the reads, for a map made from pairs files) and counts the contacts
    between contigs again. A map written by metacooler 0.2.0, with no count
    per contig, has its block computed again from the pixels (at one bin
    per contig, every contact within a contig then counts as short-range),
    which metacooler warns about for a map made from pairs files. The
    assembly name and the history of the commands are those of the input
    map.

    The map is written to a temporary file next to outfile, then given
    metacooler's full validation, and only then moved to outfile, which it
    replaces: a symbolic link at outfile is replaced, not its target.
    outfile may thus be the input map itself (the same path, or a symbolic
    link to it or from it): outfile is only replaced by a valid map of the
    bins, and an input map that outfile links to is left as it was. On any
    failure, outfile is left as it was, and the temporary files are
    removed, metacooler's partial file among them.

    Parameters
    ----------
    metacool : str
        Path to the input .metacool file, or URI of an entry of a
        .mmetacool collection: ``file.mmetacool::entries/<name>``.
    contig_data : str
        Path to the contig_data_final.txt file: the final bin of the
        contigs (Name, Final_bin), their size (Size, which must be their
        length in the map), and what metaTOR measured on them.
    outfile : str
        Path to the .metacool file to write.
    bin_summary : str, optional
        Path to the bin_summary.txt file: its columns are added to the
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

    Raises
    ------
    MetacoolError
        If the outfile exists (without force) or is the .mmetacool
        collection of the input map, an input file is missing, the contig
        data cannot be read as a table or has no Name, Size or Final_bin
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
            f"{outfile} is the .mmetacool collection of the input map {metacool}: the map of the bins would replace "
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
            f"{len(absent)} binned contig(s) of {contig_data} are not in {metacool}, and not in the map of the bins "
            f"either: {_shown(absent)}."
        )
    # Every contig of the map, the unbinned ones as ND (no genome for
    # metacooler): never an empty table, which metacooler refuses.
    assignment = {name: contig_bin.get(name, "ND") for name in names}
    final_bins = set(assignment.values()) - {"ND"}
    logger.info(
        f"Writing {outfile} from the contact map {metacool}: {len(names)} contigs, "
        f"{sum(value != 'ND' for value in assignment.values())} of them in {len(final_bins)} clusters."
    )
    with _metacooler_logging():
        with _staged_map(outfile) as staged, tempfile.TemporaryDirectory() as scratch:
            if bin_summary and final_bins:
                genomes_table = bin_summary
            else:
                # Without a bin summary or without a bin, the bin names
                # alone: the genome columns that the input map carried
                # (those of an earlier run's bins) are then replaced by
                # none, where metacooler would keep them and warn that
                # they no longer describe their bins.
                genomes_table = os.path.join(scratch, "bins.txt")
                with open(genomes_table, "w") as handle:
                    handle.write("".join(f"{name}\n" for name in ["bin", *sorted(final_bins)]))
            metacooler.reassign(
                metacool,
                assignment,
                staged,
                # Grouped as write_metacool groups the contigs, whatever
                # the order of the input map: the pixels are remapped.
                layout="genome-then-length",
                # The bin names are the first column of the table,
                # unnamed in bin_summary.txt. It replaces the columns of
                # the genomes table that the input map carried.
                genomes_table=genomes_table,
                genomes_table_format="1",
                metadata={"metator": {"version": __version__, **(metadata or {})}},
                # Over the empty temporary file.
                force=True,
                # On a failure, metacooler removes the partial file it
                # writes next to the temporary file.
                keep_partial=False,
            )
            _write_back_columns(staged, contig_data, contig_columns, header)
            store_rounds(staged, rounds, name=metacool)
            _rewrite_in_a_fresh_file(staged)
        logger.info(f"{outfile} is written: the contact map of the clusters.")
    _warn_dropped_rounds(later, metacool, outfile, in_place)


@contextmanager
def _staged_map(outfile):
    """A map written next to outfile, then moved over it once it is valid.

    Yields the path of an empty temporary file next to outfile, where the
    map is written. Once it is, the map is given metacooler's full
    validation, and only then moved to outfile, which it replaces
    (os.replace): a symbolic link at outfile is replaced, not its target,
    so outfile may be the input map of the run (the same path, or a
    symbolic link to it or from it). The map keeps the mode of the file it
    replaces (that of the target of a symbolic link): a private map stays
    private. A new map gets the mode of any new file, from the umask, not
    the private mode of the temporary file (mkstemp's 0600). On any
    failure, Ctrl-C included, outfile is left as it was,
    and the temporary file is removed. A ValueError, OSError or KeyError
    (metacooler's above all) is raised as a MetacoolError, its message
    kept.
    """
    import metacooler

    with _as_metacool_errors():
        directory = os.path.dirname(outfile) or "."
        os.makedirs(directory, exist_ok=True)
        handle, staged = tempfile.mkstemp(prefix=f".{os.path.basename(outfile)}.", suffix=".tmp", dir=directory)
        os.close(handle)
    try:
        with _as_metacool_errors():
            yield staged
            metacooler.validate(staged, full=True)
            if os.path.exists(outfile):
                shutil.copymode(outfile, staged)
            else:
                os.chmod(staged, 0o666 & ~_umask())
            os.replace(staged, outfile)
    except BaseException:
        if os.path.lexists(staged):
            os.remove(staged)
        raise


def _umask():
    """The umask of the process, which can only be read by setting it: it
    is set back at once."""
    mask = os.umask(0o022)
    os.umask(mask)
    return mask


def _map_group(metacool):
    """The file of a map and the group of the map in it: "" for a map that
    is its file (given as file, file::/ or file::), the group of an entry
    of a .mmetacool collection (entries/NAME, from file::entries/NAME)."""
    path, _, group = str(metacool).partition("::")
    return path, group.strip().strip("/")


def _copy_map(metacool, target, skip=()):
    """Copy the map metacool, as stored, to the root of target, a new HDF5
    file: every object of its group (the root of its file, or an entry of
    a .mmetacool collection, which becomes a map of its own), with its
    dtype, chunks, filters and attributes, the columns of a table in their
    order, and the soft link /chroms as a link, then the attributes of the
    group, with their dtypes; but the objects of the group named in skip."""
    import h5py

    path, group = _map_group(metacool)
    with h5py.File(path, "r") as source, h5py.File(target, "w") as copy:
        node = source[group] if group else source
        for key in node:
            if key in skip:
                continue
            link = node.get(key, getlink=True)
            if isinstance(link, h5py.SoftLink):
                copy[key] = h5py.SoftLink(link.path)
            else:
                source.copy(node[key], copy, name=key)
        for key in node.attrs:
            copy.attrs.create(key, node.attrs[key], dtype=node.attrs.get_id(key).dtype)


def _rewrite_in_a_fresh_file(path):
    """Copy the map at path into a fresh file, which takes its place.

    metacooler reassign copies a map it relabels as it is, and metacooler
    annotate and the clustering rounds edit it in place: each leaves the
    space of what it replaces in the file, which HDF5 does not reclaim. A
    map updated run after run grew with every run (from 129 kB to 272 kB
    over four updates of the map of the synthetic test data). Every object
    is copied as stored (see _copy_map); the file keeps its mode.
    """
    directory = os.path.dirname(os.path.abspath(path))
    handle, fresh = tempfile.mkstemp(prefix=f"{os.path.basename(path)}.", suffix=".fresh", dir=directory)
    os.close(handle)
    try:
        _copy_map(path, fresh)
        shutil.copymode(path, fresh)
        os.replace(fresh, path)
    except BaseException:
        if os.path.lexists(fresh):
            os.remove(fresh)
        raise


def _write_back_columns(path, contig_data, contig_columns, header):
    """Give a map the metadata columns of this run, as metacooler annotate
    edits them: the map that write_metacool has just written, or the map
    written back from an input map.

    The columns of the contig data that metaTOR measured (contig_columns)
    are stored last, in their order, as write_metacool stores them in the
    map it writes from the pairs files: the columns of a contig data file
    (header: all of them) that the map carried are dropped first, then
    those of this run are added. metacooler reassign stores the columns
    that it carries from its input map before genome_id or before the QC
    columns of the contigs: replaced where they were, they would list the
    columns of the map of the bins in another order than in the map of the
    pairs. The dropped columns are also those that this run does not
    store: the bins of each step, the contig ids, Size..., which a map made
    with the contig data of an earlier run as its contigs table holds
    (metacooler create pairs --contigs-format metator), and one of
    CONTIG_COLUMNS that this run had no values for (Shotgun_coverage
    without a depth file, nor a column of the input map that gives every
    contig its coverage...).
    """
    import metacooler

    stored = metacooler.MetaCooler(path)
    # The columns of a contig data file but Name (the contigs' own name
    # column) and those of the map that metacooler derives.
    of_metator = (set(header) | set(CONTIG_COLUMNS)) - {"Name"}
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


def write_rounds(metacool, outfile, rounds, names=None, force=False):
    """Write the map metacool at outfile with the clustering rounds of a run
    in its /clusters group (see store_rounds), in place of those it held,
    and everything else as it was: its bins, contacts, genomes and contigs
    tables, QC block and metadata. A partition run from a map stores its
    overlapping clusters so, before its validation writes the map of the
    final clusters (see reassign_metacool).

    The map is copied as stored to a temporary file next to outfile, but
    for its rounds, which are then stored there: the copy is given
    metacooler's full validation and only then moved to outfile, which it
    replaces, as reassign_metacool writes a map (see _staged_map). outfile
    may thus be the input map itself (the same path, or a symbolic link to
    it or from it): a symbolic link at outfile is replaced, not its target,
    and the map keeps the mode of the file it replaces. On any failure,
    outfile is left as it was, and the temporary file is removed. An entry
    of a .mmetacool collection is written as a map of its own.

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

    Raises
    ------
    MetacoolError
        If the outfile is no file's path (a URI), exists (without force),
        or is the .mmetacool collection of the input map, the input map is
        missing or cannot be read, or the rounds cannot be stored (a
        ValueError, OSError or KeyError of metacooler, message kept).
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
    )
    _warn_dropped_rounds(later, metacool, outfile, in_place)


def write_round(metacool, outfile, rounds, entry, force=False):
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
    written as write_rounds writes it: copied as stored next to outfile,
    but for its rounds, which are then stored there, validated in full,
    then moved over outfile (see _staged_map), which may be the input map
    itself.

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

    Raises
    ------
    MetacoolError
        If the outfile is no file's path (a URI), exists (without force),
        or is the .mmetacool collection of the input map, the input map is
        missing or cannot be read, or the rounds cannot be stored (a
        ValueError, OSError or KeyError of metacooler, message kept).
    """
    from metacooler.schema import CLUSTER_GENOME_COLUMN

    clusters, members = rounds
    name = entry["name"]

    def store(staged, mc):
        held = mc.clusters
        # The columns that metacooler computes, which it refuses as given.
        computed = ["n_contigs", "length", CLUSTER_GENOME_COLUMN]
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
    )


def _stacked(first, second):
    """The rows of two tables, one after the other, with the columns of
    both: those of first, then the others of second."""
    return second.reset_index(drop=True) if first.empty else pd.concat([first, second], ignore_index=True)


def _write_with_rounds(metacool, outfile, force, store, writing, written):
    """What write_rounds and write_round share: the checks of outfile, then
    the copy of the map metacool, as stored but for its /clusters group,
    staged next to outfile, where store(staged, mc) stores the rounds (mc:
    the input map, a metacooler.MetaCooler), validated and moved over
    outfile (see _staged_map). writing and written say what the map holds,
    in the messages of its writing."""
    from metacooler.schema import CLUSTERS_GROUP

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
        with _staged_map(outfile) as staged:
            # The rounds that the map held are not copied: the new ones
            # replace them.
            _copy_map(metacool, staged, skip={CLUSTERS_GROUP})
            store(staged, mc)
        logger.info(f"{outfile} is written: the contact map, {written}.")


def lock_path(path):
    """The lock file of the map at path: <path>.lock, next to it (next to a
    symbolic link, not its target, as the link is what a run replaces). A
    map given as file::/ is the file itself."""
    return f"{map_file(path)}.lock"


def _pid_is_running(pid):
    """Whether the process pid runs on this host."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _create_lock(lock, holder):
    """Create the lock file, with what says who holds it (holder); False when
    it exists already. Any other failure (a directory that is missing or
    not writable) is raised as it is."""
    try:
        handle = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return False
    with os.fdopen(handle, "w") as stream:
        json.dump(holder, stream)
    return True


def _lock_holder(lock):
    """What a lock file says of its holder: {} when it cannot be read, or
    does not hold a JSON object."""
    try:
        with open(lock) as stream:
            holder = json.load(stream)
    except (OSError, ValueError):
        return {}
    return holder if isinstance(holder, dict) else {}


def _take_stale_lock(lock, stale, holder):
    """Replace a lock whose holder is gone (stale, what it says) by this run's
    (holder), unless another run did first: True when this run holds it.

    The stale lock is moved aside first, which only one run can do: a run
    that moved another run's lock instead (taken in the meantime) puts it
    back, and has not taken it. One case is left open: a third run that
    creates the lock while it is moved aside holds it with the run whose
    lock was moved, which cannot be put back. Both then run, the last to
    end replacing the map: each map is whole, as each is staged and
    validated before it replaces the map."""
    aside = f"{lock}.{os.getpid()}.{uuid.uuid4().hex}"
    try:
        os.rename(lock, aside)
    except FileNotFoundError:
        # Another run moved it aside first: the lock is that run's to take.
        return _create_lock(lock, holder)
    if _lock_holder(aside) != stale:
        try:
            os.link(aside, lock)
        except FileExistsError:
            pass
        os.remove(aside)
        return False
    os.remove(aside)
    return _create_lock(lock, holder)


@contextmanager
def map_lock(path, elsewhere="--output-map"):
    """Hold the lock of the map at path for the time of a run that writes it.

    Two runs that write one map, at once, would each replace it with its own
    clusters, the last one winning. The lock, <path>.lock, is created before
    any work, and removed at the end of the run, whatever its end, by the run
    that holds it only. It names that run (its pid, host, start and command,
    and a uuid of its own, which the release compares): a lock whose process
    no longer runs on this host (a run killed with SIGKILL) is removed, with
    a warning, and taken (see _take_stale_lock). The refusals say how to
    write the map elsewhere: with elsewhere, the option of the command.

    Raises
    ------
    MetacoolError
        If another run holds the lock, its holder cannot be told, or the lock
        cannot be created (a directory that is missing or not writable).
    """
    lock = lock_path(path)
    mine = {
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "started": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "command": " ".join(["metator", *sys.argv[1:]]),
        "run": uuid.uuid4().hex,
    }
    try:
        taken = _create_lock(lock, mine)
        if not taken:
            holder = _lock_holder(lock)
            pid, host = holder.get("pid"), holder.get("host")
            if isinstance(pid, int) and host == socket.gethostname() and not _pid_is_running(pid):
                logger.warning(f"Removing {lock}: the run that held it (pid {pid}, since {holder.get('started')}) is gone.")
                taken = _take_stale_lock(lock, holder, mine)
    except OSError as error:
        message = (
            f"{lock} cannot be created ({error.strerror}): a run that writes {path} locks it there, so the directory "
            f"of the map must exist and be writable; write this run's map elsewhere with {elsewhere}."
        )
        logger.error(message)
        raise MetacoolError(message) from error
    if not taken:
        holder = _lock_holder(lock)
        pid, host = holder.get("pid"), holder.get("host")
        who = f"pid {pid} on {host}, since {holder.get('started')}" if holder else "a run that cannot be told"
        message = (
            f"{path} is being written by another metaTOR run ({who}): wait for it to end, write this run's map "
            f"elsewhere with {elsewhere}, or, if no such run is going on, remove {lock}."
        )
        logger.error(message)
        raise MetacoolError(message)
    try:
        yield lock
    finally:
        # This run's lock only: not one that another run took after it was
        # removed by hand.
        if _lock_holder(lock) == mine:
            os.remove(lock)
