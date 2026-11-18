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
its bins and contacts, grouped by final bin in the same way.

Core functions to write the file:
    - write_metacool, from the pairs files
    - reassign_metacool, from a .metacool map
"""

import gzip
import logging
import os
import tempfile
from contextlib import contextmanager

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
        metaTOR measured nothing there (Shotgun_coverage without a depth
        file, Restriction_site without an enzyme).
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


def map_file(metacool):
    """The file of a .metacool map: the map, or the .mmetacool collection of
    an entry given as file.mmetacool::entries/NAME."""
    return str(metacool).split("::", 1)[0]


def is_map(path):
    """Whether path is a contact map rather than a pairs file: the URI of
    one inside a file (file::group), or an HDF5 file, which a pairs file
    never is. A .metacool map is one; so are a .mmetacool collection and a
    .cool file, which open_map then refuses with metacooler's reason."""
    if "::" in str(path):
        return True
    try:
        with open(path, "rb") as handle:
            return handle.read(len(HDF5_SIGNATURE)) == HDF5_SIGNATURE
    except OSError:
        return False


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


def _checked_inputs(contig_data, outfile, force, inputs):
    """What both writers check before any work, in this order: the outfile
    (replaced only with force), the input files (None for one not given),
    then the contig data and its names. Returns what read_contig_data
    reads."""
    if os.path.exists(outfile) and not force:
        raise MetacoolError(f"{outfile} already exists. Use --force to overwrite it.")
    missing = [path for path in inputs if path and not os.path.isfile(path)]
    if missing:
        raise MetacoolError(f"Input file(s) not found: {', '.join(map(str, missing))}.")
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
        except BaseException:
            # No map without the columns of its contigs.
            os.remove(outfile)
            raise


def reassign_metacool(metacool, contig_data, outfile, bin_summary=None, force=False, metadata=None):
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
        contigs = open_map(metacool).contigs
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
        f"{sum(value != 'ND' for value in assignment.values())} of them in {len(final_bins)} bins."
    )
    with _metacooler_logging():
        with _as_metacool_errors():
            directory = os.path.dirname(outfile) or "."
            os.makedirs(directory, exist_ok=True)
            handle, staged = tempfile.mkstemp(prefix=f".{os.path.basename(outfile)}.", suffix=".tmp", dir=directory)
            os.close(handle)
        try:
            with _as_metacool_errors(), tempfile.TemporaryDirectory() as scratch:
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
                metacooler.validate(staged, full=True)
                os.replace(staged, outfile)
            logger.info(f"{outfile} is written: the contact map of the bins.")
        except BaseException:
            if os.path.lexists(staged):
                os.remove(staged)
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
    CONTIG_COLUMNS that this run did not measure (Shotgun_coverage without
    a depth file...).
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
