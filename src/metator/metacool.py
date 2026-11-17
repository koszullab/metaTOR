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

import bz2
import gzip
import logging
import lzma
import os
import re
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

# The columns of the genomes table that metacooler derives from the contigs,
# with the qc_* columns (_GENOMES_RESERVED_COLUMNS in metacooler/io/create.py):
# the other columns are metadata, as those of the bin summary.
DERIVED_GENOME_COLUMNS = frozenset({"name", "length", "n_contigs", "n_bins"})

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
        If pandas cannot read the file as a table (empty, not UTF-8...), or
        would read it shifted or padded (a row with more or fewer fields
        than the header: the surplus leading fields of the first rows as
        their index, every value then under the name of the column to its
        left, or the fields a row lacks as empty, its contig then in no
        bin), or it has no Name, Size or Final_bin column, or a Size that is
        not a whole number.
    """
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
    that have no name, as metacooler reads such a table. A file compressed
    as gzip, bzip2 or xz is read as pandas reads it; one compressed
    otherwise is not checked."""
    name = str(path)
    if name.endswith((".zip", ".zst", ".tar", ".tgz", ".tar.gz", ".tar.bz2", ".tar.xz")):
        return
    opener = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}.get(os.path.splitext(name)[1], open)
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


def map_file(metacool):
    """The file of a .metacool map: the map, or the .mmetacool collection of
    an entry given as file.mmetacool::entries/NAME."""
    return str(metacool).split("::", 1)[0]


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
        columns cannot be stored.
    """
    pairs_files = [str(pairs_files)] if isinstance(pairs_files, (str, os.PathLike)) else [str(p) for p in pairs_files]
    if threads < 1:
        raise MetacoolError(f"The number of threads must be at least 1, got {threads}.")
    if os.path.exists(outfile) and not force:
        raise MetacoolError(f"{outfile} already exists. Use --force to overwrite it.")
    inputs = [*pairs_files, contig_data] + ([bin_summary] if bin_summary else [])
    missing = [path for path in inputs if not os.path.isfile(path)]
    if missing:
        raise MetacoolError(f"Input file(s) not found: {', '.join(map(str, missing))}.")
    bins, contig_bin, contig_columns, _ = read_contig_data(contig_data)
    check_ascii_names(bins["chrom"], contig_bin.values())
    check_bin_names(contig_bin.values())

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
        )
        if contig_columns:
            try:
                metacooler.annotate(
                    outfile,
                    contigs_table=contig_data,
                    contigs_table_format="Name",
                    # Named one by one: the other columns of the table are
                    # bins, or the unnamed pandas index.
                    columns=[f"contigs.{column}" for column in contig_columns],
                    inplace=True,
                )
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
    (see read_contig_data): they replace those that the input map
    carried, metaTOR's own as those of the bins of an earlier run. Its
    other columns of the contigs table are kept. The QC block is
    recomputed from the pixels when the input map has one. The assembly
    name and the history of the commands are those of the input map.

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
        genomes table, one row per bin, named by its first column. Not read
        when no contig of the map is in a final bin: the map then has no
        genome but ``__unbinned__``, and no column of the genomes table
        that the input map carried.
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
        ValueError, OSError or KeyError, message kept; when metacooler has
        left its partial file, which it logs as left for inspection, the
        message then says that it was removed). An error of another type is
        raised as it is; when metacooler has left its partial file, a
        warning logged before it says that the file was removed.
    """
    metacool, outfile = str(metacool), str(outfile)
    if os.path.exists(outfile) and not force:
        raise MetacoolError(f"{outfile} already exists. Use --force to overwrite it.")
    inputs = [map_file(metacool), contig_data] + ([bin_summary] if bin_summary else [])
    missing = [path for path in inputs if not os.path.isfile(path)]
    if missing:
        raise MetacoolError(f"Input file(s) not found: {', '.join(map(str, missing))}.")
    if is_entry(metacool) and is_input_map(metacool, outfile):
        raise MetacoolError(
            f"{outfile} is the .mmetacool collection of the input map {metacool}: the map of the bins would replace "
            "the whole collection. Write it to another file."
        )
    bins, contig_bin, contig_columns, header = read_contig_data(contig_data)
    check_ascii_names(bins["chrom"], contig_bin.values())
    check_bin_names(contig_bin.values())
    sizes = dict(zip(bins["chrom"], bins["end"], strict=True))

    import metacooler

    with _metacooler_logging():
        try:
            contigs = metacooler.MetaCooler(metacool).contigs
        except (OSError, KeyError, ValueError) as error:
            reason = error.args[0] if isinstance(error, KeyError) and error.args else error
            raise MetacoolError(f"Cannot read the contact map {metacool}: {reason}") from error
    names = [str(name) for name in contigs["name"]]
    # The binning of another version of the assembly (other contigs under
    # the same names) is no binning of this map.
    wrong = [
        f"{name!r} ({int(length)} bp in the map, {sizes[name]} bp in the contig data)"
        for name, length in zip(names, contigs["length"], strict=True)
        if name in sizes and int(length) != sizes[name]
    ]
    if wrong:
        raise MetacoolError(
            f"{metacool} is not a contact map of the contigs of {contig_data}: {len(wrong)} contig(s) have another "
            f"length in it: {_shown(wrong, show=str)}."
        )
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
    genomes_table = bin_summary if bin_summary and final_bins else None
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
            with _as_metacool_errors():
                metacooler.reassign(
                    metacool,
                    assignment,
                    staged,
                    # The contigs of no final bin go to __unbinned__,
                    # whatever genome the input map gave them.
                    missing="unbinned",
                    # Grouped as write_metacool groups the contigs, whatever
                    # the order of the input map: the pixels are remapped.
                    layout="genome-then-length",
                    # The bin names are the first column of bin_summary.txt,
                    # unnamed. Given, it replaces the columns of the genomes
                    # table that the input map carried.
                    genomes_table=genomes_table,
                    genomes_table_format="1",
                    metadata={"metator": {"version": __version__, **(metadata or {})}},
                    # Over the empty temporary file.
                    force=True,
                )
                _write_back_columns(staged, contig_data, contig_columns, header, genomes_table is None)
                metacooler.validate(staged, full=True)
                os.replace(staged, outfile)
            logger.info(f"{outfile} is written: the contact map of the bins.")
        except BaseException as error:
            # metacooler leaves its own partial file next to the one it
            # writes, and logs, at ERROR, that it is left there for
            # inspection: the error then says that it is removed, or, for
            # an error of another type (raised as it is), a warning.
            partial = f"{staged}.partial"
            left = os.path.lexists(partial)
            for path in [staged, partial]:
                if os.path.lexists(path):
                    os.remove(path)
            if left:
                removed = f"The temporary files of the map, {staged} and metacooler's partial file {partial}, were removed."
                if isinstance(error, MetacoolError):
                    raise MetacoolError(f"{str(error).rstrip('.')}. {removed}") from error.__cause__
                logger.warning(removed)
            raise


def _write_back_columns(path, contig_data, contig_columns, header, drop_genome_columns):
    """Give the map written back from an input map the metadata columns of
    this run, as metacooler annotate edits them.

    The columns of the contig data that metaTOR measured (contig_columns)
    replace those that the input map carried. Its other columns (header:
    all of them) that the
    input map carried are dropped: the bins of each step, the contig ids,
    Size..., which a map made with the contig data of an earlier run as its
    contigs table holds (metacooler create pairs --contigs-format metator),
    and one of CONTIG_COLUMNS that this run did not measure
    (Shotgun_coverage without a depth file...). So are, with
    drop_genome_columns (no bin summary given), the metadata columns of the
    genomes table: they describe the bins of another run.
    """
    import metacooler

    stored = metacooler.MetaCooler(path)
    # The columns of a contig data file but Name (the contigs' own name
    # column) and those of the map that metacooler derives.
    of_metator = (set(header) | set(CONTIG_COLUMNS)) - {"Name"}
    drop = [
        f"contigs.{column}"
        for column in stored.contigs.columns
        if column in of_metator - set(contig_columns) - {"name", "length", "genome_id"}
    ]
    if drop_genome_columns:
        drop += [
            f"genomes.{column}"
            for column in stored.genomes.columns
            if column not in DERIVED_GENOME_COLUMNS and not column.startswith("qc_")
        ]
    if not contig_columns and not drop:
        return
    with _rows_of_other_contigs_ignored(contig_data):
        metacooler.annotate(
            path,
            contigs_table=contig_data if contig_columns else None,
            contigs_table_format="Name",
            # Named one by one: the other columns of the table are bins, or
            # the unnamed pandas index.
            columns=[f"contigs.{column}" for column in contig_columns] or None,
            drop=drop or None,
            # The values of this run, on every contig of the map: a contig
            # that the contig data lacks has none.
            replace=True,
            missing="na",
            inplace=True,
        )


@contextmanager
def _rows_of_other_contigs_ignored(contig_data):
    """Drop the warning of metacooler annotate about the rows of contig_data
    that name no contig of the map.

    The contig data of a run lists every contig of the assembly, and a map
    may lack some of them (a map made with a contig filter): annotate would
    warn about its rows that name them, "Contigs table (contig_data): N of M
    row(s) match no contig of the map". Its other messages are kept, and
    annotate reads contig_data itself: they name it, as does the history of
    the map (metadata["annotate"]).
    """
    # The wording of metacooler 0.2.0 (metacooler.io.collect, _resolve_contig_metadata).
    pattern = re.compile(rf"Contigs table \({re.escape(str(contig_data))}\): \d+ of \d+ row\(s\) match no contig\b")

    def kept(record):
        return not (record.levelno == logging.WARNING and pattern.match(record.getMessage()))

    mc_logger = logging.getLogger("metacooler")
    mc_logger.addFilter(kept)
    try:
        yield
    finally:
        mc_logger.removeFilter(kept)
