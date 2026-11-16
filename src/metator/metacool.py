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

Core function to write the file:
    - write_metacool
"""

import logging
import os
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


def _read_contig_table(contig_data, **kwargs):
    """A contig data file as a table of text (kwargs: more keyword arguments
    of pandas.read_csv). A file that pandas cannot read as a table (empty,
    not UTF-8, a row with more fields than the header...) is a
    MetacoolError naming the file, as are the other refusals of the inputs.
    """
    try:
        return pd.read_csv(contig_data, sep="\t", dtype=str, keep_default_na=False, **kwargs)
    except (ValueError, OSError) as error:
        raise MetacoolError(f"{contig_data} cannot be read as a table: {str(error).strip().rstrip('.')}.") from error


def read_contig_bins(contig_data):
    """Bins (one per contig) and final bin of the contigs of a contig data
    file.

    Parameters
    ----------
    contig_data : str
        Path to the contig_data_final.txt file of a validation or pipeline
        run: its Name, Size and Final_bin columns are read.

    Returns
    -------
    pandas.DataFrame
        One bin per contig, from 0 to its size: columns chrom, start and
        end.
    dict
        Final bin of each binned contig, by contig name. The unbinned
        contigs (Final_bin ND, - or empty) are left out.

    Raises
    ------
    MetacoolError
        If pandas cannot read the file as a table, or it has no Name, Size
        or Final_bin column, or a Size that is not a whole number.
    """
    table = _read_contig_table(contig_data)
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
    return bins, dict(zip(table["Name"][binned], final_bin[binned], strict=True))


def measured_contig_columns(contig_data):
    """The columns of CONTIG_COLUMNS that a contig data file has and fills.

    A column that is ``-`` on every contig is left out: metaTOR measured
    nothing there (Shotgun_coverage without a depth file, Restriction_site
    without an enzyme).

    Parameters
    ----------
    contig_data : str
        Path to the contig_data_final.txt file.

    Returns
    -------
    list of str
        The columns, in the order of CONTIG_COLUMNS.

    Raises
    ------
    MetacoolError
        If pandas cannot read the file as a table.
    """
    header = _read_contig_table(contig_data, nrows=0).columns
    present = [column for column in CONTIG_COLUMNS if column in header]
    table = _read_contig_table(contig_data, usecols=present)
    return [column for column in present if not table[column].str.strip().eq("-").all()]


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
            "A .metacool file cannot have bins named as metacooler marks a contig in no bin (NA, NaN, N/A, ND, N/D, -1 "
            f"or __unbinned__, in upper or lower case): it would put their contigs in __unbinned__. Rename the {len(bad)} "
            f"bin(s): {_shown(bad)}."
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
        measured_contig_columns).
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
    bins, contig_bin = read_contig_bins(contig_data)
    check_ascii_names(bins["chrom"], contig_bin.values())
    check_bin_names(contig_bin.values())
    contig_columns = measured_contig_columns(contig_data)

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
