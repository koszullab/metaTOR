#!/usr/bin/env python3
# coding: utf-8

"""The marker genes of the contigs of a .metacool map, as miComplete finds
those of a bin: metator markers.

miComplete gives a bin four numbers from the single-copy marker genes found
in it: its completeness, redundancy, weighted completeness and weighted
redundancy. Each comes from the marker genes of the bin, and each marker
gene lies on one contig: stored per contig, the markers of a genome give its
four numbers again, with no gene to find (see marker_quality), exactly those
that miComplete gives it as the markers were found (trained in single genome
mode on its contigs, in the order of the map, see below). Those of any other
group of contigs (two genomes merged, a contig moved from one to another,
contigs of no genome recruited) give an estimate of the numbers of a bin of
them, from the genes found in the genomes of the map, which miComplete's own
run on the group can differ from: it trains pyrodigal on the group itself,
and counts the E-values of the hits over its proteins. Of 28 curated groups
of the mock community MM133 (halves of MAGs, merges of two MAGs, MAGs with
contigs of no genome), 23 came out exact; the merge of the MAGs of two
organisms (GC 66.5% and 39.0%) has a redundancy of 1.6667 from the stored
markers (175 copies of the 105 markers of Bact105), 1.1827 from miComplete,
which, trained on the two together, finds 123 copies of 104 of them. This is
what ``metator markers`` stores in a map (see metator.metacool.write_map), so
that a viewer recomputes the numbers of a bin as a curation edits it: those
of a genome of the map, and an estimate of those of a bin the curation
made.

The genes and their markers are found as metaTOR's miComplete finds them in
a bin (see metator._micomplete): pyrodigal in its single genome mode,
trained on all the contigs of the bin together, in their order, then run on
each contig; then pyhmmer's hmmsearch of the marker HMMs of a set (Bact105:
105 bacterial markers, Arch131: 131 archaeal ones) in the proteins of the
genes, with miComplete's E-value threshold (--evalue, 4e-10), its E-values
counted over the proteins of the bin (Z, the number of proteins). miComplete
reads the hits from the table that hmmsearch writes (--tblout), and drops a
hit whose bias is too high for its score (--bias, 0.3) or whose best domain
is much less significant than the whole protein (--domain-cutoff, 1e-5):
read here from the same table, written by pyhmmer as metaTOR's miComplete
writes it, so that the hits are those of miComplete, its rounding of the
table's numbers included (see read_hits).

On a map, each genome of the map (not ``__unbinned__``) is a bin: its contigs
are trained on in the order of the map, the order in which metaTOR writes
the FASTA file of a genome (``metacooler fasta``, ``fasta/mags/``). Summed
over the contigs of such a genome, the markers give the four numbers that
miComplete gives that FASTA file. miComplete finds no gene in a bin of less
than 20 kb, which prodigal cannot train on: the contigs of a genome that
small, and those of no genome, have their genes found in metagenome mode
(pyrodigal's pre-trained models, contig by contig), and their markers
searched with Z the median number of proteins of the genomes trained in
single mode: every genome of the map of 20 kb or more, whatever it is, an
MGE genome of metator mge --reassign too, whose fewer proteins lower it.
The genomes, and groups of the other contigs, are searched side by side in
threads: pyrodigal and pyhmmer release the GIL, and each search runs on one
CPU, as miComplete runs hmmsearch in each of its processes.

The markers of a contig, in a set, are written as a text (see
encode_markers): its markers sorted by name, comma-separated, each with the
number of its proteins that miComplete keeps when above one
(``Ribosomal_L22,SecY*2``); ``*0`` for a marker that miComplete sees on the
contig (a hit in the table of hmmsearch) but whose every hit it drops, which
still counts in its weighted completeness, as miComplete counts the markers
it sees there, and in none of its other numbers. Both columns hold text, an
empty one for a contig without a marker.
"""

import functools
import io
import os
import statistics
import threading
import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np
import pandas as pd

import metator.io as mio
import metator.metacool as mtmc
from metator.log import logger

#: The marker sets of miComplete, in the order in which metaTOR keeps them on a
#: tie (see kept_quality): bacterial, then archaeal.
MARKER_SETS = ("Bact105", "Arch131")

#: What metaTOR calls a bin of each marker set, in the column Markers of the
#: results of miComplete (see metator.validation.merge_micomplete).
KINGDOMS = {"Bact105": "Bacteria", "Arch131": "Archaea"}

#: miComplete's thresholds, as metaTOR runs it (the defaults of its command
#: line: metator validation gives no --evalue, --bias or --domain-cutoff).
EVALUE = 4e-10
BIAS = 0.3
DOMAIN_CUTOFF = 1e-5

#: The length that pyrodigal needs to train on (prodigal's 20 kb): the
#: contigs of a bin joined as prodigal joins them, a linker of 12 bp after
#: each when there are several. A bin under it gets no gene from miComplete.
MIN_TRAINING = 20_000

#: How the genes of the contigs of a map are found (see contig_markers), in
#: the words of the markers entry of the metadata of the map (genes, see
#: metator.cli.markers_entry, beside the numbers that it names).
GENE_RULES = {
    "single": "Each genome but __unbinned__ whose contigs make min_training bp or more, joined as prodigal joins them (a "
    "12 bp linker after each when there are several): pyrodigal trained on its contigs in map order, the E-values of "
    "their hits counted over the proteins of the genome, as miComplete finds the markers of a bin.",
    "metagenome": "The contigs of __unbinned__ and of the shorter genomes: pyrodigal's metagenome mode, contig by "
    "contig, the E-values of their hits counted over metagenome_z proteins: the median number of proteins of the "
    "genomes trained in single mode, every one of them, MGE genomes included (without one, the proteins of these "
    "contigs, as one bin).",
}

#: How the markers of a contig are written in the contigs table of a map (see
#: encode_markers): the markers apart, and a marker from its count.
SEPARATOR = ","
COUNT = "*"

#: How many groups of contigs per thread the contigs in metagenome mode are
#: searched in: smaller groups share the threads better.
GROUPS_PER_THREAD = 4


def marker_column(name):
    """The column of the contigs table of a map that holds the markers of each
    contig in the marker set name (see encode_markers): markers_<name>."""
    return f"markers_{name}"


@dataclass(frozen=True)
class MarkerSet:
    """A marker set of miComplete (see marker_set).

    Attributes
    ----------
    name : str
        Its name: Bact105 or Arch131.
    hmms : str
        The path to its HMMs (share/<name>.hmm of miComplete).
    markers : tuple of str
        The names of its markers, in the order of its HMMs.
    weights : tuple of float
        The weight of each marker, in the same order, as miComplete reads them
        (share/<name>.weights)."""

    name: str
    hmms: str
    markers: tuple
    weights: tuple

    def entry(self):
        """The set as the markers entry of the metadata of a map records it: its
        markers in the order of its HMMs, and their weights."""
        return {"markers": list(self.markers), "weights": list(self.weights)}


@functools.cache
def marker_set(name):
    """The marker set name of miComplete (Bact105 or Arch131), from the files
    that miComplete installs: its markers named as miComplete names them (the
    NAME lines of the HMMs), and their weights (the lines name<TAB>weight of
    the weights file, but its Standard deviation line)."""
    import importlib.util

    # Found, not imported: miComplete imports names that recent Biopython
    # lacks, which only metator._micomplete puts back, in its own process.
    spec = importlib.util.find_spec("micomplete")
    if spec is None or spec.origin is None:
        raise ImportError("miComplete is not installed: metaTOR reads its marker sets.")
    share = os.path.join(os.path.dirname(os.path.abspath(spec.origin)), "share")
    hmms = os.path.join(share, f"{name}.hmm")
    with open(hmms) as handle:
        markers = tuple(line.split(" ")[2].strip() for line in handle if line.startswith("NAME"))
    weights = {}
    with open(os.path.join(share, f"{name}.weights")) as handle:
        for line in handle:
            fields = line.split("\t")
            if fields[0] != "Standard deviation:":
                weights[fields[0]] = float(fields[1])
    return MarkerSet(name, hmms, markers, tuple(weights[marker] for marker in markers))


# The HMMs of each marker set, read once by each thread that searches them
# (see _hmms).
_LOCAL = threading.local()


def _hmms(markers):
    """The profile HMMs of a marker set (a MarkerSet), as pyhmmer reads them: read
    once by each thread, which searches them on its own."""
    import pyhmmer

    held = getattr(_LOCAL, "hmms", None)
    if held is None:
        held = _LOCAL.hmms = {}
    if markers.hmms not in held:
        with pyhmmer.plan7.HMMFile(markers.hmms) as hmm_file:
            held[markers.hmms] = tuple(hmm_file)
    return held[markers.hmms]


# ---------------------------------------------------------------------------
# The genes and the markers of a group of contigs
# ---------------------------------------------------------------------------


@dataclass
class Proteins:
    """The proteins of the genes of a group of contigs (see find_proteins).

    Attributes
    ----------
    fasta : str
        The proteins, as pyrodigal writes them, and prodigal (-a): one record
        per gene, ``<contig>_<number of the gene in the contig>``, its
        translation ending with ``*`` for a complete gene.
    contigs : list of str
        The contig of each protein, in the order of the records."""

    fasta: str
    contigs: list

    def __len__(self):
        return len(self.contigs)


def find_proteins(contigs, meta=False):
    """The proteins of the genes of a group of contigs, as metaTOR's miComplete
    finds those of a bin (see metator._micomplete.create_proteome): pyrodigal
    in its single genome mode, trained on all the contigs together, joined in
    their order as prodigal joins them, then run on each contig; or, with meta,
    in its metagenome mode, contig by contig.

    Parameters
    ----------
    contigs : list of (str, bytes)
        The name and the sequence of each contig, in order.
    meta : bool
        Whether to find the genes in metagenome mode.

    Returns
    -------
    Proteins, or None
        None in single genome mode when the contigs, joined, are too short to
        train on (MIN_TRAINING): miComplete finds no gene in such a bin."""
    import pyrodigal

    finder = pyrodigal.GeneFinder(meta=meta)
    if not meta:
        if not contigs:
            return None
        try:
            # Under 100 kb, pyrodigal warns (see contig_markers, which silences
            # it): prodigal warns too, on a stderr that miComplete discards.
            finder.train(*(sequence for _, sequence in contigs))
        except ValueError:
            return None
    fasta = io.StringIO()
    owners = []
    for name, sequence in contigs:
        genes = finder.find_genes(sequence)
        genes.write_translations(fasta, sequence_id=name)
        owners.extend([name] * len(genes))
    return Proteins(fasta.getvalue(), owners)


def search_table(proteins, markers, evalue=EVALUE, z=None):
    """The hits of the HMMs of a marker set (a MarkerSet) in proteins (a
    Proteins), as metaTOR's miComplete runs hmmsearch on a bin (see
    metator._micomplete.hmm_search): ``hmmsearch -E <evalue> --tblout``, the
    table written by pyhmmer as hmmsearch writes it (its target E-values per
    protein of the search, Z: by default, the number of proteins searched).

    Returns
    -------
    str
        The table (--tblout), its header once: empty without a protein."""
    import pyhmmer

    if not len(proteins):
        return ""
    hmms = _hmms(markers)
    with pyhmmer.easel.SequenceFile(
        io.BytesIO(proteins.fasta.encode()), format="fasta", digital=True, alphabet=hmms[0].alphabet
    ) as sequences:
        block = sequences.read_block()
    # One CPU: the genomes are searched side by side, in threads.
    if z is None:
        searches = pyhmmer.hmmsearch(hmms, block, cpus=1, E=evalue)
    else:
        searches = pyhmmer.hmmsearch(hmms, block, cpus=1, E=evalue, Z=z)
    table = io.BytesIO()
    for i, hits in enumerate(searches):
        hits.write(table, format="targets", header=i == 0)
    return table.getvalue().decode()


def read_hits(table, bias=BIAS, domain_cutoff=DOMAIN_CUTOFF):
    """The hits of a table of hmmsearch (--tblout) as miComplete reads them
    (micomplete.completeness.calcCompleteness.get_completeness): the target
    protein, the marker (the query) and whether miComplete keeps the hit, from
    the numbers as the table writes them (two significant digits for an
    E-value, one decimal for a score): it drops a hit whose score times bias
    is no more than its bias, or whose best domain has an E-value more than
    domain_cutoff above that of the whole protein (suspicion_check). Every hit
    of the table counts in the markers miComplete sees, kept or dropped.

    Returns
    -------
    list of (str, str, bool)
        The protein, the marker and whether the hit is kept, in the order of
        the table."""
    hits = []
    for line in table.splitlines():
        # miComplete stops at a line "#", hmmsearch's footer (pyhmmer writes
        # none); the header's lines name no marker.
        if line == "#":
            break
        if line.startswith("#"):
            continue
        fields = line.split()
        full_evalue, score, hit_bias, domain_evalue = (float(value) for value in fields[4:8])
        dropped = score * bias <= hit_bias or domain_evalue - full_evalue > domain_cutoff
        hits.append((fields[0], fields[2], not dropped))
    return hits


def count_markers(proteins, hits):
    """The markers of each contig of a group (see find_proteins), from the hits
    of a marker set in their proteins (see read_hits): by contig, the number of
    its proteins that miComplete keeps for each marker that it sees on the
    contig, 0 for a marker whose every hit there it drops. A protein counts for
    each marker that it is kept for.

    Returns
    -------
    dict
        By contig (those with a marker only): a dict of the count of each
        marker."""
    owner = {}
    for line in proteins.fasta.splitlines():
        if line.startswith(">"):
            owner[line[1:].split(None, 1)[0]] = len(owner)
    found = {}
    for protein, marker, kept in hits:
        markers = found.setdefault(proteins.contigs[owner[protein]], {})
        markers[marker] = markers.get(marker, 0) + int(kept)
    return found


@dataclass
class GroupMarkers:
    """The markers of the contigs of a group (see group_markers).

    Attributes
    ----------
    proteins : int, or None
        The number of proteins of the group; None when it is too short to train
        on (in single genome mode).
    markers : dict
        By marker set, the markers of each contig that has some (see
        count_markers)."""

    proteins: int | None
    markers: dict = field(default_factory=dict)


def group_markers(
    contigs, sets=MARKER_SETS, meta=False, z=None, evalue=EVALUE, bias=BIAS, domain_cutoff=DOMAIN_CUTOFF, proteins=None
):
    """The markers of each contig of a group of contigs, found as metaTOR's
    miComplete finds those of a bin of these contigs, in their order (see
    find_proteins, search_table and read_hits): in single genome mode, or in
    metagenome mode with meta; the E-values counted over z proteins (by
    default, the proteins of the group: those of the bin).

    Parameters
    ----------
    contigs : list of (str, bytes)
        The name and the sequence of each contig, in order.
    sets : sequence of str
        The marker sets (see MARKER_SETS).
    meta : bool
        Whether to find the genes in metagenome mode.
    z : float, optional
        The number of proteins that the E-values are counted over.
    evalue, bias, domain_cutoff : float
        miComplete's thresholds (--evalue, --bias, --domain-cutoff).
    proteins : Proteins, optional
        The proteins of the group, when found already.

    Returns
    -------
    GroupMarkers
        Its proteins and the markers of its contigs; without markers, and
        proteins None, for a group too short to train on (in single genome
        mode)."""
    if proteins is None:
        proteins = find_proteins(contigs, meta=meta)
    if proteins is None:
        return GroupMarkers(None)
    found = GroupMarkers(len(proteins))
    for name in sets:
        table = search_table(proteins, marker_set(name), evalue, z)
        found.markers[name] = count_markers(proteins, read_hits(table, bias, domain_cutoff))
    return found


# ---------------------------------------------------------------------------
# The markers of the contigs of a map
# ---------------------------------------------------------------------------


@dataclass
class ContigMarkers:
    """The markers of the contigs of a map (see contig_markers).

    Attributes
    ----------
    table : pandas.DataFrame
        One row per contig of the map, in its order: its name (Name) and its
        markers in each set (markers_Bact105, markers_Arch131; see
        encode_markers), an empty text for none.
    trained : dict
        The genomes whose genes were found in single genome mode, trained on
        their contigs: the number of proteins of each, by name, in the order of
        the genomes table.
    small : list of str
        The genomes too short to train on, whose contigs are in metagenome
        mode, in the order of the genomes table.
    meta : int
        The number of contigs in metagenome mode: those of the small genomes and
        of no genome (__unbinned__).
    z : float, or None
        The number of proteins over which the E-values of the contigs in
        metagenome mode are counted (see metagenome_z); None without such
        contigs."""

    table: pd.DataFrame
    trained: dict
    small: list
    meta: int
    z: float | None


def metagenome_z(trained, proteins):
    """The number of proteins that the E-values of the contigs in metagenome mode
    are counted over: the median number of proteins of the genomes trained in
    single mode (trained, by genome), every genome of the map that pyrodigal
    trains on, whatever it is, an MGE genome of metator mge --reassign too,
    whose fewer proteins lower it; without a trained genome, proteins, the
    number of proteins of the contigs in metagenome mode, searched as one
    bin."""
    # At least one protein: HMMER counts no E-value over none.
    return max(1.0, float(statistics.median(trained.values())) if trained else float(proteins))


def read_contigs(metacool, mc, assembly):
    """The sequences of the contigs of a map (mc, opened) from its assembly, as
    metaTOR's miComplete reads those of a bin (Biopython: a record's name is
    the first word of its header, and its sequence its bases as they are, in
    upper or lower case).

    Parameters
    ----------
    metacool : str
        The map, for the messages.
    mc : metacooler.MetaCooler
        The map, opened.
    assembly : str
        The assembly that the map was made from: a FASTA file, plain or
        compressed with gzip.

    Returns
    -------
    dict
        The sequence of each contig of the map, as bytes, by name.

    Raises
    ------
    metator.metacool.MetacoolError
        If the assembly names a contig by several records, or the map is not a
        map of the contigs of the assembly: a contig that the assembly lacks,
        or of another length there."""
    from Bio import SeqIO

    with mtmc.read_errors(metacool):
        contigs = mc.contigs
    wanted = set(contigs["name"].astype(str))
    sequences, lengths, named_again = {}, {}, []
    with mio.read_compressed(assembly) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id in lengths:
                named_again.append(record.id)
                continue
            lengths[record.id] = len(record.seq)
            if record.id in wanted:
                sequences[record.id] = bytes(record.seq)
    if named_again:
        named_again = list(dict.fromkeys(named_again))
        raise mtmc.MetacoolError(
            f"{len(named_again)} contig(s) are named by several records of the assembly {assembly}: "
            f"{mtmc._shown(named_again)}. A contig is read by its name, the first word of its header: give each record "
            "a name of its own."
        )
    missing = [str(name) for name in contigs["name"] if str(name) not in lengths]
    if missing:
        raise mtmc.MetacoolError(
            f"{metacool} is not a contact map of the contigs of the assembly: {len(missing)} of its contig(s) are not "
            f"in the assembly {assembly}: {mtmc._shown(missing)}. Give the assembly that the map was made from."
        )
    mtmc.check_map_lengths(metacool, contigs, lengths, f"the assembly {assembly}", "the assembly")
    return sequences


def _in_threads(function, jobs, threads):
    """function on each of jobs, in threads (threads at most), the results in the
    order of jobs. Threads, not processes: pyrodigal and pyhmmer find the genes
    and search the HMMs without holding the GIL, so the threads run side by
    side, and no process is forked from this one, which runs threads, nor
    spawned, which would import the script of a caller of the Python API
    again (see metator.partition.get_hamming_distance). On a map of 20,566
    contigs and 66 genomes to train on, 4 threads find the markers 3.5 times
    faster than one, and faster than 4 processes, which the sequences are
    pickled to."""
    if threads <= 1 or len(jobs) <= 1:
        return [function(job) for job in jobs]
    with ThreadPoolExecutor(max_workers=min(threads, len(jobs))) as executor:
        return list(executor.map(function, jobs))


def _by_size(groups):
    """The keys of groups (lists of (name, sequence)), the longest group first:
    the threads that take them in turn end at about the same time."""
    return sorted(groups, key=lambda key: -sum(len(sequence) for _, sequence in groups[key]))


def _split(contigs, parts):
    """contigs (a list of (name, sequence)) in at most parts groups of about the
    same length, each in the order of contigs."""
    total = sum(len(sequence) for _, sequence in contigs)
    size = max(1, -(-total // max(1, parts)))
    groups, group, length = [], [], 0
    for contig in contigs:
        group.append(contig)
        length += len(contig[1])
        if length >= size:
            groups.append(group)
            group, length = [], 0
    if group:
        groups.append(group)
    return groups


def contig_markers(metacool, mc, sequences, threads=1, evalue=EVALUE, bias=BIAS, domain_cutoff=DOMAIN_CUTOFF):
    """The markers of every contig of a map, in each marker set of miComplete
    (see the module's docstring): each genome of the map (not __unbinned__)
    as a bin of its contigs, in the order of the map, trained on in single
    genome mode, its E-values counted over its proteins; the contigs of the
    genomes too short to train on, and of no genome, in metagenome mode, their
    E-values counted over the median number of proteins of the genomes
    trained in single mode, MGE genomes included (see metagenome_z).

    Parameters
    ----------
    metacool : str
        The map, for the messages.
    mc : metacooler.MetaCooler
        The map, opened: its genomes and contigs.
    sequences : dict
        The sequence of each contig of the map, by name (see read_contigs).
    threads : int
        The number of genomes, or groups of contigs, searched side by side.
    evalue, bias, domain_cutoff : float
        miComplete's thresholds (--evalue, --bias, --domain-cutoff).

    Returns
    -------
    ContigMarkers
        The markers of the contigs, and how their genes were found."""
    with mtmc.read_errors(metacool):
        contigs = mc.contigs
        genome_names = mc.genomes["name"].astype(str).to_numpy()
    names = contigs["name"].astype(str).to_numpy()
    of_genome = genome_names[contigs["genome_id"].to_numpy()]
    # The contigs of each genome, in the order of the map, the genomes in the
    # order of the genomes table.
    groups = {name: [] for name in genome_names if name != mtmc.UNBINNED_GENOME}
    for name, genome in zip(names, of_genome, strict=True):
        if genome != mtmc.UNBINNED_GENOME:
            groups[genome].append((name, sequences[name]))
    groups = {genome: group for genome, group in groups.items() if group}
    options = {"evalue": evalue, "bias": bias, "domain_cutoff": domain_cutoff}
    with warnings.catch_warnings():
        # A genome under 100 kb, as many are: pyrodigal warns, and prodigal,
        # on a stderr that miComplete discards. Set here, once, for the
        # threads: catch_warnings is not safe in threads.
        warnings.filterwarnings("ignore", message="sequence should be at least", category=UserWarning)
        order = _by_size(groups)
        found = dict(
            zip(order, _in_threads(lambda genome: group_markers(groups[genome], **options), order, threads), strict=True)
        )
        trained = {genome: found[genome].proteins for genome in groups if found[genome].proteins is not None}
        small = [genome for genome in groups if found[genome].proteins is None]
        # The contigs in metagenome mode, in the order of the map: their genes
        # first, then their markers, once the number of proteins that their
        # E-values are counted over is known.
        small_genomes = set(small)
        meta = [
            (name, sequences[name])
            for name, genome in zip(names, of_genome, strict=True)
            if genome == mtmc.UNBINNED_GENOME or genome in small_genomes
        ]
        parts = _split(meta, max(1, threads) * GROUPS_PER_THREAD) if meta else []
        proteins = _in_threads(lambda part: find_proteins(part, meta=True), parts, threads)
        z = metagenome_z(trained, sum(len(found_proteins) for found_proteins in proteins)) if meta else None
        jobs = list(zip(parts, proteins, strict=True))
        searched = _in_threads(lambda job: group_markers(job[0], meta=True, z=z, proteins=job[1], **options), jobs, threads)
    markers = {name: {} for name in MARKER_SETS}
    for group in [*found.values(), *searched]:
        for name, by_contig in group.markers.items():
            markers[name].update(by_contig)
    table = pd.DataFrame({"Name": names})
    for name in MARKER_SETS:
        table[marker_column(name)] = [encode_markers(markers[name].get(contig, {})) for contig in names]
    result = ContigMarkers(table, trained, small, len(meta), z)
    _log_markers(result)
    return result


def _log_markers(found):
    """Log how the genes of the contigs of a map were found (found, a
    ContigMarkers), and how many contigs have markers of each set."""
    logger.info(
        f"The genes of {len(found.trained)} genome(s) were found as miComplete finds those of a bin, by pyrodigal trained "
        f"on their contigs; those of {found.meta} contig(s), in no genome or in {len(found.small)} genome(s) too short to "
        "train on (under 20 kb), in metagenome mode"
        + (f", their marker genes searched as in a genome of {found.z:g} proteins." if found.z is not None else ".")
    )
    for name in MARKER_SETS:
        column = found.table[marker_column(name)]
        with_markers = column != ""
        genes = sum(sum(decode_markers(text).values()) for text in column[with_markers])
        logger.info(f"{name}: {genes} marker gene(s) on {int(with_markers.sum())} of the {len(column)} contigs.")


# ---------------------------------------------------------------------------
# The markers of a contig as a text, and the numbers of miComplete
# ---------------------------------------------------------------------------


def encode_markers(counts):
    """The markers of a contig (counts: the number of proteins of each, see
    count_markers) as the contigs table of a map stores them: sorted by name
    (in the order of their bytes, ASCII), comma-separated, each followed by
    ``*`` and its count when that is not 1 (``Ribosomal_L22,SecY*2``; ``*0``
    for a marker that miComplete sees but drops); an empty text for none."""
    return SEPARATOR.join(name if count == 1 else f"{name}{COUNT}{count}" for name, count in sorted(counts.items()))


def decode_markers(text):
    """The markers of a contig from their text (see encode_markers): the count of
    each, by name. An empty text is none, and so is a missing value (None or
    NaN), which a contigs table edited by another tool may hold: metaTOR
    stores both columns as text (see metator.metacool._store_markers)."""
    if text is None or (isinstance(text, float) and np.isnan(text)) or text == "":
        return {}
    counts = {}
    for item in str(text).split(SEPARATOR):
        name, star, count = item.partition(COUNT)
        counts[name] = int(count) if star else 1
    return counts


def sum_markers(counts):
    """The markers of a group of contigs: those of its contigs (counts, an
    iterable of dicts, see decode_markers), summed by marker. A marker that
    miComplete sees on a contig but drops (0) stays, as seen, at 0."""
    total = {}
    for markers in counts:
        for name, count in markers.items():
            total[name] = total.get(name, 0) + count
    return total


class Quality(NamedTuple):
    """The four numbers that miComplete gives a bin, as it writes them (see
    marker_quality)."""

    completeness: float
    redundancy: float
    weighted_completeness: float
    weighted_redundancy: float


def marker_quality(counts, markers, weights):
    """The completeness, redundancy, weighted completeness and weighted
    redundancy that miComplete gives a bin whose markers are counts, in a
    marker set, as miComplete computes and rounds them
    (micomplete.micomplete._compile_results and
    micomplete.completeness.calcCompleteness.attribute_weights).

    With P the markers present (a count of one or more) and N the markers of
    the set: completeness = P / N, and redundancy = the sum of the counts of
    the present markers / P, each rounded to 4 decimals (by Python's round,
    as miComplete: an exact tie to the even digit, 1.03125 to 1.0312). The weighted
    completeness WC is the sum of the weights of the markers that miComplete
    sees (present, or seen and dropped: a count of 0), rounded to 4 decimals,
    0.0001 if that is 0; the weighted redundancy is (the sum of the weights of
    the markers of more than one copy, each once, + WC) / WC, WC not rounded,
    the quotient rounded to 4 decimals. A bin without a present marker has 0
    for all four.

    Parameters
    ----------
    counts : dict
        The count of each marker of the bin, by name (see sum_markers): the
        markers that miComplete sees in it only.
    markers : sequence of str
        The markers of the set (see MarkerSet.markers).
    weights : sequence of float
        The weight of each marker of the set, in the order of markers.

    Returns
    -------
    Quality
        The four numbers.

    Raises
    ------
    ValueError
        If counts names a marker that is not in the set."""
    weight = dict(zip(markers, weights, strict=True))
    unknown = [name for name in counts if name not in weight]
    if unknown:
        raise ValueError(f"Markers of no marker set of {len(markers)} markers: {', '.join(sorted(unknown))}.")
    present = [name for name in markers if counts.get(name, 0) > 0]
    if not present:
        return Quality(0.0, 0.0, 0.0, 0.0)
    completeness = round(len(present) / len(markers), 4)
    redundancy = round(sum(counts[name] for name in present) / len(present), 4)
    seen = sum(weight[name] for name in markers if name in counts)
    duplicated = sum(weight[name] for name in present if counts[name] > 1)
    weighted_redundancy = round((duplicated + seen) / seen, 4)
    weighted_completeness = round(seen, 4) or 0.0001
    return Quality(completeness, redundancy, weighted_completeness, weighted_redundancy)


def kept_quality(qualities):
    """The marker set that metaTOR keeps for a bin, and its numbers, from the
    numbers of each set (qualities, by set name; see marker_quality), as
    metator.validation.merge_micomplete keeps miComplete's results: Bact105,
    unless Arch131 has the higher weighted completeness.

    Returns
    -------
    tuple
        The set's name for metaTOR (Bacteria or Archaea, see KINGDOMS), and its
        Quality."""
    bacteria, archaea = qualities["Bact105"], qualities["Arch131"]
    if bacteria.weighted_completeness >= archaea.weighted_completeness:
        return KINGDOMS["Bact105"], bacteria
    return KINGDOMS["Arch131"], archaea


def quality_of(texts):
    """The numbers that miComplete gives a group of contigs, from their markers
    as a map stores them (texts: by marker set, the markers of each contig, see
    encode_markers), as metaTOR keeps them (see kept_quality).

    Returns
    -------
    tuple
        Bacteria or Archaea, and the Quality of that set."""
    qualities = {}
    for name in MARKER_SETS:
        markers = marker_set(name)
        counts = sum_markers(decode_markers(text) for text in texts[name])
        qualities[name] = marker_quality(counts, markers.markers, markers.weights)
    return kept_quality(qualities)


def entry_sets():
    """The marker sets as the markers entry of the metadata of a map records
    them (see metator.cli.markers_entry): by set, its column in the contigs
    table, its markers in the order of its HMMs and their weights."""
    return {name: {"column": marker_column(name), **marker_set(name).entry()} for name in MARKER_SETS}
