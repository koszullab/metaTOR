#!/usr/bin/env python3
# coding: utf-8

"""The inputs of metator network and metator pipeline.

A run starts from the .metacool contact map of the contigs, or makes it from
the pairs files or the BAM files of its Hi-C libraries: one or several of
them, each argument a path, a comma-separated list or a glob pattern (see
expand), and what each file holds is read from its first bytes (see
file_kind). make_map then writes the map of the contigs, one bin per contig
of the assembly, all of them in no cluster yet, as ``metacooler create
pairs`` makes it from the pairs files: the BAM files give their pairs first
(see bam_pairs), into the temporary directory of the run.

A BAM file holds both mates of its reads (paired-end alignments, as ``bwa mem
-5SP`` writes them), or one mate of each (as bowtie2 aligns them, one mate at
a time, in two files per library). A pair is the primary alignments of the
two mates of a read, both aligned with a mapping quality of at least
--min-quality, each at the leftmost base of its alignment, as metaTOR 1.x
wrote them.
"""

import glob
import gzip
import os
import re
import shutil
import tempfile
from collections import OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass, field

import pandas as pd
import pysam

import metator.io as mio
import metator.metacool as mtmc
from metator.log import logger

#: What a file of the inputs of a run holds: the .metacool contact map of the
#: contigs, pairs files, or alignments (BAM, or SAM, files).
MAP, PAIRS, BAM = "map", "pairs", "BAM"

#: The two layouts of a BAM file: both mates of each read (paired-end
#: alignments), or one mate (the alignments of one mate at a time, in two files
#: per library, the R1 file and the R2 file).
BOTH_MATES, ONE_MATE = "both mates", "one mate"

#: The characters that make an item of an argument a glob pattern.
_GLOB_CHARS = re.compile(r"[*?[]")

#: The first line of a SAM file: a header line.
_SAM_HEADER = re.compile(rb"@(HD|SQ|RG|PG|CO)\t")

#: The first bytes of the compressions that a pairs file may have and that its
#: header is not read through (metacooler reads them): zstd and LZ4.
_UNREAD_COMPRESSIONS = (b"\x28\xb5\x2f\xfd", b"\x04\x22\x4d\x18")

#: How many alignments bam_layout reads to find the layout of a BAM file.
LAYOUT_PROBE = 10_000

#: How many names of reads with one mate a pass over a BAM file of both mates
#: remembers (the latest), to find that the file is not grouped by read name
#: (see _GroupingCheck).
GROUPING_WINDOW = 100_000

#: A mate suffix at the end of a read name, which some aligners keep.
_MATE_SUFFIX = re.compile(r"/[12]$")

#: The runs of digits of a read name, which samtools sort -n compares by value.
_DIGIT_RUNS = re.compile(r"([0-9]+)")


class _NameOrderError(Exception):
    """A BAM file is not sorted in the read name order it is walked in; its
    message is the path of the file."""


class _NotGrouped(Exception):
    """A BAM file of both mates is not grouped by read name: a read comes back
    after reads of other names."""


# ---------------------------------------------------------------------------
# The arguments
# ---------------------------------------------------------------------------


def expand(argument):
    """The files that an argument of a run names, in order.

    A path that exists is that file, whatever its name holds (a comma, a glob
    character). Otherwise the argument is a comma-separated list, each item a
    path or a glob pattern (``*``, ``?``, ``[...]``), which stands for the
    files it matches, sorted by name, as hicstuff 4.0 reads its inputs. A
    pattern that matches no file, and an empty item, are refused.

    Parameters
    ----------
    argument : str
        A path, comma-separated paths or glob patterns.

    Returns
    -------
    list of str
        The paths, which may not exist (see resolve).

    Raises
    ------
    metator.metacool.MetacoolError
        If a pattern matches no file, or the list holds an empty item.

    Examples
    --------
    >>> expand("L1.pairs.gz,L2.pairs.gz")
    ['L1.pairs.gz', 'L2.pairs.gz']
    """
    if os.path.exists(argument):
        return [argument]
    files = []
    for item in argument.split(","):
        if not item:
            raise mtmc.MetacoolError(f"The list of inputs {argument!r} holds an empty path.")
        if _GLOB_CHARS.search(item) is None:
            files.append(item)
            continue
        matches = sorted(glob.glob(item))
        if not matches:
            raise mtmc.MetacoolError(f"No file matches the input pattern {item!r}.")
        files.extend(matches)
    return files


def _first_bytes(path, size=64):
    """The first bytes of a file, of its content when it is compressed with gzip
    (BGZF, which BAM files are compressed with, is gzip): b"" for an empty or
    truncated one."""
    with open(path, "rb") as handle:
        head = handle.read(size)
    if head.startswith(b"\x1f\x8b"):
        try:
            with gzip.open(path, "rb") as handle:
                return handle.read(size)
        except (OSError, EOFError):
            return b""
    return head


def file_kind(path):
    """What a file of the inputs of a run holds, from its first bytes, whatever its
    name: MAP, PAIRS or BAM.

    An HDF5 file is a map (open_map refuses one that is not a .metacool map,
    with metacooler's reason), and so is file::group, a map inside a file,
    which is no path. A BAM file (compressed with BGZF) and a SAM file (text
    that starts with its header) are alignments. Any other file is a pairs
    file, plain or compressed, which metacooler reads. Reads (FASTQ), an
    assembly (FASTA) and a CRAM file are refused.

    Raises
    ------
    metator.metacool.MetacoolError
        If the file is not there, is a directory, or holds reads, an assembly
        or CRAM alignments.
    """
    path = str(path)
    if mtmc.is_map(path) or ("::" in path and not os.path.exists(path)):
        return MAP
    if not os.path.exists(path):
        raise mtmc.MetacoolError(f"Cannot read the input {path}: no such file.")
    if os.path.isdir(path):
        raise mtmc.MetacoolError(f"Cannot read the input {path}: it is a directory.")
    head = _first_bytes(path)
    if head.startswith(b"BAM\x01") or _SAM_HEADER.match(head):
        return BAM
    if head.startswith(b"CRAM"):
        raise mtmc.MetacoolError(
            f"{path} is a CRAM file, which metaTOR does not read: convert it to BAM (samtools view -b -T ASSEMBLY)."
        )
    if head.startswith(b"@"):
        raise mtmc.MetacoolError(
            f"{path} holds reads (FASTQ): metaTOR aligns no reads. Give the BAM files or the pairs files of the run, or "
            "the map of its contigs, which hicstuff makes from the reads (hicstuff pipeline -M metacool)."
        )
    if head.startswith(b">"):
        raise mtmc.MetacoolError(f"{path} is a FASTA file: give the assembly with -a, and the pairs or BAM files of the run.")
    return PAIRS


@dataclass
class Inputs:
    """The inputs of a run, as resolve finds them.

    Attributes
    ----------
    kind : str
        MAP, PAIRS or BAM.
    files : list of str
        The files, in the order of the arguments: the map alone, for MAP.
    libraries : list of tuple of str
        The files of each Hi-C library: a pairs file, a BAM file of both
        mates, or the BAM files of its two mates (R1, R2).
    layout : str or None
        For BAM files, BOTH_MATES or ONE_MATE.
    """

    kind: str
    files: list = field(default_factory=list)
    libraries: list = field(default_factory=list)
    layout: str | None = None

    @property
    def map(self):
        """The map, for MAP; None otherwise."""
        return self.files[0] if self.kind == MAP else None

    def describe(self):
        """What the inputs are, for the log: "a map", "a pairs file", "3 pairs
        files", "2 BAM files of both mates", "4 BAM files of one mate each"."""
        if self.kind == MAP:
            return "a map"
        what = "pairs file" if self.kind == PAIRS else "BAM file"
        files = f"a {what}" if len(self.files) == 1 else f"{len(self.files)} {what}s"
        if self.kind == BAM:
            files += " of both mates" if self.layout == BOTH_MATES else " of one mate each"
        return files


def resolve(arguments):
    """The inputs of a run (see Inputs) from its arguments: one .metacool map, or
    pairs files, or BAM files (or SAM files), each argument a path, a
    comma-separated list or a glob pattern (see expand), each file found to be
    what it is from its first bytes (see file_kind).

    Each pairs file, and each BAM file of both mates, is a library. BAM files
    of one mate each come in two halves, as hicstuff takes them: two
    arguments are the R1 files and the R2 files; more arguments (as the shell
    expands unquoted patterns) are as many R1 arguments, then as many R2
    arguments. The i-th R1 file and the i-th R2 file are the mates of one
    library.

    Raises
    ------
    metator.metacool.MetacoolError
        If an argument names no file, a file is given twice, a map is given
        with other inputs, kinds or layouts are mixed, or BAM files of one mate
        each do not come in two halves of as many files.
    """
    names = [expand(argument) for argument in arguments]
    files = [path for paths in names for path in paths]
    if not files:
        raise mtmc.MetacoolError("A run needs its inputs: the map of its contigs, or its pairs files or BAM files.")
    kinds = {path: file_kind(path) for path in files}
    seen = {}
    for path in files:
        real = os.path.realpath(mtmc.map_file(path)) if "::" in path else os.path.realpath(path)
        if real in seen:
            raise mtmc.MetacoolError(f"The input {path} is given twice: {seen[real]} is the same file.")
        seen[real] = path
    found = set(kinds.values())
    if MAP in found:
        if len(files) > 1:
            raise mtmc.MetacoolError(
                f"A .metacool map is the one input of a run, given alone, and this run has {len(files)}: {', '.join(files)}."
            )
        return Inputs(MAP, files, [tuple(files)])
    if found == {PAIRS, BAM}:
        bams = [path for path in files if kinds[path] == BAM]
        pairs = [path for path in files if kinds[path] == PAIRS]
        raise mtmc.MetacoolError(
            f"The inputs of a run are pairs files or BAM files, not both: BAM {mtmc._shown(bams)}; pairs {mtmc._shown(pairs)}."
        )
    (kind,) = found
    if kind == PAIRS:
        return Inputs(PAIRS, files, [(path,) for path in files])
    layouts = {path: bam_layout(path) for path in files}
    if len(set(layouts.values())) > 1:
        both = [path for path in files if layouts[path] == BOTH_MATES]
        one = [path for path in files if layouts[path] == ONE_MATE]
        raise mtmc.MetacoolError(
            f"The BAM files of a run hold both mates of each read, or one mate each (an R1 and an R2 file per "
            f"library), not both: both mates in {mtmc._shown(both)}; one mate in {mtmc._shown(one)}."
        )
    if layouts[files[0]] == BOTH_MATES:
        return Inputs(BAM, files, [(path,) for path in files], BOTH_MATES)
    if len(names) < 2 or len(names) % 2:
        raise mtmc.MetacoolError(
            f"One mate of each read in {', '.join(files)}: give the R1 files, then the R2 files, as two arguments (two "
            "comma-separated lists, or two patterns), or as many R1 arguments as R2 arguments."
        )
    half = len(names) // 2
    first = [path for paths in names[:half] for path in paths]
    second = [path for paths in names[half:] for path in paths]
    if len(first) != len(second):
        raise mtmc.MetacoolError(
            f"{len(first)} R1 file(s) and {len(second)} R2 file(s): each BAM file of the first mates needs the BAM file "
            f"of the second mates of its library. R1: {', '.join(first)}; R2: {', '.join(second)}."
        )
    return Inputs(BAM, files, list(zip(first, second, strict=True)), ONE_MATE)


# ---------------------------------------------------------------------------
# BAM files
# ---------------------------------------------------------------------------


@contextmanager
def _alignments(path, threads=1):
    """A BAM (or SAM) file, opened with pysam, which htslib's messages do not
    interrupt (a BAM file without index is read from its start): a file that
    pysam cannot read is a MetacoolError naming it."""
    save = pysam.set_verbosity(0)
    try:
        try:
            bam = pysam.AlignmentFile(str(path), "r", check_sq=False, threads=threads)
        except (OSError, ValueError) as error:
            raise mtmc.MetacoolError(f"Cannot read the BAM file {path}: {error}") from error
        with bam:
            yield bam
    finally:
        pysam.set_verbosity(save)


def bam_layout(path, threads=1):
    """BOTH_MATES or ONE_MATE: whether a BAM (or SAM) file holds the two mates of
    its reads (paired-end alignments), or one mate of each, from its first
    LAYOUT_PROBE primary alignments: a read flagged first or second of its
    pair is a mate, an unpaired read one mate aligned on its own.

    Raises
    ------
    metator.metacool.MetacoolError
        If the file cannot be read, or holds no read.
    """
    mates = set()
    with _alignments(path, threads) as bam:
        try:
            for i, read in enumerate(bam):
                if i == LAYOUT_PROBE:
                    break
                if read.is_secondary or read.is_supplementary:
                    continue
                mates.add(2 if read.is_read2 else 1 if read.is_read1 else 0)
                if {1, 2} <= mates:
                    return BOTH_MATES
        except (OSError, ValueError) as error:
            raise mtmc.MetacoolError(f"Cannot read the BAM file {path}: {error}") from error
    if not mates:
        raise mtmc.MetacoolError(f"{path} holds no read.")
    return ONE_MATE


def _references(path):
    """The reference sequences of a BAM file, as its header lists them: (name,
    length) pairs."""
    with _alignments(path) as bam:
        return list(zip(bam.references, bam.lengths, strict=True))


def check_references(path, references, lengths, what="reference sequences"):
    """Refuse the pairs or the alignments of a file that is not of the assembly:
    sequences (references: (name, length) pairs) that the assembly gives
    another length (lengths: by contig name), or none of which is a contig of
    the assembly. The sequences that are not contigs of the assembly (a
    spike-in, a host genome) are left out, with a warning: so are the pairs
    with a read on them.

    Parameters
    ----------
    path : str
        The BAM or pairs file, for the messages.
    references : list of tuple
        Its sequences, as (name, length).
    lengths : dict
        The length of each contig of the assembly, by name.
    what : str
        What the sequences are, in the messages.

    Raises
    ------
    metator.metacool.MetacoolError
        If a sequence has another length in the assembly, or none is a contig
        of the assembly.
    """
    other = [
        f"{name!r} ({length} bp; {lengths[name]} bp in the assembly)"
        for name, length in references
        if name in lengths and lengths[name] != length
    ]
    if other:
        raise mtmc.MetacoolError(
            f"{path} is not of the contigs of the assembly: {len(other)} of its {what} have another length in the "
            f"assembly: {mtmc._shown(other, show=str)}."
        )
    unknown = [name for name, _ in references if name not in lengths]
    if references and len(unknown) == len(references):
        raise mtmc.MetacoolError(
            f"{path} is not of the contigs of the assembly: none of its {len(references)} {what} is a contig of the "
            f"assembly: {mtmc._shown(unknown)}."
        )
    if unknown:
        logger.warning(
            f"{len(unknown)} of the {what} of {path} are not contigs of the assembly: the pairs with a read on them are "
            f"left out: {mtmc._shown(unknown)}."
        )


class _PairsWriter:
    """The pairs file of a library, made from its alignments: a pair per line,
    with the 7 columns of a pairs file, read name, the contig and position of
    each read, and their strands (1-based leftmost base of each alignment, the
    contig first in the order of the assembly, as an upper triangle). The read
    name is ".": the map does not keep it.

    The pairs with a read on a sequence that is not a contig of the assembly
    are left out (counted as other).
    """

    def __init__(self, path, lengths):
        self.path = path
        self.order = {name: i for i, name in enumerate(lengths)}
        self.pairs = self.other = 0

    def __enter__(self):
        self.handle = open(self.path, "w", buffering=1 << 20)
        # No #chromsize line: metacooler takes the contigs from the bins of the
        # map, those of the assembly.
        self.handle.write("## pairs format v1.0\n#shape: upper triangle\n")
        self.handle.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        self.pairs = self.other = 0
        return self

    def __exit__(self, *exc):
        self.handle.close()

    def write(self, read1, names1, read2, names2):
        """Write the pair of two mates (pysam.AlignedSegment), whose contigs are
        named by names1 and names2 (the reference names of their files, by id)."""
        contig1, contig2 = names1[read1.reference_id], names2[read2.reference_id]
        rank1, rank2 = self.order.get(contig1), self.order.get(contig2)
        if rank1 is None or rank2 is None:
            self.other += 1
            return
        pos1, pos2 = read1.reference_start + 1, read2.reference_start + 1
        strand1, strand2 = "-" if read1.is_reverse else "+", "-" if read2.is_reverse else "+"
        if (rank1, pos1) > (rank2, pos2):
            contig1, pos1, strand1, contig2, pos2, strand2 = contig2, pos2, strand2, contig1, pos1, strand1
        self.handle.write(f".\t{contig1}\t{pos1}\t{contig2}\t{pos2}\t{strand1}\t{strand2}\n")
        self.pairs += 1


def _passes(read, min_quality):
    """Whether a mate is aligned, with a mapping quality of at least min_quality."""
    return not read.is_unmapped and read.mapping_quality >= min_quality


class _GroupingCheck:
    """The names of the latest reads with one mate in a pass over a BAM file of
    both mates (GROUPING_WINDOW of them): a read that comes back among them is
    split, and the file not grouped by read name (_NotGrouped)."""

    def __init__(self, size=GROUPING_WINDOW):
        self.size = size
        self.names = OrderedDict()

    def lone(self, name):
        if name in self.names:
            raise _NotGrouped(name)
        self.names[name] = None
        if len(self.names) > self.size:
            self.names.popitem(last=False)


def _both_mates_pass(path, writer, min_quality, threads, check=None):
    """Write the pairs of a BAM file of both mates, its alignments grouped by read
    name (as bwa writes them, or sorted by name): the primary alignments of the
    first and the second mate of each read. With check (a _GroupingCheck),
    _NotGrouped is raised as soon as a read is found split in two groups.

    Returns
    -------
    dict
        The number of reads ("reads"), of reads with one mate in the file
        ("lone"), and of reads with both mates aligned with a mapping
        quality of at least min_quality on contigs of the assembly
        ("pairs").
    """
    counts = {"reads": 0, "lone": 0}
    with _alignments(path, threads) as bam:
        names = list(bam.references)
        current, mates = None, [None, None]

        def finish():
            counts["reads"] += 1
            first, second = mates
            if first is None or second is None:
                counts["lone"] += 1
                if check is not None:
                    check.lone(current)
            elif _passes(first, min_quality) and _passes(second, min_quality):
                writer.write(first, names, second, names)

        try:
            for read in bam:
                if read.is_secondary or read.is_supplementary or not (read.is_read1 or read.is_read2):
                    continue
                if read.query_name != current:
                    if current is not None:
                        finish()
                    current, mates = read.query_name, [None, None]
                index = 1 if read.is_read2 else 0
                if mates[index] is None:
                    mates[index] = read
            if current is not None:
                finish()
        except (OSError, ValueError) as error:
            raise mtmc.MetacoolError(f"Cannot read the BAM file {path}: {error}") from error
    counts["pairs"] = writer.pairs
    return counts


def _both_mates_pairs(path, out, lengths, min_quality, tmp_dir, threads):
    """Write the pairs file of a BAM file of both mates (see _both_mates_pass):
    read as it is, its reads grouped by read name, unless its header says it is
    sorted by coordinate, or a read turns out to be split: then grouped by name
    first (samtools collate, in pysam), into tmp_dir."""
    with _alignments(path) as bam:
        sorted_by = bam.header.to_dict().get("HD", {}).get("SO")
    if sorted_by != "coordinate":
        try:
            with _PairsWriter(out, lengths) as writer:
                return _both_mates_pass(path, writer, min_quality, threads, _GroupingCheck()), writer.other
        except _NotGrouped as split:
            logger.info(f"{path} is not grouped by read name (read {split} comes back after other reads): grouping it.")
    collated = tempfile.mkdtemp(prefix="collate_", dir=tmp_dir)
    try:
        grouped = os.path.join(collated, "grouped.bam")
        logger.info(f"Grouping the alignments of {path} by read name (samtools collate), in {collated}.")
        try:
            pysam.collate("-@", str(threads), "-o", grouped, str(path), os.path.join(collated, "part"))
        except pysam.SamtoolsError as error:
            raise mtmc.MetacoolError(f"Cannot group the alignments of {path} by read name: {error}") from error
        with _PairsWriter(out, lengths) as writer:
            return _both_mates_pass(grouped, writer, min_quality, threads), writer.other
    finally:
        shutil.rmtree(collated)


def natural_key(name):
    """The sort key of a read name in the order of ``samtools sort -n``: its
    characters one by one, but its runs of digits found at the same place in
    two names compared by value (strnum_cmp of samtools: r9 before r10), and
    their leading zeros ignored (samtools 1.17 and later). Adapted from
    hicstuff (BSD-3-Clause), which matches the mates of two BAM files so.

    Returns
    -------
    tuple
        The text between runs of digits, each ended with "0" when digits follow
        (samtools compares a digit and another character as characters, and
        every digit compares alike with a character that is no digit), and the
        value of each run of digits.

    Examples
    --------
    >>> sorted(["r10", "r9", "r1"], key=natural_key)
    ['r1', 'r9', 'r10']
    >>> natural_key("r010") == natural_key("r10")
    True
    """
    parts = _DIGIT_RUNS.split(name)
    key = []
    for i in range(1, len(parts), 2):
        key.append(parts[i - 1] + "0")
        key.append(int(parts[i]))
    key.append(parts[-1])
    return tuple(key)


#: How samtools sorts read names in each order.
_ORDERS = {natural_key: "samtools sort -n", str: "samtools sort -N"}


def _declared_order(path):
    """The read name order that a BAM file sorted by name says it follows:
    lexicographical (str: ``samtools sort -N``, from samtools 1.19, records it
    in its header) or, by default, natural_key (``samtools sort -n``)."""
    with _alignments(path) as bam:
        sub_sort = bam.header.to_dict().get("HD", {}).get("SS", "")
    return str if str(sub_sort).endswith("lexicographical") else natural_key


class _Mate:
    """The primary alignments of the BAM file of one mate, read by read, in the
    order of the file: the first primary alignment of each read, its name
    without a mate suffix (/1 or /2), and the sort key of that name (see
    order_key)."""

    def __init__(self, bam, path, name_key):
        self.records = iter(bam)
        self.path = path
        self.names = list(bam.references)
        self.name_key = name_key
        self.key_suffix = ""
        self.last_key = None
        self.reads = 0
        self.read = self.name = self.key = None
        self.advance()

    def advance(self):
        """Move to the next read: read and name are None at the end."""
        previous = self.name
        for read in self.records:
            if read.is_secondary or read.is_supplementary:
                continue
            name = _MATE_SUFFIX.sub("", read.query_name)
            if name == previous:
                # Another primary alignment of the same read: the first is kept.
                continue
            self.read, self.name, self.key = read, name, None
            self.reads += 1
            return
        self.read = self.name = self.key = None

    def order_key(self):
        """The sort key of the name of the current read, which must not come
        before the last one computed in this file: _NameOrderError otherwise."""
        if self.key is None:
            self.key = self.name_key(f"{self.name}{self.key_suffix}")
            if self.last_key is not None and self.key < self.last_key:
                raise _NameOrderError(self.path)
            self.last_key = self.key
        return self.key


def _one_mate_pass(path1, path2, writer, name_key, min_quality, threads):
    """Write the pairs of the BAM files of the two mates of a library, both
    sorted by read name in the order of name_key, walking both files together:
    the mates are matched by name, without a mate suffix. When the names
    differ, the read whose name comes first has no mate in the other file. A
    name that comes before the previous one of its file raises
    _NameOrderError: had a read been passed while its mate was still ahead,
    its mate would break the order. Names equal in that order but different
    (numbers padded with zeros: r010 and r10) are matched by name. Adapted
    from hicstuff's bam2pairs (BSD-3-Clause).

    Returns
    -------
    dict
        The reads of each file ("reads1", "reads2"), the reads with both mates
        ("matched") and the pairs written ("pairs").
    """
    matched = 0
    with _alignments(path1, threads) as bam1, _alignments(path2, threads) as bam2:
        try:
            first, second = _Mate(bam1, path1, name_key), _Mate(bam2, path2, name_key)
            # Names compared as both files sorted them: with the /1 suffix when
            # both have suffixes (r1.5/1 comes before r1/1, but r1 before r1.5).
            if all(mate.read is not None and mate.read.query_name != mate.name for mate in (first, second)):
                first.key_suffix = second.key_suffix = "/1"

            def pair(read1, read2):
                nonlocal matched
                matched += 1
                if _passes(read1, min_quality) and _passes(read2, min_quality):
                    writer.write(read1, first.names, read2, second.names)

            while first.read is not None or second.read is not None:
                if first.name == second.name:
                    pair(first.read, second.read)
                    first.advance()
                    second.advance()
                    continue
                key1 = None if first.read is None else first.order_key()
                key2 = None if second.read is None else second.order_key()
                if key1 is not None and key1 == key2:
                    tied = []
                    for mate in (first, second):
                        reads = {}
                        while mate.read is not None and mate.order_key() == key1:
                            reads.setdefault(mate.name, []).append(mate.read)
                            mate.advance()
                        tied.append(reads)
                    for name, reads1 in tied[0].items():
                        for read1, read2 in zip(reads1, tied[1].pop(name, []), strict=False):
                            pair(read1, read2)
                    continue
                # The read that comes first has no mate in the other file.
                if key2 is None or (key1 is not None and key1 < key2):
                    first.advance()
                else:
                    second.advance()
        except (OSError, ValueError) as error:
            raise mtmc.MetacoolError(f"Cannot read the BAM files {path1} and {path2}: {error}") from error
    return {"reads1": first.reads, "reads2": second.reads, "matched": matched, "pairs": writer.pairs}


def _one_mate_pairs(path1, path2, out, lengths, min_quality, tmp_dir, threads):
    """Write the pairs file of the BAM files of the two mates of a library (see
    _one_mate_pass), in the read name order that they say they are sorted in
    (natural, or lexicographical), else in the other one, as some tools sort
    names lexicographically without saying so. Files in neither order are
    sorted by name first (samtools sort -n, in pysam), into tmp_dir."""
    orders = {_declared_order(path1), _declared_order(path2)}
    declared = orders.pop() if len(orders) == 1 else natural_key
    unsorted = path1
    for name_key in (declared, str if declared is natural_key else natural_key):
        try:
            with _PairsWriter(out, lengths) as writer:
                return _one_mate_pass(path1, path2, writer, name_key, min_quality, threads), writer.other
        except _NameOrderError as error:
            unsorted = str(error)
            logger.debug(f"{unsorted} is not sorted by read name as {_ORDERS[name_key]} sorts them.")
    sort_dir = tempfile.mkdtemp(prefix="sort_", dir=tmp_dir)
    try:
        logger.info(
            f"{unsorted} is not sorted by read name: sorting the BAM files of both mates by read name (samtools sort -n), "
            f"in {sort_dir}."
        )
        sorted_bams = [os.path.join(sort_dir, "mate1.bam"), os.path.join(sort_dir, "mate2.bam")]
        for i, (bam, sorted_bam) in enumerate(zip((path1, path2), sorted_bams, strict=True)):
            try:
                pysam.sort("-n", "-@", str(threads), "-T", os.path.join(sort_dir, f"part{i}"), "-o", sorted_bam, str(bam))
            except pysam.SamtoolsError as error:
                raise mtmc.MetacoolError(f"Cannot sort the alignments of {bam} by read name: {error}") from error
        try:
            with _PairsWriter(out, lengths) as writer:
                counts = _one_mate_pass(*sorted_bams, writer, natural_key, min_quality, threads)
                return counts, writer.other
        except _NameOrderError:
            raise mtmc.MetacoolError(
                f"Cannot match the mates of {path1} and {path2}: sorted by read name with samtools sort -n (samtools "
                f"{pysam.__samtools_version__}), their read names are still not in the order expected."
            ) from None
    finally:
        shutil.rmtree(sort_dir)


def bam_pairs(library, out, lengths, layout, min_quality=30, tmp_dir=".", threads=1):
    """Write the pairs file of a library of BAM files: the pairs of the reads
    whose two mates are aligned, each with a mapping quality of at least
    min_quality, on contigs of the assembly (their primary alignments; the
    secondary and supplementary ones are left out).

    Parameters
    ----------
    library : tuple of str
        The BAM file of both mates, or the BAM files of the two mates.
    out : str
        The pairs file to write.
    lengths : dict
        The length of each contig of the assembly, by name, in its order.
    layout : str
        BOTH_MATES or ONE_MATE.
    min_quality : int
        The mapping quality that both mates of a pair must reach.
    tmp_dir : str
        Where the files are grouped or sorted by read name, when they need it.
    threads : int
        Threads decompressing the BAM files, and grouping or sorting them.

    Returns
    -------
    int
        The number of pairs written.

    Raises
    ------
    metator.metacool.MetacoolError
        If a file cannot be read, holds no read with both mates (BAM files of
        the mates of other libraries), or the mates cannot be matched.
    """
    if layout == BOTH_MATES:
        (path,) = library
        counts, other = _both_mates_pairs(path, out, lengths, min_quality, tmp_dir, threads)
        with_both, what = counts["reads"] - counts["lone"], path
        if counts["lone"]:
            logger.info(f"{counts['lone']} of the {counts['reads']} reads of {path} have one mate in it: no pair.")
    else:
        path1, path2 = library
        counts, other = _one_mate_pairs(path1, path2, out, lengths, min_quality, tmp_dir, threads)
        with_both, what = counts["matched"], f"{path1} and {path2}"
        alone = counts["reads1"] + counts["reads2"] - 2 * counts["matched"]
        if alone:
            logger.info(f"{alone} reads have one mate only in {what} (an unaligned mate that the aligner left out): no pair.")
    if not with_both:
        mates = "the two mates of a read" if layout == BOTH_MATES else "the R1 and R2 files of one library"
        raise mtmc.MetacoolError(f"{what} hold no read with both mates: they are not {mates}.")
    if other:
        logger.info(f"{other} pairs of {what} have a read on a sequence that is not a contig of the assembly: left out.")
    share = 100 * counts["pairs"] / with_both
    logger.info(
        f"{what}: {counts['pairs']} pairs of {with_both} reads with both mates ({share:.2f}%), both aligned with a "
        f"mapping quality of at least {min_quality} on contigs of the assembly."
    )
    return counts["pairs"]


# ---------------------------------------------------------------------------
# Pairs files
# ---------------------------------------------------------------------------


def pairs_chromsizes(path):
    """The sequences that the header of a pairs file lists (#chromsize lines), as
    (name, length) pairs: none for a file without them, or compressed with
    zstd or LZ4, whose header is not read here (metacooler reads the file)."""
    with open(path, "rb") as handle:
        head = handle.read(4)
    if head.startswith(_UNREAD_COMPRESSIONS):
        return []
    sizes = []
    if mio.file_compression(path) == "zip":
        return sizes
    try:
        with mio.read_compressed(path) as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                if line.startswith("#chromsize:"):
                    fields = line[len("#chromsize:") :].split()
                    if len(fields) == 2 and fields[1].isdigit():
                        sizes.append((fields[0], int(fields[1])))
    except (UnicodeDecodeError, OSError, EOFError):
        # No text: metacooler says what it is not.
        return []
    return sizes


# ---------------------------------------------------------------------------
# The map of the contigs
# ---------------------------------------------------------------------------


def map_name(prefix, assembly):
    """The name of the map that a run makes from pairs or BAM files:
    <prefix>.metacool, or <assembly>.metacool without a prefix, the name of the
    FASTA file of the assembly without its extensions (.gz, .fa, .fasta,
    .fna...)."""
    return f"{prefix or assembly_name(assembly)}.metacool"


#: The extensions of a FASTA file, which assembly_name drops.
_FASTA_EXTENSIONS = (".fasta", ".fas", ".fa", ".fna", ".ffn", ".fsa", ".seq")


def assembly_name(assembly):
    """The name of an assembly: that of its FASTA file, without the extensions of
    a FASTA file and of its compression (assembly from assembly.fa.gz)."""
    name = os.path.basename(str(assembly))
    for extension in (".gz", ".bz2", ".zip"):
        if name.endswith(extension) and len(name) > len(extension):
            name = name[: -len(extension)]
            break
    for extension in _FASTA_EXTENSIONS:
        if name.lower().endswith(extension) and len(name) > len(extension):
            return name[: -len(extension)]
    return name


def make_map(inputs, contig_data, outfile, assembly=None, min_quality=30, threads=1, tmp_dir=".", command=None):
    """Write the map of the contigs at outfile, from the pairs files or the BAM
    files of a run (inputs, see resolve), as ``metacooler create pairs`` makes
    it: one bin per contig of the assembly (contig_data, by name, in its
    order: see metator.network.create_contig_data), every contig in
    ``__unbinned__``, the longest first, and the contacts of the pairs of
    every library, summed, with the contact statistics that metacooler counts
    on the reads (the QC block, which metator mge reads).

    The BAM files give the pairs of each library (see bam_pairs), written in
    tmp_dir. The headers of the BAM files, and the #chromsize lines of the
    pairs files, are checked against the assembly first (see
    check_references).

    Parameters
    ----------
    inputs : Inputs
        Pairs files or BAM files.
    contig_data : dict
        The contigs of the assembly, by name, with their "length".
    outfile : str
        The map to write. An earlier one there is replaced once the new map is
        written.
    assembly : str, optional
        The FASTA file of the assembly, whose name the map stores (see
        assembly_name).
    min_quality : int
        The mapping quality that both mates of a pair of a BAM file must reach.
    threads : int
        Threads reading the BAM files and binning the pairs.
    tmp_dir : str
        Where the pairs of the BAM files are written.
    command : str, optional
        The command line, stored in the history of the map.

    Raises
    ------
    metator.metacool.MetacoolError
        If an input is not of the assembly, cannot be read, or no pair of the
        inputs is on the contigs of the assembly. outfile is then left as it
        was: the map is written next to it, as a hidden file, and renamed onto
        it once it holds contacts.
    """
    import metacooler

    lengths = {name: data["length"] for name, data in contig_data.items() if data["length"] > 0}
    mtmc.check_ascii_names(list(lengths), [])
    if inputs.kind == PAIRS:
        for path in inputs.files:
            check_references(path, pairs_chromsizes(path), lengths, "sequences (#chromsize)")
        pairs = list(inputs.files)
    else:
        for path in inputs.files:
            references = _references(path)
            if not references:
                raise mtmc.MetacoolError(f"{path} lists no reference sequence in its header: its reads are not aligned.")
            check_references(path, references, lengths)
        pairs = []
        for i, library in enumerate(inputs.libraries):
            pairs.append(os.path.join(tmp_dir, f"library_{i + 1}.pairs"))
            bam_pairs(library, pairs[-1], lengths, inputs.layout, min_quality, tmp_dir, threads)
    logger.info(f"Writing {outfile}: the contact map of the {len(lengths)} contigs of the assembly, in no cluster yet.")
    # Written next to outfile, under a hidden name, and renamed onto it once it
    # holds contacts: an earlier map there is left as it was until then.
    directory, name = os.path.split(os.path.abspath(outfile))
    staged = os.path.join(directory, f".{name}.{os.getpid()}.{os.urandom(4).hex()}.tmp")
    try:
        with mtmc._metacooler_logging(), mtmc._as_metacool_errors():
            metacooler.from_pairs(
                staged,
                bins=pd.DataFrame({"chrom": list(lengths), "start": 0, "end": list(lengths.values())}),
                pairs=pairs,
                # No clustering yet: every contig in __unbinned__, the longest
                # first.
                contig_genome_map=None,
                assembly_name=assembly_name(assembly) if assembly else None,
                metadata={"commands": [command]} if command else None,
                # The QC block, counted on the reads, which only the pairs give.
                qc=True,
                # No P(s) table next to the map: a run writes the map only.
                ps=False,
                nproc=threads,
                keep_partial=False,
            )
            contacts = mtmc.open_map(staged).info.get("nnz", 0)
        if not contacts:
            raise mtmc.MetacoolError(
                f"No pair from {inputs.describe()} is on the contigs of the assembly: give the pairs or the alignments "
                "of the reads of this assembly."
            )
        os.replace(staged, outfile)
    finally:
        if os.path.lexists(staged):
            os.remove(staged)
    logger.info(f"{outfile} is written: the contact map of the contigs, from {inputs.describe()}.")
