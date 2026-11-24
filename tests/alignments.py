"""BAM files of pairs, as aligners write them, for the tests of metator.inputs:
the alignments of the two mates of each pair, which metaTOR reads back as those
pairs (each mate at the 1-based 5' end of its alignment: its first aligned base
on the forward strand, its last on the reverse strand)."""

import os
import tempfile

import pandas as pd
import pysam

#: The length of the alignment of a mate (shorter on the last bases of a
#: contig).
READ_LENGTH = 50

#: The columns of a pairs table.
PAIRS_COLUMNS = ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]


def read_pairs(path):
    """The pairs of a pairs file, as a table of PAIRS_COLUMNS."""
    return pd.read_csv(
        path, sep="\t", comment="#", header=None, usecols=range(7), names=PAIRS_COLUMNS, dtype={0: str, 1: str, 3: str}
    )


def contacts(pairs):
    """The pairs of a table as a sorted list of contacts, each mate (contig,
    position, strand), the two mates in a set order: what a map counts."""
    found = []
    for row in pairs[PAIRS_COLUMNS[1:]].itertuples(index=False):
        first, second = (row[0], int(row[1]), row[4]), (row[2], int(row[3]), row[5])
        found.append(tuple(sorted([first, second])))
    return sorted(found)


def header(lengths, sort_order=None, sub_sort=None):
    """The header of a BAM file of the contigs of lengths (by name), sorted as
    sort_order says (SO), and sub_sort (SS)."""
    hd = {"VN": "1.6"}
    if sort_order:
        hd["SO"] = sort_order
    if sub_sort:
        hd["SS"] = sub_sort
    return pysam.AlignmentHeader.from_dict({"HD": hd, "SQ": [{"SN": name, "LN": int(n)} for name, n in lengths.items()]})


def segment(bam_header, name, contig=None, pos=None, strand="+", flag=0, mapq=60):
    """The alignment of a mate whose 5' end is at the 1-based position pos of
    contig, on strand: from pos on the forward strand, up to pos on the reverse
    strand, READ_LENGTH long (shorter at either end of the contig); unmapped
    without a contig."""
    read = pysam.AlignedSegment(bam_header)
    read.query_name = name
    if contig is None:
        read.flag = flag | 0x4
        read.query_sequence = "N" * READ_LENGTH
        return read
    pos = int(pos)
    if strand == "-":
        start = max(0, pos - READ_LENGTH)
        length = pos - start
    else:
        start = pos - 1
        length = min(READ_LENGTH, bam_header.get_reference_length(contig) - start)
    read.flag = flag | (0x10 if strand == "-" else 0)
    read.reference_id = bam_header.get_tid(contig)
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigartuples = [(0, length)]
    read.query_sequence = "A" * length
    return read


def _mates(bam_header, pairs, paired, mapq=60):
    """The alignments of the two mates of each pair, in the order of the pairs:
    flagged as the first and second mates of paired-end alignments, or
    unpaired (one mate at a time)."""
    for row in pairs[PAIRS_COLUMNS].itertuples(index=False):
        name, chr1, pos1, chr2, pos2, strand1, strand2 = row
        one = segment(bam_header, name, chr1, pos1, strand1, flag=0x41 if paired else 0, mapq=mapq)
        two = segment(bam_header, name, chr2, pos2, strand2, flag=0x81 if paired else 0, mapq=mapq)
        if paired:
            for read, mate in ((one, two), (two, one)):
                read.next_reference_id, read.next_reference_start = mate.reference_id, mate.reference_start
                read.mate_is_reverse = mate.is_reverse
        yield one, two


def _sorted(path, order):
    """Sort a BAM file in place: by read name ("name", samtools sort -n, or
    "lexicographic", -N) or by coordinate ("coordinate"). "aligner" leaves it
    as it is."""
    if order == "aligner":
        return path
    flag = {"name": ["-n"], "lexicographic": ["-N"], "coordinate": []}[order]
    directory = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(path)))
    sorted_path = os.path.join(directory, "sorted.bam")
    pysam.sort(*flag, "-T", os.path.join(directory, "part"), "-o", sorted_path, str(path))
    os.replace(sorted_path, path)
    os.rmdir(directory)
    return path


def both_mates(path, pairs, lengths, order="aligner", mapq=60, extra=()):
    """Write a BAM file of both mates of each pair (paired-end alignments, as bwa
    mem -5SP writes them): in the order of the pairs, the two mates of a read in
    a row, with no sort order in the header ("aligner"), or sorted by read name
    or by coordinate (see _sorted). extra: more alignments (see segment), each
    after the mates of its read."""
    bam_header = header(lengths)
    extra_of = {}
    for read in extra:
        extra_of.setdefault(read.query_name, []).append(read)
    with pysam.AlignmentFile(str(path), "wb", header=bam_header) as bam:
        for one, two in _mates(bam_header, pairs, True, mapq):
            bam.write(one)
            bam.write(two)
            for read in extra_of.pop(one.query_name, []):
                bam.write(read)
        for reads in extra_of.values():
            for read in reads:
                bam.write(read)
    return _sorted(path, order)


def one_mate(paths, pairs, lengths, order="name", mapq=60):
    """Write the two BAM files of the first and the second mates of each pair
    (unpaired alignments, as bowtie2 writes them, one mate at a time), each
    sorted by read name (or as order says: see _sorted)."""
    bam_header = header(lengths)
    with pysam.AlignmentFile(str(paths[0]), "wb", header=bam_header) as first:
        with pysam.AlignmentFile(str(paths[1]), "wb", header=bam_header) as second:
            for one, two in _mates(bam_header, pairs, False, mapq):
                first.write(one)
                second.write(two)
    for path in paths:
        _sorted(path, order)
    return paths
