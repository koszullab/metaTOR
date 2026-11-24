"""The inputs of network and pipeline other than a map: pairs files and BAM files,
which a run makes the map of its contigs from (metator.inputs)."""

import bz2
import gzip
import os
import re
import shlex
import shutil
from pathlib import Path

import cooler
import metacooler
import pandas as pd
import pysam
import pytest

import metator.inputs as mti
import metator.metacool as mtmc
import metator.network as mtn
from tests import alignments
from tests.maps import contig_map

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"


@pytest.fixture(scope="module")
def community(metagenome, metagenome_files):
    """The pairs of the mock community, the lengths of its contigs, and its
    contig data (see metator.network.create_contig_data)."""
    return {
        "pairs": metagenome.pairs,
        "lengths": {contig.name: contig.length for contig in metagenome.contigs},
        "contig_data": mtn.create_contig_data(str(metagenome_files["assembly"])),
    }


def _split_pairs(source, paths, compress=gzip.open):
    """Write the pairs of a pairs file into several, each with its header, as
    libraries of one sample: their contacts sum to those of source."""
    lines = Path(source).read_text().splitlines(keepends=True)
    head = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    size = -(-len(body) // len(paths))
    for i, path in enumerate(paths):
        with compress(path, "wt") as handle:
            handle.writelines(head + body[i * size : (i + 1) * size])
    return paths


def _same_map(found, expected):
    """Two maps of the same contacts: the same bins, pixels, contigs and genomes
    (their contact statistics too), and QC block."""
    clr, ref = cooler.Cooler(str(found)), cooler.Cooler(str(expected))
    for table in ["bins", "pixels", "chroms"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(found)).genomes, metacooler.MetaCooler(str(expected)).genomes)
    assert clr.info["metadata"]["qc"] == ref.info["metadata"]["qc"]


# ---------------------------------------------------------------------------
# The arguments
# ---------------------------------------------------------------------------


def test_expand_paths_lists_and_patterns(tmp_path, monkeypatch):
    """An argument is a path, a comma-separated list or glob patterns, whose
    matches are sorted by name; a path that exists is that file, whatever its
    name holds."""
    monkeypatch.chdir(tmp_path)
    for name in ["L2.pairs", "L10.pairs", "L1.pairs", "a,b.pairs", "x[1].pairs"]:
        Path(name).write_text("")
    assert mti.expand("L1.pairs") == ["L1.pairs"]
    assert mti.expand("L2.pairs,L1.pairs") == ["L2.pairs", "L1.pairs"]
    assert mti.expand("L*.pairs") == ["L1.pairs", "L10.pairs", "L2.pairs"]
    assert mti.expand("L?.pairs,x*") == ["L1.pairs", "L2.pairs", "x[1].pairs"]
    assert mti.expand("a,b.pairs") == ["a,b.pairs"]
    assert mti.expand("x[1].pairs") == ["x[1].pairs"]
    # A path that is not there stays as it is: file_kind refuses it.
    assert mti.expand("missing.pairs") == ["missing.pairs"]


@pytest.mark.parametrize(
    "argument, message",
    [
        ("none*.bam", "No file matches the input pattern 'none*.bam'."),
        ("L1.pairs,,L2.pairs", "The list of inputs 'L1.pairs,,L2.pairs' holds an empty path."),
        ("L1.pairs,", "The list of inputs 'L1.pairs,' holds an empty path."),
    ],
)
def test_expand_refusals(argument, message, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path("L1.pairs").write_text("")
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        mti.expand(argument)


def _write_kinds(tmp_path, lengths):
    """One file of each kind that file_kind tells apart, by name."""
    pairs = "## pairs format v1.0\n.\tc1\t5\tc2\t9\t+\t-\n"
    files = {
        "map": contig_map(
            TESTS_DATA / "alignment.pairs", {"NODE_522": {"length": 22786}}, tmp_path / "x.metacool", columns=False
        ),
        "plain pairs": tmp_path / "plain.pairs",
        "gzip pairs": tmp_path / "gz.pairs",
        "bzip2 pairs": tmp_path / "bz2.pairs",
        "pairs without header": tmp_path / "bare.txt",
        "BAM": tmp_path / "aligned.data",
        "SAM": tmp_path / "aligned.sam",
    }
    files["plain pairs"].write_text(pairs)
    with gzip.open(files["gzip pairs"], "wt") as handle:
        handle.write(pairs)
    with bz2.open(files["bzip2 pairs"], "wt") as handle:
        handle.write(pairs)
    files["pairs without header"].write_text(".\tc1\t5\tc2\t9\t+\t-\n")
    bam_header = alignments.header(lengths)
    with pysam.AlignmentFile(str(files["BAM"]), "wb", header=bam_header) as bam:
        bam.write(alignments.segment(bam_header, "r1", "c1", 5))
    with pysam.AlignmentFile(str(files["SAM"]), "w", header=bam_header) as sam:
        sam.write(alignments.segment(bam_header, "r1", "c1", 5))
    return files


def test_file_kind(tmp_path):
    """What a file holds is read from its first bytes, whatever its name: a map
    (and file::group, a map inside a file), alignments (BAM, SAM), or pairs,
    plain or compressed."""
    files = _write_kinds(tmp_path, {"c1": 100, "c2": 100})
    kinds = {name: mti.file_kind(path) for name, path in files.items()}
    assert kinds == {
        "map": mti.MAP,
        "plain pairs": mti.PAIRS,
        "gzip pairs": mti.PAIRS,
        "bzip2 pairs": mti.PAIRS,
        "pairs without header": mti.PAIRS,
        "BAM": mti.BAM,
        "SAM": mti.BAM,
    }
    assert mti.file_kind(f"{tmp_path / 'maps.mmetacool'}::entries/rep1") == mti.MAP


@pytest.mark.parametrize("kind", ["reads", "gzipped reads", "assembly", "CRAM", "missing", "directory"])
def test_file_kind_refusals(kind, tmp_path):
    """Reads, an assembly, CRAM alignments, a file that is not there and a
    directory are refused, with what to give instead."""
    path = tmp_path / "input"
    message = {
        "reads": f"{path} holds reads (FASTQ): metaTOR aligns no reads.",
        "gzipped reads": f"{path} holds reads (FASTQ): metaTOR aligns no reads.",
        "assembly": f"{path} is a FASTA file: give the assembly with -a, and the pairs or BAM files of the run.",
        "CRAM": f"{path} is a CRAM file, which metaTOR does not read: convert it to BAM",
        "missing": f"Cannot read the input {path}: no such file.",
        "directory": f"Cannot read the input {path}: it is a directory.",
    }[kind]
    if kind == "reads":
        path.write_text("@r1\nACGT\n+\nIIII\n")
    elif kind == "gzipped reads":
        with gzip.open(path, "wt") as handle:
            handle.write("@r1\nACGT\n+\nIIII\n")
    elif kind == "assembly":
        path.write_text(">c1\nACGT\n")
    elif kind == "CRAM":
        path.write_bytes(b"CRAM\x03\x00" + bytes(20))
    elif kind == "directory":
        path.mkdir()
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}"):
        mti.file_kind(path)


def _bams(tmp_path, names, layout, lengths=None):
    """Small BAM files of a pair each, named as names, of both mates or of one."""
    lengths = lengths or {"c1": 1000, "c2": 1000}
    pairs = pd.DataFrame([[f"r{i}", "c1", 10 + i, "c2", 20 + i, "+", "-"] for i in range(3)], columns=alignments.PAIRS_COLUMNS)
    paths = []
    if layout == mti.BOTH_MATES:
        for name in names:
            paths.append(alignments.both_mates(tmp_path / name, pairs, lengths))
    else:
        for first, second in zip(names[::2], names[1::2], strict=True):
            paths.extend(alignments.one_mate([tmp_path / first, tmp_path / second], pairs, lengths))
    return [str(path) for path in paths]


def test_resolve_pairs_and_maps(tmp_path, monkeypatch):
    """Pairs files are each a library, in the order of the arguments; a map is
    the one input."""
    monkeypatch.chdir(tmp_path)
    for name in ["L1.pairs", "L2.pairs", "L3.pairs"]:
        Path(name).write_text("## pairs format v1.0\n")
    inputs = mti.resolve(["L3.pairs", "L[12].pairs"])
    assert (inputs.kind, inputs.files, inputs.libraries) == (
        mti.PAIRS,
        ["L3.pairs", "L1.pairs", "L2.pairs"],
        [("L3.pairs",), ("L1.pairs",), ("L2.pairs",)],
    )
    assert inputs.describe() == "3 pairs files" and inputs.map is None
    assert mti.resolve(["L1.pairs"]).describe() == "a pairs file"
    metacool = contig_map(TESTS_DATA / "alignment.pairs", {"NODE_522": {"length": 22786}}, tmp_path / "x.metacool", columns=False)
    inputs = mti.resolve([str(metacool)])
    assert (inputs.kind, inputs.map, inputs.describe()) == (mti.MAP, str(metacool), "a map")
    with pytest.raises(mtmc.MetacoolError, match=r"^A run needs its inputs: the map of its contigs, or its pairs files"):
        mti.resolve([])


def test_resolve_bam_files(tmp_path, monkeypatch):
    """BAM files of both mates are each a library; BAM files of one mate come
    as the R1 files then the R2 files: two arguments (lists, or patterns), or
    as many R1 arguments as R2 arguments, as the shell expands patterns."""
    monkeypatch.chdir(tmp_path)
    both = _bams(tmp_path, ["A.bam", "B.bam"], mti.BOTH_MATES)
    inputs = mti.resolve([f"{both[0]},{both[1]}"])
    assert (inputs.kind, inputs.layout, inputs.libraries) == (mti.BAM, mti.BOTH_MATES, [(both[0],), (both[1],)])
    assert inputs.describe() == "2 BAM files of both mates"
    _bams(tmp_path, ["A_R1.bam", "A_R2.bam", "B_R1.bam", "B_R2.bam"], mti.ONE_MATE)
    expected = [("A_R1.bam", "A_R2.bam"), ("B_R1.bam", "B_R2.bam")]
    for arguments in (
        ["*_R1.bam", "*_R2.bam"],
        ["A_R1.bam,B_R1.bam", "A_R2.bam,B_R2.bam"],
        ["A_R1.bam", "B_R1.bam", "A_R2.bam", "B_R2.bam"],
    ):
        inputs = mti.resolve(arguments)
        assert (inputs.kind, inputs.layout, inputs.libraries) == (mti.BAM, mti.ONE_MATE, expected), arguments
    assert inputs.describe() == "4 BAM files of one mate each"
    assert mti.resolve(["A_R1.bam", "A_R2.bam"]).describe() == "2 BAM files of one mate each"
    assert mti.resolve(["A.bam"]).describe() == "a BAM file of both mates"


@pytest.mark.parametrize("case", ["map and pairs", "pairs and BAM", "layouts", "one argument", "odd", "halves", "twice"])
def test_resolve_refusals(case, tmp_path, monkeypatch):
    """What the inputs of a run cannot be, with the reason."""
    monkeypatch.chdir(tmp_path)
    Path("L1.pairs").write_text("## pairs format v1.0\n")
    _bams(tmp_path, ["A.bam"], mti.BOTH_MATES)
    _bams(tmp_path, ["A_R1.bam", "A_R2.bam", "B_R1.bam", "B_R2.bam"], mti.ONE_MATE)
    metacool = contig_map(TESTS_DATA / "alignment.pairs", {"NODE_522": {"length": 22786}}, tmp_path / "x.metacool", columns=False)
    arguments, message = {
        "map and pairs": (
            [str(metacool), "L1.pairs"],
            f"A .metacool map is the one input of a run, given alone, and this run has 2: {metacool}, L1.pairs.",
        ),
        "pairs and BAM": (
            ["A.bam", "L1.pairs"],
            "The inputs of a run are pairs files or BAM files, not both: BAM 'A.bam'; pairs 'L1.pairs'.",
        ),
        "layouts": (
            ["A.bam", "A_R1.bam", "A_R2.bam"],
            "The BAM files of a run hold both mates of each read, or one mate each (an R1 and an R2 file per library), "
            "not both: both mates in 'A.bam'; one mate in 'A_R1.bam', 'A_R2.bam'.",
        ),
        "one argument": (
            ["A_R1.bam,A_R2.bam"],
            "One mate of each read in A_R1.bam, A_R2.bam: give the R1 files, then the R2 files, as two arguments (two "
            "comma-separated lists, or two patterns), or as many R1 arguments as R2 arguments.",
        ),
        "odd": (
            ["A_R1.bam", "B_R1.bam", "A_R2.bam"],
            "One mate of each read in A_R1.bam, B_R1.bam, A_R2.bam: give the R1 files, then the R2 files, as two "
            "arguments (two comma-separated lists, or two patterns), or as many R1 arguments as R2 arguments.",
        ),
        "halves": (
            ["A_R1.bam,B_R1.bam", "A_R2.bam"],
            "2 R1 file(s) and 1 R2 file(s): each BAM file of the first mates needs the BAM file of the second mates of its "
            "library. R1: A_R1.bam, B_R1.bam; R2: A_R2.bam.",
        ),
        "twice": (["L1.pairs", "./L1.pairs"], "The input ./L1.pairs is given twice: L1.pairs is the same file."),
    }[case]
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        mti.resolve(arguments)


def test_bam_layout(tmp_path):
    """A BAM file of both mates has reads flagged first and second of their
    pair; one of one mate has unpaired reads, or the first (or second) mates
    only, as hicstuff writes the alignments of bwa mem -5SPM split by mate."""
    lengths = {"c1": 1000, "c2": 1000}
    bam_header = alignments.header(lengths)
    cases = {
        "both": [(0x41, "c1"), (0x81, "c2")],
        "unpaired": [(0, "c1"), (0x10, "c2")],
        "first mates": [(0x41, "c1"), (0x41, "c2")],
        "secondary second mate": [(0x41, "c1"), (0x181, "c2")],
    }
    for name, reads in cases.items():
        path = tmp_path / f"{name}.bam"
        with pysam.AlignmentFile(str(path), "wb", header=bam_header) as bam:
            for i, (flag, contig) in enumerate(reads):
                bam.write(alignments.segment(bam_header, f"r{i}", contig, 10, flag=flag))
        expected = mti.BOTH_MATES if name == "both" else mti.ONE_MATE
        assert mti.bam_layout(path) == expected, name
    empty = tmp_path / "empty.bam"
    pysam.AlignmentFile(str(empty), "wb", header=bam_header).close()
    with pytest.raises(mtmc.MetacoolError, match=f"^{empty} holds no read.$"):
        mti.bam_layout(empty)


# ---------------------------------------------------------------------------
# The pairs of BAM files
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("order", ["aligner", "name", "coordinate", "coordinate, unsaid"])
def test_both_mates_in_any_order(order, community, tmp_path):
    """A BAM file of both mates gives the pairs it was written from, whatever
    its order: as bwa writes it (the mates of a read in a row), sorted by name,
    or sorted by coordinate, which its header says (grouped by name first,
    with samtools collate), or does not say (found when a read comes back)."""
    work = tmp_path / "work"
    work.mkdir()
    path = work / "lib.bam"
    alignments.both_mates(path, community["pairs"], community["lengths"], order=order.split(",")[0])
    if order == "coordinate, unsaid":
        with pysam.AlignmentFile(str(path)) as bam:
            reads = list(bam)
            unsaid = alignments.header(community["lengths"])
        with pysam.AlignmentFile(str(path), "wb", header=unsaid) as bam:
            for read in reads:
                bam.write(pysam.AlignedSegment.from_dict(read.to_dict(), unsaid))
    out = work / "lib.pairs"
    assert mti.bam_pairs((str(path),), str(out), community["lengths"], mti.BOTH_MATES, tmp_dir=str(work)) == 60000
    assert alignments.contacts(alignments.read_pairs(out)) == alignments.contacts(community["pairs"])
    # Grouped or sorted in a directory of its own, removed.
    assert sorted(path.name for path in work.iterdir()) == ["lib.bam", "lib.pairs"]


@pytest.mark.parametrize("order", ["name", "lexicographic", "lexicographic, unsaid", "aligner", "coordinate"])
def test_one_mate_in_any_order(order, community, tmp_path):
    """The BAM files of the two mates give the pairs they were written from,
    whatever their order: sorted by name, naturally (samtools sort -n) or
    lexicographically (samtools sort -N, which the header says, or another
    tool that does not say it), both in the order of the reads (as the aligner
    wrote them), or in another order (sorted by name first)."""
    work = tmp_path / "work"
    work.mkdir()
    paths = [work / "R1.bam", work / "R2.bam"]
    alignments.one_mate(paths, community["pairs"], community["lengths"], order=order.split(",")[0])
    if order.endswith("unsaid"):
        for path in paths:
            with pysam.AlignmentFile(str(path)) as bam:
                reads = list(bam)
            unsaid = alignments.header(community["lengths"], sort_order="queryname")
            with pysam.AlignmentFile(str(path), "wb", header=unsaid) as bam:
                for read in reads:
                    bam.write(pysam.AlignedSegment.from_dict(read.to_dict(), unsaid))
    out = work / "lib.pairs"
    library = tuple(str(path) for path in paths)
    assert mti.bam_pairs(library, str(out), community["lengths"], mti.ONE_MATE, tmp_dir=str(work)) == 60000
    assert alignments.contacts(alignments.read_pairs(out)) == alignments.contacts(community["pairs"])
    assert sorted(path.name for path in work.iterdir()) == ["R1.bam", "R2.bam", "lib.pairs"]


def _pairs_table(rows):
    return pd.DataFrame(rows, columns=alignments.PAIRS_COLUMNS)


@pytest.mark.parametrize("order", ["name", "lexicographic", "aligner"])
def test_read_names_with_numbers(order, tmp_path):
    """Read names whose numbers sort as numbers (r9 before r10), or are padded
    with zeros (r010 and r10 alike, for samtools sort -n), and that end with a
    mate suffix (/1, /2), match their mates whatever the order of the files;
    a read of one file only has no pair."""
    lengths = {"c1": 1000, "c2": 1000}
    names = ["r9", "r10", "r010", "r1.5", "r1", "x_1", "x1", "rA", "ra"]
    rows = [[name, "c1", 10 + i, "c2", 100 + i, "+", "-"] for i, name in enumerate(names)]
    pairs = _pairs_table(rows)
    paths = [tmp_path / "R1.bam", tmp_path / "R2.bam"]
    alignments.one_mate(paths, pairs, lengths, order=order)
    # The mate suffixes that bowtie2 keeps in the names it writes, and a read
    # of the first file only.
    for path, suffix in zip(paths, ["/1", "/2"], strict=True):
        with pysam.AlignmentFile(str(path)) as bam:
            reads, bam_header = list(bam), bam.header
        if path == paths[0]:
            reads.append(alignments.segment(bam_header, "zzz", "c1", 500))
        with pysam.AlignmentFile(str(path), "wb", header=bam_header) as bam:
            for read in reads:
                read.query_name += suffix
                bam.write(read)
    out = tmp_path / "lib.pairs"
    library = tuple(str(path) for path in paths)
    assert mti.bam_pairs(library, str(out), lengths, mti.ONE_MATE, tmp_dir=str(tmp_path)) == len(names)
    assert alignments.contacts(alignments.read_pairs(out)) == alignments.contacts(pairs)


def test_what_the_pairs_of_both_mates_leave_out(tmp_path, logged):
    """The pairs of the primary alignments of the two mates of each read, both
    aligned with a mapping quality of at least --min-quality on contigs of the
    assembly: not the secondary and supplementary alignments, a read with an
    unaligned mate, or one mate in the file, nor the pairs on a sequence that
    is not a contig of the assembly."""
    lengths = {"c1": 1000, "c2": 1000, "spike": 500}
    bam_header = alignments.header(lengths)
    rows = [
        ["kept", "c1", 10, "c2", 20, "+", "-"],
        ["low", "c1", 30, "c2", 40, "+", "+"],
        ["spiked", "c1", 50, "spike", 60, "-", "+"],
        ["other", "c2", 70, "c1", 80, "-", "-"],
    ]
    path = tmp_path / "lib.bam"
    extra = [
        # A secondary and a supplementary alignment, elsewhere: left out.
        alignments.segment(bam_header, "kept", "c2", 900, flag=0x141),
        alignments.segment(bam_header, "kept", "c1", 900, flag=0x881),
    ]
    alignments.both_mates(path, _pairs_table(rows), lengths, extra=extra)
    with pysam.AlignmentFile(str(path)) as bam:
        reads = list(bam)
    with pysam.AlignmentFile(str(path), "wb", header=bam_header) as bam:
        for read in reads:
            if read.query_name == "low" and read.is_read2:
                read.mapping_quality = 29
            if read.query_name == "other" and read.is_read2:
                continue
            bam.write(read)
        # A read whose second mate is not aligned.
        bam.write(alignments.segment(bam_header, "unaligned", "c1", 300, flag=0x41))
        bam.write(alignments.segment(bam_header, "unaligned", flag=0x81))
    out = tmp_path / "lib.pairs"
    assembly = {"c1": 1000, "c2": 1000}
    assert mti.bam_pairs((str(path),), str(out), assembly, mti.BOTH_MATES, min_quality=30, tmp_dir=str(tmp_path)) == 1
    assert alignments.contacts(alignments.read_pairs(out)) == alignments.contacts(_pairs_table(rows[:1]))
    assert mti.bam_pairs((str(path),), str(out), assembly, mti.BOTH_MATES, min_quality=29, tmp_dir=str(tmp_path)) == 2
    messages = [message for _, message in logged()]
    assert f"1 of the 5 reads of {path} have one mate in it: no pair." in messages
    assert f"1 pairs of {path} have a read on a sequence that is not a contig of the assembly: left out." in messages
    assert (
        f"{path}: 2 pairs of 4 reads with both mates (50.00%), both aligned with a mapping quality of at least 29 on "
        "contigs of the assembly." in messages
    )


def test_the_mates_of_other_libraries_are_refused(tmp_path):
    """BAM files of one mate whose reads match none of the other file's (the
    R1 file of a library with the R2 file of another) are refused."""
    lengths = {"c1": 1000, "c2": 1000}
    first = _pairs_table([[f"a{i}", "c1", 10, "c2", 20, "+", "-"] for i in range(3)])
    second = _pairs_table([[f"b{i}", "c1", 10, "c2", 20, "+", "-"] for i in range(3)])
    alignments.one_mate([tmp_path / "A_R1.bam", tmp_path / "A_R2.bam"], first, lengths)
    alignments.one_mate([tmp_path / "B_R1.bam", tmp_path / "B_R2.bam"], second, lengths)
    library = (str(tmp_path / "A_R1.bam"), str(tmp_path / "B_R2.bam"))
    message = f"{library[0]} and {library[1]} hold no read with both mates: they are not the R1 and R2 files of one library."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        mti.bam_pairs(library, str(tmp_path / "lib.pairs"), lengths, mti.ONE_MATE, tmp_dir=str(tmp_path))


def test_natural_key():
    """The order of samtools sort -n: numbers by value, their leading zeros
    ignored, and digits before letters as characters."""
    names = ["r10", "r9", "r1", "r1a", "r", "r_1", "R2", "r1.5"]
    assert sorted(names, key=mti.natural_key) == ["R2", "r", "r1", "r1.5", "r1a", "r9", "r10", "r_1"]
    assert mti.natural_key("r010") == mti.natural_key("r10")


# ---------------------------------------------------------------------------
# The headers, and the map
# ---------------------------------------------------------------------------


def test_check_references(tmp_path, logged):
    """The sequences of a BAM or pairs file must be the contigs of the assembly,
    with their lengths: another length, or none of them, is refused; others
    (a spike-in) are left out, with a warning."""
    lengths = {"c1": 1000, "c2": 2000}
    mti.check_references("lib.bam", [("c1", 1000), ("c2", 2000)], lengths)
    mti.check_references("lib.bam", [], lengths)
    mti.check_references("lib.bam", [("c1", 1000), ("spike", 10)], lengths)
    assert logged() == [
        (
            "WARNING",
            "1 of the reference sequences of lib.bam are not contigs of the assembly: the pairs with a read on them are "
            "left out: 'spike'.",
        )
    ]
    other = "lib.bam is not of the contigs of the assembly: 1 of its reference sequences have another length in the "
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(other)}assembly: 'c2' \\(1999 bp; 2000 bp in"):
        mti.check_references("lib.bam", [("c1", 1000), ("c2", 1999)], lengths)
    none = "lib.bam is not of the contigs of the assembly: none of its 2 reference sequences is a contig of the assembly: "
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(none)}'x1', 'x2'.$"):
        mti.check_references("lib.bam", [("x1", 1000), ("x2", 2000)], lengths)


def test_pairs_chromsizes(tmp_path):
    """The #chromsize lines of the header of a pairs file, plain or compressed
    with gzip or bzip2; none for a file compressed otherwise, which metacooler
    reads."""
    text = "## pairs format v1.0\n#chromsize: c1 1000\n#chromsize: c2 2000\n#columns: readID\n.\tc1\t5\tc2\t9\t+\t-\n"
    plain, zipped, bzipped, zstd = (tmp_path / name for name in ["a.pairs", "b.pairs.gz", "c.pairs.bz2", "d.pairs.zst"])
    plain.write_text(text)
    with gzip.open(zipped, "wt") as handle:
        handle.write(text)
    with bz2.open(bzipped, "wt") as handle:
        handle.write(text)
    zstd.write_bytes(b"\x28\xb5\x2f\xfd" + bytes(20))
    for path in (plain, zipped, bzipped):
        assert mti.pairs_chromsizes(path) == [("c1", 1000), ("c2", 2000)]
    assert mti.pairs_chromsizes(zstd) == []


@pytest.mark.parametrize(
    "assembly, name",
    [
        ("data/assembly.fa", "assembly"),
        ("assembly.fa.gz", "assembly"),
        ("sample.contigs.fasta", "sample.contigs"),
        ("MAGS.FNA", "MAGS"),
        ("assembly", "assembly"),
        ("assembly.gz", "assembly"),
        (".fa", ".fa"),
    ],
)
def test_map_name(assembly, name):
    """The map that a run makes is named after --prefix, or after the FASTA file
    of the assembly, without its extensions."""
    assert mti.assembly_name(assembly) == name
    assert mti.map_name(None, assembly) == f"{name}.metacool"
    assert mti.map_name("sample_1", assembly) == "sample_1.metacool"


def test_make_map_from_pairs_files(community, metagenome_files, tmp_path):
    """The map that a run makes from pairs files is the map that metacooler
    create pairs makes from them, one bin per contig (see tests.maps):
    several files, the libraries of a sample, sum their contacts, and their
    contact statistics."""
    libraries = _split_pairs(metagenome_files["pairs"], [tmp_path / "L1.pairs.gz", tmp_path / "L2.pairs.gz"])
    outfile = tmp_path / "out" / "assembly.metacool"
    outfile.parent.mkdir()
    inputs = mti.resolve([str(path) for path in libraries])
    mti.make_map(inputs, community["contig_data"], str(outfile), assembly="assembly.fa", command="metator network x")
    expected = contig_map(metagenome_files["pairs"], community["contig_data"], tmp_path / "expected.metacool", columns=False)
    _same_map(outfile, expected)
    mc = metacooler.MetaCooler(str(outfile))
    assert mc.info["assembly"] == "assembly" and mc.metadata["commands"] == ["metator network x"]
    assert set(mc.genomes["name"]) == {"__unbinned__"} and mc.info["sum"] == 60000
    # No file of the run left next to the map.
    assert sorted(path.name for path in outfile.parent.iterdir()) == ["assembly.metacool"]


@pytest.mark.parametrize("layout", [mti.BOTH_MATES, mti.ONE_MATE])
def test_make_map_from_bam_files(layout, community, metagenome_files, tmp_path):
    """The map that a run makes from BAM files is the map of the pairs that
    they were written from: their pairs, written in the temporary directory of
    the run, each mate at the 5' end of its alignment."""
    pairs = community["pairs"]
    half = len(pairs) // 2
    if layout == mti.BOTH_MATES:
        alignments.both_mates(tmp_path / "L1.bam", pairs[:half], community["lengths"])
        alignments.both_mates(tmp_path / "L2.bam", pairs[half:], community["lengths"], order="coordinate")
        arguments = [str(tmp_path / "L1.bam"), str(tmp_path / "L2.bam")]
    else:
        alignments.one_mate([tmp_path / "L1_R1.bam", tmp_path / "L1_R2.bam"], pairs[:half], community["lengths"])
        alignments.one_mate([tmp_path / "L2_R1.bam", tmp_path / "L2_R2.bam"], pairs[half:], community["lengths"])
        arguments = [str(tmp_path / "L*_R1.bam"), str(tmp_path / "L*_R2.bam")]
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    outfile = tmp_path / "assembly.metacool"
    mti.make_map(mti.resolve(arguments), community["contig_data"], str(outfile), tmp_dir=str(tmp_dir))
    expected = contig_map(metagenome_files["pairs"], community["contig_data"], tmp_path / "expected.metacool", columns=False)
    _same_map(outfile, expected)
    # The pairs of each library, in the temporary directory of the run.
    assert sorted(path.name for path in tmp_dir.iterdir()) == ["library_1.pairs", "library_2.pairs"]


@pytest.mark.parametrize("case", ["pairs of another assembly", "no pair on the contigs", "BAM of another length"])
def test_make_map_refusals_leave_the_map(case, community, tmp_path):
    """Inputs that are not of the assembly are refused, and leave the map that
    an earlier run left where the run writes its map as it was, with no file
    of the run next to it."""
    outfile = tmp_path / "out" / "assembly.metacool"
    outfile.parent.mkdir()
    outfile.write_text("an earlier map")
    lengths = community["lengths"]
    if case == "pairs of another assembly":
        path = shutil.copyfile(TESTS_DATA / "subset.pairs.gz", tmp_path / "other.pairs.gz")
        message = f"{path} is not of the contigs of the assembly: none of its 12507 sequences (#chromsize) is a contig"
    elif case == "no pair on the contigs":
        path = tmp_path / "bare.pairs"
        path.write_text("## pairs format v1.0\n.\tx1\t5\tx2\t9\t+\t-\n")
        message = "No pair from a pairs file is on the contigs of the assembly"
    else:
        path = tmp_path / "other.bam"
        other = {name: length + (name == "contig_001") for name, length in lengths.items()}
        alignments.both_mates(path, community["pairs"][:10], other)
        message = f"{path} is not of the contigs of the assembly: 1 of its reference sequences have another length"
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}"):
        mti.make_map(mti.resolve([str(path)]), community["contig_data"], str(outfile), tmp_dir=str(tmp_path))
    assert outfile.read_text() == "an earlier map"
    assert sorted(os.listdir(outfile.parent)) == ["assembly.metacool"]


# ---------------------------------------------------------------------------
# The command line: network and pipeline from pairs or BAM files
# ---------------------------------------------------------------------------


def _open_by_suffix(path, mode):
    """A file opened to be written, compressed as its suffix says (.gz, .bz2)."""
    return {".gz": gzip.open, ".bz2": bz2.open}.get(Path(path).suffix, open)(path, mode)


#: The inputs that a run takes, besides a map (see _inputs).
INPUT_KINDS = ["pairs as a pattern", "pairs as a list", "pairs as two arguments", "BAM of both mates", "BAM of one mate"]


def _inputs(kind, community, metagenome_files, directory):
    """The arguments of a run from the pairs of the mock community, split into
    two libraries: pairs files, compressed with gzip and bzip2 (as a pattern, a
    list, or two arguments), or BAM files of both mates (one sorted by
    coordinate), or of one mate each (the R1 files, then the R2 files; one
    library in the order of its reads, as the aligner writes it)."""
    directory.mkdir()
    pairs, lengths = community["pairs"], community["lengths"]
    half = len(pairs) // 2
    if kind.startswith("pairs"):
        paths = [directory / "L1.pairs.gz", directory / "L2.pairs.bz2"]
        libraries = [str(path) for path in _split_pairs(metagenome_files["pairs"], paths, compress=_open_by_suffix)]
        return {
            "pairs as a pattern": [str(directory / "L*.pairs.*")],
            "pairs as a list": [",".join(libraries)],
            "pairs as two arguments": libraries,
        }[kind]
    if kind == "BAM of both mates":
        alignments.both_mates(directory / "L1.bam", pairs[:half], lengths)
        alignments.both_mates(directory / "L2.bam", pairs[half:], lengths, order="coordinate")
        return [str(directory / "L1.bam"), str(directory / "L2.bam")]
    alignments.one_mate([directory / "L1_R1.bam", directory / "L1_R2.bam"], pairs[:half], lengths)
    alignments.one_mate([directory / "L2_R1.bam", directory / "L2_R2.bam"], pairs[half:], lengths, order="aligner")
    return [f"{directory / 'L1_R1.bam'},{directory / 'L2_R1.bam'}", f"{directory / 'L1_R2.bam'},{directory / 'L2_R2.bam'}"]


@pytest.mark.parametrize("kind", INPUT_KINDS)
def test_network_from_pairs_or_bam_files(kind, community, metagenome_files, tmp_path, run_cli):
    """metator network from the pairs files or the BAM files of a run makes the
    map of the contigs in its output directory, named after the assembly, and
    stores in it the network of metator network on the map of these pairs:
    the same map, network and contig data. Its metadata give the inputs as
    given, and the map that the run made, as the map it started from."""
    common = ["-a", str(metagenome_files["assembly"]), "-e", "DpnII", "-T", str(tmp_path / "tmp")]
    expected = contig_map(metagenome_files["pairs"], community["contig_data"], tmp_path / "assembly.metacool", columns=False)
    run_cli("network", str(expected), "-o", str(tmp_path / "from_map"), *common)
    arguments = _inputs(kind, community, metagenome_files, tmp_path / "inputs")
    options = ["-q", "30", "-t", "2"] if kind.startswith("BAM") else []
    out = tmp_path / "out"
    run_cli("network", *arguments, "-o", str(out), *options, *common)

    written, reference = out / "assembly.metacool", tmp_path / "from_map" / "assembly.metacool"
    _same_map(written, reference)
    pd.testing.assert_frame_equal(mtn.map_network(str(written)), mtn.map_network(str(reference)))
    pd.testing.assert_frame_equal(mtn.map_contig_table(str(written)), mtn.map_contig_table(str(reference)))
    entry = metacooler.MetaCooler(str(written)).metadata["metator"]["network"]
    assert entry["parameters"]["input"] == arguments and entry["parameters"]["map"] == str(written)
    assert entry["parameters"]["min_quality"] == 30 and entry["command"].startswith("metator network ")
    (log,) = out.glob("metator_network_*.log")
    logged = log.read_text()
    assert f"## input: {shlex.join(arguments)}\n" in logged
    what = "2 pairs files" if kind.startswith("pairs") else f"{4 if kind.endswith('one mate') else 2} BAM files"
    assert f"INFO :: The run makes its map, {written}, from {what}" in logged
    assert sorted(path.name for path in out.iterdir()) == ["assembly.metacool", log.name]
    # The temporary directory of the run, with the pairs of the BAM files, is
    # removed.
    assert not (tmp_path / "tmp").exists()


def test_a_run_from_bam_files_takes_its_pairs_from_the_bam_files(community, metagenome_files, tmp_path, run_cli):
    """The pairs of the BAM files of a run are written in its temporary
    directory, one file per library, which -N keeps."""
    arguments = _inputs("BAM of one mate", community, metagenome_files, tmp_path / "inputs")
    run_cli("network", *arguments, "-a", str(metagenome_files["assembly"]), "-o", str(tmp_path / "out"), "-T",
            str(tmp_path / "tmp"), "-N")  # fmt: skip
    (kept,) = (tmp_path / "tmp").iterdir()
    found = [alignments.read_pairs(kept / f"library_{i}.pairs") for i in (1, 2)]
    assert alignments.contacts(pd.concat(found)) == alignments.contacts(community["pairs"])


@pytest.mark.parametrize(
    "case",
    [
        "a map there",
        "--inplace",
        "--min-quality with pairs",
        "--min-quality with a map",
        "a normalization that needs the restriction sites",
        "an input where the map goes",
        "an output directory that holds ::",
        "pairs of another assembly",
    ],
)
@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_refusals_of_a_run_from_pairs_or_bam_files(command, case, community, metagenome_files, tmp_path, run_cli):
    """What a run from pairs or BAM files refuses, before it reads them (but
    the pairs of another assembly, refused once the headers of the inputs are
    read), the reason in its log: a map left where the run writes its map, as
    it was, and no other file of the run than its log."""
    arguments = _inputs("pairs as two arguments", community, metagenome_files, tmp_path / "inputs")
    out = tmp_path / "out"
    out.mkdir()
    options = []
    map_there = out / "assembly.metacool"
    message = None
    if case == "a map there":
        map_there.write_text("an earlier map")
        message = (
            f"{map_there} exists, where the run writes the map that it makes from 2 pairs files: give -F to replace it, "
            "or another output directory with -o; or start the run from that map."
        )
    elif case == "--inplace":
        options = ["--inplace"]
        message = (
            "--inplace updates the map that a run starts from, and this run starts from 2 pairs files, which it makes "
            f"its map from, in its output directory ({map_there}): leave out --inplace."
        )
    elif case == "--min-quality with pairs":
        options = ["-q", "20"]
        message = (
            "--min-quality filters the alignments of BAM files, and the inputs of the run are 2 pairs files: leave out "
            "--min-quality."
        )
    elif case == "--min-quality with a map":
        arguments = [str(contig_map(metagenome_files["pairs"], community["contig_data"], tmp_path / "x.metacool", columns=False))]
        options = ["--min-quality", "30"]
        message = (
            "--min-quality filters the alignments of BAM files, and the inputs of the run are a map: leave out --min-quality."
        )
    elif case == "a normalization that needs the restriction sites":
        options = ["-n" if command == "network" else "-m", "RS"]
        message = (
            'For "RS" and "theoritical_hit" normalization, enzyme is required, as the map that the run makes from its '
            "pairs or BAM files gives no restriction sites."
        )
    elif case == "an input where the map goes":
        shutil.copyfile(arguments[0], map_there)
        arguments = [str(map_there)]
        options = ["-F"]
        message = (
            f"{map_there}, where the run writes its map, is {map_there}, an input of the run: give the run another "
            "output directory with -o."
        )
    elif case == "an output directory that holds ::":
        out = tmp_path / "a::b"
        message = (
            f"{out / 'assembly.metacool'}, the map that the run makes in its output directory, holds ::, which metacooler "
            "reads as file::group: give the run another output directory with -o."
        )
    else:
        arguments = [str(shutil.copyfile(TESTS_DATA / "subset.pairs.gz", tmp_path / "other.pairs.gz"))]
        map_there.write_text("an earlier map")
        options = ["-F"]
        message = f"{arguments[0]} is not of the contigs of the assembly: none of its 12507 sequences (#chromsize) is a contig"
    with pytest.raises(ValueError):
        run_cli(command, *arguments, "-a", str(metagenome_files["assembly"]), "-o", str(out), "-T", str(tmp_path / "tmp"),
                *options)  # fmt: skip
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}" in log.read_text()
    kept = [log.name] + (["assembly.metacool"] if map_there.exists() and out == map_there.parent else [])
    assert sorted(path.name for path in out.iterdir()) == sorted(kept)
    if case in ("a map there", "pairs of another assembly"):
        assert map_there.read_text() == "an earlier map"


def test_pipeline_from_pairs_refuses_mge_contigs_of_another_assembly(community, metagenome_files, tmp_path, run_cli):
    """pipeline --mges from pairs or BAM files: a list none of whose contigs is a
    contig of the assembly, which the map that the run makes would hold none
    of, is refused once the run has read the assembly, before it makes its map
    from them: no map, and no other file of the run than its log."""
    arguments = _inputs("pairs as two arguments", community, metagenome_files, tmp_path / "inputs")
    out = tmp_path / "out"
    mges = tmp_path / "mges.txt"
    mges.write_text("not_a_contig\n")
    assembly = str(metagenome_files["assembly"])
    with pytest.raises(ValueError):
        run_cli("pipeline", *arguments, "-a", assembly, "-o", str(out), "-T", str(tmp_path / "tmp"), "--mges", str(mges))
    (log,) = out.glob("metator_pipeline_*.log")
    message = f"None of the MGE contigs is a contig of {assembly}: there is no MGE contig to cluster."
    assert f"ERROR :: {message}\n" in log.read_text()
    assert sorted(path.name for path in out.iterdir()) == [log.name]


def test_a_map_there_is_replaced_with_force(community, metagenome_files, tmp_path, run_cli):
    """With -F, the map of an earlier run where the run writes the map that it
    makes is replaced (a symbolic link there, not its target)."""
    arguments = _inputs("pairs as two arguments", community, metagenome_files, tmp_path / "inputs")
    out = tmp_path / "out"
    out.mkdir()
    target = tmp_path / "earlier.metacool"
    target.write_text("an earlier map")
    (out / "assembly.metacool").symlink_to(target)
    run_cli("network", *arguments, "-a", str(metagenome_files["assembly"]), "-o", str(out), "-F", "-T", str(tmp_path / "tmp"))
    written = out / "assembly.metacool"
    assert not written.is_symlink() and target.read_text() == "an earlier map"
    assert metacooler.MetaCooler(str(written)).info["sum"] == 60000


# ---------------------------------------------------------------------------
# The pairs of BAM files: details
# ---------------------------------------------------------------------------


def test_pairs_as_an_upper_triangle(tmp_path):
    """Each pair is written with the contig that comes first in the assembly
    first, or the mate at the lower position first on one contig, whichever
    mate is the first of the read."""
    lengths = {"c1": 1000, "c2": 1000}
    rows = [["a", "c2", 10, "c1", 20, "+", "-"], ["b", "c1", 90, "c1", 30, "-", "+"], ["c", "c1", 5, "c2", 5, "+", "+"]]
    alignments.both_mates(tmp_path / "lib.bam", _pairs_table(rows), lengths)
    out = tmp_path / "lib.pairs"
    mti.bam_pairs((str(tmp_path / "lib.bam"),), str(out), lengths, mti.BOTH_MATES, tmp_dir=str(tmp_path))
    written = alignments.read_pairs(out)[alignments.PAIRS_COLUMNS[1:]].to_numpy().tolist()
    assert written == [["c1", 20, "c2", 10, "-", "+"], ["c1", 30, "c1", 90, "+", "-"], ["c1", 5, "c2", 5, "+", "+"]]


def test_the_grouping_check_remembers_the_latest_reads():
    """The check of a BAM file of both mates remembers the latest reads with one
    mate only: one that comes back among them is split (the file is not
    grouped by name); an older one is forgotten."""
    check = mti._GroupingCheck(size=2)
    for name in ["r1", "r2", "r3"]:
        check.lone(name)
    check.lone("r1")
    with pytest.raises(mti._NotGrouped):
        check.lone("r3")


def test_one_mate_files_keep_the_first_primary_alignment(tmp_path):
    """In the BAM file of one mate, the secondary and supplementary alignments
    of a read are left out, and so is any primary alignment of it after the
    first, as hicstuff's iterative alignment writes some."""
    lengths = {"c1": 1000, "c2": 1000}
    pairs = _pairs_table([["r1", "c1", 10, "c2", 20, "+", "-"], ["r2", "c1", 30, "c2", 40, "+", "+"]])
    paths = [tmp_path / "R1.bam", tmp_path / "R2.bam"]
    alignments.one_mate(paths, pairs, lengths, order="aligner")
    with pysam.AlignmentFile(str(paths[0])) as bam:
        reads, bam_header = list(bam), bam.header
    with pysam.AlignmentFile(str(paths[0]), "wb", header=bam_header) as bam:
        for read in reads:
            bam.write(read)
            if read.query_name == "r1":
                bam.write(alignments.segment(bam_header, "r1", "c2", 500, flag=0x100))
                bam.write(alignments.segment(bam_header, "r1", "c2", 600, flag=0x800))
                bam.write(alignments.segment(bam_header, "r1", "c2", 700))
    out = tmp_path / "lib.pairs"
    library = tuple(str(path) for path in paths)
    assert mti.bam_pairs(library, str(out), lengths, mti.ONE_MATE, tmp_dir=str(tmp_path)) == 2
    assert alignments.contacts(alignments.read_pairs(out)) == alignments.contacts(pairs)


def test_names_tied_in_the_order_of_samtools(tmp_path):
    """Names that samtools sort -n puts at one place (r010 and r10, their
    numbers equal) come in no set order in each file: the reads of that place
    are matched by name, in whichever order the two files have them."""
    lengths = {"c1": 1000, "c2": 1000}
    pairs = _pairs_table([["r010", "c1", 10, "c2", 20, "+", "-"], ["r10", "c1", 30, "c2", 40, "+", "+"]])
    bam_header = alignments.header(lengths, sort_order="queryname")
    paths = [tmp_path / "R1.bam", tmp_path / "R2.bam"]
    mates = {name: (one, two) for (one, two), name in zip(alignments._mates(bam_header, pairs, False), pairs.readID, strict=True)}
    for path, mate, names in [(paths[0], 0, ["r010", "r10"]), (paths[1], 1, ["r10", "r010"])]:
        with pysam.AlignmentFile(str(path), "wb", header=bam_header) as bam:
            for name in names:
                bam.write(mates[name][mate])
    out = tmp_path / "lib.pairs"
    library = tuple(str(path) for path in paths)
    assert mti.bam_pairs(library, str(out), lengths, mti.ONE_MATE, tmp_dir=str(tmp_path)) == 2
    assert alignments.contacts(alignments.read_pairs(out)) == alignments.contacts(pairs)


@pytest.mark.parametrize("case", ["no reference", "corrupt"])
def test_bam_files_that_cannot_give_pairs(case, community, tmp_path):
    """A BAM file that lists no reference sequence (its reads are not
    aligned), or that pysam cannot read, is refused, with the reason."""
    path = tmp_path / "lib.bam"
    if case == "no reference":
        bam_header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
        with pysam.AlignmentFile(str(path), "wb", header=bam_header) as bam:
            for flag in (0x45, 0x85):
                read = pysam.AlignedSegment(bam_header)
                read.query_name, read.flag, read.query_sequence = "r1", flag, "ACGT"
                bam.write(read)
        message = f"{path} lists no reference sequence in its header: its reads are not aligned."
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            mti.make_map(mti.resolve([str(path)]), community["contig_data"], str(tmp_path / "x.metacool"), tmp_dir=str(tmp_path))
    else:
        with gzip.open(path, "wb") as handle:
            handle.write(b"BAM\x01" + bytes(range(64)))
        with pytest.raises(mtmc.MetacoolError, match=f"^Cannot read the BAM file {re.escape(str(path))}: "):
            mti.resolve([str(path)])


def test_mates_at_their_five_prime_end(tmp_path):
    """Each mate is at the 5' end of its alignment, soft-clipped bases left out:
    its first aligned base on the forward strand, its last on the reverse
    strand, as pairtools and hicstuff's mapping mode 5spm record it (metaTOR
    1.x recorded the leftmost base)."""
    lengths = {"c1": 1000, "c2": 1000}
    bam_header = alignments.header(lengths)
    path = tmp_path / "lib.bam"
    with pysam.AlignmentFile(str(path), "wb", header=bam_header) as bam:
        for name, (contig1, start1, cigar1, reverse1), (contig2, start2, cigar2, reverse2) in [
            ("a", ("c1", 100, "50M", False), ("c2", 200, "50M", True)),
            ("b", ("c1", 300, "10S40M", True), ("c1", 600, "5S45M", False)),
        ]:
            for flag, contig, start, cigar, reverse in [
                (0x41, contig1, start1, cigar1, reverse1),
                (0x81, contig2, start2, cigar2, reverse2),
            ]:
                read = pysam.AlignedSegment(bam_header)
                read.query_name, read.flag = name, flag | (0x10 if reverse else 0)
                read.reference_id, read.reference_start = bam_header.get_tid(contig), start
                read.mapping_quality, read.cigarstring, read.query_sequence = 60, cigar, "A" * 50
                bam.write(read)
    out = tmp_path / "lib.pairs"
    mti.bam_pairs((str(path),), str(out), lengths, mti.BOTH_MATES, tmp_dir=str(tmp_path))
    written = alignments.read_pairs(out)[alignments.PAIRS_COLUMNS[1:]].to_numpy().tolist()
    # a: 101 (0-based 100, forward), 250 (200 + 50, reverse); b: 340 (300 + 40
    # aligned bases, reverse), 601 (forward, after its soft clip).
    assert written == [["c1", 101, "c2", 250, "+", "-"], ["c1", 340, "c1", 601, "-", "+"]]


@pytest.mark.parametrize("layout", [mti.BOTH_MATES, mti.ONE_MATE])
def test_duplicates_and_reads_failing_qc_are_left_out(layout, tmp_path, logged):
    """A read with a mate flagged as a PCR or optical duplicate (0x400, as
    Picard MarkDuplicates or samtools markdup flag them) or as failing QC
    (0x200) gives no pair, and the log counts them."""
    lengths = {"c1": 1000, "c2": 1000}
    rows = [[name, "c1", 10 * i + 10, "c2", 10 * i + 500, "+", "-"] for i, name in enumerate(["kept", "dup", "qc"])]
    pairs = _pairs_table(rows)
    flags = {("dup", 1): 0x400, ("dup", 2): 0x400, ("qc", 2): 0x200}
    if layout == mti.BOTH_MATES:
        paths = [tmp_path / "lib.bam"]
        alignments.both_mates(paths[0], pairs, lengths)
    else:
        paths = [tmp_path / "R1.bam", tmp_path / "R2.bam"]
        alignments.one_mate(paths, pairs, lengths)
    for i, path in enumerate(paths, 1):
        with pysam.AlignmentFile(str(path)) as bam:
            reads, bam_header = list(bam), bam.header
        with pysam.AlignmentFile(str(path), "wb", header=bam_header) as bam:
            for read in reads:
                mate = 2 if read.is_read2 or i == 2 else 1
                read.flag |= flags.get((read.query_name, mate), 0)
                bam.write(read)
    out = tmp_path / "lib.pairs"
    library = tuple(str(path) for path in paths)
    assert mti.bam_pairs(library, str(out), lengths, layout, tmp_dir=str(tmp_path)) == 1
    assert alignments.contacts(alignments.read_pairs(out)) == alignments.contacts(pairs[:1])
    what = library[0] if layout == mti.BOTH_MATES else f"{library[0]} and {library[1]}"
    assert (
        "INFO",
        f"1 reads of {what} with a mate flagged as a PCR or optical duplicate, and 1 with a mate flagged as failing "
        "QC: left out.",
    ) in logged()
