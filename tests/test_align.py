# Test for align module

import functools
import glob
import gzip
import itertools
import os
import random
import string

import metator.align as mta
import metator.io as mio
import pysam
import pytest

CONTIG_DATA = {"c1": {"id": 1, "length": 10_000}, "c2": {"id": 2, "length": 20_000}}


def _write_bam(path, reads):
    """Write a name-sorted BAM of 50 bp reads.

    Each read is (name, contig index, 0-based start, flag[, mapping quality]).
    """
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "c1", "LN": 10_000}, {"SN": "c2", "LN": 20_000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, ref_id, start, flag, *mapq in reads:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.flag = flag
            read.reference_id = ref_id
            read.reference_start = start
            read.mapping_quality = mapq[0] if mapq else 60
            read.cigartuples = [(0, 50)]
            bam.write(read)


def _first_reads(fastq, n_reads):
    """The lines of the first reads of a gzipped fastq."""
    with gzip.open(fastq, "rt") as reads:
        return list(itertools.islice(reads, 4 * n_reads))


def test_max_read_length(tmp_path):
    """The longest read, not the first one; gzipped or not."""
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGTACG\n+\nIIIIIII\n@r3\n\n+\n\n")
    assert mta.max_read_length(str(reads)) == 7
    assert mta.max_read_length("tests_data/for_paired.fq.gz") == 151
    (tmp_path / "empty.fq").write_text("")
    assert mta.max_read_length(str(tmp_path / "empty.fq")) == 0


@pytest.mark.parametrize("first_read_length", [None, 35], ids=["untrimmed", "trimmed-first-read"])
def test_align_iterative_extends_reads(tmp_path, first_read_length):
    """Reads that do not align once truncated are extended and aligned again,
    up to the longest read (151 bp in the test data), even when the first read
    of the file is shorter. Told that the reads were 20 bp long, hicstuff never
    aligned more than their first 40 bp; left to itself, it extends them up to
    the length of the first read only."""
    lines = _first_reads("tests_data/for_paired.fq.gz", 2000)
    if first_read_length:
        # Trim the sequence and quality of the first read, as an adapter
        # trimmer would.
        for i in (1, 3):
            lines[i] = lines[i][:first_read_length] + "\n"
    reads = tmp_path / "reads.fq"
    reads.write_text("".join(lines))
    bam = tmp_path / "alignment.bam"
    mta.align(str(reads), "tests_data/assembly", "bowtie2", str(bam), 1, str(tmp_path), iterative=True)
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as alignment:
        lengths = {read.infer_read_length() for read in alignment.fetch(until_eof=True) if not read.is_unmapped}
    assert max(lengths) == 151


@pytest.mark.parametrize("min_qual", [30, 40])
def test_align_iterative_min_quality(tmp_path, min_qual):
    """A read is kept at the first length at which it aligns with a mapping
    quality of at least min_qual, and extended otherwise. hicstuff was always
    given 40, whatever --min-quality, so a read aligned at a mapping quality of
    30 to 39 was extended even with the default --min-quality of 30."""
    lines = _first_reads("tests_data/for_paired.fq.gz", 2000)
    lengths = {lines[i][1:].split()[0]: len(lines[i + 1].rstrip()) for i in range(0, len(lines), 4)}
    reads = tmp_path / "reads.fq"
    reads.write_text("".join(lines))
    bam = tmp_path / "alignment.bam"
    mta.align(str(reads), "tests_data/assembly", "bowtie2", str(bam), 1, str(tmp_path), iterative=True, min_qual=min_qual)
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as alignment:
        # Mapping quality of the reads kept before their full length.
        truncated = [
            read.mapping_quality
            for read in alignment.fetch(until_eof=True)
            if read.flag in (0, 16) and read.infer_read_length() < lengths[read.query_name]
        ]
    assert min(truncated) >= min_qual
    assert any(mapq < 40 for mapq in truncated) == (min_qual < 40)


@pytest.mark.parametrize("aligner", ["bowtie2", "bwa"])
def test_align_iterative_index_in_parent_directory(tmp_path, monkeypatch, aligner):
    """An index given as ../idx/<name> is found in iterative mode, as in
    normal mode. hicstuff looks the index up again and strips the leading dots
    of the path it finds: it aligned the reads on /idx/<name>, which does not
    exist, and the run stopped with FileNotFoundError."""
    index_dir = tmp_path / "idx"
    index_dir.mkdir()
    if aligner == "bwa":
        index = mio.generate_fasta_index(os.path.abspath("tests_data/assembly.fa"), "bwa", str(index_dir))
    else:
        for path in glob.glob("tests_data/assembly*.bt2"):
            (index_dir / os.path.basename(path)).symlink_to(os.path.abspath(path))
        index = str(index_dir / "assembly")
    reads = tmp_path / "reads.fq"
    reads.write_text("".join(_first_reads("tests_data/for_paired.fq.gz", 200)))
    bam = tmp_path / "alignment.bam"
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    monkeypatch.chdir(run_dir)
    relative_index = os.path.join("..", "idx", os.path.basename(index))
    mta.align(str(reads), relative_index, aligner, str(bam), 1, str(tmp_path), iterative=True)
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as alignment:
        assert any(read.flag in (0, 16) for read in alignment.fetch(until_eof=True))


@pytest.mark.parametrize("aligner", ["bowtie2", "bwa"])
def test_get_contact_pairs_iterative(tmp_path, monkeypatch, aligner):
    """In iterative mode, both aligners align each mate on its own, then the
    mates are paired, one pair per read pair. With bwa, the run crashed once the
    reads were aligned (UnboundLocalError): it looked for the BAM that bwa
    writes for both mates in the other modes. hicstuff writes some bwa reads
    several times (read 2,712 of the test data, on both mates): they make one
    pair, not one per copy. The minimum mapping quality is passed on to
    hicstuff."""
    min_quals = []
    iterative_align = mta.hci.iterative_align

    def record_min_qual(*args, **kwargs):
        min_quals.append(kwargs["min_qual"])
        return iterative_align(*args, **kwargs)

    monkeypatch.setattr(mta.hci, "iterative_align", record_min_qual)
    fastqs = []
    for mate in ("for", "rev"):
        fastq = tmp_path / f"{mate}.fq"
        fastq.write_text("".join(_first_reads(f"tests_data/{mate}_paired.fq.gz", 3000)))
        fastqs.append(str(fastq))
    if aligner == "bwa":
        index = mio.generate_fasta_index("tests_data/assembly.fa", "bwa", str(tmp_path))
    else:
        index = "tests_data/assembly"
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (pairs_file,), _, _ = mta.get_contact_pairs(
        *fastqs,
        index,
        "tests_data/assembly.fa",
        aligner,
        "iterative",
        35,
        "fastq",
        None,
        None,
        str(out_dir),
        str(tmp_path),
        1,
    )
    with gzip.open(pairs_file, "rt") as pairs:
        names = [line.split("\t", 1)[0] for line in pairs if not line.startswith("#")]
    assert len(names) > 1000
    assert len(set(names)) == len(names)
    assert min_quals == [35, 35]


def test_process_bwa_bamfile(tmp_path):
    """Each mate keeps its own 1-based position; pairs are upper-triangle."""
    bam = tmp_path / "alignment.bam"
    _write_bam(
        bam,
        [
            ("r1", 0, 100, 65),  # read 1 on c1, forward
            ("r1", 1, 5000, 145),  # read 2 on c2, reverse
            ("r2", 1, 300, 65),  # both on c2, read 1 downstream:
            ("r2", 1, 200, 129),  # written swapped
        ],
    )
    out = tmp_path / "alignment.pairs"
    n_pairs = mta.process_bwa_bamfile(str(bam), 30, CONTIG_DATA, str(out))
    pairs = [line.split() for line in out.read_text().splitlines() if not line.startswith("#")]
    assert n_pairs == 2
    assert pairs == [
        ["r1", "c1", "101", "c2", "5001", "+", "-"],
        ["r2", "c2", "201", "c2", "301", "+", "+"],
    ]


def test_process_bamfile(tmp_path):
    """Keeps primary alignments (flag 0 or 16) above the quality threshold."""
    bam = tmp_path / "for.bam"
    _write_bam(
        bam,
        [
            ("r1", 0, 100, 0),  # kept, forward
            ("r2", 1, 5000, 16),  # kept, reverse
            ("r3", 0, 200, 0, 10),  # mapping quality too low
            ("r4", 0, 300, 256),  # secondary alignment
        ],
    )
    out = tmp_path / "for.tsv"
    assert mta.process_bamfile(str(bam), 30, str(out)) == 2
    assert [line.split("\t") for line in out.read_text().splitlines()] == [
        ["r1", "c1", "100", "+"],
        ["r2", "c2", "5000", "-"],
    ]


def test_process_bamfile_keeps_one_alignment_per_read(tmp_path):
    """The first alignment of a read that passes the filters is kept, its
    other copies are not: hicstuff's iterative alignment with bwa writes a read
    once per round when it also has a supplementary alignment."""
    bam = tmp_path / "for.bam"
    _write_bam(
        bam,
        [
            ("r1", 0, 100, 0),  # kept
            ("r1", 0, 100, 0),  # the same alignment, from a later round
            ("r2", 1, 5000, 16, 10),  # mapping quality too low
            ("r2", 1, 4980, 16),  # kept
            ("r2", 1, 4960, 16),  # from a later round
        ],
    )
    out = tmp_path / "for.tsv"
    assert mta.process_bamfile(str(bam), 30, str(out)) == 2
    assert [line.split("\t") for line in out.read_text().splitlines()] == [
        ["r1", "c1", "100", "+"],
        ["r2", "c2", "4980", "-"],
    ]


def test_merge_alignement(tmp_path):
    """Only reads aligned on both sides become pairs, in upper-triangle order."""
    forward = tmp_path / "for.tsv"
    reverse = tmp_path / "rev.tsv"
    forward.write_text("r1\tc1\t100\t+\nr2\tc1\t50\t+\nr4\tc2\t300\t+\n")
    reverse.write_text("r1\tc2\t5000\t-\nr3\tc2\t60\t-\nr4\tc2\t200\t+\n")
    out = tmp_path / "alignment.pairs"
    assert mta.merge_alignment(str(forward), str(reverse), CONTIG_DATA, str(out)) == 2
    pairs = [line.split() for line in out.read_text().splitlines() if not line.startswith("#")]
    assert pairs == [
        ["r1", "c1", "101", "c2", "5001", "+", "-"],
        ["r4", "c2", "201", "c2", "301", "+", "+"],
    ]


@pytest.mark.parametrize(
    ("forward_names", "reverse_names", "paired"),
    [
        # samtools sort -n puts r.5 before r5, and r5 before r_1.
        (["r5"], ["r.5", "r5"], ["r5"]),
        (["r5", "r_1"], ["r_1"], ["r_1"]),
        (["r1", "r5", "r_1"], ["r.5", "r5", "r_1"], ["r5", "r_1"]),
        (["r07", "r8"], ["r7", "r8"], ["r8"]),
    ],
)
def test_merge_alignement_follows_samtools_name_order(tmp_path, forward_names, reverse_names, paired):
    """Reads are paired in the order of samtools sort -n."""
    forward = tmp_path / "for.tsv"
    reverse = tmp_path / "rev.tsv"
    forward.write_text("".join(f"{name}\tc1\t10\t+\n" for name in forward_names))
    reverse.write_text("".join(f"{name}\tc2\t20\t-\n" for name in reverse_names))
    out = tmp_path / "alignment.pairs"
    assert mta.merge_alignment(str(forward), str(reverse), CONTIG_DATA, str(out)) == len(paired)
    pairs = [line.split() for line in out.read_text().splitlines() if not line.startswith("#")]
    assert pairs == [[name, "c1", "11", "c2", "21", "+", "-"] for name in paired]


@pytest.mark.parametrize(
    ("first", "second", "expected"),
    [
        ("r9", "r10", -1),  # digit runs compare as numbers
        ("r10", "r9", 1),
        ("r12", "r13", ord("2") - ord("3")),  # same length: first digit that differs
        ("r07", "r7", 0),  # leading zeros are ignored
        ("r007", "r8", -1),
        ("r0", "r00", 0),
        ("r01", "r001", 0),
        ("r0b", "r00", ord("b")),
        ("a12b", "a7b", 1),
        ("r1a", "r1", ord("a")),
        ("r", "r1", -1),  # a prefix comes first
        ("", "r", -1),
        ("", "", 0),
        ("r_1", "r1", ord("_") - ord("1")),  # any other character by its code
        ("r.5", "r5", ord(".") - ord("5")),
        ("rAb", "rAB", ord("b") - ord("B")),
        ("r1.10", "r1.9", 1),
        ("x12y03", "x012y3", 0),
        ("x12y3", "x12y30", -1),
        ("A00123:45:HXXX:1:1101:2:1000", "A00123:45:HXXX:1:1101:10:999", -1),
    ],
)
def test_strnum_cmp(first, second, expected):
    """The values returned by samtools' strnum_cmp."""
    assert mta.strnum_cmp(first, second) == expected
    # Antisymmetric.
    assert mta.strnum_cmp(second, first) == -expected


def _random_read_name(rng):
    """A read name mixing digit runs (some with leading zeros), letters of both
    cases and punctuation."""
    parts = []
    for _ in range(rng.randrange(1, 6)):
        kind = rng.random()
        if kind < 0.4:
            parts.append("0" * rng.randrange(3) + str(rng.randrange(200)))
        elif kind < 0.7:
            parts.append("".join(rng.choice(string.ascii_letters) for _ in range(rng.randrange(1, 3))))
        else:
            parts.append(rng.choice("_-.:/|#"))
    return "".join(parts)


def test_strnum_cmp_sorts_as_samtools(tmp_path):
    """Sorting with strnum_cmp gives the order of pysam.sort("-n"), ties
    (r7, r07) included: both keep them in their input order."""
    rng = random.Random(42)
    names = [_random_read_name(rng) for _ in range(2000)]
    names += [f"SRR123.{i}" for i in rng.sample(range(5000), 200)]
    names += [f"A00123:45:HXXX:1:{rng.choice([1101, 2101])}:{rng.randrange(5000)}:{rng.randrange(5000)}" for _ in range(200)]
    names += ["r7", "r07", "r007", "r0", "r00", "r", "r1", "r_1", "r.5", "r5", "rAb", "rAB"] * 2
    rng.shuffle(names)
    unsorted_bam = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_bam), "wb", header={"HD": {"VN": "1.6"}}) as bam:
        for name in names:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "ACGT"
            read.query_qualities = pysam.qualitystring_to_array("IIII")
            read.flag = 4
            bam.write(read)
    sorted_bam = tmp_path / "sorted.bam"
    pysam.sort("-n", "-o", str(sorted_bam), str(unsorted_bam))
    with pysam.AlignmentFile(str(sorted_bam), "rb", check_sq=False) as bam:
        samtools_order = [read.query_name for read in bam.fetch(until_eof=True)]
    assert samtools_order == sorted(names, key=functools.cmp_to_key(mta.strnum_cmp))
