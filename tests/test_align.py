# Test for align module

import metator.align as mta
import pysam

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
