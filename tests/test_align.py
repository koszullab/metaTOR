# Test for align module

import metator.align as mta
import pysam
import pytest

assembly = "tests_data/assembly.fa"


def test_align():
    # Not tested here as it's tested in the commands tests.
    ...


def test_get_contact_pairs():
    ...


def test_merge_alignement():
    ...


def test_process_bamfile():
    ...


def _write_bam(path, reads):
    """Write a name-sorted BAM of 50 bp reads: (name, contig index, 0-based start, flag)."""
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "c1", "LN": 10_000}, {"SN": "c2", "LN": 20_000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, ref_id, start, flag in reads:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.flag = flag
            read.reference_id = ref_id
            read.reference_start = start
            read.mapping_quality = 60
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
    contig_data = {"c1": {"id": 1, "length": 10_000}, "c2": {"id": 2, "length": 20_000}}
    out = tmp_path / "alignment.pairs"
    n_pairs = mta.process_bwa_bamfile(str(bam), 30, contig_data, str(out))
    pairs = [line.split() for line in out.read_text().splitlines() if not line.startswith("#")]
    assert n_pairs == 2
    assert pairs == [
        ["r1", "c1", "101", "c2", "5001", "+", "-"],
        ["r2", "c2", "201", "c2", "301", "+", "+"],
    ]
