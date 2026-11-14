# Tests of the genomic region structures used by the scaffolder.

import numpy as np
import pyfastx
import pytest

from metator.regions import Bin, BreakPoint, Contig, Position, Scaffold
from tests import synthetic

SEQS = {"c1": 100, "c2": 50, "c3": 80}


@pytest.fixture
def bin_fasta(tmp_path):
    """A bin of three contigs (100, 50 and 80 bp) and their sequences."""
    rng = np.random.default_rng(3)
    sequences = {name: synthetic.random_sequence(rng, length, 0.5) for name, length in SEQS.items()}
    with open(tmp_path / "bin.fa", "w") as fasta:
        for name, seq in sequences.items():
            fasta.write(f">{name}\n{seq}\n")
    return pyfastx.Fasta(str(tmp_path / "bin.fa")), sequences


def _layout(contigs):
    return [repr(contig) for contig in contigs]


# -- Position ----------------------------------------------------------------


def test_position():
    assert repr(Position("c", 10, True)) == "c:10:+"
    assert repr(Position("c", 10, False)) == "c:10:-"
    assert repr(Position("c", 10)) == "c:10:"
    assert Position("c", 10, True).has_sign()
    assert not Position("c", 10).has_sign()
    # Ordered by contig then coordinate; the sign is not compared.
    assert Position("a", 5) < Position("a", 7) < Position("b", 1)
    assert Position("a", 5, True) == Position("a", 5, False)


# -- Contig ------------------------------------------------------------------


def test_contig():
    contig = Contig("c", 10, 110)
    assert len(contig) == 100
    assert repr(contig) == "c:10-110:+"
    assert contig.middle() == 60
    assert hash(contig) == hash(Contig("c", 10, 110, bin="other"))
    contig.flip()
    assert contig.is_reverse and repr(contig) == "c:10-110:-"
    with pytest.raises(ValueError, match="smaller"):
        Contig("c", 10, 5)
    with pytest.raises(ValueError, match="negative"):
        Contig("c", -1, 5)


def test_contig_intersect():
    assert Contig("c", 0, 50).intersect(Contig("c", 40, 100)) == 10
    assert Contig("c", 40, 100).intersect(Contig("c", 0, 50)) == 10
    assert Contig("c", 0, 50).intersect(Contig("d", 0, 50)) == 0
    assert Contig("c", 0, 50).intersect(Contig("c", 60, 100)) == 0


def test_contig_split():
    left, right = Contig("c", 10, 110).split(30)
    assert _layout([left, right]) == ["c:10-40:+", "c:40-110:+"]
    # On the reverse strand the first 30 bp read are the last 30 of the region.
    left, right = Contig("c", 10, 110, True).split(30)
    assert _layout([left, right]) == ["c:80-110:-", "c:10-80:-"]
    with pytest.raises(ValueError, match="negative"):
        Contig("c", 10, 110).split(-1)
    with pytest.raises(ValueError, match="Beyond"):
        Contig("c", 10, 110).split(101)


# -- BreakPoint --------------------------------------------------------------


def _breakpoint(left, right):
    """Breakpoint from the 3' end of ``left`` to the 5' start of ``right``."""
    bp = BreakPoint(Position(left.name, left.end, True), Position(right.name, right.start, False))
    bp.contig1, bp.contig2 = left, right
    return bp


def test_breakpoint():
    bp = BreakPoint(Position("a", 100, True), Position("b", 0, False))
    assert repr(bp) == "a:100:+|b:0:- UNK"
    assert bp.signs == ("+", "-")
    assert bp.has_signs() and not bp.has_contigs()
    bp.contig1 = Contig("a", 0, 100)
    bp.contig2 = Contig("b", 0, 50)
    assert bp.has_contigs()
    assert repr(bp) == "[a:0-100:+]a:100:+|b:0:-[b:0-50:+] UNK"
    # The contigs must end where the breakpoint is.
    with pytest.raises(ValueError, match="contig1"):
        bp.contig1 = Contig("a", 0, 90)
    with pytest.raises(ValueError, match="contig2"):
        bp.contig2 = Contig("other", 0, 50)
    # On the reverse strand, the 3' end is the start.
    reverse = BreakPoint(Position("a", 0, True), Position("b", 50, False), "INV")
    reverse.contig1 = Contig("a", 0, 100, True)
    reverse.contig2 = Contig("b", 0, 50, True)
    # Without sign, either end fits.
    unsigned = BreakPoint(Position("a", 0), Position("b", 50))
    unsigned.contig1 = Contig("a", 0, 100)
    unsigned.contig2 = Contig("b", 0, 50)
    assert unsigned.signs is None


def test_breakpoint_intersect_and_connect():
    first = _breakpoint(Contig("a", 0, 100), Contig("b", 0, 5000))
    second = _breakpoint(Contig("b", 0, 5000), Contig("c", 0, 100))
    assert first.intersect(second) == 5000
    assert first.can_connect(second, min_intersect=1000)
    assert not first.can_connect(second, min_intersect=10_000)
    assert not second.can_connect(first)
    bare = BreakPoint(Position("b", 5000), Position("c", 0))
    with pytest.raises(ValueError, match="sign"):
        first.can_connect(bare)
    with pytest.raises(ValueError, match="overlap"):
        first.intersect(bare)


# -- Bin ---------------------------------------------------------------------


def test_bin(bin_fasta):
    fasta = bin_fasta[0]
    bin = Bin("bin", fasta)
    assert len(bin) == 230
    assert bin.get_contigs_counts() == 3
    assert _layout(bin.contigs) == ["c1:0-100:+", "c2:0-50:+", "c3:0-80:+"]
    assert all(contig.bin == "bin" for contig in bin.contigs)
    assert list(bin.boundaries) == [0, 100, 150, 230]
    assert bin.boundaries is bin.boundaries  # memorised...
    bin.contigs.pop()
    assert list(bin.boundaries) == [0, 100, 150]  # ...until the contigs change
    assert bin.get_contig_bounds(120) == (1, (100, 150))
    with pytest.raises(ValueError, match="out of bounds"):
        bin.get_contig_bounds(150)


def test_bin_contigs_between(bin_fasta):
    bin = Bin("bin", bin_fasta[0])
    inside = bin.get_contigs_between(10, 60)
    assert [(c.start, c.end) for c in inside] == [(10, 60)]
    across = bin.get_contigs_between(50, 180)
    assert [(c.start, c.end) for c in across] == [(50, 100), (0, 50), (0, 30)]
    for start, end in [(-1, 10), (0, 231)]:
        with pytest.raises(ValueError, match="out of bounds"):
            bin.get_contigs_between(start, end)


@pytest.mark.xfail(
    strict=True, raises=AssertionError, reason="Bin.get_contigs_between names every piece after the bin instead of its contig"
)
def test_bin_contigs_between_keep_contig_names(bin_fasta):
    bin = Bin("bin", bin_fasta[0])
    assert [c.name for c in bin.get_contigs_between(50, 180)] == ["c1", "c2", "c3"]


def test_bin_sequence(bin_fasta):
    fasta, sequences = bin_fasta
    bin = Bin("bin", fasta)
    bin.contigs[1].flip()
    bin.contigs.append(Contig("c3", 5, 5))  # empty: purged
    seq = "".join(bin.get_seq(fasta))
    assert seq == sequences["c1"] + synthetic.reverse_complement(sequences["c2"]) + sequences["c3"]
    assert len(bin.contigs) == 3


def test_bin_single_contig_has_no_breakpoint(tmp_path):
    (tmp_path / "one.fa").write_text(">c1\nACGTACGTAC\n")
    assert list(Bin("bin", pyfastx.Fasta(str(tmp_path / "one.fa"))).get_breakpoints()) == []


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="Bin.get_breakpoints positions breakpoints on the bin name, which the contig setters reject",
)
def test_bin_breakpoints(bin_fasta):
    breakpoints = list(Bin("bin", bin_fasta[0]).get_breakpoints())
    assert [(bp.contig1.name, bp.contig2.name) for bp in breakpoints] == [("c1", "c2"), ("c2", "c3")]


EDITS = {
    "invert across contigs": (
        lambda b: b.invert(50, 120),
        ["c1:0-50:+", "c2:0-20:-", "c1:50-100:-", "c2:20-50:+", "c3:0-80:+"],
    ),
    "invert inside a contig": (
        lambda b: b.invert(10, 20),
        ["c1:0-10:+", "c1:10-20:-", "c1:20-100:+", "c2:0-50:+", "c3:0-80:+"],
    ),
    "delete across contigs": (lambda b: b.delete(50, 120), ["c1:0-50:+", "c2:20-50:+", "c3:0-80:+"]),
    "delete inside a contig": (lambda b: b.delete(10, 20), ["c1:0-10:+", "c1:20-100:+", "c2:0-50:+", "c3:0-80:+"]),
    "insert between contigs": (
        lambda b: b.insert_contigs(100, [Contig("ins", 0, 10)]),
        ["c1:0-100:+", "ins:0-10:+", "c2:0-50:+", "c3:0-80:+"],
    ),
    "insert inside a contig": (
        lambda b: b.insert(30, Contig("ins", 0, 10)),
        ["c1:0-30:+", "ins:0-10:+", "c1:30-100:+", "c2:0-50:+", "c3:0-80:+"],
    ),
    "insert at the end": (
        lambda b: b.insert(230, Contig("ins", 0, 10)),
        ["c1:0-100:+", "c2:0-50:+", "c3:0-80:+", "ins:0-10:+"],
    ),
}


@pytest.mark.xfail(
    strict=True,
    raises=AttributeError,
    reason="Bin.insert/invert/delete append to Bin.breakpoints, a numpy array without append",
)
@pytest.mark.parametrize("edit", list(EDITS))
def test_bin_edits(bin_fasta, edit):
    bin = Bin("bin", bin_fasta[0])
    change, layout = EDITS[edit]
    change(bin)
    assert _layout(bin.contigs) == layout


# -- Scaffold ----------------------------------------------------------------


def test_scaffold(bin_fasta):
    fasta, sequences = bin_fasta
    contigs = [Contig("c1", 0, 100), Contig("c2", 0, 50, True), Contig("c3", 10, 10)]
    scaffold = Scaffold("scaffold_1", contigs, circular=True)
    assert len(scaffold) == 150
    assert scaffold.get_fasta_entry() == ">scaffold_1 circular length: 150"
    assert "".join(scaffold.get_seq(fasta)) == sequences["c1"] + synthetic.reverse_complement(sequences["c2"])
    assert len(scaffold.contigs) == 2  # the empty contig was purged
    assert Scaffold("s", [Contig("c1", 0, 100)]).get_fasta_entry() == ">s linear length: 100"
