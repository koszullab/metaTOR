# Scaffolding tests: contigs cut from one synthetic genome, with a known
# order and orientation, are put back together from their Hi-C contacts.

import shutil

import numpy as np
import pyfastx
import pytest

import metator.commands as mtc
import metator.io as mio
import metator.scaffold as mts
from metator.regions import Bin
from tests import synthetic

pytestmark = pytest.mark.skipif(shutil.which("pairtools") is None, reason="pairtools is not installed")


def _genome(seed, circular, length=100_000, strands="forward", shuffle=True):
    """One genome cut into 12-25 kb contigs, with a steep contact decay
    (P(s) ~ s^-1.5, no long-range plateau) so that the two ends of a 100 kb
    linear genome touch as little as those of a Mb-long one."""
    spec = synthetic.GenomeSpec("genome", gc=0.45, abundance=1.0, length=length, circular=circular)
    return synthetic.make_metagenome(
        seed,
        genomes=(spec,),
        mges=(),
        n_pairs=40_000,
        contig_median=16_000,
        contig_bounds=(12_000, 25_000),
        p_noise=0.01,
        decay=1.5,
        intra_uniform=0,
        strands=strands,
        shuffle=shuffle,
        prefix="frag",
    )


def _inputs(community, directory):
    """Assembly and sorted, indexed pairs of a community."""
    directory.mkdir()
    files = community.write_all(directory)
    pairs = mio.sort_pairs_pairtools(str(files["pairs"]), tmp_dir=str(directory))
    return files["assembly"], pairs


def _canonical(order, circular):
    """An order up to its reading direction (and rotation if circular)."""
    candidates = []
    for oriented in (list(order), list(order)[::-1]):
        shifts = range(len(oriented)) if circular else [0]
        candidates += [tuple(oriented[k:] + oriented[:k]) for k in shifts]
    return min(candidates)


def _true_order(community):
    return [contig.name for contig in community.molecules[0].contigs]


def _scaffold(assembly, pairs, window_size=5000, threshold=0.05):
    scaffolder = mts.Scaffolder(
        Bin("bin", pyfastx.Fasta(str(assembly))), [mio.get_pairs_data(str(pairs))], window_size, threshold
    )
    return scaffolder.solve_graph()


@pytest.mark.parametrize("circular", [False, True], ids=["linear", "circular"])
def test_scaffold_order(tmp_path, circular):
    community = _genome(seed=4, circular=circular)
    scaffolds = _scaffold(*_inputs(community, tmp_path / "inputs"))
    assert len(scaffolds) == 1
    assert scaffolds[0].circular == circular
    order = [contig.name for contig in scaffolds[0].contigs]
    assert _canonical(order, circular) == _canonical(_true_order(community), circular)


def test_scaffold_single_contig(tmp_path):
    community = _genome(seed=5, circular=False, length=20_000)
    assert len(community.contigs) == 1
    # A window larger than half the contig: no intra signal to normalise
    # with, the default normalisation is used.
    scaffolds = _scaffold(*_inputs(community, tmp_path / "inputs"), window_size=12_000)
    assert len(scaffolds) == 1
    assert [repr(c) for c in scaffolds[0].contigs] == ["frag_001:0-20000:+"]
    assert not scaffolds[0].circular


def test_scaffolds_partition_the_bin(tmp_path):
    """Whatever the window, every contig ends in exactly one scaffold."""
    community = _genome(seed=6, circular=False)
    scaffolds = _scaffold(*_inputs(community, tmp_path / "inputs"), window_size=13_000, threshold=0.5)
    names = [contig.name for scaffold in scaffolds for contig in scaffold.contigs]
    assert sorted(names) == sorted(c.name for c in community.contigs)
    assert len({scaffold.name for scaffold in scaffolds}) == len(scaffolds)


def _assert_scaffold_is_the_genome(fasta, community, junction):
    """The scaffold sequence is the genome, read on either strand and, the
    genome being circular, from any origin."""
    sequences = synthetic.read_fasta(fasta)
    assert len(sequences) == 1
    (seq,) = sequences.values()
    assert seq.count(junction) == len(community.contigs) - 1
    seq = seq.replace(junction, "")
    genome = community.molecules[0].seq
    assert seq in genome + genome or seq in synthetic.reverse_complement(genome + genome)


def test_scaffold_command(tmp_path):
    community = _genome(seed=7, circular=True)
    assembly, pairs = _inputs(community, tmp_path / "inputs")
    out_fasta, out_frags = tmp_path / "scaffold.fa", tmp_path / "frags.txt"
    args = ["-b", "bin", "-i", str(assembly), "-o", str(out_fasta), "-O", str(out_frags), "-j", "NNNNNNNNNN", str(pairs)]
    mtc.Scaffold(args, {}).execute()
    _assert_scaffold_is_the_genome(out_fasta, community, "NNNNNNNNNN")
    with open(out_frags) as frags:
        lines = frags.read().splitlines()
    assert lines[0] == f">metator_scaffold_0001 circular length: {community.molecules[0].length}"
    rows = [line.split("\t") for line in lines[1:]]
    assert sorted(row[0] for row in rows) == sorted(c.name for c in community.contigs)
    assert all(row[1] == "0" and row[3] in ("1", "-1") for row in rows)


def test_parallel_scaffold_default_outputs(tmp_path):
    """The pipeline worker reads <bin>.fa and writes <bin>_scaffolded.fa and
    <bin>_info_frags.txt."""
    community = _genome(seed=8, circular=True)
    assembly, pairs = _inputs(community, tmp_path / "inputs")
    final_dir, out_dir = tmp_path / "final_bin", tmp_path / "scaffold_bin"
    final_dir.mkdir()
    out_dir.mkdir()
    shutil.copyfile(assembly, final_dir / "bin_1.fa")
    mts.parallel_scaffold("bin_1", str(final_dir), [str(pairs)], str(out_dir), "NNNNN")
    _assert_scaffold_is_the_genome(out_dir / "bin_1_scaffolded.fa", community, "NNNNN")
    assert (out_dir / "bin_1_info_frags.txt").stat().st_size > 0


def test_scaffold_joins_contig_starts(tmp_path):
    """Contigs alternately on the forward and reverse strand, in genome
    order: the junction of contigs 2 and 3 joins their starts."""
    community = _genome(seed=9, circular=True, strands="alternate", shuffle=False)
    assembly, pairs = _inputs(community, tmp_path / "inputs")
    out_fasta = tmp_path / "scaffold.fa"
    mts.get_scaffolds("bin", str(assembly), [str(pairs)], str(out_fasta), str(tmp_path / "frags.txt"), junctions="NNNNN")
    _assert_scaffold_is_the_genome(out_fasta, community, "NNNNN")


def test_scaffold_joins_contig_end_to_next_start(tmp_path):
    """Contigs on the forward strand, in genome order: each junction joins
    the end of a contig to the start of the next one in the fasta."""
    community = _genome(seed=9, circular=True, shuffle=False)
    assembly, pairs = _inputs(community, tmp_path / "inputs")
    out_fasta = tmp_path / "scaffold.fa"
    mts.get_scaffolds("bin", str(assembly), [str(pairs)], str(out_fasta), str(tmp_path / "frags.txt"), junctions="NNNNN")
    _assert_scaffold_is_the_genome(out_fasta, community, "NNNNN")


def test_normalize_matrix_scores_every_junction(tmp_path):
    """Contig k has the nodes 2k (start) and 2k + 1 (end): every pair of
    nodes is scored but the two ends of one contig."""
    (tmp_path / "bin.fa").write_text("".join(f">c{k}\n{'ACGT' * 10}\n" for k in range(4)))
    scaffolder = mts.Scaffolder(Bin("bin", pyfastx.Fasta(str(tmp_path / "bin.fa"))), [])
    scaffolder.score_matrice = np.ones((8, 8))
    scaffolder.normalize = np.ones(8)
    contig = np.arange(8) // 2
    assert np.array_equal(scaffolder.normalize_matrix() > 0, contig[:, None] != contig[None, :])


def _oriented(contigs):
    return [(contig.name, contig.is_reverse) for contig in contigs]


def _canonical_oriented(oriented, circular):
    """Oriented contigs up to the reading direction (and rotation if
    circular): read backwards, every contig changes strand."""
    flipped = [(name, not reverse) for name, reverse in oriented[::-1]]
    candidates = []
    for order in (list(oriented), flipped):
        shifts = range(len(order)) if circular else [0]
        candidates += [tuple(order[k:] + order[:k]) for k in shifts]
    return min(candidates)


def test_scaffold_terminal_contig_keeps_its_orientation(tmp_path):
    """Contigs A, B, C of a linear genome. The outer end of A (node 0) has
    no partner and contacts the start of B (node 2) only through A; its
    single association must not take node 2 from the end of A (node 1)."""
    (tmp_path / "bin.fa").write_text("".join(f">{name}\n{'ACGT' * 10}\n" for name in "ABC"))
    scaffolder = mts.Scaffolder(Bin("bin", pyfastx.Fasta(str(tmp_path / "bin.fa"))), [], threshold=0.05)
    scores = np.zeros((6, 6))
    scores[1, 2], scores[3, 4] = 2.0, 2.1  # the junctions A-B and B-C
    scores[0, 2] = 0.3  # the start of A sees the start of B through A
    scores[1, 5], scores[0, 5] = 0.06, 0.01  # weak contacts A-C
    scaffolder.score_matrice = scores
    scaffolder.normalize = np.ones(6)
    (scaffold,) = scaffolder.solve_graph()
    assert not scaffold.circular
    expected = [("A", False), ("B", False), ("C", False)]
    assert _canonical_oriented(_oriented(scaffold.contigs), False) == _canonical_oriented(expected, False)


@pytest.mark.parametrize("seed", [5, 85, 122])
def test_scaffold_linear_genome_random_strands(tmp_path, seed):
    """Seeds where the outer end of a terminal contig used to take the
    partner of its other end, flipping the contig."""
    community = _genome(seed=seed, circular=False, strands="random")
    scaffolds = _scaffold(*_inputs(community, tmp_path / "inputs"))
    assert len(scaffolds) == 1
    assert not scaffolds[0].circular
    truth = [(contig.name, contig.reverse) for contig in community.molecules[0].contigs]
    assert _canonical_oriented(_oriented(scaffolds[0].contigs), False) == _canonical_oriented(truth, False)
