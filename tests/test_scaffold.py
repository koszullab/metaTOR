# Scaffolding tests: contigs cut from one synthetic genome, with a known
# order and orientation, are put back together from their Hi-C contacts.

import shutil

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


@pytest.mark.xfail(
    strict=True,
    reason="Scaffolder.normalize_matrix keeps scores for j - i > 2 instead of > 1: the junction between the starts of "
    "two contigs adjacent in the fasta is never scored",
)
def test_scaffold_joins_contig_starts(tmp_path):
    """Contigs alternately on the forward and reverse strand, in genome
    order: the junction of contigs 2 and 3 joins their starts."""
    community = _genome(seed=9, circular=True, strands="alternate", shuffle=False)
    assembly, pairs = _inputs(community, tmp_path / "inputs")
    out_fasta = tmp_path / "scaffold.fa"
    mts.get_scaffolds("bin", str(assembly), [str(pairs)], str(out_fasta), str(tmp_path / "frags.txt"), junctions="NNNNN")
    _assert_scaffold_is_the_genome(out_fasta, community, "NNNNN")
