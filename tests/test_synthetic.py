# Tests of the mock metagenome generator itself (tests/synthetic.py): the
# other synthetic-data tests rely on these properties.

import filecmp

import numpy as np
import pandas as pd
import pytest

from tests import synthetic


def test_same_seed_same_files(tmp_path):
    first = synthetic.make_metagenome(seed=7).write_all(tmp_path / "a")
    second = synthetic.make_metagenome(seed=7).write_all(tmp_path / "b")
    third = synthetic.make_metagenome(seed=8).write_all(tmp_path / "c")
    for key in first:
        assert filecmp.cmp(first[key], second[key], shallow=False)
    assert not filecmp.cmp(first["pairs"], third["pairs"], shallow=False)


def test_community(metagenome):
    truth = metagenome.truth
    genomes = truth[truth.Kind == "genome"]
    # Distinct GC contents, from 30% to 70%.
    gc = genomes.groupby("Genome")["GC"].mean()
    assert np.allclose(sorted(gc), [30, 40, 50, 60, 70], atol=1.5)
    # Contigs tile each molecule, with lengths of roughly 2-40 kb.
    for molecule in metagenome.molecules:
        assert sum(c.length for c in molecule.contigs) == molecule.length
    assert genomes.Length.min() >= 2_000
    assert genomes.Length.max() <= 42_000
    # Both strands are represented, and contigs are shuffled in the assembly.
    assert truth.Reverse.any() and not truth.Reverse.all()
    assert list(truth.Genome) != sorted(truth.Genome)
    # Every MGE travels with a host genome.
    mges = truth[truth.Kind == "mge"]
    assert set(mges.Host) <= set(genomes.Genome)


def test_contig_sequences_match_their_molecule(metagenome):
    for molecule in metagenome.molecules:
        for contig in molecule.contigs:
            piece = molecule.seq[contig.start : contig.end]
            expected = synthetic.reverse_complement(piece) if contig.reverse else piece
            assert contig.seq == expected


def test_files(metagenome, metagenome_files):
    names = [c.name for c in metagenome.contigs]
    fasta = synthetic.read_fasta(metagenome_files["assembly"])
    assert list(fasta) == names
    depth = pd.read_csv(metagenome_files["depth"], sep="\t")
    assert list(depth.columns[:3]) == ["contigName", "contigLen", "totalAvgDepth"]
    assert list(depth.contigName) == names
    expected = 10 * metagenome.truth.Abundance.to_numpy()
    assert np.allclose(depth.totalAvgDepth, expected, rtol=0.25)
    with open(metagenome_files["pairs"]) as pairs:
        header = [line.rstrip("\n") for line in pairs if line.startswith("#")]
    assert header[:3] == [
        "## pairs format v1.0",
        "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2",
        "#sorted: readID",
    ]
    assert [line.split()[1] for line in header if line.startswith("#chromsize")] == names


def test_contacts(metagenome):
    """Contacts decay with the distance and stay mostly inside genomes."""
    pairs = metagenome.pairs
    truth = metagenome.truth
    sizes = truth.Length
    assert (pairs.pos1 >= 1).all() and (pairs.pos1 <= pairs.chr1.map(sizes)).all()
    assert (pairs.pos2 >= 1).all() and (pairs.pos2 <= pairs.chr2.map(sizes)).all()
    # Upper triangle, as written by metator.align.
    rank = {name: i for i, name in enumerate(truth.Name)}
    r1, r2 = pairs.chr1.map(rank), pairs.chr2.map(rank)
    assert ((r1 < r2) | ((r1 == r2) & (pairs.pos1 <= pairs.pos2))).all()
    # Decay: inside contigs, short distances outnumber long ones.
    intra = pairs[pairs.chr1 == pairs.chr2]
    distance = intra.pos2 - intra.pos1
    assert (distance < 1_000).sum() > 3 * ((distance >= 1_000) & (distance < 2_000)).sum()
    # Few contacts between genomes, but MGEs touch their host.
    host = truth.Host
    same_host = pairs.chr1.map(host) == pairs.chr2.map(host)
    assert 0 < (~same_host).mean() < 0.03
    molecule = truth.Genome
    mge_pairs = pairs[(pairs.chr1.map(truth.Kind) == "mge") ^ (pairs.chr2.map(truth.Kind) == "mge")]
    mge_pairs = mge_pairs[mge_pairs.chr1.map(molecule) != mge_pairs.chr2.map(molecule)]
    assert (mge_pairs.chr1.map(host) == mge_pairs.chr2.map(host)).mean() > 0.9


def test_linear_genome_ends_do_not_touch():
    """On a linear genome the two ends contact each other least."""
    spec = synthetic.GenomeSpec("line", gc=0.5, abundance=1, length=100_000, circular=False)
    line = synthetic.make_metagenome(seed=1, genomes=(spec,), mges=(), n_pairs=20_000, p_noise=0, intra_uniform=0, shuffle=False)
    first, last = line.contigs[0].name, line.contigs[-1].name
    pairs = line.pairs
    ends = ((pairs.chr1 == first) & (pairs.chr2 == last)).sum()
    neighbours = ((pairs.chr1 == first) & (pairs.chr2 == line.contigs[1].name)).sum()
    assert ends < neighbours / 5


@pytest.mark.parametrize("strands", ["forward", "alternate"])
def test_strand_modes(strands):
    community = synthetic.make_metagenome(seed=2, n_pairs=100, strands=strands, shuffle=False)
    reverse = [c.reverse for m in community.molecules for c in m.contigs]
    if strands == "forward":
        assert not any(reverse)
    else:
        assert reverse[:4] == [False, True, False, True]
    with pytest.raises(ValueError, match="strands"):
        synthetic.make_metagenome(seed=2, n_pairs=100, strands="both")


def test_oracle_quality(metagenome):
    genome_1 = [c.name for c in metagenome.molecule("genome_1").contigs]
    genome_2 = [c.name for c in metagenome.molecule("genome_2").contigs]
    phage = [c.name for c in metagenome.molecule("phage_1").contigs]
    assert synthetic.bin_quality(genome_1, metagenome) == (1.0, 1.0)
    assert synthetic.bin_quality(genome_1 + genome_2, metagenome) == (1.0, 2.0)
    assert synthetic.bin_quality(phage, metagenome) == (0.0, 1.0)
    completeness, redundancy = synthetic.bin_quality(genome_1[:3], metagenome)
    assert 0 < completeness < 1 and redundancy == 1.0


def test_fake_micomplete_writes_the_merged_format(metagenome, tmp_path):
    import metator.io as mio

    bins = tmp_path / "bins"
    bins.mkdir()
    metagenome.write_fasta(bins / "b1.fa", [c.name for c in metagenome.molecule("genome_3").contigs])
    (bins / "notes.txt").write_text("not a bin\n")
    calls = []
    synthetic.fake_micomplete(metagenome, calls)(str(bins), str(tmp_path / "out.tsv"), 1)
    assert calls == [(str(bins), str(tmp_path / "out.tsv"))]
    assert (tmp_path / "out.tsv").read_text().startswith("## miComplete\n")
    summary = mio.micomplete_results_to_dict(str(tmp_path / "out.tsv"))
    assert list(summary) == ["b1"]
    assert summary["b1"]["Weighted completeness"] == 1.0
    assert summary["b1"]["Length"] == metagenome.molecule("genome_3").length


def test_comparison_helpers():
    assert synthetic.adjusted_rand([1, 1, 2, 2], ["a", "a", "b", "b"]) == pytest.approx(1)
    assert synthetic.adjusted_rand([1, 1, 2, 2], ["a", "b", "a", "b"]) < 0
    assert synthetic.purity(["x", "x", "y", "y"], [1, 1, 1, 2]) == 0.75
