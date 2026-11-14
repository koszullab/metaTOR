# Test for CLI tools of metator, MGE related


import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import metator.commands as mtc
import metator.io as mio
import metator.mge as mtm
from tests import synthetic


# Use global variables for temporary tmp/out folders
@pytest.fixture(scope="session")
def tmp_dir(tmp_path_factory):
    p = tmp_path_factory.mktemp("tmp")
    os.makedirs(p, exist_ok=True)
    return str(p)


global_args = {
    "NETWORK": "tests_data/network.txt",
    "CONTIGS": "tests_data/contig_data_final.txt",
    "BINNING": "tests_data/binning.txt",
    "SUBSET_FASTA": "tests_data/subset.assembly.fa",
    "SUBSET_PAIRS": "tests_data/subset.pairs.gz",
    "SUBSET_MGES": "tests_data/subset.mges.txt",
}


def test_mge(tmp_dir):
    """Test the metator mge command with subset data from pipeline output."""
    args = (
        "-c {CONTIGS} -a {SUBSET_FASTA} -b {BINNING} -m {SUBSET_MGES} "
        "-o {OUT_TEST} --no-clean-up --threshold-bin 0.5 {SUBSET_PAIRS}"
    ).format(
        OUT_TEST=Path(tmp_dir, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator mge", *args.split(" "))
    proc = mtc.Mge(args.split(" "), {})
    proc.execute()


def test_mge_random(tmp_path):
    """--random shuffles the MGE bins and still writes the bins and summary."""
    args = (
        "-c {CONTIGS} -a {SUBSET_FASTA} -b {BINNING} -m {SUBSET_MGES} -o {OUT_TEST} --random --threshold-bin 0.5 {SUBSET_PAIRS}"
    ).format(OUT_TEST=tmp_path / "out", **global_args)
    mtc.Mge(args.split(" "), {}).execute()
    assert len(os.listdir(tmp_path / "out")) > 0


# Tests on the mock metagenome (tests/synthetic.py): a phage assembled in two
# contigs and a plasmid in one, each bound to a different host.


def _mge_contigs(metagenome):
    """Names and sizes of the MGE contigs, in assembly order."""
    truth = metagenome.truth
    mges = truth[truth.Kind == "mge"]
    return list(mges.Name), list(mges.Length)


def _expected_matrix(metagenome, names):
    """Contacts between the given contigs, counted from the simulated pairs:
    upper triangle, and intra-contig contacts only beyond 1 kb."""
    pairs = metagenome.pairs
    mat = np.zeros((len(names), len(names)))
    for i, first in enumerate(names):
        for j in range(i, len(names)):
            second = names[j]
            if i == j:
                sel = (pairs.chr1 == first) & (pairs.chr2 == first)
                mat[i, i] = (np.abs(pairs.pos1[sel] - pairs.pos2[sel]) > 1000).sum()
            else:
                sel = ((pairs.chr1 == first) & (pairs.chr2 == second)) | ((pairs.chr1 == second) & (pairs.chr2 == first))
                mat[i, j] = sel.sum()
    return mat


def test_build_matrix(metagenome, metagenome_files, metagenome_indexed_pairs, tmp_path):
    names, sizes = _mge_contigs(metagenome)
    expected = _expected_matrix(metagenome, names)
    assert np.array_equal(mtm.build_matrix(names, sizes, [str(metagenome_indexed_pairs)]), expected)
    # Without pypairix index the pairs file is read line by line.
    plain = tmp_path / "plain.pairs"
    shutil.copyfile(metagenome_files["pairs"], plain)
    assert np.array_equal(mtm.build_matrix(names, sizes, [str(plain)]), expected)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="without a pypairix index, build_matrix stores the contacts of a pair whose chr1 comes after chr2 in the "
    "contig list in the lower triangle, which resolve_matrix ignores",
)
def test_build_matrix_contig_order(metagenome, metagenome_files, metagenome_indexed_pairs, tmp_path):
    names, sizes = _mge_contigs(metagenome)
    names, sizes = names[::-1], sizes[::-1]
    plain = tmp_path / "plain.pairs"
    shutil.copyfile(metagenome_files["pairs"], plain)
    indexed = mtm.build_matrix(names, sizes, [str(metagenome_indexed_pairs)])
    assert np.array_equal(mtm.build_matrix(names, sizes, [str(plain)]), indexed)


def test_resolve_matrix():
    mat = np.array([[100, 90, 0], [0, 100, 5], [0, 0, 100]], dtype=float)
    assert mtm.resolve_matrix(mat, 0.8) == [[0, 1, pytest.approx(0.9)]]
    assert mtm.resolve_matrix(np.array([[100, 90], [0, 100]], dtype=float), 0.95) == []
    # Ties: the first maximum is merged first.
    tied = np.array([[100, 90, 0, 0], [0, 100, 0, 0], [0, 0, 100, 90], [0, 0, 0, 100]], dtype=float)
    assert [b[:2] for b in mtm.resolve_matrix(tied, 0.8)] == [[0, 1], [2, 3]]


@pytest.mark.parametrize(
    "pairs, bins, score",
    [
        # Two new bins fused by a third association.
        ([[0, 1, 0.9], [2, 3, 0.85], [1, 2, 0.82]], [["c0", "c1", "c2", "c3"], ["c4"]], 0.82),
        # The same, the second bin coming first.
        ([[0, 1, 0.9], [2, 3, 0.85], [3, 0, 0.8]], [["c0", "c1", "c2", "c3"], ["c4"]], 0.8),
        # A contig joining an existing bin.
        ([[0, 1, 0.9], [1, 2, 0.85]], [["c0", "c1", "c2"], ["c3"], ["c4"]], 0.85),
    ],
)
def test_update_mge_data(pairs, bins, score):
    data = pd.DataFrame({"Name": [f"c{i}" for i in range(5)]}, index=[10, 11, 12, 13, 14])
    data, mge_bins = mtm.update_mge_data(data, pairs)
    assert [mge_bins[k]["Contigs"] for k in sorted(mge_bins)] == bins
    assert mge_bins[1]["Score"] == pytest.approx(score)
    assert all(np.isnan(mge_bins[k]["Score"]) for k in sorted(mge_bins)[1:] if len(mge_bins[k]["Contigs"]) == 1)
    assert list(data.MetaTOR_MGE_bin) == [k for k in sorted(mge_bins) for _ in mge_bins[k]["Contigs"]]


def test_shuffle_mge_bins():
    np.random.seed(0)
    data = pd.DataFrame({"Name": ["a", "b", "c", "d"], "MetaTOR_MGE_bin": [1, 1, 2, 3]})
    data, mge_bins = mtm.shuffle_mge_bins(data)
    assert sorted(len(b["Contigs"]) for b in mge_bins.values()) == [1, 1, 2]
    assert sorted(c for b in mge_bins.values() for c in b["Contigs"]) == ["a", "b", "c", "d"]


def _mge_binning(metagenome, metagenome_inputs, metagenome_network, pairs, tmp_path, random=False):
    names, _ = _mge_contigs(metagenome)
    binning = tmp_path / "binning.txt"
    truth = metagenome.truth
    genomes = truth[truth.Kind == "genome"]
    genomes[["Name", "Genome"]].to_csv(binning, sep="\t", header=False, index=False)
    contigs_data, mges_list_id = mio.import_contig_data_mges(
        str(metagenome_network["contig_data"]), mio.import_anvio_binning(str(binning)), names
    )
    out, tmp = tmp_path / "out", tmp_path / "tmp"
    out.mkdir()
    tmp.mkdir()
    summary = mtm.mge_binning(
        str(metagenome_inputs["assembly"]), contigs_data, mges_list_id, str(out), [str(pairs)], str(tmp), 0.8, random
    )
    return summary, out


def test_mge_binning(metagenome, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path):
    """The two contigs of the phage form one MGE bin, the plasmid another."""
    summary, out = _mge_binning(metagenome, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path)
    found = {frozenset(contigs.split(",")) for contigs in summary.Contigs}
    assert found == {frozenset(c.name for c in mge.contigs) for mge in metagenome.mges}
    lengths = dict(zip(summary.Contigs.map(lambda c: frozenset(c.split(","))), summary.BinLength))
    for mge in metagenome.mges:
        assert lengths[frozenset(c.name for c in mge.contigs)] == mge.length
    assert summary.BinningScore.notna().sum() == 1  # only the phage was associated
    written = pd.read_csv(out / "mges_bin_summary.tsv", sep="\t")
    assert list(written.BinName) == list(summary.BinName)
    # One fasta entry per MGE bin, its contigs separated by N spacers.
    sequences = synthetic.read_fasta(out / "mges_binned.fa")
    assert sorted(sequences) == sorted(summary.BinName)
    contig_seqs = {c.name: c.seq for c in metagenome.contigs}
    for name, contigs in zip(summary.BinName, summary.Contigs):
        pieces = [piece for piece in sequences[name].split("N") if piece]
        assert sorted(pieces) == sorted(contig_seqs[c] for c in contigs.split(","))


def test_mge_binning_random(metagenome, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path):
    """--random shuffles the contigs between MGE bins of the same sizes."""
    np.random.seed(1)
    summary, _ = _mge_binning(metagenome, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path, random=True)
    assert sorted(summary.ContigsNumber) == sorted(len(mge.contigs) for mge in metagenome.mges)
    names, _ = _mge_contigs(metagenome)
    assert sorted(",".join(summary.Contigs).split(",")) == sorted(names)
