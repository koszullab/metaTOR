# Tests of metator mge: the MGE contigs of the map of the clusters of a run,
# clustered from the contacts of the map, and stored in its clustering
# rounds. The mock metagenome of tests/synthetic.py: a phage assembled in two
# contigs and a plasmid in one, each bound to a different host.

import json
import os
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import metator.io as mio
import metator.mge as mtm
from tests import synthetic
from tests.hdf5 import walk
from tests.logs import logs


def _mge_contigs(metagenome):
    """Names of the MGE contigs of the mock metagenome, in assembly order."""
    truth = metagenome.truth
    return list(truth.Name[truth.Kind == "mge"])


def _expected_matrix(metagenome, names):
    """Contacts between the given contigs, counted from the simulated pairs, as
    metator mge counted them in the pairs files: upper triangle, and
    intra-contig contacts only beyond 1 kb."""
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


def _bins(metagenome, binsize=None):
    """The bins of a map of the mock metagenome: one per contig (binsize None), or
    of binsize bp, the last of each contig shorter."""
    truth = metagenome.truth
    if binsize is None:
        return pd.DataFrame({"chrom": truth.Name, "start": 0, "end": truth.Length})
    rows = [
        (name, start, min(start + binsize, length))
        for name, length in zip(truth.Name, truth.Length, strict=True)
        for start in range(0, length, binsize)
    ]
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def _map_of_pairs(metagenome, metagenome_files, path, binsize=None, final=None):
    """The map of the pairs of the mock metagenome, as metacooler create pairs
    makes it: with the contact statistics of its reads (the QC block, at 1001
    bp, with the counts of each contig), one bin per contig (binsize None) or
    of binsize bp, its contigs in the genomes of final (a contig name: its
    genome; none by default)."""
    import metacooler

    metacooler.from_pairs(
        str(path),
        pairs=str(metagenome_files["pairs"]),
        bins=_bins(metagenome, binsize),
        contig_genome_map=final,
        ps=False,
        keep_partial=False,
    )
    return path


@pytest.mark.parametrize("binsize", [None, 1000, 5000], ids=["one bin per contig", "1 kb", "5 kb"])
def test_the_contacts_of_a_map_are_those_of_its_pairs(binsize, metagenome, metagenome_files, tmp_path):
    """The contacts of the MGE contigs that mge reads from the map made from the
    pairs (between two contigs, from its pixels; within a contig, from its
    contact statistics at 1001 bp) are those that it read from the pairs files:
    every pair between two MGE contigs, and, within a contig, the pairs of
    reads more than 1 kb apart."""
    import metacooler

    metacool = _map_of_pairs(metagenome, metagenome_files, tmp_path / "map.metacool", binsize)
    mc = metacooler.MetaCooler(str(metacool))
    assert {key: mc.metadata["qc"][key] for key in ["source", "distance_threshold"]} == {
        "source": "pairs",
        "distance_threshold": 1001,
    }
    names = mc.contigs["name"].astype(str)
    positions = np.flatnonzero(names.isin(set(_mge_contigs(metagenome))).to_numpy())
    found = mtm.map_matrix(str(metacool), len(names), positions, mc.contigs[mtm.INTRA_LONG].to_numpy()[positions])
    expected = _expected_matrix(metagenome, list(names.iloc[positions]))
    np.testing.assert_array_equal(found, expected)
    # Contacts within and between the contigs: the phage is merged.
    assert np.diag(expected).all() and np.triu(expected, 1).any()
    clusters = mtm.resolve_matrix(found.copy(), 0.8)
    assert clusters == mtm.resolve_matrix(expected.copy(), 0.8) and len(clusters) == 1


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


def _clusters_map(metagenome, metagenome_files, directory, hosted=False):
    """The map of the clusters of a binning run of the mock metagenome (see
    clusters_map), and the list of its MGE contigs. hosted: the contigs of the
    first MGE (the phage) are in the final cluster of their host, and those of
    the others (the plasmid) in a final cluster of their own, and the genomes
    have the values that a validation gives them (decimals and text, a
    recursive column of True and False among them), where the MGE contigs are
    in none and the genomes have no values."""
    import metacooler

    import metator.clusters as mtcl
    import metator.metacool as mtmc
    import metator.network as mtn

    truth = metagenome.truth
    clusters = {genome.name: f"metator_{k:05d}_00000" for k, genome in enumerate(metagenome.genomes, start=1)}
    if hosted:
        mges = sorted(set(truth.Genome[truth.Kind == "mge"]))
        clusters |= {mge: f"metator_{len(clusters) + k:05d}_00000" for k, mge in enumerate(mges[1:], start=1)}
        hosts = {mges[0]: clusters[str(truth.Host[truth.Genome == mges[0]].iloc[0])]}
    else:
        hosts = {}
    final = {name: hosts.get(genome, clusters.get(genome, "ND")) for name, genome in zip(truth.Name, truth.Genome, strict=True)}
    metacool = _map_of_pairs(metagenome, metagenome_files, directory / "clusters.metacool", final=final)
    contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
    contacts = mtn.contig_contacts(str(metacool))
    names = metacooler.MetaCooler(str(metacool)).contigs["name"]
    table = mtn.contig_data_table(contig_data, names, mtn.contig_hits(contacts, len(names)))
    mtmc.write_map(str(metacool), str(metacool), network=mtmc.Network(table, "empirical_hit"), force=True)
    members = {cluster: [name for name, held in final.items() if held == cluster] for cluster in clusters.values()}
    rounds = mtcl.Rounds()
    for cluster, contigs in members.items():
        columns = {"Completeness": 0.9, "Lineage": "Bacteria", "Markers": 40, "Circular": True}
        rounds.add(mtcl.OVERLAPPING, cluster, contigs, columns=columns)
    rounds.close_round(mtcl.OVERLAPPING)
    for cluster, contigs in members.items():
        columns = {"Markers": 41, "Circular": False}
        rounds.add(mtcl.FINAL, cluster, contigs, parent=cluster, columns=columns)
    mtmc.store_rounds(str(metacool), rounds.tables())
    if hosted:
        values = directory / "values.txt"
        values.write_text(
            "\tCompleteness\tLineage\trecursive\n"
            + "".join(f"{cluster}\t{0.9 + k / 300}\tBacteria {k}\t{k % 2 == 0}\n" for k, cluster in enumerate(members, start=1))
        )
        metacooler.annotate(str(metacool), genomes_table=str(values), genomes_table_format="1", inplace=True)
    mges = directory / "mges.txt"
    mges.write_text("".join(f"{name}\n" for name in _mge_contigs(metagenome)))
    return metacool, mges


@pytest.fixture(scope="module")
def clusters_map(metagenome, metagenome_files, tmp_path_factory):
    """The .metacool map of the clusters of a binning run of the mock metagenome,
    from its pairs, one bin per contig, as metator validation writes it: one
    genome per genome of the community, its final clusters, the MGE contigs in
    none; the contact statistics of its reads (the QC block, from the pairs, at
    1001 bp); the network that metator network stores (GC_content, Hit, the
    weights of empirical_hit); and the clustering rounds of the run,
    overlapping and final (with columns of text, decimals, whole numbers and
    booleans)."""
    return _clusters_map(metagenome, metagenome_files, tmp_path_factory.mktemp("clusters_map"))


@pytest.fixture(scope="module")
def hosted_map(metagenome, metagenome_files, tmp_path_factory):
    """clusters_map, but the phage in the final cluster of its host, the plasmid
    in a final cluster of its own, and the genomes with values (see
    _clusters_map)."""
    return _clusters_map(metagenome, metagenome_files, tmp_path_factory.mktemp("hosted_map"), hosted=True)


def _run_on(clusters, metagenome_inputs, tmp_path):
    """A private copy of a map of clusters (see _clusters_map), and the arguments
    of an mge run on it but the map, its outputs in tmp_path/out."""
    source, mges = clusters
    metacool = tmp_path / "clusters.metacool"
    shutil.copyfile(source, metacool)
    args = ["-a", metagenome_inputs["assembly"], "-m", mges, "-o", tmp_path / "out", "-T", tmp_path / "tmp"]
    return metacool, [str(arg) for arg in args]


@pytest.fixture
def mge_run(clusters_map, metagenome_inputs, tmp_path):
    """A private copy of clusters_map, and the arguments of an mge run on it but
    the map, its outputs in tmp_path/out."""
    return _run_on(clusters_map, metagenome_inputs, tmp_path)


@pytest.fixture
def hosted_run(hosted_map, metagenome_inputs, tmp_path):
    """A private copy of hosted_map, and the arguments of an mge run on it but the
    map, its outputs in tmp_path/out."""
    return _run_on(hosted_map, metagenome_inputs, tmp_path)


def test_mge_binning(metagenome, clusters_map):
    """The two contigs of the phage form one MGE cluster, the plasmid another."""
    import metacooler

    metacool, _ = clusters_map
    summary, bins = mtm.mge_binning(str(metacool), _mge_contigs(metagenome))
    found = {frozenset(contigs.split(",")) for contigs in summary.Contigs}
    assert found == {frozenset(c.name for c in mge.contigs) for mge in metagenome.mges}
    assert list(summary.BinName) == [f"MetaTOR_MGE_{key:05d}" for key in bins]
    contigs = metacooler.MetaCooler(str(metacool)).contigs.set_index("name")
    assert list(summary.columns) == ["BinName", "BinLength", "GC", "Hit", "BinningScore", "ContigsNumber", "Contigs"]
    for row in summary.to_dict("records"):
        members = str(row["Contigs"]).split(",")
        expected = [str(c) for c in bins[int(str(row["BinName"])[-5:])]["Contigs"]]
        assert members == sorted(members, key=lambda name: (-contigs.loc[name, "length"], name)) == expected
        assert row["BinLength"] == contigs.loc[members, "length"].sum() and row["ContigsNumber"] == len(members)
        assert row["Hit"] == contigs.loc[members, "Hit"].sum()
        weighted = (contigs.loc[members, "GC_content"] * contigs.loc[members, "length"]).sum()
        assert row["GC"] == pytest.approx(weighted / row["BinLength"])
    assert summary.BinningScore.notna().sum() == 1  # only the phage was associated
    assert summary.BinningScore.max() > 0.8


def _tie_community(directory):
    """A small community whose MGE scores tie: MGE contigs A (8 kb), B (6 kb) and
    C (4 kb), with 4, 1 and 16 contacts within them between reads more than 1
    kb apart, 2 between A and B and 8 between A and C: the scores of A-B and
    A-C are both 1 (2 / sqrt(4 x 1), 8 / sqrt(4 x 16))."""
    rng = np.random.default_rng(1)
    lengths = {"G1": 60_000, "G2": 50_000, "A": 8_000, "B": 6_000, "C": 4_000}
    assembly = directory / "assembly.fa"
    with open(assembly, "w") as fasta:
        for name, length in lengths.items():
            seq = "".join(rng.choice(list("ACGT"), size=length))
            fasta.write(f">{name}\n" + "".join(seq[i : i + 60] + "\n" for i in range(0, length, 60)))
    pairs = []
    for genome in ("G1", "G2"):
        for _ in range(200):
            start = int(rng.integers(1, lengths[genome] - 3000))
            pairs.append((genome, start, genome, start + int(rng.integers(10, 3000))))
    counts = [("A", 100, "A", 600, 5), ("B", 100, "B", 600, 5), ("C", 100, "C", 600, 5), ("A", 100, "A", 5000, 4),
              ("B", 100, "B", 5000, 1), ("C", 100, "C", 3000, 16), ("A", 2000, "B", 2000, 2), ("A", 3000, "C", 1000, 8),
              ("G1", 1000, "A", 1500, 30), ("G1", 2000, "C", 1500, 30), ("G2", 1000, "B", 1500, 30),
              ("G1", 5000, "G2", 5000, 3)]  # fmt: skip
    for first, pos1, second, pos2, n in counts:
        pairs += [(first, pos1, second, pos2)] * n
    path = directory / "community.pairs"
    order = {name: i for i, name in enumerate(lengths)}
    with open(path, "w") as handle:
        handle.write("## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        for i, (first, pos1, second, pos2) in enumerate(pairs):
            if order[first] > order[second]:
                first, pos1, second, pos2 = second, pos2, first, pos1
            handle.write(f"r{i}\t{first}\t{pos1}\t{second}\t{pos2}\t+\t-\n")
    return assembly, path, lengths


def _network_map(path, pairs, contig_data, bins, final=None):
    """The map of pairs, its network stored as metator network stores it, its
    contigs in the genomes of final (none by default): grouped by genome, then
    the longest first, as the map of the clusters of a run."""
    import metacooler

    import metator.metacool as mtmc
    import metator.network as mtn

    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=final, ps=False, keep_partial=False)
    contacts = mtn.contig_contacts(str(path))
    names = metacooler.MetaCooler(str(path)).contigs["name"]
    table = mtn.contig_data_table(contig_data, names, mtn.contig_hits(contacts, len(names)))
    mtmc.write_map(str(path), str(path), network=mtmc.Network(table, "empirical_hit"), force=True)
    return path


def _clusters(summary):
    """The MGE clusters of a summary, by name: their contigs and score."""
    return {
        name: (contigs, None if pd.isna(score) else round(float(score), 12))
        for name, contigs, score in zip(summary.BinName, summary.Contigs, summary.BinningScore, strict=True)
    }


def test_the_order_of_the_map_changes_no_tie(tmp_path):
    """The MGE contigs are taken the longest first, ties by name, whatever the
    order of the map: a tie of scores is broken as on the map of the contigs,
    longest first (as hicstuff and metacooler make it, and as the contig data
    tables of a run followed it), where the map of the clusters, grouped by
    cluster, put C before B and merged A with C."""
    import metacooler

    import metator.network as mtn

    assembly, pairs, lengths = _tie_community(tmp_path)
    contig_data = mtn.create_contig_data(str(assembly))
    bins = pd.DataFrame({"chrom": list(lengths), "start": 0, "end": list(lengths.values())})
    contigs = _network_map(tmp_path / "contigs.metacool", pairs, contig_data, bins)
    final = {"G1": "X", "A": "X", "C": "X", "G2": "Y", "B": "Y"}
    clusters = _network_map(tmp_path / "clusters.metacool", pairs, contig_data, bins, final=final)
    assert list(metacooler.MetaCooler(str(contigs)).contigs["name"]) == ["G1", "G2", "A", "B", "C"]
    assert list(metacooler.MetaCooler(str(clusters)).contigs["name"]) == ["G1", "A", "C", "G2", "B"]
    expected = {"MetaTOR_MGE_00001": ("A,B", 1.0), "MetaTOR_MGE_00002": ("C", None)}
    for metacool in (contigs, clusters):
        assert _clusters(mtm.mge_binning(str(metacool), ["C", "B", "A"])[0]) == expected


def test_the_order_of_the_map_changes_no_name(metagenome, metagenome_files, clusters_map, tmp_path):
    """The clusters of a single MGE contig are numbered the longest contig first,
    whatever the order of the map: the map of the clusters of the mock
    community, grouped by cluster, and the map of its contigs, the longest
    first, give the MGE clusters the same names (here with two genome contigs
    in the list of MGE contigs, each a cluster of its own)."""
    import metacooler

    import metator.network as mtn

    metacool, _ = clusters_map
    contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
    contigs = _network_map(tmp_path / "contigs.metacool", metagenome_files["pairs"], contig_data, _bins(metagenome))
    in_contigs = {name: i for i, name in enumerate(metacooler.MetaCooler(str(contigs)).contigs["name"])}
    in_clusters = {name: i for i, name in enumerate(metacooler.MetaCooler(str(metacool)).contigs["name"])}
    # Two contigs of a genome that hosts no MGE, in another order in each map.
    truth = metagenome.truth
    hosts = {mge.host for mge in metagenome.mges}
    genomes = list(truth.Name[(truth.Kind == "genome") & ~truth.Genome.isin(hosts)])
    pair = next((x, y) for x in genomes for y in genomes if in_contigs[x] < in_contigs[y] and in_clusters[x] > in_clusters[y])
    mges = [*_mge_contigs(metagenome), *pair]
    on_contigs = _clusters(mtm.mge_binning(str(contigs), mges)[0])
    assert on_contigs == _clusters(mtm.mge_binning(str(metacool), mges)[0])
    assert len(on_contigs) == 4 and sum(len(contigs.split(",")) == 1 for contigs, _ in on_contigs.values()) == 3


def test_mge_on_a_map_without_clusters(metagenome, metagenome_files, metagenome_inputs, tmp_path, run_cli):
    """A map that holds the network of metator network and no cluster (the map of
    the contigs of a run, before partition) gets the round mge alone: mge reads
    no cluster of the map."""
    import metacooler

    import metator.network as mtn

    contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
    metacool = _network_map(tmp_path / "contigs.metacool", metagenome_files["pairs"], contig_data, _bins(metagenome))
    mges = tmp_path / "mges.txt"
    mges.write_text("".join(f"{name}\n" for name in _mge_contigs(metagenome)))
    summary = mtm.mge_binning(str(metacool), _mge_contigs(metagenome))[0]
    before = _outside_clusters(metacool)
    run_cli("mge", metacool, "-a", metagenome_inputs["assembly"], "-m", mges, "-o", tmp_path / "out")
    metacool = tmp_path / "out" / "contigs.metacool"
    mc = metacooler.MetaCooler(str(metacool))
    assert [entry["name"] for entry in mc.cluster_rounds] == ["mge"]
    _check_mge_round(mc, summary)
    # genome is metacooler's, computed from the genomes of the map.
    assert list(mc.clusters.columns) == ["round", "name", "parent_round", "parent", "n_contigs", "length", "genome",
                                         "GC", "Hit", "BinningScore"]  # fmt: skip
    assert _outside_clusters(metacool, "mge") == before


def test_mge_binning_random(metagenome, clusters_map):
    """--random shuffles the contigs between MGE clusters of the same sizes."""
    metacool, _ = clusters_map
    np.random.seed(1)
    summary, _ = mtm.mge_binning(str(metacool), _mge_contigs(metagenome), random=True)
    assert sorted(summary.ContigsNumber) == sorted(len(mge.contigs) for mge in metagenome.mges)
    assert sorted(",".join(summary.Contigs).split(",")) == sorted(_mge_contigs(metagenome))
    assert summary.BinningScore.isna().all()


def test_mge_contigs_that_the_map_lacks(metagenome, clusters_map, logged):
    """An MGE contig that the map lacks (a map made with a contig filter) is left
    out of the MGE clusters, with a warning that names it; none of the MGE
    contigs in the map is refused: there is nothing to cluster."""
    import metator.metacool as mtmc

    metacool, _ = clusters_map
    summary, _ = mtm.mge_binning(str(metacool), [*_mge_contigs(metagenome), "absent_contig"])
    assert "absent_contig" not in ",".join(summary.Contigs) and len(summary) == 2
    assert (
        "WARNING",
        f"1 MGE contig(s) are not in {metacool}: they are left out of the MGE clusters: 'absent_contig'.",
    ) in logged()
    message = f"None of the MGE contigs is in {metacool}: there is no MGE contig to cluster."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        mtm.mge_binning(str(metacool), ["absent_contig"])


def test_write_mges_fasta(metagenome, metagenome_files, clusters_map, tmp_path):
    """Two FASTA files per MGE cluster, named after it, side by side: its contigs,
    as records of their own (MGE.fa), and one record (MGE.stitched.fa), its
    contigs joined by 200 N, in lines of 80 bases; both in the order given
    (that of the cluster). A file of their names is replaced, another file is
    left as it is, and nothing else is left in the directory."""
    metacool, _ = clusters_map
    sequences = synthetic.read_fasta(metagenome_files["assembly"])
    phage, plasmid = ([c.name for c in mge.contigs] for mge in metagenome.mges)
    clusters = {1: {"Contigs": phage[::-1]}, 3: {"Contigs": plasmid}}
    out = tmp_path / "mges"
    out.mkdir()
    (out / "MetaTOR_MGE_00001.stitched.fa").write_text("an older file")
    (out / "notes.txt").write_text("the user's")
    mtm.write_mges_fasta(str(metacool), str(metagenome_files["assembly"]), clusters, str(out))
    names = [("MetaTOR_MGE_00001", phage[::-1]), ("MetaTOR_MGE_00003", plasmid)]
    for name, contigs in names:
        record = ("N" * 200).join(sequences[contig] for contig in contigs)
        stitched = f">{name}\n" + "".join(f"{record[i : i + 80]}\n" for i in range(0, len(record), 80))
        assert (out / f"{name}.stitched.fa").read_text() == stitched
        records = synthetic.read_fasta(out / f"{name}.fa")
        assert list(records) == contigs and all(records[contig] == sequences[contig] for contig in contigs)
    assert sorted(path.name for path in out.iterdir()) == sorted(
        ["notes.txt", *(f"{name}{suffix}" for name, _ in names for suffix in (".fa", ".stitched.fa"))]
    )
    assert (out / "notes.txt").read_text() == "the user's"


# The MGE clusters stored in the map of the run: the round mge of its
# clustering rounds.


def _round(mc, name="mge"):
    """The clusters of the round name of the map mc."""
    return mc.clusters[(mc.clusters["round"] == name).to_numpy()].reset_index(drop=True)


def _check_mge_round(mc, summary, name="mge"):
    """The round name of the map mc is the summary of the MGE clusters of the run:
    one cluster per row, in its order, named BinName, without a parent, whose
    length and number of contigs are BinLength and ContigsNumber, whose contigs
    are Contigs (in the order of the map, as metacooler stores the members of a
    cluster), with its other columns: GC, Hit and BinningScore (no host column:
    mge finds no host)."""
    stored = _round(mc, name)
    assert list(stored["name"]) == list(summary["BinName"]) and stored["parent"].isna().all()
    assert list(stored["length"]) == list(summary["BinLength"]) and list(stored["n_contigs"]) == list(summary["ContigsNumber"])
    for column in ["GC", "Hit", "BinningScore"]:
        np.testing.assert_array_equal(stored[column].to_numpy(float), summary[column].to_numpy(float))
    assert not {"MetagenomicBin", "AssociationScore", "AssociatedBins"} & set(stored.dropna(axis=1, how="all").columns)
    members = mc.cluster_members(name)
    found = {cluster: list(contigs) for cluster, contigs in members.groupby("cluster", sort=False)["contig"]}
    position = {contig: i for i, contig in enumerate(mc.contigs["name"])}
    assert found == {
        cluster: sorted(contigs.split(","), key=position.__getitem__)
        for cluster, contigs in zip(summary["BinName"], summary["Contigs"], strict=True)
    }


def _expected_summary(metacool, args, random=False):
    """The summary of the MGE clusters of an mge run on the map metacool (args:
    see mge_run), computed before the run."""
    mges = mio.import_mges_contigs(args[args.index("-m") + 1])
    if random:
        np.random.seed(7)
    return mtm.mge_binning(str(metacool), mges, random=random)[0]


def _binning_rounds(mc):
    """The rounds of the binning run that the map mc holds: their entries, their
    clusters and their contigs."""
    members = mc.cluster_members()
    binning = members["round"].isin(["overlapping", "final"]).to_numpy()
    return (
        [entry for entry in mc.cluster_rounds if entry["name"] in ("overlapping", "final")],
        mc.clusters[mc.clusters["round"].isin(["overlapping", "final"]).to_numpy()].reset_index(drop=True),
        members[binning].reset_index(drop=True),
    )


def _check_binning_rounds(held, mc, without=()):
    """The rounds of the binning run that the map mc holds are those that it held
    (held, see _binning_rounds): the same entries and contigs, and the same
    values in the columns of their clusters, but those of without. The
    contigs of a cluster come in the order of the map: with without, which
    a reassignment changes, in any order."""
    entries, clusters, members = _binning_rounds(mc)
    assert entries == held[0]
    expected = held[1].drop(columns=list(without))
    assert expected["Markers"].dtype == "int64" and expected["Circular"].dtype == "bool"
    for column in expected.columns:
        if pd.api.types.is_bool_dtype(expected[column]):
            expected[column] = expected[column].map({True: "True", False: "False"}).astype("str")
        elif pd.api.types.is_integer_dtype(expected[column]) and column not in ("n_contigs", "length"):
            expected[column] = expected[column].astype("float64")
    pd.testing.assert_frame_equal(clusters[list(expected.columns)], expected)
    if without:
        members, expected = (table.sort_values(list(table.columns), ignore_index=True) for table in (members, held[2]))
        pd.testing.assert_frame_equal(members, expected)
    else:
        pd.testing.assert_frame_equal(members, held[2])


def _outside_clusters(path, *entries):
    """Every group, dataset and attribute of a map, as stored, but its clustering
    rounds, and the entries of the metator key of its metadata (that of an mge
    run, by the name of its round)."""
    found = {name: value for name, value in walk(path).items() if not name.startswith("/clusters")}
    kind, keys, attributes = found["/"]
    kept = []
    for key, dtype, value in attributes:
        if key == "metadata":
            metadata = json.loads(value)
            metadata["metator"] = {name: entry for name, entry in metadata.get("metator", {}).items() if name not in entries}
            value = json.dumps(metadata, sort_keys=True)
        kept.append((key, dtype, value))
    found["/"] = (kind, [key for key in keys if key != "clusters"], kept)
    return found


@pytest.mark.parametrize("output", ["--inplace", "the output directory", "-F the output directory", "an entry"])
def test_mge_stores_its_clusters_in_the_map(output, mge_run, tmp_path, monkeypatch, run_cli):
    """mge MAP stores the MGE clusters of the run in the map of its output
    directory (MAP itself, with --inplace; for an entry of a collection, a map
    of its own, named after the entry): the round mge of its clustering
    rounds, with the contigs of each cluster, after the rounds of the binning
    run, which are kept as they were; the round records the version, the
    command line and the parameters of the run, and so does the entry mge of
    the metator metadata of the map, with the columns of the map that the run
    read."""
    import metacooler

    import metator.metacool as mtmc
    from metator.version import __version__

    metacool, args = mge_run
    summary = _expected_summary(metacool, args)
    held = _binning_rounds(metacooler.MetaCooler(str(metacool)))
    before = _outside_clusters(metacool)
    source = metacool
    if output.startswith("an entry"):
        source = f"{tmp_path / 'maps.mmetacool'}::entries/lib1"
        metacooler.create_mmetacool(tmp_path / "maps.mmetacool", {"lib1": metacool})
    metacool.chmod(0o640)
    unchanged = {path: path.read_bytes() for path in [metacool, tmp_path / "maps.mmetacool"] if path.exists()}
    out = tmp_path / "out"
    output_map = out / ("lib1.metacool" if output == "an entry" else metacool.name)
    options = ["--inplace"] if output == "--inplace" else []
    if output.startswith("-F"):
        out.mkdir()
        output_map.write_text("an older map")
        # The FASTA files of an earlier run, which -F replaces.
        (out / "fasta" / "mges").mkdir(parents=True)
        (out / "fasta" / "mges" / "MetaTOR_MGE_09999.fa").write_text(">an older cluster\nACGT\n")
        options.append("-F")
    locks = []
    write_map = mtmc.write_map

    def locked(*args, **kwargs):
        locks.append(sorted(str(path.relative_to(tmp_path)) for path in [*tmp_path.glob("*.lock"), *out.glob("*.lock")]))
        return write_map(*args, **kwargs)

    monkeypatch.setattr("metator.metacool.write_map", locked)
    run_cli("mge", source, *args, *options)
    written = metacool if output == "--inplace" else output_map
    assert locks == [[f"{written.relative_to(tmp_path)}.lock"]]
    assert not list(tmp_path.glob("*.lock")) and not list(out.glob("*.lock"))
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    _check_mge_round(mc, summary)
    assert len(_round(mc)) == 2
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping", "final", "mge"]
    entry = mc.cluster_rounds[-1]
    assert {key: entry[key] for key in ["name", "tool", "version", "step"]} == {"name": "mge", "tool": "metaTOR",
                                                                                 "version": __version__, "step": "mge"}  # fmt: skip
    assert entry["command"] == f"metator mge {' '.join([str(source), *args, *map(str, options)])}"
    assert entry["parameters"]["map"] == str(source) and entry["parameters"]["threshold_bin"] == 0.8
    run = mc.metadata["metator"]["mge"]
    read = {"contigs": ["GC_content", "Hit", "qc_intra_contig_long"], "metadata": ["qc"]}
    assert run == {key: entry[key] for key in ["version", "command", "parameters"]} | {"read": read} and entry["read"] == read
    _check_binning_rounds(held, mc)
    assert _outside_clusters(written, "mge") == before
    if output == "--inplace":
        assert metacool.stat().st_mode & 0o777 == 0o640
        del unchanged[metacool]
    assert {path: path.read_bytes() for path in unchanged} == unchanged
    assert not [path.name for path in [*tmp_path.iterdir(), *out.iterdir()] if path.name.startswith(".")]
    maps = [] if output == "--inplace" else [output_map.name]
    assert sorted(path.name for path in out.iterdir()) == sorted(["fasta", *maps])
    fasta = out / "fasta" / "mges"
    assert sorted(path.name for path in fasta.iterdir()) == sorted(
        f"{name}{suffix}" for name in summary.BinName for suffix in (".fa", ".stitched.fa")
    )
    for name, contigs in zip(summary.BinName, summary.Contigs, strict=True):
        assert list(synthetic.read_fasta(fasta / f"{name}.stitched.fa")) == [name]
        assert sorted(synthetic.read_fasta(fasta / f"{name}.fa")) == sorted(contigs.split(","))
    assert not (tmp_path / "tmp").exists()
    (log,) = logs(tmp_path, "mge")
    logged = log.read_text()
    ended = rf"INFO :: Contact map with the MGE clusters written after 0h 0m [0-9]+s: {re.escape(str(written))}\n"
    assert re.search(ended, logged)
    assert ("Removing" in logged) == output.startswith("-F")


def test_mge_rerun_replaces_its_round(mge_run, tmp_path, run_cli):
    """A second run of mge on the map replaces its round mge, and only it: the map
    holds the MGE clusters of the second run (here with a threshold that merges
    no contig: each MGE contig is a cluster of its own), and the rounds of the
    binning run and everything else as they were."""
    import metacooler

    metacool, args = mge_run
    held = _binning_rounds(metacooler.MetaCooler(str(metacool)))
    before = _outside_clusters(metacool)
    first = _expected_summary(metacool, args)
    run_cli("mge", metacool, *args, "--inplace")
    _check_mge_round(metacooler.MetaCooler(str(metacool)), first)
    second = mtm.mge_binning(str(metacool), mio.import_mges_contigs(args[3]), threshold_bin=100)[0]
    run_cli("mge", metacool, *args, "-s", "100", "-o", tmp_path / "again", "--inplace")
    mc = metacooler.MetaCooler(str(metacool))
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping", "final", "mge"]
    assert mc.cluster_rounds[-1]["parameters"]["threshold_bin"] == 100
    assert mc.metadata["metator"]["mge"]["parameters"]["threshold_bin"] == 100
    _check_mge_round(mc, second)
    assert len(second) == 3 and len(first) == 2
    _check_binning_rounds(held, mc)
    assert _outside_clusters(metacool, "mge") == before


def test_mge_random_is_a_round_of_its_own(mge_run, tmp_path, monkeypatch, run_cli):
    """--random stores its random clustering as the round mge_random, with its
    entry in the metator metadata, and writes fasta/mges_random: the round mge,
    its entry, and fasta/mges of an earlier run are left as they were. A second
    run replaces fasta/mges_random with -F."""
    import metacooler

    metacool, args = mge_run
    run_cli("mge", metacool, *args, "--inplace")
    real = metacooler.MetaCooler(str(metacool))
    rounds, entry = _round(real), real.metadata["metator"]["mge"]
    fasta = {path.name: path.read_bytes() for path in (tmp_path / "out" / "fasta" / "mges").iterdir()}
    summary = _expected_summary(metacool, args, random=True)
    for force in [[], ["-F"]]:
        monkeypatch.setattr("numpy.random.permutation", lambda ids: np.random.RandomState(7).permutation(ids))
        np.random.seed(7)
        run_cli("mge", metacool, *args, "--random", "--inplace", *force)
    mc = metacooler.MetaCooler(str(metacool))
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping", "final", "mge", "mge_random"]
    pd.testing.assert_frame_equal(_round(mc), rounds)
    assert mc.metadata["metator"]["mge"] == entry and mc.metadata["metator"]["mge_random"]["parameters"]["random"] is True
    assert sorted(len(contigs.split(",")) for contigs in _round(mc, "mge_random")["name"].map(
        dict(zip(summary.BinName, summary.Contigs, strict=True)))) == sorted(summary.ContigsNumber)  # fmt: skip
    assert {path.name: path.read_bytes() for path in (tmp_path / "out" / "fasta" / "mges").iterdir()} == fasta
    stitched = sorted((tmp_path / "out" / "fasta" / "mges_random").glob("*.stitched.fa"))
    assert sorted(path.name[: -len(".stitched.fa")] for path in stitched) == sorted(_round(mc, "mge_random")["name"])


@pytest.mark.parametrize("random", [False, True], ids=["fasta/mges", "fasta/mges_random"])
def test_the_fasta_files_of_an_earlier_run(random, mge_run, tmp_path, run_cli):
    """The FASTA files of the MGE clusters of an earlier run (fasta/mges, or
    fasta/mges_random with --random) are refused without -F, before any work,
    and left as they are; a fasta that is no directory is refused with or
    without -F."""
    metacool, args = mge_run
    out = tmp_path / "out"
    directory = out / "fasta" / ("mges_random" if random else "mges")
    directory.mkdir(parents=True)
    (directory / "MetaTOR_MGE_00001.fa").write_text(">an earlier run's\nACGT\n")
    before = metacool.read_bytes()
    options = ["--random"] if random else []
    with pytest.raises(ValueError):
        run_cli("mge", metacool, *args, *options)
    assert metacool.read_bytes() == before
    assert (directory / "MetaTOR_MGE_00001.fa").read_text() == ">an earlier run's\nACGT\n"
    (log,) = logs(tmp_path, "mge")
    assert f"ERROR :: {directory} already existed." in log.read_text()
    shutil.rmtree(out / "fasta")
    (out / "fasta").write_text("a file\n")
    with pytest.raises(ValueError):
        run_cli("mge", metacool, *args, *options, "-F")
    assert metacool.read_bytes() == before and (out / "fasta").read_text() == "a file\n"


def test_mge_random_with_a_seed(mge_run, run_cli):
    """--seed makes the random clustering of --random the same at each run,
    whatever numpy's random state before it (np.random.permutation draws it),
    and is stored with the other parameters of the run."""
    import metacooler

    metacool, args = mge_run
    found = []
    for state in (1, 2):
        np.random.seed(state)
        run_cli("mge", metacool, *args, "--random", "--seed", "11", "--inplace", "-F")
        mc = metacooler.MetaCooler(str(metacool))
        members = mc.cluster_members("mge_random")
        found.append(sorted(zip(members["cluster"].astype(str), members["contig"].astype(str), strict=True)))
        assert mc.metadata["metator"]["mge_random"]["parameters"]["seed"] == 11
    assert found[0] == found[1] and found[0]


def test_seed_random_seeds_python_and_numpy():
    """seed_random seeds both random number generators that a run draws from:
    Python's (igraph's Louvain and Leiden) and numpy's (mge --random); None
    seeds neither."""
    import random

    import metator.cli as mtcli

    draws = []
    for _ in range(2):
        mtcli.seed_random(5)
        draws.append((random.random(), list(np.random.permutation(10))))
    assert draws[0] == draws[1]
    random.seed(1)
    expected = random.random()
    random.seed(1)
    mtcli.seed_random(None)
    assert random.random() == expected


@pytest.mark.parametrize("name", ["mges_bin_summary.tsv", "mges_cluster_summary.tsv", "mges_binned.fa", "mges_binned_random.fa"])
@pytest.mark.parametrize("force", [False, True], ids=["refused", "-F"])
def test_the_summary_of_an_earlier_version(name, force, mge_run, tmp_path, run_cli):
    """An earlier MGE summary (1.x's mges_bin_summary.tsv, a development build's
    mges_cluster_summary.tsv) or FASTA file of the MGE clusters (mges_binned.fa,
    mges_binned_random.fa, which fasta/mges and fasta/mges_random replace), or a
    user's own file of that name, is refused with or without -F, before any
    work, and left as it is."""
    metacool, args = mge_run
    out = tmp_path / "out"
    out.mkdir()
    (out / name).write_text("an earlier run's\n")
    (out / "mges_other.tsv").write_text("another tool's\n")
    before = metacool.read_bytes()
    with pytest.raises(ValueError):
        run_cli("mge", metacool, *args, *(["-F"] if force else []))
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in out.iterdir()) == sorted([name, "mges_other.tsv"])
    assert (out / name).read_text() == "an earlier run's\n"
    assert not list(tmp_path.glob("*.lock"))
    (log,) = logs(tmp_path, "mge")
    assert (
        f"ERROR :: {out} holds results of an earlier run, which this run does not write and does not replace: "
        f"{name!r}. Move them away, or give the run another output directory with -o."
    ) in log.read_text()


@pytest.mark.parametrize(("case", "reason"), [("not there", "No such file or directory"), ("a directory", "Is a directory")])
def test_mge_refuses_a_list_that_cannot_be_read(case, reason, mge_run, tmp_path, run_cli):
    """A list of MGE contigs that cannot be read is refused before anything is
    written, the reason logged as an error, as pipeline --mges refuses it."""
    metacool, args = mge_run
    mges = tmp_path / "missing.txt" if case == "not there" else tmp_path / "a directory"
    if case == "a directory":
        mges.mkdir()
    args[args.index("-m") + 1] = str(mges)
    before = metacool.read_bytes()
    message = f"Cannot read the list of MGE contigs {mges} ({reason}): give a text file, one contig name per line."
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        run_cli("mge", metacool, *args)
    assert metacool.read_bytes() == before and not (tmp_path / "out" / "fasta").exists()
    (log,) = logs(tmp_path, "mge")
    assert f"ERROR :: {message}\n" in log.read_text()


@pytest.mark.parametrize("assembly", ["bzip2", "zip", "no FASTA"])
def test_mge_refuses_an_assembly_that_metacooler_cannot_read(assembly, mge_run, tmp_path, run_cli, logged):
    """An assembly that metacooler cannot write fasta/mges from (compressed
    with bzip2 or zip, or no FASTA file) is refused before any work, as the
    network, the validation and the pipeline refuse it: the map is left as it
    was, nothing is written, and nothing is removed."""
    import bz2
    import zipfile

    metacool, args = mge_run
    fasta = tmp_path / "assembly.fa"
    if assembly == "bzip2":
        fasta.write_bytes(bz2.compress(b">c1\nACGT\n"))
    elif assembly == "zip":
        with zipfile.ZipFile(fasta, "w") as archive:
            archive.writestr("assembly.fa", ">c1\nACGT\n")
    else:
        fasta.write_text("c1\tACGT\n")
    args[args.index("-a") + 1] = str(fasta)
    out = tmp_path / "out"
    out.mkdir()
    (out / "mges_cluster_summary.tsv").write_text("an earlier run's\n")
    before = metacool.read_bytes()
    with pytest.raises(ValueError):
        run_cli("mge", metacool, *args, "-F")
    assert metacool.read_bytes() == before and not list(tmp_path.glob("*.lock"))
    assert [path.name for path in out.iterdir()] == ["mges_cluster_summary.tsv"]
    message = {
        "bzip2": f"{fasta} is compressed with bzip2: give the assembly as a FASTA file, plain or compressed with gzip.",
        "zip": f"{fasta} is compressed with zip: give the assembly as a FASTA file, plain or compressed with gzip.",
        "no FASTA": f"{fasta} is not a FASTA file: give the assembly that the map was made from.",
    }[assembly]
    assert ("ERROR", message) in logged()


def _edit_metadata(path, change):
    """Change the metadata of the map at path, as a program that takes no lock of
    metacooler would (change: a function of the metadata)."""
    import h5py

    with h5py.File(path, "r+") as handle:
        metadata = json.loads(handle.attrs["metadata"])
        change(metadata)
        handle.attrs["metadata"] = json.dumps(metadata)


def _cannot_cluster(metacool, case):
    """Make the map metacool one that mge cannot cluster the MGE contigs of (case:
    see MGE_MAP_ERRORS)."""
    import metacooler

    if case == "a map without contact statistics":
        _edit_metadata(metacool, lambda metadata: metadata.pop("qc"))
    elif case == "contact statistics of the pixels":
        _edit_metadata(metacool, lambda metadata: metadata["qc"].update(source="pixels"))
    elif case == "contact statistics at 10 kb":
        _edit_metadata(metacool, lambda metadata: metadata["qc"].update(distance_threshold=10_000))
    elif case == "a map without the network":
        metacooler.annotate(str(metacool), drop=["contigs.Hit"], inplace=True)


# How to make a map that mge reads, as its refusals say.
HOW = (
    "Make the map from the pairs files of the run, with the contact statistics of its reads at 1001 bp: hicstuff "
    "pipeline -M metacool, or metacooler create pairs, with metacooler 0.4.0 or later (metacooler create pairs "
    "--qc-threshold 1001 with an earlier default), then run metator network on it (or pipeline)."
)

# What mge refuses of its map, before anything is written, the reason logged
# as an error in the log file of the run: what partition MAP refuses of the
# map that it writes (see tests/test_cli.py), its inputs as the map of its
# output directory included, even with -F; and a map that it cannot cluster
# the MGE contigs of (see metator.mge.check_mge_map): without the contact
# statistics of its pairs at 1001 bp or the network of metator network, or
# without any of the MGE contigs.
MGE_MAP_ERRORS = {
    **{
        f"{what} as the map of the output directory": (
            ["{map}", "-F"],
            f"{{out_map}}, where the run writes its map, is {{{name}}}, an input of the run: give the run another "
            "output directory with -o.",
        )
        for what, name in {"the assembly": "assembly", "the MGE contigs": "mges"}.items()
    },
    "a map that is not there": (["{missing}"], "Cannot read the contact map {missing}: no such file."),
    "a URI of a directory": (["{directory}::/"], "Cannot read the contact map {directory}::/: {directory} is not a file."),
    "no map": (
        ["{pairs}"],
        "{pairs} is not a .metacool contact map: give the map of the run as the first argument of mge.",
    ),
    "a collection": (["{collection}"], "Cannot read the contact map {collection}: "),
    "a map whose name holds ::": (
        ["{colon_map}"],
        "Cannot read the contact map {colon_map}: metacooler reads a path that holds :: as file::group. Rename the map "
        "without ::.",
    ),
    "an entry with --inplace": (
        ["{entry}", "--inplace"],
        "{entry} is an entry of a .mmetacool collection, which a run does not update in place: leave out --inplace, "
        "and the run writes its map, a map of its own, in its output directory.",
    ),
    "a directory as the map of the output directory": (
        ["{map}"],
        "{out_map} is a directory, where the run writes its map: move it away, or give the run another output directory with -o.",
    ),
    "an existing map of the output directory": (
        ["{map}"],
        "{out_map} exists: give -F to replace it with the map of this run, or another output directory with -o.",
    ),
    "a lock held by another run": (
        ["{map}", "--inplace"],
        # metacooler's reason (see tests.locks.refusal).
        "{refusal} Or write this run's map elsewhere: leave out --inplace.",
    ),
    "a map without contact statistics": (
        ["{map}"],
        "{map} holds no contact statistics of the reads of its contigs (metacooler's QC block, with "
        "qc_intra_contig_long in its contigs table): metator mge counts the contacts within each MGE contig between "
        f"reads more than 1 kb apart, which the pixels of a map do not tell. {HOW}",
    ),
    "contact statistics of the pixels": (
        ["{map}"],
        "{map} has the contact statistics of its contigs measured from its pixels (QC source 'pixels'), not from the "
        "pairs: metator mge counts the contacts within each MGE contig between reads more than 1 kb apart, which the "
        f"pixels of a map do not tell. {HOW}",
    ),
    "contact statistics at 10 kb": (
        ["{map}"],
        "{map} splits the contacts within a contig at 10000 bp (its QC distance_threshold): metator mge counts those "
        f"between reads more than 1 kb apart, at 1001 bp and more. {HOW}",
    ),
    "a map without the network": (
        ["{map}"],
        "{map} lacks the Hit column(s) of its contigs table that metator network stores: run metator network on the "
        "map first (or pipeline), or give the map that such a run wrote in its output directory (network, pipeline "
        "or validation).",
    ),
    "none of the MGE contigs in the map": (
        ["{map}"],
        "None of the MGE contigs is in {map}: there is no MGE contig to cluster.",
    ),
}


@pytest.mark.parametrize("case", list(MGE_MAP_ERRORS))
def test_mge_refuses_a_map(case, mge_run, metagenome_files, tmp_path, run_cli):
    """mge refuses, before anything is written (no MGE fasta, no map): what
    partition MAP refuses of the map that it writes (a map that is not there,
    no map, one that metacooler cannot open, an entry with --inplace, a map of
    the output directory that is a directory, an input of the run, even with
    -F, or that exists, without -F, a map that another run holds); and a map
    that it cannot cluster the MGE contigs of: without the contact statistics
    of its pairs (or of its pixels, or at another distance than 1001 bp),
    without the network of metator network, or without any of the MGE
    contigs."""
    from contextlib import ExitStack

    import metacooler

    from tests.locks import held_lock, lock_file, refusal

    import metator.metacool as mtmc

    metacool, args = mge_run
    names = {
        "map": metacool,
        "missing": tmp_path / "missing.metacool",
        "directory": tmp_path / "maps",
        "pairs": metagenome_files["pairs"],
        "collection": tmp_path / "maps.mmetacool",
        "colon_map": tmp_path / "run::1.metacool",
        # The map that the run writes in its output directory.
        "out_map": tmp_path / "out" / metacool.name,
        "assembly": args[1],
        "mges": args[3],
    }
    names["directory"].mkdir()
    metacooler.create_mmetacool(names["collection"], {"lib1": metacool})
    names["entry"] = f"{names['collection']}::entries/lib1"
    shutil.copyfile(metacool, names["colon_map"])
    # What the case puts in the output directory, as the map of the run.
    out_map = names["out_map"]
    if case.startswith(("the assembly", "the MGE contigs")):
        out_map.parent.mkdir()
        out_map.symlink_to(Path(names["assembly" if case.startswith("the assembly") else "mges"]).resolve())
    elif case == "a directory as the map of the output directory":
        out_map.mkdir(parents=True)
    elif case == "an existing map of the output directory":
        out_map.parent.mkdir()
        out_map.write_text("an older map")
    _cannot_cluster(metacool, case)
    if case == "none of the MGE contigs in the map":
        mges = tmp_path / "absent.txt"
        mges.write_text("absent_contig\n")
        args = [str(mges) if arg == args[3] else arg for arg in args]
    inputs = {name: Path(names[name]).read_bytes() for name in ["assembly", "mges"]}
    before = metacool.read_bytes()
    options, message = MGE_MAP_ERRORS[case]
    with ExitStack() as held:
        if case == "a lock held by another run":
            # Another run holds the map: another process holds its lock.
            holder = held.enter_context(held_lock(metacool, command="metator mge"))
            names["refusal"] = refusal(metacool, holder)
        message = message.format(**names)
        pattern = f"^{re.escape(message)}" + (".*entries/lib1" if case == "a collection" else "$")
        with pytest.raises(mtmc.MetacoolError, match=pattern):
            run_cli("mge", *[option.format(**names) for option in options], *args)
        assert (case == "a lock held by another run") == lock_file(metacool).exists()
    # Nothing in the output directory but what the case put there.
    assert [path.name for path in (tmp_path / "out").glob("*")] == ([out_map.name] if os.path.lexists(out_map) else [])
    if case == "an existing map of the output directory":
        assert out_map.read_text() == "an older map"
    assert metacool.read_bytes() == before
    assert {name: Path(names[name]).read_bytes() for name in ["assembly", "mges"]} == inputs
    assert not lock_file(metacool).exists()
    (log,) = logs(tmp_path, "mge")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("error", [ValueError("not valid"), KeyboardInterrupt()], ids=["ValueError", "Ctrl-C"])
def test_mge_on_a_map_that_cannot_be_written(error, mge_run, tmp_path, monkeypatch, run_cli):
    """A map with the MGE clusters that cannot be written (here, its validation
    fails) stops the run with an error that says why and that the map is left
    as it was, as partition MAP says it; Ctrl-C stops it too."""
    import metator.metacool as mtmc

    metacool, args = mge_run
    (tmp_path / "maps").mkdir()
    metacool = metacool.rename(tmp_path / "maps" / metacool.name)
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise error

    # The full validation of the staged map, before metacooler renames it
    # over the map.
    monkeypatch.setattr("metacooler.schema.validate", validate)
    message = (
        f"{metacool} is not written: not valid. It is still the input map of the run, as it was. The other outputs of "
        "the run are kept."
    )
    if isinstance(error, ValueError):
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("mge", metacool, *args, "--inplace")
        (log,) = logs(tmp_path, "mge")
        assert f"ERROR :: {message}\n" in log.read_text()
    else:
        assert "Aborted" in run_cli("mge", metacool, *args, "--inplace", exit_code=1).output
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in (tmp_path / "maps").iterdir()) == [metacool.name]
    assert list((tmp_path / "out" / "fasta" / "mges").glob("*.stitched.fa"))


def _dropped(where, rounds):
    """The warning about the clustering rounds that a map no longer holds (where:
    the map, or the file written and the map), each named with what made it."""
    return (
        "WARNING",
        f"{where} the clustering round(s) {rounds}, which this run does not make: a run that replaces the clusters of "
        "a map stores its own rounds only. Run what made them again, on the new clusters, to store them.",
    )


def _final_clusters(metagenome):
    """The final clusters of a validation of the mock metagenome (see
    metator.metacool.FinalClusters): its first contig in a cluster, the others
    in none, no summary and no round."""
    import metator.metacool as mtmc

    return mtmc.FinalClusters({str(metagenome.truth.Name.iloc[0]): "metator_00001_00000"}, None, None, {})


@pytest.mark.parametrize("with_mge", [True, False], ids=["with the round mge", "without"])
@pytest.mark.parametrize("writer", ["rounds", "clusters"])
def test_new_clusters_drop_the_round_mge_with_a_warning(
    writer, with_mge, metagenome, mge_run, tmp_path, run_cli, logged, debug_log
):
    """The map of the clusters of a new binning run on the map does not keep the
    round mge, which mge made from the map of the clusters that it replaces:
    once the map is written, a warning says so, and to run mge again on the new
    clusters."""
    import metacooler

    import metator.clusters as mtcl
    import metator.metacool as mtmc

    metacool, args = mge_run
    if with_mge:
        run_cli("mge", metacool, *args, "--inplace")
    if writer == "rounds":
        outfile = metacool
        mtmc.write_map(str(metacool), str(outfile), rounds=mtcl.Rounds().tables(), round_names=[mtcl.OVERLAPPING], force=True)
        where = f"{metacool} no longer holds"
    else:
        outfile = tmp_path / "new.metacool"
        mtmc.write_map(str(metacool), str(outfile), clusters=_final_clusters(metagenome))
        where = f"{outfile} does not keep, of {metacool},"
    warning = _dropped(where, "mge (by metaTOR mge)")
    assert (warning in logged()) == with_mge
    assert not [message for level, message in logged() if level == "WARNING" and message != warning[1]]
    written = metacooler.MetaCooler(str(outfile))
    assert "mge" not in [entry["name"] for entry in written.cluster_rounds]
    # Nor what made it: the entry mge of the metator metadata of the map.
    assert "mge" not in written.metadata.get("metator", {})
    assert ("mge" in metacooler.MetaCooler(str(metacool)).metadata["metator"]) == (with_mge and writer != "rounds")
    # Once the map is written.
    if with_mge:
        messages = [message for _, message in logged()]
        assert messages.index(warning[1]) > max(i for i, message in enumerate(messages) if " is written: " in message)


def test_new_clusters_name_what_made_each_round_they_drop(mge_run, tmp_path, logged):
    """A round that another tool stored in the map, which the binning did not make
    either, is dropped with a warning that names its tool, or says that its
    entry names none: the warning does not take every such round for the round
    mge."""
    import metator.clusters as mtcl
    import metator.metacool as mtmc

    metacool, _ = mge_run
    contig = Path(mge_run[1][3]).read_text().split()[0]
    for entry in [{"name": "dastool", "tool": "DAS_Tool"}, {"name": "custom"}]:
        rounds = mtcl.Rounds()
        rounds.add(entry["name"], f"{entry['name']}_1", [contig])
        mtmc.write_map(str(metacool), str(metacool), mge=(rounds.tables(), entry), force=True)
    mtmc.write_map(str(metacool), str(metacool), rounds=mtcl.Rounds().tables(), round_names=[mtcl.OVERLAPPING], force=True)
    rounds = "dastool (by DAS_Tool), custom (by a tool it does not name)"
    assert _dropped(f"{metacool} no longer holds", rounds) in logged()


@pytest.mark.parametrize("writer", ["rounds", "clusters"])
def test_a_refused_map_of_new_clusters_drops_no_round(writer, metagenome, mge_run, tmp_path, run_cli, logged):
    """A map of new clusters that is refused (here, a cluster named as metacooler
    marks a contig in no cluster, or an existing file without force) is not
    written: the map keeps its round mge, and no warning says that it is
    dropped."""
    import metator.clusters as mtcl
    import metator.metacool as mtmc

    metacool, args = mge_run
    run_cli("mge", metacool, *args)
    before = metacool.read_bytes()
    outfile = tmp_path / "new.metacool"
    outfile.write_text("an older map")
    with pytest.raises(mtmc.MetacoolError):
        if writer == "rounds":
            mtmc.write_map(str(metacool), str(outfile), rounds=mtcl.Rounds().tables(), round_names=[mtcl.OVERLAPPING])
        else:
            clusters = _final_clusters(metagenome)
            clusters.assignment[str(metagenome.truth.Name.iloc[0])] = "NA"
            mtmc.write_map(str(metacool), str(outfile), clusters=clusters, force=True)
    assert metacool.read_bytes() == before
    assert not [message for level, message in logged() if level == "WARNING" and "clustering round" in message]


def _assignment(mc):
    """The genome of each contig of the map mc, by name."""
    names = mc.genomes["name"].astype(str).to_numpy()
    return dict(zip(mc.contigs["name"].astype(str), names[mc.contigs["genome_id"].to_numpy()], strict=True))


def _values(mc):
    """The values of the genomes of the map mc, indexed by their names: its
    genomes table without the columns that metacooler computes (their
    length, contigs, bins and contacts)."""
    genomes = mc.genomes.set_index(mc.genomes["name"].astype(str).rename(None))
    computed = [column for column in genomes.columns if column in ("name", "length", "n_contigs", "n_bins")]
    return genomes.drop(columns=computed + [column for column in genomes.columns if column.startswith("qc_")])


def _with_mge_genomes(assignment, summary):
    """assignment (the genome of each contig, by name), each contig of an MGE
    cluster of summary in the genome of its cluster."""
    pairs = zip(summary.BinName, summary.Contigs, strict=True)
    mges = {contig: name for name, contigs in pairs for contig in contigs.split(",")}
    return {contig: mges.get(contig, genome) for contig, genome in assignment.items()}


def test_mge_reassign_makes_genomes_of_the_mge_clusters(hosted_run, run_cli, logged):
    """--reassign makes each MGE cluster a genome of the map, named after it,
    with its GC content, hits and clustering score: its contigs leave the
    genomes that held them, the final cluster of the host of the phage, which
    keeps its other contigs and its values, and that of the plasmid, which is
    no genome then. The column mge says which genomes are MGE genomes (True;
    False for the others, none for __unbinned__). The clustering rounds, the
    contigs table and the network are kept, and the map is valid; the log says
    which genomes lost contigs."""
    import metacooler

    import metator.metacool as mtmc
    import metator.network as mtn

    metacool, args = hosted_run
    before = metacooler.MetaCooler(str(metacool))
    held, values, contigs = _assignment(before), _values(before), before.contigs.set_index("name")
    rounds = _binning_rounds(before)
    network = mtn.stored_network(str(metacool)).digest
    summary = _expected_summary(metacool, args)
    run_cli("mge", metacool, *args, "--reassign", "--inplace")
    metacooler.validate(str(metacool), full=True)
    mc = metacooler.MetaCooler(str(metacool))
    assert _assignment(mc) == _with_mge_genomes(held, summary)
    after = _values(mc)
    plasmid = held[summary.Contigs[summary.ContigsNumber == 1].iloc[0]]
    kept = [genome for genome in values.index if genome != plasmid]
    assert list(after.index) == sorted([*kept, *summary.BinName], key=lambda name: (name != "__unbinned__", name))
    pd.testing.assert_frame_equal(after.loc[kept, list(values.columns)], values.loc[kept])
    assert list(after.loc[kept, "mge"]) == ["" if genome == "__unbinned__" else "False" for genome in kept]
    assert (after.loc[list(summary.BinName), "mge"] == "True").all()
    for column in ["GC", "Hit", "BinningScore"]:
        # From a genomes table of text, which metacooler reads with pandas: to
        # the bit with pandas 3, to about one part in 10^13 with pandas 2.
        mges = after.loc[list(summary.BinName), column].to_numpy(float)
        np.testing.assert_allclose(mges, summary[column].to_numpy(float), rtol=1e-12, atol=0)
        assert after.loc[kept, column].isna().all()
    missing = after.loc[list(summary.BinName), list(values.columns)].map(lambda value: value == "" or pd.isna(value))
    assert missing.all().all()
    _check_mge_round(mc, summary)
    # The clusters of the rounds are kept, but the genome of each, which
    # metacooler computes: none for the clusters of the host, whose contigs
    # are no genome any more, and the MGE genome of the plasmid for its own.
    _check_binning_rounds(rounds, mc, without=["genome"])
    binning = _binning_rounds(mc)[1].set_index("name")["genome"]
    host = held[summary.Contigs[summary.ContigsNumber == 2].iloc[0].split(",")[0]]
    plasmid_genome = summary.BinName[summary.ContigsNumber == 1].iloc[0]
    expected = {name: None if name == host else plasmid_genome if name == plasmid else name for name in binning.index}
    assert {name: None if pd.isna(genome) else genome for name, genome in binning.items()} == expected
    # Every column of the contigs table, with its values (metacooler puts its
    # QC columns last).
    pd.testing.assert_frame_equal(
        mc.contigs.set_index("name").loc[contigs.index].drop(columns="genome_id"),
        contigs.drop(columns="genome_id"),
        check_like=True,
    )
    assert mtn.stored_network(str(metacool)).digest == network
    assert ("INFO", "The 2 MGE clusters of the run are genomes of the map, of 3 contig(s): 3 of them leave 2 "
            "genome(s), 0 were in none.") in logged()  # fmt: skip
    assert ("INFO", f"1 genome(s) lose contigs to MGE genomes, and keep the values of their validation, which describe "
            f"them with those contigs: {host} (2).") in logged()  # fmt: skip
    assert ("WARNING", f"1 genome(s) of the map are no genome any more, as all their contigs are in MGE genomes: "
            f"'{plasmid}'. Their values leave the genomes table with them.") in logged()  # fmt: skip
    assert mtmc.mge_genomes(mc) == list(summary.BinName)


def test_mge_reruns_replace_or_undo_the_mge_genomes(hosted_run, tmp_path, run_cli, logged):
    """A run of mge replaces the MGE genomes of an earlier run with --reassign,
    with its own (here each MGE contig a cluster of its own), and does not add
    them; --random leaves them as they are; a run without --reassign undoes
    them: their contigs go back to their final clusters, the host of the phage
    with its values, the cluster of the plasmid without the values that it
    lost with its contig, with a warning, and the columns of the MGE genomes
    are gone."""
    import metacooler

    import metator.metacool as mtmc

    metacool, args = hosted_run
    before = metacooler.MetaCooler(str(metacool))
    held, values = _assignment(before), _values(before)
    first = _expected_summary(metacool, args)
    plasmid = held[first.Contigs[first.ContigsNumber == 1].iloc[0]]
    run_cli("mge", metacool, *args, "--reassign", "--inplace")
    second = mtm.mge_binning(str(metacool), mio.import_mges_contigs(args[3]), threshold_bin=100)[0]
    run_cli("mge", metacool, *args, "-s", "100", "--reassign", "--inplace", "-o", tmp_path / "again")
    mc = metacooler.MetaCooler(str(metacool))
    assert len(second) == 3 and _assignment(mc) == _with_mge_genomes(held, second)
    assert sorted(mtmc.mge_genomes(mc)) == sorted(second.BinName)
    kept = [genome for genome in values.index if genome in set(_assignment(mc).values())]
    pd.testing.assert_frame_equal(_values(mc).loc[kept, list(values.columns)], values.loc[kept])
    reassigned, genomes = _assignment(mc), _values(mc)
    np.random.seed(7)
    run_cli("mge", metacool, *args, "--random", "--inplace", "-o", tmp_path / "random")
    mc = metacooler.MetaCooler(str(metacool))
    assert _assignment(mc) == reassigned
    pd.testing.assert_frame_equal(_values(mc), genomes)
    run_cli("mge", metacool, *args, "--inplace", "-o", tmp_path / "undone")
    mc = metacooler.MetaCooler(str(metacool))
    assert _assignment(mc) == held and not mtmc.mge_genomes(mc)
    after = _values(mc)
    assert list(after.columns) == list(values.columns)
    pd.testing.assert_frame_equal(after.drop(index=plasmid), values.drop(index=plasmid))
    assert after.loc[plasmid].map(lambda value: value == "" or pd.isna(value)).all()
    assert ("INFO", "The 3 MGE genome(s) of the map are undone (" + ", ".join(f"'{name}'" for name in sorted(second.BinName))
            + "): of their 3 contig(s), 3 go back to their final cluster, and 0 to none.") in logged()  # fmt: skip
    assert ("WARNING", f"1 final cluster(s) are genomes of the map again, without values: '{plasmid}'. Their values left "
            "the genomes table with their last contigs, which an MGE genome took: run metator validation again to "
            "measure them.") in logged()  # fmt: skip


def test_mge_undoes_the_mge_genomes_made_by_hand(hosted_run, tmp_path, run_cli):
    """Genomes named after the clusters of the round mge of the map, as
    metacooler reassign --round mge makes them, are MGE genomes too: the next
    run of mge undoes them, though the column mge does not say so."""
    import metacooler

    metacool, args = hosted_run
    held = _assignment(metacooler.MetaCooler(str(metacool)))
    run_cli("mge", metacool, *args, "--inplace")
    summary = _expected_summary(metacool, args)
    metacooler.reassign(str(metacool), round="mge", missing="keep", inplace=True)
    assert _assignment(metacooler.MetaCooler(str(metacool))) == _with_mge_genomes(held, summary)
    run_cli("mge", metacool, *args, "--inplace", "-o", tmp_path / "again")
    assert _assignment(metacooler.MetaCooler(str(metacool))) == held


def test_mge_refuses_reassign_with_random(mge_run, tmp_path, run_cli):
    """--reassign with --random is refused before anything is read or written:
    a random clustering makes no genome."""
    metacool, args = mge_run
    before = metacool.read_bytes()
    result = run_cli("mge", metacool, *args, "--random", "--reassign", "--inplace", exit_code=2)
    message = (
        "Invalid value for '--reassign': --reassign makes genomes of the MGE clusters of the round mge, and --random "
        "stores a random clustering as the round mge_random, which leaves them as they were: give one or the other."
    )
    assert message in " ".join(result.output.replace("│", " ").split())
    assert metacool.read_bytes() == before
    assert not (tmp_path / "out").exists() and not list(logs(tmp_path, "mge"))


def test_mge_reassign_refuses_an_mge_cluster_named_after_a_genome(mge_run, run_cli):
    """An MGE cluster named after a genome that the run keeps (here a genome
    renamed by hand) is refused, and the map is left as it was: their contigs
    would make one genome."""
    import metacooler

    import metator.metacool as mtmc

    metacool, args = mge_run
    summary = _expected_summary(metacool, args)
    held = _assignment(metacooler.MetaCooler(str(metacool)))
    taken = next(genome for genome in held.values() if genome != "__unbinned__")
    renamed = {contig: summary.BinName[0] if genome == taken else genome for contig, genome in held.items()}
    metacooler.reassign(
        str(metacool), {contig: genome for contig, genome in renamed.items() if genome != "__unbinned__"}, inplace=True
    )
    before = metacool.read_bytes()
    named = rf"The MGE cluster\(s\) '{summary.BinName[0]}' of the run are named after genome\(s\) of the map that the run keeps"
    with pytest.raises(mtmc.MetacoolError, match=named):
        run_cli("mge", metacool, *args, "--reassign", "--inplace")
    assert metacool.read_bytes() == before


def test_mge_undoes_the_mge_genomes_of_a_map_without_values(mge_run, tmp_path, run_cli):
    """On a map whose genomes have no values (its final clusters, the MGE contigs
    in none), a run without --reassign undoes the MGE genomes of an earlier run
    with it: the genomes are the final clusters again, with the columns that
    they had, none of those of the MGE genomes."""
    import metacooler

    metacool, args = mge_run
    before = metacooler.MetaCooler(str(metacool))
    held, columns = _assignment(before), list(before.genomes.columns)
    run_cli("mge", metacool, *args, "--reassign", "--inplace")
    assert _assignment(metacooler.MetaCooler(str(metacool))) != held
    run_cli("mge", metacool, *args, "--inplace", "-o", tmp_path / "undone")
    mc = metacooler.MetaCooler(str(metacool))
    assert _assignment(mc) == held and list(mc.genomes.columns) == columns
