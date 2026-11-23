# Test for CLI tools of metator, MGE related


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


@pytest.fixture(autouse=True)
def private_subset_assembly(private_link, monkeypatch):
    """Each test reads the subset assembly through a link of its own (see
    private_link): two test sessions sharing tests_data removed each
    other's pyfastx index."""
    monkeypatch.setitem(global_args, "SUBSET_FASTA", private_link(global_args["SUBSET_FASTA"]))


def test_mge(tmp_dir, run_cli):
    """Test the metator mge command with subset data from pipeline output."""
    args = (
        "-c {CONTIGS} -a {SUBSET_FASTA} -b {BINNING} -m {SUBSET_MGES} "
        "-o {OUT_TEST} --no-clean-up --threshold-bin 0.5 {SUBSET_PAIRS}"
    ).format(
        OUT_TEST=Path(tmp_dir, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator mge", *args.split(" "))
    run_cli("mge", *args.split(" "))


def test_mge_random(tmp_path, run_cli):
    """--random shuffles the MGE bins and still writes the bins and summary."""
    args = (
        "-c {CONTIGS} -a {SUBSET_FASTA} -b {BINNING} -m {SUBSET_MGES} -o {OUT_TEST} --random --threshold-bin 0.5 {SUBSET_PAIRS}"
    ).format(OUT_TEST=tmp_path / "out", **global_args)
    run_cli("mge", *args.split(" "))
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


def test_build_matrix_contig_order(metagenome, metagenome_files, metagenome_indexed_pairs, tmp_path):
    names, sizes = _mge_contigs(metagenome)
    names, sizes = names[::-1], sizes[::-1]
    plain = tmp_path / "plain.pairs"
    shutil.copyfile(metagenome_files["pairs"], plain)
    indexed = mtm.build_matrix(names, sizes, [str(metagenome_indexed_pairs)])
    assert np.array_equal(indexed, _expected_matrix(metagenome, names))
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
    out = tmp_path / "out"
    out.mkdir()
    summary = mtm.mge_binning(str(metagenome_inputs["assembly"]), contigs_data, mges_list_id, str(out), [str(pairs)], 0.8, random)
    return summary, out


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_generate_mges_fasta(newline, tmp_path):
    """One entry per MGE bin, its contigs (named by numbers here, read by
    name) joined by 200 N, every line ended by a line feed, whatever the
    line ends of the assembly; a header holding a form feed is one line."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_bytes(newline.join([">1 desc\x0cmore", "AAAA", "AA", ">2", "CCCC", ">3", "GGGG", ""]).encode())
    out = tmp_path / "mges.fa"
    mtm.generate_mges_fasta(str(fasta), {1: {"Contigs": [3, 1]}, 2: {"Contigs": [2]}}, str(out))
    assert out.read_bytes() == (">MetaTOR_MGE_00001\nGGGG\n" + "N" * 200 + "AAAA\nAA\n>MetaTOR_MGE_00002\nCCCC\n").encode()


def test_mge_binning(metagenome, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path, monkeypatch):
    """The two contigs of the phage form one MGE bin, the plasmid another.
    The fasta of the bins is written without the pyfastx command, which is
    not on the PATH here."""
    monkeypatch.setenv("PATH", str(tmp_path))
    summary, out = _mge_binning(metagenome, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path)
    found = {frozenset(contigs.split(",")) for contigs in summary.Contigs}
    assert found == {frozenset(c.name for c in mge.contigs) for mge in metagenome.mges}
    lengths = dict(zip(summary.Contigs.map(lambda c: frozenset(c.split(","))), summary.BinLength, strict=True))
    for mge in metagenome.mges:
        assert lengths[frozenset(c.name for c in mge.contigs)] == mge.length
    assert summary.BinningScore.notna().sum() == 1  # only the phage was associated
    written = pd.read_csv(out / "mges_cluster_summary.tsv", sep="\t")
    assert list(written.BinName) == list(summary.BinName)
    # One fasta entry per MGE bin, its contigs separated by N spacers.
    sequences = synthetic.read_fasta(out / "mges_binned.fa")
    assert sorted(sequences) == sorted(summary.BinName)
    contig_seqs = {c.name: c.seq for c in metagenome.contigs}
    for name, contigs in zip(summary.BinName, summary.Contigs, strict=True):
        pieces = [piece for piece in sequences[name].split("N") if piece]
        assert sorted(pieces) == sorted(contig_seqs[c] for c in contigs.split(","))


def test_mge_binning_random(metagenome, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path):
    """--random shuffles the contigs between MGE bins of the same sizes."""
    np.random.seed(1)
    summary, _ = _mge_binning(metagenome, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path, random=True)
    assert sorted(summary.ContigsNumber) == sorted(len(mge.contigs) for mge in metagenome.mges)
    names, _ = _mge_contigs(metagenome)
    assert sorted(",".join(summary.Contigs).split(",")) == sorted(names)


# The MGE clusters stored in the map of the run (--map): the round mge of its
# clustering rounds.


@pytest.fixture(scope="module")
def clusters_map(metagenome, metagenome_files, tmp_path_factory):
    """The .metacool map of the clusters of a binning run of the mock
    metagenome, from its pairs, one bin per contig: one genome per genome
    of the community, its final clusters, the MGE contigs in none, and the
    clustering rounds of the run, overlapping and final (with columns of
    text, decimals, whole numbers and booleans). With the binning and the
    MGE contigs of an mge run on it."""
    import metacooler

    import metator.clusters as mtcl
    import metator.metacool as mtmc

    directory = tmp_path_factory.mktemp("clusters_map")
    truth = metagenome.truth
    genomes = truth[truth.Kind == "genome"]
    clusters = {genome.name: f"metator_{k:05d}_00000" for k, genome in enumerate(metagenome.genomes, start=1)}
    metacool = directory / "clusters.metacool"
    bins = pd.DataFrame({"chrom": truth.Name, "start": 0, "end": truth.Length})
    final = {name: clusters.get(genome, "ND") for name, genome in zip(truth.Name, truth.Genome, strict=True)}
    metacooler.from_pairs(str(metacool), pairs=str(metagenome_files["pairs"]), bins=bins, contig_genome_map=final, ps=False)
    rounds = mtcl.Rounds()
    for genome, cluster in clusters.items():
        contigs = genomes.Name[genomes.Genome == genome]
        columns = {"Completeness": 0.9, "Lineage": "Bacteria", "Markers": 40, "Circular": True}
        rounds.add(mtcl.OVERLAPPING, cluster, contigs, columns=columns)
    rounds.close_round(mtcl.OVERLAPPING)
    for genome, cluster in clusters.items():
        columns = {"Markers": 41, "Circular": False}
        rounds.add(mtcl.FINAL, cluster, genomes.Name[genomes.Genome == genome], parent=cluster, columns=columns)
    mtmc.store_rounds(str(metacool), rounds.tables())
    binning = directory / "binning.txt"
    genomes[["Name", "Genome"]].to_csv(binning, sep="\t", header=False, index=False)
    mges = directory / "mges.txt"
    mges.write_text("".join(f"{name}\n" for name in truth.Name[truth.Kind == "mge"]))
    return metacool, binning, mges


@pytest.fixture
def mge_run(clusters_map, metagenome_inputs, metagenome_network, metagenome_indexed_pairs, tmp_path):
    """A private copy of clusters_map, and the arguments of an mge run on
    the mock metagenome, its outputs in tmp_path/out."""
    source, binning, mges = clusters_map
    metacool = tmp_path / "clusters.metacool"
    shutil.copyfile(source, metacool)
    args = ["-a", metagenome_inputs["assembly"], "-c", metagenome_network["contig_data"], "-b", binning, "-m", mges,
            "-o", tmp_path / "out", "-T", tmp_path / "tmp", metagenome_indexed_pairs]  # fmt: skip
    return metacool, [str(arg) for arg in args]


def _mge_round(mc):
    """The clusters of the round mge of the map mc."""
    return mc.clusters[(mc.clusters["round"] == "mge").to_numpy()].reset_index(drop=True)


def _check_mge_round(mc, summary_file):
    """The round mge of the map mc is the summary of the run: one cluster per
    row, in its order, named BinName, without a parent, whose length and
    number of contigs are BinLength and ContigsNumber, whose contigs are
    Contigs, and with its other columns, as written (NA, a missing value:
    no number, or an empty text)."""
    summary = pd.read_csv(summary_file, sep="\t", keep_default_na=False, na_values=["NA"],
                          dtype=dict.fromkeys(["BinName", "Contigs", "MetagenomicBin", "AssociatedBins"], str))  # fmt: skip
    stored = _mge_round(mc)
    assert list(stored["name"]) == list(summary["BinName"]) and stored["parent"].isna().all()
    assert list(stored["length"]) == list(summary["BinLength"]) and list(stored["n_contigs"]) == list(summary["ContigsNumber"])
    for column in ["GC", "Hit", "BinningScore", "AssociationScore"]:
        np.testing.assert_array_equal(stored[column].to_numpy(float), summary[column].to_numpy(float))
    assert list(stored["MetagenomicBin"]) == list(summary["MetagenomicBin"]) == ["None"] * len(summary)
    assert list(stored["AssociatedBins"]) == list(summary["AssociatedBins"].fillna("")) == [""] * len(summary)
    members = mc.cluster_members("mge")
    found = {cluster: sorted(contigs) for cluster, contigs in members.groupby("cluster")["contig"]}
    assert found == {
        name: sorted(contigs.split(",")) for name, contigs in zip(summary["BinName"], summary["Contigs"], strict=True)
    }
    return summary


def _binning_rounds(mc):
    """The rounds of the binning run that the map mc holds: their entries,
    their clusters and their contigs."""
    members = mc.cluster_members()
    return (
        [entry for entry in mc.cluster_rounds if entry["name"] != "mge"],
        mc.clusters[(mc.clusters["round"] != "mge").to_numpy()].reset_index(drop=True),
        members[(members["round"] != "mge").to_numpy()].reset_index(drop=True),
    )


def _check_binning_rounds(held, mc):
    """The rounds of the binning run that the map mc holds are those that it
    held (held, see _binning_rounds): the same entries and contigs, and the
    same values in the columns of their clusters. A column of whole numbers
    is stored as decimals then, and one of booleans as text (True, False),
    as the round mge has none of them: a column with a missing value can be
    neither. The length and number of contigs of each cluster, which
    metacooler computes, stay whole numbers."""
    entries, clusters, members = _binning_rounds(mc)
    assert entries == held[0]
    expected = held[1].copy()
    assert expected["Markers"].dtype == "int64" and expected["Circular"].dtype == "bool"
    for column in expected.columns:
        if pd.api.types.is_bool_dtype(expected[column]):
            expected[column] = expected[column].map({True: "True", False: "False"}).astype("str")
        elif pd.api.types.is_integer_dtype(expected[column]) and column not in ("n_contigs", "length"):
            expected[column] = expected[column].astype("float64")
    pd.testing.assert_frame_equal(clusters[list(expected.columns)], expected)
    pd.testing.assert_frame_equal(members, held[2])


def _outside_clusters(path):
    """Every group, dataset and attribute of a map, as stored, but its
    clustering rounds."""
    return {name: value for name, value in walk(path).items() if not name.startswith("/clusters")}


@pytest.mark.parametrize("output", ["in place", "--output-map", "-F --output-map", "an entry, --output-map"])
def test_mge_stores_its_clusters_in_the_map(output, mge_run, tmp_path, monkeypatch, run_cli):
    """mge --map MAP stores the MGE clusters of the run in the map: the round
    mge of its clustering rounds, which is mges_cluster_summary.tsv, with the
    contigs of each cluster, after the rounds of the binning run, which are
    kept as they were; the round records the version, the command line and
    the parameters of the run. Everything else is as it was: the genomes and
    contigs tables (the final clusters), the bins, the contacts, the
    metadata. The map is valid in full. In place, the map keeps its mode;
    with --output-map, the map (or the collection of an entry) is left as
    it was, and the map is written to the file (an existing one, with -F,
    is replaced: removed first). The map is locked while the run writes
    it."""
    import metacooler

    import metator.metacool as mtmc
    from metator.version import __version__

    metacool, args = mge_run
    held = _binning_rounds(metacooler.MetaCooler(str(metacool)))
    before = _outside_clusters(metacool)
    source = metacool
    if output.startswith("an entry"):
        source = f"{tmp_path / 'maps.mmetacool'}::entries/lib1"
        metacooler.create_mmetacool(tmp_path / "maps.mmetacool", {"lib1": metacool})
    metacool.chmod(0o640)
    unchanged = {path: path.read_bytes() for path in [metacool, tmp_path / "maps.mmetacool"] if path.exists()}
    output_map = tmp_path / "copy.metacool"
    options = [] if output == "in place" else ["--output-map", output_map]
    if output.startswith("-F"):
        output_map.write_text("an older map")
        options.append("-F")
    locks = []
    write_round = mtmc.write_round

    def locked(*args, **kwargs):
        locks.append(sorted(path.name for path in tmp_path.glob("*.lock")))
        return write_round(*args, **kwargs)

    monkeypatch.setattr("metator.metacool.write_round", locked)
    run_cli("mge", *args, "--map", source, *options)
    written = metacool if output == "in place" else output_map
    assert locks == [[f"{written.name}.metator.lock"]] and not list(tmp_path.glob("*.lock"))
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    _check_mge_round(mc, tmp_path / "out" / "mges_cluster_summary.tsv")
    assert len(_mge_round(mc)) == 2
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping", "final", "mge"]
    entry = mc.cluster_rounds[-1]
    assert {key: entry[key] for key in ["name", "tool", "version", "step"]} == {"name": "mge", "tool": "metaTOR",
                                                                                 "version": __version__, "step": "mge"}  # fmt: skip
    assert entry["command"] == f"metator mge {' '.join([*args, '--map', str(source), *map(str, options)])}"
    assert entry["parameters"]["map"] == str(source) and entry["parameters"]["threshold_bin"] == 0.8
    _check_binning_rounds(held, mc)
    assert _outside_clusters(written) == before
    if output == "in place":
        assert metacool.stat().st_mode & 0o777 == 0o640
        del unchanged[metacool]
    assert {path: path.read_bytes() for path in unchanged} == unchanged
    assert not [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]
    (log,) = tmp_path.glob("metator_mge_*.log")
    logged = log.read_text()
    assert f"INFO :: {written} is written: the contact map, with the round mge of its clustering rounds." in logged
    assert ("Removing" in logged) == output.startswith("-F")


def test_mge_rerun_replaces_its_round(mge_run, tmp_path, run_cli):
    """A second run of mge on the map replaces its round mge, and only it:
    the map holds the MGE clusters of the second run (here with a threshold
    that merges no contig: each MGE contig is a cluster of its own), and the
    rounds of the binning run and everything else as they were."""
    import metacooler

    metacool, args = mge_run
    held = _binning_rounds(metacooler.MetaCooler(str(metacool)))
    before = _outside_clusters(metacool)
    run_cli("mge", *args, "--map", metacool)
    first = _check_mge_round(metacooler.MetaCooler(str(metacool)), tmp_path / "out" / "mges_cluster_summary.tsv")
    run_cli("mge", *args, "-s", "100", "-o", tmp_path / "again", "--map", metacool)
    mc = metacooler.MetaCooler(str(metacool))
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping", "final", "mge"]
    assert mc.cluster_rounds[-1]["parameters"]["threshold_bin"] == 100
    second = _check_mge_round(mc, tmp_path / "again" / "mges_cluster_summary.tsv")
    assert len(second) == 3 and len(first) == 2
    _check_binning_rounds(held, mc)
    assert _outside_clusters(metacool) == before


def _lock(metacool, pid):
    """A lock of the map metacool, as a run of pid on this host holds it."""
    import json
    import socket

    holder = {"pid": pid, "host": socket.gethostname(), "started": "2026-11-18T10:00:00+0000", "command": "metator mge"}
    lock = Path(f"{metacool}.metator.lock")
    lock.write_text(json.dumps(holder))
    return lock


def _map_of(path, lengths):
    """A .metacool map of the contigs lengths (their lengths, by name), one
    bin per contig, with one pair within the first one."""
    import metacooler

    first = next(iter(lengths))
    pairs = path.with_suffix(".pairs")
    pairs.write_text(f"## pairs format v1.0\nr1\t{first}\t10\t{first}\t90\t+\t-\n")
    bins = pd.DataFrame({"chrom": list(lengths), "start": 0, "end": list(lengths.values())})
    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    pairs.unlink()
    return path


# The inputs of an mge run (see mge_run) that --output-map must not name: the
# run would replace them with its map.
MGE_INPUTS = {
    "the assembly": "assembly",
    "a file of the index of the assembly": "index",
    "the contig data": "contigs",
    "the binning": "binning",
    "the MGE contigs": "mges",
    "a pairs file": "pairs",
    "the index of a pairs file": "pairs_index",
}


# What mge refuses of the map that it stores its MGE clusters in, before any
# work, the reason logged as an error in the log file of the run: a map given
# in place of the pairs files, --random with a map, and what partition MAP
# refuses (see tests/test_cli.py), its inputs as --output-map included, even
# with -F.
MGE_MAP_ERRORS = {
    **{
        f"{what} as --output-map": (
            ["--map", "{map}", "--output-map", f"{{{name}}}", "-F"],
            f"--output-map {{{name}}} is {{{name}}}, an input of the run: write the map elsewhere.",
        )
        for what, name in MGE_INPUTS.items()
    },
    "a map of other contigs": (
        ["--map", "{other_map}"],
        "{other_map} is not a contact map of the contigs of {contigs}: 1 of its contig(s) are not in it: 'other'.",
    ),
    "a map of other lengths": (
        ["--map", "{short_map}"],
        "{short_map} is not a contact map of the contigs of {contigs}: 1 contig(s) have another length in it: "
        "'{contig}' ({short} bp in the map, {length} bp in {contigs}).",
    ),
    "a map in place of the pairs files": (
        [],
        "{map} is a contact map, not a pairs file: metator mge counts the contacts within each MGE contig between reads "
        "more than 1 kb apart, which a map does not tell. Give the pairs files of the run, and the map to store the MGE "
        "clusters in with --map.",
    ),
    "--random": (
        ["--map", "{map}", "-r"],
        "--random makes a random clustering, to compare the MGE clusters with: it is not stored in a map. Drop --map.",
    ),
    "a map that is not there": (["--map", "{missing}"], "Cannot read the contact map {missing}: no such file."),
    "a URI of a directory": (
        ["--map", "{directory}::/"],
        "Cannot read the contact map {directory}::/: {directory} is not a file.",
    ),
    "no map": (
        ["--map", "{pairs}"],
        "{pairs} is not a .metacool contact map: give the map of the run with --map, or leave it out.",
    ),
    "a collection": (["--map", "{collection}"], "Cannot read the contact map {collection}: "),
    "a map whose name holds ::": (
        ["--map", "{colon_map}"],
        "Cannot read the contact map {colon_map}: metacooler reads a path that holds :: as file::group. Rename the map "
        "without ::.",
    ),
    "an entry without --output-map": (
        ["--map", "{entry}"],
        "{entry} is an entry of a .mmetacool collection, which a run does not update in place: write the map of its "
        "clusters to another file with --output-map.",
    ),
    "--output-map without a map": (
        ["--output-map", "{copy}"],
        "--output-map writes the map that mge --map MAP updates: give the map with --map, or drop --output-map.",
    ),
    "a directory as --output-map": (
        ["--map", "{map}", "--output-map", "{directory}"],
        "{directory} is a directory: --output-map takes the path of a .metacool file.",
    ),
    "an existing --output-map": (
        ["--map", "{map}", "--output-map", "{copy}"],
        "{copy} exists: give -F to replace it with the map of this run, or another --output-map.",
    ),
    "a lock held by another run": (
        ["--map", "{map}"],
        "{map} is being written by another metaTOR run (pid {pid} on {host}, since 2026-11-18T10:00:00+0000): wait "
        "for it to end, write this run's map elsewhere with --output-map, or, if no such run is going on, remove "
        "{map}.metator.lock.",
    ),
}


@pytest.mark.parametrize("case", list(MGE_MAP_ERRORS))
def test_mge_on_a_map_refused(case, mge_run, metagenome_indexed_pairs, tmp_path, run_cli):
    """mge refuses, before any work (no MGE binning, nothing written): a map
    given in place of its pairs files, as the contacts within a contig need
    the positions of the reads; --random with --map, a random binning; and
    what partition MAP refuses of the map that it updates (a map that is
    not there, no map, one that metacooler cannot open, a map of other
    contigs or lengths than the contig data, an entry without --output-map,
    an --output-map without a map, that is a directory, an input of the
    run, even with -F, or that exists, without -F, a map that another run
    holds). The reason is logged as an error in the log file of the run,
    then raised as a MetacoolError. The map, the inputs, an existing
    --output-map and the lock of another run are left as they were."""
    import socket

    import metacooler

    import metator.metacool as mtmc

    metacool, args = mge_run
    names = {
        "map": metacool,
        "missing": tmp_path / "missing.metacool",
        "directory": tmp_path / "maps",
        "pairs": metagenome_indexed_pairs,
        "collection": tmp_path / "maps.mmetacool",
        "colon_map": tmp_path / "run::1.metacool",
        "copy": tmp_path / "copy.metacool",
        "pid": os.getpid(),
        "host": socket.gethostname(),
    }
    names["directory"].mkdir()
    metacooler.create_mmetacool(names["collection"], {"lib1": metacool})
    names["entry"] = f"{names['collection']}::entries/lib1"
    shutil.copyfile(metacool, names["colon_map"])
    names["copy"].write_text("an older map")
    if case == "a lock held by another run":
        _lock(metacool, os.getpid())
    names.update(zip(["assembly", "contigs", "binning", "mges"], args[1:8:2], strict=True))
    # The pyfastx index of the assembly, left by an earlier run.
    names["index"] = f"{names['assembly']}.fxi"
    names["pairs_index"] = f"{names['pairs']}.px2"
    Path(names["index"]).write_text("index")
    contig_data = pd.read_csv(names["contigs"], sep="\t")
    names["contig"], names["length"] = str(contig_data.Name[0]), int(contig_data.Size[0])
    names["short"] = names["length"] - 500
    names["other_map"] = _map_of(tmp_path / "other.metacool", {names["contig"]: names["length"], "other": 1000})
    names["short_map"] = _map_of(tmp_path / "short.metacool", {names["contig"]: names["short"]})
    inputs = {name: Path(names[name]).read_bytes() for name in MGE_INPUTS.values()}
    before = metacool.read_bytes()
    options, message = MGE_MAP_ERRORS[case]
    message = message.format(**names)
    if case == "a map in place of the pairs files":
        args = [*args[:-1], str(metacool)]
    pattern = f"^{re.escape(message)}" + (".*entries/lib1" if case == "a collection" else "$")
    with pytest.raises(mtmc.MetacoolError, match=pattern):
        run_cli("mge", *args, *[option.format(**names) for option in options])
    assert not (tmp_path / "out").exists()
    assert metacool.read_bytes() == before and names["copy"].read_text() == "an older map"
    assert {name: Path(names[name]).read_bytes() for name in MGE_INPUTS.values()} == inputs
    assert (case == "a lock held by another run") == Path(f"{metacool}.metator.lock").exists()
    (log,) = tmp_path.glob("metator_mge_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("error", [ValueError("not valid"), KeyboardInterrupt()], ids=["ValueError", "Ctrl-C"])
def test_mge_on_a_map_that_cannot_be_written(error, mge_run, tmp_path, monkeypatch, run_cli):
    """A map with the MGE clusters that cannot be written (here, its
    validation fails) stops the run with an error that says why and that
    the map is left as it was, as partition MAP says it; Ctrl-C stops it
    too. The map is left byte for byte, with no file beside it (no staged
    map, no lock); the other outputs of the run are kept."""
    import metator.metacool as mtmc

    metacool, args = mge_run
    (tmp_path / "maps").mkdir()
    metacool = metacool.rename(tmp_path / "maps" / metacool.name)
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise error

    monkeypatch.setattr("metacooler.validate", validate)
    message = (
        f"{metacool} is not written: not valid. It is still the input map of the run, as it was. The other outputs of "
        "the run are kept."
    )
    if isinstance(error, ValueError):
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("mge", *args, "--map", metacool)
        (log,) = tmp_path.glob("metator_mge_*.log")
        assert f"ERROR :: {message}\n" in log.read_text()
    else:
        assert "Aborted" in run_cli("mge", *args, "--map", metacool, exit_code=1).output
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in (tmp_path / "maps").iterdir()) == [metacool.name]
    assert (tmp_path / "out" / "mges_cluster_summary.tsv").exists()


def test_mge_on_a_map_without_rounds(mge_run, tmp_path, run_cli):
    """A map that holds no clustering round (as the map that metator
    metacool writes) gets the round mge alone, with the columns of the
    summary; the rest of the map is as it was."""
    import metacooler

    metacool, args = mge_run
    metacooler.drop_clusters(str(metacool))
    before = walk(metacool)
    run_cli("mge", *args, "--map", metacool)
    mc = metacooler.MetaCooler(str(metacool))
    assert [entry["name"] for entry in mc.cluster_rounds] == ["mge"]
    _check_mge_round(mc, tmp_path / "out" / "mges_cluster_summary.tsv")
    # genome is metacooler's, computed from the genomes of the map.
    assert list(mc.clusters.columns) == ["round", "name", "parent_round", "parent", "n_contigs", "length", "genome",
                                         "GC", "Hit", "BinningScore", "MetagenomicBin", "AssociationScore",
                                         "AssociatedBins"]  # fmt: skip
    after = walk(metacool)
    kind, keys, attributes = before.pop("/")
    assert after.pop("/") == (kind, sorted([*keys, "clusters"]), attributes)
    assert {name: value for name, value in after.items() if not name.startswith("/clusters")} == before


def _dropped(where, rounds):
    """The warning about the clustering rounds that a map no longer holds
    (where: the map, or the file written and the map), each named with what
    made it."""
    return (
        "WARNING",
        f"{where} the clustering round(s) {rounds}, which this run does not make: a run that replaces the clusters of "
        "a map stores its own rounds only. Run what made them again, on the new clusters, to store them.",
    )


def _final_contig_data(metagenome, path, shorter=0):
    """A contig data file of the mock metagenome with no binned contig, as a
    validation writes it, its contigs shorter by shorter bp."""
    truth = metagenome.truth
    table = {"ID": range(1, len(truth) + 1), "Name": truth.Name, "Size": truth.Length - shorter, "Final_cluster": "ND"}
    pd.DataFrame(table).to_csv(path, sep="\t", index=False)
    return path


@pytest.mark.parametrize("with_mge", [True, False], ids=["with the round mge", "without"])
@pytest.mark.parametrize("writer", ["write_rounds", "reassign_metacool"])
def test_new_clusters_drop_the_round_mge_with_a_warning(writer, with_mge, metagenome, mge_run, tmp_path, run_cli, logged):
    """The map of the clusters of a new binning run on the map does not keep
    the round mge, which mge --map made from the clusters that it replaces:
    once the map is written, a warning says so, and to run mge again on the
    new clusters. So writes write_rounds the map of partition MAP (in
    place, here), and reassign_metacool those of validation MAP, pipeline
    MAP and metator metacool MAP (to another file, here). Without the round
    mge, the rounds of the binning are replaced with no warning: a guard,
    which passes without the warning too."""
    import metacooler

    import metator.clusters as mtcl
    import metator.metacool as mtmc

    metacool, args = mge_run
    if with_mge:
        run_cli("mge", *args, "--map", metacool)
    if writer == "write_rounds":
        outfile = metacool
        mtmc.write_rounds(str(metacool), str(outfile), mtcl.Rounds().tables(), names=[mtcl.OVERLAPPING], force=True)
        where = f"{metacool} no longer holds"
    else:
        outfile = tmp_path / "new.metacool"
        mtmc.reassign_metacool(
            str(metacool), str(_final_contig_data(metagenome, tmp_path / "contig_data_final.txt")), str(outfile)
        )
        where = f"{outfile} does not keep, of {metacool},"
    warning = _dropped(where, "mge (by metaTOR mge)")
    assert (warning in logged()) == with_mge
    assert not [message for level, message in logged() if level == "WARNING" and message != warning[1]]
    assert "mge" not in [entry["name"] for entry in metacooler.MetaCooler(str(outfile)).cluster_rounds]
    # Once the map is written.
    if with_mge:
        messages = [message for _, message in logged()]
        assert messages.index(warning[1]) > max(i for i, message in enumerate(messages) if " is written: " in message)


def test_new_clusters_name_what_made_each_round_they_drop(mge_run, tmp_path, logged):
    """A round that another tool stored in the map, which the binning did
    not make either, is dropped with a warning that names its tool, or says
    that its entry names none: the warning does not take every such round
    for the round mge."""
    import metator.clusters as mtcl
    import metator.metacool as mtmc

    metacool, _ = mge_run
    contig = str(pd.read_csv(mge_run[1][3], sep="\t").Name[0])
    for entry in [{"name": "dastool", "tool": "DAS_Tool"}, {"name": "custom"}]:
        rounds = mtcl.Rounds()
        rounds.add(entry["name"], f"{entry['name']}_1", [contig])
        mtmc.write_round(str(metacool), str(metacool), rounds.tables(), entry, force=True)
    mtmc.write_rounds(str(metacool), str(metacool), mtcl.Rounds().tables(), names=[mtcl.OVERLAPPING], force=True)
    rounds = "dastool (by DAS_Tool), custom (by a tool it does not name)"
    assert _dropped(f"{metacool} no longer holds", rounds) in logged()


@pytest.mark.parametrize("writer", ["write_rounds", "reassign_metacool"])
def test_a_refused_map_of_new_clusters_drops_no_round(writer, metagenome, mge_run, tmp_path, run_cli, logged):
    """A map of new clusters that is refused (here, a contig data of other
    contig lengths, or an existing file without force) is not written: the
    map keeps its round mge, and no warning says that it is dropped.
    reassign_metacool said it before its checks; write_rounds refused the
    file before it, as here: a guard, which passes without the fix too."""
    import metator.clusters as mtcl
    import metator.metacool as mtmc

    metacool, args = mge_run
    run_cli("mge", *args, "--map", metacool)
    before = metacool.read_bytes()
    outfile = tmp_path / "new.metacool"
    outfile.write_text("an older map")
    with pytest.raises(mtmc.MetacoolError):
        if writer == "write_rounds":
            mtmc.write_rounds(str(metacool), str(outfile), mtcl.Rounds().tables(), names=[mtcl.OVERLAPPING])
        else:
            contig_data = _final_contig_data(metagenome, tmp_path / "contig_data_final.txt", shorter=7)
            mtmc.reassign_metacool(str(metacool), str(contig_data), str(outfile), force=True)
    assert metacool.read_bytes() == before
    assert not [message for level, message in logged() if level == "WARNING" and "clustering round" in message]
