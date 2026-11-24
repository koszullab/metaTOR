"""metator dump: the tables of a run, written from the map of its clusters."""

import json
import os
import random
import shutil
import sys
from pathlib import Path

import metacooler
import numpy as np
import pandas as pd
import pytest

import metator.dump as mtd
import metator.io as mio
import metator.metacool as mtmc
import metator.network as mtn
from tests import synthetic
from tests.maps import contig_map

TABLES = list(mtd.TABLES)

#: The contigs of the run with no contact with another contig.
ISOLATED = ["contig_004", "contig_011", "contig_018"]


def _without_contacts(pairs, isolated, path):
    """The pairs file at pairs without the pairs of the isolated contigs with
    another contig, written to path."""
    with open(pairs) as source, open(path, "w") as kept:
        for line in source:
            fields = line.split("\t")
            if not line.startswith("#") and fields[1] != fields[3] and (fields[1] in isolated or fields[3] in isolated):
                continue
            kept.write(line)
    return path


@pytest.fixture(scope="module")
def run(metagenome, metagenome_files, tmp_path_factory):
    """A seeded pipeline run of the mock community, its map of the clusters (map),
    the tables it kept in its temporary directory (kept), and its assembly."""
    from click.testing import CliRunner

    from metator.cli import cli

    directory = tmp_path_factory.mktemp("dump_run")
    for key in ("assembly", "depth"):
        shutil.copyfile(metagenome_files[key], directory / metagenome_files[key].name)
    assembly, depth = (directory / metagenome_files[key].name for key in ("assembly", "depth"))
    pairs = _without_contacts(metagenome_files["pairs"], ISOLATED, directory / "isolated.pairs")
    contig_data = mtn.create_contig_data(str(assembly), str(depth), "DpnII")
    metacool = contig_map(pairs, contig_data, directory / "contigs.metacool", assembly_name="synthetic")
    args = ["pipeline", "contigs.metacool", "-a", assembly.name, "-d", depth.name, "-e", "DpnII", "-s", "20000",
            "-i", "10", "-j", "5", "-o", "out", "-T", "tmp", "-N"]  # fmt: skip
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(directory)
        monkeypatch.setattr(sys, "argv", ["metator", *args])
        monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.contaminated_twice(metagenome))
        random.seed(0)
        result = CliRunner().invoke(cli, args, catch_exceptions=False)
    assert result.exit_code == 0, result.output
    (kept,) = [path.parent for path in (directory / "tmp").glob("*/contig_data_final.txt")]
    # The map of the clusters, in the output directory of the run.
    return {"map": directory / "out" / metacool.name, "kept": kept, "assembly": assembly}


@pytest.fixture
def copied(run, tmp_path):
    """A private copy of the map of the run."""
    path = tmp_path / "clusters.metacool"
    shutil.copyfile(run["map"], path)
    return path


def _lines(path):
    return Path(path).read_text().splitlines()


def _table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def test_the_run_takes_every_path_of_the_tables(run):
    """What the other tests rely on: contigs in no core cluster and in no final
    cluster, a step stored as decimals, and two recursive steps."""
    written = _table(run["kept"] / "contig_data_final.txt").set_index("Name")
    assert (written.loc[ISOLATED, ["Core_cluster_ID", "Final_cluster"]] == ["-", "ND"]).all().all()
    assert (written["Recursive_cluster_ID"] != "00000").any()
    mc = metacooler.MetaCooler(str(run["map"]))
    assert "__unbinned__" in set(mc.genomes["name"]) and mc.genomes["step"].dtype == np.float64
    assert [entry["name"] for entry in mc.cluster_rounds] == ["core", "overlapping", "recursive_1", "recursive_2", "final"]
    # A contig of a cluster that step 2 split again: in both recursive rounds.
    members = mc.cluster_members()
    first, second = (set(members.contig[(members["round"] == name).to_numpy()]) for name in ["recursive_1", "recursive_2"])
    assert first & second


def test_dump_writes_the_tables_that_the_run_kept(run, tmp_path, run_cli):
    """The tables of the dump are those that the run wrote, but for the ID of the
    contigs and the order of the rows of the contigs: each in the order of the
    map, which the run regrouped by cluster."""
    out = tmp_path / "dump"
    run_cli("dump", run["map"], "-o", out)
    assert sorted(path.name for path in out.iterdir() if not path.name.endswith(".log")) == sorted(TABLES)
    kept = run["kept"]
    # Byte for byte: the summary of the final clusters (its step a whole
    # number, though the map stores it as decimals) and the clustering rounds,
    # both recursive rounds included.
    for name in ["cluster_summary.txt", "clusters.txt"]:
        assert (out / name).read_bytes() == (kept / name).read_bytes(), name
    rounds = {line.split("\t")[0] for line in _lines(out / "clusters.txt")[1:]}
    assert rounds == {"core", "overlapping", "recursive_1", "recursive_2", "final"}
    # The same contigs, in each of their clusters and in their final one: no
    # contig in no final cluster in binning.txt.
    for name in ["binning.txt", "cluster_members.txt"]:
        assert sorted(_lines(out / name)) == sorted(_lines(kept / name)), name
    assert not {line.split("\t")[0] for line in _lines(out / "binning.txt")} & set(ISOLATED)
    # The contig data: every column of the run's table but its index and ID,
    # each value as the run wrote it (the - of the contigs in no core cluster,
    # the IDs of the others with five digits, the recursive cluster of the last
    # step, ND).
    written = _table(kept / "contig_data_final.txt")
    dumped = _table(out / "contig_data_final.txt")
    assert list(dumped.columns) == [column for column in written.columns if column not in ("", "Unnamed: 0", "ID")]
    pd.testing.assert_frame_equal(dumped.set_index("Name"), written.set_index("Name")[list(dumped.columns[1:])].loc[dumped.Name])
    assert set(dumped.loc[dumped.Name.isin(ISOLATED), "Core_cluster_ID"]) == {"-"}
    assert dumped.loc[~dumped.Name.isin(ISOLATED), "Core_cluster_ID"].str.fullmatch(r"\d{5}").all()
    mc = metacooler.MetaCooler(str(run["map"]))
    assert list(dumped.Name) == list(mc.contigs.name)
    # The network: the stored one, by contig name, without the contacts
    # within a contig, in the order of the map.
    network = pd.read_csv(out / "network.txt", sep="\t", header=None, names=["contig1", "contig2", "weight"])
    stored = mtn.map_network(str(run["map"]))
    names = mc.contigs.name.to_numpy()
    assert list(network.contig1) == list(names[stored.contig1 - 1]) and list(network.contig2) == list(names[stored.contig2 - 1])
    assert network.weight.to_numpy() == pytest.approx(stored.weight.to_numpy(), rel=1e-15)
    assert not (network.contig1 == network.contig2).any()
    assert not (set(network.contig1) | set(network.contig2)) & set(ISOLATED)


def test_dump_gives_each_cluster_its_size(copied, metagenome, run, tmp_path, run_cli):
    """clusters.txt gives every cluster of every round its number of contigs and
    their length, after its parent (n_contigs, length): those that the map
    computes from its contigs, for the MGE clusters that metator mge stores
    there too, which 1.x's mges_bin_summary.tsv gave as ContigsNumber and
    BinLength."""
    truth = metagenome.truth
    mges = tmp_path / "mges.txt"
    mges.write_text("".join(f"{name}\n" for name in truth.Name[truth.Kind == "mge"]))
    run_cli("mge", copied, "-a", run["assembly"], "-m", mges, "-o", tmp_path / "mge", "-T", tmp_path / "tmp", "--inplace")
    run_cli("dump", copied, "-o", tmp_path / "dump")
    dumped = pd.read_csv(
        tmp_path / "dump" / "clusters.txt", sep="\t", dtype={"name": str}, keep_default_na=False, na_values=["NA"]
    )
    assert list(dumped.columns[:6]) == ["round", "name", "parent_round", "parent", "n_contigs", "length"]
    computed = metacooler.MetaCooler(str(copied)).clusters
    assert list(dumped["round"]) == list(computed["round"]) and list(dumped["name"]) == list(computed["name"])
    for column in ["n_contigs", "length"]:
        assert dumped[column].tolist() == computed[column].astype("int64").tolist(), column
    mge = dumped[dumped["round"] == "mge"]
    assert len(mge) and (mge["n_contigs"] > 0).all() and (mge["length"] > 0).all()


def test_dump_of_an_entry_is_that_of_its_map(run, tmp_path, run_cli):
    """An entry of a .mmetacool collection dumps as the map it holds."""
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"run": Path(run["map"])})
    run_cli("dump", run["map"], "-o", tmp_path / "map")
    run_cli("dump", f"{collection}::entries/run", "-o", tmp_path / "entry")
    for name in TABLES:
        assert (tmp_path / "entry" / name).read_bytes() == (tmp_path / "map" / name).read_bytes(), name


def test_dump_refuses_a_map_without_a_validation(metagenome_network, tmp_path, run_cli, logged):
    """A map of the contigs, with its network but no clusters, is refused, naming
    both what it lacks."""
    metacool = tmp_path / "contigs.metacool"
    shutil.copyfile(metagenome_network["map"], metacool)
    with pytest.raises(
        mtmc.MetacoolError, match="lacks the validation entry of its metator metadata and the clustering round final"
    ):
        run_cli("dump", metacool, "-o", tmp_path / "out")
    assert [level for level, _ in logged() if level == "ERROR"] == ["ERROR"]
    assert not [path for path in (tmp_path / "out").iterdir() if not path.name.endswith(".log")]


def test_dump_refuses_a_map_without_the_round_final(copied, tmp_path):
    """The validation entry alone is not a run to dump: the round final, which the
    validation stores with it, is named."""
    mc = metacooler.MetaCooler(str(copied))
    rounds = [entry for entry in mc.cluster_rounds if entry["name"] != "final"]
    clusters = mc.clusters
    clusters = clusters[(clusters["round"] != "final").to_numpy()].drop(columns=mtmc._computed_cluster_columns())
    members = mc.cluster_members()
    members = members[(members["round"] != "final").to_numpy()]
    metacooler.write_clusters(str(copied), rounds, clusters, members)
    with pytest.raises(mtmc.MetacoolError, match=r"it lacks the clustering round final\. Run metator validation"):
        mtd.dump(str(copied), str(tmp_path / "out"))


def _genomes(mc):
    """The genome of each contig of a map, by name."""
    names = mc.genomes["name"].astype(str).to_numpy()
    return dict(zip(mc.contigs["name"].astype(str), names[mc.contigs["genome_id"].to_numpy()], strict=True))


def _summary(path):
    """A summary of the clusters (cluster_summary.txt), as text, indexed by the
    names of its rows."""
    return pd.read_csv(path, sep="\t", index_col=0, dtype=str, keep_default_na=False)


def _binning(path):
    """The cluster of each contig of binning.txt, by name."""
    return dict(line.split("\t") for line in _lines(path))


def test_dump_follows_clusters_edited_after_the_run(copied, run, tmp_path):
    """A map whose contigs were moved between clusters by hand (metacooler
    reassign) is dumped as its genomes are: the moved contig in its new
    cluster in contig_data_final.txt and binning.txt, every other column as
    the run wrote it; the summary is the run's, the values that the genomes
    table carried over."""
    genomes = _genomes(metacooler.MetaCooler(str(copied)))
    binned = sorted(name for name, genome in genomes.items() if genome != "__unbinned__")
    moved = next(name for name in binned if list(genomes.values()).count(genomes[name]) > 1)
    other = next(genome for genome in genomes.values() if genome not in ("__unbinned__", genomes[moved]))
    assignment = {name: genome for name, genome in genomes.items() if genome != "__unbinned__"}
    assignment[moved] = other
    metacooler.reassign(str(copied), assignment, inplace=True)
    mtd.dump(str(copied), str(tmp_path / "out"))
    kept = run["kept"]
    assert (tmp_path / "out" / "cluster_summary.txt").read_bytes() == (kept / "cluster_summary.txt").read_bytes()
    assert _binning(tmp_path / "out" / "binning.txt") == assignment
    dumped = _table(tmp_path / "out" / "contig_data_final.txt").set_index("Name")
    expected = _table(kept / "contig_data_final.txt").set_index("Name")[list(dumped.columns)].loc[dumped.index]
    expected.loc[moved, "Final_cluster"] = other
    pd.testing.assert_frame_equal(dumped, expected)


def test_dump_follows_genomes_that_are_no_final_cluster(copied, run, tmp_path):
    """A cluster made by hand from contigs in no final cluster (the isolated
    contigs, in __unbinned__) is a genome that is no final cluster: the tables
    put its contigs in it, and the summary writes it after the final clusters,
    without values (the genomes table has none for it). A final cluster whose
    contigs all leave it is no genome, nor a row of the summary; a map without
    genomes writes the header of the summary alone, and an empty binning.txt."""
    genomes = _genomes(metacooler.MetaCooler(str(copied)))
    assignment = {name: genome for name, genome in genomes.items() if genome != "__unbinned__"}
    made = {**assignment, **dict.fromkeys(ISOLATED, "my_bin")}
    metacooler.reassign(str(copied), made, inplace=True)
    mtd.dump(str(copied), str(tmp_path / "out"))
    kept = _lines(run["kept"] / "cluster_summary.txt")
    empty = "my_bin" + "\t" * (len(kept[0].split("\t")) - 1)
    assert _lines(tmp_path / "out" / "cluster_summary.txt") == [*kept, empty]
    assert _binning(tmp_path / "out" / "binning.txt") == made
    gone = assignment[sorted(assignment)[0]]
    metacooler.reassign(str(copied), {name: "my_bin" for name, genome in assignment.items() if genome == gone}, inplace=True)
    mtd.dump(str(copied), str(tmp_path / "out"), force=True)
    left = [line for line in kept if not line.startswith(f"{gone}\t")]
    assert _lines(tmp_path / "out" / "cluster_summary.txt") == [*left, empty]
    metacooler.reassign(str(copied), {name: "ND" for name in genomes}, inplace=True)
    mtd.dump(str(copied), str(tmp_path / "out"), force=True)
    assert _lines(tmp_path / "out" / "cluster_summary.txt") == kept[:1]
    assert _lines(tmp_path / "out" / "binning.txt") == []


def test_dump_follows_the_mge_genomes(copied, run, metagenome, tmp_path, run_cli):
    """After metator mge --reassign, the tables follow the genomes of the map:
    each MGE contig in its MGE genome, in contig_data_final.txt (every other
    column as the run wrote it) and in binning.txt, which lists the final
    clusters and the MGE genomes; the summary writes the final clusters, with
    the values of the validation, then the MGE genomes, with theirs (GC, Hit,
    BinningScore), and the column mge; the clustering rounds keep the history:
    the final clusters of the run, and the MGE clusters."""
    truth = metagenome.truth
    mges = tmp_path / "mges.txt"
    mge_contigs = list(truth.Name[truth.Kind == "mge"])
    mges.write_text("".join(f"{name}\n" for name in mge_contigs))
    held = _genomes(metacooler.MetaCooler(str(copied)))
    assert {held[name] for name in mge_contigs} - {"__unbinned__"}, "no MGE contig is in a final cluster"
    run_cli("mge", copied, "-a", run["assembly"], "-m", mges, "-o", tmp_path / "mge", "-T", tmp_path / "tmp", "--inplace",
            "--reassign")  # fmt: skip
    run_cli("dump", copied, "-o", tmp_path / "dump")
    mc = metacooler.MetaCooler(str(copied))
    members = mc.cluster_members("mge")
    in_mge = dict(zip(members.contig.astype(str), members.cluster.astype(str), strict=True))
    expected = {name: in_mge.get(name, genome) for name, genome in held.items()}
    binned = {name: genome for name, genome in expected.items() if genome != "__unbinned__"}
    assert _binning(tmp_path / "dump" / "binning.txt") == binned
    kept = run["kept"]
    dumped = _table(tmp_path / "dump" / "contig_data_final.txt").set_index("Name")
    data = _table(kept / "contig_data_final.txt").set_index("Name")[list(dumped.columns)].loc[dumped.index]
    data["Final_cluster"] = ["ND" if expected[name] == "__unbinned__" else expected[name] for name in dumped.index]
    pd.testing.assert_frame_equal(dumped, data)
    summary, written = _summary(tmp_path / "dump" / "cluster_summary.txt"), _summary(kept / "cluster_summary.txt")
    clusters = mc.clusters
    round_mge = clusters[(clusters["round"] == "mge").to_numpy()].set_index("name")
    mags = [name for name in written.index if name in set(expected.values())]
    assert list(summary.index) == [*mags, *sorted(round_mge.index)]
    assert list(summary.columns) == [*written.columns, "GC", "Hit", "BinningScore", "mge"]
    pd.testing.assert_frame_equal(summary.loc[mags, list(written.columns)], written.loc[mags])
    assert (summary.loc[mags, ["GC", "Hit", "BinningScore"]] == "").all().all() and (summary.loc[mags, "mge"] == "False").all()
    genomes = summary.loc[sorted(round_mge.index)]
    assert (genomes[list(written.columns)] == "").all().all() and (genomes["mge"] == "True").all()
    for column, format in {"GC": "%.2f", "Hit": "%d", "BinningScore": "%.4f"}.items():
        values = round_mge.loc[genomes.index, column]
        assert list(genomes[column]) == ["" if pd.isna(value) else format % value for value in values], column
    dumped_members = _lines(tmp_path / "dump" / "cluster_members.txt")
    assert sorted(line for line in dumped_members if not line.startswith("mge\t")) == sorted(_lines(kept / "cluster_members.txt"))
    assert {line.split("\t")[0] for line in _lines(tmp_path / "dump" / "clusters.txt")[1:]} == {
        "core", "overlapping", "recursive_1", "recursive_2", "final", "mge"}  # fmt: skip


def test_dump_refuses_the_tables_of_an_earlier_dump(run, tmp_path, run_cli):
    """An output directory that holds a table of an earlier dump is refused,
    before the map is read; -F replaces the tables, and a symbolic link to a
    table is replaced, not its target."""
    out = tmp_path / "out"
    out.mkdir()
    (out / "binning.txt").write_text("earlier\n")
    with pytest.raises(mtmc.MetacoolError, match=r"holds the table\(s\) of an earlier dump: binning.txt\. Give -F"):
        run_cli("dump", run["map"], "-o", out)
    assert (out / "binning.txt").read_text() == "earlier\n"
    target = tmp_path / "elsewhere.txt"
    target.write_text("kept\n")
    (out / "network.txt").symlink_to(target)
    run_cli("dump", run["map"], "-o", out, "-F")
    assert target.read_text() == "kept\n"
    assert not (out / "network.txt").is_symlink()
    assert (out / "binning.txt").read_text() != "earlier\n"
    assert not [path for path in out.iterdir() if path.name.startswith(".")]


def test_dump_refuses_a_directory_in_the_place_of_a_table(run, tmp_path):
    """A directory where a table goes is refused before anything is read or
    written, with or without -F: the earlier tables stay as they are."""
    out = tmp_path / "out"
    out.mkdir()
    (out / "contig_data_final.txt").write_text("earlier\n")
    (out / "binning.txt").mkdir()
    with pytest.raises(mtmc.MetacoolError, match=r"holds a directory where metator dump writes binning\.txt"):
        mtd.dump(str(run["map"]), str(out), force=True)
    assert (out / "contig_data_final.txt").read_text() == "earlier\n"


def test_a_table_that_cannot_be_written_replaces_none(run, tmp_path, monkeypatch, logged):
    """A table that cannot be written (here clusters.txt) leaves the tables of an
    earlier dump as they were, none of them replaced, and no temporary file: a
    MetacoolError with the reason, which the command logs."""
    out = tmp_path / "out"
    out.mkdir()
    for name in TABLES:
        (out / name).write_text(f"earlier {name}\n")
    write_round_tables = mtd.mtcl.write_round_tables

    def full_disk(*args):
        write_round_tables(*args)
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(mtd.mtcl, "write_round_tables", full_disk)
    with pytest.raises(mtmc.MetacoolError, match=r"cannot be written in .*No space left on device; no table was replaced"):
        mtd.dump(str(run["map"]), str(out), force=True)
    assert sorted(path.name for path in out.iterdir()) == sorted(TABLES)
    assert all((out / name).read_text() == f"earlier {name}\n" for name in TABLES)


def test_dump_refuses_a_map_replaced_as_it_reads_its_contacts(copied, run, tmp_path, monkeypatch):
    """A map replaced while the dump reads it (here between its contigs table and
    its contacts, as a run that regroups it in place renames another map over
    it) is refused as replaced, with re-running the dump as the remedy, not as
    a map whose hits are stale; no table is written."""
    other = tmp_path / "other.metacool"
    shutil.copyfile(run["map"], other)
    names = list(metacooler.MetaCooler(str(other)).contigs["name"].astype(str))
    metacooler.reassign(str(other), {name: f"bin_{i % 3}" for i, name in enumerate(names)}, inplace=True)
    map_contig_table = mtn.map_contig_table

    def then_replaced(metacool, mc=None):
        table = map_contig_table(metacool, mc)
        os.replace(other, copied)
        return table

    monkeypatch.setattr(mtn, "map_contig_table", then_replaced)
    with pytest.raises(mtmc.MapReplacedError, match="replaced or written while metator dump read it") as refused:
        mtd.dump(str(copied), str(tmp_path / "out"))
    assert "run metator dump again once the run that writes the map has ended" in str(refused.value)
    assert "Hit" not in str(refused.value)
    assert not (tmp_path / "out").exists()


def test_dump_refuses_a_map_that_changes_while_it_is_read(copied, tmp_path, monkeypatch):
    """A map written while the dump reads its clusters (by a program that takes no
    lock) is refused, and no table is written: the tables are of one map."""
    contig_data_final = mtd.contig_data_final

    def written_meanwhile(metacool, mc, table):
        mtmc._store_metator_entries(str(copied), {"meanwhile": True})
        return contig_data_final(metacool, mc, table)

    monkeypatch.setattr(mtd, "contig_data_final", written_meanwhile)
    with pytest.raises(mtmc.MapReplacedError, match="run metator dump again"):
        mtd.dump(str(copied), str(tmp_path / "out"))
    assert not (tmp_path / "out").exists()


def test_dump_refuses_a_summary_the_map_lacks(copied, tmp_path):
    """The values of the summary come from the genomes table of the map: a column
    that an edit dropped is named, and so is one whose format does not write
    its values (a text written as a number)."""
    metacooler.annotate(str(copied), drop=["genomes.Completeness"], inplace=True)
    with pytest.raises(mtmc.MetacoolError, match=r"lacks the column\(s\) 'Completeness' of the summary"):
        mtd.dump(str(copied), str(tmp_path / "out"))
    validation = _validation(copied)
    validation["summary_columns"].remove("Completeness")
    validation["summary_formats"]["parent"] = "%.2f"
    mtmc._store_metator_entries(str(copied), {"validation": validation})
    with pytest.raises(mtmc.MetacoolError, match=r"The column parent of the genomes table of .* holds a value that"):
        mtd.dump(str(copied), str(tmp_path / "out"))
    assert not (tmp_path / "out").exists()


def _validation(metacool):
    """A copy of the validation entry of the metator metadata of a map."""
    return json.loads(json.dumps(metacooler.MetaCooler(str(metacool)).metadata["metator"]["validation"]))


def _without_rows_and_formats(validation):
    del validation["summary_rows"], validation["summary_formats"]


def _rows_as_text(validation):
    validation["summary_rows"] = " ".join(validation["summary_rows"])


def _columns_as_text(validation):
    validation["summary_columns"] = " ".join(validation["summary_columns"])


def _unknown_format(validation):
    validation["summary_formats"]["step"] = "%5.1f"


@pytest.mark.parametrize(
    "edit, message",
    [
        (
            _without_rows_and_formats,
            r"lacks summary_rows, summary_formats, which record how the summary of its clusters \(cluster_summary.txt\) "
            "is written: a development build of metaTOR 2.0 validated the map",
        ),
        (_rows_as_text, "summary_rows, which must be the names of the clusters"),
        (_columns_as_text, "summary_columns, which must be the names of the columns"),
        (_unknown_format, "summary_formats, which must give each column one of the formats %.2f, %.4f, %d, %s"),
    ],
)
def test_dump_refuses_a_summary_that_the_validation_did_not_record(copied, tmp_path, edit, message):
    """The map records how the validation wrote the summary of its clusters: a map
    that a development build validated, which does not, is refused, and so is
    one that does not record it as metator validation does."""
    validation = _validation(copied)
    edit(validation)
    mtmc._store_metator_entries(str(copied), {"validation": validation})
    with pytest.raises(mtmc.MetacoolError, match=message):
        mtd.dump(str(copied), str(tmp_path / "out"))
    assert not (tmp_path / "out").exists()


def test_cluster_summary_is_the_summary_of_the_run(run):
    """The fields of the summary, written, are the run's cluster_summary.txt: its
    rows and columns in its order, each value in the format of its column,
    though the map stores the step as decimals next to the missing value of
    __unbinned__."""
    mc = metacooler.MetaCooler(str(run["map"]))
    layout = mtd.summary_layout(str(run["map"]), mc.metadata["metator"]["validation"])
    summary = mtd.cluster_summary(str(run["map"]), mc, layout)
    assert summary.to_csv(sep="\t") == (run["kept"] / "cluster_summary.txt").read_text()
    assert list(summary.index) == layout.rows != sorted(layout.rows)


def test_cluster_summary_writes_each_value_in_its_format():
    """Each value in the format of its column, from the genomes table, in the
    order of the rows of the summary: a missing value (NaN, None, an empty
    text) is an empty field, a whole number stored as decimals is written
    whole, a number written with decimals keeps them."""

    class Map:
        genomes = pd.DataFrame(
            {
                "name": ["b", "a", "__unbinned__"],
                "x": [1.005, np.nan, np.nan],
                "n": [3.0, 2.0, np.nan],
                "t": pd.Series(["yes", "", None], dtype=object),
                "ab": [0.12345, 1.0, np.nan],
            }
        )

    layout = mio.SummaryLayout(["a", "b"], ["x", "n", "t", "ab"], {"x": "%.2f", "n": "%d", "t": "%s", "ab": "%.4f"})
    summary = mtd.cluster_summary("map", Map(), layout)
    assert summary.to_csv(sep="\t") == "\tx\tn\tt\tab\na\t\t2\t\t1.0000\nb\t1.00\t3\tyes\t0.1235\n"


def test_dump_refuses_a_map_whose_network_changed(copied, run, tmp_path, run_cli):
    """metator network run on the map after the run that wrote it, with another
    normalization, stores another network there and keeps the clusters: the
    dump refuses the map, naming metator network with the normalization of the
    run, which stores the network of the run again, and gives the tables of the
    run back."""
    common = ["-a", run["assembly"], "-o", tmp_path / "logs", "-T", tmp_path / "tmp", "--inplace"]
    run_cli("network", copied, "-n", "length", *common)
    with pytest.raises(mtmc.MetacoolError) as refused:
        mtd.dump(str(copied), str(tmp_path / "out"))
    assert str(refused.value) == (
        f"The network stored in {copied} is not the one that its clusters were made from: metator network stored "
        "another network in the map (normalization length) after the run that wrote them partitioned its own "
        "(normalization empirical_hit), and kept them. Its tables would give the network of the one and the clusters "
        "of the other. Run metator network -n empirical_hit on the map, which stores the network of the run again if "
        "its contacts are those that the run read, or metator partition, then metator validation, on the map again "
        "(or metator pipeline), which cluster the network stored now."
    )
    assert not (tmp_path / "out").exists()
    run_cli("network", copied, "-n", "empirical_hit", *common)
    mtd.dump(str(copied), str(tmp_path / "out"))
    mtd.dump(str(run["map"]), str(tmp_path / "dump"))
    for name in TABLES:
        assert (tmp_path / "out" / name).read_bytes() == (tmp_path / "dump" / name).read_bytes(), name


def test_dump_refuses_a_map_that_records_no_network_of_its_clusters(copied, tmp_path):
    """A partition entry without the digest of the network that it partitioned
    leaves the network of the clusters unknown: refused."""
    partition = metacooler.MetaCooler(str(copied)).metadata["metator"]["partition"]
    del partition[mtmc.NETWORK_DIGEST]
    mtmc._store_metator_entries(str(copied), {"partition": partition})
    with pytest.raises(mtmc.MetacoolError, match="does not record the network that its clusters were made from"):
        mtd.dump(str(copied), str(tmp_path / "out"))


def test_network_table_of_the_counts(metagenome_network):
    """Without weights (the normalization None), the network is the counts: whole
    numbers, by contig name, no contact within a contig."""
    table = metagenome_network["contig_table"]
    network = mtd.network_table(table, mtn.map_network(str(metagenome_network["map"]), None))
    assert network["weight"].dtype == np.int64 and (network["weight"] > 0).all()
    assert list(network.columns) == ["contig1", "contig2", "weight"]
    assert not (network.contig1 == network.contig2).any()
    positions = {name: i for i, name in enumerate(table["Name"])}
    assert (network.contig1.map(positions) < network.contig2.map(positions)).all()


def test_final_clusters_are_those_the_run_wrote(run):
    """The final cluster of each contig, ND for those in none (the isolated
    contigs, in metacooler's __unbinned__): as the run wrote them in its contig
    data."""
    mc = metacooler.MetaCooler(str(run["map"]))
    written = _table(run["kept"] / "contig_data_final.txt").set_index("Name")["Final_cluster"]
    final = dict(zip(mc.contigs["name"].astype(str), mtd.final_clusters(mc), strict=True))
    assert final == written.to_dict()
    assert {final[name] for name in ISOLATED} == {"ND"}


def test_an_output_directory_that_cannot_be_made(run, tmp_path):
    """An output directory that cannot be made (here under a file) is a
    MetacoolError naming it, once the map is read."""
    (tmp_path / "file").write_text("")
    with pytest.raises(mtmc.MetacoolError, match=r"The output directory .*file/out cannot be made"):
        mtd.dump(str(run["map"]), str(tmp_path / "file" / "out"))
