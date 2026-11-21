"""metator dump: the tables of a run, written from the map of its clusters.

One pipeline run of the mock community, seeded, with -N (its tables kept in
its temporary directory), is shared by the module: the dump of its map must
give back the tables it wrote. The run takes every path of the tables: three
contigs have no contact with another contig (in no core cluster, ``-``, and
in no final cluster, ``ND``, so the map has an ``__unbinned__`` genome, next
to which its genomes table stores the step of the clusters as decimals), and
a stub miComplete makes the validation go two recursive steps deep (the
rounds recursive_1 and recursive_2: a contig ends in the cluster of the last
step that split its cluster).
"""

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
    """A seeded pipeline run of the mock community, its map of the clusters
    (map) and the tables it kept in its temporary directory (kept)."""
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
    return {"map": metacool, "kept": kept}


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
    """What the other tests rely on: contigs in no core cluster and in no
    final cluster, a step stored as decimals, and two recursive steps."""
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
    """The tables of the dump are those that the run wrote, but for the ID of
    the contigs and the order of the rows of the contigs: each in the order of
    the map, which the run regrouped by cluster."""
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


def test_dump_of_an_entry_is_that_of_its_map(run, tmp_path, run_cli):
    """An entry of a .mmetacool collection dumps as the map it holds."""
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"run": Path(run["map"])})
    run_cli("dump", run["map"], "-o", tmp_path / "map")
    run_cli("dump", f"{collection}::entries/run", "-o", tmp_path / "entry")
    for name in TABLES:
        assert (tmp_path / "entry" / name).read_bytes() == (tmp_path / "map" / name).read_bytes(), name


def test_dump_refuses_a_map_without_a_validation(metagenome_network, tmp_path, run_cli, logged):
    """A map of the contigs, with its network but no clusters, is refused,
    naming both what it lacks."""
    metacool = tmp_path / "contigs.metacool"
    shutil.copyfile(metagenome_network["map"], metacool)
    with pytest.raises(
        mtmc.MetacoolError, match="lacks the validation entry of its metator metadata and the clustering round final"
    ):
        run_cli("dump", metacool, "-o", tmp_path / "out")
    assert [level for level, _ in logged() if level == "ERROR"] == ["ERROR"]
    assert not [path for path in (tmp_path / "out").iterdir() if not path.name.endswith(".log")]


def test_dump_refuses_a_map_without_the_round_final(copied, tmp_path):
    """The validation entry alone is not a run to dump: the round final, which
    the validation stores with it, is named."""
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


def test_dump_refuses_clusters_edited_after_the_run(copied, tmp_path):
    """A map whose contigs were moved between clusters by hand (metacooler
    reassign) is refused before any table is written: its summary and its
    rounds no longer describe its genomes. The clusters are named."""
    genomes = _genomes(metacooler.MetaCooler(str(copied)))
    binned = sorted(name for name, genome in genomes.items() if genome != "__unbinned__")
    moved = binned[0]
    other = next(genome for genome in genomes.values() if genome not in ("__unbinned__", genomes[moved]))
    assignment = {name: genome for name, genome in genomes.items() if genome != "__unbinned__"}
    assignment[moved] = other
    metacooler.reassign(str(copied), assignment, inplace=True)
    with pytest.raises(mtmc.MetacoolError, match=r"were edited after the run that wrote it") as refused:
        mtd.dump(str(copied), str(tmp_path / "out"))
    named = str(refused.value).split("no longer a genome of the map: ")[1].split(". ")[0]
    assert sorted(named.replace("'", "").split(", ")) == sorted([genomes[moved], other])
    assert "Run metator validation on the map again, or write the edited assignment with metacooler dump" in str(refused.value)
    assert not (tmp_path / "out").exists()


def test_dump_refuses_a_genome_that_is_no_final_cluster(copied, tmp_path):
    """A cluster made by hand from contigs in no final cluster (the isolated
    contigs, in __unbinned__) is a genome that is no final cluster, though
    every final cluster is still a genome: it is named. One that took contigs
    from a final cluster names that cluster too; a map whose every contig was
    unbinned by hand is refused for its final clusters."""
    genomes = _genomes(metacooler.MetaCooler(str(copied)))
    assignment = {name: genome for name, genome in genomes.items() if genome != "__unbinned__"}
    metacooler.reassign(str(copied), {**assignment, **dict.fromkeys(ISOLATED, "my_bin")}, inplace=True)
    with pytest.raises(
        mtmc.MetacoolError, match=r"edited after the run that wrote it .*: genome\(s\) that are no final cluster: 'my_bin'\. "
    ):
        mtd.dump(str(copied), str(tmp_path / "out"))
    taken = sorted(assignment)[0]
    metacooler.reassign(str(copied), {**assignment, taken: "my_bin"}, inplace=True)
    with pytest.raises(
        mtmc.MetacoolError,
        match=rf"no longer a genome of the map: '{assignment[taken]}'; and genome\(s\) that are no final cluster: 'my_bin'",
    ):
        mtd.dump(str(copied), str(tmp_path / "out"))
    metacooler.reassign(str(copied), {name: "ND" for name in genomes}, inplace=True)
    with pytest.raises(mtmc.MetacoolError, match=r"final cluster\(s\) whose contigs are no longer a genome of the map"):
        mtd.dump(str(copied), str(tmp_path / "out"))


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
    """A table that cannot be written (here clusters.txt) leaves the tables of
    an earlier dump as they were, none of them replaced, and no temporary file:
    a MetacoolError with the reason, which the command logs."""
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
    """A map replaced while the dump reads it (here between its contigs table
    and its contacts, as a run that regroups it in place renames another map
    over it) is refused as replaced, with re-running the dump as the remedy,
    not as a map whose hits are stale; no table is written."""
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
    """A map written while the dump reads its clusters (by a program that
    takes no lock) is refused, and no table is written: the tables are of one
    map."""
    contig_data_final = mtd.contig_data_final

    def written_meanwhile(metacool, mc, table):
        mtmc._store_metator_entries(str(copied), {"meanwhile": True})
        return contig_data_final(metacool, mc, table)

    monkeypatch.setattr(mtd, "contig_data_final", written_meanwhile)
    with pytest.raises(mtmc.MapReplacedError, match="run metator dump again"):
        mtd.dump(str(copied), str(tmp_path / "out"))
    assert not (tmp_path / "out").exists()


def test_dump_refuses_a_summary_the_map_lacks(copied, tmp_path):
    """The columns of the summary come from the genomes table of the map: one
    that an edit dropped is named, and so is a validation entry without the
    columns of its summary."""
    metacooler.annotate(str(copied), drop=["genomes.Completeness"], inplace=True)
    with pytest.raises(mtmc.MetacoolError, match=r"lacks the column\(s\) 'Completeness' of the summary"):
        mtd.dump(str(copied), str(tmp_path / "out"))
    mc = metacooler.MetaCooler(str(copied))
    validation = dict(mc.metadata["metator"]["validation"])
    del validation["summary_columns"]
    mtmc._store_metator_entries(str(copied), {"validation": validation})
    with pytest.raises(mtmc.MetacoolError, match="records no columns of the summary of its clusters"):
        mtd.dump(str(copied), str(tmp_path / "out"))


def test_cluster_summary_writes_the_values_as_the_validation(run):
    """The rows and the columns of the summary are those of the run's
    cluster_summary.txt, in its order; the step is a whole number, though the
    genomes table stores it as decimals next to the missing value of
    __unbinned__; the abundances have four decimals."""
    mc = metacooler.MetaCooler(str(run["map"]))
    summary = mtd.cluster_summary(str(run["map"]), mc, mc.metadata["metator"]["validation"])
    written = pd.read_csv(run["kept"] / "cluster_summary.txt", sep="\t", index_col=0, dtype=str, keep_default_na=False)
    assert list(summary.index) == list(written.index) and list(summary.columns) == list(written.columns)
    assert str(summary["step"].dtype) == "Int64" and list(summary["step"].astype(str)) == list(written["step"])
    assert list(summary["HiC_abundance"]) == list(written["HiC_abundance"])


def test_network_table_of_the_counts(metagenome_network):
    """Without weights (the normalization None), the network is the counts:
    whole numbers, by contig name, no contact within a contig."""
    table = metagenome_network["contig_table"]
    network = mtd.network_table(table, mtn.map_network(str(metagenome_network["map"]), None))
    assert network["weight"].dtype == np.int64 and (network["weight"] > 0).all()
    assert list(network.columns) == ["contig1", "contig2", "weight"]
    assert not (network.contig1 == network.contig2).any()
    positions = {name: i for i, name in enumerate(table["Name"])}
    assert (network.contig1.map(positions) < network.contig2.map(positions)).all()


def test_final_clusters_are_those_the_run_wrote(run):
    """The final cluster of each contig, ND for those in none (the isolated
    contigs, in metacooler's __unbinned__): as the run wrote them in its
    contig data."""
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
