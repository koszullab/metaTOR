"""metator dump: the tables of a run, written from the map of its clusters.

One pipeline run of the mock community, seeded, with -N (its tables kept in
its temporary directory) and a stub miComplete that finds the first cluster
contaminated (so the run makes a round recursive_1), is shared by the module:
the dump of its map must give back the tables it wrote.
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
    contig_data = mtn.create_contig_data(str(assembly), str(depth), "DpnII")
    metacool = contig_map(metagenome_files["pairs"], contig_data, directory / "contigs.metacool", assembly_name="synthetic")
    args = ["pipeline", "contigs.metacool", "-a", assembly.name, "-d", depth.name, "-e", "DpnII", "-s", "20000",
            "-i", "10", "-j", "5", "-o", "out", "-T", "tmp", "-N"]  # fmt: skip
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(directory)
        monkeypatch.setattr(sys, "argv", ["metator", *args])
        monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.first_cluster_contaminated(metagenome))
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


def test_dump_writes_the_tables_that_the_run_kept(run, tmp_path, run_cli):
    """The tables of the dump are those that the run wrote, but for the ID of
    the contigs and the order of the rows of the contigs: each in the order of
    the map, which the run regrouped by cluster."""
    out = tmp_path / "dump"
    run_cli("dump", run["map"], "-o", out)
    assert sorted(path.name for path in out.iterdir() if not path.name.endswith(".log")) == sorted(TABLES)
    kept = run["kept"]
    # Byte for byte: the summary of the final clusters and the clustering
    # rounds, the round recursive_1 included.
    for name in ["cluster_summary.txt", "clusters.txt"]:
        assert (out / name).read_bytes() == (kept / name).read_bytes(), name
    assert {line.split("\t")[0] for line in _lines(out / "clusters.txt")[1:]} == {"core", "overlapping", "recursive_1", "final"}
    # The same contigs, in each of their clusters and in their final one.
    for name in ["binning.txt", "cluster_members.txt"]:
        assert sorted(_lines(out / name)) == sorted(_lines(kept / name)), name
    # The contig data: every column of the run's table but its index and ID,
    # each value as the run wrote it.
    written = pd.read_csv(kept / "contig_data_final.txt", sep="\t", dtype=str, keep_default_na=False)
    dumped = pd.read_csv(out / "contig_data_final.txt", sep="\t", dtype=str, keep_default_na=False)
    assert list(dumped.columns) == [column for column in written.columns if column not in ("", "Unnamed: 0", "ID")]
    pd.testing.assert_frame_equal(dumped.set_index("Name"), written.set_index("Name")[list(dumped.columns[1:])].loc[dumped.Name])
    mc = metacooler.MetaCooler(str(run["map"]))
    assert list(dumped.Name) == list(mc.contigs.name)
    assert (dumped.Recursive_cluster_ID != "00000").any()
    # The network: the stored one, by contig name, without the contacts
    # within a contig, in the order of the map.
    network = pd.read_csv(out / "network.txt", sep="\t", header=None, names=["contig1", "contig2", "weight"])
    stored = mtn.map_network(str(run["map"]))
    names = mc.contigs.name.to_numpy()
    assert list(network.contig1) == list(names[stored.contig1 - 1]) and list(network.contig2) == list(names[stored.contig2 - 1])
    assert network.weight.to_numpy() == pytest.approx(stored.weight.to_numpy(), rel=1e-15)
    assert not (network.contig1 == network.contig2).any()


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
    clusters = clusters[(clusters["round"] != "final").to_numpy()].drop(columns=list(mtd.COMPUTED_CLUSTER_COLUMNS))
    members = mc.cluster_members()
    members = members[(members["round"] != "final").to_numpy()]
    metacooler.write_clusters(str(copied), rounds, clusters, members)
    with pytest.raises(mtmc.MetacoolError, match=r"it lacks the clustering round final\. Run metator validation"):
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


def test_dump_refuses_a_map_that_changes_while_it_is_read(copied, tmp_path, monkeypatch):
    """A map written while the dump reads it (by a program that takes no lock)
    is refused, and no table is written: the tables are of one map."""
    contig_data_final = mtd.contig_data_final

    def written_meanwhile(metacool, mc, table):
        mtmc._store_metator_entries(str(copied), {"meanwhile": True})
        return contig_data_final(metacool, mc, table)

    monkeypatch.setattr(mtd, "contig_data_final", written_meanwhile)
    with pytest.raises(mtmc.MetacoolError, match="written"):
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
    """The step of a cluster is a whole number even where the genomes table
    stores it as decimals (next to the missing value of __unbinned__), the
    abundances have four decimals, the other decimals two."""
    mc = metacooler.MetaCooler(str(run["map"]))
    validation = mc.metadata["metator"]["validation"]
    summary = mtd.cluster_summary(str(run["map"]), mc, validation)
    assert list(summary.columns) == validation["summary_columns"]
    assert str(summary["step"].dtype) == "Int64"
    assert all(len(value.split(".")[1]) == 4 for value in summary["HiC_abundance"])
    # The order of the round final: that of the summary of the run.
    final = mc.clusters
    assert list(summary.index) == list(final.loc[(final["round"] == "final").to_numpy(), "name"])


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


def test_final_clusters_are_the_genomes_of_the_map(run):
    """A contig's final cluster is its genome in the map, ND for
    __unbinned__."""
    mc = metacooler.MetaCooler(str(run["map"]))
    final = mtd.final_clusters(mc)
    genomes = mc.contigs.genome_id.map(mc.genomes.name).astype(str)
    assert list(final) == [("ND" if name == "__unbinned__" else name) for name in genomes]
    assert os.path.exists(run["kept"] / "binning.txt")
