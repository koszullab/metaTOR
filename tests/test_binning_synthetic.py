# End-to-end binning tests on the mock metagenome of tests/synthetic.py,
# checked against its ground truth: network, partition, validation (with a
# stub miComplete) and the whole pipeline.

import json
import os
import random
import re
import shlex
import shutil
import subprocess
from pathlib import Path

import cooler
import metacooler
import numpy as np
import pandas as pd
import pytest
import yaml

import metator.cli as mtcli
import metator.clusters as mtcl
import metator.dump as mtd
import metator.metacool as mtmc
import metator.network as mtn
import metator.outputs as mto
import metator.partition as mtp
import metator.report as mtr
import metator.validation as mtv
from metator.version import __version__
from tests import alignments, synthetic
from tests.hdf5 import walk
from tests.maps import contig_map, pairs_map

NORMALIZATIONS = ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"]


# The tables of a validation (contig_data_final.txt, cluster_summary.txt,
# binning.txt, clusters.txt, cluster_members.txt) are in the temporary
# directory of the run, which -N keeps: every run of this module keeps it,
# and the tests read them there (see _kept).
TABLES = ["contig_data_final.txt", "cluster_summary.txt", "binning.txt", "clusters.txt", "cluster_members.txt"]

# The FASTA files of the final clusters of a run, in its output directory.
MAGS_WITH_MGES = os.path.join("fasta", "mags_with-mges")


@pytest.fixture(autouse=True)
def _keep_the_tables(monkeypatch):
    """Keep the temporary directory of every run, as -N does."""
    private_tmp_dir = mtcli.private_tmp_dir
    monkeypatch.setattr(mtcli, "private_tmp_dir", lambda tmpdir, clean_up=True: private_tmp_dir(tmpdir, clean_up=False))


def _kept(out):
    """The directory of the tables of the latest validation or pipeline run whose
    output directory is out: the temporary directory of the run (in the -T of
    its command line, which its log records; a relative one from the directory
    of out, where the tests run) that holds them, the last written before the
    end of its log (two runs may share a -T)."""
    logs = sorted(Path(out).glob("metator_*.log"), key=os.path.getmtime)
    for log in reversed(logs):
        line = next(line for line in log.read_text().splitlines() if line.startswith("## cmd: "))
        command = shlex.split(line[len("## cmd: ") :])
        tmp = Path(command[command.index("-T") + 1] if "-T" in command else "./tmp")
        if not tmp.is_absolute():
            tmp = Path(out).resolve().parent / tmp
        ended = os.path.getmtime(log)
        found = sorted(
            (table for table in tmp.glob("*/contig_data_final.txt") if os.path.getmtime(table) <= ended),
            key=os.path.getmtime,
        )
        if found:
            return found[-1].parent
    raise FileNotFoundError(f"No tables of a validation of {out}")


def _summary_from_the_map(path):
    """The summary of the final clusters of the run that wrote the map of the
    clusters at path (cluster_summary.txt), as the map alone writes it again:
    its genomes table holds the values, and its entry validation how they are
    written (its rows, columns and their formats, see
    metator.cli.validation_entry), an empty field for a missing value."""
    mc = metacooler.MetaCooler(str(path))
    entry = mc.metadata["metator"]["validation"]
    genomes = mc.genomes.set_index("name")
    columns, formats = entry["summary_columns"], entry["summary_formats"]

    def field(value, written):
        return "" if value is None or value == "" or (isinstance(value, float) and np.isnan(value)) else written % value

    lines = ["\t".join(["", *columns])]
    for row in entry["summary_rows"]:
        lines.append("\t".join([row, *(field(genomes.loc[row, column], formats[column]) for column in columns)]))
    return "\n".join(lines) + "\n"


def _named_edges(network, contig_data):
    """Network edges with the names of their contigs (by their IDs in the contig
    data)."""
    names = contig_data.set_index("ID")["Name"]
    return network.assign(name1=network.contig1.map(names), name2=network.contig2.map(names))


def _genome_labels(table, metagenome, column):
    """Truth and found labels of the genome contigs (MGE contigs left out)."""
    truth = metagenome.truth
    table = table[table.Name.map(truth.Kind) == "genome"]
    return list(table.Name.map(truth.Genome)), list(table[column].astype(str))


@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_network_edges_follow_genomes(metagenome, metagenome_files, tmp_path, normalization):
    """With every normalisation, same-genome edges outweigh cross-genome ones."""
    contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
    metacool = contig_map(metagenome_files["pairs"], contig_data, tmp_path / "contigs.metacool", columns=False)
    contacts = mtn.contig_contacts(str(metacool))
    names = metacooler.MetaCooler(str(metacool)).contigs["name"]
    table = mtn.contig_data_table(contig_data, names, mtn.contig_hits(contacts, len(names)))
    network = mtn.map_network(str(metacool), mtn.network_weights(table, normalization), contacts)
    truth = metagenome.truth
    edges = _named_edges(network, table)
    genomes = (edges.name1.map(truth.Kind) == "genome") & (edges.name2.map(truth.Kind) == "genome")
    same = edges.name1.map(truth.Genome) == edges.name2.map(truth.Genome)
    within, across = edges.weight[genomes & same], edges.weight[genomes & ~same]
    assert len(within) > 0 and len(across) > 0
    assert within.mean() > 10 * across.mean()
    assert within.sum() > 0.9 * (within.sum() + across.sum())
    # Every contig is a node of the network.
    assert set(edges.name1) | set(edges.name2) == set(truth.Name)


def test_contig_data_matches_the_truth(metagenome, metagenome_network):
    """Sizes, GC contents and shotgun coverages of the contig data, which metator
    network stores in the map: its contigs, the longest first."""
    data = metagenome_network["contig_table"].set_index("Name")
    truth = metagenome.truth.set_index("Name", drop=False).loc[data.index]
    assert sorted(data.index) == sorted(metagenome.truth.Name)
    assert list(data.Size) == sorted(data.Size, reverse=True)
    assert (data.Size == truth.Length).all()
    assert data.GC_content.to_numpy() == pytest.approx(truth.GC.to_numpy(), abs=1e-6)
    assert data.Shotgun_coverage.to_numpy() == pytest.approx(10 * truth.Abundance.to_numpy(), rel=0.25)
    # Every read pair gives two hits.
    assert data.Hit.sum() == 2 * len(metagenome.pairs)


@pytest.mark.parametrize("algorithm", ["louvain", "leiden"])
def test_partition_recovers_the_genomes(metagenome, metagenome_inputs, metagenome_network, tmp_path, algorithm):
    """The overlapping clusters of the partition are the genomes, and those that
    reach 50 kb, which the validation writes for miComplete to check (see
    metator.validation.write_overlapping_fasta), one each."""
    random.seed(0)  # igraph draws from Python's random module.
    contig_table = metagenome_network["contig_table"]
    clusters = mtp.partition(algorithm, metagenome_network["network_table"], 10, 0.8, 1.0, 1, "metator")
    table = mtp.partition_table(contig_table, clusters)
    truth, found = _genome_labels(table, metagenome, "Overlapping_cluster_ID")
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    # Each core cluster in one overlapping cluster: the overlapping clusters
    # merge them.
    core = {contig: number for number, contigs in clusters.core.items() for contig in contigs}
    for contigs in clusters.overlapping.values():
        assert contigs == sorted(contigs, key=lambda contig: (core[contig], contig))
    assert sorted(core) == sorted(contig for contigs in clusters.overlapping.values() for contig in contigs)
    # One fasta per genome, each holding one genome.
    fasta_dir = tmp_path / "overlapping_clusters"
    written = mtv.write_overlapping_fasta(
        str(metagenome_network["map"]), str(metagenome_inputs["assembly"]), clusters, contig_table, 50_000, str(fasta_dir)
    )
    bins = sorted(fasta_dir.glob("*.fa"))
    assert [path.stem for path in bins] == sorted(written)
    assert len(bins) == len(metagenome.genomes)
    for fasta in bins:
        contigs = synthetic.read_fasta(fasta)
        assert synthetic.bin_quality(list(contigs), metagenome)[1] == 1.0
    # Nothing written next to the assembly: no pyfastx index.
    assert not os.path.exists(f"{metagenome_inputs['assembly']}.fxi")


def _map_with_merged_partition(metagenome, metagenome_network, path, merged):
    """A copy, at path, of the map of the mock community in which metator network
    stored its network (see _network_map), with a partition stored as metator
    partition stores it (see metator.metacool.write_rounds): the genomes
    ``merged`` share a cluster while every other genome has its own, MGE
    contigs going with their host, each cluster both a core and an overlapping
    cluster."""
    truth = metagenome.truth
    bin_id = {}
    for genome in metagenome.genomes:
        if genome.name in merged[1:]:
            bin_id[genome.name] = bin_id[merged[0]]
        else:
            bin_id[genome.name] = len(set(bin_id.values())) + 1
    data = metagenome_network["contig_table"]
    ids = data.Name.map(truth.Host).map(bin_id)
    contigs = {number: sorted(data.ID[ids == number]) for number in sorted(set(ids))}
    clusters = mtcl.Partition(
        core=dict(contigs),
        overlapping=dict(contigs),
        names={number: mtcl.overlapping_name("metator", number) for number in contigs},
    )
    metacool = _network_map(metagenome_network, path)
    entry = {
        "version": __version__,
        "normalization": "empirical_hit",
        mtmc.NETWORK_DIGEST: mtn.stored_network(str(metacool)).digest,
    }
    mtmc.write_map(
        str(metacool),
        str(metacool),
        partition=entry,
        rounds=mtcl.partition_rounds(clusters, list(data.Name)),
        round_names=list(mtmc.PARTITION_ROUNDS),
        force=True,
    )
    return metacool


def _network_map(metagenome_network, path):
    """A copy of the map of the mock community in which metator network stored its
    network, at path: validation MAP and partition MAP read it, and write the
    map of their clusters to their output directory (over it, with
    --inplace)."""
    return Path(shutil.copyfile(metagenome_network["map"], path))


def test_validation_splits_a_merged_bin(
    metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli, saved_figures
):
    """A bin holding two genomes is reported contaminated by (stub) miComplete:
    the recursive decontamination splits it back into the two genomes."""
    random.seed(0)
    calls = []
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome, calls))
    metacool = _map_with_merged_partition(
        metagenome, metagenome_network, tmp_path / "network.metacool", merged=("genome_1", "genome_2")
    )
    out = tmp_path / "out"
    # Outputs of an earlier run, replaced thanks to --force: its final
    # clusters, and its plots, but a file of the user's in plots/, left there.
    for stale, name in [(MAGS_WITH_MGES, "metator_00009_00000.fa"), ("plots", "bins_distribution.pdf")]:
        (out / stale).mkdir(parents=True)
        (out / stale / name).write_text(">old\nACGT\n")
    (out / "plots" / "notes.txt").write_text("the user's\n")
    args = [
        str(metacool),
        "-a", str(metagenome_inputs["assembly"]),
        "-o", str(out),
        "-T", str(tmp_path / "tmp"),
        "-s", "50000",
        "-i", "10",
        "-t", "1",
        "-F",
    ]  # fmt: skip
    run_cli("validation", *args)
    assert not (out / MAGS_WITH_MGES / "metator_00009_00000.fa").exists()

    # miComplete ran on the partition bins, then on the bins of the first
    # recursive step only: the split bins are clean.
    assert [Path(c[0]).name for c in calls] == ["overlapping_clusters", "step_1"]
    final = pd.read_csv(_kept(out) / "contig_data_final.txt", sep="\t")
    truth, found = _genome_labels(final, metagenome, "Final_cluster")
    assert "ND" not in found
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    assert synthetic.purity(truth, found) >= 0.9
    summary = pd.read_csv(_kept(out) / "cluster_summary.txt", sep="\t", index_col=0)
    assert len(summary) == len(metagenome.genomes)
    assert "metator_00001_00000" not in summary.index
    assert summary["Weighted redundancy"].max() == 1.0
    # Hi-C abundances are proportional to the hits of the bins.
    hits = final.groupby("Final_cluster").Hit.sum()[summary.index]
    hic = summary.HiC_abundance.astype(float)
    assert (hic / hic.sum()).to_numpy() == pytest.approx((hits / hits.sum()).to_numpy(), rel=1e-3)
    # The final fasta files are the final bins, and nothing else.
    assert sorted(p.name for p in (out / MAGS_WITH_MGES).iterdir()) == sorted(f"{name}.fa" for name in summary.index)
    binning = pd.read_csv(_kept(out) / "binning.txt", sep="\t", header=None, names=["Name", "Bin"])
    assert set(binning.Bin) == set(summary.index)
    assert sorted(binning.Name) == sorted(final.Name[final.Final_cluster != "ND"])
    # The plots of the final clusters, in plots/, as PDF files, the user's
    # file left there: a dot of completeness per cluster, and a wedge per
    # quality class and for the contigs in no cluster.
    plots = ["bins_distribution.pdf", "bins_size_distribution.pdf"]
    assert sorted(path.name for path in (out / "plots").iterdir()) == [*plots, "notes.txt"]
    assert all((out / "plots" / name).read_bytes().startswith(b"%PDF") for name in plots)
    assert (out / "plots" / "notes.txt").read_text() == "the user's\n"
    completion, _ = saved_figures["bins_distribution.pdf"].axes[0].collections
    dots = completion.get_offsets()[:, 1]
    assert len(dots) == len(summary) and list(dots) == sorted(dots, reverse=True)
    # miComplete's numbers, which the summary gives at two decimals.
    np.testing.assert_allclose(sorted(dots), sorted(100 * summary["Weighted completeness"]), atol=0.5)
    assert len(saved_figures["bins_size_distribution.pdf"].axes[0].patches) == 6
    # The bins of the split are recursive bins, of step 1; the others are not.
    assert list(summary.recursive) == list(summary.step > 0)
    assert summary.recursive.any() and not summary.recursive.all()
    # The map alone writes the summary again, byte for byte: the map of the
    # clusters, in the output directory.
    assert _summary_from_the_map(out / metacool.name) == (_kept(out) / "cluster_summary.txt").read_text()
    # metacooler wrote the FASTA files of the clusters of the recursive step,
    # and indexed nothing beside the assembly.
    assert not os.path.exists(f"{metagenome_inputs['assembly']}.fxi")


def test_validation_decompresses_the_assembly_once(
    metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli
):
    """An assembly compressed with gzip is decompressed once for the run, in its
    temporary directory: every FASTA file of the clusters (the overlapping
    clusters, each recursive step) is written from that copy, where metacooler
    copied it again for each write."""
    import gzip

    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    written = []
    write_fasta = metacooler.write_fasta

    def recorded(uri, fasta, out, **kwargs):
        written.append(Path(fasta))
        return write_fasta(uri, fasta, out, **kwargs)

    monkeypatch.setattr(metacooler, "write_fasta", recorded)
    compressed = tmp_path / "assembly.fa.gz"
    compressed.write_bytes(gzip.compress(Path(metagenome_inputs["assembly"]).read_bytes()))
    finals = {}
    for name, assembly in [("plain", metagenome_inputs["assembly"]), ("gz", compressed)]:
        metacool = _map_with_merged_partition(
            metagenome, metagenome_network, tmp_path / f"{name}.metacool", merged=("genome_1", "genome_2")
        )
        written.clear()
        random.seed(0)
        run_cli("validation", metacool, "-a", assembly, "-o", tmp_path / name, "-T", tmp_path / f"tmp_{name}",
                "-s", "50000", "-i", "10", "-t", "1")  # fmt: skip
        finals[name] = {path.name: path.read_bytes() for path in (tmp_path / name / MAGS_WITH_MGES).iterdir()}
    # The overlapping clusters, then the clusters of the recursive step.
    assert len(written) == 2 and len(set(written)) == 1
    (plain,) = set(written)
    assert plain.name == "assembly.fa" and plain.parent.parent == tmp_path / "tmp_gz"
    assert finals["gz"] == finals["plain"]


def test_validation_checks_the_clusters_in_the_order_of_the_partition(
    metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli
):
    """The FASTA file of each overlapping cluster that miComplete checks lists its
    contigs as metator partition lists them (see metator.partition.partition):
    by core cluster, in the order of their numbers, then by ID, as the files of
    overlapping_clusters/ did, not in the order of the map."""
    data = metagenome_network["contig_table"]
    host = data.Name.map(metagenome.truth.Host)
    # genome_1 and genome_2 in one overlapping cluster, the genome whose
    # contigs come last in the map in the first core cluster.
    late, early = sorted(["genome_1", "genome_2"], key=lambda genome: -data.ID[host == genome].min())
    others = [genome.name for genome in metagenome.genomes if genome.name not in (late, early)]
    core = {number + 1: sorted(data.ID[host == genome]) for number, genome in enumerate([late, early, *others])}
    overlapping = {1: core[1] + core[2], **{number: core[number + 1] for number in range(2, len(core))}}
    clusters = mtcl.Partition(
        core=core, overlapping=overlapping, names={number: mtcl.overlapping_name("metator", number) for number in overlapping}
    )
    assert overlapping[1] != sorted(overlapping[1])
    metacool = _network_map(metagenome_network, tmp_path / "network.metacool")
    entry = {mtmc.NETWORK_DIGEST: mtn.stored_network(str(metacool)).digest}
    mtmc.write_map(
        str(metacool),
        str(metacool),
        partition=entry,
        rounds=mtcl.partition_rounds(clusters, list(data.Name)),
        round_names=list(mtmc.PARTITION_ROUNDS),
        force=True,
    )
    assert mtmc.read_partition(str(metacool)) == clusters
    checked = {}
    fake = synthetic.fake_micomplete(metagenome)

    def micomplete_quality(fasta_dir, outfile, threads):
        if not checked:
            checked.update({path.stem: list(synthetic.read_fasta(path)) for path in Path(fasta_dir).glob("*.fa")})
        fake(fasta_dir, outfile, threads)

    monkeypatch.setattr("metator.validation.micomplete_quality", micomplete_quality)
    run_cli("validation", metacool, "-a", metagenome_inputs["assembly"], "-o", tmp_path / "out", "-T", tmp_path / "tmp",
            "-s", "50000", "-i", "5")  # fmt: skip
    names = list(data.Name)
    expected = {clusters.names[number]: [names[contig - 1] for contig in contigs] for number, contigs in overlapping.items()}
    assert checked and checked == {name: expected[name] for name in checked}
    assert checked[clusters.names[1]] != [names[contig - 1] for contig in sorted(overlapping[1])]


def test_validation_again_on_the_map_of_its_clusters(
    metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli
):
    """validation MAP runs again on the map of the clusters that it wrote, in
    which it grouped the contigs by cluster: the partition stored in the map is
    still that of its network, which the regroup left as it was (see
    metator.network.network_digest), the same contig for contig, and the run
    reads it with the IDs of the contigs in the regrouped map, their positions
    in it, which number the rows of its tables."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    metacool = _map_with_merged_partition(
        metagenome, metagenome_network, tmp_path / "network.metacool", merged=("genome_1", "genome_2")
    )

    def named(path):
        """The partition of the map at path, by contig name."""
        names = list(metacooler.MetaCooler(str(path)).contigs.name)
        partition = mtmc.read_partition(str(path), network=mtn.stored_network(str(path)).digest)
        return [
            {number: sorted(names[contig - 1] for contig in contigs) for number, contigs in clusters.items()}
            for clusters in (partition.core, partition.overlapping)
        ], names

    before, names = named(metacool)
    # In place: the second run validates the map of the clusters of the first.
    common = ["-a", metagenome_inputs["assembly"], "-T", tmp_path / "tmp", "-s", "50000", "-i", "10", "-t", "1", "--inplace"]
    finals = {}
    for run in ["first", "second"]:
        random.seed(0)
        run_cli("validation", metacool, "-o", tmp_path / run, *common)
        finals[run] = pd.read_csv(_kept(tmp_path / run) / "contig_data_final.txt", sep="\t")
    after, regrouped = named(metacool)
    assert regrouped != names and sorted(regrouped) == sorted(names)
    assert after == before
    # The second run numbers the contigs as the regrouped map orders them.
    assert list(finals["first"].Name) == names and list(finals["second"].Name) == regrouped
    assert list(finals["second"].ID) == list(range(1, len(regrouped) + 1))
    found = {run: final.set_index("Name").Final_cluster.sort_index() for run, final in finals.items()}
    assert found["second"].equals(found["first"])


def test_contigs_named_by_numbers(tmp_path, monkeypatch, run_cli):
    """Contigs named by numbers (000017, 1e5, 2.50), which pandas reads as numbers
    unless told otherwise: a run reads every table of names as text, and each
    contig goes by its name, as the map names it, to the final cluster of the
    FASTA file that lists it (read as numbers, the names lost their leading
    zeros, or their form, and the clusters their contigs)."""
    from tests.conftest import SYNTHETIC_SEED

    metagenome = synthetic.make_metagenome(seed=SYNTHETIC_SEED)
    renamed = {}
    for i, contig in enumerate(metagenome.contigs):
        renamed[contig.name] = {0: "1e5", 1: "2.50"}.get(i, f"{i:06d}")
        contig.name = renamed[contig.name]
    for column in ["chr1", "chr2"]:
        metagenome.pairs[column] = metagenome.pairs[column].map(renamed)
    files = metagenome.write_all(tmp_path / "inputs")
    metacool = contig_map(
        files["pairs"], mtn.create_contig_data(str(files["assembly"])), tmp_path / "input.metacool", columns=False
    )
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    out = tmp_path / "out"
    common = ["-o", out, "-T", tmp_path / "tmp"]
    random.seed(0)
    run_cli("network", metacool, "-a", files["assembly"], "-e", "DpnII", *common)
    # The map that network wrote in the output directory, which the next
    # steps update in place.
    metacool = common[1] / metacool.name
    run_cli("partition", metacool, "-i", "10", *common)
    run_cli("validation", metacool, "-a", files["assembly"], "-s", "50000", "-i", "5", *common)
    final = pd.read_csv(_kept(out) / "contig_data_final.txt", sep="\t", dtype={"Name": str}, keep_default_na=False)
    assert sorted(final.Name) == sorted(renamed.values())
    mc = metacooler.MetaCooler(str(metacool))
    genome_of = dict(zip(mc.contigs["name"], mc.genomes["name"].to_numpy()[mc.contigs["genome_id"]], strict=True))
    binned = final[final.Final_cluster != "ND"]
    assert len(binned) and {name: genome_of[name] for name in binned.Name} == dict(
        zip(binned.Name, binned.Final_cluster, strict=True)
    )
    # Each FASTA file of a final cluster lists the contigs of its cluster in
    # the map.
    files = sorted((out / MAGS_WITH_MGES).glob("*.fa"))
    assert sorted(path.stem for path in files) == sorted(set(binned.Final_cluster))
    for path in files:
        listed = [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]
        assert sorted(listed) == sorted(name for name, genome in genome_of.items() if genome == path.stem)


def test_a_contig_named_na(tmp_path, monkeypatch, run_cli):
    """A contig named NA, which pandas reads as a missing value unless told
    otherwise, is a contig as any other: validation reads it by name from the
    partition table of the run, and writes it into the map of the clusters, in
    its final cluster, with its size (a missing name made the Size column of
    contig_data_final.txt a float, which the map refused; the run keeps its
    tables in its temporary directory: see _keep_the_tables)."""
    from tests.conftest import SYNTHETIC_SEED

    metagenome = synthetic.make_metagenome(seed=SYNTHETIC_SEED)
    contig = max(metagenome.molecule("genome_1").contigs, key=lambda contig: contig.length)
    name, contig.name = contig.name, "NA"
    for column in ["chr1", "chr2"]:
        metagenome.pairs[column] = metagenome.pairs[column].replace(name, "NA")
    files = metagenome.write_all(tmp_path / "inputs")
    metacool = contig_map(
        files["pairs"], mtn.create_contig_data(str(files["assembly"])), tmp_path / "input.metacool", columns=False
    )
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    common = ["-o", tmp_path / "out", "-T", tmp_path / "tmp"]
    random.seed(0)
    run_cli("network", metacool, "-a", files["assembly"], "-e", "DpnII", *common)
    # The map that network wrote in the output directory, which the next
    # steps update in place.
    metacool = common[1] / metacool.name
    run_cli("partition", metacool, "-i", "10", *common)
    run_cli("validation", metacool, "-a", files["assembly"], "-s", "50000", "-i", "5", *common)
    final = pd.read_csv(_kept(tmp_path / "out") / "contig_data_final.txt", sep="\t", keep_default_na=False)
    row = final[final.Name == "NA"]
    assert len(row) == 1 and pd.api.types.is_integer_dtype(final.Size) and int(row.Size.iloc[0]) == contig.length
    mc = metacooler.MetaCooler(str(metacool))
    members = mc.cluster_members("final")
    assert list(members.cluster[members.contig == "NA"]) == [row.Final_cluster.iloc[0]]


def test_validation_without_contaminated_bin(metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli):
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    metacool = _map_with_merged_partition(metagenome, metagenome_network, tmp_path / "network.metacool", merged=())
    out = tmp_path / "out"
    args = [
        str(metacool),
        "-a", str(metagenome_inputs["assembly"]),
        "-o", str(out),
        "-T", str(tmp_path / "tmp"),
        "-s", "50000",
    ]  # fmt: skip
    run_cli("validation", *args)
    assert len(pd.read_csv(_kept(out) / "cluster_summary.txt", sep="\t", index_col=0)) == len(metagenome.genomes)


def test_the_plots_count_the_clusters_as_the_log(
    metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli, saved_figures
):
    """validation draws the plots of its final clusters from miComplete's numbers,
    which cluster_summary.txt and the map round to two decimals: a cluster of
    a weighted completeness of 0.8996 is a MQ MAG in the pie chart, as in the
    log, where its summary says 0.90, the completeness of a HQ MAG."""

    def micomplete_quality(fasta_dir, outfile, threads):
        synthetic.write_micomplete_summary(fasta_dir, outfile, metagenome)
        # Three comment lines, the header, then a row per cluster.
        with open(outfile) as table:
            lines = [line.split("\t") for line in table.read().splitlines()]
        lines[4][lines[3].index("Weighted completeness")] = "0.8996"
        with open(outfile, "w") as table:
            table.write("".join("\t".join(line) + "\n" for line in lines))

    monkeypatch.setattr("metator.validation.micomplete_quality", micomplete_quality)
    metacool = _map_with_merged_partition(metagenome, metagenome_network, tmp_path / "network.metacool", merged=())
    out = tmp_path / "out"
    run_cli("validation", metacool, "-a", metagenome_inputs["assembly"], "-o", out, "-T", tmp_path / "tmp", "-s", "50000")
    summary = pd.read_csv(_kept(out) / "cluster_summary.txt", sep="\t", index_col=0, dtype=str)
    (almost,) = summary.index[summary["Weighted completeness"] == "0.90"]
    length = int(float(summary.loc[almost, "Length"]))
    others = summary.drop(index=almost)
    assert (others["Weighted completeness"] == "1.00").all() and (summary["Weighted redundancy"] == "1.00").all()
    legend = [text.get_text() for text in saved_figures["bins_size_distribution.pdf"].axes[0].get_legend().get_texts()]
    assert legend[0].startswith(f"HQ MAGs: {len(others)} - ")
    assert legend[1] == f"MQ MAGs: 1 - {round(length / 1_000_000, 2)}Mb"
    (log,) = out.glob("metator_validation_*.log")
    assert f"MQ MAGs: 1\tTotal Size: {length}" in log.read_text()
    completion, _ = saved_figures["bins_distribution.pdf"].axes[0].collections
    assert completion.get_offsets()[-1, 1] == pytest.approx(89.96)


def test_a_run_goes_on_when_a_plot_cannot_be_drawn(
    metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli
):
    """A plot that cannot be drawn is a warning in the log, with its error and its
    traceback: the run draws the other one, writes the map of its clusters,
    and exits 0."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))

    def fails(*args):
        raise ValueError("no pie")

    monkeypatch.setattr("metator.figures.figures_bins_size_distribution", fails)
    metacool = _map_with_merged_partition(metagenome, metagenome_network, tmp_path / "network.metacool", merged=())
    out = tmp_path / "out"
    run_cli("validation", metacool, "-a", metagenome_inputs["assembly"], "-o", out, "-T", tmp_path / "tmp", "-s", "50000")
    assert [path.name for path in (out / "plots").iterdir()] == ["bins_distribution.pdf"]
    metacooler.validate(str(out / metacool.name))
    (log,) = out.glob("metator_validation_*.log")
    logged = log.read_text()
    assert f"WARNING :: {out / 'plots' / 'bins_size_distribution.pdf'} is not drawn: ValueError: no pie" in logged
    assert "Traceback (most recent call last)" in logged
    assert "Final contact map written after " in logged and logged.rstrip().endswith(str(out / metacool.name))


def test_validation_without_bin_above_the_size_threshold(metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli):
    """partition MAP stores every overlapping cluster, whatever its size; a
    validation run of which none reaches --size stops before miComplete, saying
    why, and leaves its map as it was."""
    calls = []
    monkeypatch.setattr("metator.validation.micomplete_quality", lambda *args: calls.append(args))
    metacool = _network_map(metagenome_network, tmp_path / "network.metacool")
    run_cli("partition", str(metacool), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), "-i", "2")
    # The map of the output directory, with the partition.
    metacool = tmp_path / "out" / metacool.name
    mc = metacooler.MetaCooler(str(metacool))
    assert len(mc.clusters[mc.clusters["round"] == "overlapping"]) > 0
    before = metacool.read_bytes()
    args = [str(metacool), "-a", str(metagenome_inputs["assembly"]), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"),
            "-s", "10000000"]  # fmt: skip
    with pytest.raises(ValueError):
        run_cli("validation", *args)
    assert not calls
    (log,) = (tmp_path / "out").glob("metator_validation_*.log")
    assert "ERROR :: No overlapping cluster of the partition reaches the minimum size" in log.read_text()
    assert metacool.read_bytes() == before


def test_pipeline_without_bin_above_the_size_threshold(metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """No bin of the partition reaches --size: the run stops before miComplete
    (which found no bin, then failed to read its empty results), saying why in
    its log file."""
    calls = []
    monkeypatch.setattr("metator.validation.micomplete_quality", lambda *args: calls.append(args))
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    before = metacool.read_bytes()
    out = tmp_path / "out"
    args = [str(metacool), "-a", str(metagenome_inputs["assembly"]), "-s", "10000000", "-i", "2", "-o", str(out),
            "-T", str(tmp_path / "tmp")]  # fmt: skip
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)
    assert not calls
    (log,) = out.glob("metator_pipeline_*.log")
    message = (
        "No overlapping cluster of the partition reaches the minimum size of a cluster to check (--size of validation or "
        "pipeline): there is no cluster to validate."
    )
    assert f"ERROR :: {message}" in log.read_text()
    assert metacool.read_bytes() == before


def test_pipeline_recovers_the_genomes(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """The whole pipeline from the map of the pairs, with stub miComplete."""
    random.seed(0)
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    metacool = _contig_map(metagenome_inputs, tmp_path / "inputs" / "contigs.metacool")
    out = tmp_path / "out"
    args = [
        str(metacool),
        "-a", str(metagenome_inputs["assembly"]),
        "-d", str(metagenome_inputs["depth"]),
        "-e", "DpnII",
        "-m", "abundance",
        "-s", "50000",
        "-i", "10",
        "-j", "5",
        "-o", str(out),
        "-T", str(tmp_path / "tmp"),
    ]  # fmt: skip
    run_cli("pipeline", *args)

    final = pd.read_csv(_kept(out) / "contig_data_final.txt", sep="\t")
    truth, found = _genome_labels(final, metagenome, "Final_cluster")
    assert "ND" not in found
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    summary = pd.read_csv(_kept(out) / "cluster_summary.txt", sep="\t", index_col=0)
    assert len(summary) == len(metagenome.genomes)
    # MGE contigs are binned with their host.
    host_bin = final.set_index("Name").Final_cluster
    for mge in metagenome.mges:
        host_contigs = [c.name for c in metagenome.molecule(mge.host).contigs]
        assert {host_bin[c.name] for c in mge.contigs} <= set(host_bin[host_contigs])
    # The shotgun abundance of a bin is its share of the sequenced bases.
    truth = metagenome.truth
    bases = (truth.Length * truth.Abundance).groupby(final.set_index("Name").Final_cluster).sum()
    expected = 100 * bases / (truth.Length * truth.Abundance).sum()
    assert summary.SG_abundance.astype(float).to_numpy() == pytest.approx(expected[summary.index].to_numpy(), rel=0.1)
    # The contact map of the clusters, in the output directory: one bin per
    # contig, one genome per cluster (the cluster summary columns with it),
    # and every pair a contact.
    written = out / metacool.name
    metacooler.validate(str(written))
    clr = cooler.Cooler(str(written))
    assert clr.pixels()[:]["count"].sum() == len(metagenome.pairs)
    assert clr.info["metadata"]["metator"]["validation"]["parameters"]["map"] == str(metacool)
    mc = metacooler.MetaCooler(str(written))
    assert sorted(mc.contigs.name) == sorted(final.Name)
    assert dict(zip(mc.contigs.name, mc.contigs.genome_id.map(mc.genomes.name), strict=True)) == dict(host_bin)
    assert sorted(mc.genomes.name) == sorted(summary.index)
    assert mc.genomes.set_index("name").Completeness[summary.index].to_numpy() == pytest.approx(summary.Completeness.to_numpy())
    # What metaTOR measured on each contig (with a depth file and an
    # enzyme: the four columns), not its bins.
    columns = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    assert [column for column in mc.contigs.columns if not column.startswith("qc_")] == ["name", "length", "genome_id", *columns]
    measured = final.set_index("Name").loc[mc.contigs.name, columns]
    for column in columns:
        assert mc.contigs[column].tolist() == measured[column].tolist(), column
    # Clean-up: no pyfastx index, no intermediate table in the output
    # directory (the tables of the validation stay in the temporary
    # directory, kept here as -N keeps it: see _keep_the_tables, and metator
    # dump writes them in tables/); the FASTA files of the final clusters,
    # with their MGE contigs, and of the MAGs, the cool files of the MAGs,
    # the plots of the clusters, as PDF files, in plots/, and the report.
    assert not os.path.exists(f"{metagenome_inputs['assembly']}.fxi")
    outputs = sorted(path.name for path in out.iterdir() if not path.name.endswith(".log"))
    assert outputs == sorted([metacool.name, "cool", "fasta", "plots", "tables", *mtr.REPORT_FILES])
    assert sorted(path.name for path in (out / "fasta").iterdir()) == ["mags", "mags_with-mges"]
    assert sorted(path.name for path in (out / "plots").iterdir()) == ["bins_distribution.pdf", "bins_size_distribution.pdf"]
    assert sorted(path.name for path in (out / "tables").iterdir()) == sorted(mtd.TABLES)
    assert not (out / "contig_data_partition.txt").exists()


def test_pipeline_map_is_the_map_of_the_pairs(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """The map that a run leaves is the map of the pairs of its input map and its
    final tables (see tests.maps.pairs_map): the same bins, contacts, contigs
    table (column for column), genomes table (the columns of the cluster
    summary and the contact statistics of each cluster) and QC block."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    from_pairs = metacooler.from_pairs
    maps = []

    def counted_from_pairs(uri, *args, **kwargs):
        maps.append(str(uri))
        return from_pairs(uri, *args, **kwargs)

    metacool = _contig_map(metagenome_inputs, tmp_path / "inputs" / "contigs.metacool")
    monkeypatch.setattr(metacooler, "from_pairs", counted_from_pairs)
    out = tmp_path / "out"
    args = [str(metacool), "-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]),
            "-e", "DpnII", "-s", "50000", "-i", "10", "-j", "5", "-o", str(out), "-T", str(tmp_path / "tmp")]  # fmt: skip
    random.seed(0)
    run_cli("pipeline", *args)
    assert maps == []
    # The map of the clusters, in the output directory.
    written = out / metacool.name
    metacooler.validate(str(written), full=True)
    clr = cooler.Cooler(str(written))
    expected = tmp_path / "pairs.metacool"
    pairs_map(
        metagenome_inputs["pairs"],
        _kept(out) / "contig_data_final.txt",
        expected,
        bin_summary=str(_kept(out) / "cluster_summary.txt"),
        assembly_name=clr.info["genome-assembly"],
        metadata={key: value for key, value in clr.info["metadata"]["metator"].items() if key != "version"},
    )

    ref = cooler.Cooler(str(expected))
    bins = clr.bins()[:]
    # The run's map also stores the normalization of its network, which the
    # map of the pairs does not hold.
    assert "metator_empirical_hit" in bins.columns
    pd.testing.assert_frame_equal(bins.drop(columns="metator_empirical_hit", errors="ignore"), ref.bins()[:])
    for table in ["pixels", "chroms"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    _assert_normalized_as_the_network(written, tmp_path)
    genomes = metacooler.MetaCooler(str(written)).genomes
    pd.testing.assert_frame_equal(genomes, metacooler.MetaCooler(str(expected)).genomes)
    assert {"Completeness", "qc_inter_genome"} <= set(genomes.columns)
    info, ref_info = dict(clr.info), dict(ref.info)
    metadata, ref_metadata = info.pop("metadata"), ref_info.pop("metadata")
    assert metadata["qc"] == ref_metadata["qc"] and metadata["qc"]["n_pairs_read"] == len(metagenome.pairs)
    assert {key for key in info.keys() | ref_info.keys() if info.get(key) != ref_info.get(key)} <= {"creation-date"}
    differ = {key for key in metadata.keys() | ref_metadata.keys() if metadata.get(key) != ref_metadata.get(key)}
    assert differ == {"annotate", "metator", "reassign"}
    assert [entry["source"] for entry in metadata["reassign"]] == [str(metacool)]
    # The columns of the input map (tests.maps.contig_map), then those that
    # the run stores, once, as metator network stores them (dropped, then
    # added).
    assert len(metadata["annotate"]) == 3 and len(ref_metadata["annotate"]) == 1


def _assert_normalized_as_the_network(written, tmp_path):
    """The contacts of the map of a run normalized by the weights it stores and
    summed by pair of contigs (metacooler network --balance) are the edges of
    its network (see metator.network.map_network: but those of weight 0),
    whatever the order of its contigs."""
    names = [str(name) for name in metacooler.MetaCooler(str(written)).contigs["name"]]
    network = mtn.map_network(str(written))
    expected_edges = {
        tuple(sorted((names[id1 - 1], names[id2 - 1]))): weight for id1, id2, weight in network.itertuples(index=False)
    }
    edges = tmp_path / "edges.tsv"
    metacooler.network(str(written), level="contig", balance="metator_empirical_hit", out=str(edges), header=True)
    normalized = {
        tuple(sorted((first, second))): value
        for first, second, value in pd.read_csv(edges, sep="\t").itertuples(index=False)
        if first != second
    }
    assert normalized.keys() == expected_edges.keys()
    assert all(normalized[pair] == pytest.approx(expected_edges[pair], rel=1e-9) for pair in normalized)


def test_pipeline_restarts_from_the_map_of_a_stopped_run(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """A run that stops with an error (miComplete fails in the validation, here)
    leaves its map as it was, in no cluster, even one that it updates in place
    (--inplace), from which a run starts again."""
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-T", str(tmp_path / "tmp"), "--inplace"]  # fmt: skip
    clean, out = tmp_path / "clean", tmp_path / "out"
    clean_map = _contig_map(metagenome_inputs, tmp_path / "clean.metacool")
    metacool = Path(shutil.copyfile(clean_map, tmp_path / "contigs.metacool"))
    before = metacool.read_bytes()
    micomplete = synthetic.fake_micomplete(metagenome)
    monkeypatch.setattr("metator.validation.micomplete_quality", micomplete)
    random.seed(0)
    run_cli("pipeline", str(clean_map), "-o", str(clean), *common)

    def micomplete_fails(*args):
        raise RuntimeError("miComplete failed")

    monkeypatch.setattr("metator.validation.micomplete_quality", micomplete_fails)
    random.seed(0)
    with pytest.raises(RuntimeError, match=r"^miComplete failed$"):
        run_cli("pipeline", str(metacool), "-o", str(out), *common)
    assert metacool.read_bytes() == before
    monkeypatch.setattr("metator.validation.micomplete_quality", micomplete)
    random.seed(0)
    run_cli("pipeline", str(metacool), "-o", str(out), "-F", *common)

    for name in ["contig_data_final.txt", "cluster_summary.txt", "binning.txt"]:
        assert (_kept(out) / name).read_text() == (_kept(clean) / name).read_text(), name
    final_bins = [{path.name: path.read_bytes() for path in (run / MAGS_WITH_MGES).iterdir()} for run in [out, clean]]
    assert final_bins[0] == final_bins[1] and final_bins[0]
    # The map of the clusters, written over the map that the stopped run
    # left as it was.
    metacooler.validate(str(metacool), full=True)
    clr, ref = cooler.Cooler(str(metacool)), cooler.Cooler(str(clean_map))
    for table in ["bins", "pixels", "chroms"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(metacool)).genomes, metacooler.MetaCooler(str(clean_map)).genomes)
    info, ref_info = dict(clr.info), dict(ref.info)
    metadata, ref_metadata = info.pop("metadata"), ref_info.pop("metadata")
    assert metadata["qc"] == ref_metadata["qc"]
    assert {key for key in info.keys() | ref_info.keys() if info.get(key) != ref_info.get(key)} <= {"creation-date"}
    # The parameters of the run (its map, its output directory) and the
    # paths of the history.
    differ = {key for key in metadata.keys() | ref_metadata.keys() if metadata.get(key) != ref_metadata.get(key)}
    assert differ == {"metator", "annotate", "reassign"}
    assert [entry["source"] for entry in metadata["reassign"]] == [str(metacool)]


def _assert_same_clusters(out, ref, contigs=None):
    """The outputs of two runs (out, ref: their output directories) of the same
    clusters from maps whose contigs are in another order, or that lack some of
    them (contigs: the contigs of the map of out, those of ref by default): the
    IDs of the contigs are their positions in the map of the run, and so are
    the rows of its tables, and the sums over the contigs of a cluster (its
    abundances) are made in that order, so equal to the last bits."""

    def table(run, name, **kwargs):
        return pd.read_csv(_kept(run) / name, sep="\t", **kwargs)

    binning = [table(run, "binning.txt", header=None, names=["Name", "Bin"]) for run in (out, ref)]
    if contigs is not None:
        binning[1] = binning[1][binning[1].Name.isin(contigs)]
    assert sorted(map(tuple, binning[0].to_numpy())) == sorted(map(tuple, binning[1].to_numpy()))
    final = [
        table(run, "contig_data_final.txt", index_col=0, dtype=str).drop(columns="ID").set_index("Name") for run in (out, ref)
    ]
    if contigs is not None:
        final[1] = final[1].loc[sorted(contigs)]

    def values(frame):
        # The values as numbers where they all are: the IDs of the clusters
        # of a column in which no contig is "-" are written as numbers
        # (00001 as 1), in a table of the contigs of the map.
        frame = frame.sort_index()
        for column in frame.columns:
            numbers = pd.to_numeric(frame[column], errors="coerce")
            if numbers.notna().all():
                frame[column] = numbers
        return frame

    pd.testing.assert_frame_equal(values(final[0]), values(final[1]))
    for name, index in [("cluster_summary.txt", 0), ("clusters.txt", ["round", "name"])]:
        found = [table(run, name, index_col=index).sort_index() for run in (out, ref)]
        if contigs is not None:
            found = [frame.drop(columns="SG_abundance") for frame in found]
        pd.testing.assert_frame_equal(found[0], found[1], check_exact=False, rtol=1e-12)
    members = [table(run, "cluster_members.txt").sort_values(["round", "cluster", "contig"]) for run in (out, ref)]
    pd.testing.assert_frame_equal(members[0].reset_index(drop=True), members[1].reset_index(drop=True))
    fasta = [{path.name: synthetic.read_fasta(path) for path in (run / MAGS_WITH_MGES).iterdir()} for run in (out, ref)]
    assert fasta[0] == fasta[1] and fasta[0]


def test_pipeline_from_a_metacool(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """From the contact map of the clusters of a first run (one bin per contig,
    its contigs grouped by cluster), the pipeline makes the same network and
    the same clusters, but for the order of the contigs, which is that of the
    map (see _assert_same_clusters)."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    # In place (--inplace): the first run makes its map the map of its
    # clusters, which the second starts from.
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-T", str(tmp_path / "tmp"), "--inplace"]  # fmt: skip
    expected = _contig_map(metagenome_inputs, tmp_path / "first.metacool")
    random.seed(0)
    run_cli("pipeline", str(expected), "-o", str(tmp_path / "pair"), *common)
    metacool = Path(shutil.copyfile(expected, tmp_path / "input.metacool"))
    random.seed(0)
    run_cli("pipeline", str(metacool), "-o", str(tmp_path / "metacool"), *common)

    pair, mcool = tmp_path / "pair", tmp_path / "metacool"
    _assert_same_clusters(mcool, pair)
    final = pd.read_csv(_kept(mcool) / "contig_data_final.txt", sep="\t")
    truth, found = _genome_labels(final, metagenome, "Final_cluster")
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    # The same network, stored in both maps: the same edges, by contig name.
    edges = []
    for path in [expected, metacool]:
        names = [str(name) for name in metacooler.MetaCooler(str(path)).contigs["name"]]
        network = mtn.map_network(str(path))
        edges.append({tuple(sorted((names[a - 1], names[b - 1]))): w for a, b, w in network.itertuples(index=False)})
    assert edges[0] == edges[1] and edges[0]
    # The input map, updated in place: no map goes to the output directory.
    assert not list(mcool.glob("*.metacool*")) and not list(tmp_path.glob("*.lock"))
    assert not list(pair.glob("*.metacool*"))
    written = metacool
    metacooler.validate(str(written), full=True)
    _read_by_cooler_and_metacooler(written)
    clr, ref = cooler.Cooler(str(written)), cooler.Cooler(str(expected))
    for table in ["bins", "pixels"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    # The contigs table of the input map: metacooler reassign relabels its
    # contigs in place (the same bins), and writes genome_id again, after the
    # other columns; metaTOR then stores its own columns last, as in the map
    # of the first run.
    chroms, ref_chroms = clr.chroms()[:], ref.chroms()[:]
    columns = ["name", "length", "genome_id", "GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    for table in [chroms, ref_chroms]:
        assert [column for column in table.columns if not column.startswith("qc_")] == columns
    pd.testing.assert_frame_equal(chroms[ref_chroms.columns], ref_chroms)
    # The QC block and columns of the first run, carried over.
    genomes, ref_genomes = metacooler.MetaCooler(str(written)).genomes, metacooler.MetaCooler(str(expected)).genomes
    pd.testing.assert_frame_equal(genomes, ref_genomes)
    metadata = clr.info["metadata"]
    parameters = metadata["metator"]["validation"]["parameters"]
    assert parameters["map"] == str(metacool) and "start" not in parameters and "forward" not in parameters
    assert metadata["qc"] == ref.info["metadata"]["qc"] and metadata["qc"]["source"] == "pairs"
    # The history of the input map: the map of the clusters of the first
    # run, written over its map of the contigs, then that of this run.
    assert [entry["source"] for entry in metadata["reassign"]] == [str(expected), str(metacool)]
    # The clustering rounds of the run, those of the first run: the core and
    # overlapping clusters of the partition, each recursive step, and the
    # final clusters, the map's genomes.
    mc, ref_mc = metacooler.MetaCooler(str(written)), metacooler.MetaCooler(str(expected))
    assert mc.cluster_rounds == ref_mc.cluster_rounds
    names = [entry["name"] for entry in mc.cluster_rounds]
    assert names[:2] == ["core", "overlapping"] and names[-1] == "final"
    assert all(name == f"recursive_{step}" for step, name in enumerate(names[2:-1], start=1))
    pd.testing.assert_frame_equal(mc.clusters, ref_mc.clusters)
    pd.testing.assert_frame_equal(_members(mc), _members(ref_mc))
    # Every overlapping cluster, those that reach --size with the columns of
    # miComplete, which checked them; no FASTA file of them in the output
    # directory.
    final = pd.read_csv(_kept(mcool) / "contig_data_final.txt", sep="\t", dtype=str)
    overlapping = mc.clusters[mc.clusters["round"] == "overlapping"]
    assert set(overlapping.name) == {f"metator_{number}_00000" for number in set(final.Overlapping_cluster_ID) - {"-"}}
    checked = overlapping["Weighted completeness"].notna()
    assert checked.any() and (checked == (overlapping["length"] >= 50000)).all()
    assert not (mcool / "overlapping_clusters").exists()
    genome_of = dict(zip(mc.contigs.name, mc.contigs.genome_id.map(mc.genomes.name), strict=True))
    final_members = mc.cluster_members("final")
    assert set(final_members.cluster) == set(mc.genomes.name) - {"__unbinned__"}
    assert all(genome_of[contig] == cluster for cluster, contig in zip(final_members.cluster, final_members.contig, strict=True))

    # The same map, written from the map and the tables of the run, without
    # rounds (metator.metacool.write_map, with the columns of the contig
    # data, as metator metacool of metaTOR 1.4.0 wrote it from a map).
    command = tmp_path / "command.metacool"
    table = pd.read_csv(
        _kept(mcool) / "contig_data_final.txt", sep="\t", dtype={"Name": str, "Final_cluster": str}, keep_default_na=False
    )
    final = table["Final_cluster"].str.strip()
    binned = ~final.isin(mtmc.UNBINNED)
    summary = str(_kept(mcool) / "cluster_summary.txt")
    clusters = mtmc.FinalClusters(dict(zip(table["Name"][binned], final[binned], strict=True)), summary, None, {})
    mtmc.write_map(str(expected), str(command), network=mtmc.Network(table, "None"), clusters=clusters)
    metacooler.validate(str(command), full=True)
    clr = cooler.Cooler(str(command))
    run_map = cooler.Cooler(str(written))
    # But the normalization of the network, which the run stores in its map.
    pd.testing.assert_frame_equal(clr.bins()[:], run_map.bins()[:].drop(columns="metator_empirical_hit", errors="ignore"))
    for table in ["pixels", "chroms"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(run_map, table)()[:])
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(command)).genomes, genomes)
    assert clr.info["metadata"]["qc"] == metadata["qc"]


def test_pipeline_from_the_columns_of_a_metacool(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """From a map of the contigs whose contigs table holds what metaTOR measures
    with a depth file and an enzyme (see tests.maps.contig_map), a run without
    --depth and --enzyme reads the restriction sites and the coverage of the
    contigs from its Restriction_site and Shotgun_coverage columns
    (theoritical_hit needs both): the same network, tables and clusters as with
    them, byte for byte, and the map of its clusters stores them back."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    common = ["-a", str(metagenome_inputs["assembly"]), "-m", "theoritical_hit", "-s", "50000", "-i", "10", "-j", "5",
              "-T", str(tmp_path / "tmp")]  # fmt: skip
    options = ["-d", str(metagenome_inputs["depth"]), "-e", "DpnII"]
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    columns = metacooler.MetaCooler(str(metacool)).contigs[["name", "Restriction_site", "Shotgun_coverage"]]
    runs = {"options": options, "map": []}
    for name, given in runs.items():
        random.seed(0)
        # Each run writes its own map, in its output directory: both start
        # from the same map.
        run_cli("pipeline", str(metacool), "-o", str(tmp_path / name), *common, *given)  # fmt: skip

    by_options, by_map = tmp_path / "options", tmp_path / "map"
    for name in ["contig_data_final.txt", "cluster_summary.txt", "binning.txt"]:
        assert (_kept(by_map) / name).read_bytes() == (_kept(by_options) / name).read_bytes(), name
    # The same network, stored in each map.
    networks = [mtn.map_network(str(run / metacool.name)) for run in (by_options, by_map)]
    pd.testing.assert_frame_equal(networks[0], networks[1])
    assert "SG_abundance" in pd.read_csv(_kept(by_map) / "cluster_summary.txt", sep="\t").columns
    (log,) = by_map.glob("metator_pipeline_*.log")
    assert f"INFO :: The shotgun coverage of the contigs: the Shotgun_coverage column of {metacool}." in log.read_text()
    # Stored back: the map of the bins holds the values that the run read.
    for written in [by_map / metacool.name, by_options / metacool.name]:
        metacooler.validate(str(written), full=True)
        stored = metacooler.MetaCooler(str(written)).contigs.set_index("name")
        pd.testing.assert_frame_equal(stored.loc[columns.name, columns.columns[1:]], columns.set_index("name"))


def test_pipeline_from_a_filtered_metacool(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """From a map that lacks contigs of the assembly (made with a contig filter,
    which leaves out the pairs with a read on the contigs it drops), the
    pipeline completes and makes the clusters of a run from the map of the
    pairs that the map holds, every contig in it: the contigs it lacks have no
    hits, and are in no cluster and in no table of the run (see
    _assert_same_clusters)."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    truth = metagenome.truth
    shortest = int(truth.Length.median())
    kept = set(truth.Name[truth.Length >= shortest])
    # The map, at 1 kb, of the contigs of shortest bp or more. metacooler
    # indexes the assembly (.fai): it reads a private copy.
    metacool = tmp_path / "filtered.metacool"
    metacooler.from_pairs(str(metacool), pairs=str(metagenome_inputs["pairs"]), binsize=1000, contig_genome_map=None,
                          fasta=str(shutil.copyfile(metagenome_inputs["assembly"], tmp_path / "map.fa")), ps=False,
                          filters=[f"contigs.length >= {shortest}"])  # fmt: skip
    assert set(metacooler.MetaCooler(str(metacool)).contigs.name) == kept
    # The pairs it holds: those whose two reads lie on the contigs it keeps.
    pairs = tmp_path / "kept.pairs"
    with open(metagenome_inputs["pairs"]) as all_pairs, open(pairs, "w") as kept_pairs:
        for line in all_pairs:
            fields = line.split("\t")
            if line.startswith("#") or (fields[1] in kept and fields[3] in kept):
                kept_pairs.write(line)
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-T", str(tmp_path / "tmp")]  # fmt: skip
    random.seed(0)
    kept_map = _contig_map(metagenome_inputs, tmp_path / "kept.metacool", pairs=pairs)
    run_cli("pipeline", str(kept_map), "-o", str(tmp_path / "pair"), *common)
    contacts = cooler.Cooler(str(metacool)).pixels()[:]["count"].sum()
    random.seed(0)
    run_cli("pipeline", str(metacool), "-o", str(tmp_path / "metacool"), *common)

    pair, mcool = tmp_path / "pair", tmp_path / "metacool"
    _assert_same_clusters(mcool, pair, contigs=kept)
    final = pd.read_csv(_kept(mcool) / "contig_data_final.txt", sep="\t")
    assert set(final.Name) == kept
    # The contigs that the map lacks: no hit, in no cluster, in the run from
    # the map of all the contigs.
    lacking = set(truth.Name) - kept
    pair_final = pd.read_csv(_kept(pair) / "contig_data_final.txt", sep="\t").set_index("Name")
    assert 0 < len(lacking) < len(truth)
    assert (pair_final.Hit[sorted(lacking)] == 0).all() and (pair_final.Final_cluster[sorted(lacking)] == "ND").all()
    assert len(pd.read_csv(_kept(mcool) / "cluster_summary.txt", sep="\t")) > 0
    (log,) = mcool.glob("metator_pipeline_*.log")
    assert f"INFO :: {len(lacking)} contig(s) of the assembly are not in {metacool}: they get no contacts" in log.read_text()
    # The map of the clusters, in the output directory: its contigs, at 1 kb,
    # in their clusters.
    written = mcool / metacool.name
    metacooler.validate(str(written), full=True)
    _read_by_cooler_and_metacooler(written)
    mc = metacooler.MetaCooler(str(written))
    assert mc.binsize == 1000 and set(mc.contigs.name) == kept
    final_bin = final.set_index("Name").Final_cluster.replace("ND", "__unbinned__")
    assert dict(zip(mc.contigs.name, mc.contigs.genome_id.map(mc.genomes.name), strict=True)) == final_bin.to_dict()
    assert cooler.Cooler(str(written)).pixels()[:]["count"].sum() == contacts
    # What metaTOR measured on the contigs of the map, which
    # contig_data_final.txt lists: metacooler annotate read the network's
    # table, as the history of the map says.
    columns = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    assert [column for column in mc.contigs.columns if not column.startswith("qc_")] == ["name", "length", "genome_id", *columns]
    measured = final.set_index("Name").loc[mc.contigs.name, columns]
    for column in columns:
        assert mc.contigs[column].tolist() == measured[column].tolist(), column
    assert "Contigs table" not in log.read_text()
    assert mc.metadata["annotate"][-1]["contigs_table"].endswith("contig_data_network.txt")
    assert not [path.name for path in mcool.iterdir() if path.name.startswith(".")]


def test_pipeline_from_pairs_or_bam_files(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """From the pairs files or the BAM files of the libraries of a run (two
    here), the pipeline makes the map of the contigs, one bin per contig, in
    its output directory, named after --prefix, in a step more, first; then
    gives the clusters, the tables and the map of the clusters of a run from
    the map of those pairs, as metacooler create pairs makes it: the same bins,
    contacts, contigs, clusters and contact statistics."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-p", "sample", "-T", str(tmp_path / "tmp")]  # fmt: skip
    contig_data = mtn.create_contig_data(str(metagenome_inputs["assembly"]), str(metagenome_inputs["depth"]), "DpnII")
    reference = contig_map(metagenome_inputs["pairs"], contig_data, tmp_path / "sample.metacool", columns=False)
    random.seed(0)
    run_cli("pipeline", str(reference), "-o", str(tmp_path / "map"), *common)
    expected = tmp_path / "map" / "sample.metacool"

    # (tmp_path / "inputs" holds the copies of the files of the community.)
    inputs = tmp_path / "libraries"
    inputs.mkdir()
    lengths = {contig.name: contig.length for contig in metagenome.contigs}
    pairs = metagenome.pairs
    half = len(pairs) // 2
    lines = Path(metagenome_inputs["pairs"]).read_text().splitlines(keepends=True)
    head, body = [line for line in lines if line.startswith("#")], [line for line in lines if not line.startswith("#")]
    for name, rows in [("L1.pairs", body[:half]), ("L2.pairs", body[half:])]:
        (inputs / name).write_text("".join(head + rows))
    alignments.both_mates(inputs / "L1.bam", pairs[:half], lengths)
    alignments.both_mates(inputs / "L2.bam", pairs[half:], lengths, order="coordinate")
    alignments.one_mate([inputs / "L1_R1.bam", inputs / "L1_R2.bam"], pairs[:half], lengths)
    alignments.one_mate([inputs / "L2_R1.bam", inputs / "L2_R2.bam"], pairs[half:], lengths)
    runs = {
        "pairs": ([str(inputs / "L*.pairs")], "2 pairs files"),
        "both": ([str(inputs / "L1.bam"), str(inputs / "L2.bam")], "2 BAM files of both mates"),
        "one": ([str(inputs / "L*_R1.bam"), str(inputs / "L*_R2.bam")], "4 BAM files of one mate each"),
    }
    for name, (arguments, what) in runs.items():
        out = tmp_path / name
        random.seed(0)
        run_cli("pipeline", *arguments, "-o", str(out), *common)
        _assert_same_clusters(out, tmp_path / "map")
        written = out / "sample.metacool"
        metacooler.validate(str(written), full=True)
        clr, ref = cooler.Cooler(str(written)), cooler.Cooler(str(expected))
        for table in ["bins", "pixels", "chroms"]:
            pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
        mc, ref_mc = metacooler.MetaCooler(str(written)), metacooler.MetaCooler(str(expected))
        pd.testing.assert_frame_equal(mc.genomes, ref_mc.genomes)
        assert mc.cluster_rounds == ref_mc.cluster_rounds
        pd.testing.assert_frame_equal(mc.clusters, ref_mc.clusters)
        assert clr.info["metadata"]["qc"] == ref.info["metadata"]["qc"]
        assert clr.info["metadata"]["commands"][0].startswith("metator pipeline ")
        parameters = clr.info["metadata"]["metator"]["validation"]["parameters"]
        assert parameters["input"] == arguments and parameters["map"] == str(written)
        # The output directory: the map, the log and the outputs of the run
        # (the FASTA files of its clusters and of its MAGs, their cool files,
        # the tables, the plots and the report).
        (log,) = out.glob("metator_pipeline_*.log")
        outputs = ["cool", "fasta", log.name, "plots", "sample.metacool", "tables", *mtr.REPORT_FILES]
        assert sorted(path.name for path in out.iterdir()) == sorted(outputs)
        steps = [line.split(" :: ", 2)[2] for line in log.read_text().splitlines() if " :: INFO :: [" in line]
        assert steps == [
            f"[1/5] map: the contact map of the contigs, from {what}",
            f"[2/5] network: the contigs of {written} and their contacts",
            "[3/5] partition: the clusters of the network",
            "[4/5] validation: miComplete checks the clusters, and the contaminated ones are partitioned again",
            "[5/5] export: the outputs of the run, from its map",
        ]
        # Without --mges, the run has no steps mge and host, and says so.
        assert "INFO :: No list of MGE contigs (--mges): the run has no steps mge and host, " in log.read_text()


def _clustering(path):
    """What a map says of its clustering: the contigs of each cluster of each of
    its rounds, and the genome of each contig."""
    mc = metacooler.MetaCooler(str(path))
    members = {}
    for entry in mc.cluster_rounds:
        table = mc.cluster_members(entry["name"])
        members[entry["name"]] = sorted(zip(table["cluster"].astype(str), table["contig"].astype(str), strict=True))
    genomes = dict(zip(mc.contigs["name"].astype(str), mc.contigs["genome_id"].map(mc.genomes["name"]), strict=True))
    return members, genomes


def test_a_seed_gives_the_same_clusters(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """--seed seeds the random numbers that Louvain draws (Python's random
    module, which igraph draws from) where a run partitions: from the same map,
    the same seed gives the same clusters, every round of them, whatever the
    random state before the run; and partition --seed S then validation --seed S,
    step by step, give the clusters of pipeline --seed S, which seeds its
    partition and its recursive partition alike. The seed is stored with the
    other parameters of the run."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.first_cluster_contaminated(metagenome))
    assembly, depth = str(metagenome_inputs["assembly"]), str(metagenome_inputs["depth"])
    contigs = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    tmp = ["-T", str(tmp_path / "tmp")]
    maps = {}
    for name, state in [("first", 1), ("second", 2)]:
        out = tmp_path / name
        random.seed(state)
        run_cli("pipeline", str(contigs), "-a", assembly, "-d", depth, "-e", "DpnII", "-s", "20000", "-i", "10", "-j", "5",
                "--seed", "13", "-o", str(out), *tmp)  # fmt: skip
        maps[name] = out / "contigs.metacool"
    steps = tmp_path / "steps"
    random.seed(3)
    run_cli("network", str(contigs), "-a", assembly, "-d", depth, "-e", "DpnII", "-o", str(steps), *tmp)
    metacool = steps / "contigs.metacool"
    random.seed(4)
    run_cli("partition", str(metacool), "-i", "10", "--seed", "13", "-o", str(steps), *tmp)
    random.seed(5)
    run_cli("validation", str(metacool), "-a", assembly, "-s", "20000", "-i", "5", "--seed", "13", "-o", str(steps), *tmp)
    maps["steps"] = metacool
    found = {name: _clustering(path) for name, path in maps.items()}
    rounds = [entry for entry in found["first"][0] if entry.startswith("recursive_")]
    assert rounds, "no recursive partition: the test checks none of its seed"
    assert found["first"] == found["second"] == found["steps"]
    for name, path in maps.items():
        entries = metacooler.MetaCooler(str(path)).metadata["metator"]
        assert entries["partition"]["parameters"]["seed"] == 13 and entries["validation"]["parameters"]["seed"] == 13, name


def _contig_map(metagenome_inputs, path, pairs=None):
    """The map of the contigs of the mock community at path, in no cluster yet
    (see tests.maps.contig_map): one bin per contig, what metaTOR measures on
    the assembly (with a depth file and an enzyme) and the contact statistics
    of the pairs (those of the community, or pairs)."""
    contig_data = mtn.create_contig_data(str(metagenome_inputs["assembly"]), str(metagenome_inputs["depth"]), "DpnII")
    return contig_map(pairs or metagenome_inputs["pairs"], contig_data, path, assembly_name="synthetic")


def _annotated(objects):
    """The history of the metacooler annotate edits of a map (see
    tests.hdf5.walk), the temporary table that the network was stored from
    masked (the contigs_table of the entries that name one in the temporary
    directory of the run, whose name is drawn at random)."""
    (metadata,) = [json.loads(value) for key, _, value in objects["/"][2] if key == "metadata"]
    entries = metadata.get("annotate", [])
    for entry in entries:
        if str(entry.get("contigs_table")).endswith("/contig_data_network.txt"):
            entry["contigs_table"] = "<the temporary contig data of the network>"
    return entries


def _without_run(objects):
    """The objects of a map (see tests.hdf5.walk) but what differs between two
    runs of the same clusters: the date the map was written (creation-date, an
    attribute of the file), the command line and the parameters of the run that
    wrote it and of the runs that stored its network, its partition and its
    MGE clusters, and found their hosts (the command and parameters of the
    entries of the metator key of its metadata, a JSON text: network,
    partition, validation, and mge and host after them; and of the entry of the
    round mge, which records the run that made it, in the attribute rounds of
    the /clusters group, a JSON text), and the history of its metacooler
    annotate edits, which records each store of the network (see _annotated):
    the steps store it twice, with network and validation, pipeline once."""
    kind, members, attributes = objects["/"]
    kept = []
    for key, dtype, value in attributes:
        if key == "metadata":
            metadata = json.loads(value)
            for entry in metadata["metator"].values():
                del entry["command"], entry["parameters"]
            del metadata["annotate"]
            value = metadata
        if key != "creation-date":
            kept.append((key, dtype, value))
    objects = {**objects, "/": (kind, members, kept)}
    if "/clusters" in objects:
        kind, members, attributes = objects["/clusters"]
        rounds = [(key, dtype, json.loads(value) if key == "rounds" else value) for key, dtype, value in attributes]
        for key, _, value in rounds:
            for entry in value if key == "rounds" else []:
                for run in ["command", "parameters"]:
                    entry.pop(run, None)
        objects["/clusters"] = (kind, members, rounds)
    return objects


def _metator_metadata(objects):
    """The metator key of the metadata of a map (see tests.hdf5.walk)."""
    (metadata,) = [value for key, _, value in objects["/"][2] if key == "metadata"]
    return json.loads(metadata)["metator"]


def test_a_cluster_partitioned_again_but_kept_whole(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """A contaminated cluster that the recursive step partitions again, into
    clusters too incomplete to replace it, is kept whole: its contigs keep the
    IDs of their recursive clusters in contig_data_final.txt, and are in it, in
    binning.txt (which left them out) and in the map of the clusters as in the
    round final, its FASTA file listing them."""
    monkeypatch.setattr(
        "metator.validation.micomplete_quality", synthetic.first_cluster_contaminated(metagenome, kept_whole=True)
    )
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    out = tmp_path / "out"
    random.seed(0)
    run_cli("pipeline", metacool, "-a", metagenome_inputs["assembly"], "-o", out, "-T", tmp_path / "tmp", "-s", "20000",
            "-i", "10", "-j", "5")  # fmt: skip
    final = pd.read_csv(_kept(out) / "contig_data_final.txt", sep="\t", dtype=str, keep_default_na=False)
    kept_whole = final[(final.Recursive_cluster_ID != "00000") & (final.Final_cluster != "ND")]
    assert len(kept_whole) and set(kept_whole.Final_cluster) == {"metator_00001_00000"}
    binning = pd.read_csv(_kept(out) / "binning.txt", sep="\t", header=None, names=["Name", "Bin"], dtype=str)
    assert list(binning.Name) == list(final.Name[final.Final_cluster != "ND"])
    mc = metacooler.MetaCooler(str(out / metacool.name))
    assert "recursive_1" in [entry["name"] for entry in mc.cluster_rounds]
    members = mc.cluster_members("final")
    in_cluster = set(members.contig[members.cluster == "metator_00001_00000"])
    assert set(kept_whole.Name) <= in_cluster
    listed = (out / MAGS_WITH_MGES / "metator_00001_00000.fa").read_text().splitlines()
    assert in_cluster == {line[1:].split()[0] for line in listed if line.startswith(">")}


def test_steps_on_a_map_give_the_map_of_pipeline(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """The steps run one by one on a map (network MAP, partition MAP, validation
    MAP) give the map that pipeline MAP gives, from copies of the same map of
    the contigs, and the same tables: every group, soft link, dataset (its
    dtype, shape, chunks, compression and values) and attribute of the HDF5
    file, the network that both store included (the columns of the contigs
    table and the weights of the bins), but the date it was written
    (creation-date), the command line and parameters of the runs that wrote it
    (the command and parameters of the metator key of its metadata and of its
    entries network and partition: pipeline's, network's, partition's and
    validation's), and the history of its metacooler annotate edits (see
    _without_run), in which the steps store the network twice, network MAP then
    validation MAP, which stores it again after its regroup, and pipeline once."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.first_cluster_contaminated(metagenome))
    assembly = str(metagenome_inputs["assembly"])
    contig_map = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    runs = {name: tmp_path / name for name in ["pipeline", "steps"]}
    for directory in runs.values():
        directory.mkdir()
        shutil.copyfile(contig_map, directory / "input.metacool")
    # In place (--inplace): both runs update their input.metacool, whose
    # history of reassignments names it.
    common = ["-a", assembly, "-o", "out", "-T", "tmp", "-s", "20000", "--inplace"]
    # partition reads no assembly, and keeps every cluster.
    partition_args = ["-o", "out", "-T", "tmp", "--inplace"]

    monkeypatch.chdir(runs["pipeline"])
    random.seed(0)
    run_cli("pipeline", "input.metacool", *common, "-d", metagenome_inputs["depth"], "-e", "DpnII", "-i", "10", "-j", "5")

    monkeypatch.chdir(runs["steps"])
    metacool = runs["steps"] / "input.metacool"
    random.seed(0)
    run_cli("network", "input.metacool", "-a", assembly, "-o", "out", "-T", "tmp", "-d", metagenome_inputs["depth"],
            "-e", "DpnII", "--inplace")  # fmt: skip
    metacooler.validate(str(metacool), full=True)
    # The network, and nothing else.
    before, after = walk(contig_map), walk(metacool)
    changed = {name for name in before.keys() | after.keys() if before.get(name) != after.get(name)}
    columns = {f"/contigs/{column}" for column in ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]}
    assert {"/", "/bins", "/bins/metator_empirical_hit", "/contigs", "/contigs/Hit"} <= changed
    assert changed <= {"/", "/bins", "/bins/metator_empirical_hit", "/contigs", *columns}
    mc = metacooler.MetaCooler(str(metacool))
    assert mc.metadata["metator"]["network"]["parameters"]["normalization"] == "empirical_hit"
    assert mc.cluster_rounds == [] and set(mc.genomes.name) == {"__unbinned__"}
    after_network = after
    run_cli("partition", "input.metacool", *partition_args, "-i", "10")
    metacooler.validate(str(metacool), full=True)
    mc = metacooler.MetaCooler(str(metacool))
    assert mc.cluster_rounds == [
        {"name": name, "tool": "metaTOR", "version": __version__, "step": "partition"} for name in ["core", "overlapping"]
    ]
    # Every cluster, named as the validation reads it back: the core clusters
    # core_<number>, the overlapping ones metator_<number>_00000.
    clusters = mtmc.read_partition(str(metacool))
    assert list(mc.clusters.name) == [mtcl.core_name(number) for number in clusters.core] + list(clusters.names.values())
    assert not (runs["steps"] / "out" / "overlapping_clusters").exists()
    assert not (runs["steps"] / "out" / "contig_data_partition.txt").exists()
    # genome, which metacooler computes: no cluster is a genome of the map yet.
    assert list(mc.clusters.columns) == ["round", "name", "parent_round", "parent", "n_contigs", "length", "genome"]
    assert mc.clusters.genome.isna().all()
    # Every contig of the network in one core and one overlapping cluster.
    members = mc.cluster_members()
    assert members.groupby("round").contig.apply(lambda contigs: contigs.is_unique).all()
    assert set(members.contig[members["round"] == "core"]) == set(members.contig[members["round"] == "overlapping"])
    partition = mc.metadata["metator"]["partition"]
    assert partition["normalization"] == "empirical_hit" and partition["command"].startswith("metator partition input.metacool")
    # Nothing else of the map changed: its genomes and contigs tables, its
    # contacts, its network, the rest of its metadata.
    before, after = after_network, walk(metacool)
    kind, stored, attributes = after.pop("/")
    before_kind, before_stored, before_attributes = before.pop("/")
    assert (kind, [name for name in stored if name != "clusters"]) == (before_kind, before_stored)
    metadata = {key: (json.loads(value) if key == "metadata" else value) for key, _, value in attributes}
    before_metadata = {key: (json.loads(value) if key == "metadata" else value) for key, _, value in before_attributes}
    del metadata["metadata"]["metator"]["partition"]
    assert metadata == before_metadata
    assert {name: value for name, value in after.items() if not name.startswith("/clusters")} == before
    run_cli("validation", "input.metacool", *common, "-i", "5")

    pipeline, steps = runs["pipeline"] / "out", runs["steps"] / "out"
    for name in ["cluster_summary.txt", "contig_data_final.txt", "binning.txt", "clusters.txt", "cluster_members.txt"]:
        assert (_kept(steps) / name).read_bytes() == (_kept(pipeline) / name).read_bytes(), name
    written = {name: directory / "input.metacool" for name, directory in runs.items()}
    for path in written.values():
        metacooler.validate(str(path), full=True)
    found = {name: walk(path) for name, path in written.items()}
    assert "/bins/metator_empirical_hit" in found["steps"]
    assert _without_run(found["steps"]) == _without_run(found["pipeline"])
    # The history of the edits: the validation stores again the network that
    # network MAP stored, once it regrouped the map, which put the network
    # before the QC columns (see metator.metacool.write_map); pipeline
    # stores it once, after its regroup. The network comes last in both.
    network = [f"contigs.{column}" for column in ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]]
    annotated = {name: _annotated(objects) for name, objects in found.items()}
    assert annotated["steps"][:-2] == annotated["pipeline"]
    assert [(entry["dropped"], entry["added"]) for entry in annotated["steps"][-2:]] == [(network, []), ([], network)]
    assert annotated["pipeline"][-1]["added"] == network
    assert list(metacooler.MetaCooler(str(written["steps"])).contigs.columns)[-4:] == [name[8:] for name in network]
    # What differs: the command lines and their parameters, those of
    # pipeline, and of network and validation.
    metator = {name: _metator_metadata(objects) for name, objects in found.items()}
    for step in ["network", "partition", "validation"]:
        assert metator["pipeline"][step]["command"].startswith("metator pipeline input.metacool")
        assert metator["steps"][step]["command"].startswith(f"metator {step} input.metacool")
    validation = {name: entries["validation"] for name, entries in metator.items()}
    assert validation["steps"]["parameters"]["map"] == validation["pipeline"]["parameters"]["map"] == "input.metacool"
    assert "start" not in validation["steps"]["parameters"]
    assert validation["steps"]["summary_columns"] == validation["pipeline"]["summary_columns"]
    # Both partitioned the network that both maps store, which the
    # reassignment of the map of the clusters left as it was: a validation
    # reads it.
    digests = {name: metator[name]["partition"][mtmc.NETWORK_DIGEST] for name in runs}
    assert digests["steps"] == digests["pipeline"] == mtn.stored_network(str(written["pipeline"])).digest
    # The rounds of the run, the final clusters those of the genomes of the
    # map.
    mc = metacooler.MetaCooler(str(written["steps"]))
    assert [entry["name"] for entry in mc.cluster_rounds] == ["core", "overlapping", "recursive_1", "final"]
    assert set(mc.cluster_members("final").cluster) == set(mc.genomes.name) - {"__unbinned__"}
    assert not list(tmp_path.rglob("*.lock")) and not list(tmp_path.rglob(".input.metacool.*"))
    # The tables that metator dump writes from the two maps are the same.
    for name, path in written.items():
        run_cli("dump", path, "-o", tmp_path / f"dump_{name}")
    for table in [
        "contig_data_final.txt",
        "cluster_summary.txt",
        "binning.txt",
        "clusters.txt",
        "cluster_members.txt",
        "network.txt",
    ]:
        assert (tmp_path / "dump_steps" / table).read_bytes() == (tmp_path / "dump_pipeline" / table).read_bytes(), table


def _members(mc):
    """The contigs of the clusters of a map, by name, in a fixed order."""
    return mc.cluster_members().sort_values(["round", "cluster", "contig"]).reset_index(drop=True)


def _read_by_cooler_and_metacooler(path):
    """cooler and metacooler read every table of a map written in place: its
    pixels and bins through cooler dump, every pixel through metacooler check
    --deep (their command lines, as a user runs them)."""
    for table in ("chroms", "bins", "pixels"):
        subprocess.run(["cooler", "dump", "-t", table, str(path)], check=True, capture_output=True)
    subprocess.run(["metacooler", "check", "--deep", str(path)], check=True, capture_output=True)


# The lines of the log of the recursive steps of a validation, for each way
# they end (see metator.validation.recursive_decontamination): the step, of
# RECURSIVE_STEPS at most, what it partitions again, what replaces it, and when
# no further step is needed.
STEP_TIME = r"\(0h 0m [0-9]+s\)"
RECURSIVE_LOGS = {
    "no cluster contaminated": (
        "fake",
        0,
        [r"Recursive step 1/9: no cluster is contaminated, so no recursive step is needed\."],
    ),
    "contaminated twice": (
        "twice",
        1,
        [
            r"Recursive step 1/9: 1 contaminated cluster\(s\) to partition again\.",
            rf"Recursive step 1/9: 1 of them replaced by the [0-9]+ cluster\(s\) of more than 20000 bp that they split into {STEP_TIME}\.",
            r"Recursive step 2/9: 1 contaminated cluster\(s\) to partition again\.",
            rf"Recursive step 2/9: 1 of them replaced by the [0-9]+ cluster\(s\) of more than 20000 bp that they split into {STEP_TIME}\.",
            r"Recursive step 3/9: no contaminated cluster left, so no further recursive step is needed\.",
        ],
    ),
    "kept whole": (
        "kept whole",
        1,
        [
            r"Recursive step 1/9: 1 contaminated cluster\(s\) to partition again\.",
            rf"Recursive step 1/9: none of them is replaced, as the [0-9]+ cluster\(s\) they split into lose too much of their completeness {STEP_TIME}, so no further recursive step is needed\.",
        ],
    ),
    "no part larger than --size": (
        "twice",
        1,
        [
            r"Recursive step 1/9: 1 contaminated cluster\(s\) to partition again\.",
            r"Recursive step 1/9: none of them splits into clusters of more than 20000 bp, so no further recursive step is needed\.",
        ],
    ),
    "the last step": (
        "twice",
        1,
        [
            r"Recursive step 1/2: 1 contaminated cluster\(s\) to partition again\.",
            rf"Recursive step 1/2: 1 of them replaced by the [0-9]+ cluster\(s\) of more than 20000 bp that they split into {STEP_TIME}\.",
            r"Recursive step 2/2: 1 contaminated cluster\(s\) to partition again\.",
            rf"Recursive step 2/2: 1 of them replaced by the [0-9]+ cluster\(s\) of more than 20000 bp that they split into {STEP_TIME}\.",
            r"Recursive step 2/2 was the last: its clusters are kept as they are\.",
        ],
    ),
}


@pytest.mark.parametrize("case", list(RECURSIVE_LOGS))
def test_the_log_follows_the_steps_of_the_run(case, metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli, logged):
    """The log of a pipeline run marks its four steps ([1/4] network... [4/4]
    export), says
    how many clusters miComplete found contaminated, then follows each
    recursive step of the validation (Recursive step 1/9...), what it
    partitions again and replaces, and says when no further step is needed:
    no contaminated cluster, none left, parts too incomplete to replace their
    cluster, or none larger than --size (here a partition into none, as the
    stub recursive clustering gives it); or that the last step allowed was
    the last (here 2). It ends with the map that the run wrote, and how long
    it took."""
    stub = {
        "fake": synthetic.fake_micomplete(metagenome),
        "twice": synthetic.contaminated_twice(metagenome),
        "kept whole": synthetic.first_cluster_contaminated(metagenome, kept_whole=True),
    }[RECURSIVE_LOGS[case][0]]
    monkeypatch.setattr("metator.validation.micomplete_quality", stub)
    if case == "no part larger than --size":
        clustering = mtv.recursive_clustering

        def into_none(*args, **kwargs):
            _, contigs_data, _, _ = clustering(*args, **kwargs)
            return False, contigs_data, {}, {}

        monkeypatch.setattr(mtv, "recursive_clustering", into_none)
    if case == "the last step":
        monkeypatch.setattr(mtv, "RECURSIVE_STEPS", 2)
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    out = tmp_path / "out"
    random.seed(0)
    run_cli("pipeline", metacool, "-a", metagenome_inputs["assembly"], "-o", out, "-T", tmp_path / "tmp", "-s", "20000",
            "-i", "10", "-j", "5")  # fmt: skip
    messages = [message for level, message in logged() if level == "INFO"]
    stages = [message for message in messages if message.startswith("[")]
    assert [stage[:5] for stage in stages] == ["[1/4]", "[2/4]", "[3/4]", "[4/4]"]
    assert messages.index(stages[0]) < messages.index(stages[1]) < messages.index(stages[2])
    _, contaminated, expected = RECURSIVE_LOGS[case]
    checked = [message for message in messages if message.startswith("miComplete checked")]
    assert len(checked) == 1 and re.fullmatch(
        rf"miComplete checked the 5 clusters in 0h 0m [0-9]+s, and found {contaminated} of them contaminated, which the "
        r"recursive steps partition again\.",
        checked[0],
    )
    steps = [message for message in messages if message.startswith("Recursive step")]
    assert len(steps) == len(expected) and all(re.fullmatch(pattern, step) for pattern, step in zip(expected, steps, strict=True))
    assert re.fullmatch(rf"Final contact map written after 0h 0m [0-9]+s: {re.escape(str(out / metacool.name))}", messages[-1])


# The MGE steps of the pipeline (--mges): metator mge --reassign and metator
# host on the map of its clusters.


def _mge_list(metagenome, path):
    """The list of the MGE contigs of the mock community at path, one name per
    line, as metator mge and pipeline --mges read it."""
    truth = metagenome.truth
    path.write_text("".join(f"{name}\n" for name in truth.Name[truth.Kind == "mge"]))
    return path


def test_pipeline_with_mges_is_the_steps_one_by_one(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """pipeline MAP --mges gives the map that the steps give one by one, network
    MAP, partition MAP, validation MAP, then mge MAP --reassign and host MAP,
    from copies of the same map of the contigs (see
    test_steps_on_a_map_give_the_map_of_pipeline: but the command lines of the
    runs), and, with export MAP -a ASSEMBLY last, the same outputs (see
    _assert_same_outputs): its steps mge and host, 4 and 5 here, are mge
    --reassign and host, on the map of its clusters, and its step export, 6,
    is export. Its log marks its six steps. Each MGE cluster is a genome of
    the map, and its host is the final cluster of the genome that hosts it.
    The outputs: fasta/mags, the MAGs without their MGE contigs, which
    fasta/mags_with-mges keeps; cool, a cool file per MAG with the MGE
    clusters that it hosts, balanced; and tables."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    assembly = str(metagenome_inputs["assembly"])
    contig_map = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    mges = _mge_list(metagenome, tmp_path / "mges.txt")
    runs = {name: tmp_path / name for name in ["pipeline", "steps"]}
    for directory in runs.values():
        directory.mkdir()
        shutil.copyfile(contig_map, directory / "input.metacool")
    common = ["-a", assembly, "-o", "out", "-T", "tmp", "-s", "20000", "--inplace"]

    monkeypatch.chdir(runs["pipeline"])
    random.seed(0)
    run_cli("pipeline", "input.metacool", *common, "-i", "10", "-j", "5", "--mges", mges)

    monkeypatch.chdir(runs["steps"])
    random.seed(0)
    run_cli("network", "input.metacool", "-a", assembly, "-o", "out", "-T", "tmp", "--inplace")
    run_cli("partition", "input.metacool", "-o", "out", "-T", "tmp", "--inplace", "-i", "10")
    run_cli("validation", "input.metacool", *common, "-i", "5")
    run_cli("mge", "input.metacool", "-a", assembly, "-m", mges, "--reassign", "-o", "out", "-T", "tmp", "--inplace")
    run_cli("host", "input.metacool", "-o", "out", "--inplace")
    run_cli("export", "input.metacool", "-a", assembly, "-o", "out", "-T", "tmp", "-F")

    pipeline, steps = runs["pipeline"] / "out", runs["steps"] / "out"
    _assert_same_outputs(steps, pipeline)
    for directory in [MAGS_WITH_MGES, os.path.join("fasta", "mges")]:
        files = [{path.name: path.read_bytes() for path in (out / directory).iterdir()} for out in (steps, pipeline)]
        assert files[0] == files[1] and files[1], directory
    assert sorted(path.name for path in pipeline.iterdir() if path.is_dir()) == ["cool", "fasta", "plots", "tables"]
    assert sorted(path.name for path in (pipeline / "tables").iterdir()) == sorted([*mtd.TABLES, "mge_hosts.txt"])
    written = {name: directory / "input.metacool" for name, directory in runs.items()}
    for path in written.values():
        metacooler.validate(str(path), full=True)
    found = {name: walk(path) for name, path in written.items()}
    assert _without_run(found["steps"]) == _without_run(found["pipeline"])
    # What differs: the command lines of the runs, and their parameters.
    metator = {name: _metator_metadata(objects) for name, objects in found.items()}
    assert set(metator["pipeline"]) == set(metator["steps"]) == {"network", "partition", "validation", "mge", "host"}
    for step in ["mge", "host"]:
        assert metator["pipeline"][step]["command"].startswith("metator pipeline input.metacool")
        assert metator["steps"][step]["command"].startswith(f"metator {step} input.metacool")
    # The six steps of the run, in its log, and no word of a stop.
    (log,) = pipeline.glob("metator_pipeline_*.log")
    logged = log.read_text()
    names = ["network", "partition", "validation", "mge", "host", "export"]
    assert re.findall(r" :: INFO :: \[([0-9])/([0-9])\] ([a-z]+): ", logged) == [
        (str(number), "6", name) for number, name in enumerate(names, start=1)
    ]
    assert "No list of MGE contigs" not in logged and "The run stopped" not in logged
    # Each MGE cluster, a genome of the map, and its host, the final cluster
    # of the genome that hosts it.
    mc = metacooler.MetaCooler(str(written["pipeline"]))
    members = mc.cluster_members(mtcl.MGE)
    mge_of = dict(zip(members["contig"], members["cluster"], strict=True))
    genomes = mc.genomes
    assert sorted(genomes.name[genomes.mge == "True"]) == sorted(set(mge_of.values()))
    final = mc.cluster_members(mtcl.FINAL)
    final_of = dict(zip(final["contig"], final["cluster"], strict=True))
    clusters = mc.clusters[(mc.clusters["round"] == mtcl.MGE).to_numpy()].set_index("name")
    for mge in metagenome.mges:
        (cluster,) = {mge_of[contig.name] for contig in mge.contigs}
        (host,) = {final_of[contig.name] for contig in metagenome.molecule(mge.host).contigs}
        assert (clusters.loc[cluster, "n_hosts"], clusters.loc[cluster, "hosts"]) == (1, host)
        # The MGE contigs: in the FASTA file of their host's final cluster
        # (the validation clusters them with it), not in that of its MAG.
        contigs = {contig.name for contig in mge.contigs}
        assert contigs <= set(synthetic.read_fasta(pipeline / MAGS_WITH_MGES / f"{host}.fa"))
        assert not contigs & set(synthetic.read_fasta(pipeline / "fasta" / "mags" / f"{host}.fa"))
    # fasta/mags: each MAG, a genome of the map but the MGE genomes, its
    # contigs and stitched; cool: a cool file per MAG, its contigs and those of
    # the MGE clusters that it hosts, balanced (weights in its bins).
    genome_of = dict(zip(mc.contigs["name"], genomes["name"].to_numpy()[mc.contigs["genome_id"].to_numpy()], strict=True))
    mags = list(genomes.name[genomes.mge == "False"])
    assert sorted(path.name for path in (pipeline / "fasta" / "mags").iterdir()) == sorted(
        f"{mag}{suffix}" for mag in mags for suffix in (".fa", ".stitched.fa")
    )
    assert sorted(path.name for path in (pipeline / "cool").iterdir()) == sorted(f"{mag}.cool" for mag in mags)
    hosted = {
        name: mges.split(",") for name, mges in zip(genomes.name, genomes.mges, strict=True) if isinstance(mges, str) and mges
    }
    assert sorted(hosted) == sorted(clusters.hosts[clusters.n_hosts > 0])
    for mag in mags:
        contigs = sorted(name for name, genome in genome_of.items() if genome == mag)
        assert sorted(synthetic.read_fasta(pipeline / "fasta" / "mags" / f"{mag}.fa")) == contigs
        assert list(synthetic.read_fasta(pipeline / "fasta" / "mags" / f"{mag}.stitched.fa")) == [mag]
        clr = cooler.Cooler(str(pipeline / "cool" / f"{mag}.cool"))
        mge_contigs = [contig for contig, cluster in mge_of.items() if cluster in hosted.get(mag, [])]
        assert sorted(clr.chromnames) == sorted(contigs + mge_contigs) and "weight" in clr.bins().columns


def test_pipeline_with_mges_again_on_the_map_it_wrote(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli, logged):
    """pipeline MAP --mges, run again on the map it updated in place, warns before
    any work that it replaces the clusters, the rounds, the round mge among
    them, and the hosts of the MGE clusters that the map holds; but not, once
    it has written the map of its clusters, that the map no longer holds the
    round mge and the hosts, which its steps mge and host store again."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    mges = _mge_list(metagenome, tmp_path / "mges.txt")
    args = [metacool, "-a", metagenome_inputs["assembly"], "-o", tmp_path / "out", "-T", tmp_path / "tmp", "-s", "20000",
            "-i", "10", "-j", "5", "--mges", mges, "--inplace"]  # fmt: skip
    random.seed(0)
    run_cli("pipeline", *args)
    first = len(logged())
    random.seed(0)
    run_cli("pipeline", *args, "-F")
    # The cool files of the smallest MAGs, which cooler balance cannot
    # balance on this map of one bin per contig, aside (see
    # metator.outputs.write_cool_files).
    warnings = [message for level, message in logged()[first:] if level == "WARNING" and not message.startswith("Every bin of ")]
    replaced = (
        rf"{re.escape(str(metacool))} holds [0-9]+ cluster\(s\) and the clustering round\(s\) core, overlapping, "
        r"(recursive_[0-9]+, )*final, mge and the hosts of its MGE clusters \(metator host\), which this run replaces, "
        r"in place: leave out --inplace to write the map of this run to another file, and leave .* as it is\."
    )
    assert len(warnings) == 1 and re.fullmatch(replaced, warnings[0]), warnings
    mc = metacooler.MetaCooler(str(metacool))
    assert [entry["name"] for entry in mc.cluster_rounds][-1] == mtcl.MGE
    assert set(mc.metadata["metator"]) == {"network", "partition", "validation", "mge", "host"}
    assert (mc.genomes.mge == "True").any() and mc.genomes.n_mges.notna().any()


@pytest.mark.parametrize("step", [4, 5], ids=["mge", "host"])
def test_pipeline_stopped_in_an_mge_step_keeps_its_clusters(step, metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """A pipeline run with --mges that stops in its step 4, mge, leaves the map of
    its clusters that step 3, the validation, wrote, and one that stops in its
    step 5, host, the map with its MGE clusters too, which mge wrote; the log
    says so, and what to run on the map to complete the run."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))

    def fails(*args, **kwargs):
        raise mtmc.MetacoolError("cannot")

    monkeypatch.setattr("metator.mge.mge_binning" if step == 4 else "metator.host.find_hosts", fails)
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    mges = _mge_list(metagenome, tmp_path / "mges.txt")
    out = tmp_path / "out"
    random.seed(0)
    with pytest.raises(mtmc.MetacoolError, match=r"^cannot$"):
        run_cli("pipeline", metacool, "-a", metagenome_inputs["assembly"], "-o", out, "-T", tmp_path / "tmp", "-s",
                "20000", "-i", "10", "-j", "5", "--mges", mges)  # fmt: skip
    written = out / metacool.name
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    held = ["network", "partition", "validation", *(["mge"] if step == 5 else [])]
    assert set(mc.metadata["metator"]) == set(held)
    assert (mtcl.MGE in [entry["name"] for entry in mc.cluster_rounds]) == (step == 5)
    # fasta/mges, made anew, holds the MGE clusters of mge; the hosts and
    # the outputs of the run (fasta/mags, cool, tables) are not written.
    assert bool(list((out / "fasta" / "mges").iterdir())) == (step == 5)
    assert not (out / "tables").exists() and not list((out / "fasta" / "mags").iterdir())
    assert not list((out / "cool").iterdir())
    (log,) = out.glob("metator_pipeline_*.log")
    name, what, finish = {
        4: ("mge", "the clusters of the run", "metator mge --reassign then metator host"),
        5: ("host", "the clusters of the run and their MGE clusters", "metator host"),
    }[step]
    logged = log.read_text()
    assert "ERROR :: cannot\n" in logged
    stopped = (
        f"INFO :: The run stopped in its step {step}, {name}: {written} holds {what}; run {finish} on it to complete "
        "the map, then metator export to write the outputs of the run from it.\n"
    )
    assert stopped in logged


def test_pipeline_with_mges_from_pairs_files(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """pipeline --mges from the pairs files of a run makes the map of the contigs
    first, in a step more, then gives what it gives from that map, as
    metacooler create pairs makes it: the same clusters, MGE clusters and
    hosts, and the same files (the FASTA files of the clusters, of the MAGs and
    of the MGE clusters, and the hosts). Its log marks its six steps: map, the
    three of the binning, then mge and host, 5 and 6 here."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    assembly, depth = str(metagenome_inputs["assembly"]), str(metagenome_inputs["depth"])
    mges = _mge_list(metagenome, tmp_path / "mges.txt")
    common = ["-a", assembly, "-d", depth, "-e", "DpnII", "-s", "20000", "-i", "10", "-j", "5", "-p", "sample",
              "-T", str(tmp_path / "tmp"), "--mges", str(mges)]  # fmt: skip
    contig_data = mtn.create_contig_data(assembly, depth, "DpnII")
    reference = contig_map(metagenome_inputs["pairs"], contig_data, tmp_path / "sample.metacool", columns=False)
    runs = {"map": str(reference), "pairs": str(metagenome_inputs["pairs"])}
    for name, argument in runs.items():
        random.seed(0)
        run_cli("pipeline", argument, "-o", str(tmp_path / name), *common)
    maps = {name: metacooler.MetaCooler(str(tmp_path / name / "sample.metacool")) for name in runs}
    pd.testing.assert_frame_equal(maps["pairs"].genomes, maps["map"].genomes)
    pd.testing.assert_frame_equal(maps["pairs"].clusters, maps["map"].clusters)

    def rounds(mc):
        # The round mge records the run that made it, whose input differs.
        return [
            {key: value for key, value in entry.items() if key not in ("command", "parameters")} for entry in mc.cluster_rounds
        ]

    assert rounds(maps["pairs"]) == rounds(maps["map"]) and rounds(maps["map"])[-1]["name"] == mtcl.MGE
    genomes = maps["map"].genomes
    assert (genomes.mge == "True").any() and genomes.n_mges.notna().any()
    for directory in [MAGS_WITH_MGES, os.path.join("fasta", "mags"), os.path.join("fasta", "mges")]:
        files = [{path.name: path.read_bytes() for path in (tmp_path / name / directory).iterdir()} for name in runs]
        assert files[0] == files[1] and files[0], directory
    hosts = os.path.join("tables", "mge_hosts.txt")
    assert (tmp_path / "pairs" / hosts).read_bytes() == (tmp_path / "map" / hosts).read_bytes()
    (log,) = (tmp_path / "pairs").glob("metator_pipeline_*.log")
    names = ["map", "network", "partition", "validation", "mge", "host", "export"]
    assert re.findall(r" :: INFO :: \[([0-9])/([0-9])\] ([a-z]+): ", log.read_text()) == [
        (str(number), "7", name) for number, name in enumerate(names, start=1)
    ]


def _output_files(directory):
    """The outputs of a run in its output directory, by path relative to it: its
    files but its logs and its map."""
    return {
        str(path.relative_to(directory)): path
        for path in sorted(Path(directory).rglob("*"))
        if path.is_file() and not path.name.endswith((".log", ".metacool"))
    }


def _records(path):
    """The records of a FASTA file, as (header, sequence), in no order."""
    blocks = [block.split("\n", 1) for block in Path(path).read_text().split(">")[1:]]
    return sorted((header, sequence.replace("\n", "")) for header, sequence in blocks)


def _assert_same_outputs(out, ref, without=None):
    """The outputs of two runs (out, ref: their output directories; but those of
    ref under without): the same files, of the same bytes, but the FASTA files
    of the final clusters, the same records (metator export writes them in the
    order of the map), the plots, PDF files (which record when they were
    drawn), the cool files, the same bins (weights included) and pixels, and
    the report, the same run (see _assert_same_report)."""
    files = _output_files(out)
    expected = {name: path for name, path in _output_files(ref).items() if without is None or not name.startswith(without)}
    assert sorted(files) == sorted(expected)
    for name, path in files.items():
        if name.endswith(".cool"):
            clr, ref_clr = cooler.Cooler(str(path)), cooler.Cooler(str(expected[name]))
            pd.testing.assert_frame_equal(clr.bins()[:], ref_clr.bins()[:])
            pd.testing.assert_frame_equal(clr.pixels()[:], ref_clr.pixels()[:])
        elif name.endswith(".pdf"):
            assert path.read_bytes().startswith(b"%PDF"), name
        elif name.startswith(MAGS_WITH_MGES):
            assert _records(path) == _records(expected[name]), name
        elif name in mtr.REPORT_FILES:
            _assert_same_report(name, path, expected[name], without)
        else:
            assert path.read_bytes() == expected[name].read_bytes(), name


def _report_page(path, without=None):
    """What report.html shows of a run: its figures, and the outputs of the run
    that it links to (but those under without)."""
    page = Path(path).read_text()
    (figures,) = re.findall(r'<script id="metator-figures" type="application/json">(.*?)</script>', page, re.S)
    links = [link for link in re.findall(r'href="([^"#][^"]*)"', page) if not link.startswith("https://")]
    return json.loads(figures), sorted(link for link in links if without is None or not link.startswith(without))


def _assert_same_report(name, path, ref, without=None):
    """A file of the report of two runs (path, and ref, its file in the other
    run): the same run, but when each was written, its map, and the command
    lines of its steps (report.txt: its sections report and commands; and
    report.html, which shows them): report.txt, the same values in its other
    sections, and the same tables; report_mqc.yaml, the same; report.html, the
    same figures, and the same links to the outputs of the run (those of ref
    but under without)."""
    if name == mtr.REPORT_TXT:
        (sections, tables), (ref_sections, ref_tables) = mtr.read_txt(path), mtr.read_txt(ref)
        assert list(sections) == list(ref_sections) and list(tables) == list(ref_tables)
        for section, values in sections.items():
            if section == "report":
                volatile = {"map", "written", "map_written"}
                values = {key: value for key, value in values.items() if key not in volatile}
                assert values == {key: value for key, value in ref_sections[section].items() if key not in volatile}
            elif section != "commands":
                assert values == ref_sections[section], section
        for table, rows in tables.items():
            pd.testing.assert_frame_equal(rows, ref_tables[table], obj=table)
    elif name == mtr.REPORT_MQC:
        assert yaml.safe_load(path.read_text()) == yaml.safe_load(ref.read_text())
    else:
        assert _report_page(path) == _report_page(ref, without)


def test_export_writes_the_outputs_of_a_run_again(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli, logged):
    """metator export MAP -a ASSEMBLY writes the outputs of a run again from the
    map that pipeline --mges wrote, which it leaves as it was: the same files
    (see _assert_same_outputs), the hosts found again with the thresholds of
    the run. Without --assembly, the same but the FASTA files, which the log
    says how to write. In an output directory that holds them, refused without
    -F, before anything is written, and replaced with it. A map written while
    export reads it is refused, at its end."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    assembly = str(metagenome_inputs["assembly"])
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    mges = _mge_list(metagenome, tmp_path / "mges.txt")
    out, tmp = tmp_path / "out", tmp_path / "tmp"
    random.seed(0)
    run_cli("pipeline", metacool, "-a", assembly, "-o", out, "-T", tmp, "-s", "20000", "-i", "10", "-j", "5", "--mges", mges)
    written = out / metacool.name
    before = written.read_bytes()

    exported = tmp_path / "exported"
    run_cli("export", written, "-a", assembly, "-o", exported, "-T", tmp)
    _assert_same_outputs(exported, out)
    assert _output_files(exported) and written.read_bytes() == before
    (log,) = exported.glob("metator_export_*.log")
    assert (
        "INFO :: The hosts of the MGE clusters, found again from the network of the map with --min-share 10 and "
        f"--min-contigs 5, as the host run recorded them: {exported / 'tables' / 'mge_hosts.txt'}, "
        f"{exported / 'plots' / 'mge_hosts.pdf'}\n"
    ) in log.read_text()

    partial = tmp_path / "partial"
    run_cli("export", written, "-o", partial, "-T", tmp)
    _assert_same_outputs(partial, out, without="fasta")
    (log,) = partial.glob("metator_export_*.log")
    assert (
        "INFO :: No assembly (--assembly): the FASTA files of the run (fasta) are not written. Give the assembly that "
        f"the map was made from to write them: {assembly}, as the run recorded it.\n"
    ) in log.read_text()

    # Again, in the same output directory: refused without -F, then replaced.
    (exported / "cool" / "mine.txt").write_text("an earlier run's\n")
    with pytest.raises(ValueError):
        run_cli("export", written, "-a", assembly, "-o", exported, "-T", tmp)
    refused = f"{exported / MAGS_WITH_MGES} already existed. Remove directory or use -F argument to overwrite it."
    assert ("ERROR", refused) in logged()
    assert (exported / "cool" / "mine.txt").exists()
    run_cli("export", written, "-a", assembly, "-o", exported, "-T", tmp, "-F")
    _assert_same_outputs(exported, out)

    # An assembly that lacks a contig of a final cluster: the run stops as it
    # writes the FASTA files, and says so.
    missing = str(metacooler.MetaCooler(str(written)).cluster_members(mtcl.FINAL)["contig"].iloc[0])
    lacking = tmp_path / "lacking.fa"
    records = Path(assembly).read_text().split(">")[1:]
    lacking.write_text("".join(f">{record}" for record in records if record.split()[0] != missing))
    with pytest.raises(mtmc.MetacoolError, match=missing):
        run_cli("export", written, "-a", lacking, "-o", tmp_path / "lacking", "-T", tmp)
    stopped = (
        f"The run stopped as it wrote the outputs of {written}, which is complete: metator export writes them from it "
        "(with -F, in place of what this run wrote)."
    )
    assert ("INFO", stopped) in logged()

    # The map, written as export reads it (here, as it writes its tables).
    write_map_outputs = mtcli.write_map_outputs

    def writes_the_map(metacool, *args):
        write_map_outputs(metacool, *args)
        status = os.stat(metacool)
        os.utime(metacool, ns=(status.st_atime_ns, status.st_mtime_ns + 10**9))

    monkeypatch.setattr("metator.cli.write_map_outputs", writes_the_map)
    with pytest.raises(mtmc.MapReplacedError):
        run_cli("export", written, "-o", tmp_path / "replaced", "-T", tmp)
    message = (
        f"{written} was replaced or written while metator export read it, by a run that updated it: what it wrote may "
        "mix two maps. Run metator export again, with -F for what it wrote, once the run that writes the map has ended."
    )
    assert ("ERROR", message) in logged()


@pytest.mark.parametrize("case", ["no clusters", "the assembly in cool", "a tmpdir in cool"])
def test_export_refusals(case, metagenome_inputs, tmp_path, run_cli):
    """What metator export refuses before it writes anything: a map without the
    clusters of a validation, as metator dump refuses it; an input of the run
    in a folder that it makes anew, which -F would replace with it (with -F
    here); and a --tmpdir there, as a usage error."""
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    out = tmp_path / "out"
    assembly = metagenome_inputs["assembly"]
    if case == "a tmpdir in cool":
        result = run_cli("export", metacool, "-o", out, "-T", out / "cool" / "tmp", exit_code=2)
        assert f"is in {out / 'cool'}, which the run makes anew" in " ".join(result.output.replace("│", " ").split())
        assert not out.exists()
        return
    kept = []
    if case == "the assembly in cool":
        (out / "cool").mkdir(parents=True)
        assembly = shutil.copyfile(assembly, out / "cool" / "assembly.fa")
        kept = ["cool"]
        message = (
            f"{assembly}, an input of the run, is in {out / 'cool'}, which the run makes anew: move it away, or give "
            "the run another output directory with -o."
        )
    else:
        message = (
            f"{metacool} holds no clusters of a validation: it lacks the validation entry of its metator metadata and "
            "the clustering round final. Run metator validation (or metator pipeline) on it first, which writes the "
            "map of the clusters."
        )
    with pytest.raises(ValueError):
        run_cli("export", metacool, "-a", assembly, "-o", out, "-T", tmp_path / "tmp", "-F")
    (log,) = out.glob("metator_export_*.log")
    assert f"ERROR :: {message}\n" in log.read_text()
    assert sorted(path.name for path in out.iterdir()) == sorted([log.name, *kept]) and Path(assembly).exists()


def test_the_map_holds_the_values_of_micomplete(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli, logged):
    """The map of the clusters holds the summary of the validation at full
    precision (cluster_summary_full.txt, beside cluster_summary.txt, whose
    floats have 2 decimals): a cluster of a weighted completeness of 0.8951
    is an MQ MAG in the log of the run and in what export reads from the map
    for its plots, not an HQ MAG of 0.90; metator dump writes
    cluster_summary.txt as the run wrote it, its floats rounded."""

    def just_below(fasta_dir, outfile, threads):
        synthetic.write_micomplete_summary(fasta_dir, outfile, metagenome)
        # Three comment lines, the header, then a row per cluster.
        with open(outfile) as summary:
            lines = [line.split("\t") for line in summary.read().splitlines()]
        lines[4][lines[3].index("Weighted completeness")] = "0.8951"
        with open(outfile, "w") as summary:
            summary.write("".join("\t".join(line) + "\n" for line in lines))

    monkeypatch.setattr("metator.validation.micomplete_quality", just_below)
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    out = tmp_path / "out"
    random.seed(0)
    run_cli("pipeline", metacool, "-a", metagenome_inputs["assembly"], "-o", out, "-T", tmp_path / "tmp", "-s", "20000",
            "-i", "10", "-j", "5", "-N")  # fmt: skip
    messages = [message for level, message in logged() if level == "INFO"]
    assert any(message.startswith("HQ MAGs: 4\t") for message in messages)
    assert any(message.startswith("MQ MAGs: 1\t") for message in messages)
    mc = metacooler.MetaCooler(str(out / metacool.name))
    genomes = mc.genomes.set_index("name")
    (name,) = genomes.index[genomes["Weighted completeness"] == 0.8951]
    assert mtv.bin_quality(0.8951, genomes.at[name, "Weighted redundancy"]) == "MQ"
    rows = mc.metadata["metator"]["validation"]["summary_rows"]
    assert mto.final_summary(mc, rows)[name]["Weighted completeness"] == 0.8951
    # cluster_summary.txt, as the run wrote it and as dump writes it: 0.90.
    tables = _kept(out)
    summary = pd.read_csv(tables / "cluster_summary.txt", sep="\t", index_col=0, dtype=str)
    assert summary.at[name, "Weighted completeness"] == "0.90"
    assert (out / "tables" / "cluster_summary.txt").read_bytes() == (tables / "cluster_summary.txt").read_bytes()
    full = pd.read_csv(tables / "cluster_summary_full.txt", sep="\t", index_col=0)
    assert full.at[name, "Weighted completeness"] == 0.8951
