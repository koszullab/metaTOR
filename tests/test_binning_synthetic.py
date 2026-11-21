# End-to-end binning tests on the mock metagenome of tests/synthetic.py,
# checked against its ground truth: network, partition, validation (with a
# stub miComplete) and the whole pipeline.

import json
import os
import random
import shlex
import shutil
import subprocess
from pathlib import Path

import cooler
import metacooler
import numpy as np
import pandas as pd
import pytest

import metator.cli as mtcli
import metator.clusters as mtcl
import metator.metacool as mtmc
import metator.network as mtn
import metator.partition as mtp
import metator.validation as mtv
from metator.version import __version__
from tests import synthetic
from tests.hdf5 import walk
from tests.maps import contig_map, pairs_map

NORMALIZATIONS = ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"]


# The tables of a validation (contig_data_final.txt, cluster_summary.txt,
# binning.txt, clusters.txt, cluster_members.txt) are in the temporary
# directory of the run, which -N keeps: every run of this module keeps it,
# and the tests read them there (see _kept).
TABLES = ["contig_data_final.txt", "cluster_summary.txt", "binning.txt", "clusters.txt", "cluster_members.txt"]


@pytest.fixture(autouse=True)
def _keep_the_tables(monkeypatch):
    """Keep the temporary directory of every run, as -N does."""
    private_tmp_dir = mtcli.private_tmp_dir
    monkeypatch.setattr(mtcli, "private_tmp_dir", lambda tmpdir, clean_up=True: private_tmp_dir(tmpdir, clean_up=False))


def _kept(out):
    """The directory of the tables of the latest validation or pipeline run
    whose output directory is out: the temporary directory of the run (in
    the -T of its command line, which its log records; a relative one from
    the directory of out, where the tests run) that holds them, the last
    written before the end of its log (two runs may share a -T)."""
    logs = sorted(Path(out).glob("metator_*.log"), key=os.path.getmtime)
    for log in reversed(logs):
        line = next(line for line in log.read_text().splitlines() if line.startswith("## Command used: "))
        command = shlex.split(line[len("## Command used: ") :])
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
    clusters at path (cluster_summary.txt), as the map alone writes it
    again: its genomes table holds the values, and its entry validation how
    they are written (its rows, columns and their formats, see
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
    """Network edges with the names of their contigs (by their IDs in the
    contig data)."""
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
    """Sizes, GC contents and shotgun coverages of the contig data, which
    metator network stores in the map: its contigs, the longest first."""
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
    """The overlapping clusters of the partition are the genomes, and those
    that reach 50 kb, which the validation writes for miComplete to check
    (see metator.validation.write_overlapping_fasta), one each."""
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
    """A copy, at path, of the map of the mock community in which metator
    network stored its network (see _network_map), with a partition stored
    as metator partition stores it (see metator.metacool.write_rounds):
    the genomes ``merged`` share a cluster while every other genome has its
    own, MGE contigs going with their host, each cluster both a core and an
    overlapping cluster."""
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
    """A copy of the map of the mock community in which metator network
    stored its network, at path: validation MAP and partition MAP read it,
    and write the map of their clusters over it."""
    return Path(shutil.copyfile(metagenome_network["map"], path))


def test_validation_splits_a_merged_bin(metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli):
    """A bin holding two genomes is reported contaminated by (stub) miComplete:
    the recursive decontamination splits it back into the two genomes."""
    random.seed(0)
    calls = []
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome, calls))
    metacool = _map_with_merged_partition(
        metagenome, metagenome_network, tmp_path / "network.metacool", merged=("genome_1", "genome_2")
    )
    out = tmp_path / "out"
    # Outputs of an earlier run, replaced thanks to --force.
    for stale in ["final_clusters"]:
        (out / stale).mkdir(parents=True)
        (out / stale / "metator_00009_00000.fa").write_text(">old\nACGT\n")
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
    assert not (out / "final_clusters" / "metator_00009_00000.fa").exists()

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
    assert sorted(p.name for p in (out / "final_clusters").iterdir()) == sorted(f"{name}.fa" for name in summary.index)
    binning = pd.read_csv(_kept(out) / "binning.txt", sep="\t", header=None, names=["Name", "Bin"])
    assert set(binning.Bin) == set(summary.index)
    assert sorted(binning.Name) == sorted(final.Name[final.Final_cluster != "ND"])
    # No plots: metaTOR draws none since 2.0.
    assert not (out / "plot").exists()
    # The bins of the split are recursive bins, of step 1; the others are not.
    assert list(summary.recursive) == list(summary.step > 0)
    assert summary.recursive.any() and not summary.recursive.all()
    # The map alone writes the summary again, byte for byte.
    assert _summary_from_the_map(metacool) == (_kept(out) / "cluster_summary.txt").read_text()
    # metacooler wrote the FASTA files of the clusters of the recursive step,
    # and indexed nothing beside the assembly.
    assert not os.path.exists(f"{metagenome_inputs['assembly']}.fxi")


def test_validation_decompresses_the_assembly_once(
    metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli
):
    """An assembly compressed with gzip is decompressed once for the run, in
    its temporary directory: every FASTA file of the clusters (the
    overlapping clusters, each recursive step) is written from that copy,
    where metacooler copied it again for each write. The final clusters are
    those of the run from the plain assembly."""
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
        finals[name] = {path.name: path.read_bytes() for path in (tmp_path / name / "final_clusters").iterdir()}
    # The overlapping clusters, then the clusters of the recursive step.
    assert len(written) == 2 and len(set(written)) == 1
    (plain,) = set(written)
    assert plain.name == "assembly.fa" and plain.parent.parent == tmp_path / "tmp_gz"
    assert finals["gz"] == finals["plain"]


def test_validation_checks_the_clusters_in_the_order_of_the_partition(
    metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli
):
    """The FASTA file of each overlapping cluster that miComplete checks lists
    its contigs as metator partition lists them (see
    metator.partition.partition): by core cluster, in the order of their
    numbers, then by ID, as the files of overlapping_clusters/ did, not in
    the order of the map. Here, an overlapping cluster of two core clusters,
    the first of which holds the later contigs of the map."""
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
    which it grouped the contigs by cluster: the partition stored in the map
    is still that of its network, which the regroup left as it was (see
    metator.network.network_digest), the same contig for contig, and the
    run reads it with the IDs of the contigs in the regrouped map, their
    positions in it, which number the rows of its tables."""
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
    common = ["-a", metagenome_inputs["assembly"], "-T", tmp_path / "tmp", "-s", "50000", "-i", "10", "-t", "1"]
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
    """Contigs named by numbers (000017, 1e5, 2.50), which pandas reads as
    numbers unless told otherwise: a run reads every table of names as
    text, and each contig goes by its name, as the map names it, to the
    final cluster of the FASTA file that lists it (read as numbers, the
    names lost their leading zeros, or their form, and the clusters their
    contigs)."""
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
    files = sorted((out / "final_clusters").glob("*.fa"))
    assert sorted(path.stem for path in files) == sorted(set(binned.Final_cluster))
    for path in files:
        listed = [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]
        assert sorted(listed) == sorted(name for name, genome in genome_of.items() if genome == path.stem)


def test_a_contig_named_na(tmp_path, monkeypatch, run_cli):
    """A contig named NA, which pandas reads as a missing value unless told
    otherwise, is a contig as any other: validation reads it by name from
    the partition table of the run, and writes it into the map of the
    clusters, in its final cluster, with its size (a missing name made the
    Size column of contig_data_final.txt a float, which the map refused;
    the run keeps its tables in its temporary directory: see
    _keep_the_tables)."""
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


def test_validation_without_bin_above_the_size_threshold(metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli):
    """partition MAP stores every overlapping cluster, whatever its size; a
    validation run of which none reaches --size stops before miComplete,
    saying why, and leaves its map as it was."""
    calls = []
    monkeypatch.setattr("metator.validation.micomplete_quality", lambda *args: calls.append(args))
    metacool = _network_map(metagenome_network, tmp_path / "network.metacool")
    run_cli("partition", str(metacool), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), "-i", "2")
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
    """No bin of the partition reaches --size: the run stops before
    miComplete (which found no bin, then failed to read its empty results),
    saying why in its log file. It leaves its map as it was, in no cluster."""
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
    # The contact map, updated in place: one bin per contig, one genome per
    # cluster (the cluster summary columns with it), and every pair a
    # contact.
    metacooler.validate(str(metacool))
    clr = cooler.Cooler(str(metacool))
    assert clr.pixels()[:]["count"].sum() == len(metagenome.pairs)
    assert clr.info["metadata"]["metator"]["validation"]["parameters"]["map"] == str(metacool)
    mc = metacooler.MetaCooler(str(metacool))
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
    # directory, kept here as -N keeps it: see _keep_the_tables).
    assert not os.path.exists(f"{metagenome_inputs['assembly']}.fxi")
    assert sorted(path.name for path in out.iterdir() if not path.name.endswith(".log")) == ["final_clusters"]
    assert not (out / "contig_data_partition.txt").exists()


def test_pipeline_map_is_the_map_of_the_pairs(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """The map that a run leaves is the map of the pairs of its input map
    and its final tables (see tests.maps.pairs_map): the same bins,
    contacts, contigs table (column for column), genomes table (the columns
    of the cluster summary and the contact statistics of each cluster) and
    QC block. The run reads no pairs file: it writes the map of its
    clusters over its map at the end with metacooler reassign. Only the
    history of the map differs, the annotate and reassign keys of its
    metadata (and the date of the file)."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    from_pairs = metacooler.from_pairs
    maps = []

    def counted_from_pairs(uri, *args, **kwargs):
        maps.append(str(uri))
        return from_pairs(uri, *args, **kwargs)

    written = _contig_map(metagenome_inputs, tmp_path / "inputs" / "contigs.metacool")
    monkeypatch.setattr(metacooler, "from_pairs", counted_from_pairs)
    out = tmp_path / "out"
    args = [str(written), "-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]),
            "-e", "DpnII", "-s", "50000", "-i", "10", "-j", "5", "-o", str(out), "-T", str(tmp_path / "tmp")]  # fmt: skip
    random.seed(0)
    run_cli("pipeline", *args)
    assert maps == []
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
    assert [entry["source"] for entry in metadata["reassign"]] == [str(written)]
    # The columns of the input map (tests.maps.contig_map), then those that
    # the run stores, once, as metator network stores them (dropped, then
    # added).
    assert len(metadata["annotate"]) == 3 and len(ref_metadata["annotate"]) == 1


def _assert_normalized_as_the_network(written, tmp_path):
    """The contacts of the map of a run normalized by the weights it stores
    and summed by pair of contigs (metacooler network --balance) are the
    edges of its network (see metator.network.map_network: but those of
    weight 0), whatever the order of its contigs."""
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
    """A run that stops with an error (miComplete fails in the validation,
    here) leaves its map as it was, in no cluster. A new run from it (with
    -F in the same output directory, whose cluster directories the stopped
    run left) gives the outputs of a run that does not stop, and its map,
    but for the provenance of the map."""
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-T", str(tmp_path / "tmp")]  # fmt: skip
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
    final_bins = [{path.name: path.read_bytes() for path in (run / "final_clusters").iterdir()} for run in [out, clean]]
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
    """The outputs of two runs (out, ref: their output directories) of the
    same clusters from maps whose contigs are in another order, or that
    lack some of them (contigs: the contigs of the map of out, those of
    ref by default): the IDs of the contigs are their positions in the map
    of the run, and so are the rows of its tables, and the sums over the
    contigs of a cluster (its abundances) are made in that order, so equal
    to the last bits. The same clusters, of the same contigs, under the same
    names, with the same values; the same final FASTA files, but for the
    order of their records. With contigs, the shotgun abundance of the
    clusters (SG_abundance), a share of the sequenced bases of the contigs
    of the map, is not compared: those that the map lacks are not
    counted."""

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
    fasta = [{path.name: synthetic.read_fasta(path) for path in (run / "final_clusters").iterdir()} for run in (out, ref)]
    assert fasta[0] == fasta[1] and fasta[0]


def test_pipeline_from_a_metacool(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """From the contact map of the clusters of a first run (one bin per
    contig, its contigs grouped by cluster), the pipeline makes the same
    network and the same clusters, but for the order of the contigs, which
    is that of the map (see _assert_same_clusters). It updates the input map
    in place (a copy of the map of the first run): the map of the first
    run, QC and clustering rounds included, but for its provenance."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-T", str(tmp_path / "tmp")]  # fmt: skip
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
    """From a map of the contigs whose contigs table holds what metaTOR
    measures with a depth file and an enzyme (see tests.maps.contig_map), a
    run without --depth and --enzyme reads the restriction sites and the
    coverage of the contigs from its Restriction_site and Shotgun_coverage
    columns (theoritical_hit needs both): the same network, tables and
    clusters as with them, byte for byte, and the map of its clusters
    stores them back."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    common = ["-a", str(metagenome_inputs["assembly"]), "-m", "theoritical_hit", "-s", "50000", "-i", "10", "-j", "5",
              "-T", str(tmp_path / "tmp")]  # fmt: skip
    options = ["-d", str(metagenome_inputs["depth"]), "-e", "DpnII"]
    metacool = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    columns = metacooler.MetaCooler(str(metacool)).contigs[["name", "Restriction_site", "Shotgun_coverage"]]
    runs = {"options": options, "map": []}
    for name, given in runs.items():
        random.seed(0)
        # Each run writes its own map (--output-map): a run from a map
        # updates it in place otherwise, and the second would start from
        # the first's.
        run_cli("pipeline", str(metacool), "-o", str(tmp_path / name),
                "--output-map", str(tmp_path / name / "metator.metacool"), *common, *given)  # fmt: skip

    by_options, by_map = tmp_path / "options", tmp_path / "map"
    for name in ["contig_data_final.txt", "cluster_summary.txt", "binning.txt"]:
        assert (_kept(by_map) / name).read_bytes() == (_kept(by_options) / name).read_bytes(), name
    # The same network, stored in each map.
    networks = [mtn.map_network(str(run / "metator.metacool")) for run in (by_options, by_map)]
    pd.testing.assert_frame_equal(networks[0], networks[1])
    assert "SG_abundance" in pd.read_csv(_kept(by_map) / "cluster_summary.txt", sep="\t").columns
    (log,) = by_map.glob("metator_pipeline_*.log")
    assert f"INFO :: The shotgun coverage of the contigs: the Shotgun_coverage column of {metacool}." in log.read_text()
    # Stored back: the map of the bins holds the values that the run read.
    for written in [by_map / "metator.metacool", by_options / "metator.metacool"]:
        metacooler.validate(str(written), full=True)
        stored = metacooler.MetaCooler(str(written)).contigs.set_index("name")
        pd.testing.assert_frame_equal(stored.loc[columns.name, columns.columns[1:]], columns.set_index("name"))


def test_pipeline_from_a_filtered_metacool(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """From a map that lacks contigs of the assembly (made with a contig
    filter, which leaves out the pairs with a read on the contigs it drops),
    the pipeline completes and makes the clusters of a run from the map of
    the pairs that the map holds, every contig in it: the contigs it lacks
    have no hits, and are in no cluster and in no table of the run (see
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
    run_cli("pipeline", str(kept_map), "-o", str(tmp_path / "pair"), "--output-map", str(tmp_path / "pair.metacool"), *common)
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
    # The map, updated in place: its contigs, at 1 kb, in their clusters.
    written = metacool
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


def _contig_map(metagenome_inputs, path, pairs=None):
    """The map of the contigs of the mock community at path, in no cluster
    yet (see tests.maps.contig_map): one bin per contig, what metaTOR
    measures on the assembly (with a depth file and an enzyme) and the
    contact statistics of the pairs (those of the community, or pairs)."""
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
    runs of the same clusters: the date the map was written
    (creation-date, an attribute of the file), the command line and the
    parameters of the run that wrote it and of the runs that stored its
    network and its partition (the command and parameters of the metator
    key of its metadata, a JSON text, and of its entries network and
    partition), and the history of its metacooler annotate edits, which
    records each store of the network (see _annotated): the steps store it
    twice, with network and validation, pipeline once."""
    kind, members, attributes = objects["/"]
    kept = []
    for key, dtype, value in attributes:
        if key == "metadata":
            metadata = json.loads(value)
            for run in ["command", "parameters"]:
                for step in ["network", "partition", "validation"]:
                    del metadata["metator"][step][run]
            del metadata["annotate"]
            value = metadata
        if key != "creation-date":
            kept.append((key, dtype, value))
    return {**objects, "/": (kind, members, kept)}


def _metator_metadata(objects):
    """The metator key of the metadata of a map (see tests.hdf5.walk)."""
    (metadata,) = [value for key, _, value in objects["/"][2] if key == "metadata"]
    return json.loads(metadata)["metator"]


def _first_cluster_contaminated(metagenome, kept_whole=False):
    """A stub miComplete (see synthetic.fake_micomplete) that finds the first
    overlapping cluster contaminated (a weighted redundancy of 1.5): the
    validation partitions it again, which makes the round recursive_1, whose
    clusters replace it; with kept_whole, they are too incomplete for that
    (a weighted completeness of 0.1), and the cluster is kept whole."""

    def micomplete_quality(fasta_dir, outfile, threads):
        synthetic.write_micomplete_summary(fasta_dir, outfile, metagenome)
        # Three comment lines, the header, then a row per cluster.
        lines = [line.split("\t") for line in Path(outfile).read_text().splitlines()]
        if Path(fasta_dir).name == "overlapping_clusters":
            lines[4][lines[3].index("Weighted redundancy")] = "1.5"
        elif kept_whole:
            for line in lines[4:]:
                line[lines[3].index("Weighted completeness")] = "0.1"
        Path(outfile).write_text("".join("\t".join(line) + "\n" for line in lines))

    return micomplete_quality


def test_a_cluster_partitioned_again_but_kept_whole(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """A contaminated cluster that the recursive step partitions again, into
    clusters too incomplete to replace it, is kept whole: its contigs keep
    the IDs of their recursive clusters in contig_data_final.txt, and are in
    it, in binning.txt (which left them out) and in the map of the clusters
    as in the round final, its FASTA file listing them."""
    monkeypatch.setattr("metator.validation.micomplete_quality", _first_cluster_contaminated(metagenome, kept_whole=True))
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
    mc = metacooler.MetaCooler(str(metacool))
    assert "recursive_1" in [entry["name"] for entry in mc.cluster_rounds]
    members = mc.cluster_members("final")
    in_cluster = set(members.contig[members.cluster == "metator_00001_00000"])
    assert set(kept_whole.Name) <= in_cluster
    listed = (out / "final_clusters" / "metator_00001_00000.fa").read_text().splitlines()
    assert in_cluster == {line[1:].split()[0] for line in listed if line.startswith(">")}


def test_steps_on_a_map_give_the_map_of_pipeline(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """The steps run one by one on a map (network MAP, partition MAP,
    validation MAP) give the map that pipeline MAP gives, from copies of the
    same map of the contigs, and the same tables: every group, soft link,
    dataset (its dtype, shape, chunks, compression and values) and attribute
    of the HDF5 file, the network that both store included (the columns of
    the contigs table and the weights of the bins), but the date it was
    written (creation-date), the command line and parameters of the runs
    that wrote it (the command and parameters of the metator key of its
    metadata and of its entries network and partition: pipeline's,
    network's, partition's and validation's), and the history of its
    metacooler annotate edits (see _without_run), in which the steps store
    the network twice, network MAP then validation MAP, which stores it
    again after its regroup, and pipeline once. Each run is made in a
    directory of its own, with the same relative paths, which the history of
    the map records (the input map of metacooler reassign, the tables of
    metacooler annotate). The rounds are
    core, overlapping, recursive_1 (a cluster that the stub miComplete finds
    contaminated, partitioned again) and final.

    Run in one process after one random.seed(0), as igraph draws from
    Python's random module: the steps draw the stream that pipeline draws
    (the temporary directory of each command drew from it too, so a run
    step by step, with three, drew another one than pipeline, with one).

    Each step leaves a valid map: network MAP stores the network of the
    contigs (the GC_content, Hit, Shotgun_coverage and Restriction_site
    columns of the contigs table, the weights metator_empirical_hit of the
    bins and the entry network of the metadata), and changes nothing else;
    partition MAP stores the partition, every core and overlapping cluster
    with its contigs, in its clustering rounds (the rounds core and
    overlapping, without miComplete's columns, which the validation adds to
    the overlapping clusters that it checks), and the entry partition of
    the metadata, and changes nothing else; validation MAP writes the map
    of the clusters, as pipeline MAP does at its end."""
    monkeypatch.setattr("metator.validation.micomplete_quality", _first_cluster_contaminated(metagenome))
    assembly = str(metagenome_inputs["assembly"])
    contig_map = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    runs = {name: tmp_path / name for name in ["pipeline", "steps"]}
    for directory in runs.values():
        directory.mkdir()
        shutil.copyfile(contig_map, directory / "input.metacool")
    common = ["-a", assembly, "-o", "out", "-T", "tmp", "-s", "20000"]
    # partition reads no assembly, and keeps every cluster.
    partition_args = ["-o", "out", "-T", "tmp"]

    monkeypatch.chdir(runs["pipeline"])
    random.seed(0)
    run_cli("pipeline", "input.metacool", *common, "-d", metagenome_inputs["depth"], "-e", "DpnII", "-i", "10", "-j", "5")

    monkeypatch.chdir(runs["steps"])
    metacool = runs["steps"] / "input.metacool"
    random.seed(0)
    run_cli("network", "input.metacool", "-a", assembly, "-o", "out", "-T", "tmp", "-d", metagenome_inputs["depth"],
            "-e", "DpnII")  # fmt: skip
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
