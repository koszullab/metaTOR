# End-to-end binning tests on the mock metagenome of tests/synthetic.py,
# checked against its ground truth: network, partition, validation (with a
# stub miComplete) and the whole pipeline.

import os
import random
import shutil
from pathlib import Path

import pandas as pd
import pytest

import metator.commands as mtc
import metator.network as mtn
import metator.partition as mtp
from tests import synthetic

NORMALIZATIONS = ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"]


def _named_edges(network_file, contig_data_file):
    """Network edges with the names of their contigs."""
    edges = pd.read_csv(network_file, sep="\t", header=None, names=["id1", "id2", "weight"])
    names = pd.read_csv(contig_data_file, sep="\t").set_index("ID")["Name"]
    edges["name1"] = edges.id1.map(names)
    edges["name2"] = edges.id2.map(names)
    return edges


def _genome_labels(table, metagenome, column):
    """Truth and found labels of the genome contigs (MGE contigs left out)."""
    truth = metagenome.truth
    table = table[table.Name.map(truth.Kind) == "genome"]
    return list(table.Name.map(truth.Genome)), list(table[column].astype(str))


@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_network_edges_follow_genomes(metagenome, metagenome_files, tmp_path, normalization):
    """With every normalisation, same-genome edges outweigh cross-genome ones."""
    contig_data, hit_data = mtn.create_contig_data(str(metagenome_files["assembly"]), 1, str(metagenome_files["depth"]), "DpnII")
    network, contig_data_file = mtn.alignment_to_contacts(
        [str(metagenome_files["pairs"])],
        contig_data,
        0,
        hit_data,
        str(tmp_path),
        "network.txt",
        "contig_data_network.txt",
        str(tmp_path),
        1,
        normalization,
        False,
    )
    truth = metagenome.truth
    edges = _named_edges(network, contig_data_file)
    genomes = (edges.name1.map(truth.Kind) == "genome") & (edges.name2.map(truth.Kind) == "genome")
    same = edges.name1.map(truth.Genome) == edges.name2.map(truth.Genome)
    within, across = edges.weight[genomes & same], edges.weight[genomes & ~same]
    assert len(within) > 0 and len(across) > 0
    assert within.mean() > 10 * across.mean()
    assert within.sum() > 0.9 * (within.sum() + across.sum())
    # Every contig is a node of the network.
    assert set(edges.name1) | set(edges.name2) == set(truth.Name)


def test_contig_data_matches_the_truth(metagenome, metagenome_network):
    """Sizes, GC contents and shotgun coverages of the contig data table."""
    data = pd.read_csv(metagenome_network["contig_data"], sep="\t").set_index("Name")
    truth = metagenome.truth
    assert list(data.index) == list(truth.Name)
    assert (data.Size == truth.Length).all()
    assert data.GC_content.to_numpy() == pytest.approx(truth.GC.to_numpy(), abs=1e-6)
    assert data.Shotgun_coverage.to_numpy() == pytest.approx(10 * truth.Abundance.to_numpy(), rel=0.25)
    # Every read pair gives two hits.
    assert data.Hit.sum() == 2 * len(metagenome.pairs)


@pytest.mark.parametrize("algorithm", ["louvain", "leiden"])
def test_partition_recovers_the_genomes(metagenome, metagenome_files, metagenome_network, tmp_path, algorithm):
    random.seed(0)  # igraph draws from Python's random module.
    fasta_dir = tmp_path / "overlapping_bin"
    fasta_dir.mkdir()
    _, contig_data_file = mtp.partition(
        algorithm,
        str(metagenome_files["assembly"]),
        False,
        str(metagenome_network["contig_data"]),
        10,
        str(metagenome_network["network"]),
        str(tmp_path),
        str(fasta_dir),
        0.8,
        1.0,
        50_000,
        str(tmp_path),
        1,
        "metator",
    )
    table = pd.read_csv(contig_data_file, sep="\t")
    truth, found = _genome_labels(table, metagenome, "Overlapping_bin_ID")
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    # One fasta per genome, each holding one genome.
    bins = sorted(fasta_dir.glob("*.fa"))
    assert len(bins) == len(metagenome.genomes)
    for fasta in bins:
        contigs = synthetic.read_fasta(fasta)
        assert synthetic.bin_quality(list(contigs), metagenome)[1] == 1.0


def _partition_with_merged_bin(metagenome, contig_data_file, out_dir, merged):
    """A partition table and its bins where the genomes ``merged`` share a bin
    while every other genome has its own; MGE contigs go with their host.
    """
    truth = metagenome.truth
    bin_id = {}
    for genome in metagenome.genomes:
        if genome.name in merged[1:]:
            bin_id[genome.name] = bin_id[merged[0]]
        else:
            bin_id[genome.name] = len(set(bin_id.values())) + 1
    data = pd.read_csv(contig_data_file, sep="\t")
    ids = data.Name.map(truth.Host).map(bin_id)
    sizes = data.groupby(ids)["Size"].transform("sum")
    counts = data.groupby(ids)["Size"].transform("size")
    for level in ["Core", "Overlapping"]:
        data[f"{level}_bin_ID"] = [f"{i:05d}" for i in ids]
        data[f"{level}_bin_contigs"] = counts
        data[f"{level}_bin_size"] = sizes
    partition_file = out_dir / "contig_data_partition.txt"
    data.to_csv(partition_file, sep="\t", index=False)
    fasta_dir = out_dir / "overlapping_bin"
    fasta_dir.mkdir()
    for i in sorted(set(ids)):
        metagenome.write_fasta(fasta_dir / f"metator_{i:05d}_00000.fa", list(data.Name[ids == i]))
    return partition_file, fasta_dir


def test_validation_splits_a_merged_bin(metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch):
    """A bin holding two genomes is reported contaminated by (stub) miComplete:
    the recursive decontamination splits it back into the two genomes."""
    random.seed(0)
    calls = []
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome, calls))
    partition_file, fasta_dir = _partition_with_merged_bin(
        metagenome, metagenome_network["contig_data"], tmp_path, merged=("genome_1", "genome_2")
    )
    out = tmp_path / "out"
    args = [
        "-a", str(metagenome_inputs["assembly"]),
        "-c", str(partition_file),
        "-f", str(fasta_dir),
        "-n", str(metagenome_network["network"]),
        "-o", str(out),
        "-T", str(tmp_path / "tmp"),
        "-s", "50000",
        "-i", "10",
        "-t", "1",
    ]  # fmt: skip
    mtc.Validation(args, {}).execute()

    # miComplete ran on the partition bins, then on the bins of the first
    # recursive step only: the split bins are clean.
    assert [Path(c[0]).name for c in calls] == ["overlapping_bin", "step_1"]
    final = pd.read_csv(out / "contig_data_final.txt", sep="\t")
    truth, found = _genome_labels(final, metagenome, "Final_bin")
    assert "ND" not in found
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    assert synthetic.purity(truth, found) >= 0.9
    summary = pd.read_csv(out / "bin_summary.txt", sep="\t", index_col=0)
    assert len(summary) == len(metagenome.genomes)
    assert "metator_00001_00000" not in summary.index
    assert summary["Weighted redundancy"].max() == 1.0
    # Hi-C abundances are proportional to the hits of the bins.
    hits = final.groupby("Final_bin").Hit.sum()[summary.index]
    hic = summary.HiC_abundance.astype(float)
    assert (hic / hic.sum()).to_numpy() == pytest.approx((hits / hits.sum()).to_numpy(), rel=1e-3)
    # The final fasta files are the final bins.
    assert sorted(p.stem for p in (out / "final_bin").glob("*.fa")) == sorted(summary.index)
    binning = pd.read_csv(out / "binning.txt", sep="\t", header=None, names=["Name", "Bin"])
    assert set(binning.Bin) == set(summary.index)
    assert (out / "plot" / "bins_distribution.png").stat().st_size > 0


@pytest.mark.xfail(
    strict=True,
    raises=FileNotFoundError,
    reason="Validation removes <assembly>.fxi, which only exists if a contaminated bin was re-partitioned",
)
def test_validation_without_contaminated_bin(metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch):
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    partition_file, fasta_dir = _partition_with_merged_bin(metagenome, metagenome_network["contig_data"], tmp_path, merged=())
    out = tmp_path / "out"
    args = [
        "-a", str(metagenome_inputs["assembly"]),
        "-c", str(partition_file),
        "-f", str(fasta_dir),
        "-n", str(metagenome_network["network"]),
        "-o", str(out),
        "-T", str(tmp_path / "tmp"),
        "-s", "50000",
    ]  # fmt: skip
    mtc.Validation(args, {}).execute()
    assert len(pd.read_csv(out / "bin_summary.txt", sep="\t", index_col=0)) == len(metagenome.genomes)


@pytest.mark.xfail(
    strict=True,
    raises=FileNotFoundError,
    reason="Partition removes <assembly>.fxi, which only exists if a bin passed the size threshold",
)
def test_partition_without_bin_above_the_size_threshold(metagenome_inputs, metagenome_network, tmp_path):
    args = [
        "-a", str(metagenome_inputs["assembly"]),
        "-c", str(metagenome_network["contig_data"]),
        "-n", str(metagenome_network["network"]),
        "-o", str(tmp_path / "out"),
        "-T", str(tmp_path / "tmp"),
        "-i", "2",
        "-s", "10000000",
    ]  # fmt: skip
    mtc.Partition(args, {}).execute()
    assert list((tmp_path / "out" / "overlapping_bin").iterdir()) == []


def test_pipeline_recovers_the_genomes(metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch):
    """The whole pipeline from pairs, with stub miComplete and scaffolding."""
    random.seed(0)
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    pairs = tmp_path / "inputs" / metagenome_indexed_pairs.name
    shutil.copyfile(metagenome_indexed_pairs, pairs)
    shutil.copyfile(f"{metagenome_indexed_pairs}.px2", f"{pairs}.px2")
    out = tmp_path / "out"
    args = [
        "-1", str(pairs),
        "-a", str(metagenome_inputs["assembly"]),
        "-d", str(metagenome_inputs["depth"]),
        "-S", "pair",
        "-e", "DpnII",
        "-m", "abundance",
        "-s", "50000",
        "-i", "10",
        "-j", "5",
        "-o", str(out),
        "-T", str(tmp_path / "tmp"),
        "-C",
        "--scaffold",
    ]  # fmt: skip
    mtc.Pipeline(args, {}).execute()

    final = pd.read_csv(out / "contig_data_final.txt", sep="\t")
    truth, found = _genome_labels(final, metagenome, "Final_bin")
    assert "ND" not in found
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    summary = pd.read_csv(out / "bin_summary.txt", sep="\t", index_col=0)
    assert len(summary) == len(metagenome.genomes)
    # MGE contigs are binned with their host.
    host_bin = final.set_index("Name").Final_bin
    for mge in metagenome.mges:
        host_contigs = [c.name for c in metagenome.molecule(mge.host).contigs]
        assert {host_bin[c.name] for c in mge.contigs} <= set(host_bin[host_contigs])
    # The shotgun abundance of a bin is its share of the sequenced bases.
    truth = metagenome.truth
    bases = (truth.Length * truth.Abundance).groupby(final.set_index("Name").Final_bin).sum()
    expected = 100 * bases / (truth.Length * truth.Abundance).sum()
    assert summary.SG_abundance.astype(float).to_numpy() == pytest.approx(expected[summary.index].to_numpy(), rel=0.1)
    # Outputs of the optional steps.
    assert (out / "clustering_matrix.npz").stat().st_size > 0
    scaffolds = sorted(p.name for p in (out / "scaffold_bin").glob("*_scaffolded.fa"))
    assert scaffolds == sorted(f"{name}_scaffolded.fa" for name in summary.index)
    # Clean-up: no temporary directory, no pyfastx index, no intermediate table.
    assert not (tmp_path / "tmp").exists()
    assert not os.path.exists(f"{metagenome_inputs['assembly']}.fxi")
    assert not (out / "contig_data_partition.txt").exists()
