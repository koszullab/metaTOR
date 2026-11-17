# End-to-end binning tests on the mock metagenome of tests/synthetic.py,
# checked against its ground truth: network, partition, validation (with a
# stub miComplete) and the whole pipeline.

import os
import random
import shutil
from pathlib import Path

import cooler
import metacooler
import pandas as pd
import pytest

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
def test_partition_recovers_the_genomes(metagenome, metagenome_inputs, metagenome_network, tmp_path, algorithm):
    random.seed(0)  # igraph draws from Python's random module.
    fasta_dir = tmp_path / "overlapping_bin"
    fasta_dir.mkdir()
    # A private copy of the assembly: pyfastx indexes it (.fxi) in place.
    _, contig_data_file = mtp.partition(
        algorithm,
        str(metagenome_inputs["assembly"]),
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


def test_validation_splits_a_merged_bin(metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli):
    """A bin holding two genomes is reported contaminated by (stub) miComplete:
    the recursive decontamination splits it back into the two genomes."""
    random.seed(0)
    calls = []
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome, calls))
    partition_file, fasta_dir = _partition_with_merged_bin(
        metagenome, metagenome_network["contig_data"], tmp_path, merged=("genome_1", "genome_2")
    )
    out = tmp_path / "out"
    # Outputs of an earlier run, replaced thanks to --force.
    for stale in ["recursive_bin", "final_bin"]:
        (out / stale).mkdir(parents=True)
        (out / stale / "metator_00009_00000.fa").write_text(">old\nACGT\n")
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
        "-F",
    ]  # fmt: skip
    run_cli("validation", *args)
    assert not (out / "recursive_bin" / "metator_00009_00000.fa").exists()

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


def test_validation_without_contaminated_bin(metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli):
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
    run_cli("validation", *args)
    assert len(pd.read_csv(out / "bin_summary.txt", sep="\t", index_col=0)) == len(metagenome.genomes)


def test_partition_without_bin_above_the_size_threshold(metagenome_inputs, metagenome_network, tmp_path, run_cli):
    args = [
        "-a", str(metagenome_inputs["assembly"]),
        "-c", str(metagenome_network["contig_data"]),
        "-n", str(metagenome_network["network"]),
        "-o", str(tmp_path / "out"),
        "-T", str(tmp_path / "tmp"),
        "-i", "2",
        "-s", "10000000",
    ]  # fmt: skip
    run_cli("partition", *args)
    assert list((tmp_path / "out" / "overlapping_bin").iterdir()) == []


def test_pipeline_recovers_the_genomes(metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch, run_cli):
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
    run_cli("pipeline", *args)

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
    # The worker processes that scaffold the bins write to the log file too,
    # whether they are forked or not (forkserver, Linux's default from Python
    # 3.14 on, and spawn start them without the parent's log handlers).
    (log_file,) = out.glob("metator_pipeline_*.log")
    logged = log_file.read_text()
    assert [name for name in summary.index if f"Scaffolding {name}..." not in logged] == []
    # The contact map, written by default: one bin per contig, one genome
    # per bin (the bin summary columns with it), and every pair a contact.
    metacool = out / "metator.metacool"
    metacooler.validate(str(metacool))
    clr = cooler.Cooler(str(metacool))
    assert clr.pixels()[:]["count"].sum() == len(metagenome.pairs)
    assert clr.info["metadata"]["metator"]["parameters"]["start"] == "pair"
    mc = metacooler.MetaCooler(str(metacool))
    assert sorted(mc.contigs.name) == sorted(final.Name)
    assert dict(zip(mc.contigs.name, mc.contigs.genome_id.map(mc.genomes.name), strict=True)) == dict(host_bin)
    assert sorted(mc.genomes.name) == sorted(summary.index)
    assert mc.genomes.set_index("name").Completeness[summary.index].to_numpy() == pytest.approx(summary.Completeness.to_numpy())
    # What metaTOR measured on each contig (with a depth file and an
    # enzyme: the four columns), not its bins.
    columns = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    assert list(mc.contigs.columns) == ["name", "length", "genome_id", *columns]
    measured = final.set_index("Name").loc[mc.contigs.name, columns]
    for column in columns:
        assert mc.contigs[column].tolist() == measured[column].tolist(), column
    # Clean-up: no temporary directory, no pyfastx index, no intermediate table.
    assert not (tmp_path / "tmp").exists()
    assert not os.path.exists(f"{metagenome_inputs['assembly']}.fxi")
    assert not (out / "contig_data_partition.txt").exists()


def test_pipeline_without_metacool(metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch, run_cli):
    """--no-metacool: the bins, and no contact map."""
    random.seed(0)
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    out = tmp_path / "out"
    args = ["-1", str(metagenome_indexed_pairs), "-a", str(metagenome_inputs["assembly"]), "-S", "pair", "-s", "50000",
            "-i", "5", "-j", "2", "-p", "sample", "-o", str(out), "-T", str(tmp_path / "tmp"), "--no-metacool"]  # fmt: skip
    run_cli("pipeline", *args)
    assert len(pd.read_csv(out / "bin_summary.txt", sep="\t", index_col=0)) == len(metagenome.genomes)
    assert not list(out.glob("*.metacool*"))


def test_pipeline_from_a_metacool(metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch, run_cli):
    """From the contact map written by a run from pairs (one bin per contig,
    its contigs grouped by bin), the pipeline makes the same bins: its
    network differs by the pairs with a read on the last base of a contig,
    which the map counts. It writes the map of its bins back from the input
    map: the map of the run from pairs, but for its QC, recomputed from the
    pixels, and its provenance."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-T", str(tmp_path / "tmp")]  # fmt: skip
    random.seed(0)
    run_cli("pipeline", "-1", str(metagenome_indexed_pairs), "-S", "pair", "-o", str(tmp_path / "pair"), *common)
    metacool = tmp_path / "pair" / "metator.metacool"
    random.seed(0)
    run_cli("pipeline", "-1", str(metacool), "-S", "metacool", "-o", str(tmp_path / "metacool"), *common)

    pair, mcool = tmp_path / "pair", tmp_path / "metacool"
    assert (mcool / "network.txt").read_text() != (pair / "network.txt").read_text()
    for name in ["contig_data_network.txt", "contig_data_final.txt", "bin_summary.txt", "binning.txt"]:
        assert (mcool / name).read_text() == (pair / name).read_text(), name
    assert sorted(p.name for p in (mcool / "final_bin_unscaffold").iterdir()) == sorted(
        p.name for p in (pair / "final_bin_unscaffold").iterdir()
    )
    final = pd.read_csv(mcool / "contig_data_final.txt", sep="\t")
    truth, found = _genome_labels(final, metagenome, "Final_bin")
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    # The map of the bins, written back from the input map.
    assert sorted(path.name for path in mcool.glob("*.metacool*")) == ["metator.metacool"]
    written, expected = mcool / "metator.metacool", pair / "metator.metacool"
    metacooler.validate(str(written), full=True)
    clr, ref = cooler.Cooler(str(written)), cooler.Cooler(str(expected))
    for table in ["bins", "pixels"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    # The contigs table that metacooler reassign rebuilds stores genome_id
    # after the columns that it carries from the input map.
    chroms, ref_chroms = clr.chroms()[:], ref.chroms()[:]
    assert list(chroms.columns) == ["name", "length", "GC_content", "Hit", "Shotgun_coverage", "Restriction_site", "genome_id"]
    pd.testing.assert_frame_equal(chroms[ref_chroms.columns], ref_chroms)
    genomes, ref_genomes = metacooler.MetaCooler(str(written)).genomes, metacooler.MetaCooler(str(expected)).genomes
    columns = [column for column in ref_genomes.columns if not column.startswith("qc_")]
    pd.testing.assert_frame_equal(genomes[columns], ref_genomes[columns])
    metadata = clr.info["metadata"]
    assert metadata["metator"]["parameters"]["start"] == "metacool"
    assert metadata["qc"]["source"] == "pixels"
    assert [entry["source"] for entry in metadata["reassign"]] == [str(metacool)]


def test_pipeline_from_a_filtered_metacool(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """From a map that lacks contigs of the assembly (made with a contig
    filter, which leaves out the pairs with a read on the contigs it drops),
    the pipeline completes and makes the bins of a run from the pairs that
    the map holds: the contigs it lacks have no hits and are in no bin."""
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
    run_cli("pipeline", "-1", str(pairs), "-S", "pair", "-o", str(tmp_path / "pair"), "--no-metacool", *common)
    random.seed(0)
    run_cli("pipeline", "-1", str(metacool), "-S", "metacool", "-o", str(tmp_path / "metacool"), *common)

    pair, mcool = tmp_path / "pair", tmp_path / "metacool"
    for name in ["contig_data_network.txt", "contig_data_final.txt", "bin_summary.txt", "binning.txt"]:
        assert (mcool / name).read_text() == (pair / name).read_text(), name
    assert sorted(p.name for p in (mcool / "final_bin_unscaffold").iterdir()) == sorted(
        p.name for p in (pair / "final_bin_unscaffold").iterdir()
    )
    final = pd.read_csv(mcool / "contig_data_final.txt", sep="\t")
    lacking = ~final.Name.isin(kept)
    assert 0 < lacking.sum() < len(final)
    assert (final.Hit[lacking] == 0).all() and (final.Final_bin[lacking] == "ND").all()
    assert len(pd.read_csv(mcool / "bin_summary.txt", sep="\t")) > 0
    (log,) = mcool.glob("metator_pipeline_*.log")
    assert f"INFO :: {lacking.sum()} contig(s) of the assembly are not in {metacool}: they get no contacts" in log.read_text()
    # The map of the bins, from the map: its contigs, at 1 kb, in their bins.
    written = mcool / "metator.metacool"
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    assert mc.binsize == 1000 and set(mc.contigs.name) == kept
    final_bin = final.set_index("Name").Final_bin.replace("ND", "__unbinned__")
    assert (
        dict(zip(mc.contigs.name, mc.contigs.genome_id.map(mc.genomes.name), strict=True))
        == final_bin[~lacking.to_numpy()].to_dict()
    )
    assert cooler.Cooler(str(written)).pixels()[:]["count"].sum() == cooler.Cooler(str(metacool)).pixels()[:]["count"].sum()
    # What metaTOR measured on the contigs of the map. contig_data_final.txt
    # lists the contigs that the map lacks too: metacooler annotate reads it,
    # as the history of the map says, and its warning about their rows is
    # dropped.
    columns = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    assert list(mc.contigs.columns) == ["name", "length", "genome_id", *columns]
    measured = final.set_index("Name").loc[mc.contigs.name, columns]
    for column in columns:
        assert mc.contigs[column].tolist() == measured[column].tolist(), column
    assert "Contigs table" not in log.read_text()
    assert mc.metadata["annotate"][-1]["contigs_table"] == str(mcool / "contig_data_final.txt")
    assert not [path.name for path in mcool.iterdir() if path.name.startswith(".")]
