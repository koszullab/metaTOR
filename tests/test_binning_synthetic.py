# End-to-end binning tests on the mock metagenome of tests/synthetic.py,
# checked against its ground truth: network, partition, validation (with a
# stub miComplete) and the whole pipeline.

import json
import os
import random
import shutil
import subprocess
from pathlib import Path

import cooler
import metacooler
import pandas as pd
import pytest
from scipy import sparse

import metator.metacool as mtmc
import metator.network as mtn
import metator.partition as mtp
from metator.version import __version__
from tests import synthetic
from tests.hdf5 import walk

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
    fasta_dir = tmp_path / "overlapping_clusters"
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
        1,
        "metator",
    )
    table = pd.read_csv(contig_data_file, sep="\t")
    truth, found = _genome_labels(table, metagenome, "Overlapping_cluster_ID")
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
        data[f"{level}_cluster_ID"] = [f"{i:05d}" for i in ids]
        data[f"{level}_cluster_contigs"] = counts
        data[f"{level}_cluster_size"] = sizes
    partition_file = out_dir / "contig_data_partition.txt"
    data.to_csv(partition_file, sep="\t", index=False)
    fasta_dir = out_dir / "overlapping_clusters"
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
    for stale in ["recursive_clusters", "final_clusters"]:
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
        "-C",
    ]  # fmt: skip
    run_cli("validation", *args)
    assert not (out / "recursive_clusters" / "metator_00009_00000.fa").exists()

    # miComplete ran on the partition bins, then on the bins of the first
    # recursive step only: the split bins are clean.
    assert [Path(c[0]).name for c in calls] == ["overlapping_clusters", "step_1"]
    final = pd.read_csv(out / "contig_data_final.txt", sep="\t")
    truth, found = _genome_labels(final, metagenome, "Final_cluster")
    assert "ND" not in found
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    assert synthetic.purity(truth, found) >= 0.9
    summary = pd.read_csv(out / "cluster_summary.txt", sep="\t", index_col=0)
    assert len(summary) == len(metagenome.genomes)
    assert "metator_00001_00000" not in summary.index
    assert summary["Weighted redundancy"].max() == 1.0
    # Hi-C abundances are proportional to the hits of the bins.
    hits = final.groupby("Final_cluster").Hit.sum()[summary.index]
    hic = summary.HiC_abundance.astype(float)
    assert (hic / hic.sum()).to_numpy() == pytest.approx((hits / hits.sum()).to_numpy(), rel=1e-3)
    # The final fasta files are the final bins, and nothing else.
    assert sorted(p.name for p in (out / "final_clusters").iterdir()) == sorted(f"{name}.fa" for name in summary.index)
    binning = pd.read_csv(out / "binning.txt", sep="\t", header=None, names=["Name", "Bin"])
    assert set(binning.Bin) == set(summary.index)
    assert sorted(binning.Name) == sorted(final.Name[final.Final_cluster != "ND"])
    assert (out / "plot" / "bins_distribution.png").stat().st_size > 0
    # The bins of the split are recursive bins, of step 1; the others are not.
    assert list(summary.recursive) == list(summary.step > 0)
    assert summary.recursive.any() and not summary.recursive.all()
    # The clustering matrix of the recursive partition: a row and a column per
    # contig ID, and values only between contigs of the bin partitioned again.
    # It was 2 columns wide, and validation failed with --cluster-matrix.
    matrix = sparse.load_npz(out / "clustering_matrix_recursive.npz").tocoo()
    assert matrix.shape == (len(final) + 1, len(final) + 1)
    split = set(final.ID[final.Overlapping_cluster_ID == 1])
    assert matrix.nnz > 0 and set(matrix.row) | set(matrix.col) <= split
    assert 0 < matrix.data.min() and matrix.data.max() <= 1


def test_validation_of_a_partition_of_1x(metagenome, metagenome_inputs, metagenome_network, tmp_path, monkeypatch, run_cli):
    """validation reads the contig data of a partition of metaTOR 1.x, its
    columns named after bins (Overlapping_bin_ID...): it writes what it
    writes from the same table of 2.0, byte for byte."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    partition_file, fasta_dir = _partition_with_merged_bin(
        metagenome, metagenome_network["contig_data"], tmp_path, merged=("genome_1", "genome_2")
    )
    tables = {"2.0": str(partition_file), "1.x": synthetic.as_1x(partition_file, tmp_path / "contig_data_partition_1x.txt")}
    written = {}
    for version, table in tables.items():
        random.seed(0)
        fasta = shutil.copytree(fasta_dir, tmp_path / f"overlapping_{version}")
        out = tmp_path / f"out_{version}"
        args = [
            "-a", str(metagenome_inputs["assembly"]),
            "-c", table,
            "-f", str(fasta),
            "-n", str(metagenome_network["network"]),
            "-o", str(out),
            "-T", str(tmp_path / "tmp"),
            "-s", "50000",
            "-i", "10",
            "-t", "1",
        ]  # fmt: skip
        run_cli("validation", *args)
        names = ["contig_data_final.txt", "cluster_summary.txt", "binning.txt", "clusters.txt", "cluster_members.txt"]
        written[version] = {name: (out / name).read_bytes() for name in names}
        written[version] |= {f"final_clusters/{p.name}": p.read_bytes() for p in (out / "final_clusters").glob("*.fa")}
    assert written["1.x"] == written["2.0"]
    assert b"Final_cluster" in written["2.0"]["contig_data_final.txt"].split(b"\n")[0]


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
    assert len(pd.read_csv(out / "cluster_summary.txt", sep="\t", index_col=0)) == len(metagenome.genomes)


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
    assert list((tmp_path / "out" / "overlapping_clusters").iterdir()) == []


def test_pipeline_without_bin_above_the_size_threshold(
    metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch, run_cli
):
    """No bin of the partition reaches --size: the run stops before
    miComplete (which found no bin, then failed to read its empty results),
    saying why in its log file. It leaves the map of the contigs that it
    built from the pairs files, in no bin: every pair a contact."""
    calls = []
    monkeypatch.setattr("metator.validation.micomplete_quality", lambda *args: calls.append(args))
    out = tmp_path / "out"
    args = ["-1", str(metagenome_indexed_pairs), "-a", str(metagenome_inputs["assembly"]), "-S", "pair", "-s", "10000000",
            "-i", "2", "-o", str(out), "-T", str(tmp_path / "tmp")]  # fmt: skip
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)
    assert not calls
    (log,) = out.glob("metator_pipeline_*.log")
    message = f"{out / 'overlapping_clusters'} holds no cluster to validate (no .fa, .fasta or .fna file): partition writes none"
    assert f"ERROR :: {message}" in log.read_text()
    metacool = out / "metator.metacool"
    metacooler.validate(str(metacool), full=True)
    assert list(metacooler.MetaCooler(str(metacool)).genomes.name) == ["__unbinned__"]
    assert cooler.Cooler(str(metacool)).pixels()[:]["count"].sum() == len(metagenome.pairs)


def test_pipeline_network_is_that_of_its_pairs(
    metagenome_inputs, metagenome_indexed_pairs, metagenome_network, network_from_a_map, tmp_path, run_cli
):
    """A run from the pairs of one library reads its network from the map
    of its contigs, not from its pairs files: network.txt and
    contig_data_network.txt are, byte for byte, those of its pairs at
    edge 0 (metagenome_network, from alignment_to_contacts), the pairs
    with a read on the last base of a contig included. network_0.txt, the
    network of its library, is still written, with the bytes of the pairs
    path's."""
    out = tmp_path / "out"
    args = ["-1", str(metagenome_indexed_pairs), "-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]),
            "-S", "pair", "-e", "DpnII", "-o", str(out), "-T", str(tmp_path / "tmp")]  # fmt: skip
    with pytest.raises(RuntimeError, match=r"^network built$"):
        run_cli("pipeline", *args)
    assert network_from_a_map == [str(out / "metator.metacool")]
    assert (out / "network.txt").read_bytes() == metagenome_network["network"].read_bytes()
    assert (out / "contig_data_network.txt").read_bytes() == metagenome_network["contig_data"].read_bytes()
    assert sorted(path.name for path in out.glob("network*.txt")) == ["network.txt", "network_0.txt"]
    assert (out / "network_0.txt").read_bytes() == (metagenome_network["network"].parent / "network_0.txt").read_bytes()


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
    truth, found = _genome_labels(final, metagenome, "Final_cluster")
    assert "ND" not in found
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    summary = pd.read_csv(out / "cluster_summary.txt", sep="\t", index_col=0)
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
    # Outputs of the optional steps.
    assert (out / "clustering_matrix.npz").stat().st_size > 0
    scaffolds = sorted(p.name for p in (out / "scaffolded_clusters").glob("*_scaffolded.fa"))
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
    assert [column for column in mc.contigs.columns if not column.startswith("qc_")] == ["name", "length", "genome_id", *columns]
    measured = final.set_index("Name").loc[mc.contigs.name, columns]
    for column in columns:
        assert mc.contigs[column].tolist() == measured[column].tolist(), column
    # Clean-up: no temporary directory, no pyfastx index, no intermediate table.
    assert not (tmp_path / "tmp").exists()
    assert not os.path.exists(f"{metagenome_inputs['assembly']}.fxi")
    assert not (out / "contig_data_partition.txt").exists()


def test_pipeline_map_is_the_map_of_the_pairs(
    metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch, run_cli
):
    """The map that a run from pairs leaves is the map that metator metacool
    (write_metacool) writes from its pairs files and final tables: the same
    bins, contacts, contigs table (column for column), genomes table (the
    columns of the bin summary and the contact statistics of each bin) and
    QC block. The run reads its pairs files into a map once, as soon as it
    has them, and writes the map of its bins over it at the end with
    metacooler reassign, without reading them again: only the history of
    the map differs, the annotate and reassign keys of its metadata (and
    the date of the file)."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    write_metacool, from_pairs = mtmc.write_metacool, metacooler.from_pairs
    monkeypatch.setattr("metator.metacool.write_metacool", lambda *args, **kwargs: pytest.fail("pairs read again"))
    maps = []

    def counted_from_pairs(uri, *args, **kwargs):
        maps.append(str(uri))
        return from_pairs(uri, *args, **kwargs)

    monkeypatch.setattr(metacooler, "from_pairs", counted_from_pairs)
    out = tmp_path / "out"
    args = ["-1", str(metagenome_indexed_pairs), "-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]),
            "-S", "pair", "-e", "DpnII", "-s", "50000", "-i", "10", "-j", "5", "-o", str(out), "-T", str(tmp_path / "tmp")]  # fmt: skip
    random.seed(0)
    run_cli("pipeline", *args)
    written = out / "metator.metacool"
    assert maps == [str(written)]
    metacooler.validate(str(written), full=True)
    clr = cooler.Cooler(str(written))
    expected = tmp_path / "pairs.metacool"
    write_metacool(
        [str(metagenome_indexed_pairs)],
        str(out / "contig_data_final.txt"),
        str(expected),
        bin_summary=str(out / "cluster_summary.txt"),
        assembly_name=clr.info["genome-assembly"],
        metadata={key: value for key, value in clr.info["metadata"]["metator"].items() if key != "version"},
    )

    ref = cooler.Cooler(str(expected))
    bins = clr.bins()[:]
    # The run's map also stores the normalization of its network, which
    # write_metacool does not compute: with metacooler.store_weights, from
    # metacooler's next release.
    assert ("metator_empirical_hit" in bins.columns) == mtmc.stores_weights()
    pd.testing.assert_frame_equal(bins.drop(columns="metator_empirical_hit", errors="ignore"), ref.bins()[:])
    for table in ["pixels", "chroms"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    if mtmc.stores_weights():
        _assert_normalized_as_the_network(written, out, tmp_path)
    genomes = metacooler.MetaCooler(str(written)).genomes
    pd.testing.assert_frame_equal(genomes, metacooler.MetaCooler(str(expected)).genomes)
    assert {"Completeness", "qc_inter_genome"} <= set(genomes.columns)
    info, ref_info = dict(clr.info), dict(ref.info)
    metadata, ref_metadata = info.pop("metadata"), ref_info.pop("metadata")
    assert metadata["qc"] == ref_metadata["qc"] and metadata["qc"]["n_pairs_read"] == len(metagenome.pairs)
    assert {key for key in info.keys() | ref_info.keys() if info.get(key) != ref_info.get(key)} <= {"creation-date"}
    differ = {key for key in metadata.keys() | ref_metadata.keys() if metadata.get(key) != ref_metadata.get(key)}
    assert differ == {"annotate", "reassign"}
    assert [entry["source"] for entry in metadata["reassign"]] == [str(written)]
    assert len(metadata["annotate"]) == 3 and len(ref_metadata["annotate"]) == 1


def _assert_normalized_as_the_network(written, out, tmp_path):
    """The contacts of the map of a run normalized by the weights it stores
    and summed by pair of contigs are the edges of its network.txt (but
    those of weight 0)."""
    names = pd.read_csv(out / "contig_data_final.txt", sep="\t", usecols=["ID", "Name"]).set_index("ID").Name
    network = pd.read_csv(out / "network.txt", sep="\t", header=None, names=["id1", "id2", "weight"])
    expected_edges = {
        tuple(sorted((names[id1], names[id2]))): weight for id1, id2, weight in network.itertuples(index=False) if weight
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


def test_pipeline_restarts_from_the_map_of_a_stopped_run(
    metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch, run_cli
):
    """A run from pairs that stops with an error once it has built the map of
    its contigs (miComplete fails in the validation, here) leaves that map,
    in no bin. A new run from it (the map as its first argument, with -F in
    the same output directory, whose bin directories the stopped run left)
    gives the
    outputs of a run from the pairs that does not stop, and its map, but
    for the provenance of the map."""
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-T", str(tmp_path / "tmp")]  # fmt: skip
    pairs, clean, out = str(metagenome_indexed_pairs), tmp_path / "clean", tmp_path / "out"
    micomplete = synthetic.fake_micomplete(metagenome)
    monkeypatch.setattr("metator.validation.micomplete_quality", micomplete)
    random.seed(0)
    run_cli("pipeline", "-1", pairs, "-S", "pair", "-o", str(clean), *common)

    def micomplete_fails(*args):
        raise RuntimeError("miComplete failed")

    monkeypatch.setattr("metator.validation.micomplete_quality", micomplete_fails)
    random.seed(0)
    with pytest.raises(RuntimeError, match=r"^miComplete failed$"):
        run_cli("pipeline", "-1", pairs, "-S", "pair", "-o", str(out), *common)
    metacool = out / "metator.metacool"
    metacooler.validate(str(metacool), full=True)
    assert list(metacooler.MetaCooler(str(metacool)).genomes.name) == ["__unbinned__"]
    monkeypatch.setattr("metator.validation.micomplete_quality", micomplete)
    random.seed(0)
    run_cli("pipeline", str(metacool), "-o", str(out), "-F", *common)

    for name in ["network.txt", "contig_data_network.txt", "contig_data_final.txt", "cluster_summary.txt", "binning.txt"]:
        assert (out / name).read_text() == (clean / name).read_text(), name
    final_bins = [{path.name: path.read_bytes() for path in (run / "final_clusters").iterdir()} for run in [out, clean]]
    assert final_bins[0] == final_bins[1] and final_bins[0]
    # The map of the bins, written over the map that the stopped run left.
    metacooler.validate(str(metacool), full=True)
    clr, ref = cooler.Cooler(str(metacool)), cooler.Cooler(str(clean / "metator.metacool"))
    for table in ["bins", "pixels", "chroms"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    pd.testing.assert_frame_equal(
        metacooler.MetaCooler(str(metacool)).genomes, metacooler.MetaCooler(str(clean / "metator.metacool")).genomes
    )
    info, ref_info = dict(clr.info), dict(ref.info)
    metadata, ref_metadata = info.pop("metadata"), ref_info.pop("metadata")
    assert metadata["qc"] == ref_metadata["qc"]
    assert {key for key in info.keys() | ref_info.keys() if info.get(key) != ref_info.get(key)} <= {"creation-date"}
    # The parameters of the run (its start, its input) and the paths of the
    # history.
    differ = {key for key in metadata.keys() | ref_metadata.keys() if metadata.get(key) != ref_metadata.get(key)}
    assert differ == {"metator", "annotate", "reassign"}
    assert [entry["source"] for entry in metadata["reassign"]] == [str(metacool)]


def test_pipeline_from_a_metacool(metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch, run_cli):
    """From the contact map written by a run from pairs (one bin per contig,
    its contigs grouped by cluster), the pipeline makes the same network,
    byte for byte, and the same clusters. It updates the input map in place
    (a copy of the map of the run from pairs): the map of the run from
    pairs, QC and clustering rounds included, but for its provenance."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    common = ["-a", str(metagenome_inputs["assembly"]), "-d", str(metagenome_inputs["depth"]), "-e", "DpnII", "-s", "50000",
              "-i", "10", "-j", "5", "-T", str(tmp_path / "tmp")]  # fmt: skip
    random.seed(0)
    run_cli("pipeline", "-1", str(metagenome_indexed_pairs), "-S", "pair", "-o", str(tmp_path / "pair"), *common)
    expected = tmp_path / "pair" / "metator.metacool"
    metacool = Path(shutil.copyfile(expected, tmp_path / "input.metacool"))
    random.seed(0)
    run_cli("pipeline", str(metacool), "-o", str(tmp_path / "metacool"), *common)

    pair, mcool = tmp_path / "pair", tmp_path / "metacool"
    for name in ["network.txt", "contig_data_network.txt", "contig_data_final.txt", "cluster_summary.txt", "binning.txt",
                 "clusters.txt", "cluster_members.txt"]:  # fmt: skip
        assert (mcool / name).read_text() == (pair / name).read_text(), name
    assert sorted(p.name for p in (mcool / "final_clusters").iterdir()) == sorted(
        p.name for p in (pair / "final_clusters").iterdir()
    )
    final = pd.read_csv(mcool / "contig_data_final.txt", sep="\t")
    truth, found = _genome_labels(final, metagenome, "Final_cluster")
    assert synthetic.adjusted_rand(truth, found) >= 0.9
    # The heatmap of the clusters, drawn from the map of the run, from pairs
    # or from a map (a run from a map drew none): the same map, the same
    # image.
    heatmaps = [run / "plot" / "network_heatmap.png" for run in [pair, mcool]]
    assert heatmaps[0].read_bytes() == heatmaps[1].read_bytes()
    # The input map, updated in place: no map goes to the output directory.
    assert not list(mcool.glob("*.metacool*")) and not list(tmp_path.glob("*.lock"))
    written = metacool
    metacooler.validate(str(written), full=True)
    _read_by_cooler_and_metacooler(written)
    clr, ref = cooler.Cooler(str(written)), cooler.Cooler(str(expected))
    for table in ["bins", "pixels"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    # The contigs table of the input map: metacooler reassign relabels its
    # contigs in place (the same bins), and writes genome_id again, after the
    # other columns; metaTOR then stores its own columns last, as in the map
    # of the pairs.
    chroms, ref_chroms = clr.chroms()[:], ref.chroms()[:]
    columns = ["name", "length", "genome_id", "GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    for table in [chroms, ref_chroms]:
        assert [column for column in table.columns if not column.startswith("qc_")] == columns
    pd.testing.assert_frame_equal(chroms[ref_chroms.columns], ref_chroms)
    # The QC block and columns of the run from pairs, carried over.
    genomes, ref_genomes = metacooler.MetaCooler(str(written)).genomes, metacooler.MetaCooler(str(expected)).genomes
    pd.testing.assert_frame_equal(genomes, ref_genomes)
    metadata = clr.info["metadata"]
    parameters = metadata["metator"]["parameters"]
    assert parameters["map"] == str(metacool) and parameters["start"] is None and parameters["forward"] is None
    assert metadata["qc"] == ref.info["metadata"]["qc"] and metadata["qc"]["source"] == "pairs"
    # The history of the input map, whose map of the clusters was written
    # over the map of its contigs, then that of this run.
    assert [entry["source"] for entry in metadata["reassign"]] == [str(expected), str(metacool)]
    # The clustering rounds of the run, those of the run from pairs: the
    # overlapping clusters that miComplete checked, each recursive step, and
    # the final clusters, the map's genomes.
    mc, ref_mc = metacooler.MetaCooler(str(written)), metacooler.MetaCooler(str(expected))
    assert mc.cluster_rounds == ref_mc.cluster_rounds
    names = [entry["name"] for entry in mc.cluster_rounds]
    assert names[0] == "overlapping" and names[-1] == "final"
    assert all(name == f"recursive_{step}" for step, name in enumerate(names[1:-1], start=1))
    pd.testing.assert_frame_equal(mc.clusters, ref_mc.clusters)
    pd.testing.assert_frame_equal(_members(mc), _members(ref_mc))
    overlapping = mc.clusters[mc.clusters["round"] == "overlapping"]
    assert set(overlapping.name) == {path.stem for path in (mcool / "overlapping_clusters").glob("*.fa")}
    genome_of = dict(zip(mc.contigs.name, mc.contigs.genome_id.map(mc.genomes.name), strict=True))
    final_members = mc.cluster_members("final")
    assert set(final_members.cluster) == set(mc.genomes.name) - {"__unbinned__"}
    assert all(genome_of[contig] == cluster for cluster, contig in zip(final_members.cluster, final_members.contig, strict=True))

    # The same map, written by metator metacool from the map and the tables
    # of the run: a run step by step can write it too (but for the rounds,
    # which metator metacool does not store).
    command = tmp_path / "command.metacool"
    run_cli("metacool", "-c", mcool / "contig_data_final.txt", "-b", mcool / "cluster_summary.txt", "-o", command, expected)
    metacooler.validate(str(command), full=True)
    clr = cooler.Cooler(str(command))
    run_map = cooler.Cooler(str(written))
    # But the normalization of the network, which the run stores in its map.
    pd.testing.assert_frame_equal(clr.bins()[:], run_map.bins()[:].drop(columns="metator_empirical_hit", errors="ignore"))
    for table in ["pixels", "chroms"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(run_map, table)()[:])
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(command)).genomes, genomes)
    assert clr.info["metadata"]["qc"] == metadata["qc"]


def test_pipeline_from_the_columns_of_a_metacool(
    metagenome, metagenome_inputs, metagenome_indexed_pairs, tmp_path, monkeypatch, run_cli
):
    """From the map that a run from pairs wrote with a depth file and an
    enzyme, a run without --depth and --enzyme reads the restriction sites
    and the coverage of the contigs from its Restriction_site and
    Shotgun_coverage columns (theoritical_hit needs both): the same network,
    tables and bins as with them, byte for byte, and the map of its bins
    stores them back."""
    monkeypatch.setattr("metator.validation.micomplete_quality", synthetic.fake_micomplete(metagenome))
    common = ["-a", str(metagenome_inputs["assembly"]), "-m", "theoritical_hit", "-s", "50000", "-i", "10", "-j", "5",
              "-T", str(tmp_path / "tmp")]  # fmt: skip
    options = ["-d", str(metagenome_inputs["depth"]), "-e", "DpnII"]
    random.seed(0)
    run_cli("pipeline", "-1", str(metagenome_indexed_pairs), "-S", "pair", "-o", str(tmp_path / "pair"), *common, *options)
    metacool = tmp_path / "pair" / "metator.metacool"
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
    for name in ["network.txt", "contig_data_network.txt", "contig_data_final.txt", "cluster_summary.txt", "binning.txt"]:
        assert (by_map / name).read_bytes() == (by_options / name).read_bytes(), name
    assert "SG_abundance" in pd.read_csv(by_map / "cluster_summary.txt", sep="\t").columns
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
    contacts = cooler.Cooler(str(metacool)).pixels()[:]["count"].sum()
    random.seed(0)
    run_cli("pipeline", str(metacool), "-o", str(tmp_path / "metacool"), *common)

    pair, mcool = tmp_path / "pair", tmp_path / "metacool"
    for name in ["contig_data_network.txt", "contig_data_final.txt", "cluster_summary.txt", "binning.txt"]:
        assert (mcool / name).read_text() == (pair / name).read_text(), name
    assert sorted(p.name for p in (mcool / "final_clusters").iterdir()) == sorted(
        p.name for p in (pair / "final_clusters").iterdir()
    )
    final = pd.read_csv(mcool / "contig_data_final.txt", sep="\t")
    lacking = ~final.Name.isin(kept)
    assert 0 < lacking.sum() < len(final)
    assert (final.Hit[lacking] == 0).all() and (final.Final_cluster[lacking] == "ND").all()
    assert len(pd.read_csv(mcool / "cluster_summary.txt", sep="\t")) > 0
    (log,) = mcool.glob("metator_pipeline_*.log")
    assert f"INFO :: {lacking.sum()} contig(s) of the assembly are not in {metacool}: they get no contacts" in log.read_text()
    # The map, updated in place: its contigs, at 1 kb, in their clusters.
    written = metacool
    metacooler.validate(str(written), full=True)
    _read_by_cooler_and_metacooler(written)
    mc = metacooler.MetaCooler(str(written))
    assert mc.binsize == 1000 and set(mc.contigs.name) == kept
    final_bin = final.set_index("Name").Final_cluster.replace("ND", "__unbinned__")
    assert (
        dict(zip(mc.contigs.name, mc.contigs.genome_id.map(mc.genomes.name), strict=True))
        == final_bin[~lacking.to_numpy()].to_dict()
    )
    assert cooler.Cooler(str(written)).pixels()[:]["count"].sum() == contacts
    # What metaTOR measured on the contigs of the map. contig_data_final.txt
    # lists the contigs that the map lacks too: metacooler annotate reads it,
    # as the history of the map says, and ignores their rows without a
    # warning.
    columns = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    assert [column for column in mc.contigs.columns if not column.startswith("qc_")] == ["name", "length", "genome_id", *columns]
    measured = final.set_index("Name").loc[mc.contigs.name, columns]
    for column in columns:
        assert mc.contigs[column].tolist() == measured[column].tolist(), column
    assert "Contigs table" not in log.read_text()
    assert mc.metadata["annotate"][-1]["contigs_table"] == str(mcool / "contig_data_final.txt")
    assert not [path.name for path in mcool.iterdir() if path.name.startswith(".")]


def _contig_map(metagenome_inputs, path):
    """The map of the contigs of the mock community at path, in no cluster
    yet, as a pipeline run from its pairs builds it: one bin per contig,
    what metaTOR measures on the assembly (with a depth file and an enzyme)
    and the contact statistics of the pairs."""
    contig_data, _ = mtn.create_contig_data(str(metagenome_inputs["assembly"]), 1, str(metagenome_inputs["depth"]), "DpnII")
    mtmc.write_contig_map([str(metagenome_inputs["pairs"])], contig_data, str(path), assembly_name="synthetic")
    return path


def _without_run(objects):
    """The objects of a map (see tests.hdf5.walk) but what differs between two
    runs of the same clusters: the date the map was written
    (creation-date, an attribute of the file), and the command line and
    the parameters of the run that wrote it (the command and parameters
    of the metator key of its metadata, a JSON text)."""
    kind, members, attributes = objects["/"]
    kept = []
    for key, dtype, value in attributes:
        if key == "metadata":
            metadata = json.loads(value)
            for run in ["command", "parameters"]:
                del metadata["metator"][run]
            value = metadata
        if key != "creation-date":
            kept.append((key, dtype, value))
    return {**objects, "/": (kind, members, kept)}


def _metator_metadata(objects):
    """The metator key of the metadata of a map (see tests.hdf5.walk)."""
    (metadata,) = [value for key, _, value in objects["/"][2] if key == "metadata"]
    return json.loads(metadata)["metator"]


def _first_cluster_contaminated(metagenome):
    """A stub miComplete (see synthetic.fake_micomplete) that finds the first
    overlapping cluster contaminated (a weighted redundancy of 1.5): the
    validation partitions it again, which makes the round recursive_1 (its
    parts, too incomplete for one of them to replace it)."""

    def micomplete_quality(fasta_dir, outfile, threads):
        synthetic.write_micomplete_summary(fasta_dir, outfile, metagenome)
        if Path(fasta_dir).name == "overlapping_clusters":
            # Three comment lines, the header, then a row per cluster.
            lines = [line.split("\t") for line in Path(outfile).read_text().splitlines()]
            lines[4][lines[3].index("Weighted redundancy")] = "1.5"
            Path(outfile).write_text("".join("\t".join(line) + "\n" for line in lines))

    return micomplete_quality


def test_steps_on_a_map_give_the_map_of_pipeline(metagenome, metagenome_inputs, tmp_path, monkeypatch, run_cli):
    """The steps run one by one on a map (network MAP, partition MAP,
    validation MAP) give the map that pipeline MAP gives, from copies of the
    same map of the contigs, and the same tables: every group, soft link,
    dataset (its dtype, shape, chunks, compression and values) and attribute
    of the HDF5 file, but the date it was written (creation-date) and the
    command line and parameters of the run that wrote it (the command and
    parameters of the metator key of its metadata: pipeline's, validation's).
    Each run is made in a directory of its own, with the same relative
    paths, which the history of the map records (the input map of
    metacooler reassign, the tables of metacooler annotate). The rounds are
    overlapping, recursive_1 (a cluster that the stub miComplete finds
    contaminated, partitioned again) and final.

    Run in one process after one random.seed(0), as igraph draws from
    Python's random module: the steps draw the stream that pipeline draws
    (the temporary directory of each command drew from it too, so a run
    step by step, with three, drew another one than pipeline, with one).

    Each step leaves a valid map: network MAP reads it only; partition MAP
    stores the overlapping clusters of the partition, those that reach
    --size (one per FASTA file of overlapping_clusters), in its clustering rounds
    (the round overlapping, without miComplete's columns, which the
    validation adds), and changes nothing else; validation MAP writes the
    map of the clusters, as pipeline MAP does at its end."""
    monkeypatch.setattr("metator.validation.micomplete_quality", _first_cluster_contaminated(metagenome))
    assembly = str(metagenome_inputs["assembly"])
    contig_map = _contig_map(metagenome_inputs, tmp_path / "contigs.metacool")
    runs = {name: tmp_path / name for name in ["pipeline", "steps"]}
    for directory in runs.values():
        directory.mkdir()
        shutil.copyfile(contig_map, directory / "input.metacool")
    common = ["-a", assembly, "-o", "out", "-T", "tmp", "-s", "20000"]

    monkeypatch.chdir(runs["pipeline"])
    random.seed(0)
    run_cli("pipeline", "input.metacool", *common, "-d", metagenome_inputs["depth"], "-e", "DpnII", "-i", "10", "-j", "5")

    monkeypatch.chdir(runs["steps"])
    metacool = runs["steps"] / "input.metacool"
    random.seed(0)
    run_cli("network", "input.metacool", "-a", assembly, "-o", "out", "-T", "tmp", "-d", metagenome_inputs["depth"],
            "-e", "DpnII")  # fmt: skip
    assert metacool.read_bytes() == contig_map.read_bytes()
    run_cli("partition", "input.metacool", *common, "-c", "out/contig_data_network.txt", "-n", "out/network.txt", "-i", "10")
    metacooler.validate(str(metacool), full=True)
    mc = metacooler.MetaCooler(str(metacool))
    assert mc.cluster_rounds == [{"name": "overlapping", "tool": "metaTOR", "version": __version__, "step": "partition"}]
    fasta_files = sorted((runs["steps"] / "out" / "overlapping_clusters").glob("*.fa"))
    assert list(mc.clusters.name) == [path.stem for path in fasta_files] and fasta_files
    # genome, which metacooler computes: no cluster is a genome of the map yet.
    assert list(mc.clusters.columns) == ["round", "name", "parent_round", "parent", "n_contigs", "length", "genome"]
    assert mc.clusters.genome.isna().all()
    members = mc.cluster_members("overlapping")
    for path in fasta_files:
        assert sorted(members.contig[members.cluster == path.stem]) == sorted(synthetic.read_fasta(path)), path.name
    # Nothing else of the map changed: its genomes and contigs tables, its
    # contacts, its metadata.
    before, after = walk(contig_map), walk(metacool)
    kind, stored, attributes = after.pop("/")
    assert (kind, [name for name in stored if name != "clusters"], attributes) == before.pop("/")
    assert {name: value for name, value in after.items() if not name.startswith("/clusters")} == before
    run_cli("validation", "input.metacool", *common, "-c", "out/contig_data_partition.txt", "-n", "out/network.txt",
            "-f", "out/overlapping_clusters", "-i", "5")  # fmt: skip

    pipeline, steps = runs["pipeline"] / "out", runs["steps"] / "out"
    for name in ["cluster_summary.txt", "contig_data_final.txt", "binning.txt", "clusters.txt", "cluster_members.txt",
                 "plot/network_heatmap.png"]:  # fmt: skip
        assert (steps / name).read_bytes() == (pipeline / name).read_bytes(), name
    written = {name: directory / "input.metacool" for name, directory in runs.items()}
    for path in written.values():
        metacooler.validate(str(path), full=True)
    found = {name: walk(path) for name, path in written.items()}
    assert _without_run(found["steps"]) == _without_run(found["pipeline"])
    # What differs: the command lines and their parameters, those of
    # pipeline and of validation.
    metator = {name: _metator_metadata(objects) for name, objects in found.items()}
    assert metator["pipeline"]["command"].startswith("metator pipeline input.metacool")
    assert metator["steps"]["command"].startswith("metator validation input.metacool")
    assert metator["steps"]["parameters"]["map"] == metator["pipeline"]["parameters"]["map"] == "input.metacool"
    assert "start" not in metator["steps"]["parameters"]
    # The rounds of the run, the final clusters those of the genomes of the
    # map.
    mc = metacooler.MetaCooler(str(written["steps"]))
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping", "recursive_1", "final"]
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
