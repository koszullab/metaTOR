# Test for validation module

import metator.io as mio
import metator.validation as mtv
import numpy as np
import pandas as pd
import pickle
import pytest
import os

threads = 8
fasta_dir = "tests_data/outdir_validation/overlapping_clusters"
arch_output = "tests_data/outdir_validation/tmp_micomplete/arch131.tsv"
bact_output = "tests_data/outdir_validation/tmp_micomplete/bact105.tsv"
assembly = "tests_data/outdir_validation/assembly_val.fa"
bin_summary = mio.micomplete_results_to_dict("tests_data/outdir_validation/overlapping_micomplete_results.txt")
iterations = 5
resolution_parameter = 1.0
contigs_data = pd.read_csv(
    "tests_data/outdir_validation/contig_data_final.txt",
    sep="\t",
    header=0,
    converters={"Overlapping_cluster_ID": str},
)
contigs_data["index"] = contigs_data["ID"] - 1
contigs_data = contigs_data.set_index("index")
contigs_data["Recursive_cluster_ID"] = f"{0:05d}"
contigs_data["Recursive_cluster_contigs"] = pd.Series("-", index=contigs_data.index, dtype=object)
contigs_data["Recursive_cluster_size"] = pd.Series("-", index=contigs_data.index, dtype=object)
network = pd.read_csv("tests_data/outdir_validation/network.txt", sep="\t", header=None, names=["contig1", "contig2", "weight"])
fileObj = open("tests_data/outdir_validation/recursive_bins.obj", "rb")
recursive_bins = pickle.load(fileObj)
fileObj.close()


def test_get_bin_coverage():
    bin_info = mtv.get_bin_coverage(bin_summary, contigs_data)
    # Every contig is binned: the two bins hold all the Hi-C hits.
    assert bin_info["MetaTOR_00002_00000"]["HiC_abundance"] == pytest.approx(99.45, abs=1e-2)
    assert sum(info["HiC_abundance"] for info in bin_info.values()) == pytest.approx(100)


def test_merge_micomplete():
    out_file = "tmp_micomplete.csv"
    mtv.merge_micomplete(bact_output, arch_output, out_file)
    data = pd.read_csv(out_file, sep="\t", comment="#", index_col=0)
    assert data.loc["MetaTOR_00001_00000", "Markers"] == "Archaea"
    assert data.loc["MetaTOR_00002_00000", "Markers"] == "Bacteria"
    os.remove(out_file)


COLS = ["Completeness", "Weighted completeness", "Redundancy", "Weighted redundancy", "CDs"]


def _recursive_summary(path, completeness):
    """Write a miComplete-like table: one row per recursive bin."""
    with open(path, "w") as table:
        table.write("Name\t" + "\t".join(COLS) + "\n")
        for name, value in completeness.items():
            table.write(f"{name}\t{value}\t{value}\t1.0\t1.0\t100\n")
    return str(path)


def _parent(completeness):
    values = [completeness, completeness, 1.6, 1.6, 300]
    return {"P_00001_00000": dict(zip(COLS, values, strict=True), recursive=True, step=0, parent=None)}


def test_micomplete_compare_bins_keeps_every_child(tmp_path):
    """A parent split into several good bins is replaced by all of them."""
    children = ["P_00001_00001", "P_00001_00002", "P_00001_00003"]
    table = _recursive_summary(tmp_path / "rec.tsv", dict.fromkeys(children, 0.9))
    parents = dict.fromkeys(children, "P_00001_00000")
    summary, contamination = mtv.micomplete_compare_bins(table, _parent(0.9), parents, 1)
    assert sorted(summary) == children
    assert contamination
    assert all(summary[c]["parent"] == "P_00001_00000" for c in children)
    assert all(summary[c]["step"] == 1 for c in children)


def test_micomplete_compare_bins_keeps_parent_over_worse_children(tmp_path):
    """Children that lose too much completeness are discarded."""
    children = ["P_00001_00001", "P_00001_00002"]
    table = _recursive_summary(tmp_path / "rec.tsv", dict.fromkeys(children, 0.3))
    parents = dict.fromkeys(children, "P_00001_00000")
    summary, contamination = mtv.micomplete_compare_bins(table, _parent(0.9), parents, 1)
    assert list(summary) == ["P_00001_00000"]
    assert summary["P_00001_00000"]["recursive"] is False
    assert not contamination


def test_micomplete_quality(tmp_path, monkeypatch):
    bins = os.path.abspath(fasta_dir)
    before = sorted(os.listdir(bins))
    monkeypatch.chdir(tmp_path)
    mtv.micomplete_quality(bins, "results.tsv", threads)
    data = pd.read_csv("results.tsv", sep="\t", comment="#", index_col=0)
    assert len(data) == 2
    assert len(data.columns) == 14
    # miComplete's intermediate files are left neither here nor next to the bins.
    assert sorted(os.listdir(tmp_path)) == ["results.tsv", "tests_data"]
    assert sorted(os.listdir(bins)) == before


def test_micomplete_runs_in_the_process_group_of_the_run(tmp_path, monkeypatch):
    """miComplete runs in the process group of the run, not a session of its
    own: a Ctrl-C in the terminal, or the signal that a job scheduler sends
    to the job, reaches it too. A stand-in for miComplete writes its process
    group, which Python itself reads (no ps)."""
    import shlex
    import sys

    stand_in = tmp_path / "python"
    groups = tmp_path / "groups"
    group_of_the_stand_in = "import os; print(os.getpgrp())"
    stand_in.write_text(f"#!/bin/sh\n{shlex.quote(sys.executable)} -c {shlex.quote(group_of_the_stand_in)} >> {groups}\nexit 3\n")
    stand_in.chmod(0o755)
    monkeypatch.setattr(mtv.sys, "executable", str(stand_in))
    bins = tmp_path / "bins"
    bins.mkdir()
    with pytest.raises(Exception, match="returned non-zero exit status 3"):
        mtv.micomplete_quality(str(bins), str(tmp_path / "results.tsv"), 1)
    assert [int(line) for line in groups.read_text().split()] == [os.getpgrp()]


def test_micomplete_failure_is_an_error(tmp_path, monkeypatch):
    """A guard (subprocess.run(check=True) stopped it so too): miComplete
    exiting with an error stops the validation, with its exit status and
    its command line."""
    import subprocess

    stand_in = tmp_path / "python"
    stand_in.write_text("#!/bin/sh\nexit 3\n")
    stand_in.chmod(0o755)
    monkeypatch.setattr(mtv.sys, "executable", str(stand_in))
    bins = tmp_path / "bins"
    bins.mkdir()
    with pytest.raises(subprocess.CalledProcessError) as failed:
        mtv.micomplete_quality(str(bins), str(tmp_path / "results.tsv"), 1)
    assert failed.value.returncode == 3 and failed.value.cmd[:3] == [str(stand_in), "-m", "metator._micomplete"]


def test_recursive_clustering_worker():
    partition = mtv.recursive_clustering_worker(
        "MetaTOR_00002_00000",
        network,
        "louvain",
        iterations,
        resolution_parameter,
        contigs_data,
    )
    _val = int(partition[716].split(";")[0])
    assert len(partition) == 192
    assert len(partition[716].split(";")) == iterations


def test_recursive_clustering_worker_subnetwork(monkeypatch):
    """The subnetwork of a bin is the network rows between two of its
    contigs, in the order of the network, which the partition reads in
    memory: an edge listed in two rows (network.txt of 1.3.10 from two
    libraries) keeps both, as in partition. No file is written."""
    network = pd.DataFrame(
        [(5, 2, 1.5), (1, 2, 0.25), (2, 9, 3.0), (3, 3, 2.0), (1, 3, 1e-05), (9, 8, 1.0), (2, 5, 0.5)],
        columns=["contig1", "contig2", "weight"],
    )
    data = pd.DataFrame(
        {"ID": [1, 2, 3, 5, 8, 9], "Overlapping_cluster_ID": ["00001"] * 4 + ["00002"] * 2, "Recursive_cluster_ID": "00000"}
    )
    partitioned = []
    algo_partition = mtv.mtp.algo_partition

    def recorded(algorithm, subnetwork, *args):
        partitioned.append(subnetwork)
        return algo_partition(algorithm, subnetwork, *args)

    monkeypatch.setattr(mtv.mtp, "algo_partition", recorded)
    partition = mtv.recursive_clustering_worker("bin_00001_00000", network, "louvain", 2, 1.0, data)
    (subnetwork,) = partitioned
    pd.testing.assert_frame_equal(subnetwork, network.iloc[[0, 1, 3, 4, 6]])
    assert sorted(partition) == [1, 2, 3, 5]


def test_name_recursive_clusters():
    """The recursive bins larger than the size get their IDs in the contig
    data, their parent, and their contigs, by name, for their FASTA files
    (which metacooler writes: see write_recursive_fasta); the others none."""
    data = contigs_data.copy()
    parents, groups = {}, {}
    contamination = mtv.name_recursive_clusters("MetaTOR_00002_00000", data, recursive_bins, 500000, parents, "MetaTOR", groups)
    assert contamination
    assert len(np.unique(data.Recursive_cluster_ID)) > 1
    assert sorted(groups) == ["MetaTOR_00002_00001", "MetaTOR_00002_00002"]
    assert parents == dict.fromkeys(groups, "MetaTOR_00002_00000")
    for name, contigs in groups.items():
        rows = data[data.Recursive_cluster_ID == name.rsplit("_", 1)[1]]
        assert contigs == list(rows.Name) and all(isinstance(contig, str) for contig in contigs)
    # None reaches an enormous size: no bin, no change.
    data, parents, groups = contigs_data.copy(), {}, {}
    assert not mtv.name_recursive_clusters("MetaTOR_00002_00000", data, recursive_bins, 10**12, parents, "MetaTOR", groups)
    assert parents == {} and groups == {} and set(data.Recursive_cluster_ID) == {"00000"}


def test_correct_final_bin():
    """The final bin of each contig: the bin of the FASTA file it was written
    in (written, the contigs of each file by bin name). The contigs of a
    bin partitioned again but kept whole, which keep the ID of their
    recursive bin, are in it as the others are: binning.txt, which lists the
    contigs by final bin, left them out."""
    data = pd.DataFrame(
        {"Name": ["c1", "c2", "c3", "c4"], "Recursive_cluster_ID": ["00000", "00001", "00003", "00000"], "Final_cluster": "ND"}
    )
    written = {"B_00001_00000": ["c1", "c3"], "B_00002_00001": ["c2"], "B_00003_00000": ["c4"]}
    summary = {"B_00001_00000": {}, "B_00002_00001": {}}
    corrected = mtv.correct_final_bin(data, written, summary)
    assert list(corrected.Final_cluster) == ["B_00001_00000", "B_00002_00001", "B_00001_00000", "ND"]


def test_correct_final_bin_refuses_a_contig_that_the_table_lacks():
    """A contig of a final FASTA file that the contig data lacks (a name
    that the table read as another: 00123 as the number 123) is refused,
    named, rather than added to the table without a value."""
    data = pd.DataFrame({"Name": ["123", "c2"], "Final_cluster": "ND"})
    with pytest.raises(
        ValueError,
        match=r"^1 contig\(s\) of the FASTA files of the final clusters are not in the contig data of the validation: '00123'\.$",
    ):
        mtv.correct_final_bin(data, {"B_00001_00000": ["00123", "c2"]}, {"B_00001_00000": {}})


@pytest.mark.parametrize(
    ("completeness", "redundancy", "quality"),
    [
        # One bin of each class (the classes of the binning of the mock
        # metagenome that tests/test_figures.py drew, before the plots went).
        (0.95, 1.02, "HQ"),
        (0.80, 1.05, "MQ"),
        (0.60, 1.05, "LQ"),
        (0.90, 1.30, "Contaminated"),
        (0.30, 1.00, "Other"),
        # The edges: >= 90% and <= 5% contamination is HQ, >= 70% MQ, >= 50%
        # LQ, more than 10% contamination is contaminated, below 50% other.
        (0.90, 1.05, "HQ"),
        (0.90, 1.06, "MQ"),
        (0.89, 1.0, "MQ"),
        (0.70, 1.10, "MQ"),
        (0.69, 1.0, "LQ"),
        (0.50, 1.10, "LQ"),
        (0.50, 1.11, "Contaminated"),
        (0.49, 2.0, "Other"),
        # As miComplete writes them: text.
        ("0.95", "1.02", "HQ"),
    ],
)
def test_bin_quality(completeness, redundancy, quality):
    assert mtv.bin_quality(completeness, redundancy) == quality


def test_give_results_info(caplog):
    """The log counts the bins in the quality classes of bin_quality: a bin of
    95% completeness with a redundancy of 1.08 is a MQ MAG, which the log
    counted as HQ, and a bin of 95% with 1.02 is a HQ MAG, which it counted
    as MQ."""
    summary = {
        "hq": {"Weighted completeness": 0.95, "Weighted redundancy": 1.02, "Length": 3_000_000},
        "mq": {"Weighted completeness": 0.95, "Weighted redundancy": 1.08, "Length": 2_000_000},
        "lq": {"Weighted completeness": 0.6, "Weighted redundancy": 1.08, "Length": 1_000},
        "conta": {"Weighted completeness": 0.6, "Weighted redundancy": 1.2, "Length": 10},
        "other": {"Weighted completeness": 0.2, "Weighted redundancy": 1.0, "Length": 1},
    }
    with caplog.at_level("INFO", logger="metator_logger"):
        mtv.give_results_info(summary)
    assert "5 clusters have been kept after the recursive iterations." in caplog.text
    assert "Total size of the extracted clusters: 5001011" in caplog.text
    for line in [
        "HQ MAGs: 1\tTotal Size: 3000000",
        "MQ MAGs: 1\tTotal Size: 2000000",
        "LQ MAGs: 1\tTotal Size: 1000",
        "Contaminated potential MAGs: 1\tTotal Size: 10",
        "Other clusters: 1\tTotal Size: 1",
    ]:
        assert line in caplog.text


def test_bin_abundances_are_percentages():
    contigs = pd.DataFrame(
        {
            "Final_cluster": ["A", "A", "B"],
            "Size": [100, 200, 300],
            "Shotgun_coverage": [1.0, 1.0, 2.0],
            "Hit": [10, 20, 30],
        }
    )
    summary = mtv.get_bin_coverage({"A": {}, "B": {}}, contigs)
    assert summary["A"]["SG_abundance"] == pytest.approx(100 * 300 / 900)
    assert sum(b["SG_abundance"] for b in summary.values()) == pytest.approx(100)
    assert sum(b["HiC_abundance"] for b in summary.values()) == pytest.approx(100)


def test_write_overlapping_fasta_needs_the_assembly_of_the_map(metagenome_network, tmp_path, logged):
    """The FASTA files of the overlapping clusters to check come from the
    assembly that the map was made from: one that lacks a contig of a
    cluster is refused, with metacooler's reason, logged as an error, and no
    file is written. No cluster that reaches the size: no file, and none
    needed."""
    import metator.clusters as mtcl

    table = metagenome_network["contig_table"]
    ids = list(table.ID)
    clusters = mtcl.Partition(core={1: ids}, overlapping={1: ids}, names={1: "metator_00001_00000"})
    other = tmp_path / "other.fa"
    other.write_text(">x\nACGT\n")
    with pytest.raises(ValueError, match=r"^The clusters to check cannot be written from "):
        mtv.write_overlapping_fasta(str(metagenome_network["map"]), str(other), clusters, table, 1, str(tmp_path / "fasta"))
    assert [level for level, message in logged() if message.startswith("The clusters to check")] == ["ERROR"]
    assert not list((tmp_path / "fasta").iterdir())
    written = mtv.write_overlapping_fasta(
        str(metagenome_network["map"]), str(other), clusters, table, int(table.Size.sum()) + 1, str(tmp_path / "none")
    )
    assert written == {} and not list((tmp_path / "none").iterdir())


def test_write_recursive_fasta(metagenome_network, metagenome_files, tmp_path, logged):
    """The FASTA files of the clusters of a recursive step: those of their
    contigs, through metacooler, each record as it is in the assembly, in
    the order given; none for no cluster. An assembly that lacks a contig is
    refused, with metacooler's reason, logged as an error."""
    table = metagenome_network["contig_table"]
    names = [str(name) for name in table.Name[:3]]
    groups = {"metator_00001_00001": names[::-1], "metator_00001_00002": names[:1]}
    out = tmp_path / "step_1"
    mtv.write_recursive_fasta(str(metagenome_network["map"]), str(metagenome_files["assembly"]), groups, str(out))
    assert sorted(path.name for path in out.iterdir()) == ["metator_00001_00001.fa", "metator_00001_00002.fa"]
    records = [line[1:].split()[0] for line in (out / "metator_00001_00001.fa").read_text().splitlines() if line.startswith(">")]
    assert records == names[::-1]
    mtv.write_recursive_fasta(str(metagenome_network["map"]), str(metagenome_files["assembly"]), {}, str(tmp_path / "none"))
    assert not list((tmp_path / "none").iterdir())
    other = tmp_path / "other.fa"
    other.write_text(">x\nACGT\n")
    with pytest.raises(ValueError, match=r"^The clusters of the recursive step cannot be written from "):
        mtv.write_recursive_fasta(str(metagenome_network["map"]), str(other), groups, str(tmp_path / "refused"))
    assert [level for level, message in logged() if message.startswith("The clusters of the recursive")] == ["ERROR"]
