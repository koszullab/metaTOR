# Test for validation module

import metator.io as mio
import metator.validation as mtv
import numpy as np
import pandas as pd
import pickle
import pytest
import os
import shutil

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
network = pd.read_csv(
    "tests_data/outdir_validation/network.txt",
    sep="\t",
    header=None,
    names=["contig1", "contig2", "weight"],
    dtype={"weight": str},
)
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


def _alive(pid):
    """Whether the process pid runs: a zombie, whose status no parent has
    read yet, has ended, and so has one reaped between the open of its
    /proc file and the read (ProcessLookupError)."""
    try:
        with open(f"/proc/{pid}/stat") as handle:
            return handle.read().rsplit(")", 1)[1].split()[0] not in ("Z", "X")
    except (FileNotFoundError, ProcessLookupError):
        return False


@pytest.mark.skipif(not os.path.isdir("/proc"), reason="reads the state of the processes in /proc")
def test_a_stop_kills_what_micomplete_started(tmp_path, monkeypatch):
    """A run stopped while miComplete runs (Ctrl-C, or SIGTERM: an exception
    raised as the run waits for it) kills miComplete and the processes that
    it started, the workers of its pool: subprocess.run killed miComplete
    only, and its workers lived on, orphaned and blocked for good. Here a
    stand-in for miComplete starts a process that would sleep a minute,
    then waits for it."""
    import signal
    import threading
    import time

    started = tmp_path / "started"
    stand_in = tmp_path / "python"
    stand_in.write_text(f"#!/bin/sh\nsleep 60 &\necho $! > {started}.tmp\nmv {started}.tmp {started}\nwait\n")
    stand_in.chmod(0o755)
    # miComplete is run as sys.executable -m metator._micomplete.
    monkeypatch.setattr(mtv.sys, "executable", str(stand_in))
    bins = tmp_path / "bins"
    bins.mkdir()

    def interrupt_once_started():
        deadline = time.monotonic() + 30
        while not started.exists() and time.monotonic() < deadline:
            time.sleep(0.05)
        if started.exists():
            signal.pthread_kill(threading.main_thread().ident, signal.SIGINT)

    previous = signal.signal(signal.SIGINT, signal.default_int_handler)
    interrupter = threading.Thread(target=interrupt_once_started)
    interrupter.start()
    try:
        with pytest.raises(KeyboardInterrupt):
            mtv.micomplete_quality(str(bins), str(tmp_path / "results.tsv"), 1)
    finally:
        interrupter.join()
        signal.signal(signal.SIGINT, previous)
    child = int(started.read_text())
    try:
        deadline = time.monotonic() + 10
        while _alive(child) and time.monotonic() < deadline:
            time.sleep(0.05)
        assert not _alive(child)
    finally:
        if _alive(child):
            os.kill(child, signal.SIGKILL)


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
    tmp_dir = "tmp_partition_validation_1"
    os.makedirs(tmp_dir, exist_ok=True)
    partition = mtv.recursive_clustering_worker(
        "MetaTOR_00002_00000",
        tmp_dir,
        network,
        "louvain",
        iterations,
        resolution_parameter,
        contigs_data,
    )
    _val = int(partition[716].split(";")[0])
    assert len(partition) == 192
    assert len(partition[716].split(";")) == iterations
    shutil.rmtree(tmp_dir)


def test_recursive_clustering_worker_subnetwork(tmp_path):
    """The subnetwork of a bin is the network rows between two of its
    contigs, verbatim and in the order of the network: an edge listed in
    two rows (network.txt of 1.3.10 from two libraries) keeps both, as in
    partition."""
    rows = ["5\t2\t1.5", "1\t2\t0.25", "2\t9\t3.0", "3\t3\t2", "1\t3\t1e-05", "9\t8\t1.0", "2\t5\t0.5"]
    (tmp_path / "network.txt").write_text("\n".join(rows) + "\n")
    network = pd.read_csv(
        tmp_path / "network.txt", sep="\t", header=None, names=["contig1", "contig2", "weight"], dtype={"weight": str}
    )
    data = pd.DataFrame(
        {"ID": [1, 2, 3, 5, 8, 9], "Overlapping_cluster_ID": ["00001"] * 4 + ["00002"] * 2, "Recursive_cluster_ID": "00000"}
    )
    partition = mtv.recursive_clustering_worker("bin_00001_00000", str(tmp_path), network, "louvain", 2, 1.0, data)
    subnetwork = tmp_path / "recursive_clusters" / "bin_00001_00000" / "subnetwork_bin_00001_00000.txt"
    assert subnetwork.read_text().splitlines() == [rows[0], rows[1], rows[3], rows[4], rows[6]]
    assert sorted(partition) == [1, 2, 3, 5]


def test_update_contigs_data_recursive(tmp_path, monkeypatch):
    """The fasta of the recursive bins are written without the pyfastx
    command, which is not on the PATH here (the files were empty)."""
    monkeypatch.setenv("PATH", str(tmp_path))
    # A private copy of the assembly: pyfastx indexes it (.fxi) in place.
    shutil.copyfile(assembly, tmp_path / "assembly.fa")
    out = tmp_path / "bins"
    out.mkdir()
    contamination, contig_data, _ = mtv.update_contigs_data_recursive(
        "MetaTOR_00002_00000",
        contigs_data,
        recursive_bins,
        str(tmp_path / "assembly.fa"),
        str(out),
        500000,
        False,
        dict(),
        "MetaTOR",
    )
    assert contamination
    assert len(np.unique(contig_data.Recursive_cluster_ID)) > 1
    assert sorted(os.listdir(out)) == ["MetaTOR_00002_00001.fa", "MetaTOR_00002_00002.fa"]
    assert all(path.read_text().startswith(">") for path in out.iterdir())


def test_write_bins_contigs(tmp_path):
    """binning.txt holds the contigs of the final bins, with the bin that
    Final_cluster gives them, in the order of the contigs data."""
    binning_file = tmp_path / "binning.txt"
    data = contigs_data.copy()
    data.loc[data.index[0], "Final_cluster"] = "ND"
    # The contigs of a bin partitioned again but kept whole keep the ID of
    # their recursive bin: they were left out of binning.txt.
    kept_whole = data.index[data.Final_cluster == "MetaTOR_00002_00000"][:5]
    data.loc[kept_whole, "Recursive_cluster_ID"] = "00003"
    returned = mtv.write_bins_contigs(bin_summary, data, str(binning_file))
    assert returned is data
    binning = pd.read_csv(binning_file, sep="\t", header=None, names=["Name", "Bin"])
    expected = data[data.Final_cluster != "ND"]
    assert list(binning.Name) == list(expected.Name) and list(binning.Bin) == list(expected.Final_cluster)
    assert set(data.loc[kept_whole, "Name"]) <= set(binning.Name)


def test_correct_final_bin(tmp_path):
    """The final bin of each contig is read in the fasta files of the bins,
    without leaving the index of pyfastx next to them."""
    final_dir = tmp_path / "final_clusters"
    final_dir.mkdir()
    (final_dir / "B_00001_00000.fa").write_text(">c1 description\nACGT\n>c3\nGG\n")
    (final_dir / "B_00002_00001.fa").write_text(">c2\nTT\n")
    data = pd.DataFrame({"Name": ["c1", "c2", "c3", "c4"], "Final_cluster": "ND"})
    summary = {"B_00001_00000": {}, "B_00002_00001": {}}
    corrected = mtv.correct_final_bin(data, str(final_dir), summary)
    assert list(corrected.Final_cluster) == ["B_00001_00000", "B_00002_00001", "B_00001_00000", "ND"]
    assert sorted(path.name for path in final_dir.iterdir()) == ["B_00001_00000.fa", "B_00002_00001.fa"]


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
