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
fasta_dir = "tests_data/outdir_validation/overlapping_bin"
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
    converters={"Overlapping_bin_ID": str},
)
contigs_data["index"] = contigs_data["ID"] - 1
contigs_data = contigs_data.set_index("index")
contigs_data["Recursive_bin_ID"] = f"{0:05d}"
contigs_data["Recursive_bin_contigs"] = pd.Series("-", index=contigs_data.index, dtype=object)
contigs_data["Recursive_bin_size"] = pd.Series("-", index=contigs_data.index, dtype=object)
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
        {"ID": [1, 2, 3, 5, 8, 9], "Overlapping_bin_ID": ["00001"] * 4 + ["00002"] * 2, "Recursive_bin_ID": "00000"}
    )
    partition = mtv.recursive_clustering_worker("bin_00001_00000", str(tmp_path), network, "louvain", 2, 1.0, data)
    subnetwork = tmp_path / "recursive_bins" / "bin_00001_00000" / "subnetwork_bin_00001_00000.txt"
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
    assert len(np.unique(contig_data.Recursive_bin_ID)) > 1
    assert sorted(os.listdir(out)) == ["MetaTOR_00002_00001.fa", "MetaTOR_00002_00002.fa"]
    assert all(path.read_text().startswith(">") for path in out.iterdir())


def test_write_bins_contigs():
    binning_file = "tmp_binning.txt"
    contig_data = mtv.write_bins_contigs(bin_summary, contigs_data, binning_file, "MetaTOR")
    print(np.unique(list(contig_data["Final_bin"])))
    assert len(np.unique(list(contig_data["Final_bin"]))) == 2
    os.remove(binning_file)


def test_bin_abundances_are_percentages():
    contigs = pd.DataFrame(
        {
            "Final_bin": ["A", "A", "B"],
            "Size": [100, 200, 300],
            "Shotgun_coverage": [1.0, 1.0, 2.0],
            "Hit": [10, 20, 30],
        }
    )
    summary = mtv.get_bin_coverage({"A": {}, "B": {}}, contigs)
    assert summary["A"]["SG_abundance"] == pytest.approx(100 * 300 / 900)
    assert sum(b["SG_abundance"] for b in summary.values()) == pytest.approx(100)
    assert sum(b["HiC_abundance"] for b in summary.values()) == pytest.approx(100)
