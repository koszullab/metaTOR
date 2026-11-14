# Test for validation module

import metator.io as mio
import metator.validation as mtv
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import pytest
import os
import re
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
network = nx.read_edgelist(
    "tests_data/outdir_validation/network.txt",
    nodetype=int,
    data=(("weight", float),),
)
fileObj = open("tests_data/outdir_validation/recursive_bins.obj", "rb")
recursive_bins = pickle.load(fileObj)
fileObj.close()


def test_get_bin_coverage():
    bin_info = mtv.get_bin_coverage(bin_summary, contigs_data)
    assert bin_info["MetaTOR_00002_00000"]["HiC_abundance"] == pytest.approx(49.72, abs=1e-2)


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
    return {"P_00001_00000": dict(zip(COLS, values), recursive=True, step=0, parent=None)}


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


def test_update_contigs_data_recursive():
    tmp_dir = "tmp_partition_validation_2"
    contamination = False
    parent_dict = dict()
    os.makedirs(tmp_dir, exist_ok=True)
    (
        contamination,
        contig_data,
        parent_dict,
    ) = mtv.update_contigs_data_recursive(
        "MetaTOR_00002_00000",
        contigs_data,
        recursive_bins,
        assembly,
        tmp_dir,
        tmp_dir,
        500000,
        contamination,
        parent_dict,
        "MetaTOR",
    )
    assert contamination
    assert len(np.unique(contig_data.Recursive_bin_ID)) > 1
    assert len(os.listdir(tmp_dir)) == 4
    shutil.rmtree(tmp_dir)


def test_write_bins_contigs():
    binning_file = "tmp_binning.txt"
    contig_data = mtv.write_bins_contigs(bin_summary, contigs_data, binning_file, "MetaTOR")
    print(np.unique(list(contig_data["Final_bin"])))
    assert len(np.unique(list(contig_data["Final_bin"]))) == 2
    os.remove(binning_file)


# ChekM deprecated functions.


def test_checkm_compare_bins():
    # Test checkM compare bins deprecated function.
    checkm_summary = mtv.checkm_compare_bins(
        overlapping_checkm_file="tests_data/outdir_checkM/overlapping_checkm_results.txt",
        overlapping_taxonomy_file="tests_data/outdir_checkM/overlapping_checkm_taxonomy.txt",
        recursive_checkm_file="tests_data/outdir_checkM/recursif_checkm_results.txt",
        recursive_taxonomy_file="tests_data/outdir_checkM/recursif_checkm_taxonomy.txt",
        prefix="MetaTOR",
    )
    over = mio.read_results_checkm(
        "tests_data/outdir_checkM/overlapping_checkm_results.txt",
        "tests_data/outdir_checkM/overlapping_checkm_taxonomy.txt",
    )
    rec = mio.read_results_checkm(
        "tests_data/outdir_checkM/recursif_checkm_results.txt",
        "tests_data/outdir_checkM/recursif_checkm_taxonomy.txt",
    )

    # Assert there is the rigth numbers of keys for a bin:
    keys = list(checkm_summary.keys())
    assert len(checkm_summary[keys[0]]) == 10

    # Assert there are no more keys then available:
    test_over = [re.search("_0", key) is not None for key in keys]
    assert sum(test_over) <= len(over)
    assert len(checkm_summary) - sum(test_over) <= len(rec)

    # Assert there are no duplicates  (one bin from overlapping and from
    # recursive)
    over_id = []
    rec_id = []
    for i in range(len(test_over)):
        if test_over[i]:
            over_id.append(int(keys[i].split("_")[1]))
        else:
            rec_id.append(int(keys[i].split("_")[1]))
    assert not (set(over_id) & set(rec_id))


@pytest.mark.xfail(
    strict=True,
    reason="HiC_abundance is halved twice: a fully binned metagenome sums to 50%, while SG_abundance sums to 100%",
)
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
