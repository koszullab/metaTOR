# Test for figures module

import gzip

import metator.figures as mtf
import numpy as np
import pandas as pd
import pytest


def test_reindex_df():
    data = pd.DataFrame({"name": ["a", "b"], "size_weight": [1, 3]})
    expanded = mtf.reindex_df(data, "size_weight")
    assert expanded["name"].tolist() == ["a", "b", "b", "b"]
    assert expanded.index.tolist() == [0, 1, 2, 3]


SIZES = {"c1": 120_000, "c2": 80_000, "c3": 100_000, "c4": 60_000}


@pytest.fixture
def binned_output(tmp_path):
    """A tiny but coherent metator output: contigs, bins and a gzipped pairs file."""
    contigs_data = pd.DataFrame(
        {
            "ID": [1, 2, 3, 4],
            "Name": list(SIZES),
            "Size": list(SIZES.values()),
            "GC_content": [40.0, 42.0, 55.0, 56.0],
            "Hit": [500, 400, 300, 200],
            "Shotgun_coverage": [10.0, 12.0, 30.0, 28.0],
            "Final_bin": ["MetaTOR_00001_00000"] * 2 + ["MetaTOR_00002_00000"] * 2,
        }
    )
    bin_summary = {
        "MetaTOR_00001_00000": {
            "Length": 200_000,
            "GC-content": 41.0,
            "Weighted completeness": 0.95,
            "Weighted redundancy": 1.02,
        },
        "MetaTOR_00002_00000": {
            "Length": 160_000,
            "GC-content": 55.5,
            "Weighted completeness": 0.60,
            "Weighted redundancy": 1.01,
        },
    }
    rng = np.random.default_rng(0)
    names = list(SIZES)
    with gzip.open(tmp_path / "alignment_0_sorted.pairs.gz", "wt") as pairs:
        pairs.write("## pairs format v1.0\n")
        for i in range(2_000):
            c1, c2 = rng.choice(names, 2)
            pos1, pos2 = rng.integers(1, SIZES[c1]), rng.integers(1, SIZES[c2])
            pairs.write(f"r{i}\t{c1}\t{pos1}\t{c2}\t{pos2}\t+\t-\n")
    return tmp_path, contigs_data, bin_summary


def test_plot_figures_reads_compressed_pairs(binned_output):
    """The pipeline leaves *_sorted.pairs.gz in its outdir: the heatmap must read them."""
    out_dir, contigs_data, bin_summary = binned_output
    contigs_data["Shotgun_coverage"] = "-"
    mtf.plot_figures(str(out_dir), contigs_data, bin_summary, 50_000)
    assert (out_dir / "plot" / "network_heatmap.png").stat().st_size > 0


def test_plot_figures_with_shotgun_coverage(binned_output):
    """With --depth, the shotgun coverage boxplot is drawn too."""
    out_dir, contigs_data, bin_summary = binned_output
    (out_dir / "alignment_0_sorted.pairs.gz").unlink()
    mtf.plot_figures(str(out_dir), contigs_data, bin_summary, 50_000)
    assert (out_dir / "plot" / "MAGs_SG_cov_distribution.png").stat().st_size > 0
