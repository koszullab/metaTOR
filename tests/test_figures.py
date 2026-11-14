# Test for figures module

import gzip
import io

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


# Figures of a binning of the mock metagenome (tests/synthetic.py), with
# bins of every quality class or of only some of them.

QUALITIES = {
    "HQ": (0.95, 1.02),
    "MQ": (0.80, 1.05),
    "LQ": (0.60, 1.05),
    "Contaminated": (0.90, 1.30),
    "Other": (0.30, 1.00),
}


def _binning(metagenome, qualities, with_depth=True):
    """Contig data and bin summary: one bin per genome, of the given
    qualities (in genome order); MGE contigs are left unbinned."""
    truth = metagenome.truth
    bins = {g.name: f"metator_{k:05d}_00000" for k, g in enumerate(metagenome.genomes, start=1)}
    contigs_data = pd.DataFrame(
        {
            "ID": range(1, len(truth) + 1),
            "Name": truth.Name,
            "Size": truth.Length,
            "GC_content": truth.GC,
            "Hit": (truth.Length * truth.Abundance / 100).astype(int),
            "Shotgun_coverage": 10 * truth.Abundance if with_depth else "-",
            "Final_bin": truth.Genome.map(bins).fillna("ND"),
        }
    ).reset_index(drop=True)
    bin_summary = {}
    for genome, quality in zip(metagenome.genomes, qualities):
        completeness, redundancy = QUALITIES[quality]
        bin_summary[bins[genome.name]] = {
            "Length": genome.length,
            "GC-content": 100 * genome.gc,
            "Weighted completeness": completeness,
            "Weighted redundancy": redundancy,
        }
    return contigs_data, bin_summary


@pytest.mark.parametrize(
    "qualities",
    [list(QUALITIES), ["HQ", "Other", "HQ", "Other", "HQ"]],
    ids=["every-class", "some-classes"],
)
def test_plot_figures_quality_classes(metagenome, tmp_path, qualities):
    contigs_data, bin_summary = _binning(metagenome, qualities)
    mtf.plot_figures(str(tmp_path), contigs_data, bin_summary, 50_000)
    for name in [
        "bins_distribution",
        "bins_size_distribution",
        "MAGs_GC_distribution",
        "MAGs_HiC_cov_distribution",
        "MAGs_SG_cov_distribution",
    ]:
        assert (tmp_path / "plot" / f"{name}.png").stat().st_size > 0
    # No pairs in the output directory: no heatmap.
    assert not (tmp_path / "plot" / "network_heatmap.png").exists()


def test_bins_size_distribution_classes(metagenome, tmp_path):
    _, bin_summary = _binning(metagenome, list(QUALITIES))
    table = pd.DataFrame.from_dict(bin_summary, orient="index").reset_index(drop=True)
    classified = mtf.figures_bins_size_distribution(table, 700_000, 50_000, str(tmp_path / "pie.png"))
    assert list(classified.MAG_quality) == list(QUALITIES)
    assert (tmp_path / "pie.png").stat().st_size > 0


def test_frags_network_and_heatmap(metagenome, tmp_path):
    """Binned contigs are laid out bin after bin; pairs on unbinned contigs
    are skipped."""
    contigs_data, bin_summary = _binning(metagenome, list(QUALITIES))
    info_contigs, mag_starts, n_bins = mtf.generates_frags_network(contigs_data, bin_summary, 20_000)
    binned = contigs_data[contigs_data.Final_bin != "ND"]
    assert set(info_contigs) == set(binned.Name)
    assert len(mag_starts) == len(bin_summary)
    assert n_bins == sum(g.length // 20_000 + 1 for g in metagenome.genomes)
    pairs_file = tmp_path / "pairs.txt"
    metagenome.write_pairs(pairs_file)
    matrix = mtf.build_matrix_network([str(pairs_file)], info_contigs, n_bins, 20_000)
    assert matrix.shape == (n_bins, n_bins)
    assert np.allclose(matrix, matrix.T)
    mtf.network_heatmap(matrix, mag_starts=mag_starts)  # shown, not saved


def _checkv_summary():
    """A CheckV quality summary (completeness of the viral contigs), read
    from its tab-separated file."""
    table = io.StringIO(
        "contig_id\tcontig_length\tprovirus\tcompleteness\n"
        "v0\t50000\tNo\t95.0\n"
        "v1\t40000\tNo\t80.0\n"
        "v2\t30000\tNo\t60.0\n"
        "v3\t20000\tNo\t40.0\n"
        "v4\t10000\tNo\tNA\n"
        "v5\t5000\tYes\t99.0\n"
    )
    return pd.read_csv(table, sep="\t")


@pytest.mark.xfail(
    strict=True,
    reason="build_vmags_summary counts a vMAG of >= 90% completeness in both the >90 and the >70 classes",
)
def test_build_vmags_summary():
    summary = mtf.build_vmags_summary(_checkv_summary())
    assert list(summary.bins) == [1, 1, 1, 1, 1]
    assert summary.bins.sum() == 5  # proviruses are left out


def test_vmags_figures(tmp_path):
    labels = ["sample_1", "sample_2"]
    mtf.barplot_bins_number(labels, [_checkv_summary(), _checkv_summary()], str(tmp_path / "number.png"))
    mtf.barplot_bins_size(labels, [_checkv_summary(), _checkv_summary()], str(tmp_path / "size.png"))
    mtf.pie_bins_size_distribution(_checkv_summary(), str(tmp_path / "pie.png"))
    for name in ["number.png", "size.png", "pie.png"]:
        assert (tmp_path / name).stat().st_size > 0
