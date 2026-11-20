# Test for figures module

import gzip
import shutil

import hicstuff.hicstuff as hcs
import metacooler
import metator.figures as mtf
import scipy.sparse as sp
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_hex


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
            "Final_cluster": ["MetaTOR_00001_00000"] * 2 + ["MetaTOR_00002_00000"] * 2,
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
    """metator network leaves *_sorted.pairs.gz in its outdir, which the
    heatmap of metator validation, run step by step without a map, reads:
    the heatmap must read them."""
    out_dir, contigs_data, bin_summary = binned_output
    contigs_data["Shotgun_coverage"] = "-"
    mtf.plot_figures(str(out_dir), contigs_data, bin_summary, 50_000)
    assert (out_dir / "plot" / "network_heatmap.png").stat().st_size > 0


def test_plot_figures_takes_only_pairs_files(binned_output):
    """A file whose name only looks like that of a pairs file
    (alignment_0_sorted_pairs.gz) is not taken for one: no heatmap."""
    out_dir, contigs_data, bin_summary = binned_output
    (out_dir / "alignment_0_sorted.pairs.gz").rename(out_dir / "alignment_0_sorted_pairs.gz")
    mtf.plot_figures(str(out_dir), contigs_data, bin_summary, 50_000)
    assert not (out_dir / "plot" / "network_heatmap.png").exists()


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
            "Final_cluster": truth.Genome.map(bins).fillna("ND"),
        }
    ).reset_index(drop=True)
    bin_summary = {}
    for genome, quality in zip(metagenome.genomes, qualities, strict=True):
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


def _boxes(ax):
    """(x, q1, median, q3, fill colour) of each box, from left to right."""
    medians = {}
    for line in ax.lines:
        x = np.asarray(line.get_xdata(), dtype=float)
        # Medians span the whole box; caps span half of it.
        if len(x) == 2 and np.isclose(x[1] - x[0], 1.0):
            medians[round(x.mean())] = line.get_ydata()[0]
    boxes = []
    for box in ax.patches:
        extents = box.get_path().get_extents()
        assert extents.width == pytest.approx(1.0)
        x = round((extents.x0 + extents.x1) / 2)
        boxes.append((x, extents.y0, medians[x], extents.y1, to_hex(box.get_facecolor())))
    return sorted(boxes)


def _check_boxes(ax, expected):
    """Compare the boxes of ax with (x, values, colour) triples."""
    boxes = _boxes(ax)
    assert [(x, colour) for x, *_, colour in boxes] == [(x, colour) for x, _, colour in expected]
    for (_, *stats, _), (_, values, _) in zip(boxes, expected, strict=True):
        assert stats == pytest.approx(np.percentile(values, [25, 50, 75]))


BOXPLOTS = {
    "MAGs_GC_distribution": ("GC_content", "GC content distribution", "linear"),
    "MAGs_HiC_cov_distribution": ("HiC_cov", "HiC coverage distribution", "log"),
    "MAGs_SG_cov_distribution": ("Shotgun_coverage", "Assembly coverage distribution", "log"),
}


@pytest.mark.parametrize(
    "qualities",
    [list(QUALITIES), ["HQ", "Other", "HQ", "Other", "HQ"]],
    ids=["every-class", "some-classes"],
)
def test_quality_boxplots(metagenome, tmp_path, saved_figures, qualities):
    """One box per bin of the size-weighted contig values, coloured by the
    quality of the bin, and a legend of every quality class."""
    contigs_data, bin_summary = _binning(metagenome, qualities)
    mtf.plot_figures(str(tmp_path), contigs_data.copy(), bin_summary, 50_000)
    binned = contigs_data[contigs_data.Final_cluster != "ND"].copy()
    binned["HiC_cov"] = 100 * binned.Hit / binned.Size
    weights = (binned.Size / binned.Size.min()).astype(int)
    quality = dict(zip(bin_summary, qualities, strict=True))
    for name, (column, ylabel, scale) in BOXPLOTS.items():
        values = binned[column].repeat(weights).groupby(binned.Final_cluster.repeat(weights))
        # Bins by decreasing value: the GC content of the bin, or the
        # highest coverage of its contigs.
        if column == "GC_content":
            order = sorted(bin_summary, key=lambda b: bin_summary[b]["GC-content"], reverse=True)
        else:
            order = list(values.max().sort_values(ascending=False).index)
        ax = saved_figures[f"{name}.png"].axes[0]
        _check_boxes(
            ax,
            [(x, values.get_group(b), mtf.MAG_QUALITY_COLOURS[quality[b]]) for x, b in enumerate(order)],
        )
        assert ax.get_xlabel() == "" and ax.get_ylabel() == ylabel
        assert ax.get_yscale() == scale
        assert ax.get_xlim() == (-0.5, len(order) - 0.5)
        assert [label.get_text() for label in ax.get_xticklabels()] == [""] * len(order)
        legend = ax.get_legend()
        assert [text.get_text() for text in legend.get_texts()] == list(mtf.MAG_QUALITY_LABELS.values())
        assert [to_hex(handle.get_facecolor()) for handle in legend.legend_handles] == list(mtf.MAG_QUALITY_COLOURS.values())


def test_quality_boxplots_skip_unknown_classes_and_missing_values(tmp_path, saved_figures):
    """Bins of no known quality class are left out, and so are missing and
    infinite values (seaborn treated both as missing)."""
    data = pd.DataFrame(
        {
            "Final_cluster": ["b1", "b1", "b1", "b1", "b2", "b3", "b3"],
            "MAG_quality": ["MQ", "MQ", "MQ", "MQ", "ND", "LQ", "LQ"],
            "size_weight": [1, 2, 1, 1, 1, 1, 1],
            "value": [1.0, 3.0, np.nan, np.inf, 5.0, 7.0, -np.inf],
        }
    )
    mtf._quality_boxplots(data, "value", "Value", str(tmp_path / "boxes.png"))
    ax = saved_figures["boxes.png"].axes[0]
    _check_boxes(ax, [(0, [1.0, 3.0, 3.0], "#6666ff"), (1, [7.0], "#ccccff")])
    assert ax.get_ylabel() == "Value" and ax.get_yscale() == "linear"
    assert (tmp_path / "boxes.png").stat().st_size > 0


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
    binned = contigs_data[contigs_data.Final_cluster != "ND"]
    assert set(info_contigs) == set(binned.Name)
    assert len(mag_starts) == len(bin_summary)
    assert n_bins == sum(g.length // 20_000 + 1 for g in metagenome.genomes)
    pairs_file = tmp_path / "pairs.txt"
    metagenome.write_pairs(pairs_file)
    matrix = mtf.build_matrix_network([str(pairs_file)], info_contigs, n_bins, 20_000)
    assert matrix.shape == (n_bins, n_bins)
    assert np.allclose(matrix, matrix.T)
    mtf.network_heatmap(matrix, mag_starts=mag_starts)  # shown, not saved


# The heatmap drawn from the map of a run (build_matrix_map).


def test_normalize_scn_is_that_of_hicstuff():
    """The in-place SCN normalization gives the numbers of hicstuff's, a
    window without contacts included."""
    rng = np.random.default_rng(0)
    for size in [1, 7, 64, 150]:
        matrix = (rng.random((size, size)) < 0.3) * rng.poisson(4, (size, size)).astype(float)
        matrix[rng.integers(0, size)] = 0
        matrix = np.minimum(matrix, matrix.T) + np.minimum(matrix, matrix.T).T
        np.testing.assert_allclose(
            mtf.normalize_scn(matrix), hcs.normalize_dense(matrix, norm="SCN", iterations=100), rtol=1e-12, atol=0
        )


def test_heatmap_windows():
    """The windows of build_matrix_network, MAG after MAG, the contigs of
    each end to end: each bin of a map covers its windows by share; a contig
    in no MAG covers none."""
    contigs_data = pd.DataFrame(
        {"Name": ["a", "b", "c", "d"], "Size": [30, 25, 10, 40], "Final_cluster": ["m1", "m1", "ND", "m2"]}
    )
    bin_summary = {"m1": {"Length": 55}, "m2": {"Length": 40}}
    bins = pd.DataFrame({"chrom": ["a", "b", "c", "d", "d"], "start": [0, 0, 0, 0, 25], "end": [30, 25, 10, 25, 40]})
    windows, starts, size = mtf.heatmap_windows(bins, contigs_data, bin_summary, 20)
    # m1: 55 bp, windows 0 to 2; m2: 40 bp, windows 3 to 5.
    assert (starts, size) == ([0, 3], 6)
    expected = np.zeros((5, 6))
    expected[0, [0, 1]] = [20 / 30, 10 / 30]
    expected[1, [1, 2]] = [10 / 25, 15 / 25]
    expected[3, [3, 4]] = [20 / 25, 5 / 25]
    expected[4, 4] = 1
    np.testing.assert_allclose(windows.toarray(), expected)


def _pairs_matrix(pairs_file, bins, windows):
    """The matrix of the heatmap from a pairs file, as a map of these bins
    holds its pairs: each read in its bin (one per contig, or of a fixed
    size from the first base), the contacts of each bin, with another or
    within itself, spread over its windows, those within a window left
    out, then normalized."""
    pairs = pd.read_csv(pairs_file, sep="\t", comment="#", header=None, usecols=[1, 2, 3, 4],
                        names=["chr1", "pos1", "chr2", "pos2"])  # fmt: skip
    first = bins.reset_index().groupby(bins["chrom"].astype(str), sort=False)["index"].min()
    width = (bins["end"] - bins["start"]).max() if len(bins.chrom.unique()) < len(bins) else None

    def bin_of(chrom, pos):
        index = first.loc[chrom.astype(str)].to_numpy()
        return index if width is None else index + np.maximum(pos.to_numpy() - 1, 0) // width

    b1, b2 = bin_of(pairs.chr1, pairs.pos1), bin_of(pairs.chr2, pairs.pos2)
    counts = sp.coo_matrix((np.ones(len(b1)), (np.minimum(b1, b2), np.maximum(b1, b2))), shape=(len(bins),) * 2)
    matrix = (windows.T @ counts.tocsr() @ windows).toarray()
    matrix += matrix.T
    np.fill_diagonal(matrix, 0)
    return mtf.normalize_scn(matrix)


@pytest.mark.parametrize("binsize", [None, 1000], ids=["one bin per contig", "1 kb bins"])
def test_heatmap_of_a_map_is_that_of_its_pairs(metagenome, metagenome_files, binsize, tmp_path):
    """The matrix that build_matrix_map builds from the pixels of a map is
    the matrix that a helper builds from the pairs the map holds, each pair
    counted between the bins of its two reads. The helper spreads the bins
    over the windows with heatmap_windows, and normalizes the matrix with
    normalize_scn, as build_matrix_map does: this test checks how the
    pixels are read and summed, not those two, which test_heatmap_windows
    and test_normalize_scn_is_that_of_hicstuff check (and
    test_heatmap_of_a_map_is_that_of_the_network, without them)."""
    contigs_data, bin_summary = _binning(metagenome, list(QUALITIES))
    # metacooler indexes the assembly (.fai): a private copy.
    assembly = shutil.copyfile(metagenome_files["assembly"], tmp_path / "assembly.fa")
    metacool = tmp_path / "map.metacool"
    if binsize:
        kwargs = {"binsize": binsize}
    else:
        kwargs = {"bins": pd.DataFrame({"chrom": contigs_data.Name, "start": 0, "end": contigs_data.Size})}
    metacooler.from_pairs(str(metacool), pairs=str(metagenome_files["pairs"]), fasta=str(assembly),
                          contig_genome_map=None, contig_ordering="fasta", ps=False, **kwargs)  # fmt: skip
    matrix, starts = mtf.build_matrix_map(str(metacool), contigs_data, bin_summary, 20_000, chunksize=1000)
    bins = metacooler.MetaCooler(str(metacool)).bins()[:]
    windows, expected_starts, size = mtf.heatmap_windows(bins, contigs_data, bin_summary, 20_000)
    assert starts == expected_starts and matrix.shape == (size, size)
    assert size == sum(genome.length // 20_000 + 1 for genome in metagenome.genomes)
    np.testing.assert_allclose(matrix, _pairs_matrix(metagenome_files["pairs"], bins, windows), rtol=1e-9, atol=1e-15)
    assert np.allclose(matrix, matrix.T) and matrix.max() > 0


def test_heatmap_of_a_map_is_that_of_the_network(tmp_path):
    """A guard, independent of heatmap_windows and normalize_scn, which
    test_heatmap_of_a_map_is_that_of_its_pairs shares with the code: where
    no contig spans two windows, the heatmap drawn from a map of one bin per
    contig is the heatmap drawn from its pairs by build_matrix_network,
    which a validation without a map draws (a contig in no MAG, d, in
    neither). The reads lie within their window: a read at the last base of
    a contig (position 100, 1-based, of a contig of 100 bp) falls in the
    next window for build_matrix_network."""
    pairs = tmp_path / "pairs.pairs"
    rows = [("a", 10, "b", 20), ("a", 90, "b", 50), ("a", 50, "c", 50), ("b", 10, "c", 90), ("b", 50, "c", 10),
            ("b", 90, "c", 40), ("a", 10, "a", 90), ("c", 20, "d", 30), ("d", 10, "d", 70)]  # fmt: skip
    pairs.write_text(
        "## pairs format v1.0\n" + "".join(f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t-\n" for i, (c1, p1, c2, p2) in enumerate(rows))
    )
    lengths = {"a": 100, "b": 100, "c": 100, "d": 80}
    metacool = tmp_path / "map.metacool"
    bins = pd.DataFrame({"chrom": list(lengths), "start": 0, "end": list(lengths.values())})
    metacooler.from_pairs(str(metacool), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    contigs_data = pd.DataFrame(
        {"Name": list(lengths), "Size": list(lengths.values()), "Final_cluster": ["m1", "m1", "m2", "ND"]}
    )
    bin_summary = {"m1": {"Length": 200}, "m2": {"Length": 100}}
    matrix, starts = mtf.build_matrix_map(
        str(metacool), contigs_data, {name: dict(mag) for name, mag in bin_summary.items()}, 100
    )
    info_contigs, mag_starts, size = mtf.generates_frags_network(contigs_data, bin_summary, 100)
    network = mtf.build_matrix_network([str(pairs)], info_contigs, size, 100)
    assert starts == mag_starts == [0, 3] and matrix.shape == network.shape == (5, 5)
    np.testing.assert_allclose(matrix, network, rtol=1e-12, atol=1e-15)
    assert (matrix > 0).sum() == 6


@pytest.mark.parametrize("windows", [40, 3000])
def test_heatmap_drawn_at_a_pixel_per_window(windows, tmp_path):
    """The heatmap is saved at about one pixel per window, from 100 to 2000
    dpi: its side has at least as many pixels as the matrix has windows,
    and less than 1.5 times as many (600 for a small matrix, drawn at 100
    dpi). It was saved at 2000 dpi whatever its size, 8879 pixels a side,
    which took gigabytes of memory for a few thousand windows."""
    from PIL import Image

    matrix = np.random.default_rng(0).random((windows, windows))
    out_file = tmp_path / "heatmap.png"
    mtf.network_heatmap(matrix, str(out_file), [0, windows // 2])
    width, height = Image.open(out_file).size
    assert windows <= min(width, height) <= max(width, height) < 1.5 * max(windows, 400)


def test_the_contacts_within_a_bin_are_spread_over_its_windows(tmp_path):
    """The contacts within a bin of the map are spread over its windows, as
    those between two bins are: a contig longer than a window, in one bin
    (as in the map of the contigs of a run), is no dark square on the
    heatmap. They were left out: every contig of 50 kb or more was a dark
    square, and a MAG of one contig an empty block."""
    pairs = tmp_path / "c1.pairs"
    pairs.write_text("## pairs format v1.0\n" + "".join(f"r{i}\tc1\t{1 + 37 * i}\tc1\t{290 - 41 * i}\t+\t-\n" for i in range(7)))
    metacool = tmp_path / "c1.metacool"
    bins = pd.DataFrame({"chrom": ["c1"], "start": [0], "end": [300]})
    metacooler.from_pairs(str(metacool), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    contigs_data = pd.DataFrame({"Name": ["c1"], "Size": [300], "Final_cluster": ["m1"]})
    matrix, starts = mtf.build_matrix_map(str(metacool), contigs_data, {"m1": {"Length": 300}}, 100)
    # c1 covers windows 0 to 2 of the 4 of m1: its contacts within itself
    # fill them, but for the diagonal.
    assert starts == [0] and matrix.shape == (4, 4)
    assert (matrix[:3, :3] > 0).sum() == 6 and not matrix.diagonal().any() and not matrix[3].any()


def test_heatmap_windows_without_a_mag():
    """No MAG gives no window: heatmap_windows failed on an empty bin
    summary (zip() argument 2 is longer than argument 1)."""
    contigs_data = pd.DataFrame({"Name": ["a"], "Size": [30], "Final_cluster": ["ND"]})
    bins = pd.DataFrame({"chrom": ["a"], "start": [0], "end": [30]})
    windows, starts, size = mtf.heatmap_windows(bins, contigs_data, {}, 20)
    assert (starts, size, windows.shape) == ([], 0, (1, 0))


@pytest.mark.parametrize("error", ["MetacoolError", "ValueError"])
def test_a_heatmap_that_cannot_be_drawn_is_a_warning(error, binned_output, monkeypatch, caplog):
    """The heatmap is a figure: when it cannot be drawn, the reason is a
    warning, and the other figures are drawn, where the run stopped (and
    left the map of its contigs unbinned). The reason of a map that cannot
    be read is its message; that of another error, which is unexpected,
    its type, message and traceback, as the map of a run says it. The case
    of a map that cannot be read is a guard: its warning was that already.
    """
    import metator.metacool as mtmc

    out_dir, contigs_data, bin_summary = binned_output

    def build_matrix_map(*args, **kwargs):
        raise {"MetacoolError": mtmc.MetacoolError, "ValueError": ValueError}[error]("no map")

    monkeypatch.setattr(mtf, "build_matrix_map", build_matrix_map)
    with caplog.at_level("WARNING", logger="metator"):
        mtf.plot_figures(str(out_dir), contigs_data, bin_summary, 50_000, metacool="map.metacool")
    (record,) = [record for record in caplog.records if "heatmap" in record.getMessage()]
    reason = "no map" if error == "MetacoolError" else "ValueError: no map"
    assert (record.levelname, record.getMessage()) == ("WARNING", f"The heatmap of the clusters is not drawn: {reason}")
    assert bool(record.exc_info) == (error == "ValueError")
    assert not (out_dir / "plot" / "network_heatmap.png").exists()
    assert (out_dir / "plot" / "bins_distribution.png").exists()


def test_no_heatmap_without_a_mag(binned_output, caplog):
    """With no MAG, no heatmap is drawn, and the log says so: it failed
    with an index error on an empty matrix, given as the reason."""
    out_dir, contigs_data, _ = binned_output
    (out_dir / "plot").mkdir()
    with caplog.at_level("INFO", logger="metator"):
        mtf.plot_heatmap(str(out_dir), str(out_dir / "plot"), contigs_data, {}, metacool="map.metacool")
    assert "No MAG: the heatmap of the clusters is not drawn." in caplog.text and "not drawn:" not in caplog.text
    assert not (out_dir / "plot" / "network_heatmap.png").exists()


def test_the_lines_between_mags_are_a_pixel_wide(tmp_path):
    """The lines between the MAGs are a pixel wide at the dpi of the image,
    which can be 100 dpi: they were 0.1 pt, which no pixel shows below 720
    dpi (the brightest of them barely above the colour of the matrix).
    Before the heatmap was saved at about a pixel per window, it was saved
    at 2000 dpi, where 0.1 pt is a few pixels: this passes there."""
    from PIL import Image

    out_file = tmp_path / "heatmap.png"
    mtf.network_heatmap(np.zeros((40, 40)), str(out_file), [0, 20])
    image = np.asarray(Image.open(out_file).convert("RGB")).astype(int)
    # The matrix, all 0, is the most frequent colour; the lines are white,
    # dotted, at alpha 0.8.
    colours, counts = np.unique(image.reshape(-1, 3), axis=0, return_counts=True)
    matrix = (image == colours[counts.argmax()]).all(axis=2)
    rows, columns = np.flatnonzero(matrix.any(axis=1)), np.flatnonzero(matrix.any(axis=0))
    inside = image[rows[0] : rows[-1] + 1, columns[0] : columns[-1] + 1]
    bright = inside.sum(axis=2) - colours[counts.argmax()].sum() > 150
    assert bright.sum() >= 40
