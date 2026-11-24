"""The plots of the final clusters of a validation (metator.figures)."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_hex

import metator.figures as mtf
import metator.validation as mtv

# Final clusters of every quality class (see metator.validation.bin_quality),
# by name, as the validation holds them: their weighted completeness and
# redundancy, and their length. Out of the order of their completeness, two
# of the same completeness.
CLUSTERS = {
    "metator_00001_00000": (0.60, 1.05, 1_000_000),  # LQ
    "metator_00002_00000": (0.95, 1.02, 3_000_000),  # HQ
    "metator_00003_00001": (0.30, 1.00, 500_000),  # Other
    "metator_00003_00002": (0.90, 1.30, 4_000_000),  # Contaminated
    "metator_00004_00000": (0.80, 1.05, 2_000_000),  # MQ
    "metator_00005_00000": (0.80, 1.00, 1_500_000),  # MQ
}
# The contigs of the run: 15 Mb, 3 of them in no cluster.
CONTIGS = pd.DataFrame({"Size": [10_000_000, 5_000_000]})


def _summary(clusters):
    """The summary of the final clusters, as recursive_decontamination holds it:
    by name, miComplete's columns (the length a float, as miComplete gives it)
    and those of the validation."""
    return {
        name: {
            "Length": float(length),
            "GC-content": 50.0,
            "Weighted completeness": completeness,
            "Weighted redundancy": redundancy,
            "recursive": False,
            "step": 0,
            "HiC_abundance": 1.0,
        }
        for name, (completeness, redundancy, length) in clusters.items()
    }


def _table(clusters):
    """The final clusters as the plots read them: a row each."""
    return pd.DataFrame(
        [[length, completeness, redundancy] for completeness, redundancy, length in clusters.values()],
        index=list(clusters),
        columns=["Length", "Weighted completeness", "Weighted redundancy"],
    )


def test_plot_figures(tmp_path, saved_figures):
    """Both plots, in plots/ of the output directory, made if missing (as 1.4.0
    drew them, plot_figures takes the output directory): the total size is
    that of the contigs of the run."""
    out = tmp_path / "out"
    mtf.plot_figures(str(out), CONTIGS, _summary(CLUSTERS), 500_000)
    plots = [out / "plots" / "bins_distribution.pdf", out / "plots" / "bins_size_distribution.pdf"]
    assert sorted(out.rglob("*")) == [out / "plots", *plots] and mtf.plot_paths(str(out)) == [str(path) for path in plots]
    assert all(path.stat().st_size > 0 for path in plots)
    texts = [text.get_text() for text in saved_figures["bins_size_distribution.pdf"].axes[0].texts]
    assert "Total size of the contigs: 15.0Mb" in texts


def test_plot_figures_replaces_its_plots_only(tmp_path):
    """The plots replace files of their names in plots/, a symbolic link there
    too, not its target, and the other files there are left as they are."""
    plot_dir = tmp_path / "plots"
    plot_dir.mkdir()
    (plot_dir / "bins_distribution.pdf").write_text("an earlier plot\n")
    (tmp_path / "kept.pdf").write_text("the target of a link\n")
    (plot_dir / "bins_size_distribution.pdf").symlink_to(tmp_path / "kept.pdf")
    (plot_dir / "notes.txt").write_text("the user's\n")
    mtf.plot_figures(str(tmp_path), CONTIGS, _summary(CLUSTERS), 500_000)
    for name in ["bins_distribution.pdf", "bins_size_distribution.pdf"]:
        assert not (plot_dir / name).is_symlink() and (plot_dir / name).read_bytes().startswith(b"%PDF")
    assert (tmp_path / "kept.pdf").read_text() == "the target of a link\n"
    assert (plot_dir / "notes.txt").read_text() == "the user's\n"


def test_bins_distribution(tmp_path, saved_figures):
    """A dot per cluster for its completeness (red) and its contamination (black),
    in percent, by decreasing completeness (clusters of the same completeness
    in the order of the summary), and dashed lines at the thresholds of the
    quality classes."""
    out_file = tmp_path / "dots.pdf"
    mtf.figures_bins_distribution(_table(CLUSTERS), str(out_file))
    assert out_file.stat().st_size > 0
    ax = saved_figures["dots.pdf"].axes[0]
    completion, contamination = ax.collections
    assert completion.get_label() == "Completion" and to_hex(completion.get_facecolor()[0]) == "#ff0000"
    assert contamination.get_label() == "Contamination" and to_hex(contamination.get_facecolor()[0]) == "#000000"
    # 00002 (95%), 00003_00002 (90%), 00004 and 00005 (80%), 00001 (60%),
    # 00003_00001 (30%).
    np.testing.assert_allclose(completion.get_offsets(), [[0, 95], [1, 90], [2, 80], [3, 80], [4, 60], [5, 30]])
    np.testing.assert_allclose(contamination.get_offsets(), [[0, 2], [1, 30], [2, 5], [3, 0], [4, 5], [5, 0]], atol=1e-9)
    assert sorted(line.get_ydata()[0] for line in ax.lines) == [5, 10, 50, 70, 90, 100]
    assert all(line.get_linestyle() == "--" for line in ax.lines)
    assert list(ax.get_yticks()) == [5, 10, 50, 70, 90, 100]
    assert (ax.get_xlabel(), ax.get_ylabel()) == ("Rank index", "Completeness/Contamination (%)")
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["Completion", "Contamination"]


def test_bins_distribution_draws_the_contamination_from_0_to_105(tmp_path, saved_figures):
    """The contamination is drawn from 0 to 105%: a cluster in which miComplete
    finds no marker gene (a redundancy of 0) at 0, where it was drawn at -100%,
    and a contamination of more than 105% at 105%, as the plot of metaTOR 1.x
    meant to (it compared the redundancy, a ratio, with 105)."""
    clusters = {"no_marker": (0.0, 0.0, 30_000), "twice": (0.62, 2.6, 900_000), "a_third": (0.5, 1.35, 600_000)}
    mtf.figures_bins_distribution(_table(clusters), str(tmp_path / "dots.pdf"))
    _, contamination = saved_figures["dots.pdf"].axes[0].collections
    np.testing.assert_allclose(contamination.get_offsets()[:, 1], [105, 35, 0])


def test_bins_size_distribution(tmp_path, saved_figures):
    """A wedge per quality class, by the size of its clusters, then the contigs
    in no cluster; the legend counts the clusters of each class and their
    size, and the text below gives the total size, the share of the MAGs and
    the size threshold."""
    out_file = tmp_path / "pie.pdf"
    classified = mtf.figures_bins_size_distribution(_table(CLUSTERS), 15_000_000, 500_000, str(out_file))
    assert out_file.stat().st_size > 0
    assert dict(classified["MAG_quality"]) == {
        "metator_00001_00000": "LQ",
        "metator_00002_00000": "HQ",
        "metator_00003_00001": "Other",
        "metator_00003_00002": "Contaminated",
        "metator_00004_00000": "MQ",
        "metator_00005_00000": "MQ",
    }
    ax = saved_figures["pie.pdf"].axes[0]
    wedges = ax.patches
    # HQ, MQ, LQ, Contaminated, Other: 3, 3.5, 1, 4 and 0.5 Mb; and 3 Mb in
    # no cluster.
    spans = [(wedge.theta2 - wedge.theta1) / 360 for wedge in wedges]
    np.testing.assert_allclose(spans, np.array([3, 3.5, 1, 4, 0.5, 3]) / 15)
    assert [to_hex(wedge.get_facecolor()) for wedge in wedges] == [
        "#103b6f",
        "#6666ff",
        "#ccccff",
        "#ed3139",
        "#c0bcbc",
        "#eaeded",
    ]
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [
        "HQ MAGs: 1 - 3.0Mb",
        "MQ MAGs: 2 - 3.5Mb",
        "LQ MAGs: 1 - 1.0Mb",
        "Contaminated clusters: 1 - 4.0Mb",
        "Other clusters: 1 - 0.5Mb",
        "Other contigs - 3.0Mb",
    ]
    # The pie draws an empty label per wedge.
    assert [text.get_text() for text in ax.texts if text.get_text()] == [
        "Total size of the contigs: 15.0Mb",
        "Percentage of MAGs: 50.0%",
        "Clusters threshold: 500.0kb",
    ]
    assert ax.get_title() == "Size proportion of clusters depending on their quality"


def test_the_quality_classes_of_the_pie_are_those_of_the_log(tmp_path, saved_figures, caplog):
    """The pie chart counts the clusters of each class as the log of the validation
    does (see metator.validation.give_results_info), from miComplete's
    numbers: a cluster of a weighted completeness of 0.8996 is a MQ MAG, which
    cluster_summary.txt, at two decimals (0.90), would make a HQ MAG."""
    summary = _summary({"almost": (0.8996, 1.0, 2_000_000), "complete": (0.97, 1.01, 3_000_000)})
    with caplog.at_level("INFO", logger="metator_logger"):
        mtv.give_results_info(summary)
    mtf.plot_figures(str(tmp_path), CONTIGS, summary, 500_000)
    legend = [text.get_text() for text in saved_figures["bins_size_distribution.pdf"].axes[0].get_legend().get_texts()]
    assert legend[:2] == ["HQ MAGs: 1 - 3.0Mb", "MQ MAGs: 1 - 2.0Mb"]
    assert "HQ MAGs: 1\tTotal Size: 3000000" in caplog.text and "MQ MAGs: 1\tTotal Size: 2000000" in caplog.text


def test_the_plots_without_a_cluster(tmp_path, saved_figures):
    """Without a final cluster, the dot plot has no dot, and the pie chart is that
    of the contigs in no cluster."""
    mtf.plot_figures(str(tmp_path), CONTIGS, {}, 500_000)
    completion, contamination = saved_figures["bins_distribution.pdf"].axes[0].collections
    assert len(completion.get_offsets()) == len(contamination.get_offsets()) == 0
    ax = saved_figures["bins_size_distribution.pdf"].axes[0]
    assert [round((wedge.theta2 - wedge.theta1) / 360, 9) for wedge in ax.patches] == [0, 0, 0, 0, 0, 1]
    assert ax.get_legend().get_texts()[-1].get_text() == "Other contigs - 15.0Mb"


@pytest.mark.parametrize("broken", ["figures_bins_distribution", "figures_bins_size_distribution"])
def test_a_plot_that_cannot_be_drawn_is_a_warning(broken, tmp_path, monkeypatch, caplog):
    """The plots are figures: one that cannot be drawn is a warning, with the type,
    message and traceback of the error, and the other one is drawn (a run goes
    on to write the map of its clusters: see
    tests/test_binning_synthetic.py::test_a_run_goes_on_when_a_plot_cannot_be_drawn)."""

    def fails(*args):
        raise ValueError("no plot")

    monkeypatch.setattr(mtf, broken, fails)
    plot_dir = tmp_path / "plots"
    mtf.plot_figures(str(tmp_path), CONTIGS, _summary(CLUSTERS), 500_000)
    name = {"figures_bins_distribution": "bins_distribution.pdf", "figures_bins_size_distribution": "bins_size_distribution.pdf"}
    (warning,) = [record for record in caplog.records if record.levelname == "WARNING"]
    assert warning.getMessage() == f"{plot_dir / name[broken]} is not drawn: ValueError: no plot"
    assert warning.exc_info and warning.exc_info[0] is ValueError
    assert [path.name for path in plot_dir.iterdir()] == [other for other in name.values() if other != name[broken]]


def test_a_summary_that_the_plots_cannot_read_is_a_warning(tmp_path, caplog):
    """A summary of the clusters without a column that the plots read (here, the
    Length of a cluster) is a warning for each plot, as any plot that cannot
    be drawn: plot_figures raises nothing, and draws none."""
    summary = _summary(CLUSTERS)
    del summary["metator_00001_00000"]["Length"]
    mtf.plot_figures(str(tmp_path), CONTIGS, summary, 500_000)
    warnings = [record.getMessage() for record in caplog.records if record.levelname == "WARNING"]
    assert warnings == [f"{path} is not drawn: KeyError: 'Length'" for path in mtf.plot_paths(str(tmp_path))]
    assert not list((tmp_path / "plots").iterdir())


def test_a_plot_directory_that_cannot_be_made_is_a_warning(tmp_path, caplog):
    """A plot directory that cannot be made (a file there) is a warning, which gives
    the reason of the system, and no plot is drawn."""
    plot_dir = tmp_path / "plots"
    plot_dir.write_text("a file")
    mtf.plot_figures(str(tmp_path), CONTIGS, _summary(CLUSTERS), 500_000)
    (warning,) = [record for record in caplog.records if record.levelname == "WARNING"]
    assert warning.getMessage() == f"The plots of the clusters are not drawn: {plot_dir} cannot be made (File exists)."
    assert plot_dir.read_text() == "a file"


def test_the_plots_open_no_window(tmp_path, monkeypatch, caplog):
    """The plots are drawn without pyplot, whose figures an interactive backend
    shows in a window: a figure of pyplot (plt.figure, which plt.subplots and
    plt.pie open theirs with) is an error here, and none is left open."""

    def figure(*args, **kwargs):
        raise AssertionError("a figure of pyplot")

    monkeypatch.setattr(plt, "figure", figure)
    plot_dir = tmp_path / "plots"
    mtf.plot_figures(str(tmp_path), CONTIGS, _summary(CLUSTERS), 500_000)
    assert not [record for record in caplog.records if record.levelname == "WARNING"]
    assert sorted(path.name for path in plot_dir.iterdir()) == ["bins_distribution.pdf", "bins_size_distribution.pdf"]
    assert plt.get_fignums() == []
