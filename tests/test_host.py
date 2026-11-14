# Tests of the MGE host detection (metator.host) on the mock metagenome:
# each MGE contacts one host genome much more than the others.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import to_hex

import metator.host as mth


@pytest.fixture
def binned(metagenome, metagenome_network):
    """Contig data with the true bins as Final_bin, the network as a table,
    and the bin names of the genomes and MGEs."""
    data = pd.read_csv(metagenome_network["contig_data"], sep="\t")
    mags = {g.name: f"metator_{k:05d}_00000" for k, g in enumerate(metagenome.genomes, start=1)}
    mges = {m.name: f"MetaTOR_MGE_{k:05d}" for k, m in enumerate(metagenome.mges, start=1)}
    data["Final_bin"] = data.Name.map(metagenome.truth.Genome).map({**mags, **mges})
    network = pd.read_csv(metagenome_network["network"], sep="\t", names=["contig1", "contig2", "signal"])
    return data, network, mags, mges


def _bin_summary(mag_names, completeness, redundancy):
    return pd.DataFrame({"MAG": mag_names, "Completeness": completeness, "Redundancy": redundancy})


def test_create_bins(metagenome, binned):
    data, _, mags, mges = binned
    unbinned = data.index[data.Name.map(metagenome.truth.Kind) == "genome"][0]
    data.loc[unbinned, "Final_bin"] = "ND"  # unbinned contigs are ignored
    mag_bins, mge_bins = mth.create_bins(data)
    assert set(mag_bins) == set(mags.values())
    assert set(mge_bins) == set(mges.values())
    for genome in metagenome.genomes:
        expected = {c.name for c in genome.contigs} - {data.Name[unbinned]}
        assert set(mag_bins[mags[genome.name]].contigs) == expected
    ids = data.set_index("Name").ID
    for mge in metagenome.mges:
        contigs = mge_bins[mges[mge.name]].contigs
        assert contigs == {c.name: ids[c.name] for c in mge.contigs}
    contig_to_mag = mth.map_contigs_to_mags(mag_bins)
    assert len(contig_to_mag) == sum(len(mag.contigs) for mag in mag_bins.values())
    last = metagenome.genomes[0].contigs[-1].name
    assert contig_to_mag.get(ids[last]) == (mags[metagenome.genomes[0].name] if last != data.Name[unbinned] else None)


def test_bin_classes():
    mag, mge = mth.Mag("metator_00001_00000"), mth.MgeMag("MetaTOR_MGE_00001")
    mag.add_contig("contig_1", 1)
    mge.add_contig("contig_2", 2)
    mge.add_host(mag.name, mag)
    mag.add_mge(mge.name, mge)
    assert repr(mag) == "Mag('metator_00001_00000', 1 contigs, 1 MGEs)"
    assert repr(mge) == "MgeMag('MetaTOR_MGE_00001', 1 contigs, 1 hosts)"
    assert repr(mth.Bin("bin")) == "Bin('bin', 0 contigs)"


def test_mges_find_their_host(metagenome, binned, tmp_path):
    data, network, mags, mges = binned
    mag_bins, mge_bins = mth.create_bins(data)
    mag_bins, mge_bins = mth.compute_mge_mag_interactions(
        network,
        mag_bins,
        mge_bins,
        interaction_threshold=10.0,
        min_interacting_contigs=5,
        output_file=str(tmp_path / "interactions.txt"),
        image_file=str(tmp_path / "histogram.png"),
    )
    for mge in metagenome.mges:
        host = mags[mge.host]
        assert list(mge_bins[mges[mge.name]].hosts) == [host]
        assert mges[mge.name] in mag_bins[host].mges
    hosts = {mags[mge.host] for mge in metagenome.mges}
    assert all(not mag.mges for name, mag in mag_bins.items() if name not in hosts)
    table = pd.read_csv(tmp_path / "interactions.txt", sep="\t")
    assert list(table.columns) == ["mgeMAG", "MAG", "Total Signal", "% of signal", "Interacting Contigs"]
    assert len(table) == len(metagenome.mges)
    assert (table["% of signal"] > 80).all()
    assert (tmp_path / "histogram.png").stat().st_size > 0


def test_no_host_below_the_thresholds(binned, tmp_path):
    data, network, _, _ = binned
    mag_bins, mge_bins = mth.create_bins(data)
    mth.compute_mge_mag_interactions(
        network,
        mag_bins,
        mge_bins,
        min_interacting_contigs=1000,
        output_file=str(tmp_path / "interactions.txt"),
        image_file=str(tmp_path / "histogram.png"),
    )
    assert all(not mge.hosts for mge in mge_bins.values())
    assert len(pd.read_csv(tmp_path / "interactions.txt", sep="\t")) == 0
    assert not (tmp_path / "histogram.png").exists()


def test_classify_mags():
    names = [f"m{i}" for i in range(7)]
    mags = {name: mth.Mag(name) for name in [*names, "absent"]}
    summary = _bin_summary(
        names,
        [np.nan, 0.95, 0.95, 0.95, 0.8, 0.6, 0.3],
        [1.0, 1.2, 1.02, 1.08, 1.05, 1.05, 1.0],
    )
    mags = mth.classify_mags(mags, summary)
    qualities = [mags[name].quality for name in names]
    assert qualities == ["Unknown", "Contaminated", "Complete", "HQ", "MQ", "LQ", "PQ"]
    assert mags["m4"].completeness == 0.8 and mags["m4"].contamination == 1.05
    assert mags["m4"].info["MAG"] == "m4"
    assert mags["absent"].quality is None


def test_estimate_noise():
    """Mean intra- and inter-MAG signals, the latter padded to one per other MAG."""
    mags = {name: mth.Mag(name) for name in ["A", "B", "C"]}
    for name, contigs in {"A": [1, 2], "B": [3, 4], "C": [5]}.items():
        for contig in contigs:
            mags[name].add_contig(f"contig_{contig}", contig)
    network = pd.DataFrame(
        {"contig1": [1, 3, 1, 2], "contig2": [2, 4, 3, 4], "signal": [10.0, 8.0, 1.0, 3.0]},
    )
    mags = mth.estimate_noise(mags, network)
    assert mags["A"].intra_signal == 10.0 and mags["B"].intra_signal == 8.0
    assert mags["C"].intra_signal == 0.0
    assert mags["A"].inter_signals == [2.0, 0]
    assert mags["B"].inter_signals == [2.0, 0]
    assert mags["C"].inter_signals == 0


def test_background_noise(metagenome, binned):
    """Contigs of a MAG contact each other more than other MAGs' contigs."""
    data, network, _, _ = binned
    mag_bins, _ = mth.create_bins(data)
    mag_bins = mth.estimate_noise(mag_bins, network)
    for mag in mag_bins.values():
        assert len(mag.inter_signals) == len(mag_bins) - 1
        assert mag.intra_signal > 10 * max(mag.inter_signals)


def _half_violins(ax):
    """What the split violins of ax show, by (x, side): the range and the
    quartiles of each half violin, and the value of each line drawn in
    place of a half violin."""
    violins, quartiles, lines = {}, {}, {}
    for body in [c for c in ax.collections if isinstance(c, PolyCollection)]:
        extents = body.get_paths()[0].get_extents()
        assert extents.width == pytest.approx(0.4)
        side = "high" if np.isclose(extents.x0, round(extents.x0)) else "low"
        x = round(extents.x0 if side == "high" else extents.x1)
        colour = {"low": "C0", "high": "C1"}[side]
        assert to_hex(body.get_facecolor()[0]) == to_hex(colour)
        violins[x, side] = (extents.y0, extents.y1)
    for collection in [c for c in ax.collections if isinstance(c, LineCollection)]:
        for (x0, y), (x1, _) in collection.get_segments():
            # Quartile lines span half of the half violin from its axis.
            side = "high" if np.isclose(x0, round(x0)) else "low"
            quartiles.setdefault((round(x0 if side == "high" else x1), side), []).append(y)
    for line in ax.lines:
        (x0, x1), (y, _) = line.get_xdata(), line.get_ydata()
        side = "high" if np.isclose(x0, round(x0)) else "low"
        lines[round(x0 if side == "high" else x1), side] = y
    return violins, {key: sorted(values) for key, values in quartiles.items()}, lines


def test_plot_mag_noise(binned, tmp_path, saved_figures):
    """One split violin per quality: intra-MAG signals on the left half,
    inter-MAG signals on the right half, with their quartiles."""
    data, network, _, _ = binned
    mag_bins, _ = mth.create_bins(data)
    mag_bins = mth.estimate_noise(mag_bins, network)
    qualities = ["HQ", "HQ", "MQ", "Contaminated", "MQ"]
    for mag, quality in zip(mag_bins.values(), qualities):
        mag.quality = quality
    mth.plot_mag_noise(mag_bins, image_file=str(tmp_path / "noise.png"))
    assert (tmp_path / "noise.png").stat().st_size > 0
    assert not plt.get_fignums()  # closed once saved

    ax = saved_figures["noise.png"].axes[0]
    order = ["HQ", "MQ", "Contaminated"]  # order of first appearance
    assert [label.get_text() for label in ax.get_xticklabels()] == order
    assert list(ax.get_xticks()) == [0, 1, 2]
    assert ax.get_xlabel() == "MAG Quality"
    assert ax.get_ylabel() == "Interaction Signal (log scale)"
    assert ax.get_title() == "Distribution of Intra- and Inter-MAG Interaction Signals"
    legend = ax.get_legend()
    assert legend.get_title().get_text() == "Interaction Type"
    assert [text.get_text() for text in legend.get_texts()] == ["intra", "inter"]
    assert [to_hex(handle.get_facecolor()) for handle in legend.legend_handles] == [to_hex("C0"), to_hex("C1")]

    signals = {}
    for mag, quality in zip(mag_bins.values(), qualities):
        signals.setdefault((order.index(quality), "low"), []).append(mag.intra_signal)
        signals.setdefault((order.index(quality), "high"), []).extend(mag.inter_signals)
    signals = {key: np.log10(np.array(values) + 1e-10) for key, values in signals.items()}
    violins, quartiles, lines = _half_violins(ax)
    # The only Contaminated MAG has a single intra-MAG signal: a line.
    assert lines == {(2, "low"): pytest.approx(signals[2, "low"][0])}
    del signals[2, "low"]
    assert violins.keys() == signals.keys() == quartiles.keys()
    for key, values in signals.items():
        assert violins[key] == pytest.approx((values.min(), values.max()))
        assert quartiles[key] == pytest.approx(np.percentile(values, [25, 50, 75]))


def test_plot_mag_noise_leaves_out_missing_values(tmp_path, saved_figures):
    """MAGs without a quality are left out, and so are missing signals."""
    mags = {}
    for name, quality, intra, inter in [
        ("m1", "HQ", 1.0, [0.01, 0.02]),
        ("m2", None, 5.0, [0.5, 0.6]),
        ("m3", "HQ", 2.0, [0.03, 0.01]),
        ("m4", "LQ", 3.0, []),
    ]:
        mags[name] = mth.Mag(name)
        mags[name].quality, mags[name].intra_signal, mags[name].inter_signals = quality, intra, inter
    mth.plot_mag_noise(mags, image_file=str(tmp_path / "noise.png"))
    ax = saved_figures["noise.png"].axes[0]
    assert [label.get_text() for label in ax.get_xticklabels()] == ["HQ", "LQ"]
    violins, quartiles, lines = _half_violins(ax)

    def log(*signals):
        return np.log10(np.array(signals) + 1e-10)

    assert violins == {(0, "low"): pytest.approx(log(1.0, 2.0)), (0, "high"): pytest.approx(log(0.01, 0.03))}
    assert quartiles.keys() == violins.keys()
    # LQ: one intra-MAG signal and no inter-MAG signal.
    assert lines == {(1, "low"): pytest.approx(log(3.0)[0])}


def test_half_violins_without_density():
    """Groups with one distinct value are lines, empty groups are skipped."""
    _, ax = plt.subplots()
    mth._half_violins(ax, [np.array([1.0]), np.array([]), np.array([2.0, 2.0])], "high", "C1")
    assert not ax.collections
    lines = sorted((tuple(line.get_xdata()), tuple(line.get_ydata())) for line in ax.lines)
    assert lines == [((0, 0.4), (1.0, 1.0)), ((2, 2.4), (2.0, 2.0))]


def test_annotate_hosts(metagenome, binned, tmp_path):
    data, network, mags, mges = binned
    summary = _bin_summary(list(mags.values()), [0.95, 0.95, 0.8, 0.6, 0.3], [1.0, 1.2, 1.05, 1.05, 1.0])
    mag_bins, mge_bins = mth.annotate_hosts(data, network, summary, 10, 5)
    for mge in metagenome.mges:
        assert list(mge_bins[mges[mge.name]].hosts) == [mags[mge.host]]
    qualities = {name: mag.quality for name, mag in mag_bins.items()}
    assert qualities == dict(zip(mags.values(), ["Complete", "Contaminated", "MQ", "LQ", "PQ"]))
    for name in ["background_noise.png", "mge_mag_histogram.png", "mge_mag_interactions.txt"]:
        assert (tmp_path / name).stat().st_size > 0
