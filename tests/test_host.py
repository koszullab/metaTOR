# Tests of the MGE host detection (metator.host) on the mock metagenome:
# each MGE contacts one host genome much more than the others.

import sys

import numpy as np
import pandas as pd
import pytest

import metator.commands as mtc
import metator.host as mth

needs_py310 = pytest.mark.xfail(
    sys.version_info < (3, 10),
    raises=TypeError,
    strict=True,
    reason="host.classify_mags annotates a nested function with Literal[...] | Literal[...], which needs Python 3.10",
)


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


@needs_py310
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


def test_plot_mag_noise(binned, tmp_path):
    data, network, _, _ = binned
    mag_bins, _ = mth.create_bins(data)
    mag_bins = mth.estimate_noise(mag_bins, network)
    for mag, quality in zip(mag_bins.values(), ["HQ", "HQ", "MQ", "Contaminated", "MQ"]):
        mag.quality = quality
    mth.plot_mag_noise(mag_bins, image_file=str(tmp_path / "noise.png"))
    assert (tmp_path / "noise.png").stat().st_size > 0


@needs_py310
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


@pytest.mark.xfail(
    strict=True,
    raises=AttributeError,
    reason="metator host calls metator.host.host_detection, which no longer exists",
)
def test_host_command(metagenome, binned, metagenome_network, tmp_path):
    data, _, _, _ = binned
    data.to_csv(tmp_path / "contig_data.txt", sep="\t", index=False)
    data[["Name", "Final_bin"]].to_csv(tmp_path / "binning.txt", sep="\t", index=False, header=False)
    mge_contigs = [c.name for mge in metagenome.mges for c in mge.contigs]
    (tmp_path / "mges.txt").write_text("\n".join(mge_contigs) + "\n")
    args = [
        "-n", str(metagenome_network["network"]),
        "-b", str(tmp_path / "binning.txt"),
        "-m", str(tmp_path / "mges.txt"),
        "-c", str(tmp_path / "contig_data.txt"),
        "-o", str(tmp_path / "hosts.tsv"),
    ]  # fmt: skip
    mtc.Host(args, {}).execute()
    assert (tmp_path / "hosts.tsv").stat().st_size > 0
