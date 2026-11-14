# Test for quality_check module

import shutil

import metator.commands as mtc
import metator.io as mio
import metator.quality_check as mtq
import pandas as pd
import pytest
from tests import synthetic

CONTIGS = ["NODE_522", "NODE_1404", "NODE_1814", "NODE_2100", "NODE_2210", "NODE_2398"]


@pytest.fixture
def large_contigs_pairs():
    """Pairs between the large contigs, extracted for the quality metrics."""
    contigs_data = mio.read_contig_data("tests_data/outdir/contig_data_network.txt")
    out_file = "alignment_large.pairs"
    n_pairs = mtq.extract_pairs(["tests_data/outdir/alignment.pairs"], out_file, CONTIGS, contigs_data)
    return out_file, n_pairs


def test_extract_pairs(large_contigs_pairs):
    _, n_pairs = large_contigs_pairs
    assert n_pairs == 1310


def test_hic_quality(large_contigs_pairs):
    pairs_file, _ = large_contigs_pairs
    contigs = dict.fromkeys(["NODE_522", "NODE_1814", "NODE_2398"], "MetaTOR_00001_00001")
    contigs |= dict.fromkeys(["NODE_1404", "NODE_2100", "NODE_2210"], "MetaTOR_00001_00002")
    (
        n_religated,
        n_loops,
        n_weirds,
        n_informative_intra,
        n_informative_inter,
        n_intra_mags,
        n_inter_mags,
    ) = mtq.hic_quality(
        contigs,
        2,
        pairs_file,
        pairs_file + ".idx",
        "tests_data/assembly.fa",
        ["DpnII", "HpaII"],
        prefix="test",
        plot=True,
        plot_event="event.pdf",
        plot_cam="cam.pdf",
        threshold=["10", "10"],
    )
    assert n_religated == 1
    assert n_loops == 1
    assert n_weirds == 394
    assert n_informative_intra == 882
    assert n_informative_inter == 15
    assert n_intra_mags == 1293
    assert n_inter_mags == 17


# Quality check on a mock metagenome (tests/synthetic.py) of two genomes
# assembled in contigs larger than 100 kb, the size the quality check uses.


@pytest.fixture(scope="module")
def large_contigs(tmp_path_factory):
    """Assembly, sorted pairs and metator outputs (contig data, bin summary)
    of two genomes in large contigs; a third, low-quality bin has none."""
    if shutil.which("pairtools") is None:
        pytest.skip("pairtools is not installed")
    genomes = (
        synthetic.GenomeSpec("big_1", gc=0.40, abundance=1.0, length=260_000),
        synthetic.GenomeSpec("big_2", gc=0.60, abundance=2.0, length=240_000),
    )
    community = synthetic.make_metagenome(
        seed=11, genomes=genomes, mges=(), n_pairs=30_000, contig_median=120_000, contig_bounds=(105_000, 150_000)
    )
    directory = tmp_path_factory.mktemp("large_contigs")
    files = community.write_all(directory)
    pairs = mio.sort_pairs_pairtools(str(files["pairs"]), tmp_dir=str(directory))
    truth = community.truth
    bins = {"big_1": "metator_00001_00000", "big_2": "metator_00002_00000"}
    contig_data = pd.DataFrame(
        {"ID": range(1, len(truth) + 1), "Name": truth.Name, "Size": truth.Length, "Final_bin": truth.Genome.map(bins)}
    )
    contig_data.to_csv(directory / "contig_data_final.txt", sep="\t", index=False)
    summary = pd.DataFrame(
        {"Weighted completeness": [0.95, 0.92, 0.5], "Weighted redundancy": [1.01, 1.05, 1.0]},
        index=[*bins.values(), "metator_00003_00000"],
    )
    summary.to_csv(directory / "bin_summary.txt", sep="\t")
    return community, files, pairs, directory


def test_extract_hq_contigs(large_contigs):
    directory = large_contigs[3]
    contigs_data = mio.read_contig_data(str(directory / "contig_data_final.txt"))
    summary = mio.read_bin_summary(str(directory / "bin_summary.txt"))
    hq_contigs, n_mags = mtq.extract_hq_contigs(summary, contigs_data)
    assert n_mags == 2
    assert hq_contigs == dict(zip(contigs_data.Name, contigs_data.Final_bin))
    # A contaminated MAG is not high quality.
    summary.loc["metator_00002_00000", "Weighted redundancy"] = 1.3
    hq_contigs, n_mags = mtq.extract_hq_contigs(summary, contigs_data)
    assert n_mags == 1
    assert set(hq_contigs.values()) == {"metator_00001_00000"}


def test_hic_quality_counts(large_contigs, tmp_path):
    """Every extracted pair is counted once, in the right category."""
    community, files, pairs, directory = large_contigs
    contigs_data = mio.read_contig_data(str(directory / "contig_data_final.txt"))
    names = list(contigs_data.Name)
    n_pairs = mtq.extract_pairs([str(pairs)], str(tmp_path / "hq.pairs"), names, contigs_data)
    assert n_pairs == len(community.pairs)
    contigs = dict(zip(names, contigs_data.Final_bin))
    counts = mtq.hic_quality(
        contigs,
        2,
        str(tmp_path / "hq.pairs"),
        str(tmp_path / "hq_idx.pairs"),
        str(files["assembly"]),
        ["DpnII"],
        threshold=["3", "3"],
    )
    religated, loops, weirds, informative_intra, informative_inter, intra_mags, inter_mags = counts
    simulated = community.pairs
    genome = community.truth.Genome
    same_contig = simulated.chr1 == simulated.chr2
    same_genome = simulated.chr1.map(genome) == simulated.chr2.map(genome)
    assert intra_mags + inter_mags == n_pairs
    assert inter_mags == (~same_genome).sum()
    assert informative_inter == (same_genome & ~same_contig).sum()
    assert religated + loops + weirds + informative_intra == same_contig.sum()
    assert religated > 0 and loops > 0 and weirds > 0
    # With a single MAG there is no noise score to compute, the counts stay.
    assert (
        mtq.hic_quality(
            contigs,
            1,
            str(tmp_path / "hq.pairs"),
            str(tmp_path / "idx2.pairs"),
            str(files["assembly"]),
            ["DpnII"],
            threshold=["3", "3"],
        )
        == counts
    )


def test_qc_command(large_contigs, tmp_path):
    """From the metator output directory, with automatic thresholds and plots."""
    _, files, pairs, directory = large_contigs
    out = tmp_path / "qc"
    args = ["-a", str(files["assembly"]), "-e", "DpnII", "-O", str(directory), "-o", str(out), "-p", "sample", "-P",
            "-T", str(tmp_path / "tmp"), str(pairs)]  # fmt: skip
    mtc.Qc(args, {}).execute()
    assert (out / "sample_camembert_plot.pdf").stat().st_size > 0
    assert (out / "sample_event.pdf").stat().st_size > 0
    assert not (tmp_path / "tmp").exists()


def test_qc_command_without_hq_mags(large_contigs, tmp_path):
    """Explicit tables; no high-quality MAG, so nothing is plotted."""
    _, files, pairs, directory = large_contigs
    summary = mio.read_bin_summary(str(directory / "bin_summary.txt"))
    summary["Weighted completeness"] = 0.4
    summary.to_csv(tmp_path / "bin_summary.txt", sep="\t")
    args = ["-a", str(files["assembly"]), "-e", "DpnII", "-b", str(tmp_path / "bin_summary.txt"),
            "-c", str(directory / "contig_data_final.txt"), "-t", "3,3", "-T", str(tmp_path / "tmp"), "-n",
            str(pairs)]  # fmt: skip
    mtc.Qc(args, {}).execute()
    assert not list((tmp_path / "metator_qc").glob("*.pdf"))
    assert (tmp_path / "tmp").exists()  # --no-clean-up


@pytest.mark.parametrize("missing", ["bin_summary.txt", "contig_data_final.txt"])
def test_qc_command_missing_tables(large_contigs, tmp_path, missing):
    _, files, pairs, directory = large_contigs
    metator_dir = tmp_path / "metator"
    metator_dir.mkdir()
    for name in {"bin_summary.txt", "contig_data_final.txt"} - {missing}:
        shutil.copyfile(directory / name, metator_dir / name)
    args = ["-a", str(files["assembly"]), "-e", "DpnII", "-O", str(metator_dir), "-T", str(tmp_path / "tmp"), str(pairs)]
    with pytest.raises(FileNotFoundError):
        mtc.Qc(args, {}).execute()
