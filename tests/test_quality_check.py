# Test for quality_check module

import metator.io as mio
import metator.quality_check as mtq
import pytest

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
