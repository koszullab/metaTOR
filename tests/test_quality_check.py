# Test for validation module

import metator.io as mio
import metator.quality_check as mtq
import pytest
import os


# def test_extract_hq_contigs():
#     bin_summary = mio.read_bin_summary("tests_data/outdir/bin_summary.txt")
#     contigs_data = mio.read_contig_data(
#         "tests_data/outdir/contig_data_final.txt"
#     )
#     hq_contigs = mtq.extract_hq_contigs(bin_summary, contigs_data)
#     assert hq_contigs == {"NODE_3_length_540571_cov_9.303447": "MetaTOR_3_1"}


def test_extract_pairs():
    contigs_data = mio.read_contig_data(
        "tests_data/outdir/contig_data_network.txt"
    )
    pairs_files = ["tests_data/outdir/alignment.pairs"]
    out_file = "tests_data/outdir/alignment_large.pairs"
    contigs = [
        "NODE_522",
        "NODE_1404",
        "NODE_1814",
        "NODE_2100",
        "NODE_2210",
        "NODE_2398",
    ]
    n_pairs = mtq.extract_pairs(pairs_files, out_file, contigs, contigs_data)
    assert n_pairs == 1310


def test_hic_quality():
    contigs = {
        "NODE_522": "MetaTOR_00001_00001",
        "NODE_1814": "MetaTOR_00001_00001",
        "NODE_2398": "MetaTOR_00001_00001",
        "NODE_1404": "MetaTOR_00001_00002",
        "NODE_2100": "MetaTOR_00001_00002",
        "NODE_2210": "MetaTOR_00001_00002",
    }
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
        "tests_data/outdir/alignment_large.pairs",
        "tests_data/outdir/alignment_large.pairs.idx",
        "tests_data/assembly.fa",
        ["DpnII", "HpaII"],
        prefix="test",
        plot=True,
        plot_event="tests_data/outdir/event.pdf",
        plot_cam="tests_data/outdir/cam.pdf",
        threshold=["10", "10"],
    )
    assert n_religated == 1
    assert n_loops == 1
    assert n_weirds == 394
    assert n_informative_intra == 882
    assert n_informative_inter == 15
    assert n_intra_mags == 1293
    assert n_inter_mags == 17
    os.remove("tests_data/outdir/alignment_large.pairs")
    os.remove("tests_data/outdir/alignment_large.pairs.idx")


def test_quality_check():
    ...
