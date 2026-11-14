# Test for quality_check module

import gzip
import locale
import re

import Bio
import metator.io as mio
import metator.quality_check as mtq
import pytest
from metator.log import logger
from packaging.version import Version

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


@pytest.mark.skipif(
    Version(Bio.__version__) < Version("1.87"),
    reason="Biopython < 1.87 skips text before the first record (with a warning from 1.85).",
)
@pytest.mark.parametrize("before", ["\n", "# assembly\n", "not a fasta\n"])
def test_hic_quality_rejects_text_before_the_first_record(tmp_path, caplog, monkeypatch, before):
    """metator qc names the assembly instead of failing in Biopython's parser."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(before + ">NODE_522\nGATCGATC\n")
    monkeypatch.setattr(logger, "propagate", True)
    with pytest.raises(ValueError, match=f"^Not a FASTA file: {re.escape(str(fasta))}$"):
        mtq.hic_quality({}, 0, "absent.pairs", "absent_idx.pairs", str(fasta), ["DpnII"])
    assert f"{fasta} is not a FASTA file: it must start with a '>' header line" in caplog.text


@pytest.mark.skipif(
    locale.getpreferredencoding(False).lower().replace("-", "") != "utf8",
    reason="Needs a UTF-8 locale to make the byte below undecodable.",
)
def test_parse_fasta_passes_decoding_errors_on(tmp_path, caplog, monkeypatch):
    """A byte that is not UTF-8 is not reported as a missing '>' header."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_bytes(b">NODE_522\nGATC\xe9GATC\n")
    monkeypatch.setattr(logger, "propagate", True)
    with pytest.raises(UnicodeDecodeError):
        list(mtq.parse_fasta(str(fasta)))
    assert "is not a FASTA file" not in caplog.text


def test_parse_fasta_reads_compressed_files(tmp_path):
    fasta = tmp_path / "assembly.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">NODE_522\nGATC\n>NODE_1404\nGG\n")
    assert [(record.id, str(record.seq)) for record in mtq.parse_fasta(str(fasta))] == [("NODE_522", "GATC"), ("NODE_1404", "GG")]
