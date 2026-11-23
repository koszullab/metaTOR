# Test for network module

import gzip
import metator.io as mio
import numpy as np
import metator.network as mtn
import pytest
import os
import shutil
from pathlib import Path


assembly = "tests_data/assembly.fa"
depth_file = "tests_data/depth.txt"
contigs_data = {
    "NODE_1": {
        "id": 1,
        "Name": "NODE_1",
        "length": 642311,
        "GC": 38.6876450815882,
        "hit": 30795,
        "coverage": 41.1565,
        "RS": 0,
    },
    "NODE_6": {
        "id": 6,
        "Name": "NODE_6",
        "length": 505338,
        "GC": 53.0124391991103,
        "hit": 5533054,
        "coverage": 1330.17,
        "RS": 2505,
    },
    "NODE_8": {
        "id": 8,
        "Name": "NODE_8",
        "length": 416773,
        "GC": 26.977035460550468,
        "hit": 206333,
        "coverage": 97.6355,
        "RS": 3496,
    },
    "NODE_10": {
        "id": 10,
        "Name": "NODE_10",
        "length": 396325,
        "GC": 36.0124391991103,
        "hit": 306312,
        "coverage": 0,
        "RS": 2578,
    },
}


def test_create_contig_data():
    # Test contig data builder.
    # Case without depth file.
    contig_data = mtn.create_contig_data(assembly)
    assert len(contig_data) == 1219
    assert contig_data["NODE_522"] == {
        "id": 1,
        "length": 22786,
        "GC": 62.0775915035548,
        "hit": 0,
        "coverage": "-",
        "RS": "-",
    }
    # Case with depth file.
    contig_data = mtn.create_contig_data(assembly, depth_file, "DpnII")
    assert len(contig_data) == 1219
    assert contig_data["NODE_522"] == {
        "id": 1,
        "length": 22786,
        "GC": 62.0775915035548,
        "hit": 0,
        "coverage": 4.76595,
        "RS": 162,
    }


def test_normalize_pair():
    # Test for edges normalisation methods.
    length = mtn.normalize_pair(contigs_data, ("NODE_1", "NODE_6"), 52, "length")
    abundance = mtn.normalize_pair(contigs_data, ("NODE_1", "NODE_8"), 2004, "abundance")
    abundance_null = mtn.normalize_pair(contigs_data, ("NODE_1", "NODE_10"), 275, "abundance")
    RS = mtn.normalize_pair(contigs_data, ("NODE_6", "NODE_8"), 800, "RS")
    EH = mtn.normalize_pair(contigs_data, ("NODE_8", "NODE_1"), 2022, "empirical_hit")
    TH = mtn.normalize_pair(contigs_data, ("NODE_6", "NODE_8"), 3333, "theoritical_hit")
    TH_null = mtn.normalize_pair(contigs_data, ("NODE_1", "NODE_8"), 2004, "theoritical_hit")
    assert length == pytest.approx(71.77, abs=1e-2)
    assert abundance == pytest.approx(31.61, abs=1e-2)
    assert abundance_null == 0
    assert RS == pytest.approx(221.57, abs=1e-2)
    assert EH == pytest.approx(2.54e-2, abs=1e-4)
    assert TH == pytest.approx(2.51e-05, abs=1e-7)
    assert TH_null == 0


# Contig data with a contig that has no hit and one that has no coverage, for
# which normalize_pair returns 0. The restriction sites are counted plus one,
# as create_contig_data counts them: never 0.
NORMALIZED_CONTIGS = {
    "a": {"id": 1, "length": 5000, "hit": 120, "coverage": 12.5, "RS": 31},
    "b": {"id": 2, "length": 1200, "hit": 7, "coverage": 3.0, "RS": 6},
    "c": {"id": 3, "length": 800, "hit": 0, "coverage": 4.5, "RS": 2},
    "d": {"id": 4, "length": 2500, "hit": 40, "coverage": 0, "RS": 14},
}


@pytest.mark.parametrize("normalization", ["abundance", "length", "RS", "empirical_hit", "theoritical_hit"])
def test_normalization_weights_factorise_normalize_pair(normalization):
    """One weight per contig, whose products times a contact count are the
    contacts that normalize_pair normalizes: NaN where it gives 0, for a
    contig without hits or without coverage."""
    names = list(NORMALIZED_CONTIGS)
    values = {
        key: np.array([NORMALIZED_CONTIGS[name][key] for name in names], dtype=float)
        for key in ("length", "hit", "coverage", "RS")
    }
    weights = mtn.normalization_weights(values["length"], values["hit"], values["coverage"], values["RS"], normalization)
    for i, first in enumerate(names):
        for j, second in enumerate(names):
            if i == j:
                continue
            expected = mtn.normalize_pair(NORMALIZED_CONTIGS, (first, second), 2004, normalization)
            normalized = 2004 * weights[i] * weights[j]
            assert (0.0 if np.isnan(normalized) else normalized) == pytest.approx(expected, rel=1e-12, abs=0)
    assert np.isnan(weights).any()


def test_normalization_weights_refusals():
    lengths = np.array([10.0, 20.0])
    with pytest.raises(ValueError, match="No weights for the normalization 'None'"):
        mtn.normalization_weights(lengths, lengths, lengths, lengths, "None")
    with pytest.raises(ValueError, match='"abundance" normalization uses the shotgun coverage'):
        mtn.normalization_weights(lengths, lengths, None, None, "abundance")


def test_create_contig_data_reads_the_depth_file_by_name(tmp_path):
    """The depth file is read by contig name: in another order, it gives the
    same lengths and coverages. It was read line by line, in the order of
    the assembly, and gave each contig the coverage of another."""
    header, *lines = Path(depth_file).read_text().splitlines(keepends=True)
    reversed_depth = tmp_path / "depth.txt"
    reversed_depth.write_text(header + "".join(reversed(lines)))
    expected = mtn.create_contig_data(assembly, depth_file, "DpnII")
    found = mtn.create_contig_data(assembly, str(reversed_depth), "DpnII")
    assert found == expected
    assert found["NODE_522"]["coverage"] == 4.76595


@pytest.mark.parametrize(
    "case, message",
    [
        ("missing", "1 contigs of the assembly are not in the depth file"),
        ("other length", "1 contigs of the assembly have another length in the depth file"),
    ],
)
def test_create_contig_data_refuses_the_depth_file_of_another_assembly(case, message, tmp_path, caplog):
    header, first, *lines = Path(depth_file).read_text().splitlines(keepends=True)
    if case == "other length":
        name, length, *rest = first.split("\t")
        lines.append("\t".join([name, str(int(length) + 1), *rest]))
    other = tmp_path / "depth.txt"
    other.write_text(header + "".join(lines))
    with caplog.at_level("INFO", logger="metator_logger"), pytest.raises(ValueError, match="is not that of the assembly"):
        mtn.create_contig_data(assembly, str(other))
    assert f"{message} {other}, such as NODE_522" in caplog.text


def test_create_contig_data_from_a_compressed_assembly(tmp_path):
    """A gzipped assembly gives the same contig data as the uncompressed one."""
    compressed = tmp_path / "assembly.fa.gz"
    with open(assembly, "rb") as fasta, gzip.open(compressed, "wb") as gz:
        shutil.copyfileobj(fasta, gz)
    assert mtn.create_contig_data(str(compressed), depth_file, "DpnII") == mtn.create_contig_data(assembly, depth_file, "DpnII")


def test_write_contig_data():
    # Test contig data writer.
    tmp_file = "tmp_contig_data_test.tsv"
    mtn.write_contig_data(contigs_data, tmp_file)
    data = mio.read_contig_data(tmp_file)
    assert data.loc["NODE_1", "Size"] == 642311
    assert data.loc["NODE_8", "Hit"] == 206333
    assert data.loc["NODE_6", "Restriction_site"] == 2505
    os.remove(tmp_file)
