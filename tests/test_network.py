# Test for network module

import gzip
import io
import pandas as pd
import numpy as np
import metator.network as mtn
import pytest
from tests.oracles import normalize_pair
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
    # Test contig data builder: what it measures on the assembly, by contig
    # name, in the order of the assembly. The hits are counted on the map,
    # and the IDs are positions in the map: neither is a key.
    # Case without depth file.
    contig_data = mtn.create_contig_data(assembly)
    assert len(contig_data) == 1219
    assert next(iter(contig_data)) == "NODE_522"
    assert contig_data["NODE_522"] == {"length": 22786, "GC": 62.0775915035548, "coverage": "-", "RS": "-"}
    # Case with depth file.
    contig_data = mtn.create_contig_data(assembly, depth_file, "DpnII")
    assert len(contig_data) == 1219
    assert contig_data["NODE_522"] == {"length": 22786, "GC": 62.0775915035548, "coverage": 4.76595, "RS": 162}


def _normalized(contigs, pair, n_occ, normalization):
    """The normalized contact of a pair of contigs (by name, in contigs: their
    data, as create_contig_data gives it) through the weights of
    normalization_weights: n_occ * w[i] * w[j], 0 for a contig without a weight
    (NaN), whose contacts are no edge."""
    names = list(contigs)
    values = {key: np.array([contigs[name][key] for name in names], dtype=float) for key in ("length", "hit", "coverage", "RS")}
    weights = mtn.normalization_weights(values["length"], values["hit"], values["coverage"], values["RS"], normalization)
    normalized = n_occ * weights[names.index(pair[0])] * weights[names.index(pair[1])]
    return 0.0 if np.isnan(normalized) else normalized


def test_normalized_contacts():
    # The normalizations of the contacts of metaTOR 1.x (its normalize_pair),
    # through the weights of the contigs.
    length = _normalized(contigs_data, ("NODE_1", "NODE_6"), 52, "length")
    abundance = _normalized(contigs_data, ("NODE_1", "NODE_8"), 2004, "abundance")
    abundance_null = _normalized(contigs_data, ("NODE_1", "NODE_10"), 275, "abundance")
    RS = _normalized(contigs_data, ("NODE_6", "NODE_8"), 800, "RS")
    EH = _normalized(contigs_data, ("NODE_8", "NODE_1"), 2022, "empirical_hit")
    TH = _normalized(contigs_data, ("NODE_6", "NODE_8"), 3333, "theoritical_hit")
    TH_null = _normalized(contigs_data, ("NODE_1", "NODE_8"), 2004, "theoritical_hit")
    assert length == pytest.approx(71.77, abs=1e-2)
    assert abundance == pytest.approx(31.61, abs=1e-2)
    assert abundance_null == 0
    assert RS == pytest.approx(221.57, abs=1e-2)
    assert EH == pytest.approx(2.54e-2, abs=1e-4)
    assert TH == pytest.approx(2.51e-05, abs=1e-7)
    assert TH_null == 0


# Contig data with a contig that has no hit and one that has no coverage, for
# which the normalization of metaTOR 1.x gives 0. The restriction sites are
# counted plus one, as create_contig_data counts them: never 0.
NORMALIZED_CONTIGS = {
    "a": {"id": 1, "length": 5000, "hit": 120, "coverage": 12.5, "RS": 31},
    "b": {"id": 2, "length": 1200, "hit": 7, "coverage": 3.0, "RS": 6},
    "c": {"id": 3, "length": 800, "hit": 0, "coverage": 4.5, "RS": 2},
    "d": {"id": 4, "length": 2500, "hit": 40, "coverage": 0, "RS": 14},
}


@pytest.mark.parametrize("normalization", ["abundance", "length", "RS", "empirical_hit", "theoritical_hit"])
def test_normalization_weights_factorise_the_normalizations_of_1x(normalization):
    """One weight per contig, whose products times a contact count are the
    contacts that metaTOR 1.x normalized edge by edge: NaN where it gave 0, for
    a contig without hits or without coverage."""
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
            expected = normalize_pair(NORMALIZED_CONTIGS, (first, second), 2004, normalization)
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
    """The depth file is read by contig name: in another order, it gives the same
    lengths and coverages."""
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


def one_contig_assembly(tmp_path, length):
    """An assembly of one contig, "big", of `length` bases."""
    fasta = tmp_path / "big.fa"
    bases = np.random.default_rng(0).choice(list("ACGT"), size=length)
    fasta.write_text(">big\n" + "".join(bases) + "\n")
    return str(fasta)


def test_read_depth_reads_lengths_in_scientific_notation(tmp_path):
    """jgi_summarize_bam_contig_depths writes the lengths of a million bases or
    more in scientific notation (#26): read as floats, the others as ints."""
    depth = tmp_path / "depth.txt"
    depth.write_text("contigName\tcontigLen\ttotalAvgDepth\nbig\t1.21754e+06\t3.5\nNODE_522\t22786\t4.76595\n")
    depths = mtn.read_depth(str(depth))
    assert depths == {"big": (1217540.0, 3.5), "NODE_522": (22786, 4.76595)}
    assert isinstance(depths["big"][0], float)
    assert isinstance(depths["NODE_522"][0], int)


@pytest.mark.parametrize(
    "depth_length, length, matches",
    [
        (22786, 22786, True),
        (22787, 22786, False),
        (22786.0, 22786, True),
        # 6 significant digits: 1217543 is written 1.21754e+06.
        (1217540.0, 1217543, True),
        (1217543.0, 1217543, True),
        (1217550.0, 1217543, False),
        # Written as an integer, a length is the length.
        (1217540, 1217543, False),
    ],
)
def test_depth_length_matches(depth_length, length, matches):
    assert mtn.depth_length_matches(depth_length, length) is matches


def test_create_contig_data_with_a_length_in_scientific_notation(tmp_path):
    """A contig of 1,217,543 bases, whose length the depth file writes
    1.21754e+06, gets its coverage (#26: int() refused the length)."""
    fasta = one_contig_assembly(tmp_path, 1217543)
    depth = tmp_path / "depth.txt"
    depth.write_text("contigName\tcontigLen\ttotalAvgDepth\nbig\t1.21754e+06\t3.5\n")
    contig_data = mtn.create_contig_data(fasta, str(depth))
    assert contig_data["big"]["length"] == 1217543
    assert contig_data["big"]["coverage"] == 3.5


@pytest.mark.parametrize("written", ["1.21755e+06", "1217540"])
def test_create_contig_data_refuses_another_length_in_scientific_notation(written, tmp_path, caplog):
    """A length that is not that of the contig, rounded to 6 significant digits
    or written whole, is refused as another length."""
    fasta = one_contig_assembly(tmp_path, 1217543)
    depth = tmp_path / "depth.txt"
    depth.write_text(f"contigName\tcontigLen\ttotalAvgDepth\nbig\t{written}\t3.5\n")
    with caplog.at_level("INFO", logger="metator_logger"), pytest.raises(ValueError, match="is not that of the assembly"):
        mtn.create_contig_data(fasta, str(depth))
    assert f"1 contigs of the assembly have another length in the depth file {depth}, such as big" in caplog.text


def test_create_contig_data_from_a_compressed_assembly(tmp_path):
    """A gzipped assembly gives the same contig data as the uncompressed one."""
    compressed = tmp_path / "assembly.fa.gz"
    with open(assembly, "rb") as fasta, gzip.open(compressed, "wb") as gz:
        shutil.copyfileobj(fasta, gz)
    assert mtn.create_contig_data(str(compressed), depth_file, "DpnII") == mtn.create_contig_data(assembly, depth_file, "DpnII")


def test_contig_data_table():
    """The contig data of the contigs of a map, in its order: their IDs are their
    positions in the map, from 1, and what was measured on them their columns,
    as contig_data_network.txt held them; a column that the run did not measure
    is "-"."""
    names = ["NODE_8", "NODE_1", "NODE_6"]
    table = mtn.contig_data_table(contigs_data, names, [206333, 30795, 5533054])
    assert list(table.columns) == mtn.CONTIG_TABLE_COLUMNS
    assert list(table.ID) == [1, 2, 3] and list(table.Name) == names
    assert list(table.Size) == [416773, 642311, 505338] and list(table.Hit) == [206333, 30795, 5533054]
    assert list(table.Restriction_site) == [3496, 0, 2505] and table.Shotgun_coverage[1] == 41.1565
    numeric = [pd.api.types.is_numeric_dtype(dtype) for dtype in table.dtypes]
    assert numeric == [True, False, True, True, True, True, True]
    unmeasured = {name: {**data, "coverage": "-", "RS": "-"} for name, data in contigs_data.items()}
    table = mtn.contig_data_table(unmeasured, names, [0, 0, 0])
    assert set(table.Shotgun_coverage) == set(table.Restriction_site) == {"-"}
    # As partition writes it, and reads it back.
    assert pd.read_csv(io.StringIO(table.to_csv(sep="\t", index=False)), sep="\t").equals(table)


@pytest.mark.parametrize("normalization", ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"])
def test_network_weights(normalization):
    """The weights of a contig data table: those of normalization_weights, none
    for None."""
    names = list(contigs_data)
    table = mtn.contig_data_table(contigs_data, names, [contigs_data[name]["hit"] for name in names])
    weights = mtn.network_weights(table, normalization)
    if normalization == "None":
        assert weights is None
        return
    values = {
        key: np.array([contigs_data[name][key] for name in names], dtype=float) for key in ("length", "hit", "coverage", "RS")
    }
    expected = mtn.normalization_weights(values["length"], values["hit"], values["coverage"], values["RS"], normalization)
    np.testing.assert_array_equal(weights, expected)


def test_network_weights_without_a_column():
    """A normalization that needs a column that the run did not measure."""
    unmeasured = {name: {**data, "coverage": "-"} for name, data in contigs_data.items()}
    table = mtn.contig_data_table(unmeasured, list(contigs_data), [1, 2, 3, 4])
    with pytest.raises(ValueError, match='"abundance" normalization uses the shotgun coverage'):
        mtn.network_weights(table, "abundance")
