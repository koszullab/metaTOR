# Tests of metator.metacool: the contact map of a binning as a .metacool
# file, on the mock metagenome of tests/synthetic.py. Its genome contigs are
# binned by genome and its MGE contigs left unbinned.

import bz2
import gzip
import logging
import lzma
import re
import shlex
import shutil
import tarfile
import zipfile
from pathlib import Path

import cooler
import metacooler
import numpy as np
import pandas as pd
import pytest

import metator.io as mio
import metator.log as mtl
import metator.metacool as mtmc
import metator.network as mtn
from metator.version import __version__

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"

PAIRS_HEADER = "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"


def _bin_names(metagenome):
    """Final bin of each contig: one bin per genome, none for the MGEs."""
    names = {genome.name: f"metator_{i:05d}_00000" for i, genome in enumerate(metagenome.genomes, start=1)}
    return {contig.name: names.get(contig.molecule, "ND") for contig in metagenome.contigs}


# The contig columns of the map: those of contig_data_final.txt, as
# write_contig_data writes it, but Restriction_site, which is - everywhere.
CONTIG_COLUMNS = ["GC_content", "Hit", "Shotgun_coverage"]


def write_contig_data(metagenome, path, final_bin=None):
    """A contig_data_final.txt file, laid out as metaTOR writes it (an
    unnamed index first, Final_bin last), of a run with a depth file and no
    enzyme: Restriction_site is - on every contig."""
    final_bin = final_bin or _bin_names(metagenome)
    truth = metagenome.truth
    table = pd.DataFrame(
        {
            "ID": range(1, len(truth) + 1),
            "Name": truth.Name.to_numpy(),
            "Size": truth.Length.to_numpy(),
            "GC_content": truth.GC.to_numpy(),
            "Hit": range(10, len(truth) + 10),
            "Shotgun_coverage": (truth.Abundance * 3.5).to_numpy(),
            "Restriction_site": "-",
            "Overlapping_bin_ID": "00001",
            "Recursive_bin_contigs": "-",
            "Final_bin": truth.Name.map(final_bin).to_numpy(),
        }
    )
    table.to_csv(path, sep="\t", index=True)
    return path


# The QC counts of each contig, which metacooler (0.3.0) stores with the
# QC block.
QC_CONTIG_COLUMNS = ["qc_intra_contig_short", "qc_intra_contig_long"]


def without_qc(columns):
    """The columns of a contigs table but metacooler's QC columns."""
    return [column for column in columns if not column.startswith("qc_")]


def assert_contig_columns(path, contig_data, columns=CONTIG_COLUMNS):
    """The contigs table of the map (cooler's chroms) holds these columns of
    the contig data, with their values and types, and no other: no bin, no
    pandas index. Besides, metacooler's QC counts of each contig, when the
    map has a QC block."""
    contigs = metacooler.MetaCooler(str(path)).contigs
    assert without_qc(contigs.columns) == ["name", "length", "genome_id", *columns]
    has_qc = "qc" in cooler.Cooler(str(path)).info["metadata"]
    assert [column for column in contigs.columns if column.startswith("qc_")] == (QC_CONTIG_COLUMNS if has_qc else [])
    chroms = cooler.Cooler(str(path)).chroms()[:]
    assert list(chroms.columns) == list(contigs.columns)
    table = pd.read_csv(contig_data, sep="\t", index_col="Name").loc[contigs.name]
    for column in columns:
        assert contigs[column].dtype == table[column].dtype, column
        assert contigs[column].tolist() == table[column].tolist(), column
        assert chroms[column].tolist() == table[column].tolist(), column


def write_bin_summary(metagenome, path):
    """A bin_summary.txt file, written as validation writes it."""
    final_bin = _bin_names(metagenome)
    summary = {}
    for i, genome in enumerate(metagenome.genomes):
        summary[final_bin[genome.contigs[0].name]] = {
            "Length": float(genome.length),
            "GC-content": 40.0 + i,
            "Present Markers": 100.0 + i,
            "Completeness": 0.9 - i / 100,
            "Weighted redundancy": 1.0 + i / 100,
            "recursive": False,
            "step": 0,
            "parent": None if i else "metator_00009_00000",
            "HiC_abundance": 10.0 + i / 3,
            "SG_abundance": 20.0 + i / 7,
        }
    mio.write_bin_summary(summary, path)
    return path


def write_pairs(pairs, path):
    with open(path, "w") as out:
        out.write(PAIRS_HEADER)
        pairs.to_csv(out, sep="\t", header=False, index=False)
    return path


@pytest.fixture
def binning(metagenome, metagenome_files, tmp_path):
    """Pairs, contig data and bin summary of a binning of the mock community."""
    return {
        "pairs": metagenome_files["pairs"],
        "contig_data": write_contig_data(metagenome, tmp_path / "contig_data_final.txt"),
        "bin_summary": write_bin_summary(metagenome, tmp_path / "bin_summary.txt"),
    }


def expected_pixels(pairs, contig_order):
    """Contacts between contigs, counted from the pairs: one per pair, both
    ends on the contigs, in the upper triangle of the map."""
    index = {name: i for i, name in enumerate(contig_order)}
    keep = pairs.chr1.isin(index) & pairs.chr2.isin(index)
    bin1 = pairs.chr1[keep].map(index)
    bin2 = pairs.chr2[keep].map(index)
    pixels = pd.DataFrame({"bin1_id": np.minimum(bin1, bin2), "bin2_id": np.maximum(bin1, bin2)})
    return pixels.value_counts().sort_index()


def map_pixels(path):
    pixels = cooler.Cooler(str(path)).pixels()[:]
    return pixels.set_index(["bin1_id", "bin2_id"])["count"].sort_index()


def test_metacool_of_a_binning(metagenome, binning, tmp_path):
    out = tmp_path / "out" / "sample.metacool"
    mtmc.write_metacool(
        [str(binning["pairs"])],
        str(binning["contig_data"]),
        str(out),
        bin_summary=str(binning["bin_summary"]),
        assembly_name="synthetic",
        metadata={"command": "metator metacool", "parameters": {"threads": 1}},
    )
    metacooler.validate(str(out), full=True)

    # cooler reads it: one bin per contig, the whole contig.
    clr = cooler.Cooler(str(out))
    truth = metagenome.truth
    bins = clr.bins()[:]
    assert sorted(bins.chrom) == sorted(truth.Name)
    assert (bins.start == 0).all()
    assert (bins.end.to_numpy() == truth.Length[bins.chrom].to_numpy()).all()
    assert clr.info["genome-assembly"] == "synthetic"

    # One genome per bin, plus the unbinned contigs (the MGEs).
    mc = metacooler.MetaCooler(str(out))
    genomes = mc.genomes.set_index("name")
    summary = mio.read_bin_summary(binning["bin_summary"])
    assert list(genomes.index) == ["__unbinned__", *sorted(summary.index)]
    contigs = mc.contigs
    found = contigs.genome_id.map(mc.genomes.name)
    final_bin = _bin_names(metagenome)
    assert list(found) == [final_bin[name] if final_bin[name] != "ND" else "__unbinned__" for name in contigs.name]
    assert genomes.loc["__unbinned__", "n_contigs"] == sum(len(mge.contigs) for mge in metagenome.mges)
    # Grouped by genome, the longest contigs first.
    for _, group in contigs.groupby("genome_id"):
        assert list(group.length) == sorted(group.length, reverse=True)

    # The bin summary columns, one row per bin.
    for column in summary.columns:
        assert column in genomes.columns, column
    for column in ["Length", "GC-content", "Present Markers", "Completeness", "HiC_abundance", "SG_abundance"]:
        assert genomes.loc[summary.index, column].to_numpy() == pytest.approx(summary[column].to_numpy(), abs=1e-4)
    assert list(genomes.loc[summary.index, "parent"]) == ["metator_00009_00000"] + [""] * (len(summary) - 1)

    # Every pair is a contact between its two contigs.
    pixels = map_pixels(out)
    assert pixels.sum() == len(metagenome.pairs)
    expected = expected_pixels(metagenome.pairs, list(mc.contigs.name))
    assert pixels.to_dict() == expected.to_dict()

    # What metaTOR measured on each contig, in the contigs table.
    assert_contig_columns(out, binning["contig_data"])
    assert contigs.GC_content.dtype == "float64" and contigs.Hit.dtype == "int64"

    metadata = clr.info["metadata"]
    assert metadata["metator"] == {
        "version": __version__,
        "command": "metator metacool",
        "parameters": {"threads": 1},
    }
    assert metadata["annotate"][0]["added"] == [f"contigs.{column}" for column in CONTIG_COLUMNS]
    assert metadata["qc"]["n_pairs_read"] == len(metagenome.pairs)
    # Only the map: no P(s) table next to it.
    assert sorted(path.name for path in out.parent.iterdir()) == ["sample.metacool"]


def test_unbinned_contigs_keep_their_contacts(metagenome, metagenome_files, tmp_path):
    """Contigs with Final_bin ND, - or empty are in __unbinned__, with all
    their contacts; a map without bin summary or assembly name."""
    final_bin = _bin_names(metagenome)
    names = [contig.name for contig in metagenome.contigs]
    for name, value in zip(names[:3], ["ND", "-", ""], strict=True):
        final_bin[name] = value
    contig_data = write_contig_data(metagenome, tmp_path / "contig_data_final.txt", final_bin)
    out = tmp_path / "map.metacool"
    mtmc.write_metacool(str(metagenome_files["pairs"]), str(contig_data), str(out))

    mc = metacooler.MetaCooler(str(out))
    unbinned = set(mc.contigs.name[mc.contigs.genome_id == 0])
    assert mc.genomes.name[0] == "__unbinned__"
    assert unbinned == set(names[:3]) | {c.name for mge in metagenome.mges for c in mge.contigs}
    assert "-" not in set(mc.genomes.name) and "" not in set(mc.genomes.name)
    pixels = map_pixels(out)
    assert pixels.to_dict() == expected_pixels(metagenome.pairs, list(mc.contigs.name)).to_dict()
    # Every contact of an unbinned contig is there.
    pairs = metagenome.pairs
    touching = pairs.chr1.isin(unbinned) | pairs.chr2.isin(unbinned)
    ids = set(np.flatnonzero(mc.contigs.name.isin(unbinned)))
    on_unbinned = pixels[[b1 in ids or b2 in ids for b1, b2 in pixels.index]]
    assert on_unbinned.sum() == touching.sum() > 0
    assert cooler.Cooler(str(out)).info["genome-assembly"] == "unknown"


def test_libraries_are_summed(metagenome, binning, metagenome_indexed_pairs, tmp_path):
    """Two libraries, one plain and one sorted and bgzipped by metaTOR: their
    contacts add up."""
    pairs = metagenome.pairs
    half = len(pairs) // 2
    first = write_pairs(pairs.iloc[:half], tmp_path / "library_1.pairs")
    second = write_pairs(pairs.iloc[half:], tmp_path / "library_2.pairs")
    second = mio.sort_and_index_pairs(str(second), tmp_dir=str(tmp_path))
    out = tmp_path / "two.metacool"
    mtmc.write_metacool([str(first), second], str(binning["contig_data"]), str(out), threads=2)
    whole = tmp_path / "whole.metacool"
    mtmc.write_metacool([str(metagenome_indexed_pairs)], str(binning["contig_data"]), str(whole))
    assert map_pixels(out).sum() == len(pairs)
    assert map_pixels(out).to_dict() == map_pixels(whole).to_dict()


def test_pairs_off_the_contigs_are_left_out(metagenome, binning, tmp_path):
    """A pair with an end on a contig the contig data does not list, or past
    the end of its contig, is no contact."""
    pairs = metagenome.pairs.iloc[:100].copy()
    first, last = pairs.chr1.iloc[0], pairs.chr1.iloc[1]
    extra = pd.DataFrame(
        [
            ["r_unknown", "contig_unknown", 10, first, 10, "+", "-"],
            ["r_past_end", first, 10, last, metagenome.truth.Length[last] + 5, "+", "-"],
        ],
        columns=pairs.columns,
    )
    path = write_pairs(pd.concat([pairs, extra]), tmp_path / "extra.pairs")
    out = tmp_path / "extra.metacool"
    mtmc.write_metacool(str(path), str(binning["contig_data"]), str(out))
    pixels = map_pixels(out)
    assert pixels.sum() == 100
    assert pixels.to_dict() == expected_pixels(pairs, list(metacooler.MetaCooler(str(out)).contigs.name)).to_dict()


def test_missing_input_files(binning, tmp_path):
    out = tmp_path / "map.metacool"
    with pytest.raises(mtmc.MetacoolError, match=re.escape("not found: missing.pairs, missing_summary.txt")):
        mtmc.write_metacool(
            [str(binning["pairs"]), "missing.pairs"], str(binning["contig_data"]), str(out), bin_summary="missing_summary.txt"
        )
    assert not out.exists()


@pytest.mark.parametrize("kind", ["contig", "bin"])
def test_non_ascii_names_are_refused(metagenome, metagenome_files, kind, tmp_path, monkeypatch):
    """Before the pairs are read."""
    final_bin = _bin_names(metagenome)
    names = list(final_bin)
    if kind == "bin":
        final_bin[names[0]] = final_bin[names[1]] = "bin_é"
    contig_data = write_contig_data(metagenome, tmp_path / "contig_data_final.txt", final_bin)
    if kind == "contig":
        text = contig_data.read_text()
        for name, new in [(names[0], "contig_é"), (names[1], "contig_ü")]:
            text = text.replace(f"\t{name}\t", f"\t{new}\t")
        contig_data.write_text(text)
    monkeypatch.setattr(metacooler, "from_pairs", lambda *args, **kwargs: pytest.fail("pairs read"))
    expected = "'contig_é', 'contig_ü'" if kind == "contig" else "'bin_é'."
    with pytest.raises(mtmc.MetacoolError, match=f"needs ASCII {kind} names .* are not: {expected}"):
        mtmc.write_metacool(str(metagenome_files["pairs"]), str(contig_data), str(tmp_path / "map.metacool"))


def test_non_ascii_names_shown():
    """The first names, then how many more: the message ends with one period."""
    names = r"'é0', 'é1', 'é2', 'é3', 'é4', \.\.\. \(\+2 more\)"
    with pytest.raises(mtmc.MetacoolError, match=rf"7 contig name\(s\) are not: {names}\.$"):
        mtmc.check_ascii_names(["ok", *[f"é{i}" for i in range(7)], "é0"], [])
    mtmc.check_ascii_names(["contig_1"], ["metator_00001_00000"])


NO_BIN_NAMES = ["NA", "nan", "N/A", "nd", "N/d", "-1", "__unbinned__", "__UNBINNED__"]


def test_bin_names_read_as_no_bin_are_refused(metagenome, metagenome_files, tmp_path, monkeypatch):
    """Bins that metacooler would silently empty into __unbinned__: refused,
    named, before the pairs are read. ND is metaTOR's own mark of a contig
    in no bin, not a bin name."""
    final_bin = _bin_names(metagenome)
    names = list(final_bin)
    for name, bin_name in zip(names, [*NO_BIN_NAMES, "NA", "ND"], strict=False):
        final_bin[name] = bin_name
    contig_data = write_contig_data(metagenome, tmp_path / "contig_data_final.txt", final_bin)
    monkeypatch.setattr(metacooler, "from_pairs", lambda *args, **kwargs: pytest.fail("pairs read"))
    shown = r"'NA', 'nan', 'N/A', 'nd', 'N/d', \.\.\. \(\+3 more\)"
    message = (
        "A .metacool file cannot have bins named as metacooler marks a contig in no bin (NA, NaN, N/A, N/D, -1 or "
        "__unbinned__ in any case, and nd in any case other than ND, which marks a contig in no bin): it would put "
        "their contigs in __unbinned__. Rename the 8 bin(s): "
    )
    with pytest.raises(mtmc.MetacoolError, match=rf"^{re.escape(message)}{shown}\.$"):
        mtmc.write_metacool(str(metagenome_files["pairs"]), str(contig_data), str(tmp_path / "map.metacool"))


def test_no_bin_names_are_those_of_metacooler():
    """The bin names refused are exactly the genome names that metacooler
    reads as no genome (-1 among them from metacooler 0.2.0)."""
    from metacooler.io.create import _UNBINNED_SENTINELS

    assert set(_UNBINNED_SENTINELS) == mtmc.METACOOLER_UNBINNED
    mtmc.check_bin_names(["metator_00001_00000", "MetaTOR_MGE_00001", "ND_1", "na1"])
    for name in mtmc.METACOOLER_UNBINNED:
        with pytest.raises(mtmc.MetacoolError, match="Rename the 2 bin"):
            mtmc.check_bin_names(["metator_00001_00000", name.upper(), f" {name} "])


def test_contig_columns_the_table_has(metagenome, binning, tmp_path):
    """A contig data file with some of the contig columns: those it has; one
    with none of them: none, and the map all the same."""
    table = pd.read_csv(binning["contig_data"], sep="\t", index_col=0)
    some = tmp_path / "some_columns.txt"
    table.drop(columns=["GC_content", "Shotgun_coverage"]).to_csv(some, sep="\t")
    out = tmp_path / "some.metacool"
    mtmc.write_metacool(str(binning["pairs"]), str(some), str(out))
    metacooler.validate(str(out), full=True)
    assert_contig_columns(out, some, ["Hit"])

    none = tmp_path / "no_columns.txt"
    table[["Name", "Size", "Final_bin"]].to_csv(none, sep="\t")
    out = tmp_path / "none.metacool"
    mtmc.write_metacool(str(binning["pairs"]), str(none), str(out))
    metacooler.validate(str(out), full=True)
    assert_contig_columns(out, none, [])
    assert "annotate" not in cooler.Cooler(str(out)).info["metadata"]


@pytest.mark.parametrize("force", [False, True])
def test_contig_column_metacooler_refuses(force, binning, tmp_path):
    """A contig column that cooler could not read (not ASCII): metacooler
    refuses it, and no map is left, rather than a map without its contig
    columns."""
    contig_data = Path(binning["contig_data"])
    lines = contig_data.read_text().splitlines(keepends=True)
    fields = lines[1].split("\t")
    fields[4] = "38.5é"
    contig_data.write_text("".join([lines[0], "\t".join(fields), *lines[2:]]))
    out = tmp_path / "map.metacool"
    if force:
        out.write_text("old")
    with pytest.raises(mtmc.MetacoolError, match="GC_content") as error:
        mtmc.write_metacool(str(binning["pairs"]), str(contig_data), str(out), force=force)
    assert isinstance(error.value.__cause__, ValueError)
    assert not out.exists()


@pytest.mark.parametrize("suffix", [".gz", ".bz2"])
def test_maps_of_a_compressed_contig_data(suffix, metagenome, binning, tmp_path):
    """A contig data file compressed with gzip as .gz: both maps, with its
    contig columns, as metacooler reads it too (this case passed before the
    refusal of other compressions, which it guards). Compressed with bzip2:
    refused before any work, rather than once the map is written, when
    metacooler failed to read it."""
    contig_data = tmp_path / f"contig_data_final.txt{suffix}"
    opener = {".gz": gzip.open, ".bz2": bz2.open}[suffix]
    with opener(contig_data, "wb") as handle:
        handle.write(Path(binning["contig_data"]).read_bytes())
    from_pairs = tmp_path / "pairs.metacool"
    from_map = tmp_path / "map.metacool"
    unbinned = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    if suffix == ".bz2":
        with pytest.raises(mtmc.MetacoolError, match="compressed in a way that metacooler"):
            mtmc.write_metacool(str(binning["pairs"]), str(contig_data), str(from_pairs))
        with pytest.raises(mtmc.MetacoolError, match="compressed in a way that metacooler"):
            mtmc.reassign_metacool(str(unbinned), str(contig_data), str(from_map))
        assert not from_pairs.exists() and not from_map.exists()
        return
    mtmc.write_metacool(str(binning["pairs"]), str(contig_data), str(from_pairs))
    mtmc.reassign_metacool(str(unbinned), str(contig_data), str(from_map))
    for out in [from_pairs, from_map]:
        metacooler.validate(str(out), full=True)
        assert_contig_columns(out, binning["contig_data"])
        assert assignment(out) == final_bins(metagenome)


def test_read_contig_data(tmp_path):
    """The bins and final bins of the contigs, the contig columns the table
    has but those that are - on every contig (no depth file, no enzyme),
    and all its columns; not the bins, nor the other columns."""
    table = tmp_path / "contig_data_final.txt"
    table.write_text(
        "\tID\tName\tSize\tGC_content\tShotgun_coverage\tRestriction_site\tHit\tFinal_bin\n"
        "0\t1\t01\t100\t40.5\t-\t-\t3\tb1\n"
        "1\t2\tc2\t200\t41.5\t-\t12\t4\tND\n"
        "2\t3\tc3\t300\t41.5\t-\t12\t4\t-\n"
        "3\t4\tc4\t400\t41.5\t-\t12\t4\t\n"
    )
    bins, final_bin, columns, header = mtmc.read_contig_data(str(table))
    assert bins.to_dict("list") == {"chrom": ["01", "c2", "c3", "c4"], "start": [0] * 4, "end": [100, 200, 300, 400]}
    assert final_bin == {"01": "b1"}
    assert columns == ["GC_content", "Hit", "Restriction_site"]
    assert header == ["Unnamed: 0", "ID", "Name", "Size", "GC_content", "Shotgun_coverage", "Restriction_site", "Hit",
                      "Final_bin"]  # fmt: skip
    table.write_text("Name\tSize\tShotgun_coverage\tFinal_bin\nc1\t100\t - \tb1\n")
    assert mtmc.read_contig_data(str(table))[2] == []
    table.write_text("Name\tSize\tFinal_bin\nc1\t100\tb1\n")
    assert mtmc.read_contig_data(str(table))[2] == []


@pytest.mark.parametrize(
    "content, reason",
    [
        (b"", "No columns to parse from file."),
        (
            b"Name\tSize\tFinal_bin\nc1\t100\tb1\nc2\t200\tb1\tb2\n",
            ": 1 row(s) do not have the 3 tab-separated fields of the header: line 3 (4).",
        ),
        (
            b"Name\tSize\tFinal_bin\n0\tc1\t100\tb1\n1\tc2\t200\tND\n",
            ": 2 row(s) do not have the 3 tab-separated fields of the header: line 2 (4), line 3 (4).",
        ),
        (
            b"\nName\tSize\tFinal_bin\n\nc1\t100\tb1\n  \nc2\t200\n\t\t\t\n",
            ": 2 row(s) do not have the 3 tab-separated fields of the header: line 6 (2), line 7 (4).",
        ),
        (b"Name\tSize\tFinal_bin\nc\xe91\t100\tb1\n", "'utf-8' codec can't decode byte 0xe9"),
    ],
    ids=["empty", "row of more fields", "rows of more fields from the first", "rows of other fields", "not UTF-8"],
)
def test_read_contig_data_unreadable_table(content, reason, tmp_path):
    """A contig data file that pandas cannot read as a table, or would read
    shifted or padded (a row with more or fewer fields than the header,
    blank lines left out and counted): a MetacoolError naming the file,
    the cause chained."""
    table = tmp_path / "contig_data_final.txt"
    table.write_bytes(content)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(str(table))} cannot be read as a table: ") as error:
        mtmc.read_contig_data(str(table))
    assert reason in str(error.value)
    assert isinstance(error.value.__cause__, ValueError)


def test_read_contig_data_of_a_gzipped_file(tmp_path):
    """A contig data file compressed with gzip is read, and its rows
    checked, as pandas reads it."""
    table = tmp_path / "contig_data_final.txt.gz"
    with gzip.open(table, "wt") as handle:
        handle.write("Name\tSize\tFinal_bin\nc1\t100\tb1\nc2\t200\tND\n")
    bins, final_bin, *_ = mtmc.read_contig_data(str(table))
    assert bins.to_dict("list") == {"chrom": ["c1", "c2"], "start": [0, 0], "end": [100, 200]}
    assert final_bin == {"c1": "b1"}
    with gzip.open(table, "wt") as handle:
        handle.write("Name\tSize\tFinal_bin\nc1\t100\tb1\tb2\n")
    message = r"1 row\(s\) do not have the 3 tab-separated fields of the header: line 2 \(4\)\.$"
    with pytest.raises(mtmc.MetacoolError, match=message):
        mtmc.read_contig_data(str(table))


@pytest.mark.parametrize("suffix", [".bz2", ".xz", ".zip", ".tar.gz", ".GZ"])
def test_read_contig_data_compressed_otherwise(suffix, tmp_path):
    """A contig data file that pandas reads compressed otherwise than with
    gzip as .gz is refused: metacooler, which reads it again to store the
    contig columns, cannot (the map failed there, once written)."""
    content = "Name\tSize\tFinal_bin\nc1\t100\tb1\n"
    table = tmp_path / f"contig_data_final.txt{suffix}"
    source = tmp_path / "contig_data_final.txt"
    source.write_text(content)
    if suffix == ".zip":
        with zipfile.ZipFile(table, "w") as archive:
            archive.write(source, arcname=source.name)
    elif suffix == ".tar.gz":
        with tarfile.open(table, "w:gz") as archive:
            archive.add(source, arcname=source.name)
    else:
        opener = {".bz2": bz2.open, ".xz": lzma.open, ".GZ": gzip.open}[suffix]
        with opener(table, "wt") as handle:
            handle.write(content)
    assert pd.read_csv(table, sep="\t", dtype=str).to_dict("list") == {
        "Name": ["c1"],
        "Size": ["100"],
        "Final_bin": ["b1"],
    }
    message = (
        f"{table} is compressed in a way that metacooler, which reads it again to store the contig columns, cannot "
        "read: give it plain, or compressed with gzip as .gz."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        mtmc.read_contig_data(str(table))


def test_read_contig_data_header_with_a_trailing_tab(tmp_path):
    """A header that ends in a tab names an empty last column: the rows may
    leave it out, as pandas and metacooler read them. A row with more fields
    than the header is still refused."""
    table = tmp_path / "contig_data_final.txt"
    table.write_text("Name\tSize\tFinal_bin\t\nc1\t100\tb1\nc2\t200\tND\t\n")
    bins, final_bin, *_ = mtmc.read_contig_data(str(table))
    assert list(bins.chrom) == ["c1", "c2"] and final_bin == {"c1": "b1"}
    table.write_text("Name\tSize\tFinal_bin\t\nc1\t100\tb1\t\tx\n")
    message = r"1 row\(s\) do not have the 3 to 4 tab-separated fields of the header: line 2 \(5\)\.$"
    with pytest.raises(mtmc.MetacoolError, match=message):
        mtmc.read_contig_data(str(table))


@pytest.mark.parametrize("size", ["1x00", "", "100.0"])
def test_read_contig_data_size_not_a_whole_number(size, tmp_path):
    table = tmp_path / "contig_data_final.txt"
    table.write_text(f"Name\tSize\tFinal_bin\nc1\t100\tb1\nc2\t{size}\tND\n")
    message = f"The Size column of {table} must give the length of each contig in bp, as a whole number: "
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}.*'{re.escape(size)}'\\.$") as error:
        mtmc.read_contig_data(str(table))
    assert isinstance(error.value.__cause__, ValueError)


def test_outfile_in_a_file(binning, tmp_path):
    """An outfile whose directory cannot be made: a MetacoolError, with the
    message of the OSError."""
    (tmp_path / "file").write_text("not a directory")
    with pytest.raises(mtmc.MetacoolError, match=r"^\[Errno 17\] File exists: ") as error:
        mtmc.write_metacool(str(binning["pairs"]), str(binning["contig_data"]), str(tmp_path / "file" / "map.metacool"))
    assert isinstance(error.value.__cause__, FileExistsError)


def test_metacooler_messages_go_to_the_log_file(metagenome, binning, tmp_path):
    """metacooler's messages reach the log file of the run; its logger is
    left as it was."""
    log = tmp_path / "run.log"
    mc_logger = logging.getLogger("metacooler")
    before = (mc_logger.level, mc_logger.propagate, list(mc_logger.handlers))
    mtl.set_file_handler(str(log))
    try:
        mtmc.write_metacool(str(binning["pairs"]), str(binning["contig_data"]), str(tmp_path / "map.metacool"))
    finally:
        for handler in mtl.logger.handlers[:]:
            if isinstance(handler, logging.FileHandler):
                mtl.logger.removeHandler(handler)
                handler.close()
    logged = log.read_text()
    n_binned = sum(len(genome.contigs) for genome in metagenome.genomes)
    n_bins = len(metagenome.genomes)
    assert f"contact map of {len(metagenome.contigs)} contigs, {n_binned} of them in {n_bins} bins." in logged
    assert "INFO :: Created metacool file" in logged
    assert (mc_logger.level, mc_logger.propagate, list(mc_logger.handlers)) == before


def test_one_pairs_file_as_a_path(binning, tmp_path):
    """One pairs file, given as a path; it is only read."""
    pairs = tmp_path / "copy.pairs"
    shutil.copyfile(binning["pairs"], pairs)
    mtmc.write_metacool(Path(pairs), str(binning["contig_data"]), str(tmp_path / "map.metacool"))
    assert pairs.read_bytes() == Path(binning["pairs"]).read_bytes()


# The metacool command.


def test_is_map(binning, tmp_path):
    """A map is an HDF5 file, or the URI of a map inside one (file::group),
    whatever the name: a pairs file may have :: in its name, and a URI is
    no map when its file is no HDF5 file (or no file)."""
    metacool = tmp_path / "map.metacool"
    mtmc.write_metacool(str(binning["pairs"]), str(binning["contig_data"]), str(metacool))
    colons = shutil.copyfile(binning["pairs"], tmp_path / "lib::L1.pairs")
    named = shutil.copyfile(metacool, tmp_path / "lib::L1.metacool")
    assert [mtmc.is_map(path) for path in [metacool, named, f"{metacool}::/", f"{metacool}::entries/a"]] == [True] * 4
    pairs = binning["pairs"]
    others = [pairs, colons, f"{pairs}::/", tmp_path / "missing.metacool", f"{tmp_path}/missing::x", tmp_path]
    assert [mtmc.is_map(path) for path in others] == [False] * 6


def test_metacool_command_from_a_pairs_file_named_with_colons(metagenome, binning, tmp_path, run_cli):
    """A pairs file whose name holds :: is a pairs file, not the URI of a
    map inside a file: the map is written from its pairs."""
    pairs = shutil.copyfile(binning["pairs"], tmp_path / "lib::L1.pairs")
    out = tmp_path / "out.metacool"
    run_cli("metacool", "-c", binning["contig_data"], "-b", binning["bin_summary"], "-o", out, pairs)
    metacooler.validate(str(out), full=True)
    assert map_pixels(out).sum() == len(metagenome.pairs)
    assert "reassign" not in cooler.Cooler(str(out)).info["metadata"]


def test_metacool_command(metagenome, binning, tmp_path, run_cli):
    pairs = metagenome.pairs
    half = len(pairs) // 2
    first = write_pairs(pairs.iloc[:half], tmp_path / "library_1.pairs")
    second = write_pairs(pairs.iloc[half:], tmp_path / "library_2.pairs")
    args = ["-c", binning["contig_data"], "-b", binning["bin_summary"], "-o", "maps/sample.metacool", "-n", "synthetic",
            "-t", "2", first, second]  # fmt: skip
    run_cli("metacool", *args)

    out = tmp_path / "maps" / "sample.metacool"
    metacooler.validate(str(out))
    clr = cooler.Cooler(str(out))
    assert clr.info["genome-assembly"] == "synthetic"
    assert clr.pixels()[:]["count"].sum() == len(pairs)
    genomes = metacooler.MetaCooler(str(out)).genomes
    assert list(genomes.name) == ["__unbinned__", *sorted(mio.read_bin_summary(binning["bin_summary"]).index)]
    assert "Completeness" in genomes.columns
    assert_contig_columns(out, binning["contig_data"])
    assert clr.info["metadata"]["metator"] == {
        "version": __version__,
        "command": shlex.join(["metator", "metacool", *map(str, args)]),
        "parameters": {
            "pairsfile": [str(first), str(second)],
            "contig_data": str(binning["contig_data"]),
            "bin_summary": str(binning["bin_summary"]),
            "outfile": "maps/sample.metacool",
            "assembly_name": "synthetic",
            "threads": 2,
            "force": False,
        },
    }
    # The log of the run is in the working directory.
    (log,) = tmp_path.glob("metator_metacool_*.log")
    logged = log.read_text()
    assert "## Command used: metator metacool -c" in logged
    assert "Created metacool file" in logged


@pytest.mark.parametrize("case", ["map and pairs", "assembly name"])
def test_metacool_command_from_a_map_refused(case, metagenome, binning, tmp_path, run_cli):
    """From a .metacool map: alone (a map sums its libraries already), and
    without --assembly-name (the map keeps that of the input map). Refused
    before anything is written, logged, then raised."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    args = {"map and pairs": [metacool, binning["pairs"]], "assembly name": ["-n", "assembly", metacool]}[case]
    message = {
        "map and pairs": "metator metacool writes the map from pairs files, or from one .metacool map: got 2 files, "
        "1 of them .metacool maps.",
        "assembly name": f"--assembly-name names the assembly of a map written from pairs files: the map written "
        f"from {metacool} keeps the assembly name of that map.",
    }[case]
    out = tmp_path / "map.metacool"
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("metacool", "-c", binning["contig_data"], "-o", out, *args)
    assert f"ERROR :: {message}" in _logged(tmp_path)
    assert not out.exists()


def test_metacool_command_from_an_mmetacool_entry(metagenome, binning, tmp_path, run_cli):
    """From an entry of a .mmetacool collection (file::entries/NAME): the
    map of the bins, with the contacts of the entry."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    out = tmp_path / "bins.metacool"
    run_cli("metacool", "-c", binning["contig_data"], "-b", binning["bin_summary"], "-o", out, f"{collection}::entries/lib1")

    metacooler.validate(str(out), full=True)
    assert assignment(out) == final_bins(metagenome)
    assert contig_contacts(out) == contig_contacts(metacool)


@pytest.mark.parametrize("kind", ["collection", "cooler"])
def test_metacool_command_from_another_hdf5_file(kind, metagenome, binning, tmp_path, run_cli):
    """From an HDF5 file that is no .metacool map: a .mmetacool collection,
    whose entries the error lists, or a .cool file. Taken for a map, not
    read as a pairs file, and refused with metacooler's reason, logged,
    then raised."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    if kind == "collection":
        path = tmp_path / "maps.mmetacool"
        metacooler.create_mmetacool(path, {"lib1": metacool})
        reason = "entries/lib1"
    else:
        path = tmp_path / "contigs.cool"
        clr = cooler.Cooler(str(metacool))
        cooler.create_cooler(str(path), clr.bins()[:][["chrom", "start", "end"]], clr.pixels()[:])
        reason = "metacooler-format"
    out = tmp_path / "map.metacool"
    with pytest.raises(mtmc.MetacoolError, match=f"^Cannot read the contact map {re.escape(str(path))}: .*{reason}") as error:
        run_cli("metacool", "-c", binning["contig_data"], "-o", out, path)
    assert f"ERROR :: {error.value}" in _logged(tmp_path)
    assert not out.exists()


def test_metacool_command_without_final_bin(metagenome, binning, tmp_path, run_cli):
    """A run with no final bin: the map of every contig, all in __unbinned__
    with all their contacts. The bin summary has no row and adds no
    column."""
    contig_data = write_contig_data(metagenome, tmp_path / "no_bin.txt", dict.fromkeys(_bin_names(metagenome), "ND"))
    summary = tmp_path / "no_bin_summary.txt"
    summary.write_text(Path(binning["bin_summary"]).read_text().splitlines(keepends=True)[0])
    run_cli("metacool", "-c", contig_data, "-b", summary, "-o", "map.metacool", binning["pairs"])

    out = tmp_path / "map.metacool"
    metacooler.validate(str(out), full=True)
    mc = metacooler.MetaCooler(str(out))
    assert list(mc.genomes.name) == ["__unbinned__"]
    assert mc.genomes.n_contigs[0] == len(metagenome.contigs)
    assert set(mio.read_bin_summary(binning["bin_summary"]).columns) & set(mc.genomes.columns) == set()
    assert map_pixels(out).to_dict() == expected_pixels(metagenome.pairs, list(mc.contigs.name)).to_dict()
    assert map_pixels(out).sum() == len(metagenome.pairs)
    assert_contig_columns(out, contig_data)
    assert f"contact map of {len(metagenome.contigs)} contigs, 0 of them in 0 bins." in _logged(tmp_path)


def _logged(directory):
    return "".join(log.read_text() for log in directory.glob("metator_metacool_*.log"))


def test_metacool_command_existing_outfile(binning, tmp_path, run_cli):
    out = tmp_path / "map.metacool"
    out.write_text("old")
    args = ["-c", binning["contig_data"], "-o", out, binning["pairs"]]
    message = f"{out} already exists. Use --force to overwrite it."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("metacool", *args)
    assert f"ERROR :: {message}" in _logged(tmp_path)
    assert out.read_text() == "old"
    run_cli("metacool", *args, "-F")
    metacooler.validate(str(out))


def test_metacool_command_misaligned_bin_summary(binning, tmp_path, run_cli):
    """A bin summary whose header has more fields than its rows, as
    tests_data/outdir/bin_summary.txt (18 and 15): metacooler cannot tell
    which column each value belongs to. Its reason is logged, then raised
    as a MetacoolError, and no map is left."""
    summary = TESTS_DATA / "outdir" / "bin_summary.txt"
    with pytest.raises(mtmc.MetacoolError, match=r"^Genomes table is misaligned: .*bin_summary\.txt\n") as error:
        run_cli("metacool", "-c", binning["contig_data"], "-b", summary, "-o", "map.metacool", binning["pairs"])
    assert "header names 18 column(s), but no data row has more than 15" in str(error.value)
    assert isinstance(error.value.__cause__, ValueError)
    assert f"ERROR :: {error.value}\n" in _logged(tmp_path)
    assert not (tmp_path / "map.metacool").exists()


@pytest.mark.parametrize(
    "error",
    [KeyError("genome_name"), OSError("No space left on device"), ValueError("bad")],
    ids=["KeyError", "OSError", "ValueError"],
)
def test_metacooler_errors_are_metacool_errors(error, binning, tmp_path, monkeypatch, run_cli):
    """What metacooler raises is logged, then raised again as a
    MetacoolError, its message kept (not the repr of a KeyError)."""

    def from_pairs(*args, **kwargs):
        raise error

    monkeypatch.setattr(metacooler, "from_pairs", from_pairs)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(error.args[0])}$") as raised:
        run_cli("metacool", "-c", binning["contig_data"], "-o", "map.metacool", binning["pairs"])
    assert raised.value.__cause__ is error
    assert f"ERROR :: {error.args[0]}\n" in _logged(tmp_path)


ERRORS = {
    "missing pairs file": (["missing.pairs"], "Input file(s) not found: missing.pairs."),
    "contig data without Final_bin": ([], "has no Final_bin column"),
    "non-ASCII contig names": ([], "A .metacool file needs ASCII contig names"),
    "bin named NA": ([], "Rename the 1 bin(s): 'NA'."),
    "Size not a whole number": ([], "must give the length of each contig in bp, as a whole number: "),
}


@pytest.mark.parametrize("case", list(ERRORS))
def test_metacool_command_errors(case, metagenome, metagenome_network, binning, tmp_path, run_cli):
    """Refused before any pair is read, as metaTOR refuses an input: the
    reason is logged as an error, in the log file of the run, and the
    MetacoolError stops the command."""
    pairs, message = ERRORS[case]
    contig_data = binning["contig_data"]
    if case == "contig data without Final_bin":
        contig_data = metagenome_network["contig_data"]
    elif case == "non-ASCII contig names":
        name = metagenome.contigs[0].name
        contig_data.write_text(contig_data.read_text().replace(f"\t{name}\t", "\tcontig_é\t"))
    elif case == "bin named NA":
        contig_data = write_contig_data(
            metagenome, tmp_path / "na.txt", {**_bin_names(metagenome), metagenome.contigs[0].name: "NA"}
        )
    elif case == "Size not a whole number":
        name = metagenome.contigs[0].name
        contig_data.write_text(re.sub(rf"\t{name}\t\d+\t", f"\t{name}\t1x00\t", contig_data.read_text()))
    with pytest.raises(mtmc.MetacoolError, match=re.escape(message)) as error:
        run_cli("metacool", "-c", contig_data, "-o", "map.metacool", *(pairs or [binning["pairs"]]))
    assert f"ERROR :: {error.value}\n" in _logged(tmp_path)
    assert not (tmp_path / "map.metacool").exists()


def test_metacool_command_threads_below_one(binning, tmp_path, monkeypatch, run_cli):
    """-t 0: refused with a message (not metacooler's ValueError), before
    any pair is read."""
    monkeypatch.setattr(metacooler, "from_pairs", lambda *args, **kwargs: pytest.fail("pairs read"))
    message = "The number of threads must be at least 1, got 0."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("metacool", "-c", binning["contig_data"], "-o", "map.metacool", "-t", "0", binning["pairs"])
    assert f"ERROR :: {message}" in _logged(tmp_path)
    assert not (tmp_path / "map.metacool").exists()


# The map of the bins written back from a .metacool map of the contigs
# (reassign_metacool), as a run that starts from a map writes it.


def unbinned_map(path, metagenome, pairs, binsize=None, qc=True):
    """A .metacool map of the pairs made before any binning, as metacooler
    create pairs without --contigs makes it: every contig in __unbinned__,
    in the order of the assembly (the #chromsize lines of the pairs), one
    bin per contig without binsize."""
    truth = metagenome.truth
    bins = (
        {"binsize": binsize}
        if binsize
        else {"bins": pd.DataFrame({"chrom": truth.Name.to_numpy(), "start": 0, "end": truth.Length.to_numpy()})}
    )
    metacooler.from_pairs(str(path), pairs=str(pairs), contig_genome_map=None, contig_ordering="fasta", ps=False, qc=qc, **bins)
    return Path(path)


def genomes_without_qc(path):
    genomes = metacooler.MetaCooler(str(path)).genomes
    return genomes[[column for column in genomes.columns if not column.startswith("qc_")]]


def assert_same_map(path, expected):
    """The two maps hold the same bins, contacts, contigs table and genomes
    table, but the QC columns of the genomes."""
    clr, ref = cooler.Cooler(str(path)), cooler.Cooler(str(expected))
    pd.testing.assert_frame_equal(clr.bins()[:], ref.bins()[:])
    pd.testing.assert_frame_equal(clr.pixels()[:], ref.pixels()[:])
    pd.testing.assert_frame_equal(clr.chroms()[:], ref.chroms()[:])
    pd.testing.assert_frame_equal(genomes_without_qc(path), genomes_without_qc(expected))


def contig_pairs(first, second):
    """Each pair of contigs, its two names in order."""
    first, second = np.asarray(first, dtype=object), np.asarray(second, dtype=object)
    return pd.DataFrame({"a": np.where(first <= second, first, second), "b": np.where(first <= second, second, first)})


def contig_contacts(path):
    """The contacts of a map summed by contig pair, whatever its bin size."""
    clr = cooler.Cooler(str(path))
    chrom = clr.bins()[:]["chrom"].astype(str).to_numpy(dtype=object)
    pixels = clr.pixels()[:]
    pairs = contig_pairs(chrom[pixels.bin1_id.to_numpy()], chrom[pixels.bin2_id.to_numpy()])
    pairs["count"] = pixels["count"].to_numpy()
    return pairs.groupby(["a", "b"])["count"].sum().to_dict()


def assignment(path):
    """The genome of each contig of a map, by contig name."""
    mc = metacooler.MetaCooler(str(path))
    return dict(zip(mc.contigs.name, mc.contigs.genome_id.map(mc.genomes.name), strict=True))


def final_bins(metagenome):
    return {name: (bin_name if bin_name != "ND" else "__unbinned__") for name, bin_name in _bin_names(metagenome).items()}


def test_map_of_the_bins_from_a_map(metagenome, binning, tmp_path):
    """From a map of the contigs made before any binning (one bin per
    contig, in the order of the assembly): the map that write_metacool
    writes from the pairs, QC included (metacooler 0.3.0 carries it over to
    the bins), but for its provenance. The input map is only read, and
    nothing but the map is left next to it."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    before = metacool.read_bytes()
    expected = tmp_path / "pairs.metacool"
    tables = {"bin_summary": str(binning["bin_summary"]), "metadata": {"command": "metator pipeline"}}
    mtmc.write_metacool(str(binning["pairs"]), str(binning["contig_data"]), str(expected), **tables)
    out = tmp_path / "out" / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(binning["contig_data"]), str(out), **tables)

    metacooler.validate(str(out), full=True)
    assert_same_map(out, expected)
    assert_contig_columns(out, binning["contig_data"])
    assert metacool.read_bytes() == before
    clr = cooler.Cooler(str(out))
    metadata = clr.info["metadata"]
    assert metadata["metator"] == {"version": __version__, "command": "metator pipeline"}
    assert metacooler.MetaCooler(str(out)).info["ordering"] == "genome-then-length"
    # The QC block and columns that write_metacool measures on the pairs with
    # the final bins: metacooler carries the intra-contig counts of the input
    # map over, read positions included, and measures the others again. The
    # assembly and the history are those of the input map.
    assert metadata["qc"]["source"] == "pairs" and metadata["qc"]["n_pairs_read"] == len(metagenome.pairs)
    assert metadata["qc"] == cooler.Cooler(str(expected)).info["metadata"]["qc"]
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(out)).genomes, metacooler.MetaCooler(str(expected)).genomes)
    assert [entry["source"] for entry in metadata["reassign"]] == [str(metacool)]
    assert clr.info["genome-assembly"] == cooler.Cooler(str(metacool)).info["genome-assembly"]
    assert sorted(path.name for path in out.parent.iterdir()) == ["bins.metacool"]
    # Read and written as the map of the pairs (not the private temporary file).
    assert out.stat().st_mode == expected.stat().st_mode


def test_map_of_the_bins_from_a_map_of_metacooler_0_2_0(metagenome, binning, tmp_path, logged):
    """From a map written by metacooler 0.2.0, which holds no QC count per
    contig (here, a map of 0.3.0 without its two columns): metacooler
    computes the QC block of the bins again from the pixels, and warns that
    what only the reads gave is lost. The map of the bins is written."""
    import h5py

    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    with h5py.File(metacool, "r+") as handle:
        for column in QC_CONTIG_COLUMNS:
            del handle["contigs"][column]
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    qc = cooler.Cooler(str(out)).info["metadata"]["qc"]
    assert qc["source"] == "pixels" and "n_pairs_read" not in qc
    (warning,) = [message for level, message in logged() if level == "WARNING"]
    assert warning.startswith("Recomputing the source's QC block from the pixels")
    assert "intra_contig_long drops to 0" in warning


def test_map_of_the_bins_from_a_1kb_map(metagenome, binning, tmp_path):
    """From a map of 1 kb bins: its bins, each contig in its final bin, the
    genomes in one block of bins each, and the contacts of the pairs by
    contig pair."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"], binsize=1000)
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    bins, source = cooler.Cooler(str(out)).bins()[:], cooler.Cooler(str(metacool)).bins()[:]
    columns = ["chrom", "start", "end"]
    assert sorted(map(tuple, bins[columns].astype(str).to_numpy())) == sorted(map(tuple, source[columns].astype(str).to_numpy()))
    assert assignment(out) == final_bins(metagenome)
    assert contig_contacts(out) == contig_contacts(metacool)
    pairs = metagenome.pairs
    assert contig_contacts(out) == contig_pairs(pairs.chr1, pairs.chr2).value_counts().to_dict()
    # Each genome is one block of bins.
    mc = metacooler.MetaCooler(str(out))
    for genome in mc.genomes.name:
        lo, hi = mc.genome_extent(genome)
        assert set(bins.chrom.iloc[lo:hi].astype(str)) == {name for name, found in assignment(out).items() if found == genome}
    # The genomes of the map of the pairs, but their number of bins.
    expected = genomes_without_qc(mtmc_map(binning, tmp_path))
    pd.testing.assert_frame_equal(genomes_without_qc(out).drop(columns="n_bins"), expected.drop(columns="n_bins"))
    assert (mc.genomes.n_bins > expected.n_bins).all()


def mtmc_map(binning, tmp_path):
    """The map that write_metacool writes from the pairs of the binning."""
    path = tmp_path / "pairs.metacool"
    if not path.exists():
        mtmc.write_metacool(
            str(binning["pairs"]), str(binning["contig_data"]), str(path), bin_summary=str(binning["bin_summary"])
        )
    return path


@pytest.mark.parametrize("binned", ["map", "assembly"])
def test_map_of_the_bins_from_a_cut_map(binned, metagenome, binning, tmp_path, logged):
    """From a map cut by metacooler extract (a .cool file), then made a map
    again by create cool: the contigs that it keeps, in their final bins.
    A binned contig that the map lacks (a binning made on more contigs than
    the map's) is left out, with a warning."""
    truth = metagenome.truth
    shortest = int(truth.Length.median())
    kept = set(truth.Name[truth.Length >= shortest])
    full = unbinned_map(tmp_path / "full.metacool", metagenome, binning["pairs"])
    cut = tmp_path / "cut.cool"
    metacooler.extract(str(full), filters=[f"contigs.length >= {shortest}"], output=str(cut))
    metacool = tmp_path / "cut.metacool"
    metacooler.from_cool(str(cut), str(metacool), contig_genome_map=None, contig_ordering="cool")
    contig_data = binning["contig_data"]
    final_bin = _bin_names(metagenome)
    if binned == "map":
        # As a run from the cut map bins them: the contigs it lacks are in no bin.
        final_bin = {name: (bin_name if name in kept else "ND") for name, bin_name in final_bin.items()}
        contig_data = write_contig_data(metagenome, tmp_path / "cut_contig_data.txt", final_bin)
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(contig_data), str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    assert assignment(out) == {name: found for name, found in final_bins(metagenome).items() if name in kept}
    assert contig_contacts(out) == contig_contacts(metacool)
    assert_contig_columns(out, contig_data)
    assert cooler.Cooler(str(out)).info["metadata"]["qc"]["source"] == "pixels"
    warnings = [message for level, message in logged() if level == "WARNING"]
    # The contig data lists contigs that the map lacks: annotate reads it,
    # ignoring their rows without a warning (unmatched="ignore"), and the
    # history of the map names the contig data.
    assert not [message for message in warnings if "Contigs table" in message]
    assert cooler.Cooler(str(out)).info["metadata"]["annotate"][-1]["contigs_table"] == str(contig_data)
    assert not [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]
    lacking = sorted(name for name, bin_name in final_bin.items() if bin_name != "ND" and name not in kept)
    if binned == "map":
        assert not lacking and not [message for message in warnings if "binned contig(s)" in message]
    else:
        assert lacking
        assert (
            f"{len(lacking)} binned contig(s) of {contig_data} are not in {metacool}, and not in the map of the bins either: "
        ) in "\n".join(warnings)


def test_map_of_the_bins_from_a_cut_map_other_warnings(metagenome, binning, tmp_path, logged):
    """Only the rows of the contig data that name no contig of the map are
    ignored without a warning: another warning of metacooler annotate about
    the table, a contig described by two rows here, names the contig data,
    as does the history of the map."""
    truth = metagenome.truth
    shortest = int(truth.Length.median())
    full = unbinned_map(tmp_path / "full.metacool", metagenome, binning["pairs"])
    cut = tmp_path / "cut.cool"
    metacooler.extract(str(full), filters=[f"contigs.length >= {shortest}"], output=str(cut))
    metacool = tmp_path / "cut.metacool"
    metacooler.from_cool(str(cut), str(metacool), contig_genome_map=None, contig_ordering="cool")
    table = pd.read_csv(binning["contig_data"], sep="\t", index_col=0)
    twice = truth.Name[truth.Length >= shortest].iloc[0]
    contig_data = tmp_path / "twice.txt"
    pd.concat([table, table[table.Name == twice]]).to_csv(contig_data, sep="\t")
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(contig_data), str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    warnings = [message for level, message in logged() if level == "WARNING"]
    assert [message for message in warnings if "Contigs table" in message] == [
        f"Contigs table ({contig_data}): 1 contig(s) are described by more than one row; the last row wins: {twice!r}"
    ]
    assert cooler.Cooler(str(out)).info["metadata"]["annotate"][-1]["contigs_table"] == str(contig_data)


def test_map_of_the_bins_from_an_mmetacool_entry(metagenome, binning, tmp_path):
    """From an entry of a .mmetacool collection: the map of the bins of the
    map it holds, a .metacool file of its own."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    uri = f"{collection}::entries/lib1"
    for source, out in [(str(metacool), tmp_path / "from_map.metacool"), (uri, tmp_path / "from_entry.metacool")]:
        mtmc.reassign_metacool(source, str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))
        metacooler.validate(str(out), full=True)
    assert_same_map(tmp_path / "from_entry.metacool", tmp_path / "from_map.metacool")
    metadata = cooler.Cooler(str(tmp_path / "from_entry.metacool")).info["metadata"]
    assert [entry["source"] for entry in metadata["reassign"]] == [uri]


def test_map_of_the_bins_from_a_map_without_qc(metagenome, binning, tmp_path):
    """From a map that has no QC block: the map of the bins has none either."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"], qc=False)
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))
    metacooler.validate(str(out), full=True)
    assert "qc" not in cooler.Cooler(str(out)).info["metadata"]
    assert not [column for column in metacooler.MetaCooler(str(out)).genomes.columns if column.startswith("qc_")]
    assert assignment(out) == final_bins(metagenome)


@pytest.mark.parametrize("source", ["unbinned", "binned"])
def test_map_of_no_bin_from_a_map(source, metagenome, binning, tmp_path, logged):
    """A run with no final bin: every contig of the map in __unbinned__, its
    only genome, with no column of a bin summary, even when the input map
    is that of an earlier binning, whose bin summary columns would be stale
    (and no warning of metacooler that __unbinned__ keeps them)."""
    if source == "unbinned":
        metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    else:
        metacool = mtmc_map(binning, tmp_path)
        assert "Completeness" in metacooler.MetaCooler(str(metacool)).genomes.columns
    contig_data = write_contig_data(metagenome, tmp_path / "no_bin.txt", dict.fromkeys(_bin_names(metagenome), "ND"))
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(contig_data), str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    mc = metacooler.MetaCooler(str(out))
    assert list(mc.genomes.name) == ["__unbinned__"]
    assert list(genomes_without_qc(out).columns) == ["name", "length", "n_contigs", "n_bins"]
    assert set(assignment(out).values()) == {"__unbinned__"} and len(mc.contigs) == len(metagenome.contigs)
    assert_contig_columns(out, contig_data)
    assert contig_contacts(out) == contig_contacts(metacool)
    # From a map of no bin either, no contig moves: metacooler says so at
    # INFO, as it is given a genomes table (the bin names, here none).
    assert [message for level, message in logged() if level == "WARNING"] == []


def test_map_of_the_bins_without_a_bin_summary(metagenome, binning, tmp_path, logged):
    """From the map of an earlier binning, another binning without a bin
    summary: no column of the earlier bin summary in the genomes table, and
    no warning of metacooler that a bin keeps them, which metaTOR's
    command, with no genomes table to give, could not act on. No file is
    left but the map."""
    earlier = mtmc_map(binning, tmp_path)
    final_bin = _bin_names(metagenome)
    merged = {final_bin[genome.contigs[0].name] for genome in metagenome.genomes[:2]}
    final_bin = {name: ("metator_00099_00000" if bin_name in merged else bin_name) for name, bin_name in final_bin.items()}
    contig_data = write_contig_data(metagenome, tmp_path / "new.txt", final_bin)
    out = tmp_path / "bins.metacool"
    before = set(tmp_path.iterdir())
    mtmc.reassign_metacool(str(earlier), str(contig_data), str(out))

    metacooler.validate(str(out), full=True)
    assert assignment(out) == {name: (bin_name if bin_name != "ND" else "__unbinned__") for name, bin_name in final_bin.items()}
    assert list(genomes_without_qc(out).columns) == ["name", "length", "n_contigs", "n_bins"]
    assert_contig_columns(out, contig_data)
    assert [message for level, message in logged() if level == "WARNING"] == []
    assert set(tmp_path.iterdir()) - before == {out}


def test_map_of_the_bins_replaces_the_columns_of_the_input_map(metagenome, binning, tmp_path):
    """From the map of an earlier binning: the bins, the genomes table and
    the contig columns of this run. A contig column that this run did not
    measure (Shotgun_coverage without a depth file) is dropped; a column
    that metaTOR does not write is kept."""
    earlier = mtmc_map(binning, tmp_path)
    taxa = tmp_path / "taxa.tsv"
    taxa.write_text("Name\ttaxon\n" + "".join(f"{name}\tt_{i % 3}\n" for i, name in enumerate(metagenome.truth.Name)))
    metacooler.annotate(str(earlier), contigs_table=str(taxa), contigs_table_format="Name", inplace=True)
    # Another binning: the first two genomes merged into one bin, with
    # another bin summary; another GC content, a restriction site count
    # and no shotgun coverage.
    final_bin = _bin_names(metagenome)
    merged = {final_bin[genome.contigs[0].name] for genome in metagenome.genomes[:2]}
    final_bin = {name: ("metator_00099_00000" if bin_name in merged else bin_name) for name, bin_name in final_bin.items()}
    table = pd.read_csv(write_contig_data(metagenome, tmp_path / "new.txt", final_bin), sep="\t", index_col=0)
    table["GC_content"] += 1.5
    table["Shotgun_coverage"] = "-"
    table["Restriction_site"] = range(len(table))
    contig_data = tmp_path / "contig_data_final.txt"
    table.to_csv(contig_data, sep="\t")
    summary = tmp_path / "new_summary.txt"
    names = sorted(set(final_bin.values()) - {"ND"})
    summary.write_text("\tContamination\n" + "".join(f"{name}\t0.{i}\n" for i, name in enumerate(names)))
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(earlier), str(contig_data), str(out), bin_summary=str(summary))

    metacooler.validate(str(out), full=True)
    assert assignment(out) == {name: (bin_name if bin_name != "ND" else "__unbinned__") for name, bin_name in final_bin.items()}
    assert list(genomes_without_qc(out).columns) == ["name", "length", "n_contigs", "n_bins", "Contamination"]
    # metaTOR's columns last, in their order, after the column it does not
    # write.
    contigs = metacooler.MetaCooler(str(out)).contigs
    assert without_qc(contigs.columns) == ["name", "length", "genome_id", "taxon", "GC_content", "Hit", "Restriction_site"]
    # Read back from the file, as metacooler does: the float parser of
    # pandas 2 can differ from the in-memory value in the last place.
    measured = pd.read_csv(contig_data, sep="\t", index_col="Name").loc[contigs.name]
    for column in ["GC_content", "Hit", "Restriction_site"]:
        assert contigs[column].tolist() == measured[column].tolist(), column
    assert contigs.set_index("name").taxon.to_dict() == {
        name: f"t_{i % 3}" for i, name in enumerate(metagenome.truth.Name)
    }  # fmt: skip


@pytest.mark.parametrize("earlier_bins", ["same", "other"])
def test_map_of_the_bins_from_the_map_of_an_earlier_run(earlier_bins, metagenome, binning, tmp_path):
    """From the map of an earlier run, which carries metaTOR's contig
    columns: the map that write_metacool writes from the pairs, contigs table
    included, metaTOR's columns last, as there. metacooler reassign stores
    the columns that it carries before the QC columns when it groups the
    contigs again (other bins): the contigs table is then that of the map
    of the pairs, column for column. When it relabels them in place (the
    same bins), it writes genome_id again, after the other columns, and the
    QC columns stay where they were: the contigs table is that of the map
    of the pairs but for the place of genome_id."""
    earlier = tmp_path / "earlier.metacool"
    if earlier_bins == "same":
        shutil.copyfile(mtmc_map(binning, tmp_path), earlier)
    else:
        # The first two genomes in one bin.
        final_bin = _bin_names(metagenome)
        merged = {final_bin[genome.contigs[0].name] for genome in metagenome.genomes[:2]}
        final_bin = {name: ("metator_00099_00000" if bin_name in merged else bin_name) for name, bin_name in final_bin.items()}
        contig_data = write_contig_data(metagenome, tmp_path / "earlier.txt", final_bin)
        mtmc.write_metacool(str(binning["pairs"]), str(contig_data), str(earlier))
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(earlier), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    expected = mtmc_map(binning, tmp_path)
    clr, ref = cooler.Cooler(str(out)), cooler.Cooler(str(expected))
    for table in ["bins", "pixels"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    chroms, ref_chroms = clr.chroms()[:], ref.chroms()[:]
    assert list(ref_chroms.columns) == ["name", "length", "genome_id", *QC_CONTIG_COLUMNS, *CONTIG_COLUMNS]
    if earlier_bins == "same":
        assert list(chroms.columns) == ["name", "length", *QC_CONTIG_COLUMNS, "genome_id", *CONTIG_COLUMNS]
        chroms = chroms[ref_chroms.columns]
    pd.testing.assert_frame_equal(chroms, ref_chroms)
    assert_contig_columns(out, binning["contig_data"])
    # The genomes, QC block and columns included, carried over to the bins.
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(out)).genomes, metacooler.MetaCooler(str(expected)).genomes)
    assert clr.info["metadata"]["qc"] == ref.info["metadata"]["qc"]


def test_map_of_the_bins_drops_the_columns_of_an_earlier_run(metagenome, binning, tmp_path):
    """From a map made with the contig data of an earlier run as its
    contigs table (metacooler create pairs --contigs-format metator, which
    stores all its columns): the bins and contig ids of that run, and its
    other columns, are not carried into the map of the bins; only this
    run's contig columns, and a column that metaTOR does not write."""
    table = pd.read_csv(binning["contig_data"], sep="\t", index_col=0)
    earlier = tmp_path / "earlier.metacool"
    metacooler.from_pairs(
        str(earlier),
        bins=pd.DataFrame({"chrom": table.Name, "start": 0, "end": table.Size}),
        pairs=str(binning["pairs"]),
        contig_genome_map=str(binning["contig_data"]),
        contig_table_format="metator",
        ps=False,
    )
    taxa = tmp_path / "taxa.tsv"
    taxa.write_text("Name\ttaxon\n" + "".join(f"{name}\tt_{i % 3}\n" for i, name in enumerate(table.Name)))
    metacooler.annotate(str(earlier), contigs_table=str(taxa), contigs_table_format="Name", inplace=True)
    stale = {"Unnamed: 0", "ID", "Size", "Restriction_site", "Overlapping_bin_ID", "Recursive_bin_contigs"}
    assert stale <= set(metacooler.MetaCooler(str(earlier)).contigs.columns)
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(earlier), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    assert assignment(out) == final_bins(metagenome)
    contigs = metacooler.MetaCooler(str(out)).contigs
    assert sorted(without_qc(contigs.columns)) == sorted(["name", "length", "genome_id", *CONTIG_COLUMNS, "taxon"])


@pytest.mark.parametrize("source", ["unbinned", "binned"])
def test_map_of_the_bins_from_a_contig_data_without_columns(source, metagenome, binning, tmp_path):
    """A contig data file with none of the contig columns (a hand-made one):
    no contig column in the map of the bins, and those that the input map
    carried (from an earlier run) are dropped."""
    if source == "unbinned":
        metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    else:
        metacool = mtmc_map(binning, tmp_path)
    contig_data = tmp_path / "no_columns.txt"
    pd.read_csv(binning["contig_data"], sep="\t", index_col=0)[["Name", "Size", "Final_bin"]].to_csv(contig_data, sep="\t")
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(contig_data), str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    assert_contig_columns(out, contig_data, [])
    assert assignment(out) == final_bins(metagenome)
    annotated = cooler.Cooler(str(out)).info["metadata"].get("annotate", [])
    if source == "unbinned":
        assert annotated == []
    else:
        assert annotated[-1]["dropped"] == [f"contigs.{column}" for column in CONTIG_COLUMNS]


def test_map_of_the_bins_from_a_map_of_other_lengths(metagenome, binning, tmp_path):
    """A contig data whose Size is not the length of a contig in the map (a
    binning of another version of the assembly): refused, with the contigs,
    before anything is written."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    table = pd.read_csv(binning["contig_data"], sep="\t", index_col=0)
    first, length = table.Name.iloc[0], int(table.Size.iloc[0])
    table.loc[table.index[0], "Size"] = length + 10
    contig_data = tmp_path / "other_lengths.txt"
    table.to_csv(contig_data, sep="\t")
    out = tmp_path / "bins.metacool"
    message = (
        f"{metacool} is not a contact map of the contigs of {contig_data}: 1 contig(s) have another length in it: "
        f"{first!r} ({length} bp in the map, {length + 10} bp in the contig data)."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        mtmc.reassign_metacool(str(metacool), str(contig_data), str(out), bin_summary=str(binning["bin_summary"]))
    assert sorted(path.name for path in tmp_path.iterdir() if "metacool" in path.name) == ["contigs.metacool"]


def test_map_of_the_bins_from_a_contig_data_of_other_contigs(metagenome, binning, tmp_path):
    """A contig data that names no contig of the map (the binning of
    another assembly): metacooler annotate refuses it, naming it, and no
    file is left."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    table = pd.read_csv(binning["contig_data"], sep="\t", index_col=0)
    table["Name"] = "other_" + table["Name"]
    contig_data = tmp_path / "other_contigs.txt"
    table.to_csv(contig_data, sep="\t")
    out = tmp_path / "out"
    message = f"The contigs table {str(contig_data)!r} matches no contig of "
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}"):
        mtmc.reassign_metacool(str(metacool), str(contig_data), str(out / "bins.metacool"))
    assert list(out.iterdir()) == []


@pytest.mark.parametrize("how", ["path", "root URI", "link to the outfile", "outfile a link"])
def test_map_of_the_bins_replaces_the_input_map(how, metagenome, binning, tmp_path, logged):
    """The outfile may be the input map itself (a rebinning into the same
    directory): the same path (also as the URI of its root, file::/), the
    input given as a symbolic link to the outfile, or the outfile a
    symbolic link to the input map. The map of the bins is written next to
    it, validated, and only then replaces it (the outfile, not the target of
    a link): no other file is left, and the log says it is written."""
    out = tmp_path / "out"
    out.mkdir()
    outfile = out / "metator.metacool"
    if how == "outfile a link":
        metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
        outfile.symlink_to(metacool)
    else:
        metacool = unbinned_map(outfile, metagenome, binning["pairs"])
        if how == "link to the outfile":
            metacool = tmp_path / "link.metacool"
            metacool.symlink_to(outfile)
    before = Path(mtmc.map_file(metacool)).read_bytes()
    if how == "root URI":
        metacool = f"{metacool}::/"
    assert mtmc.is_input_map(str(metacool), str(outfile))
    mtmc.reassign_metacool(
        str(metacool), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True
    )

    assert not outfile.is_symlink()
    metacooler.validate(str(outfile), full=True)
    assert_same_map(outfile, mtmc_map(binning, tmp_path))
    assert sorted(path.name for path in out.iterdir()) == ["metator.metacool"]
    assert ("INFO", f"{outfile} is written: the contact map of the bins.") in logged()
    if how == "outfile a link":
        # The input map itself is left as it was.
        assert metacool.read_bytes() == before
    else:
        assert Path(mtmc.map_file(metacool)).read_bytes() != before


def _fail_in_reassign(uri, table, output, **kwargs):
    """metacooler.reassign, failing as it does when it cannot write: told
    to remove its partial file (keep_partial=False), it leaves none."""
    assert kwargs["keep_partial"] is False
    raise OSError("No space left on device")


def _fail_in_validate(path, full=True):
    raise ValueError("pixels out of order")


def removed_partial(messages):
    """The file that metacooler was writing, from its INFO line: it removed
    its partial file (keep_partial=False)."""
    (staged,) = [
        re.match(r"Creation failed; the partial file (\S+)\.partial was removed and \1 was not written\.$", message)[1]
        for level, message in messages
        if level == "INFO" and message.startswith("Creation failed")
    ]
    return staged


@pytest.mark.parametrize("failure", ["reassign", "validate", "in metacooler"])
@pytest.mark.parametrize("how", ["path", "link to the outfile"])
def test_input_map_kept_when_the_map_of_the_bins_fails(how, failure, metagenome, binning, tmp_path, monkeypatch, logged):
    """When the map of the bins cannot be written or does not validate, the
    input map (the outfile) is left byte for byte as it was, and no
    temporary file stays next to it: metacooler removes the partial file it
    was writing, and says so at INFO, not that it is left for inspection."""
    out = tmp_path / "out"
    out.mkdir()
    outfile = unbinned_map(out / "metator.metacool", metagenome, binning["pairs"])
    metacool = outfile
    if how == "link to the outfile":
        metacool = tmp_path / "link.metacool"
        metacool.symlink_to(outfile)
    before = outfile.read_bytes()
    if failure == "reassign":
        monkeypatch.setattr(metacooler, "reassign", _fail_in_reassign)
        message = "No space left on device"
    elif failure == "validate":
        monkeypatch.setattr(metacooler, "validate", _fail_in_validate)
        message = "pixels out of order"
    else:
        # The validation of metacooler reassign, which writes the map to its
        # partial file first.
        monkeypatch.setattr("metacooler.io.reassign.validate", _fail_in_validate)
        message = "pixels out of order"
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$") as raised:
        mtmc.reassign_metacool(
            str(metacool), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True
        )
    assert str(raised.value.__cause__) == message
    assert outfile.read_bytes() == before
    assert sorted(path.name for path in out.iterdir()) == ["metator.metacool"]
    messages = logged()
    assert not [message for level, message in messages if level == "ERROR" and "for inspection" in message]
    if failure == "in metacooler":
        # The temporary file of the map, next to the outfile.
        staged = removed_partial(messages)
        assert Path(staged).parent == out and Path(staged).name.startswith(".metator.metacool.")


def _fail_unexpectedly(path, full=True):
    raise RuntimeError("unexpected failure")


def test_partial_file_removed_on_an_unexpected_error(metagenome, binning, tmp_path, monkeypatch, logged):
    """An error of another type than those metacooler is expected to raise,
    within metacooler reassign once it has started its partial file (in its
    validation, here), is raised as it is: metacooler removes its partial
    file, the temporary file of the map is removed, and the input map (the
    outfile) is left as it was, with no other file next to it."""
    out = tmp_path / "out"
    out.mkdir()
    outfile = unbinned_map(out / "metator.metacool", metagenome, binning["pairs"])
    before = outfile.read_bytes()
    monkeypatch.setattr("metacooler.io.reassign.validate", _fail_unexpectedly)
    with pytest.raises(RuntimeError, match=r"^unexpected failure$"):
        mtmc.reassign_metacool(
            str(outfile), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True
        )
    assert outfile.read_bytes() == before
    assert sorted(path.name for path in out.iterdir()) == ["metator.metacool"]
    messages = logged()
    staged = removed_partial(messages)
    assert Path(staged).parent == out and Path(staged).name.startswith(".metator.metacool.")
    assert [message for level, message in messages if level in ("WARNING", "ERROR")] == []


def test_metacool_failure_leaves_no_partial_file(metagenome, binning, tmp_path, monkeypatch, logged):
    """When metacooler cannot write the map of the pairs (its validation
    fails, here), write_metacool raises the reason, and leaves neither the
    map nor metacooler's partial file next to the outfile."""
    out = tmp_path / "out"
    out.mkdir()
    outfile = out / "metator.metacool"
    monkeypatch.setattr("metacooler.io.create.validate", _fail_in_validate)
    with pytest.raises(mtmc.MetacoolError, match=r"^pixels out of order$"):
        mtmc.write_metacool(
            str(binning["pairs"]), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"])
        )
    assert list(out.iterdir()) == []
    messages = logged()
    assert removed_partial(messages) == str(outfile)
    assert not [message for level, message in messages if level == "ERROR"]


def test_map_of_the_bins_existing_outfile(metagenome, binning, tmp_path):
    """An existing outfile, the input map or another file: refused without
    force, and left as it was."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    other = tmp_path / "other.metacool"
    other.write_text("old")
    for outfile in [metacool, other]:
        before = outfile.read_bytes()
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(f'{outfile} already exists. Use --force to overwrite it.')}$"):
            mtmc.reassign_metacool(str(metacool), str(binning["contig_data"]), str(outfile))
        assert outfile.read_bytes() == before
    mtmc.reassign_metacool(str(metacool), str(binning["contig_data"]), str(other), force=True)
    metacooler.validate(str(other), full=True)


@pytest.mark.parametrize(
    "error",
    [KeyError("genome_name"), OSError("No space left on device"), ValueError("bad")],
    ids=["KeyError", "OSError", "ValueError"],
)
def test_metacooler_errors_of_the_map_of_the_bins(error, metagenome, binning, tmp_path, monkeypatch):
    """What metacooler raises while writing the map of the bins is raised
    again as a MetacoolError, its message kept (not the repr of a
    KeyError), and leaves no file."""

    def reassign(*args, **kwargs):
        raise error

    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    monkeypatch.setattr(metacooler, "reassign", reassign)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(error.args[0])}$") as raised:
        mtmc.reassign_metacool(str(metacool), str(binning["contig_data"]), str(tmp_path / "bins.metacool"))
    assert raised.value.__cause__ is error
    assert sorted(path.name for path in tmp_path.iterdir() if "metacool" in path.name) == ["contigs.metacool"]


REASSIGN_ERRORS = {
    "missing map": "Input file(s) not found: missing.metacool.",
    "not a map": "Cannot read the contact map ",
    "missing entry": "Cannot read the contact map maps.mmetacool::entries/lib2: ",
    "collection as outfile": (
        "maps.mmetacool is the .mmetacool collection of the input map maps.mmetacool::entries/lib1: the map of the "
        "bins would replace the whole collection. Write it to another file."
    ),
    "contig data without Final_bin": "has no Final_bin column",
    "bin named NA": "Rename the 1 bin(s): 'NA'.",
    "Size not a whole number": "must give the length of each contig in bp, as a whole number: ",
}


@pytest.mark.parametrize("case", list(REASSIGN_ERRORS))
def test_map_of_the_bins_refused(case, metagenome, metagenome_network, binning, tmp_path, monkeypatch):
    """Inputs that cannot give the map of the bins: a MetacoolError, before
    anything is written."""
    monkeypatch.chdir(tmp_path)
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    metacooler.create_mmetacool(Path("maps.mmetacool"), {"lib1": metacool})
    uri, contig_data, outfile = str(metacool), binning["contig_data"], "bins.metacool"
    if case == "missing map":
        uri = "missing.metacool"
    elif case == "not a map":
        uri = str(binning["pairs"])
    elif case == "missing entry":
        uri = "maps.mmetacool::entries/lib2"
    elif case == "collection as outfile":
        uri, outfile = "maps.mmetacool::entries/lib1", "maps.mmetacool"
    elif case == "contig data without Final_bin":
        contig_data = metagenome_network["contig_data"]
    elif case == "bin named NA":
        contig_data = write_contig_data(
            metagenome, tmp_path / "na.txt", {**_bin_names(metagenome), metagenome.contigs[0].name: "NA"}
        )
    elif case == "Size not a whole number":
        name = metagenome.contigs[0].name
        contig_data.write_text(re.sub(rf"\t{name}\t\d+\t", f"\t{name}\t1x00\t", contig_data.read_text()))
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()}
    with pytest.raises(mtmc.MetacoolError, match=re.escape(REASSIGN_ERRORS[case])):
        mtmc.reassign_metacool(uri, str(contig_data), outfile, force=True)
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()} == before


# The map of the contigs, in no bin yet, that a pipeline run builds from its
# pairs files (write_contig_map), and the map of its bins written from it.


def assembly_contig_data(inputs, measured=True):
    """The contig data of the assembly of the mock community, as the
    pipeline makes it: with its depth file and an enzyme (the shotgun
    coverages and restriction sites measured), or neither."""
    depth, enzyme = (str(inputs["depth"]), "DpnII") if measured else (None, None)
    contig_data, _ = mtn.create_contig_data(str(inputs["assembly"]), 1, depth, enzyme)
    return contig_data


@pytest.mark.parametrize("measured", [True, False], ids=["depth and enzyme", "neither"])
def test_map_of_the_contigs(measured, metagenome, metagenome_inputs, tmp_path):
    """The map that a pipeline run builds from its pairs files: one bin per
    contig, every contig in __unbinned__, the longest first, every pair a
    contact, the QC block measured on the reads, and the contig columns
    measured on the assembly as the contig data files of the run hold them
    (Shotgun_coverage with a depth file, Restriction_site with an enzyme),
    but Hit, which the network counts. Only the map is written: no P(s)
    table next to it, no index next to the assembly, no table left in the
    temporary directory."""
    contig_data = assembly_contig_data(metagenome_inputs, measured)
    inputs = sorted(metagenome_inputs["assembly"].parent.iterdir())
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    out = tmp_path / "out" / "metator.metacool"
    mtmc.write_contig_map(
        str(metagenome_inputs["pairs"]),
        contig_data,
        str(out),
        assembly_name="synthetic",
        metadata={"command": "metator pipeline"},
        tmp_dir=str(tmp_dir),
    )

    metacooler.validate(str(out), full=True)
    mc, clr = metacooler.MetaCooler(str(out)), cooler.Cooler(str(out))
    assert list(mc.genomes.name) == ["__unbinned__"]
    contigs = mc.contigs
    truth = metagenome.truth
    assert sorted(contigs.name) == sorted(truth.Name)
    assert list(contigs.length) == sorted(contigs.length, reverse=True)
    bins = clr.bins()[:]
    assert list(bins.chrom.astype(str)) == list(contigs.name) and (bins.start == 0).all()
    assert (bins.end.to_numpy() == truth.Length[contigs.name].to_numpy()).all()
    assert map_pixels(out).to_dict() == expected_pixels(metagenome.pairs, list(contigs.name)).to_dict()
    metadata = clr.info["metadata"]
    assert metadata["qc"]["source"] == "pairs" and metadata["qc"]["n_pairs_read"] == len(metagenome.pairs)
    assert metadata["metator"] == {"version": __version__, "command": "metator pipeline"}
    assert clr.info["genome-assembly"] == "synthetic"
    assert mc.info["ordering"] == "genome-then-length"
    contig_data_file = tmp_path / "contig_data_network.txt"
    mtn.write_contig_data(contig_data, str(contig_data_file))
    columns = ["GC_content", "Shotgun_coverage", "Restriction_site"] if measured else ["GC_content"]
    assert_contig_columns(out, contig_data_file, columns)
    assert sorted(path.name for path in out.parent.iterdir()) == ["metator.metacool"]
    assert sorted(metagenome_inputs["assembly"].parent.iterdir()) == inputs
    # The table of the columns was written in tmp_dir, as the history of the
    # map records it, and removed.
    (annotated,) = metadata["annotate"]
    assert Path(annotated["contigs_table"]).parent.parent == tmp_dir
    assert list(tmp_dir.iterdir()) == []


def final_contig_data(contig_data, pairs, final_bin, path):
    """The contig_data_final.txt file of a run with this contig data, laid
    out as validation writes it (an unnamed index first, Final_bin last):
    its columns as write_contig_data writes them, the hits of each contig
    counted on the pairs, as the network step counts them."""
    hits = pd.concat([pairs.chr1, pairs.chr2]).value_counts()
    for name, data in contig_data.items():
        data["hit"] = int(hits.get(name, 0))
    mtn.write_contig_data(contig_data, str(path))
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    table["Final_bin"] = table.Name.map(final_bin)
    table.to_csv(path, sep="\t", index=True)
    return path


@pytest.mark.parametrize("final", ["bins", "no bin", "one bin", "longest first"])
def test_map_of_the_bins_from_the_map_of_the_contigs(final, metagenome, metagenome_inputs, tmp_path):
    """The map of the bins that reassign_metacool writes over the map of the
    contigs, as a pipeline run does at its end, is the map that
    write_metacool writes from the same pairs files and final tables: the
    same bins, contacts, contigs table (metaTOR's columns last, as there),
    genomes table (QC columns included) and QC block, but for its history
    (the annotate and reassign keys of its metadata) and its date. When the
    final bins leave the contigs in the order of the map of the contigs (for
    example, no final bin, one bin holding every contig, or the longest
    contigs in the first bin and the others in the second), metacooler
    relabels them where they are, and writes genome_id again, after the QC
    columns of the contigs."""
    pairs = str(metagenome_inputs["pairs"])
    contig_data = assembly_contig_data(metagenome_inputs)
    out = tmp_path / "out"
    outfile = out / "metator.metacool"
    run = {"metadata": {"command": "metator pipeline"}}
    mtmc.write_contig_map(pairs, contig_data, str(outfile), assembly_name="synthetic", **run)
    final_bin = _bin_names(metagenome)
    if final in ["no bin", "one bin"]:
        final_bin = dict.fromkeys(final_bin, "ND" if final == "no bin" else "metator_00001_00000")
    elif final == "longest first":
        # The 20 first contigs of the map of the contigs, the longest, in
        # the first bin, and the others in the second.
        order = metacooler.MetaCooler(str(outfile)).contigs.name
        final_bin = {name: f"metator_0000{1 + (i >= 20)}_00000" for i, name in enumerate(order)}
    contig_data_file = final_contig_data(contig_data, metagenome.pairs, final_bin, out / "contig_data_final.txt")
    tables = {"bin_summary": str(write_bin_summary(metagenome, out / "bin_summary.txt")), **run}
    mtmc.reassign_metacool(str(outfile), str(contig_data_file), str(outfile), force=True, **tables)
    expected = tmp_path / "pairs.metacool"
    mtmc.write_metacool(pairs, str(contig_data_file), str(expected), assembly_name="synthetic", **tables)

    metacooler.validate(str(outfile), full=True)
    clr, ref = cooler.Cooler(str(outfile)), cooler.Cooler(str(expected))
    for table in ["bins", "pixels"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    chroms, ref_chroms = clr.chroms()[:], ref.chroms()[:]
    columns = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    assert list(ref_chroms.columns) == ["name", "length", "genome_id", *QC_CONTIG_COLUMNS, *columns]
    if final != "bins":
        assert list(chroms.columns) == ["name", "length", *QC_CONTIG_COLUMNS, "genome_id", *columns]
        chroms = chroms[ref_chroms.columns]
    pd.testing.assert_frame_equal(chroms, ref_chroms)
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(outfile)).genomes, metacooler.MetaCooler(str(expected)).genomes)
    info, ref_info = dict(clr.info), dict(ref.info)
    metadata, ref_metadata = info.pop("metadata"), ref_info.pop("metadata")
    assert metadata["qc"] == ref_metadata["qc"] and metadata["qc"]["source"] == "pairs"
    assert {key for key in info.keys() | ref_info.keys() if info.get(key) != ref_info.get(key)} <= {"creation-date"}
    assert {key for key in metadata.keys() | ref_metadata.keys() if metadata.get(key) != ref_metadata.get(key)} == {
        "annotate",
        "reassign",
    }
    # The history: the columns measured on the assembly, added to the map of
    # the contigs, dropped from the map of the bins, then those of the
    # contig data added again, last.
    assert [entry["source"] for entry in metadata["reassign"]] == [str(outfile)]
    assert [(entry["added"], entry["dropped"]) for entry in metadata["annotate"]] == [
        (["contigs.GC_content", "contigs.Shotgun_coverage", "contigs.Restriction_site"], []),
        ([], ["contigs.GC_content", "contigs.Shotgun_coverage", "contigs.Restriction_site"]),
        ([f"contigs.{column}" for column in columns], []),
    ]
    assert sorted(path.name for path in out.iterdir()) == ["bin_summary.txt", "contig_data_final.txt", "metator.metacool"]


def test_map_of_the_contigs_of_an_assembly_with_an_empty_contig(metagenome, metagenome_inputs, tmp_path, logged):
    """An assembly that holds a contig of no length (an empty fasta record):
    metacooler leaves it out of the map, and says so; the columns of the
    other contigs are stored, without a warning about the row of the table
    that names it."""
    contig_data = assembly_contig_data(metagenome_inputs)
    contig_data["empty"] = {"id": len(contig_data) + 1, "length": 0, "GC": 0.0, "hit": 0, "coverage": 0.0, "RS": 1}
    out = tmp_path / "map.metacool"
    mtmc.write_contig_map(str(metagenome_inputs["pairs"]), contig_data, str(out))

    metacooler.validate(str(out), full=True)
    contigs = metacooler.MetaCooler(str(out)).contigs
    assert sorted(contigs.name) == sorted(metagenome.truth.Name)
    assert contigs.GC_content.notna().all()
    warnings = [message for level, message in logged() if level == "WARNING"]
    assert warnings == ["Dropping 1 contig(s) with non-positive length: ['empty']"]


CONTIG_MAP_ERRORS = {
    "threads below 1": "The number of threads must be at least 1, got 0.",
    "existing outfile": "map.metacool already exists. Use --force to overwrite it.",
    "missing pairs file": "Input file(s) not found: missing.pairs.",
    "contig name not ASCII": "A .metacool file needs ASCII contig names (cooler reads them as ASCII text), and 1 "
    "contig name(s) are not: 'contig_é'.",
}


@pytest.mark.parametrize("case", list(CONTIG_MAP_ERRORS))
def test_map_of_the_contigs_refused(case, metagenome_files, tmp_path, monkeypatch):
    """Inputs that cannot give the map of the contigs: a MetacoolError,
    before the pairs are read, and nothing written."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(metacooler, "from_pairs", lambda *args, **kwargs: pytest.fail("pairs read"))
    contig_data = {"c1": {"length": 100, "GC": 40.5, "coverage": "-", "RS": "-"}}
    pairs, kwargs = [str(metagenome_files["pairs"])], {}
    if case == "threads below 1":
        kwargs["threads"] = 0
    elif case == "existing outfile":
        Path("map.metacool").write_text("earlier")
    elif case == "missing pairs file":
        pairs.append("missing.pairs")
    else:
        contig_data["contig_é"] = contig_data["c1"]
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()}
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(CONTIG_MAP_ERRORS[case])}$"):
        mtmc.write_contig_map(pairs, contig_data, "map.metacool", **kwargs)
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()} == before


@pytest.mark.parametrize(
    "failure, error",
    [
        ("map", ValueError("pixels out of order")),
        ("columns", ValueError("bad table")),
        ("columns", RuntimeError("unexpected")),
        ("columns", KeyboardInterrupt()),
    ],
    ids=["map", "columns", "columns, unexpected", "columns, interrupted"],
)
def test_map_of_the_contigs_failure_leaves_nothing(failure, error, metagenome_inputs, tmp_path, monkeypatch, logged):
    """When metacooler cannot write the map (its validation fails, here) or
    store its contig columns, write_contig_map raises the reason (a
    MetacoolError, or an error of another type as it is, an interruption
    (Ctrl-C) included), and leaves neither the map nor metacooler's partial
    file next to the outfile."""

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr("metacooler.io.create.validate" if failure == "map" else "metacooler.annotate", fail)
    out = tmp_path / "out"
    out.mkdir()
    outfile = out / "metator.metacool"
    raised = mtmc.MetacoolError if isinstance(error, ValueError) else type(error)
    with pytest.raises(raised, match=f"^{re.escape(str(error))}$"):
        mtmc.write_contig_map(str(metagenome_inputs["pairs"]), assembly_contig_data(metagenome_inputs), str(outfile))
    assert list(out.iterdir()) == []
    if failure == "map":
        assert removed_partial(logged()) == str(outfile)
    assert not [message for level, message in logged() if level == "ERROR"]


def _split_pairs(pairs, directory, n):
    """The pairs of a file dealt into n libraries, each with its header."""
    lines = (
        Path(pairs).read_text().splitlines(keepends=True)
        if not str(pairs).endswith(".gz")
        else gzip.decompress(Path(pairs).read_bytes()).decode().splitlines(keepends=True)
    )
    header = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    paths = []
    for i in range(n):
        path = Path(directory) / f"lib{i}.pairs"
        path.write_text("".join(header + body[i::n]))
        paths.append(str(path))
    return paths


@pytest.mark.parametrize("libraries", [2, 3])
def test_the_maps_of_the_libraries_and_their_sum(metagenome_inputs, tmp_path, libraries, logged):
    """The map of each library is the map of its pairs file alone, but for
    the columns measured on the assembly, which the map of their sum holds;
    and the map of their sum, merged from them, is the map of all the pairs
    files at once, its contact statistics included, but for its history
    and the type of its counts: int64, which metacooler merge writes so
    that a sum cannot overflow, where the map of all the pairs files at
    once stores int32. The run's metadata, which merge does not carry, is
    stored in the sum; and the log names the pairs file of each library.
    An earlier run's map at the outfile is replaced (force)."""
    contig_data = assembly_contig_data(metagenome_inputs)
    pairs = _split_pairs(metagenome_inputs["pairs"], tmp_path, libraries)
    out, collection = tmp_path / "out" / "run.metacool", tmp_path / "out" / "run.mmetacool"
    out.parent.mkdir()
    out.write_text("an earlier run's")
    mtmc.write_library_maps(pairs, contig_data, str(out), str(collection), force=True, metadata={"command": "metator x"})

    expected = tmp_path / "all.metacool"
    mtmc.write_contig_map(pairs, contig_data, str(expected), metadata={"command": "metator x"})
    metacooler.validate(str(out), full=True)
    clr, ref = cooler.Cooler(str(out)), cooler.Cooler(str(expected))
    summed, at_once = clr.pixels()[:], ref.pixels()[:]
    assert (summed["count"].dtype, at_once["count"].dtype) == (np.int64, np.int32)
    pd.testing.assert_frame_equal(summed, at_once.astype({"count": np.int64}))
    pd.testing.assert_frame_equal(clr.bins()[:], ref.bins()[:])
    pd.testing.assert_frame_equal(clr.chroms()[:], ref.chroms()[:])
    pd.testing.assert_frame_equal(
        metacooler.MetaCooler(str(out)).genomes, metacooler.MetaCooler(str(expected)).genomes, check_dtype=False
    )
    assert clr.info["metadata"]["qc"] == ref.info["metadata"]["qc"]
    assert clr.info["metadata"]["metator"] == ref.info["metadata"]["metator"] == {"version": __version__, "command": "metator x"}
    for i, library in enumerate(pairs):
        assert ("INFO", f"Library {i} ({mtmc.library_entry(str(collection), i)}): {library}.") in logged()
        entry = mtmc.library_entry(str(collection), i)
        metacooler.validate(entry, full=True)
        alone = tmp_path / f"alone{i}.metacool"
        mtmc.write_contig_map(library, contig_data, str(alone))
        own, its = cooler.Cooler(entry), cooler.Cooler(str(alone))
        pd.testing.assert_frame_equal(own.pixels()[:], its.pixels()[:])
        assert own.info["metadata"]["qc"] == its.info["metadata"]["qc"]
        assert "GC_content" not in metacooler.MetaCooler(entry).contigs.columns
    assert sorted(path.name for path in out.parent.iterdir()) == ["run.metacool", "run.mmetacool"]


@pytest.mark.parametrize("error", [KeyboardInterrupt, RuntimeError])
@pytest.mark.parametrize("where", ["the second library", "the collection", "the merge", "the columns of the sum"])
def test_the_maps_of_the_libraries_are_both_written_or_neither(metagenome_inputs, tmp_path, monkeypatch, where, error):
    """A guard: a failure while the maps are written (the map of a library,
    the collection, the merge, the columns of the sum) leaves neither file,
    no file of metacooler's, and none in the run's temporary directory."""
    contig_data = assembly_contig_data(metagenome_inputs)
    pairs = _split_pairs(metagenome_inputs["pairs"], tmp_path, 2)
    out, collection = tmp_path / "out" / "run.metacool", tmp_path / "out" / "run.mmetacool"
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    target, name = {
        "the second library": (metacooler, "from_pairs"),
        "the collection": (metacooler, "create_mmetacool"),
        "the merge": (metacooler, "merge"),
        "the columns of the sum": (mtmc, "_store_assembly_columns"),
    }[where]
    original, calls = getattr(target, name), []

    def failing(*args, **kwargs):
        calls.append(args)
        if where != "the second library" or len(calls) == 2:
            raise error("stopped")
        return original(*args, **kwargs)

    monkeypatch.setattr(target, name, failing)
    with pytest.raises(error):
        mtmc.write_library_maps(pairs, contig_data, str(out), str(collection), tmp_dir=str(tmp_dir))
    assert list(out.parent.iterdir()) == [] and list(tmp_dir.iterdir()) == []


def test_the_maps_of_the_libraries_are_written_in_the_temporary_directory(metagenome_inputs, tmp_path, monkeypatch):
    """The map of each library is written in the run's temporary directory,
    then copied into the collection: not next to the collection, where a
    run killed meanwhile left them (in a directory of no run's)."""
    contig_data = assembly_contig_data(metagenome_inputs)
    pairs = _split_pairs(metagenome_inputs["pairs"], tmp_path, 2)
    out, collection = tmp_path / "out" / "run.metacool", tmp_path / "out" / "run.mmetacool"
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    from_pairs, written = metacooler.from_pairs, []

    def recorded(path, *args, **kwargs):
        written.append(Path(path))
        return from_pairs(path, *args, **kwargs)

    monkeypatch.setattr(metacooler, "from_pairs", recorded)
    mtmc.write_library_maps(pairs, contig_data, str(out), str(collection), tmp_dir=str(tmp_dir))
    assert [path.parent.parent for path in written] == [tmp_dir, tmp_dir]
    assert list(tmp_dir.iterdir()) == [] and sorted(path.name for path in out.parent.iterdir()) == [
        "run.metacool",
        "run.mmetacool",
    ]


def _empty_pairs(path):
    """A pairs file of no pair: its header only."""
    Path(path).write_text("## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
    return str(path)


def test_a_library_without_contacts(metagenome_inputs, tmp_path):
    """A guard: a library without any contact is allowed. Its map has none,
    the sum is the map of the other library, its network is empty and its
    hits are all 0 (the pairs path divided by zero there)."""
    contig_data = assembly_contig_data(metagenome_inputs)
    (pairs,) = _split_pairs(metagenome_inputs["pairs"], tmp_path, 1)
    libraries = [pairs, _empty_pairs(tmp_path / "empty.pairs")]
    out, collection = tmp_path / "run.metacool", tmp_path / "run.mmetacool"
    mtmc.write_library_maps(libraries, contig_data, str(out), str(collection))
    entries = [mtmc.library_entry(str(collection), i) for i in range(2)]
    assert [metacooler.MetaCooler(entry).info["nnz"] for entry in entries][1] == 0
    alone = tmp_path / "alone.metacool"
    mtmc.write_contig_map(pairs, contig_data, str(alone))
    pd.testing.assert_frame_equal(
        cooler.Cooler(str(out)).pixels()[:], cooler.Cooler(str(alone)).pixels()[:].astype({"count": np.int64})
    )
    net = tmp_path / "net"
    net.mkdir()
    mtn.metacool_to_contacts(str(out), contig_data, str(net), "network.txt", "contig_data_network.txt", str(net), 1,
                             "empirical_hit", libraries=entries)  # fmt: skip
    assert (net / "network_1.txt").read_text() == "" and (net / "network_0.txt").read_text() != ""
    hits = pd.read_csv(net / "hit_data_alignment.txt", sep="\t", header=None)
    assert (hits.iloc[:, -1] == 0).all() and (hits.iloc[:, -2] > 0).any()


def test_libraries_without_any_contact(metagenome_inputs, tmp_path):
    """No library holds a contact: cooler merges no map without one, so the
    map of the sum is the map of all the pairs files at once, with no
    contact either, whose network the run refuses (it stopped inside the
    merge, on "No objects to concatenate")."""
    contig_data = assembly_contig_data(metagenome_inputs)
    libraries = [_empty_pairs(tmp_path / f"empty{i}.pairs") for i in range(2)]
    out, collection = tmp_path / "run.metacool", tmp_path / "run.mmetacool"
    mtmc.write_library_maps(libraries, contig_data, str(out), str(collection))
    metacooler.validate(str(out), full=True)
    assert metacooler.MetaCooler(str(out)).info["nnz"] == 0
    assert cooler.Cooler(str(out)).info["metadata"]["metator"] == {"version": __version__}
    with pytest.raises(mtmc.MetacoolError, match="holds no contact: there is no network to build"):
        mtn.metacool_to_contacts(str(out), contig_data, str(tmp_path), "network.txt", "contig_data_network.txt",
                                 str(tmp_path), 1, "None")  # fmt: skip


def test_a_library_that_cannot_be_read_is_named(metagenome_inputs, tmp_path):
    """The map of each library is written from a temporary file: an error
    names the library and its pairs file, not the temporary file."""
    contig_data = assembly_contig_data(metagenome_inputs)
    pairs = _split_pairs(metagenome_inputs["pairs"], tmp_path, 2)
    lines = Path(pairs[1]).read_text().splitlines(keepends=True)
    # A row with an extra field, which metacooler refuses.
    lines[-1] = lines[-1].rstrip("\n") + "\textra\n"
    Path(pairs[1]).write_text("".join(lines))
    out, collection = tmp_path / "run.metacool", tmp_path / "run.mmetacool"
    with pytest.raises(mtmc.MetacoolError, match=f"^Library 1, {re.escape(pairs[1])}: "):
        mtmc.write_library_maps(pairs, contig_data, str(out), str(collection))
    assert not out.exists() and not collection.exists()


@pytest.mark.parametrize(
    "case, message",
    [
        ("threads below 1", "The number of threads must be at least 1, got 0."),
        ("existing collection", "run.mmetacool already exists. Use --force to overwrite it."),
    ],
)
def test_the_maps_of_the_libraries_refuse(metagenome_inputs, tmp_path, case, message):
    contig_data = assembly_contig_data(metagenome_inputs)
    pairs = _split_pairs(metagenome_inputs["pairs"], tmp_path, 2)
    out, collection = tmp_path / "run.metacool", tmp_path / "run.mmetacool"
    kwargs = {"threads": 0} if case == "threads below 1" else {}
    if case == "existing collection":
        collection.write_text("an earlier run's")
    with pytest.raises(mtmc.MetacoolError, match=re.escape(message)):
        mtmc.write_library_maps(pairs, contig_data, str(out), str(collection), **kwargs)


# The normalization of the network, which the map of the bins stores, with
# metacooler.store_weights: from metacooler's next release (js2264/metacooler
# pull request #19).
needs_store_weights = pytest.mark.skipif(
    not mtmc.stores_weights(), reason="metacooler.store_weights comes with metacooler's next release"
)


def _contig_map_and_final_tables(metagenome, inputs, tmp_path, measured=True):
    """The map of the contigs of a pipeline run, built from its pairs, and
    its final contig data (the hits counted on the pairs, as the network
    step counts them). Returns the map, the contig data file and the
    contig data, whose hits are those of the file."""
    contig_data = assembly_contig_data(inputs, measured)
    outfile = tmp_path / "metator.metacool"
    mtmc.write_contig_map(str(inputs["pairs"]), contig_data, str(outfile), assembly_name="synthetic")
    contig_data_file = final_contig_data(
        contig_data, metagenome.pairs, _bin_names(metagenome), tmp_path / "contig_data_final.txt"
    )
    return outfile, contig_data_file, contig_data


def _weights(path):
    """The weight columns of a map, by name: their method and attributes."""
    from metacooler.io.weights import describe_weights

    return {row["name"]: row for row in describe_weights(str(path))}


def _normalized_contacts(path, name):
    """The contacts between contigs of a map, normalized by its weights
    ``name`` and summed by pair of contigs, as metacooler network --balance
    writes them: the contacts of a contig with no weight are left out."""
    edges = Path(path).parent / f"{name}.tsv"
    metacooler.network(str(path), level="contig", balance=name, out=str(edges), header=True)
    table = pd.read_csv(edges, sep="\t")
    pairs = contig_pairs(table.node1, table.node2)
    pairs["value"] = table.value.to_numpy()
    return {(a, b): value for a, b, value in pairs.itertuples(index=False) if a != b}


@needs_store_weights
@pytest.mark.parametrize("normalization", ["empirical_hit", "length", "RS", "abundance", "theoritical_hit"])
def test_map_of_the_bins_stores_the_normalization_of_the_network(normalization, metagenome, metagenome_inputs, tmp_path, logged):
    """Written over the map of the contigs, as a pipeline run does at its
    end, the map of the bins stores the normalization of the network, as
    the weights metator_<normalization> of its bins, from the columns of its
    contigs table that the normalization reads, which it records. The
    contacts of the map normalized by them and summed by pair of contigs
    are metaTOR's network: normalize_pair of the contact count of each pair
    (none where it gives 0)."""
    outfile, contig_data_file, contig_data = _contig_map_and_final_tables(metagenome, metagenome_inputs, tmp_path)
    mtmc.reassign_metacool(str(outfile), str(contig_data_file), str(outfile), force=True, normalization=normalization)

    metacooler.validate(str(outfile), full=True)
    name = f"metator_{normalization}"
    (row,) = _weights(outfile).values()
    attrs = row["attrs"]
    assert (row["name"], row["method"], attrs["mode"], attrs["metator_version"]) == (name, "metator", normalization, __version__)
    columns = [f"contigs/{column}" for column in mtmc.NORMALIZATION_COLUMNS[normalization]]
    assert attrs["depends_on"] == (["pixels"] if "contigs/Hit" in columns else []) + columns
    expected = {
        pair: mtn.normalize_pair(contig_data, pair, count, normalization)
        for pair, count in contig_contacts(outfile).items()
        if pair[0] != pair[1]
    }
    normalized = _normalized_contacts(outfile, name)
    assert normalized.keys() == {pair for pair, value in expected.items() if value}
    assert all(normalized[pair] == pytest.approx(expected[pair], rel=1e-9) for pair in normalized)
    assert not [message for level, message in logged() if level == "WARNING"]


@needs_store_weights
def test_map_of_the_bins_at_1kb_stores_the_weight_of_each_contig_on_its_bins(metagenome, metagenome_inputs, tmp_path):
    """From a map of the contigs at a fixed bin size, made before any
    binning: every bin of a contig has the weight of its contig."""
    contig_data = assembly_contig_data(metagenome_inputs)
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, metagenome_inputs["pairs"], binsize=1000)
    contig_data_file = final_contig_data(
        contig_data, metagenome.pairs, _bin_names(metagenome), tmp_path / "contig_data_final.txt"
    )
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(metacool), str(contig_data_file), str(out), normalization="empirical_hit")

    weights = metacooler.read_weights(str(out), "metator_empirical_hit")
    bins = cooler.Cooler(str(out)).bins()[:]
    assert len(bins) > len(contig_data)
    hits = bins.chrom.astype(str).map({name: data["hit"] for name, data in contig_data.items()}).to_numpy(dtype=float)
    with np.errstate(divide="ignore"):
        expected = 1 / np.sqrt(hits)
    expected[~np.isfinite(expected)] = np.nan
    np.testing.assert_allclose(weights, expected)


@needs_store_weights
def test_map_of_the_bins_without_a_column_that_the_normalization_reads(metagenome, metagenome_inputs, tmp_path, logged):
    """A run without an enzyme has no restriction sites: the map of the bins
    is written, without the RS normalization, and a warning says why."""
    outfile, contig_data_file, _ = _contig_map_and_final_tables(metagenome, metagenome_inputs, tmp_path, measured=False)
    mtmc.reassign_metacool(str(outfile), str(contig_data_file), str(outfile), force=True, normalization="RS")

    metacooler.validate(str(outfile), full=True)
    assert _weights(outfile) == {}
    (warning,) = [message for level, message in logged() if level == "WARNING"]
    assert warning == (
        f"{outfile} stores no RS normalization of the network. The RS normalization of the network reads the "
        "Restriction_site column(s) of the contigs table, which the map lacks."
    )


@needs_store_weights
@pytest.mark.parametrize("normalization", ["abundance", "None"])
def test_map_of_the_bins_replaces_the_normalization_of_an_earlier_run(
    normalization, metagenome, metagenome_inputs, tmp_path, logged
):
    """From the map of an earlier run, which stores the normalization of its
    network: the map of the bins stores this run's alone, or none with
    "None", without a warning. A weight that another tool stored is kept."""
    earlier, contig_data_file, _ = _contig_map_and_final_tables(metagenome, metagenome_inputs, tmp_path)
    mtmc.reassign_metacool(str(earlier), str(contig_data_file), str(earlier), force=True, normalization="empirical_hit")
    n_bins = len(cooler.Cooler(str(earlier)).bins()[:])
    metacooler.store_weights(str(earlier), "other", np.ones(n_bins), method="other", depends_on=["pixels"])
    out = tmp_path / "bins.metacool"
    mtmc.reassign_metacool(str(earlier), str(contig_data_file), str(out), normalization=normalization)

    weights = _weights(out)
    assert sorted(weights) == sorted(["other", *([] if normalization == "None" else [f"metator_{normalization}"])])
    assert weights["other"]["method"] == "other"
    assert not [message for level, message in logged() if level == "WARNING"]
    assert (
        "INFO",
        "Removed the normalization(s) of the network of an earlier run from the map: 'metator_empirical_hit'.",
    ) in logged()


def test_store_normalization_hands_its_weights_to_metacooler(metagenome, metagenome_inputs, tmp_path, monkeypatch):
    """What store_normalization gives metacooler.store_weights, whatever the
    version of metacooler: the weights of normalization_weights from the
    columns of the map, on each bin of the contig, and what they depend on.
    "None" and an unknown normalization have none."""
    outfile, _, contig_data = _contig_map_and_final_tables(metagenome, metagenome_inputs, tmp_path)
    # The map of the contigs has no Hit column: add the hits of the run.
    hits = pd.DataFrame({"Name": list(contig_data), "Hit": [data["hit"] for data in contig_data.values()]})
    hits.to_csv(tmp_path / "hits.tsv", sep="\t", index=False)
    metacooler.annotate(str(outfile), contigs_table=str(tmp_path / "hits.tsv"), contigs_table_format="Name", inplace=True)
    stored = {}

    def store_weights(uri, name, weights, **kwargs):
        stored.update(uri=uri, name=name, weights=weights, **kwargs)

    monkeypatch.setattr(metacooler, "store_weights", store_weights, raising=False)
    mtmc.store_normalization(str(outfile), "length")

    names = cooler.Cooler(str(outfile)).bins()[:].chrom.astype(str)
    with np.errstate(divide="ignore"):
        expected = np.sqrt(names.map({name: data["length"] / data["hit"] for name, data in contig_data.items()}).to_numpy(float))
    expected[~np.isfinite(expected)] = np.nan
    np.testing.assert_allclose(stored.pop("weights"), expected)
    assert stored == {
        "uri": str(outfile),
        "name": "metator_length",
        "method": "metator",
        "depends_on": ["pixels", "contigs/Hit"],
        "parameters": {"mode": "length", "metator_version": __version__},
        "force": True,
    }
    for normalization in ["None", "bogus"]:
        with pytest.raises(mtmc.MetacoolError, match=f"No weights for the normalization '{normalization}'"):
            mtmc.store_normalization(str(outfile), normalization)
    metacooler.annotate(str(outfile), drop=["contigs.Shotgun_coverage", "contigs.Restriction_site"], inplace=True)
    with pytest.raises(mtmc.MetacoolError, match="reads the Shotgun_coverage and Restriction_site column"):
        mtmc.store_normalization(str(outfile), "theoritical_hit")


def test_map_of_the_bins_with_a_metacooler_that_stores_no_weights(metagenome, metagenome_inputs, tmp_path, monkeypatch, logged):
    """With metacooler 0.3.0, which has no store_weights: the map of the
    bins is written as before, without the normalization, which the log
    says, without a warning; store_normalization refuses."""
    # metacooler exports store_weights lazily, through its __getattr__.
    monkeypatch.setattr(mtmc, "stores_weights", lambda: False)
    outfile, contig_data_file, _ = _contig_map_and_final_tables(metagenome, metagenome_inputs, tmp_path)
    mtmc.reassign_metacool(str(outfile), str(contig_data_file), str(outfile), force=True, normalization="empirical_hit")

    metacooler.validate(str(outfile), full=True)
    assert list(cooler.Cooler(str(outfile)).bins()[:].columns) == ["chrom", "start", "end"]
    assert not [message for level, message in logged() if level == "WARNING"]
    assert any(
        level == "INFO" and message.startswith(f"{outfile} stores no empirical_hit normalization of the network: metacooler")
        for level, message in logged()
    )
    with pytest.raises(mtmc.MetacoolError, match="cannot store weights"):
        mtmc.store_normalization(str(outfile), "empirical_hit")
