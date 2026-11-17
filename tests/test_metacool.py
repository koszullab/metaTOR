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


def assert_contig_columns(path, contig_data, columns=CONTIG_COLUMNS):
    """The contigs table of the map (cooler's chroms) holds these columns of
    the contig data, with their values and types, and no other: no bin, no
    pandas index."""
    contigs = metacooler.MetaCooler(str(path)).contigs
    assert list(contigs.columns) == ["name", "length", "genome_id", *columns]
    chroms = cooler.Cooler(str(path)).chroms()[:]
    assert list(chroms.columns) == ["name", "length", "genome_id", *columns]
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
    writes from the pairs, but for its QC, recomputed from the pixels, and
    its provenance. The input map is only read, and nothing but the map is
    left next to it."""
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
    # The QC of the pixels; the assembly and the history of the input map.
    assert metadata["qc"]["source"] == "pixels" and metadata["qc"]["n_contacts"] == len(metagenome.pairs)
    assert [entry["source"] for entry in metadata["reassign"]] == [str(metacool)]
    assert clr.info["genome-assembly"] == cooler.Cooler(str(metacool)).info["genome-assembly"]
    assert sorted(path.name for path in out.parent.iterdir()) == ["bins.metacool"]
    # Read and written as the map of the pairs (not the private temporary file).
    assert out.stat().st_mode == expected.stat().st_mode


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
    # its warning about their rows is dropped, and the history of the map
    # names the contig data.
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
    """Only the warning of metacooler annotate about the rows of the contig
    data that name no contig of the map is dropped: another warning about
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


def test_rows_of_other_contigs_ignored(tmp_path, logged):
    """The filter drops the warning of metacooler annotate about the rows of
    the contig data that name no contig of the map, and that one only:
    neither the same warning about another table, nor another message about
    the contig data, nor this message at another level. It is removed on
    exit."""
    contig_data = str(tmp_path / "contig_data_final.txt")
    other = str(tmp_path / "other.txt")
    unmatched = "Contigs table ({}): 3 of 10 row(s) match no contig of the map and are ignored, e.g. 'c1'."
    mc_logger = logging.getLogger("metacooler")
    with mtmc._metacooler_logging(), mtmc._rows_of_other_contigs_ignored(contig_data):
        mc_logger.warning(unmatched.format(contig_data))
        mc_logger.warning(unmatched.format(other))
        mc_logger.warning(f"Contigs table ({contig_data}): 1 contig(s) are described by more than one row")
        mc_logger.error(unmatched.format(contig_data))
    with mtmc._metacooler_logging():
        mc_logger.warning(unmatched.format(contig_data))
    assert logged() == [
        ("WARNING", unmatched.format(other)),
        ("WARNING", f"Contigs table ({contig_data}): 1 contig(s) are described by more than one row"),
        ("ERROR", unmatched.format(contig_data)),
        ("WARNING", unmatched.format(contig_data)),
    ]


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
def test_map_of_no_bin_from_a_map(source, metagenome, binning, tmp_path):
    """A run with no final bin: every contig of the map in __unbinned__, its
    only genome, with no column of a bin summary, even when the input map
    is that of an earlier binning, whose bin summary columns would be stale."""
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
    contigs = metacooler.MetaCooler(str(out)).contigs
    assert list(contigs.columns) == ["name", "length", "genome_id", "GC_content", "Hit", "taxon", "Restriction_site"]
    # Read back from the file, as metacooler does: the float parser of
    # pandas 2 can differ from the in-memory value in the last place.
    measured = pd.read_csv(contig_data, sep="\t", index_col="Name").loc[contigs.name]
    for column in ["GC_content", "Hit", "Restriction_site"]:
        assert contigs[column].tolist() == measured[column].tolist(), column
    assert contigs.set_index("name").taxon.to_dict() == {
        name: f"t_{i % 3}" for i, name in enumerate(metagenome.truth.Name)
    }  # fmt: skip


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
    assert sorted(contigs.columns) == sorted(["name", "length", "genome_id", *CONTIG_COLUMNS, "taxon"])


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
    """metacooler.reassign, failing once it has started to write its partial
    file next to the output."""
    Path(f"{output}.partial").write_text("partial")
    raise OSError("No space left on device")


def _fail_in_validate(path, full=True):
    raise ValueError("pixels out of order")


@pytest.mark.parametrize("failure", ["reassign", "validate", "in metacooler"])
@pytest.mark.parametrize("how", ["path", "link to the outfile"])
def test_input_map_kept_when_the_map_of_the_bins_fails(how, failure, metagenome, binning, tmp_path, monkeypatch, logged):
    """When the map of the bins cannot be written or does not validate, the
    input map (the outfile) is left byte for byte as it was, and no
    temporary file stays next to it. When metacooler has left its partial
    file, which it logs as left for inspection, the error says that it was
    removed."""
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
    expected = message
    if failure != "validate":
        # The temporary file, next to the outfile, and the partial file that
        # metacooler (or the stub) writes next to it.
        staged = "(" + re.escape(str(out / ".metator.metacool.")) + r"\w+\.tmp)"
        expected = (
            rf"{message}\. The temporary files of the map, {staged} and metacooler's partial file \1\.partial, were removed\."
        )
    with pytest.raises(mtmc.MetacoolError, match=f"^{expected}$") as raised:
        mtmc.reassign_metacool(
            str(metacool), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True
        )
    assert str(raised.value.__cause__) == message
    assert outfile.read_bytes() == before
    assert sorted(path.name for path in out.iterdir()) == ["metator.metacool"]
    if failure == "in metacooler":
        # metacooler has logged that it left the partial file it was writing.
        staged = re.match(expected, str(raised.value)).group(1)
        assert (
            "ERROR",
            f"Creation failed. The partial file was left at {staged}.partial for inspection; {staged} was not written.",
        ) in logged()


def _fail_unexpectedly(path, full=True):
    raise RuntimeError("unexpected failure")


def test_partial_file_removed_on_an_unexpected_error(metagenome, binning, tmp_path, monkeypatch, logged):
    """An error of another type than those metacooler is expected to raise,
    within metacooler reassign once it has started its partial file (in its
    validation, here), is raised as it is: a warning says that the partial
    file, which metacooler logs as left for inspection, was removed. The
    input map (the outfile) is left as it was, and no other file is left."""
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
    left = [message for level, message in messages if level == "ERROR" and "partial file was left" in message]
    (staged,) = [
        re.match(r"Creation failed\. The partial file was left at (\S+)\.partial for inspection", message)[1] for message in left
    ]
    assert Path(staged).parent == out and Path(staged).name.startswith(".metator.metacool.")
    removed = f"The temporary files of the map, {staged} and metacooler's partial file {staged}.partial, were removed."
    assert messages[-1] == ("WARNING", removed)


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
