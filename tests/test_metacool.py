# Tests of metator.metacool: the contact map of a binning as a .metacool
# file, on the mock metagenome of tests/synthetic.py. Its genome contigs are
# binned by genome and its MGE contigs left unbinned.

import logging
import re
import shlex
import shutil
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

PAIRS_HEADER = "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"


def _bin_names(metagenome):
    """Final bin of each contig: one bin per genome, none for the MGEs."""
    names = {genome.name: f"metator_{i:05d}_00000" for i, genome in enumerate(metagenome.genomes, start=1)}
    return {contig.name: names.get(contig.molecule, "ND") for contig in metagenome.contigs}


def write_contig_data(metagenome, path, final_bin=None):
    """A contig_data_final.txt file, laid out as metaTOR writes it (an
    unnamed index first, Final_bin last)."""
    final_bin = final_bin or _bin_names(metagenome)
    truth = metagenome.truth
    table = pd.DataFrame(
        {
            "ID": range(1, len(truth) + 1),
            "Name": truth.Name.to_numpy(),
            "Size": truth.Length.to_numpy(),
            "GC_content": truth.GC.to_numpy(),
            "Hit": 1,
            "Recursive_bin_contigs": "-",
            "Final_bin": truth.Name.map(final_bin).to_numpy(),
        }
    )
    table.to_csv(path, sep="\t", index=True)
    return path


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

    metadata = clr.info["metadata"]
    assert metadata["metator"] == {
        "version": __version__,
        "command": "metator metacool",
        "parameters": {"threads": 1},
    }
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


def test_existing_outfile(binning, tmp_path):
    out = tmp_path / "map.metacool"
    out.write_text("old")
    args = ([str(binning["pairs"])], str(binning["contig_data"]), str(out))
    with pytest.raises(mtmc.MetacoolError, match=re.escape("already exists. Use --force")):
        mtmc.write_metacool(*args)
    assert out.read_text() == "old"
    mtmc.write_metacool(*args, force=True)
    metacooler.validate(str(out))


def test_missing_input_files(binning, tmp_path):
    out = tmp_path / "map.metacool"
    with pytest.raises(mtmc.MetacoolError, match=re.escape("not found: missing.pairs, missing_summary.txt")):
        mtmc.write_metacool(
            [str(binning["pairs"]), "missing.pairs"], str(binning["contig_data"]), str(out), bin_summary="missing_summary.txt"
        )
    assert not out.exists()


def test_contig_data_without_final_bins(metagenome_network, binning, tmp_path):
    """The contig data of the network step has no Final_bin column."""
    out = tmp_path / "map.metacool"
    with pytest.raises(mtmc.MetacoolError, match=re.escape("no Final_bin column: give the contig_data_final.txt")):
        mtmc.write_metacool(str(binning["pairs"]), str(metagenome_network["contig_data"]), str(out))
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


def test_read_contig_bins(tmp_path):
    table = tmp_path / "contig_data_final.txt"
    table.write_text("\tID\tName\tSize\tFinal_bin\n0\t1\t01\t100\tb1\n1\t2\tc2\t200\tND\n2\t3\tc3\t300\t-\n3\t4\tc4\t400\t\n")
    bins, final_bin = mtmc.read_contig_bins(str(table))
    assert bins.to_dict("list") == {"chrom": ["01", "c2", "c3", "c4"], "start": [0] * 4, "end": [100, 200, 300, 400]}
    assert final_bin == {"01": "b1"}


def test_write_genomes_table(tmp_path):
    """The bin names become genome_name; the values are copied verbatim."""
    summary = tmp_path / "bin_summary.txt"
    summary.write_text(
        "\tLength\tPresent Markers\tparent\tHiC_abundance\nb1\t100.00\t3.00\t\t0.0597\nb2\t20.50\t1.00\tb0\t1.0000\n"
    )
    mtmc.write_genomes_table(str(summary), str(tmp_path / "genomes.tsv"))
    assert (tmp_path / "genomes.tsv").read_text() == (
        "genome_name\tLength\tPresent Markers\tparent\tHiC_abundance\nb1\t100.00\t3.00\t\t0.0597\nb2\t20.50\t1.00\tb0\t1.0000\n"
    )


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


ERRORS = {
    "missing pairs file": (["missing.pairs"], "Input file(s) not found: missing.pairs."),
    "contig data without Final_bin": ([], "has no Final_bin column"),
    "non-ASCII contig names": ([], "A .metacool file needs ASCII contig names"),
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
