# Tests of metator markers: the marker genes of the contigs of a map, found as
# miComplete finds those of a bin, and the four numbers of miComplete
# computed again from them. miComplete itself, run as metaTOR runs it
# (python -m metator._micomplete, with pyrodigal and pyhmmer), is the oracle:
# on the bins of tests_data/micomplete, and on the genomes of a small map made
# of their contigs, written in the order of the map as metaTOR writes the
# FASTA file of a genome.

import itertools
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import metacooler
import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO

import metator.clusters as mtcl
import metator.markers as mtmk
import metator.metacool as mtmc
import metator.validation as mtv
from metator.version import __version__
from tests.hdf5 import walk
from tests.logs import logs

BINS = Path(__file__).resolve().parents[1] / "tests_data" / "micomplete"

# pyrodigal warns that a bin under 100 kb is short to train on, as prodigal
# does; contig_markers silences it, not group_markers, which the tests call.
pytestmark = pytest.mark.filterwarnings("ignore:sequence should be at least:UserWarning")

#: Two contigs of metator_00002_00000, 110 kb: miComplete sees a hit of
#: Cons_hypoth95 on Contig_413455 and drops it (its best domain is too weak),
#: and no other protein of the bin is one: the marker counts in the weighted
#: completeness of the bin, not in its completeness.
DROPPED = ("Contig_352019", "Contig_413455")

#: The genomes of the map of the tests (see marker_map): two bins of
#: tests_data, which pyrodigal trains on (38 and 37 kb); binC, two contigs of
#: metator_00009_00000 too short to train on (2 kb); and three of its contigs
#: in no genome. A contig of each of the last two holds a marker of Arch131.
MAP_GENOMES = {"binA": "metator_00063_00000", "binB": "metator_00082_00000"}
SMALL = ("Contig_23233", "Contig_241358")
UNBINNED = ("Contig_535749", "Contig_5196", "Contig_490599")


def _records(name):
    """The records of a bin of tests_data/micomplete."""
    return list(SeqIO.parse(BINS / f"{name}.fa", "fasta"))


def _write(records, path):
    SeqIO.write(records, path, "fasta")
    return path


def _first_contigs(length):
    """The first contigs of metator_00009_00000, as long as they total less than
    length bp (as tests/test_micomplete_shim.py makes them)."""
    contigs, total = [], 0
    for record in _records("metator_00009_00000"):
        if total + len(record) >= length:
            break
        contigs.append(record)
        total += len(record)
    return contigs


def _micomplete(fastas, directory):
    """miComplete's results on the FASTA files fastas, run as metaTOR runs it, with
    both marker sets: by set, then Markers (metaTOR's merge of the two, see
    metator.validation.merge_micomplete), a table of the rows of miComplete by
    bin name."""
    seq_tab = directory / "seq.tab"
    seq_tab.write_text("".join(f"{fasta}\tfna\n" for fasta in fastas))
    tables = {}
    for name in mtmk.MARKER_SETS:
        outfile = directory / f"{name}.tsv"
        cmd = [sys.executable, "-m", "metator._micomplete", str(seq_tab), "--hmms", name, "--weights", name,
               "--threads", "2", "--outfile", str(outfile)]  # fmt: skip
        subprocess.run(cmd, cwd=directory, check=True)
        tables[name] = pd.read_csv(outfile, sep="\t", comment="#", index_col=0)
    merged = directory / "merged.tsv"
    mtv.merge_micomplete(str(directory / "Bact105.tsv"), str(directory / "Arch131.tsv"), str(merged))
    tables["Markers"] = pd.read_csv(merged, sep="\t", comment="#", index_col=0)
    return tables


def marker_map(directory):
    """A map of the contigs of tests_data/micomplete (one bin per contig), with its
    assembly: the genomes of MAP_GENOMES, each of the contigs of its bin;
    binC, the contigs SMALL (under 20 kb); and the contigs UNBINNED in no
    genome. Returns the map and the assembly."""
    genome_of, records = {}, []
    for genome, name in MAP_GENOMES.items():
        for record in _records(name):
            genome_of[record.id] = genome
            records.append(record)
    others = {record.id: record for record in _records("metator_00009_00000")}
    for genome, names in [("binC", SMALL), ("ND", UNBINNED)]:
        for name in names:
            genome_of[name] = genome
            records.append(others[name])
    # The assembly in another order than the map's.
    assembly = _write(records[::-1], directory / "assembly.fa")
    names = [record.id for record in records]
    pairs = directory / "contigs.pairs"
    pairs.write_text(
        "## pairs format v1.0\n"
        + "".join(f"r{i}\t{first}\t10\t{second}\t20\t+\t-\n" for i, (first, second) in enumerate(itertools.pairwise(names)))
    )
    path = directory / "contigs.metacool"
    bins = pd.DataFrame({"chrom": names, "start": 0, "end": [len(record) for record in records]})
    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=genome_of, ps=False, keep_partial=False)
    pairs.unlink()
    return path, assembly


@pytest.fixture(scope="module")
def oracle(tmp_path_factory):
    """The bins checked against miComplete, the map of the tests (see marker_map)
    and miComplete's results on them all (see _micomplete): the bins of
    tests_data/micomplete; dropped, the two contigs of DROPPED; too_short and
    just_enough, the first contigs of metator_00009_00000 under 15 and 20 kb
    (prodigal trains on the second, its contigs joined with linkers); and the
    genomes binA and binB of the map, written in the order of the map, as
    metacooler fasta writes the FASTA file of a genome."""
    directory = tmp_path_factory.mktemp("markers_oracle")
    fastas = {path.stem: path for path in sorted(BINS.glob("*.fa"))}
    records = {record.id: record for record in _records("metator_00002_00000") if record.id in DROPPED}
    fastas["dropped"] = _write([records[name] for name in DROPPED], directory / "dropped.fa")
    fastas["too_short"] = _write(_first_contigs(15_000), directory / "too_short.fa")
    fastas["just_enough"] = _write(_first_contigs(20_000), directory / "just_enough.fa")
    metacool, assembly = marker_map(directory)
    genomes = directory / "genomes"
    metacooler.write_fasta(str(metacool), str(assembly), str(genomes), genomes=list(MAP_GENOMES), force=True)
    fastas.update({name: genomes / f"{name}.fa" for name in MAP_GENOMES})
    return {"fastas": fastas, "map": metacool, "assembly": assembly, "micomplete": _micomplete(list(fastas.values()), directory)}


def _contigs(fasta):
    """The name and the sequence of each record of a FASTA file, in order, as
    metaTOR's miComplete reads them."""
    return [(record.id, bytes(record.seq)) for record in SeqIO.parse(fasta, "fasta")]


def _row(table, name):
    """The four numbers of miComplete in its row of a bin (Completeness...)."""
    row = table.loc[name]
    return mtmk.Quality(
        float(row["Completeness"]),
        float(row["Redundancy"]),
        float(row["Weighted completeness"]),
        float(row["Weighted redundancy"]),
    )


# The markers of a bin, summed over its contigs, give the numbers of miComplete.


@pytest.mark.parametrize("name", ["metator_00002_00000", "metator_00009_00000", "metator_00063_00000",
                                  "metator_00082_00000", "dropped", "too_short", "just_enough"])  # fmt: skip
def test_the_markers_of_a_bin_give_the_numbers_of_micomplete(name, oracle):
    """The markers of each contig of a bin, found as miComplete finds those of the
    bin (trained on its contigs, in the order of its FASTA file, the E-values
    counted over its proteins), written as the map stores them and read
    back, then summed, give the four numbers that miComplete gives the bin, in
    each marker set, and the set that metaTOR keeps; and as many proteins as
    miComplete finds (none in too_short)."""
    found = mtmk.group_markers(_contigs(oracle["fastas"][name]))
    tables = oracle["micomplete"]
    for markers in mtmk.MARKER_SETS:
        marker_set = mtmk.marker_set(markers)
        by_contig = found.markers.get(markers, {})
        counts = mtmk.sum_markers(mtmk.decode_markers(mtmk.encode_markers(contig)) for contig in by_contig.values())
        assert mtmk.marker_quality(counts, marker_set.markers, marker_set.weights) == _row(tables[markers], name), markers
    texts = {markers: [mtmk.encode_markers(contig) for contig in found.markers.get(markers, {}).values()]
             for markers in mtmk.MARKER_SETS}  # fmt: skip
    kingdom, quality = mtmk.quality_of(texts)
    assert (kingdom, quality) == (tables["Markers"].loc[name, "Markers"], _row(tables["Markers"], name))
    assert (found.proteins or 0) == tables["Bact105"].loc[name, "CDs"]


def test_a_marker_seen_and_dropped_counts_in_the_weighted_completeness(oracle):
    """In dropped, miComplete sees Cons_hypoth95 (a hit in the table of hmmsearch)
    but drops it: the marker counts in its weighted completeness, as miComplete
    counts the markers that it sees, but not in its completeness. The contig
    holds it as Cons_hypoth95*0."""
    found = mtmk.group_markers(_contigs(oracle["fastas"]["dropped"]), sets=["Bact105"])
    by_contig = found.markers["Bact105"]
    assert by_contig["Contig_413455"]["Cons_hypoth95"] == 0
    assert "Cons_hypoth95*0" in mtmk.encode_markers(by_contig["Contig_413455"]).split(",")
    marker_set = mtmk.marker_set("Bact105")
    counts = mtmk.sum_markers(by_contig.values())
    seen = mtmk.marker_quality(counts, marker_set.markers, marker_set.weights)
    present = mtmk.marker_quality({name: n for name, n in counts.items() if n}, marker_set.markers, marker_set.weights)
    weight = dict(zip(marker_set.markers, marker_set.weights, strict=True))["Cons_hypoth95"]
    assert seen.completeness == present.completeness and seen.redundancy == present.redundancy
    assert seen.weighted_completeness == pytest.approx(present.weighted_completeness + weight, abs=1e-4)
    assert seen == _row(oracle["micomplete"]["Bact105"], "dropped")


def test_read_hits_as_micomplete():
    """The hits of a table of hmmsearch, as miComplete reads them: the numbers as
    the table writes them, a hit dropped when its score times the bias cutoff
    is no more than its bias (at equality too), or its best domain's E-value
    exceeds that of the protein by more than the domain cutoff; the header
    names no marker, and a line "#" (hmmsearch's footer) ends the table."""
    header = (
        "#                                                               --- full sequence ---- --- best 1 domain ----\n"
        "# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias\n"
        "#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----\n"
    )
    rows = [
        "c_1 - SecY PF00344.15 1.2e-50 170.0 0.4 1.5e-50 169.0 0.4 1.0 1 0 0 1 1 1 1 # 1 # 900 # 1 # ID=1_1",
        # 10 * 0.3 = 3.0, above the bias 2.9: kept.
        "c_2 - SecY PF00344.15 1e-20 10.0 2.9 1e-20 9.0 2.9 1.0 1 0 0 1 1 1 1 # 1 # 900 # 1 # ID=1_2",
        # 20 * 0.3 = 6.0, no more than the bias 6.0: dropped.
        "c_3 - RecR PF02132.14 1e-20 20.0 6.0 1e-20 19.0 6.0 1.0 1 0 0 1 1 1 1 # 1 # 900 # 1 # ID=1_3",
        # The best domain 2e-5 above the protein: dropped.
        "c_4 - RecR PF02132.14 1e-12 50.0 0.1 2e-05 48.0 0.1 1.0 1 0 0 1 1 1 1 # 1 # 900 # 1 # ID=1_4",
        "#",
        "c_5 - RecR PF02132.14 1e-50 50.0 0.1 1e-50 48.0 0.1 1.0 1 0 0 1 1 1 1 # 1 # 900 # 1 # ID=1_5",
    ]
    assert mtmk.read_hits(header + "\n".join(rows) + "\n") == [
        ("c_1", "SecY", True),
        ("c_2", "SecY", True),
        ("c_3", "RecR", False),
        ("c_4", "RecR", False),
    ]
    assert mtmk.read_hits(header + rows[3], domain_cutoff=1e-4) == [("c_4", "RecR", True)]


def test_a_group_without_a_protein():
    """No contig, or no protein: no gene to train on, and no hit to search."""
    assert mtmk.find_proteins([]) is None
    assert mtmk.search_table(mtmk.Proteins("", []), mtmk.marker_set("Bact105")) == ""
    assert mtmk.group_markers([]) == mtmk.GroupMarkers(None)


def test_marker_sets_need_micomplete(monkeypatch):
    """The marker sets are miComplete's files: without miComplete, none."""
    import importlib.util

    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
    with pytest.raises(ImportError, match="miComplete is not installed"):
        mtmk.marker_set.__wrapped__("Bact105")


def test_count_markers():
    """The markers of each contig, from the hits in its proteins: a protein counts
    for each marker that it is kept for, and a marker that a contig only has
    dropped hits of counts 0 there."""
    proteins = mtmk.Proteins(">a_1 # x\nMK*\n>a_2 # x\nMR\n>b_1 # x\nMA*\n", ["a", "a", "b"])
    hits = [("a_1", "SecY", True), ("a_2", "SecY", True), ("a_1", "RecR", True), ("b_1", "RecR", False)]
    assert mtmk.count_markers(proteins, hits) == {"a": {"SecY": 2, "RecR": 1}, "b": {"RecR": 0}}


# The markers as text, and the four numbers.


def test_encode_and_decode_markers():
    """Sorted by name, in the order of their bytes (upper case first),
    comma-separated, a count but 1 after *: written and read back. An empty
    text, a missing value or None is no marker."""
    counts = {"tRNA-synt_1d": 1, "SecY": 2, "Ribosomal_L22": 1, "RecR": 0, "TruB_N": 3}
    text = mtmk.encode_markers(counts)
    assert text == "RecR*0,Ribosomal_L22,SecY*2,TruB_N*3,tRNA-synt_1d"
    assert mtmk.decode_markers(text) == counts
    assert mtmk.encode_markers({}) == ""
    for none in ("", None, float("nan"), np.nan):
        assert mtmk.decode_markers(none) == {}
    assert mtmk.sum_markers([{"SecY": 1, "RecR": 0}, {"SecY": 1}, {"RecR": 0, "TruB_N": 1}, {}]) == {
        "SecY": 2,
        "RecR": 0,
        "TruB_N": 1,
    }


#: A marker set of four markers, for the numbers of marker_quality.
MARKERS = ("A", "B", "C", "D")
WEIGHTS = (0.1, 0.2, 0.3, 0.4)


@pytest.mark.parametrize(
    "counts, expected",
    [
        ({}, (0.0, 0.0, 0.0, 0.0)),
        # Seen, every hit dropped: no present marker, nothing counts.
        ({"A": 0}, (0.0, 0.0, 0.0, 0.0)),
        ({"A": 1, "B": 1}, (0.5, 1.0, 0.3, 1.0)),
        # B twice: redundancy 3/2; its weight once in the redundancy.
        ({"A": 1, "B": 2}, (0.5, 1.5, 0.3, round((0.2 + 0.1 + 0.2) / (0.1 + 0.2), 4))),
        # C seen and dropped: in the weighted completeness, nowhere else.
        ({"A": 1, "C": 0}, (0.25, 1.0, 0.4, 1.0)),
        ({"A": 3, "B": 1, "C": 2, "D": 1}, (1.0, 1.75, 1.0, round((0.1 + 0.3 + 1.0) / 1.0, 4))),
    ],
)
def test_marker_quality(counts, expected):
    assert mtmk.marker_quality(counts, MARKERS, WEIGHTS) == expected


def test_marker_quality_rounds_as_micomplete():
    """Rounded to 4 decimals, a weighted completeness of less than 0.00005 written
    0.0001, and the weighted redundancy computed from the weighted completeness
    before it is rounded. An exact tie is rounded to the even digit, as
    Python's round does: 33 copies of 32 markers, 1.03125, give 1.0312."""
    markers = [f"M{number}" for number in range(32)]
    counts = {name: 2 if name == "M0" else 1 for name in markers}
    assert mtmk.marker_quality(counts, markers, [1.0] * 32) == (1.0, 1.0312, 32.0, 1.0312)
    quality = mtmk.marker_quality({"A": 2, "B": 1}, ("A", "B", "C"), (1e-05, 2e-05, 0.5))
    assert quality == (round(2 / 3, 4), 1.5, 0.0001, round((1e-05 + 3e-05) / 3e-05, 4))
    with pytest.raises(ValueError, match="Markers of no marker set of 4 markers: X"):
        mtmk.marker_quality({"X": 1}, MARKERS, WEIGHTS)


@pytest.mark.parametrize("archaea, kept", [(0.5, "Bacteria"), (0.6, "Archaea"), (0.4, "Bacteria")])
def test_kept_quality(archaea, kept):
    """The set of the higher weighted completeness, Bact105 on a tie, as
    metator.validation.merge_micomplete keeps miComplete's results."""
    bacteria = mtmk.Quality(0.5, 1.0, 0.5, 1.0)
    other = mtmk.Quality(0.4, 1.1, archaea, 1.2)
    assert mtmk.kept_quality({"Bact105": bacteria, "Arch131": other}) == (
        kept,
        bacteria if kept == "Bacteria" else other,
    )


def test_marker_sets():
    """The sets of miComplete: their markers, as miComplete names them, in the
    order of its HMMs, with their weights; 35 markers in both."""
    bacteria, archaea = mtmk.marker_set("Bact105"), mtmk.marker_set("Arch131")
    assert (len(bacteria.markers), len(archaea.markers)) == (105, 131)
    assert len(set(bacteria.markers) & set(archaea.markers)) == 35
    assert bacteria.markers[:2] == ("CoaE", "Cons_hypoth95")
    assert dict(zip(bacteria.markers, bacteria.weights, strict=True))["Ribosomal_L22"] == 1.0227835678392199e-05
    assert all(len(name) <= 15 and name.isascii() for name in bacteria.markers + archaea.markers)
    entry = mtmk.entry_sets()
    assert list(entry) == ["Bact105", "Arch131"]
    assert entry["Bact105"] == {"column": "markers_Bact105", "markers": list(bacteria.markers), "weights": list(bacteria.weights)}


def test_metagenome_z():
    """The median number of proteins of the trained genomes; without one, the
    proteins searched; at least one."""
    assert mtmk.metagenome_z({"a": 10, "b": 30, "c": 25}, 7) == 25.0
    assert mtmk.metagenome_z({"a": 10, "b": 31}, 7) == 20.5
    assert mtmk.metagenome_z({}, 7) == 7.0
    assert mtmk.metagenome_z({}, 0) == 1.0


# metator markers: the markers stored in the map.


def _stored(metacool):
    """The markers that the map stores: the columns of its contigs table, by
    contig, and the markers entry of its metator metadata."""
    mc = metacooler.MetaCooler(str(metacool))
    table = mc.contigs.set_index(mc.contigs["name"].astype(str))
    return table[list(mtmc.MARKER_COLUMNS)], mc.metadata["metator"]["markers"]


def _genome_texts(metacool, genome):
    """The markers of the contigs of a genome of the map, by set, as stored."""
    mc = metacooler.MetaCooler(str(metacool))
    names = mc.genomes["name"].astype(str).to_numpy()[mc.contigs["genome_id"].to_numpy()]
    contigs = mc.contigs[names == genome]
    return {name: list(contigs[mtmk.marker_column(name)]) for name in mtmk.MARKER_SETS}


@pytest.fixture
def marker_run(oracle, tmp_path):
    """A private copy of the map of the tests, its assembly, and the output
    directory of a run (tmp_path/out)."""
    metacool = tmp_path / "contigs.metacool"
    shutil.copyfile(oracle["map"], metacool)
    return metacool, oracle["assembly"], tmp_path / "out"


@pytest.mark.parametrize("output", ["the output directory", "--inplace", "an entry"])
def test_markers_stores_the_markers_of_the_contigs(output, marker_run, oracle, tmp_path, run_cli):
    """markers MAP stores the markers of each contig in the map of its output
    directory (MAP itself with --inplace; for an entry of a collection, a map of
    its own): summed over a genome that pyrodigal trains on, they give the
    numbers that miComplete gives the FASTA file of the genome, in the order of
    the map; the contigs of binC (under 20 kb) and of no genome get theirs in
    metagenome mode, searched as in a genome of the median number of proteins
    of the others. The run is recorded under the markers entry of the metator
    metadata, and everything else is kept as it was."""
    metacool, assembly, out = marker_run
    before = walk(metacool)
    unchanged = metacool.read_bytes()
    source = metacool
    if output == "an entry":
        source = f"{tmp_path / 'maps.mmetacool'}::entries/lib1"
        metacooler.create_mmetacool(tmp_path / "maps.mmetacool", {"lib1": metacool})
    options = ["--inplace"] if output == "--inplace" else []
    written = metacool if output == "--inplace" else out / ("lib1.metacool" if output == "an entry" else metacool.name)
    arguments = [str(source), "-a", str(assembly), "-o", str(out), "-t", "2", "-T", str(tmp_path / "tmp"), *options]
    run_cli("markers", *arguments)
    metacooler.validate(str(written), full=True)
    tables = oracle["micomplete"]
    for genome in MAP_GENOMES:
        kingdom, quality = mtmk.quality_of(_genome_texts(written, genome))
        assert (kingdom, quality) == (tables["Markers"].loc[genome, "Markers"], _row(tables["Markers"], genome)), genome
    table, entry = _stored(written)
    # The markers of Arch131 of the contigs in metagenome mode, as a bin of
    # metator_00009_00000 has them.
    assert list(table.loc[[*SMALL, *UNBINNED], "markers_Arch131"]) == ["", "Adenylsucc_synt", "", "", "ATP-synt_F"]
    assert (table.loc[[*SMALL, *UNBINNED], "markers_Bact105"] == "").all()
    assert entry["version"] == __version__
    assert entry["command"] == "metator markers " + " ".join(arguments)
    assert entry["parameters"]["evalue"] == 4e-10 and entry["parameters"]["threads"] == 2
    assert entry["thresholds"] == {"evalue": 4e-10, "bias": 0.3, "domain_cutoff": 1e-05}
    proteins = [int(tables["Bact105"].loc[genome, "CDs"]) for genome in MAP_GENOMES]
    genes = entry["genes"]
    assert (genes["trained_genomes"], genes["metagenome_contigs"], genes["min_training"]) == (2, 5, 20_000)
    # How the genes were found, in words: the rules of the two modes.
    assert {mode: genes[mode] for mode in ("single", "metagenome")} == mtmk.GENE_RULES
    assert genes["metagenome_z"] == (proteins[0] + proteins[1]) / 2
    assert entry["sets"] == mtmk.entry_sets() and entry["kept"] == {"by": "weighted_completeness", "ties": "Bact105"}
    assert entry["read"] == {"genomes": ["name"]}
    # Everything else as it was: the contigs table but its two columns, and the
    # metadata but the markers entry and metacooler's record of the annotation.
    after = walk(written)
    for name in mtmc.MARKER_COLUMNS:
        del after[f"/contigs/{name}"]
    kind, keys, attributes = after["/contigs"]
    after["/contigs"] = (kind, [key for key in keys if key not in mtmc.MARKER_COLUMNS], attributes)
    for state in (before, after):
        kind, keys, attributes = state["/"]
        state["/"] = (kind, keys, [(key, dtype, value) for key, dtype, value in attributes if key != "metadata"])
    assert after == before
    metadata = metacooler.MetaCooler(str(written)).metadata
    assert set(metadata["metator"]) == {"markers"}
    if output != "--inplace":
        assert metacool.read_bytes() == unchanged
    (log,) = logs(out, "markers")
    logged = log.read_text()
    assert "The genes of 2 genome(s) were found as miComplete finds those of a bin" in logged
    assert "those of 5 contig(s), in no genome or in 1 genome(s) too short to train on (under 20 kb)" in logged
    assert "Contact map with the marker genes of its contigs written after" in logged
    assert not list(tmp_path.glob("*.lock")) and not list(out.glob("*.lock"))
    assert not (tmp_path / "tmp").exists()


def test_markers_in_threads_or_not(marker_run, tmp_path):
    """The markers found by one thread are those found by several."""
    metacool, assembly, _ = marker_run
    mc = metacooler.MetaCooler(str(metacool))
    sequences = mtmk.read_contigs(str(metacool), mc, str(assembly))
    one = mtmk.contig_markers(str(metacool), mc, sequences, threads=1)
    several = mtmk.contig_markers(str(metacool), mc, sequences, threads=3)
    pd.testing.assert_frame_equal(one.table, several.table)
    assert (one.trained, one.small, one.meta, one.z) == (several.trained, several.small, several.meta, several.z)
    assert list(one.trained) == ["binA", "binB"] and one.small == ["binC"]


def test_markers_of_a_map_without_genomes(oracle, tmp_path):
    """On a map whose contigs are all in __unbinned__ (before any binning), every
    contig is in metagenome mode, its E-values counted over the proteins of
    them all, searched as one bin: there is no genome to take the median of."""
    metacool = tmp_path / "unbinned.metacool"
    mc = metacooler.MetaCooler(str(oracle["map"]))
    names = mc.contigs["name"].astype(str).tolist()
    bins = pd.DataFrame({"chrom": names, "start": 0, "end": mc.contigs["length"].to_numpy()})
    pairs = tmp_path / "contigs.pairs"
    pairs.write_text(f"## pairs format v1.0\nr1\t{names[0]}\t10\t{names[1]}\t20\t+\t-\n")
    metacooler.from_pairs(str(metacool), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False, keep_partial=False)
    unbinned = metacooler.MetaCooler(str(metacool))
    sequences = mtmk.read_contigs(str(metacool), unbinned, str(oracle["assembly"]))
    found = mtmk.contig_markers(str(metacool), unbinned, sequences, threads=2)
    proteins = mtmk.find_proteins([(name, sequences[name]) for name in names], meta=True)
    assert (found.trained, found.small, found.meta, found.z) == ({}, [], len(names), float(len(proteins)))
    assert (found.table[mtmc.MARKER_COLUMNS[0]] != "").any()


def test_markers_rerun_replaces_its_markers(marker_run, tmp_path, run_cli, logged):
    """A second run replaces the markers of the first, with its thresholds: in
    place, with a warning before any work; -N keeps the table of the markers
    in the temporary directory of the run."""
    metacool, assembly, out = marker_run
    run_cli("markers", metacool, "-a", assembly, "--inplace", "-o", out, "-T", tmp_path / "tmp")
    first, _ = _stored(metacool)
    run_cli("markers", metacool, "-a", assembly, "--inplace", "-o", out, "-T", tmp_path / "tmp", "-N", "--evalue", "1e-30")
    warning = (
        f"{metacool} holds the marker genes of its contigs (metator markers), which this run replaces, in place: leave "
        f"out --inplace to write the map of this run to another file, and leave {metacool} as it is."
    )
    assert ("WARNING", warning) in logged()
    second, entry = _stored(metacool)
    assert entry["thresholds"]["evalue"] == 1e-30
    # Fewer markers with a lower E-value, and none that the first run lacked.
    for column in mtmc.MARKER_COLUMNS:
        for before, after in zip(first[column], second[column], strict=True):
            assert set(mtmk.decode_markers(after)) <= set(mtmk.decode_markers(before))
    assert (first != second).any(axis=None)
    (kept,) = (tmp_path / "tmp").iterdir()
    table = pd.read_csv(kept / mtmc.MARKERS_FILE, sep="\t", keep_default_na=False)
    assert list(table.columns) == ["Name", *mtmc.MARKER_COLUMNS]
    stored = second.reindex(table["Name"]).fillna("")
    assert list(table[mtmc.MARKER_COLUMNS[0]]) == list(stored[mtmc.MARKER_COLUMNS[0]])


def test_a_set_without_any_marker_is_stored_as_empty_texts(marker_run, tmp_path, run_cli, logged):
    """With an E-value that no hit reaches, no contig has a marker: each column
    holds an empty text on every contig, as text still, which metacooler would
    store as numbers, all missing: a row more of the table of the markers
    (which -N keeps), naming no contig, keeps the columns text, and metacooler
    leaves it out without a word. A run that finds markers again stores them."""
    metacool, assembly, out = marker_run
    tmp = tmp_path / "tmp"
    run_cli("markers", metacool, "-a", assembly, "--inplace", "-o", out, "-T", tmp, "-N", "--evalue", "1e-300")
    table, entry = _stored(metacool)
    assert entry["thresholds"]["evalue"] == 1e-300
    assert all(table[column].tolist() == [""] * len(table) for column in mtmc.MARKER_COLUMNS)
    with h5py.File(metacool, "r") as handle:
        assert [handle["contigs"][column].dtype.kind for column in mtmc.MARKER_COLUMNS] == ["S", "S"]
    assert mtmk.quality_of(_genome_texts(metacool, "binA")) == ("Bacteria", (0.0, 0.0, 0.0, 0.0))
    (kept,) = tmp.iterdir()
    rows = pd.read_csv(kept / mtmc.MARKERS_FILE, sep="\t", keep_default_na=False)
    assert list(rows["Name"]) == [*table.index, mtmc.NO_CONTIG]
    assert list(rows.iloc[-1]) == [mtmc.NO_CONTIG, mtmc.NO_CONTIG_TEXT, mtmc.NO_CONTIG_TEXT]
    assert not [message for level, message in logged() if level in ("WARNING", "ERROR")]
    run_cli("markers", metacool, "-a", assembly, "--inplace", "-o", out)
    table, _ = _stored(metacool)
    assert (table[mtmc.MARKER_COLUMNS[0]] != "").any()


def test_the_row_that_names_no_contig():
    """The row of the table of the markers that keeps their columns text names
    no contig of the map, nor one by its name in the input of a collected map
    (the name after <genome>__), which metacooler annotate matches too."""
    assert mtmc._no_contig(["c1", "c2"]) == "(no contig)"
    assert mtmc._no_contig(["c1", "(no contig)"]) == "_(no contig)"
    assert mtmc._no_contig(["g__(no contig)", "g___(no contig)"]) == "__(no contig)"


@pytest.mark.parametrize("case", ["a contig missing", "another length", "two records of a contig", "not a FASTA file"])
def test_markers_refuses_an_assembly_that_is_not_that_of_the_map(case, marker_run, tmp_path, run_cli, logged):
    """An assembly that lacks a contig of the map, gives one another length, or
    names a contig by several records (or that is no FASTA file) is refused
    before anything is written, the reason logged as an error."""
    metacool, assembly, out = marker_run
    records = list(SeqIO.parse(assembly, "fasta"))
    name = records[0].id
    if case == "a contig missing":
        records = records[1:]
    elif case == "another length":
        records[0] = records[0][:-1]
    elif case == "two records of a contig":
        records.append(records[0])
    other = tmp_path / "other.fa"
    if case == "not a FASTA file":
        other.write_text("not a FASTA file\n")
    else:
        _write(records, other)
    before = metacool.read_bytes()
    with pytest.raises((mtmc.MetacoolError, ValueError)):
        run_cli("markers", metacool, "-a", other, "-o", out, "--inplace")
    assert metacool.read_bytes() == before
    message = {
        "a contig missing": f"{metacool} is not a contact map of the contigs of the assembly: 1 of its contig(s) are not "
        f"in the assembly {other}: '{name}'. Give the assembly that the map was made from.",
        "another length": f"{metacool} is not a contact map of the contigs of the assembly {other}: 1 contig(s) have "
        f"another length in it: '{name}' ({len(records[0]) + 1} bp in the map, {len(records[0])} bp in the assembly).",
        "two records of a contig": f"1 contig(s) are named by several records of the assembly {other}: '{name}'. A contig "
        "is read by its name, the first word of its header: give each record a name of its own.",
        "not a FASTA file": f"{other} is not a FASTA file: give the assembly that the map was made from.",
    }[case]
    assert ("ERROR", message) in logged()


def test_markers_of_a_gzip_compressed_assembly(marker_run, tmp_path, run_cli):
    """A gzip-compressed assembly gives the markers of the plain one."""
    import gzip

    metacool, assembly, _ = marker_run
    compressed = tmp_path / "assembly.fa.gz"
    with open(assembly, "rb") as plain, gzip.open(compressed, "wb") as handle:
        shutil.copyfileobj(plain, handle)
    run_cli("markers", metacool, "-a", assembly, "-o", tmp_path / "plain")
    run_cli("markers", metacool, "-a", compressed, "-o", tmp_path / "gzip")
    pd.testing.assert_frame_equal(_stored(tmp_path / "plain" / metacool.name)[0], _stored(tmp_path / "gzip" / metacool.name)[0])


def test_the_runs_that_regroup_the_contigs_keep_the_markers(marker_run, tmp_path, run_cli):
    """The markers are facts about the contigs: a run that regroups them (the
    final clusters of a validation, the MGE genomes of mge --reassign) or stores
    the network keeps the columns and the entry."""
    metacool, assembly, out = marker_run
    run_cli("markers", metacool, "-a", assembly, "--inplace", "-o", out)
    table, entry = _stored(metacool)
    contigs = list(table.index)
    run_cli("network", metacool, "-a", assembly, "--inplace", "-o", out, "-n", "None")
    clusters = mtmc.FinalClusters({name: "metator_00001_00000" for name in contigs[:10]}, None, None, {})
    mtmc.write_map(str(metacool), str(metacool), clusters=clusters, force=True)
    mc = metacooler.MetaCooler(str(metacool))
    assert set(mc.genomes["name"]) == {"__unbinned__", "metator_00001_00000"}
    kept, kept_entry = _stored(metacool)
    pd.testing.assert_frame_equal(kept.reindex(contigs), table)
    assert kept_entry == entry
    assert {"network", "validation", "markers"} <= set(mc.metadata["metator"])


def test_write_map_stores_the_markers_alone(marker_run):
    """write_map stores the markers alone, as it stores the hosts."""
    metacool, _, _ = marker_run
    markers = mtmc.Markers(pd.DataFrame({"Name": []}), {})
    with pytest.raises(ValueError, match="stores the marker genes of the contigs of metator markers alone"):
        mtmc.write_map(str(metacool), str(metacool), markers=markers, rounds=mtcl.Rounds().tables(), force=True)
    with pytest.raises(ValueError, match="stores the hosts of the MGE clusters of metator host alone"):
        mtmc.write_map(str(metacool), str(metacool), markers=markers, hosts=({}, {}), force=True)


def test_holds_markers(marker_run, run_cli, tmp_path):
    """A map holds markers once metator markers stored them."""
    metacool, assembly, out = marker_run
    assert not mtmc.holds_markers(metacooler.MetaCooler(str(metacool)))
    run_cli("markers", metacool, "-a", assembly, "--inplace", "-o", out)
    assert mtmc.holds_markers(metacooler.MetaCooler(str(metacool)))


def test_the_metadata_stays_small(marker_run, run_cli):
    """The markers entry, the marker sets with it, adds about 10 KB to the
    metadata of the map, a JSON attribute: the names and weights of the
    markers, 8.7 KB, which miComplete's files give to 17 digits."""
    metacool, assembly, out = marker_run
    run_cli("markers", metacool, "-a", assembly, "--inplace", "-o", out)
    _, entry = _stored(metacool)
    assert len(json.dumps(entry)) < 12_000
    assert not any(math.isnan(weight) for spec in entry["sets"].values() for weight in spec["weights"])
    assert os.path.getsize(metacool) > 0
