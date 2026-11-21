# Tests of metator.metacool: the contact map of a binning as a .metacool
# file, on the mock metagenome of tests/synthetic.py. Its genome contigs are
# binned by genome and its MGE contigs left unbinned.

import logging
import os
import re
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
import metator.network as mtn
from metator.version import __version__
from tests import synthetic
from tests.hdf5 import walk
from tests.maps import contig_map, pairs_map
from tests.oracles import normalize_pair

PAIRS_HEADER = "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"


def reassign(
    metacool, contig_data, outfile, bin_summary=None, force=False, rounds=None, planned=None, metadata=None,
    normalization="None", **stores,
):  # fmt: skip
    """The map of the clusters of a run, from the input map metacool, written
    at outfile as metator validation writes it (see
    metator.metacool.write_map): the final clusters of a contig data file
    (its Name and Final_cluster columns, ND, - or empty for none), the
    summary of the clusters (bin_summary: none, the cluster names alone),
    the clustering rounds of the run (rounds), its entry validation (the
    version of metaTOR, and metadata), and what stores adds (the partition
    entry of a pipeline run), with the network of the contigs, as pipeline
    stores it: the columns of the contig data that metaTOR measured (those
    of CONTIG_COLUMNS with numbers), and the weights of the normalization
    (none for "None", the default); normalization None: no network."""
    table = pd.read_csv(contig_data, sep="\t", dtype={"Name": str, "Final_cluster": str}, keep_default_na=False)
    if normalization is not None:
        stores["network"] = mtmc.Network(table, normalization)
    final = table["Final_cluster"].str.strip()
    binned = ~final.isin(mtmc.UNBINNED)
    clusters = mtmc.FinalClusters(
        dict(zip(table["Name"][binned], final[binned], strict=True)),
        None if bin_summary is None else str(bin_summary),
        rounds,
        {"version": __version__, **(metadata or {})},
    )
    mtmc.write_map(str(metacool), str(outfile), clusters=clusters, force=force, planned=planned, **stores)


def _bin_names(metagenome):
    """Final bin of each contig: one bin per genome, none for the MGEs."""
    names = {genome.name: f"metator_{i:05d}_00000" for i, genome in enumerate(metagenome.genomes, start=1)}
    return {contig.name: names.get(contig.molecule, "ND") for contig in metagenome.contigs}


# The contig columns of the map: those of contig_data_final.txt, as
# write_contig_data writes it, but Restriction_site, which is - everywhere.
CONTIG_COLUMNS = ["GC_content", "Hit", "Shotgun_coverage"]


def write_contig_data(metagenome, path, final_bin=None):
    """A contig_data_final.txt file, laid out as metaTOR writes it (an
    unnamed index first, Final_cluster last), of a run with a depth file and no
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
            "Overlapping_cluster_ID": "00001",
            "Recursive_cluster_contigs": "-",
            "Final_cluster": truth.Name.map(final_bin).to_numpy(),
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
    """A cluster_summary.txt file, written as validation writes it."""
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


@pytest.fixture
def binning(metagenome, metagenome_files, tmp_path):
    """Pairs, contig data and bin summary of a binning of the mock community."""
    return {
        "pairs": metagenome_files["pairs"],
        "contig_data": write_contig_data(metagenome, tmp_path / "contig_data_final.txt"),
        "bin_summary": write_bin_summary(metagenome, tmp_path / "cluster_summary.txt"),
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


def test_unbinned_contigs_keep_their_contacts(metagenome, binning, tmp_path):
    """Contigs with Final_cluster ND, - or empty are in __unbinned__, with all
    their contacts; a map without bin summary."""
    final_bin = _bin_names(metagenome)
    names = [contig.name for contig in metagenome.contigs]
    for name, value in zip(names[:3], ["ND", "-", ""], strict=True):
        final_bin[name] = value
    contig_data = write_contig_data(metagenome, tmp_path / "contig_data_final.txt", final_bin)
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    out = tmp_path / "map.metacool"
    reassign(str(metacool), str(contig_data), str(out))

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


def test_non_ascii_cluster_names_are_refused(metagenome, binning, tmp_path, monkeypatch):
    """Before the map is read."""
    final_bin = _bin_names(metagenome)
    names = list(final_bin)
    final_bin[names[0]] = final_bin[names[1]] = "bin_é"
    contig_data = write_contig_data(metagenome, tmp_path / "contig_data_final.txt", final_bin)
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    monkeypatch.setattr(metacooler, "MetaCooler", lambda *args, **kwargs: pytest.fail("map read"))
    with pytest.raises(mtmc.MetacoolError, match=r"needs ASCII cluster names .* are not: 'bin_é'\."):
        reassign(str(metacool), str(contig_data), str(tmp_path / "map.metacool"))


def test_a_clustered_contig_that_the_map_lacks_is_refused(metagenome, binning, tmp_path):
    """A contig of the final clusters that the map lacks (here one whose name
    is not ASCII, which no map holds) is refused, named, and the map is not
    written: the clusters of a run are clusters of the contigs of its map, by
    their names in it, and a contig named otherwise (a name read as a number,
    00123 as 123) is another contig."""
    final_bin = _bin_names(metagenome)
    clustered = next(name for name, cluster in final_bin.items() if cluster != "ND")
    contig_data = write_contig_data(metagenome, tmp_path / "contig_data_final.txt", final_bin)
    text = contig_data.read_text()
    contig_data.write_text(text.replace(f"\t{clustered}\t", "\tcontig_é\t"))
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    out = tmp_path / "map.metacool"
    message = (
        f"1 contig(s) of the final clusters are not contigs of {metacool}: 'contig_é'. The clusters of a run are "
        "clusters of the contigs of its map, by their names in it."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        reassign(str(metacool), str(contig_data), str(out))
    assert not out.exists()


def test_non_ascii_names_shown():
    """The first names, then how many more: the message ends with one period."""
    names = r"'é0', 'é1', 'é2', 'é3', 'é4', \.\.\. \(\+2 more\)"
    with pytest.raises(mtmc.MetacoolError, match=rf"7 contig name\(s\) are not: {names}\.$"):
        mtmc.check_ascii_names(["ok", *[f"é{i}" for i in range(7)], "é0"], [])
    mtmc.check_ascii_names(["contig_1"], ["metator_00001_00000"])


NO_BIN_NAMES = ["NA", "nan", "N/A", "nd", "N/d", "-1", "__unbinned__", "__UNBINNED__"]


def test_bin_names_read_as_no_bin_are_refused(metagenome, binning, tmp_path, monkeypatch):
    """Bins that metacooler would silently empty into __unbinned__: refused,
    named, before the map is read. ND is metaTOR's own mark of a contig in
    no bin, not a bin name."""
    final_bin = _bin_names(metagenome)
    names = list(final_bin)
    for name, bin_name in zip(names, [*NO_BIN_NAMES, "NA", "ND"], strict=False):
        final_bin[name] = bin_name
    contig_data = write_contig_data(metagenome, tmp_path / "contig_data_final.txt", final_bin)
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    monkeypatch.setattr(metacooler, "MetaCooler", lambda *args, **kwargs: pytest.fail("map read"))
    shown = r"'NA', 'nan', 'N/A', 'nd', 'N/d', \.\.\. \(\+3 more\)"
    message = (
        "A .metacool file cannot have clusters named as metacooler marks a contig in no cluster (NA, NaN, N/A, N/D, "
        "-1 or __unbinned__ in any case, and nd in any case other than ND, which marks a contig in no cluster): it "
        "would put their contigs in __unbinned__. Rename the 8 cluster(s): "
    )
    with pytest.raises(mtmc.MetacoolError, match=rf"^{re.escape(message)}{shown}\.$"):
        reassign(str(metacool), str(contig_data), str(tmp_path / "map.metacool"))


def test_no_bin_names_are_those_of_metacooler():
    """The bin names refused are exactly the genome names that metacooler
    reads as no genome (-1 among them from metacooler 0.2.0)."""
    from metacooler.io.create import _UNBINNED_SENTINELS

    assert set(_UNBINNED_SENTINELS) == mtmc.METACOOLER_UNBINNED
    mtmc.check_bin_names(["metator_00001_00000", "MetaTOR_MGE_00001", "ND_1", "na1"])
    for name in mtmc.METACOOLER_UNBINNED:
        with pytest.raises(mtmc.MetacoolError, match="Rename the 2 cluster"):
            mtmc.check_bin_names(["metator_00001_00000", name.upper(), f" {name} "])


def test_misaligned_bin_summary(metagenome, binning, tmp_path):
    """A bin summary whose header has more fields than its rows (the names
    of miComplete's columns split at their spaces, as an old
    cluster_summary.txt had them: 18 and 15): metacooler cannot tell which
    column each value belongs to. Its reason is raised as a MetacoolError,
    and no map is left."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    columns = ["Length", "GC-content", "Present Markers", "Completeness", "Redundancy", "Weighted completeness",
               "Weighted redundancy", "Contigs", "N50", "L50", "N90", "L90", "CDs", "HiC_abundance"]  # fmt: skip
    clusters = sorted(set(_bin_names(metagenome).values()) - {"ND"})
    header = "\t" + "\t".join(" ".join(columns).split(" ")) + "\n"
    summary = tmp_path / "cluster_summary.txt"
    summary.write_text(header + "".join(name + "\t1" * len(columns) + "\n" for name in clusters))
    out = tmp_path / "map.metacool"
    with pytest.raises(mtmc.MetacoolError, match=r"^Genomes table is misaligned: .*cluster_summary\.txt\n") as error:
        reassign(str(metacool), str(binning["contig_data"]), str(out), bin_summary=str(summary))
    assert "header names 18 column(s), but no data row has more than 15" in str(error.value)
    assert isinstance(error.value.__cause__, ValueError)
    assert not out.exists()


def test_outfile_in_a_file(metagenome, binning, tmp_path):
    """An outfile whose directory cannot be made: a MetacoolError that says
    so, with the reason of the OSError."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    (tmp_path / "file").write_text("not a directory")
    outfile = tmp_path / "file" / "map.metacool"
    message = f"The directory of {outfile} cannot be made (File exists): write the map elsewhere."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$") as error:
        reassign(str(metacool), str(binning["contig_data"]), str(outfile))
    assert isinstance(error.value.__cause__, FileExistsError)


def test_metacooler_messages_go_to_the_log_file(metagenome, binning, tmp_path):
    """metacooler's messages reach the log file of the run; its logger is
    left as it was."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    log = tmp_path / "run.log"
    mc_logger = logging.getLogger("metacooler")
    before = (mc_logger.level, mc_logger.propagate, list(mc_logger.handlers))
    mtl.set_file_handler(str(log))
    try:
        reassign(str(metacool), str(binning["contig_data"]), str(tmp_path / "map.metacool"))
    finally:
        for handler in mtl.logger.handlers[:]:
            if isinstance(handler, logging.FileHandler):
                mtl.logger.removeHandler(handler)
                handler.close()
    logged = log.read_text()
    n_binned = sum(len(genome.contigs) for genome in metagenome.genomes)
    n_bins = len(metagenome.genomes)
    assert f": {len(metagenome.contigs)} contigs, {n_binned} of them in {n_bins} clusters." in logged
    assert "INFO :: Reassigned " in logged
    assert (mc_logger.level, mc_logger.propagate, list(mc_logger.handlers)) == before


def test_is_map(metagenome, binning, tmp_path):
    """A map is an HDF5 file, or the URI of a map inside one (file::group),
    whatever the name: a pairs file may have :: in its name, and a URI is
    no map when its file is no HDF5 file (or no file)."""
    metacool = unbinned_map(tmp_path / "map.metacool", metagenome, binning["pairs"])
    colons = shutil.copyfile(binning["pairs"], tmp_path / "lib::L1.pairs")
    named = shutil.copyfile(metacool, tmp_path / "lib::L1.metacool")
    assert [mtmc.is_map(path) for path in [metacool, named, f"{metacool}::/", f"{metacool}::entries/a"]] == [True] * 4
    pairs = binning["pairs"]
    others = [pairs, colons, f"{pairs}::/", tmp_path / "missing.metacool", f"{tmp_path}/missing::x", tmp_path]
    assert [mtmc.is_map(path) for path in others] == [False] * 6


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
    # The columns of the contigs table in any order: the oracle stores its
    # columns after the QC counts, which the map of a run carries before
    # them (see assert_contig_columns).
    chroms, ref_chroms = clr.chroms()[:], ref.chroms()[:]
    pd.testing.assert_frame_equal(chroms.sort_index(axis=1), ref_chroms.sort_index(axis=1))
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
    contig, in the order of the assembly): the map of the pairs (see
    pairs_map), QC included (metacooler 0.3.0 carries it over to
    the bins), but for its provenance. The input map is only read, and
    nothing but the map is left next to it."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    before = metacool.read_bytes()
    expected = tmp_path / "pairs.metacool"
    tables = {"bin_summary": str(binning["bin_summary"]), "metadata": {"command": "metator pipeline"}}
    pairs_map(binning["pairs"], binning["contig_data"], expected, **tables)
    out = tmp_path / "out" / "bins.metacool"
    reassign(str(metacool), str(binning["contig_data"]), str(out), normalization="None", **tables)

    metacooler.validate(str(out), full=True)
    assert_same_map(out, expected)
    assert_contig_columns(out, binning["contig_data"])
    assert metacool.read_bytes() == before
    clr = cooler.Cooler(str(out))
    metadata = clr.info["metadata"]
    assert metadata["metator"]["validation"] == {"version": __version__, "command": "metator pipeline"}
    assert metacooler.MetaCooler(str(out)).info["ordering"] == "genome-then-length"
    # The QC block and columns that metacooler measures on the pairs with
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
    reassign(str(metacool), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))

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
    reassign(str(metacool), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))

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
    """The map of the pairs of the binning (see pairs_map)."""
    path = tmp_path / "pairs.metacool"
    if not path.exists():
        pairs_map(binning["pairs"], binning["contig_data"], path, bin_summary=binning["bin_summary"])
    return path


@pytest.mark.parametrize("binned", ["map", "assembly"])
def test_map_of_the_bins_from_a_cut_map(binned, metagenome, binning, tmp_path, logged):
    """From a map cut by metacooler extract (a .cool file), then made a map
    again by create cool: the contigs that it keeps, in their final bins.
    A binning made on more contigs than the map's, with binned contigs that
    the map lacks, is refused: the clusters of a run are clusters of the
    contigs of its map."""
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
    lacking = sorted(name for name, bin_name in final_bin.items() if bin_name != "ND" and name not in kept)
    assert bool(lacking) == (binned == "assembly")
    if lacking:
        with pytest.raises(mtmc.MetacoolError, match=rf"^{len(lacking)} contig\(s\) of the final clusters are not contigs"):
            reassign(str(metacool), str(contig_data), str(out), bin_summary=str(binning["bin_summary"]))
        assert not out.exists()
        return
    reassign(str(metacool), str(contig_data), str(out), bin_summary=str(binning["bin_summary"]), normalization="None")

    metacooler.validate(str(out), full=True)
    assert assignment(out) == {name: found for name, found in final_bins(metagenome).items() if name in kept}
    assert contig_contacts(out) == contig_contacts(metacool)
    assert_contig_columns(out, contig_data)
    assert cooler.Cooler(str(out)).info["metadata"]["qc"]["source"] == "pixels"
    warnings = [message for level, message in logged() if level == "WARNING"]
    # The contig data lists contigs that the map lacks: annotate reads the
    # table of the network, ignoring their rows without a warning
    # (unmatched="ignore"), and the history of the map names that table.
    assert not [message for message in warnings if "Contigs table" in message]
    assert cooler.Cooler(str(out)).info["metadata"]["annotate"][-1]["contigs_table"].endswith("contig_data_network.txt")
    assert not [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]


def test_map_of_the_bins_from_an_mmetacool_entry(metagenome, binning, tmp_path):
    """From an entry of a .mmetacool collection: the map of the bins of the
    map it holds, a .metacool file of its own."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    uri = f"{collection}::entries/lib1"
    for source, out in [(str(metacool), tmp_path / "from_map.metacool"), (uri, tmp_path / "from_entry.metacool")]:
        reassign(source, str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))
        metacooler.validate(str(out), full=True)
    assert_same_map(tmp_path / "from_entry.metacool", tmp_path / "from_map.metacool")
    metadata = cooler.Cooler(str(tmp_path / "from_entry.metacool")).info["metadata"]
    assert [entry["source"] for entry in metadata["reassign"]] == [uri]


def test_map_of_the_bins_from_a_map_without_qc(metagenome, binning, tmp_path):
    """From a map that has no QC block: the map of the bins has none either."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"], qc=False)
    out = tmp_path / "bins.metacool"
    reassign(str(metacool), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]))
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
    reassign(str(metacool), str(contig_data), str(out), bin_summary=str(binning["bin_summary"]), normalization="None")

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
    reassign(str(earlier), str(contig_data), str(out), normalization="None")

    metacooler.validate(str(out), full=True)
    assert assignment(out) == {name: (bin_name if bin_name != "ND" else "__unbinned__") for name, bin_name in final_bin.items()}
    assert list(genomes_without_qc(out).columns) == ["name", "length", "n_contigs", "n_bins"]
    assert_contig_columns(out, contig_data)
    assert [message for level, message in logged() if level == "WARNING"] == []
    assert set(tmp_path.iterdir()) - before == {out}


@pytest.mark.parametrize("earlier", ["2.0", "1.x"])
def test_map_of_the_bins_drops_the_columns_of_an_earlier_run(earlier, metagenome, binning, tmp_path):
    """From a map made with the contig data of an earlier run as its
    contigs table (metacooler create pairs --contigs-format Name,Final_cluster,
    which stores all its columns): the clusters and contig ids of that run,
    and its other columns, are not carried into the map of the clusters;
    only the columns of the network of this run, and a column that metaTOR
    does not write. The same from the contig data of a run of metaTOR 1.x,
    its cluster columns named after bins (Final_bin, Overlapping_bin_ID...):
    they would contradict the clusters of this run (see
    metator.metacool.METATOR_TABLE_COLUMNS)."""
    contig_data = str(binning["contig_data"])
    earlier_table = contig_data if earlier == "2.0" else synthetic.as_1x(contig_data, tmp_path / "contig_data_1x.txt")
    final, bins = {"2.0": ("Final_cluster", "cluster"), "1.x": ("Final_bin", "bin")}[earlier]
    table = pd.read_csv(earlier_table, sep="\t", index_col=0)
    earlier_map = tmp_path / "earlier.metacool"
    metacooler.from_pairs(
        str(earlier_map),
        bins=pd.DataFrame({"chrom": table.Name, "start": 0, "end": table.Size}),
        pairs=str(binning["pairs"]),
        contig_genome_map=earlier_table,
        contig_table_format=f"Name,{final}",
        ps=False,
    )
    taxa = tmp_path / "taxa.tsv"
    taxa.write_text("Name\ttaxon\n" + "".join(f"{name}\tt_{i % 3}\n" for i, name in enumerate(table.Name)))
    metacooler.annotate(str(earlier_map), contigs_table=str(taxa), contigs_table_format="Name", inplace=True)
    stale = {"Unnamed: 0", "ID", "Size", "Restriction_site", f"Overlapping_{bins}_ID", f"Recursive_{bins}_contigs"}
    assert stale <= set(metacooler.MetaCooler(str(earlier_map)).contigs.columns)
    out = tmp_path / "clusters.metacool"
    reassign(str(earlier_map), contig_data, str(out), bin_summary=str(binning["bin_summary"]))

    metacooler.validate(str(out), full=True)
    assert assignment(out) == final_bins(metagenome)
    contigs = metacooler.MetaCooler(str(out)).contigs
    assert sorted(without_qc(contigs.columns)) == sorted(["name", "length", "genome_id", *CONTIG_COLUMNS, "taxon"])
    # The network last, after the QC columns.
    assert list(contigs.columns)[-len(CONTIG_COLUMNS) :] == CONTIG_COLUMNS


@pytest.mark.parametrize("earlier_bins", ["same", "other"])
def test_map_of_the_bins_from_the_map_of_an_earlier_run(earlier_bins, metagenome, binning, tmp_path):
    """From the map of an earlier run, which carries metaTOR's contig
    columns: the map of the pairs (see pairs_map), contigs table
    included, metaTOR's columns last, as there. metacooler reassign stores
    the columns that it carries before the QC columns when it groups the
    contigs again (other bins): the contigs table is then that of the map
    of the pairs, column for column; and so it is when it relabels them in
    place (the same bins), which keeps the order of the columns that
    create writes, genome_id third (metacooler 0.4.0)."""
    earlier = tmp_path / "earlier.metacool"
    if earlier_bins == "same":
        shutil.copyfile(mtmc_map(binning, tmp_path), earlier)
    else:
        # The first two genomes in one bin.
        final_bin = _bin_names(metagenome)
        merged = {final_bin[genome.contigs[0].name] for genome in metagenome.genomes[:2]}
        final_bin = {name: ("metator_00099_00000" if bin_name in merged else bin_name) for name, bin_name in final_bin.items()}
        contig_data = write_contig_data(metagenome, tmp_path / "earlier.txt", final_bin)
        pairs_map(binning["pairs"], contig_data, earlier)
    out = tmp_path / "bins.metacool"
    reassign(str(earlier), str(binning["contig_data"]), str(out), bin_summary=str(binning["bin_summary"]), normalization="None")

    metacooler.validate(str(out), full=True)
    expected = mtmc_map(binning, tmp_path)
    clr, ref = cooler.Cooler(str(out)), cooler.Cooler(str(expected))
    for table in ["bins", "pixels"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    chroms, ref_chroms = clr.chroms()[:], ref.chroms()[:]
    assert list(ref_chroms.columns) == ["name", "length", "genome_id", *QC_CONTIG_COLUMNS, *CONTIG_COLUMNS]
    pd.testing.assert_frame_equal(chroms.sort_index(axis=1), ref_chroms.sort_index(axis=1))
    assert_contig_columns(out, binning["contig_data"])
    # The genomes, QC block and columns included, carried over to the bins.
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(out)).genomes, metacooler.MetaCooler(str(expected)).genomes)
    assert clr.info["metadata"]["qc"] == ref.info["metadata"]["qc"]


@pytest.mark.parametrize("how", ["path", "root URI", "link to the outfile", "outfile a link"])
def test_map_of_the_bins_replaces_the_input_map(how, metagenome, binning, tmp_path, logged):
    """The outfile may be the input map itself (a rebinning into the same
    directory): the same path (also as the URI of its root, file::/), or the
    file that the input map, given as a symbolic link, points to. The map
    of the bins is written in an edit session of metacooler, validated, and
    only then replaces the input map: the file that a symbolic link points
    to, as metacooler edits the target of a link, which keeps pointing at
    it. An outfile that is a symbolic link to the input map is a map of its
    own: the link is replaced (force), and the input map left as it was, as
    before the map went through metacooler. No other file is left, and the
    log says it is written."""
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
    in_place = how != "outfile a link"
    assert mtmc.is_input_map(str(metacool), str(outfile)) == in_place
    reassign(str(metacool), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True)

    assert not outfile.is_symlink()
    metacooler.validate(str(outfile), full=True)
    assert_same_map(outfile, mtmc_map(binning, tmp_path))
    assert sorted(path.name for path in out.iterdir()) == ["metator.metacool"]
    assert ("INFO", f"{outfile} is written: the contact map of the clusters.") in logged()
    # In place, the input map is the map edited (through a link, its
    # target); else it is left as it was.
    assert (Path(mtmc.map_file(metacool)).read_bytes() == before) != in_place
    if how == "link to the outfile":
        assert Path(metacool).is_symlink() and Path(metacool).resolve() == outfile
    assert not [path.name for path in tmp_path.iterdir() if path.name.startswith(".") or path.suffix == ".lock"]


@pytest.mark.parametrize("writer", ["reassign_metacool", "write_rounds"])
def test_a_hard_link_of_the_input_map_is_refused(writer, metagenome, binning, tmp_path):
    """An outfile that is the input map under another name (a hard link) is
    refused before any work, with metaTOR's reason: written as a map of its
    own, it would replace that name only, and metacooler refuses an output
    that is the same file as an input. It updated the input map in place,
    and the hard link kept the old map."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    outfile = tmp_path / "bins.metacool"
    os.link(metacool, outfile)
    before = metacool.read_bytes()
    message = (
        f"{outfile} is a hard link of the input map {metacool}: the same file under another name. Write the map to "
        f"another file, or to {os.path.realpath(metacool)} to update the map in place."
    )
    write = {
        "reassign_metacool": lambda: reassign(str(metacool), str(binning["contig_data"]), str(outfile), force=True),
        "write_rounds": lambda: write_rounds(str(metacool), str(outfile), _one_round(metacool), force=True),
    }[writer]
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        write()
    assert metacool.read_bytes() == before and os.path.samefile(metacool, outfile)
    assert not [path.name for path in tmp_path.iterdir() if path.name.startswith(".") or path.suffix == ".lock"]


def test_map_of_the_bins_whose_only_edit_is_the_reassignment(metagenome, binning, tmp_path):
    """A map of the bins whose only edit is the reassignment (a contig data
    without a column that metaTOR measured, no rounds, no normalization):
    metacooler reassign writes the staged map by a rename of its own, which
    the edit session must count as a write, and the map is written. A guard
    for the re-pin: metacooler's revised edit session installs a map only
    when it was written."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    table = pd.read_csv(binning["contig_data"], sep="\t", dtype=str, keep_default_na=False)[["Name", "Size", "Final_cluster"]]
    contig_data = tmp_path / "contig_data_final.txt"
    table.to_csv(contig_data, sep="\t", index=False)
    reassign(str(metacool), str(contig_data), str(metacool), force=True)
    genomes = set(metacooler.MetaCooler(str(metacool)).genomes.name.astype(str))
    assert genomes - {"__unbinned__"} == set(table.Final_cluster) - mtmc.UNBINNED


def test_a_map_written_as_metacooler_copies_it_keeps_metacoolers_reason(metagenome, binning, tmp_path, monkeypatch):
    """A map that a program writes as metacooler copies it into the map of
    the bins (to another file): the map is not written, with metacooler's
    reason, kept, not a read of the map that would mix two maps."""
    import metacooler.io.staging as staging

    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    copy_fresh = staging.copy_fresh

    def copied_then_written(source, target, **kwargs):
        copy_fresh(source, target, **kwargs)
        state = os.stat(metacool)
        os.utime(metacool, ns=(state.st_atime_ns, state.st_mtime_ns + 1_000_000))

    monkeypatch.setattr(staging, "copy_fresh", copied_then_written)
    outfile = tmp_path / "bins.metacool"
    with pytest.raises(mtmc.MetacoolError) as raised:
        reassign(str(metacool), str(binning["contig_data"]), str(outfile))
    assert isinstance(raised.value.__cause__, metacooler.MapChangedError)
    assert str(raised.value) == str(raised.value.__cause__)
    assert not outfile.exists()


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
        # The full validation of the staged map, before metacooler renames
        # it over the outfile.
        monkeypatch.setattr("metacooler.schema.validate", _fail_in_validate)
        message = "pixels out of order"
    else:
        # The validation of metacooler reassign, which writes the map to its
        # partial file first.
        monkeypatch.setattr("metacooler.io.reassign.validate", _fail_in_validate)
        message = "pixels out of order"
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$") as raised:
        reassign(str(metacool), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True)
    assert str(raised.value.__cause__) == message
    assert outfile.read_bytes() == before
    assert sorted(path.name for path in out.iterdir()) == ["metator.metacool"]
    messages = logged()
    assert not [message for level, message in messages if level == "ERROR" and "for inspection" in message]
    if failure == "in metacooler":
        # The temporary file of the map, next to the outfile.
        staged = removed_partial(messages)
        # Staged by metacooler in a private directory next to the map.
        assert Path(staged).parent.parent == out and Path(staged).parent.name.startswith(".metator.metacool.")


def _replaced_while_read(monkeypatch, path, by):
    """Have the map at path replaced by the map by, as a run that updates a
    map in place replaces it (os.replace): once the first chunk of its
    pixels is read (MetaCooler.iter_pixels), or once metacooler reassign
    has read it, or once metacooler.network has summed the first chunk of
    its pixels."""

    def replace():
        copy = f"{path}.copy"
        shutil.copyfile(by, copy)
        os.replace(copy, path)

    iter_pixels = metacooler.MetaCooler.iter_pixels

    def read_then_replace(self, *args, **kwargs):
        chunks = iter_pixels(self, *args, **kwargs)
        yield next(chunks)
        replace()
        yield from chunks

    def then_replace(read):
        def replaced(*args, **kwargs):
            written = read(*args, **kwargs)
            replace()
            return written

        return replaced

    # metacooler.network, which metator.network.contig_contacts reads the map
    # with: once the first chunk of its pixels is summed.
    import metacooler.io.network as network_module

    aggregate_pixels, replaced = network_module.aggregate_pixels, []

    def aggregate_then_replace(*args, **kwargs):
        summed = aggregate_pixels(*args, **kwargs)
        if not replaced:
            replaced.append(True)
            replace()
        return summed

    monkeypatch.setattr(metacooler.MetaCooler, "iter_pixels", read_then_replace)
    monkeypatch.setattr(network_module, "aggregate_pixels", aggregate_then_replace)
    monkeypatch.setattr(metacooler, "reassign", then_replace(metacooler.reassign))


def _map_and_its_clusters(metagenome, binning, tmp_path):
    """A map of the contigs, in the order of the assembly, and the map of
    its clusters, whose bins are grouped by cluster: another order."""
    unbinned = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    clusters = tmp_path / "clusters.metacool"
    reassign(str(unbinned), str(binning["contig_data"]), str(clusters), bin_summary=str(binning["bin_summary"]))
    return unbinned, clusters


@pytest.mark.parametrize("reader", ["contig_contacts"])
def test_a_map_replaced_while_read_is_refused(reader, metagenome, metagenome_files, binning, tmp_path, monkeypatch):
    """A map replaced while it is read, as a run that updates it in place
    replaces it (here by the map of its clusters), is refused: the map that
    a run writes to --output-map is not written, the file there left as it
    was. metacooler opens the file of a map for each chunk of pixels that it
    reads: these reads (network MAP, a run from a map that writes
    --output-map) mixed the bins of one map with the pixels of the other,
    silently (metacooler's tests show it)."""
    unbinned, clusters = _map_and_its_clusters(metagenome, binning, tmp_path)
    output = tmp_path / "output"
    output.write_text("an earlier output")
    reads = {"contig_contacts": lambda: mtn.contig_contacts(str(unbinned), chunksize=100)}
    _replaced_while_read(monkeypatch, unbinned, clusters)
    message = (
        f"{unbinned} was replaced or written while it was read, by a run that updated it: what was read of it would "
        "mix two maps. Read it again once that run has ended."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        reads[reader]()
    assert unbinned.read_bytes() == clusters.read_bytes()
    assert output.read_text() == "an earlier output"
    assert sorted(path.name for path in tmp_path.iterdir() if path.name.startswith(".")) == []


def test_a_read_that_fails_on_a_replaced_map_says_why(metagenome, metagenome_files, binning, tmp_path, monkeypatch):
    """A read that fails because the map was replaced meanwhile (here by a
    file that is no map, which cannot be read on) says that the map was
    replaced, from the error of the read, which alone named no cause."""
    unbinned = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    not_a_map = tmp_path / "not_a_map.txt"
    not_a_map.write_text("no map")
    _replaced_while_read(monkeypatch, unbinned, not_a_map)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(str(unbinned))} was replaced or written") as raised:
        mtn.contig_contacts(str(unbinned), chunksize=100)
    assert isinstance(raised.value.__cause__, OSError)


def _stopped_as_cooler_writes(uri, table, output, **kwargs):
    """metacooler.reassign stopped by a Ctrl-C (or the SIGTERM of a run) as
    cooler writes the map: cooler builds its writer out of the with that
    would remove its file on a failure, and metacooler leaves its partial
    file too."""
    Path(f"{output}.partial").write_text("metacooler's partial file")
    Path(f"{output}.partial.__writing__").write_text("cooler's file")
    raise KeyboardInterrupt


def test_a_stop_as_cooler_writes_leaves_no_file(metagenome, binning, tmp_path, monkeypatch):
    """A Ctrl-C or a SIGTERM that stops cooler as it writes the map of the
    clusters leaves the input map (the outfile) as it was, and no file of
    the run next to it: cooler's (<staged>.partial.__writing__) stayed
    there. The files of another run (another temporary file) are left."""
    out = tmp_path / "out"
    out.mkdir()
    outfile = unbinned_map(out / "metator.metacool", metagenome, binning["pairs"])
    before = outfile.read_bytes()
    others = [out / ".metator.metacool.another.tmp", out / ".metator.metacool.another.tmp.partial.__writing__"]
    for path in others:
        path.write_text("another run's")
    monkeypatch.setattr(metacooler, "reassign", _stopped_as_cooler_writes)
    with pytest.raises(KeyboardInterrupt):
        reassign(str(outfile), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True)
    assert outfile.read_bytes() == before
    assert sorted(path.name for path in out.iterdir()) == sorted(["metator.metacool", *(path.name for path in others)])


def _fail_unexpectedly(path, full=True):
    raise RuntimeError("unexpected failure")


def test_partial_file_removed_on_an_unexpected_error(metagenome, binning, tmp_path, monkeypatch, logged):
    """An error of another type than those metacooler is expected to raise,
    once the staged map is written (in its validation, here), is raised as
    it is: the staged copy is removed, and the input map (the outfile) is
    left as it was, with no other file next to it."""
    out = tmp_path / "out"
    out.mkdir()
    outfile = unbinned_map(out / "metator.metacool", metagenome, binning["pairs"])
    before = outfile.read_bytes()
    # The full validation of the staged map, before metacooler renames it.
    monkeypatch.setattr("metacooler.schema.validate", _fail_unexpectedly)
    with pytest.raises(RuntimeError, match=r"^unexpected failure$"):
        reassign(str(outfile), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True)
    assert outfile.read_bytes() == before
    # The staged copy, in a private directory next to the map, is removed.
    assert sorted(path.name for path in out.iterdir()) == ["metator.metacool"]
    assert [message for level, message in logged() if level in ("WARNING", "ERROR")] == []


def test_map_of_the_bins_existing_outfile(metagenome, binning, tmp_path):
    """An existing outfile, the input map or another file: refused without
    force, and left as it was."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    other = tmp_path / "other.metacool"
    other.write_text("old")
    for outfile in [metacool, other]:
        before = outfile.read_bytes()
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(f'{outfile} already exists. Use --force to overwrite it.')}$"):
            reassign(str(metacool), str(binning["contig_data"]), str(outfile))
        assert outfile.read_bytes() == before
    reassign(str(metacool), str(binning["contig_data"]), str(other), force=True)
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

    def refused(*args, **kwargs):
        raise error

    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    monkeypatch.setattr(metacooler, "reassign", refused)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(error.args[0])}$") as raised:
        reassign(str(metacool), str(binning["contig_data"]), str(tmp_path / "bins.metacool"))
    assert raised.value.__cause__ is error
    assert sorted(path.name for path in tmp_path.iterdir() if "metacool" in path.name) == ["contigs.metacool"]


REASSIGN_ERRORS = {
    "missing map": "Input file(s) not found: missing.metacool.",
    "not a map": "Cannot read the contact map ",
    "missing entry": "Cannot read the contact map maps.mmetacool::entries/lib2: ",
    "collection as outfile": (
        "maps.mmetacool is the .mmetacool collection of the input map maps.mmetacool::entries/lib1: the map would "
        "replace the whole collection. Write it to another file."
    ),
    "bin named NA": "Rename the 1 cluster(s): 'NA'.",
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
    elif case == "bin named NA":
        contig_data = write_contig_data(
            metagenome, tmp_path / "na.txt", {**_bin_names(metagenome), metagenome.contigs[0].name: "NA"}
        )
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()}
    with pytest.raises(mtmc.MetacoolError, match=re.escape(REASSIGN_ERRORS[case])):
        reassign(uri, str(contig_data), outfile, force=True)
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()} == before


# The map of the contigs, in no cluster yet, that a run starts from (see
# tests.maps.contig_map), and the map of its clusters written from it.


def assembly_contig_data(inputs, measured=True):
    """The contig data of the assembly of the mock community, as the
    pipeline makes it: with its depth file and an enzyme (the shotgun
    coverages and restriction sites measured), or neither."""
    depth, enzyme = (str(inputs["depth"]), "DpnII") if measured else (None, None)
    return mtn.create_contig_data(str(inputs["assembly"]), depth, enzyme)


def final_contig_data(contig_data, pairs, final_bin, path):
    """The contig_data_final.txt file of a run with this contig data, laid
    out as validation writes it (an unnamed index first, Final_cluster last):
    its columns as contig_data_network.txt held them, the hits of each
    contig counted on the pairs, as the network step counts them."""
    hits = pd.concat([pairs.chr1, pairs.chr2]).value_counts()
    with open(path, "w") as table:
        table.write("ID\tName\tSize\tGC_content\tHit\tShotgun_coverage\tRestriction_site\n")
        for i, (name, data) in enumerate(contig_data.items(), start=1):
            data["hit"] = int(hits.get(name, 0))
            values = [i, name, data["length"], data["GC"], data["hit"], data["coverage"], data["RS"]]
            table.write("\t".join(map(str, values)) + "\n")
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    table["Final_cluster"] = table.Name.map(final_bin)
    table.to_csv(path, sep="\t", index=True)
    return path


@pytest.mark.parametrize("final", ["bins", "no bin", "one bin", "longest first"])
def test_map_of_the_bins_from_the_map_of_the_contigs(final, metagenome, metagenome_inputs, tmp_path):
    """The map of the bins that reassign_metacool writes over the map of the
    contigs, as a pipeline run does at its end, is the map of the same
    pairs files and final tables (see pairs_map): the
    same bins, contacts, contigs table (metaTOR's columns last, as there),
    genomes table (QC columns included) and QC block, but for its history
    (the annotate and reassign keys of its metadata) and its date. When the
    final bins leave the contigs in the order of the map of the contigs (for
    example, no final bin, one bin holding every contig, or the longest
    contigs in the first bin and the others in the second), metacooler
    relabels them where they are, in the order of the columns that create
    writes, genome_id third (metacooler 0.4.0)."""
    pairs = str(metagenome_inputs["pairs"])
    contig_data = assembly_contig_data(metagenome_inputs)
    out = tmp_path / "out"
    outfile = out / "metator.metacool"
    run = {"metadata": {"command": "metator pipeline"}}
    out.mkdir()
    contig_map(pairs, contig_data, outfile, assembly_name="synthetic")
    final_bin = _bin_names(metagenome)
    if final in ["no bin", "one bin"]:
        final_bin = dict.fromkeys(final_bin, "ND" if final == "no bin" else "metator_00001_00000")
    elif final == "longest first":
        # The 20 first contigs of the map of the contigs, the longest, in
        # the first bin, and the others in the second.
        order = metacooler.MetaCooler(str(outfile)).contigs.name
        final_bin = {name: f"metator_0000{1 + (i >= 20)}_00000" for i, name in enumerate(order)}
    contig_data_file = final_contig_data(contig_data, metagenome.pairs, final_bin, out / "contig_data_final.txt")
    tables = {"bin_summary": str(write_bin_summary(metagenome, out / "cluster_summary.txt")), **run}
    reassign(str(outfile), str(contig_data_file), str(outfile), force=True, **tables)
    expected = tmp_path / "pairs.metacool"
    pairs_map(pairs, contig_data_file, expected, assembly_name="synthetic", **tables)

    metacooler.validate(str(outfile), full=True)
    clr, ref = cooler.Cooler(str(outfile)), cooler.Cooler(str(expected))
    for table in ["bins", "pixels"]:
        pd.testing.assert_frame_equal(getattr(clr, table)()[:], getattr(ref, table)()[:])
    chroms, ref_chroms = clr.chroms()[:], ref.chroms()[:]
    columns = ["GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]
    assert list(ref_chroms.columns) == ["name", "length", "genome_id", *QC_CONTIG_COLUMNS, *columns]
    pd.testing.assert_frame_equal(chroms.sort_index(axis=1), ref_chroms.sort_index(axis=1))
    pd.testing.assert_frame_equal(metacooler.MetaCooler(str(outfile)).genomes, metacooler.MetaCooler(str(expected)).genomes)
    info, ref_info = dict(clr.info), dict(ref.info)
    metadata, ref_metadata = info.pop("metadata"), ref_info.pop("metadata")
    assert metadata["qc"] == ref_metadata["qc"] and metadata["qc"]["source"] == "pairs"
    assert {key for key in info.keys() | ref_info.keys() if info.get(key) != ref_info.get(key)} <= {"creation-date"}
    assert {key for key in metadata.keys() | ref_metadata.keys() if metadata.get(key) != ref_metadata.get(key)} == {
        "annotate",
        "metator",
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
    assert sorted(path.name for path in out.iterdir()) == ["cluster_summary.txt", "contig_data_final.txt", "metator.metacool"]


@pytest.mark.parametrize("mode", [0o600, 0o640])
def test_a_map_updated_in_place_keeps_its_mode(metagenome, binning, tmp_path, mode):
    """The map of the clusters takes the mode of the file it replaces: a
    private map stays private (it got the mode of a new file)."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    metacool.chmod(mode)
    reassign(str(metacool), str(binning["contig_data"]), str(metacool), bin_summary=str(binning["bin_summary"]), force=True)
    assert metacool.stat().st_mode & 0o777 == mode


def test_a_map_updated_again_and_again_does_not_grow(metagenome, binning, tmp_path):
    """A map updated in place run after run stays bounded: each edit session
    starts from a fresh copy of the map, so the map holds at most the space
    that the edits of one run left behind. From the second run on, it keeps
    its size, but for the history of the runs that its metadata gains (it
    grew by what each run's edits replaced)."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    sizes = []
    for _ in range(6):
        reassign(str(metacool), str(binning["contig_data"]), str(metacool), bin_summary=str(binning["bin_summary"]), force=True)
        sizes.append(metacool.stat().st_size)
    assert max(sizes[1:]) - sizes[1] < 8_000, sizes
    metacooler.validate(str(metacool), full=True)


def test_a_network_stored_again_and_again_stays_bounded(metagenome, metagenome_inputs, tmp_path):
    """The network stored in place run after run (metator network) stays
    bounded as the map of the clusters does: each edit session starts from a
    fresh copy of the map, so the map holds at most the space that the edits
    of one run left behind (about 28 kB here, which the second run adds).
    From the second run on, it keeps its size, but for the history of the
    runs that its metadata gains (two annotate entries a run: the columns of
    the network dropped, then added again), which the edits of a run
    rewrite: two pages in six runs."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, metagenome_inputs["pairs"])
    names = metacooler.MetaCooler(str(metacool)).contigs["name"]
    hits = mtn.contig_hits(mtn.contig_contacts(str(metacool)), len(names))
    table = mtn.contig_data_table(assembly_contig_data(metagenome_inputs), names, hits)
    sizes = []
    for _ in range(6):
        mtmc.write_map(str(metacool), str(metacool), network=mtmc.Network(table, "empirical_hit"), force=True)
        sizes.append(metacool.stat().st_size)
    assert max(sizes[1:]) - sizes[1] < 12_000, sizes
    assert len(metacooler.MetaCooler(str(metacool)).metadata["annotate"]) == 1 + 2 * 5
    metacooler.validate(str(metacool), full=True)


def _one_round(metacool):
    """The tables of one clustering round, overlapping: one cluster of the
    first two contigs of the map."""
    contigs = list(metacooler.MetaCooler(str(metacool)).contigs.name[:2])
    clusters = pd.DataFrame({"round": ["overlapping"], "name": ["metator_00001_00000"], "parent_round": [None], "parent": [None]})
    members = pd.DataFrame({"round": "overlapping", "cluster": "metator_00001_00000", "contig": contigs})
    return clusters, members


@pytest.mark.parametrize("umask", [0o022, 0o077])
@pytest.mark.parametrize("writer", ["reassign_metacool", "write_rounds"])
def test_a_new_map_has_the_mode_of_a_new_file(metagenome, binning, tmp_path, writer, umask):
    """A map written to a new file (--output-map) gets the mode of any new file, from the umask, not the
    private mode of the file staged for it (mkstemp's 0600): write_rounds,
    which copies the map into that file, wrote it private. reassign_metacool
    is a guard: metacooler writes its map next to the staged file and
    renames it over it, a new file; so is umask 077, under which both
    agree."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    outfile = tmp_path / "clusters.metacool"
    previous = os.umask(umask)
    try:
        if writer == "reassign_metacool":
            reassign(
                str(metacool), str(binning["contig_data"]), str(outfile), bin_summary=str(binning["bin_summary"]), force=True
            )
        else:
            write_rounds(str(metacool), str(outfile), _one_round(metacool))
    finally:
        os.umask(previous)
    assert outfile.stat().st_mode & 0o777 == 0o666 & ~umask


@pytest.mark.parametrize("table", ["no Size column", "no file"])
def test_the_contigs_of_a_map_against_a_table_that_cannot_be_read(metagenome, binning, tmp_path, table):
    """A contig data file without a Name or a Size column, or that is not
    there, is refused by name: the map cannot be checked against it."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    contig_data = tmp_path / "contig_data.txt"
    if table == "no Size column":
        contig_data.write_text("ID\tName\n1\tc1\n")
    with pytest.raises(mtmc.MetacoolError, match=f"^Cannot read the contigs of {re.escape(str(contig_data))}: "):
        mtmc.check_map_of_contigs(str(metacool), str(contig_data))


def test_rounds_name_the_map_that_lacks_their_contigs(metagenome, binning, tmp_path, caplog):
    """The contigs of the rounds that the map lacks are left out, with a
    warning that names the map, not the file staged for it."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    clusters, members = _one_round(metacool)
    members = pd.concat([members, pd.DataFrame({"round": ["overlapping"], "cluster": ["metator_00001_00000"],
                                                "contig": ["not_a_contig"]})])  # fmt: skip
    with caplog.at_level("WARNING"):
        write_rounds(str(metacool), str(metacool), (clusters, members), force=True)
    assert f"1 contig(s) of the clusters of the run are not in {metacool}, and not in its clusters either" in caplog.text
    assert ".tmp" not in caplog.text


def test_rounds_written_in_place_keep_the_size_of_the_map(metagenome, binning, tmp_path):
    """Rounds written in place, run after run, keep the size of the map from
    the second run on: each edit session starts from a fresh copy of the
    map, so the map holds the space of the rounds that the last run replaced
    at most, which HDF5 does not reclaim. Replaced in the map itself, they
    left their space in the file at every run: the map grew with each run."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    sizes = []
    for _ in range(6):
        write_rounds(str(metacool), str(metacool), _one_round(metacool), force=True)
        sizes.append(metacool.stat().st_size)
    assert max(sizes[1:]) - sizes[1] < 1024, sizes


def test_rounds_written_to_a_copy_of_a_map(metagenome, binning, tmp_path):
    """write_rounds writes the map with the rounds of a run, and the rest of
    the map as stored (an h5py walk), in place of the rounds that it held;
    from an entry of a collection, the map of its own that the map of the
    entry gives. The map it was written from is left as it was."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    metacooler.write_clusters(
        str(metacool),
        ["earlier"],
        pd.DataFrame({"round": ["earlier"], "name": ["x"]}),
        pd.DataFrame({"round": ["earlier"], "cluster": ["x"], "contig": [metacooler.MetaCooler(str(metacool)).contigs.name[0]]}),
    )
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    before = metacool.read_bytes()
    from_map, from_entry = tmp_path / "from_map.metacool", tmp_path / "from_entry.metacool"
    write_rounds(str(metacool), str(from_map), _one_round(metacool), names=["overlapping"])
    write_rounds(f"{collection}::entries/lib1", str(from_entry), _one_round(metacool))
    assert metacool.read_bytes() == before
    for path in [from_map, from_entry]:
        metacooler.validate(str(path), full=True)
        mc = metacooler.MetaCooler(str(path))
        assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping"]
        assert list(mc.clusters.name) == ["metator_00001_00000"]
    found, source = walk(from_map), walk(metacool)
    assert {name: value for name, value in found.items() if not name.startswith("/clusters")} == {
        name: value for name, value in source.items() if not name.startswith("/clusters")
    }
    assert walk(from_entry) == found


@pytest.mark.parametrize("case", ["the collection of its entry", "an existing file", "a group"])
def test_rounds_written_to_a_map_refused(case, metagenome, binning, tmp_path):
    """write_rounds refuses, before any work, to write over the collection of
    an entry (the whole collection would be replaced), over an existing
    file without force (the map itself among them), and to a group of a
    file (file::group, a URI)."""
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, binning["pairs"])
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    entry = f"{collection}::entries/lib1"
    source, outfile, message = {
        "the collection of its entry": (
            entry,
            collection,
            f"{collection} is the .mmetacool collection of the input map {entry}: the map would replace the whole "
            "collection. Write it to another file.",
        ),
        "an existing file": (metacool, metacool, f"{metacool} already exists. Use --force to overwrite it."),
        "a group": (metacool, f"{metacool}::/", f"{metacool}::/ is not the path of a file"),
    }[case]
    before = {path: path.read_bytes() for path in [metacool, collection]}
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}"):
        write_rounds(str(source), str(outfile), _one_round(metacool))
    assert {path: path.read_bytes() for path in before} == before
    assert sorted(path.name for path in tmp_path.iterdir() if path.name.startswith(".")) == []


def _contig_map_and_final_tables(metagenome, inputs, tmp_path, measured=True):
    """The map of the contigs that a pipeline run starts from, made from its
    pairs (see tests.maps.contig_map), and its final contig data (the hits
    counted on the pairs, as the network step counts them). Returns the
    map, the contig data file and the contig data, whose hits are those of
    the file."""
    contig_data = assembly_contig_data(inputs, measured)
    outfile = tmp_path / "metator.metacool"
    contig_map(inputs["pairs"], contig_data, outfile, assembly_name="synthetic")
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


@pytest.mark.parametrize("normalization", ["empirical_hit", "length", "RS", "abundance", "theoritical_hit"])
def test_map_of_the_bins_stores_the_normalization_of_the_network(normalization, metagenome, metagenome_inputs, tmp_path, logged):
    """Written over the map of the contigs, as a pipeline run does at its
    end, the map of the bins stores the normalization of the network, as
    the weights metator_<normalization> of its bins, from the columns of its
    contigs table that the normalization reads, which it records. The
    contacts of the map normalized by them and summed by pair of contigs
    are metaTOR's network: the normalization of metaTOR 1.x of the contact
    count of each pair (see tests.oracles.normalize_pair; none where it gives
    0)."""
    outfile, contig_data_file, contig_data = _contig_map_and_final_tables(metagenome, metagenome_inputs, tmp_path)
    reassign(str(outfile), str(contig_data_file), str(outfile), force=True, normalization=normalization)

    metacooler.validate(str(outfile), full=True)
    name = f"metator_{normalization}"
    (row,) = _weights(outfile).values()
    attrs = row["attrs"]
    assert (row["name"], row["method"], attrs["mode"], attrs["metator_version"]) == (name, "metator", normalization, __version__)
    columns = [f"contigs/{column}" for column in mtmc.NORMALIZATION_COLUMNS[normalization]]
    assert attrs["depends_on"] == (["pixels"] if "contigs/Hit" in columns else []) + columns
    expected = {
        pair: normalize_pair(contig_data, pair, count, normalization)
        for pair, count in contig_contacts(outfile).items()
        if pair[0] != pair[1]
    }
    normalized = _normalized_contacts(outfile, name)
    assert normalized.keys() == {pair for pair, value in expected.items() if value}
    assert all(normalized[pair] == pytest.approx(expected[pair], rel=1e-9) for pair in normalized)
    assert not [message for level, message in logged() if level == "WARNING"]


def test_map_of_the_bins_at_1kb_stores_the_weight_of_each_contig_on_its_bins(metagenome, metagenome_inputs, tmp_path):
    """From a map of the contigs at a fixed bin size, made before any
    binning: every bin of a contig has the weight of its contig."""
    contig_data = assembly_contig_data(metagenome_inputs)
    metacool = unbinned_map(tmp_path / "contigs.metacool", metagenome, metagenome_inputs["pairs"], binsize=1000)
    contig_data_file = final_contig_data(
        contig_data, metagenome.pairs, _bin_names(metagenome), tmp_path / "contig_data_final.txt"
    )
    out = tmp_path / "bins.metacool"
    reassign(str(metacool), str(contig_data_file), str(out), normalization="empirical_hit")

    weights = metacooler.read_weights(str(out), "metator_empirical_hit")
    bins = cooler.Cooler(str(out)).bins()[:]
    assert len(bins) > len(contig_data)
    hits = bins.chrom.astype(str).map({name: data["hit"] for name, data in contig_data.items()}).to_numpy(dtype=float)
    with np.errstate(divide="ignore"):
        expected = 1 / np.sqrt(hits)
    expected[~np.isfinite(expected)] = np.nan
    np.testing.assert_allclose(weights, expected)


def test_map_of_the_bins_not_written_when_its_normalization_cannot_be_stored(
    metagenome, metagenome_inputs, tmp_path, monkeypatch
):
    """A normalization that metacooler fails to store fails the map: in an
    edit session of metacooler, a write that fails fails the session, so
    the map is not written, with the reason, and the input map is left as
    it was. (A normalization whose columns the run has not is refused by
    metator network before any work.)"""
    outfile, contig_data_file, _ = _contig_map_and_final_tables(metagenome, metagenome_inputs, tmp_path)
    before = outfile.read_bytes()

    def refused(*args, **kwargs):
        raise ValueError("weights refused")

    monkeypatch.setattr(metacooler, "store_weights", refused)
    with pytest.raises(mtmc.MetacoolError, match=r"^weights refused$"):
        reassign(str(outfile), str(contig_data_file), str(outfile), force=True, normalization="length")
    assert outfile.read_bytes() == before
    # No lock, and no staging directory (.<map>.<uuid>.<n>.d) left.
    assert sorted(
        path.name for path in outfile.parent.iterdir() if path.name.startswith((outfile.name, f".{outfile.name}."))
    ) == [outfile.name]


@pytest.mark.parametrize("normalization", ["abundance", "None"])
def test_map_of_the_bins_replaces_the_normalization_of_an_earlier_run(
    normalization, metagenome, metagenome_inputs, tmp_path, logged
):
    """From the map of an earlier run, which stores the normalization of its
    network: the map of the bins stores this run's alone, or none with
    "None", without a warning. A weight that another tool stored is kept."""
    earlier, contig_data_file, _ = _contig_map_and_final_tables(metagenome, metagenome_inputs, tmp_path)
    reassign(str(earlier), str(contig_data_file), str(earlier), force=True, normalization="empirical_hit")
    n_bins = len(cooler.Cooler(str(earlier)).bins()[:])
    metacooler.store_weights(str(earlier), "other", np.ones(n_bins), method="other", depends_on=["pixels"])
    out = tmp_path / "bins.metacool"
    reassign(str(earlier), str(contig_data_file), str(out), normalization=normalization)

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
        # The map the weights were computed from: metacooler refuses them if
        # the map changed since.
        "source": stored["source"],
    }
    assert isinstance(stored["source"], metacooler.MetaCooler) and str(stored["source"].uri).startswith(str(outfile))
    for normalization in ["None", "bogus"]:
        with pytest.raises(mtmc.MetacoolError, match=f"No weights for the normalization '{normalization}'"):
            mtmc.store_normalization(str(outfile), normalization)
    metacooler.annotate(str(outfile), drop=["contigs.Shotgun_coverage", "contigs.Restriction_site"], inplace=True)
    with pytest.raises(mtmc.MetacoolError, match="reads the Shotgun_coverage and Restriction_site column"):
        mtmc.store_normalization(str(outfile), "theoritical_hit")


def _three_contigs(path):
    """A map of the contigs a (300 bp), b (200 bp) and c (100 bp), in that
    order, every one in __unbinned__, one bin each."""
    import metacooler

    pairs = path.with_suffix(".pairs")
    pairs.write_text("## pairs format v1.0\nr1\ta\t10\tb\t20\t+\t-\nr2\tb\t10\tc\t90\t+\t-\n")
    bins = pd.DataFrame({"chrom": ["a", "b", "c"], "start": 0, "end": [300, 200, 100]})
    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    pairs.unlink()
    return path


def _stored_rounds(path, clusters):
    """Store clusters in the map at path as metacooler stores clustering
    rounds: by round, each cluster's name and contigs."""
    import metacooler

    rounds = list(dict.fromkeys(round_name for round_name, _, _ in clusters))
    metacooler.write_clusters(
        str(path),
        rounds,
        pd.DataFrame({"round": [round_name for round_name, _, _ in clusters], "name": [name for _, name, _ in clusters]}),
        pd.DataFrame(
            [(round_name, name, contig) for round_name, name, contigs in clusters for contig in contigs],
            columns=["round", "cluster", "contig"],
        ),
    )


def test_read_partition_gives_the_partition_back(tmp_path):
    """The partition that a map stores, read back with the IDs of the
    contigs, their positions in the map: the core clusters by number, their
    contigs by ID; the overlapping clusters by number, named as stored, their
    contigs in the order of their core clusters, then of their IDs (that of
    the FASTA file of the cluster, which the rounds do not keep). The
    partition stored by metator partition comes back as it was."""
    import metator.clusters as mtcl

    metacool = _three_contigs(tmp_path / "map.metacool")
    _stored_rounds(
        metacool,
        [
            ("core", "core_00002", ["a", "c"]),
            ("core", "core_00001", ["b"]),
            ("overlapping", "my_sample_metator_00003_00000", ["a", "b", "c"]),
        ],
    )
    partition = mtmc.read_partition(str(metacool))
    assert partition.core == {1: [2], 2: [1, 3]}
    assert partition.overlapping == {3: [2, 1, 3]}
    assert partition.names == {3: "my_sample_metator_00003_00000"}
    # As metator partition stores it, and back.
    stored = mtcl.Partition(
        core={1: [1, 3], 2: [2]}, overlapping={1: [2], 7: [1, 3]}, names={1: "m_00001_00000", 7: "m_00007_00000"}
    )
    write_rounds(
        str(metacool),
        str(metacool),
        mtcl.partition_rounds(stored, ["a", "b", "c"]),
        names=list(mtmc.PARTITION_ROUNDS),
        force=True,
    )
    assert mtmc.read_partition(str(metacool)) == stored


@pytest.mark.parametrize(
    "case",
    ["no round", "no overlapping", "a core cluster of another name", "two of one number", "a contig in no core cluster"],
)
def test_read_partition_refused(case, tmp_path):
    """A map that holds no partition, or one that metator partition did not
    store (its clusters named otherwise, or an overlapping cluster with a
    contig in no core cluster), is refused, with what to run."""
    metacool = _three_contigs(tmp_path / "map.metacool")
    clusters = {
        "no round": [],
        "no overlapping": [("core", "core_00001", ["a"])],
        "a core cluster of another name": [("core", "core_1", ["a"]), ("overlapping", "m_00001_00000", ["a"])],
        "two of one number": [
            ("core", "core_00001", ["a"]),
            ("core", "core_000001", ["b"]),
            ("overlapping", "m_00001_00000", ["a", "b"]),
        ],
        "a contig in no core cluster": [("core", "core_00001", ["a"]), ("overlapping", "m_00001_00000", ["a", "c"])],
    }[case]
    if clusters:
        _stored_rounds(metacool, clusters)
    message = {
        "no round": "holds no partition of its network (the clustering round(s) core and overlapping): run metator partition",
        "no overlapping": "holds no partition of its network (the clustering round(s) overlapping): run metator partition",
        "a core cluster of another name": "The clustering round core of",
        "two of one number": "The clustering round core of",
        "a contig in no core cluster": "hold contig(s) in no core cluster: 'c'.",
    }[case]
    with pytest.raises(mtmc.MetacoolError, match=re.escape(message)):
        mtmc.read_partition(str(metacool))


def test_read_partition_of_the_network_of_the_map(tmp_path):
    """With the digest of the network of the map, read_partition reads the
    partition only when the map records that it was made on that network,
    as metator partition records it with the partition: a partition without
    that record, or with that of another network, is refused, with the step
    to run again. holds_partition says whether the map holds a partition at
    all, and partition_network gives the network that it records."""
    import metacooler

    import metator.clusters as mtcl

    metacool = _three_contigs(tmp_path / "map.metacool")
    assert not mtmc.holds_partition(metacooler.MetaCooler(str(metacool)))
    clusters = mtcl.Partition(core={1: [1, 2, 3]}, overlapping={1: [1, 2, 3]}, names={1: "m_00001_00000"})
    for entries, made_on in [
        (None, None),
        ({"partition": {"version": "x"}}, None),
        ({"partition": {mtmc.NETWORK_DIGEST: "f00"}}, "f00"),
    ]:
        write_rounds(
            str(metacool),
            str(metacool),
            mtcl.partition_rounds(clusters, ["a", "b", "c"]),
            names=list(mtmc.PARTITION_ROUNDS),
            force=True,
            entries=entries,
        )
        mc = metacooler.MetaCooler(str(metacool))
        assert mtmc.holds_partition(mc) and mtmc.partition_network(mc) == made_on
        assert mtmc.read_partition(str(metacool)) == clusters
    assert mtmc.read_partition(str(metacool), network="f00") == clusters
    with pytest.raises(mtmc.MetacoolError) as refused:
        mtmc.read_partition(str(metacool), network="ba5")
    assert str(refused.value) == (
        f"{metacool} holds a partition of another network than the one stored in it now: metator network stored another "
        f"network in the map since metator partition ran on it. Run metator partition {metacool} again, which "
        "partitions the network stored in the map."
    )
    write_rounds(
        str(metacool),
        str(metacool),
        mtcl.partition_rounds(clusters, ["a", "b", "c"]),
        names=list(mtmc.PARTITION_ROUNDS),
        force=True,
        entries={"partition": "not an entry"},
    )
    with pytest.raises(mtmc.MetacoolError, match=re.escape(f"{metacool} does not record the network that its partition")):
        mtmc.read_partition(str(metacool), network="f00")


def test_the_partition_entry_is_kept_in_the_map_of_the_clusters(tmp_path):
    """The map of the clusters that a validation writes keeps the entry
    partition of the metator metadata of its input map, which metator
    partition stored there, as it keeps the entry network; a pipeline run
    stores its own."""
    import metacooler

    import metator.clusters as mtcl

    metacool = _three_contigs(tmp_path / "map.metacool")
    clusters = mtcl.Partition(core={1: [1, 2]}, overlapping={1: [1, 2]}, names={1: "metator_00001_00000"})
    entry = {"version": "x", "normalization": "None"}
    write_rounds(
        str(metacool),
        str(metacool),
        mtcl.partition_rounds(clusters, ["a", "b", "c"]),
        names=list(mtmc.PARTITION_ROUNDS),
        force=True,
        entries={"partition": entry},
    )
    assert metacooler.MetaCooler(str(metacool)).metadata["metator"] == {"partition": entry}
    table = tmp_path / "contig_data_final.txt"
    table.write_text("Name\tSize\tFinal_cluster\na\t300\tmetator_00001_00000\nb\t200\tmetator_00001_00000\nc\t100\tND\n")
    for given, expected in [(None, entry), ({"version": "y"}, {"version": "y"})]:
        written = tmp_path / f"clusters_{given is None}.metacool"
        reassign(str(metacool), str(table), str(written), partition=given)
        assert metacooler.MetaCooler(str(written)).metadata["metator"]["partition"] == expected


def write_rounds(metacool, outfile, rounds, names=None, force=False, planned=None, entries=None):
    """Write the map metacool at outfile with these clustering rounds in place
    of those it held, as metator partition writes it (see
    metator.metacool.write_map), with the entry partition of entries."""
    partition = (entries or {}).get("partition")
    mtmc.write_map(
        str(metacool), str(outfile), partition=partition, rounds=rounds, round_names=names, force=force, planned=planned
    )
