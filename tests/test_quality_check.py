# Test for quality_check module

import gzip
import locale
import re
import shutil

import Bio
import metator.io as mio
import metator.metacool as mtmc
import metator.quality_check as mtq
import pandas as pd
import pytest
from metator.log import logger
from tests import synthetic

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
    tuple(int(x) for x in Bio.__version__.split(".")[:2]) < (1, 87),
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


# Quality check on a mock metagenome (tests/synthetic.py) of two genomes
# assembled in contigs larger than 100 kb, the size the quality check uses.


@pytest.fixture(scope="module")
def large_contigs(tmp_path_factory):
    """Assembly, sorted pairs and metator outputs (contig data, bin summary)
    of two genomes in large contigs; a third, low-quality bin has none."""
    genomes = (
        synthetic.GenomeSpec("big_1", gc=0.40, abundance=1.0, length=260_000),
        synthetic.GenomeSpec("big_2", gc=0.60, abundance=2.0, length=240_000),
    )
    community = synthetic.make_metagenome(
        seed=11, genomes=genomes, mges=(), n_pairs=30_000, contig_median=120_000, contig_bounds=(105_000, 150_000)
    )
    directory = tmp_path_factory.mktemp("large_contigs")
    files = community.write_all(directory)
    pairs = mio.sort_and_index_pairs(str(files["pairs"]), tmp_dir=str(directory))
    truth = community.truth
    bins = {"big_1": "metator_00001_00000", "big_2": "metator_00002_00000"}
    contig_data = pd.DataFrame(
        {"ID": range(1, len(truth) + 1), "Name": truth.Name, "Size": truth.Length, "Final_cluster": truth.Genome.map(bins)}
    )
    contig_data.to_csv(directory / "contig_data_final.txt", sep="\t", index=False)
    summary = pd.DataFrame(
        {"Weighted completeness": [0.95, 0.92, 0.5], "Weighted redundancy": [1.01, 1.05, 1.0]},
        index=[*bins.values(), "metator_00003_00000"],
    )
    summary.to_csv(directory / "cluster_summary.txt", sep="\t")
    return community, files, pairs, directory


def test_extract_hq_contigs(large_contigs):
    directory = large_contigs[3]
    contigs_data = mio.read_contig_data(str(directory / "contig_data_final.txt"))
    summary = mio.read_bin_summary(str(directory / "cluster_summary.txt"))
    hq_contigs, n_mags = mtq.extract_hq_contigs(summary, contigs_data)
    assert n_mags == 2
    assert hq_contigs == dict(zip(contigs_data.Name, contigs_data.Final_cluster, strict=True))
    # A contaminated MAG is not high quality.
    summary.loc["metator_00002_00000", "Weighted redundancy"] = 1.3
    hq_contigs, n_mags = mtq.extract_hq_contigs(summary, contigs_data)
    assert n_mags == 1
    assert set(hq_contigs.values()) == {"metator_00001_00000"}


def test_hic_quality_counts(large_contigs, tmp_path):
    """Every extracted pair is counted once, in the right category."""
    community, files, pairs, directory = large_contigs
    contigs_data = mio.read_contig_data(str(directory / "contig_data_final.txt"))
    names = list(contigs_data.Name)
    n_pairs = mtq.extract_pairs([str(pairs)], str(tmp_path / "hq.pairs"), names, contigs_data)
    assert n_pairs == len(community.pairs)
    contigs = dict(zip(names, contigs_data.Final_cluster, strict=True))
    counts = mtq.hic_quality(
        contigs,
        2,
        str(tmp_path / "hq.pairs"),
        str(tmp_path / "hq_idx.pairs"),
        str(files["assembly"]),
        ["DpnII"],
        threshold=["3", "3"],
    )
    religated, loops, weirds, informative_intra, informative_inter, intra_mags, inter_mags = counts
    simulated = community.pairs
    genome = community.truth.Genome
    same_contig = simulated.chr1 == simulated.chr2
    same_genome = simulated.chr1.map(genome) == simulated.chr2.map(genome)
    assert intra_mags + inter_mags == n_pairs
    assert inter_mags == (~same_genome).sum()
    assert informative_inter == (same_genome & ~same_contig).sum()
    assert religated + loops + weirds + informative_intra == same_contig.sum()
    assert religated > 0 and loops > 0 and weirds > 0
    # With a single MAG there is no noise score to compute, the counts stay.
    assert (
        mtq.hic_quality(
            contigs,
            1,
            str(tmp_path / "hq.pairs"),
            str(tmp_path / "idx2.pairs"),
            str(files["assembly"]),
            ["DpnII"],
            threshold=["3", "3"],
        )
        == counts
    )


def test_hic_quality_without_pairs(large_contigs, tmp_path, caplog):
    """No pair on the contigs of the high-quality MAGs: every count, score
    and ratio is 0, with a warning, not a ZeroDivisionError."""
    _, files, _, directory = large_contigs
    contigs_data = mio.read_contig_data(str(directory / "contig_data_final.txt"))
    names = list(contigs_data.Name)
    assert mtq.extract_pairs([], str(tmp_path / "hq.pairs"), names, contigs_data) == 0
    contigs = dict(zip(names, contigs_data.Final_cluster, strict=True))
    camembert = tmp_path / "camembert.pdf"
    with caplog.at_level("INFO", logger="metator"):
        counts = mtq.hic_quality(
            contigs, 2, str(tmp_path / "hq.pairs"), str(tmp_path / "hq_idx.pairs"), str(files["assembly"]), ["DpnII"],
            plot=True, plot_cam=str(camembert), threshold=["3", "3"],
        )  # fmt: skip
    assert counts == (0, 0, 0, 0, 0, 0, 0)
    assert "No pairs have been extracted" in caplog.text
    assert "Noise score: 0.00E+00" in caplog.text and "Ratio inter/intra contigs: 0.00%." in caplog.text
    assert f"No pair to plot: {camembert} is not drawn." in caplog.text and not camembert.exists()


@pytest.mark.parametrize(
    "case, message",
    [
        ("no MAG", "No high-quality MAG (weighted completeness above 70%, weighted redundancy under 1.15): no quality check."),
        ("no large contig", "The 2 high-quality MAGs have no contig longer than 100 kb: no quality check."),
        ("no pair", "No pair between the contigs longer than 100 kb of the high-quality MAGs: no quality check."),
    ],
)
def test_quality_check_without_pairs_says_why(case, message, large_contigs, tmp_path, caplog):
    """Without pairs to check, the log says why: it said "No MAGs have been
    found." whatever the reason."""
    _, files, _, directory = large_contigs
    contig_data = pd.read_csv(directory / "contig_data_final.txt", sep="\t")
    summary = pd.read_csv(directory / "cluster_summary.txt", sep="\t", index_col=0)
    if case == "no MAG":
        summary["Weighted completeness"] = 0.4
    if case == "no large contig":
        contig_data["Size"] = 50_000
    contig_data.to_csv(tmp_path / "contig_data.txt", sep="\t", index=False)
    summary.to_csv(tmp_path / "cluster_summary.txt", sep="\t")
    with caplog.at_level("INFO", logger="metator_logger"):
        mtq.quality_check(str(tmp_path / "contig_data.txt"), str(tmp_path / "cluster_summary.txt"), str(files["assembly"]),
                          [], str(tmp_path), str(tmp_path), "sample", False, ["DpnII"], None)  # fmt: skip
    assert "0 pairs have been extracted." in caplog.text
    assert message in caplog.text


def test_qc_command(large_contigs, tmp_path, run_cli):
    """From the metator output directory, with automatic thresholds and plots."""
    _, files, pairs, directory = large_contigs
    out = tmp_path / "qc"
    args = ["-a", str(files["assembly"]), "-e", "DpnII", "-O", str(directory), "-o", str(out), "-p", "sample", "-P",
            "-T", str(tmp_path / "tmp"), str(pairs)]  # fmt: skip
    run_cli("qc", *args)
    assert (out / "sample_camembert_plot.pdf").stat().st_size > 0
    assert (out / "sample_event.pdf").stat().st_size > 0
    assert not (tmp_path / "tmp").exists()


def test_qc_command_without_hq_mags(large_contigs, tmp_path, run_cli):
    """Explicit tables; no high-quality MAG, so nothing is plotted."""
    _, files, pairs, directory = large_contigs
    summary = mio.read_bin_summary(str(directory / "cluster_summary.txt"))
    summary["Weighted completeness"] = 0.4
    summary.to_csv(tmp_path / "cluster_summary.txt", sep="\t")
    args = ["-a", str(files["assembly"]), "-e", "DpnII", "-b", str(tmp_path / "cluster_summary.txt"),
            "-c", str(directory / "contig_data_final.txt"), "-t", "3,3", "-T", str(tmp_path / "tmp"), "-n",
            str(pairs)]  # fmt: skip
    run_cli("qc", *args)
    assert not list((tmp_path / "metator_qc").glob("*.pdf"))
    assert (tmp_path / "tmp").exists()  # --no-clean-up


@pytest.mark.parametrize("run", ["1.x", "both summaries"])
def test_qc_command_in_the_directory_of_a_run_of_1x(run, large_contigs, tmp_path, run_cli, logged):
    """-O finds the bin_summary.txt of a run of metaTOR 1.x, and its contig
    data, its columns named after bins, is read under their names of 2.0:
    qc measures what it measures in the directory of the same run of 2.0.
    The cluster_summary.txt of a directory that has both is read."""
    _, files, pairs, directory = large_contigs
    metator_dir = tmp_path / "metator"
    metator_dir.mkdir()
    if run == "1.x":
        shutil.copyfile(directory / "cluster_summary.txt", metator_dir / "bin_summary.txt")
        synthetic.as_1x(directory / "contig_data_final.txt", metator_dir / "contig_data_final.txt")
    else:
        shutil.copyfile(directory / "cluster_summary.txt", metator_dir / "cluster_summary.txt")
        (metator_dir / "bin_summary.txt").write_text("not a summary\n")
        shutil.copyfile(directory / "contig_data_final.txt", metator_dir / "contig_data_final.txt")
    measured = []
    for path in (directory, metator_dir):
        start = len(logged())
        args = ["-a", files["assembly"], "-e", "DpnII", "-O", path, "-o", tmp_path / path.name, "-T", tmp_path / "tmp", pairs]
        run_cli("qc", *map(str, args))
        measured.append([message for _, message in logged()[start:] if re.match(r"(Number|Contacts|Noise|\d+ pairs)", message)])
    assert measured[0] == measured[1]
    assert any(message.startswith("Noise contact ratio") for message in measured[0])


@pytest.mark.parametrize("missing", ["cluster_summary.txt", "contig_data_final.txt"])
def test_qc_command_missing_tables(large_contigs, tmp_path, missing, run_cli):
    _, files, pairs, directory = large_contigs
    metator_dir = tmp_path / "metator"
    metator_dir.mkdir()
    for name in {"cluster_summary.txt", "contig_data_final.txt"} - {missing}:
        shutil.copyfile(directory / name, metator_dir / name)
    args = ["-a", str(files["assembly"]), "-e", "DpnII", "-O", str(metator_dir), "-T", str(tmp_path / "tmp"), str(pairs)]
    with pytest.raises(FileNotFoundError):
        run_cli("qc", *args)


# The quality check from a .metacool contact map, in place of the pairs files:
# the noise of the library only.


def _edge_pairs(community):
    """Five pairs on the edges of the contigs of community (two contigs of a
    genome, a and b, and a contig of another genome, c), in the upper
    triangle, as metator.align writes them: a read at position 0 (which a
    1-based pairs file should not hold, but some do) within a and between a
    and c, a read on the last base of a with one on the first base of b,
    and a read past the end of a within a, and of b with c. The map and the
    pairs path count the first three, and neither the last two (see
    test_noise_from_a_map_is_that_of_its_pairs)."""
    rank = {contig.name: i for i, contig in enumerate(community.contigs)}
    length = {contig.name: contig.length for contig in community.contigs}
    genome = dict(zip(community.truth.Name, community.truth.Genome, strict=True))
    first = genome[community.contigs[0].name]
    a, b = [contig.name for contig in community.contigs if genome[contig.name] == first][:2]
    c = next(contig.name for contig in community.contigs if genome[contig.name] != first)

    def pair(read, contig1, pos1, contig2, pos2):
        if rank[contig1] > rank[contig2]:
            contig1, pos1, contig2, pos2 = contig2, pos2, contig1, pos1
        return f"{read}\t{contig1}\t{pos1}\t{contig2}\t{pos2}\t+\t-\n"

    return [
        pair("edge_1", a, 0, a, 2000),
        pair("edge_2", a, length[a], b, 1),
        pair("edge_3", a, length[a], c, 0),
        pair("edge_4", a, 100, a, length[a] + 1),
        pair("edge_5", b, length[b] + 50, c, 7),
    ]


@pytest.fixture(scope="module")
def large_contigs_map(large_contigs, tmp_path_factory):
    """The pairs of large_contigs with the pairs on the edges of its contigs
    (_edge_pairs), sorted and indexed for the pairs path, and the .metacool
    map that metacooler builds from them, one bin per contig, every contig
    in no genome yet."""
    import metacooler

    community, files, _, _ = large_contigs
    directory = tmp_path_factory.mktemp("large_contigs_map")
    pairs = directory / "edges.pairs"
    shutil.copyfile(files["pairs"], pairs)
    with open(pairs, "a") as handle:
        handle.writelines(_edge_pairs(community))
    indexed = mio.sort_and_index_pairs(str(pairs), tmp_dir=str(directory))
    truth = community.truth
    metacool = directory / "map.metacool"
    bins = pd.DataFrame({"chrom": truth.Name, "start": 0, "end": truth.Length})
    metacooler.from_pairs(str(metacool), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    return indexed, metacool


def _metrics(messages):
    """The messages of the noise, which a map gives: the number of MAGs, the
    pairs within contigs, within MAGs and between MAGs, the noise contact
    ratio and the noise score."""
    return [message for message in messages if message.startswith(("Number of MAGs", "Contacts ", "Noise "))]


def test_noise_from_a_map_is_that_of_its_pairs(large_contigs, large_contigs_map, tmp_path, logged):
    """The noise that qc measures from a .metacool map (extract_map_pairs,
    map_quality) is the noise that it measures from the pairs files that
    the map was made of (extract_pairs, hic_quality): the pairs within
    contigs (the uncuts, loops, weirds and informative pairs within a
    contig of the pairs path), within MAGs and between MAGs, then the noise
    contact ratio and the noise score, logged line for line as from the
    pairs. The events are not measured from the map, which the log says.

    What is compared: the pairs that the map holds, those whose two reads
    lie on the contigs, within their length. The pairs files hold five
    pairs on the edges of the contigs (see _edge_pairs): a read at position
    0, which metacooler puts on the first base and the pairs path queries
    from 0, and a read on the last base count in both; a read past the end
    of its contig counts in neither, as metacooler leaves it out and the
    pairs path queries a contig up to its length. The pairs are in the
    upper triangle, as metator.align writes them: the pairs path reads one
    order of each two contigs, where the map holds both."""
    community, files, _, directory = large_contigs
    indexed, metacool = large_contigs_map
    contigs_data = mio.read_contig_data(str(directory / "contig_data_final.txt"))
    names = list(contigs_data.Name)
    contigs = dict(zip(names, contigs_data.Final_cluster, strict=True))
    n_pairs = mtq.extract_pairs([str(indexed)], str(tmp_path / "hq.pairs"), names, contigs_data)
    events = mtq.hic_quality(contigs, 2, str(tmp_path / "hq.pairs"), str(tmp_path / "hq_idx.pairs"), str(files["assembly"]),
                             ["DpnII"], threshold=["3", "3"])  # fmt: skip
    from_pairs = [message for _, message in logged()]
    contacts = mtq.extract_map_pairs(str(metacool), names, contigs_data, "contig_data_final.txt")
    counts = mtq.map_quality(contigs, 2, contacts)
    from_map = [message for _, message in logged()][len(from_pairs) :]
    religated, loops, weirds, informative_intra, informative_inter, intra_mags, inter_mags = events
    assert n_pairs == contacts["count"].sum() == len(community.pairs) + 3
    assert counts == (religated + loops + weirds + informative_intra, intra_mags, inter_mags)
    assert informative_inter == intra_mags - counts[0] and inter_mags > 0
    assert len(_metrics(from_map)) == 6 and _metrics(from_map) == _metrics(from_pairs)
    assert f"{n_pairs} pairs have been extracted from {metacool}." in from_map and mtq.NOT_FROM_A_MAP in from_map
    assert not [message for message in from_map if message.startswith(("Religated", "Loop", "Weirds", "Informative", "Ratio"))]


def test_noise_metrics():
    """The noise contact ratio and the noise score of hic_quality, without a
    pair, with one MAG, and with two (2 * 9 pairs within, 1 between)."""
    assert mtq.noise_metrics(2, 0, 0) == (0, 0)
    assert mtq.noise_metrics(1, 9, 1) == (10, 0)
    assert mtq.noise_metrics(2, 18, 2) == (10, pytest.approx(2 / (9 + 2)))


@pytest.mark.parametrize("given", ["a map", "an entry of a collection"])
def test_qc_command_from_a_map(given, large_contigs, large_contigs_map, tmp_path, run_cli):
    """metator qc MAP, the map in place of the pairs files (a file, or an
    entry of a .mmetacool collection), logs the noise that qc logs from the
    pairs files the map was made of, line for line (see
    test_noise_from_a_map_is_that_of_its_pairs), then what a map does not
    give; it writes nothing but its log. --assembly and --enzyme are not
    used."""
    import metacooler

    community, files, _, directory = large_contigs
    indexed, metacool = large_contigs_map
    if given == "an entry of a collection":
        metacooler.create_mmetacool(tmp_path / "maps.mmetacool", {"lib1": metacool})
        metacool = f"{tmp_path / 'maps.mmetacool'}::entries/lib1"
    common = ["-e", "DpnII", "-O", directory, "-T", tmp_path / "tmp"]
    run_cli("qc", *common, "-a", files["assembly"], "-o", tmp_path / "pairs", "-t", "3,3", indexed)
    # No assembly: it is not read.
    run_cli("qc", *common, "-a", tmp_path / "absent.fa", "-o", tmp_path / "map", metacool)
    (pairs_log,) = (tmp_path / "pairs").glob("metator_qc_*.log")
    (map_log,) = (tmp_path / "map").glob("metator_qc_*.log")

    def messages(log):
        return [line.split(" :: ", 2)[2] for line in log.read_text().splitlines() if " :: INFO :: " in line]

    assert len(_metrics(messages(map_log))) == 6 and _metrics(messages(map_log)) == _metrics(messages(pairs_log))
    assert f"{len(community.pairs) + 3} pairs have been extracted from {metacool}." in messages(map_log)
    assert mtq.NOT_FROM_A_MAP in messages(map_log) and "Religated ratio" not in map_log.read_text()
    assert sorted(path.name for path in (tmp_path / "map").iterdir()) == [map_log.name]
    assert not (tmp_path / "tmp").exists()


# What qc refuses of a map in place of the pairs files, the reason logged as
# an error: before any work, how the map is given; then, once the tables are
# read, a map that cannot give the contacts of the contigs of the run.
QC_MAP_ERRORS = {
    "a map and pairs files": (
        ["{map}", "{pairs}"],
        "metator qc reads pairs files, or one .metacool map: got 2 files, 1 of them .metacool maps.",
    ),
    "two maps": (
        ["{map}", "{map}"],
        "metator qc reads pairs files, or one .metacool map: got 2 files, 2 of them .metacool maps.",
    ),
    "--threshold": (
        ["-t", "3,3", "{map}"],
        "metator qc from the contact map {map} takes no -t/--threshold: they set and draw the events of the pairs within "
        "a contig (uncuts, loops, weirds), which need the positions and strands of their reads, which only pairs files "
        "hold.",
    ),
    "--threshold and --plot": (
        ["-P", "-t", "3,3", "{map}"],
        "metator qc from the contact map {map} takes no -t/--threshold or -P/--plot: ",
    ),
    "a collection": (["{collection}"], "Cannot read the contact map {collection}: "),
    "a map whose name holds ::": (
        ["{colon_map}"],
        "Cannot read the contact map {colon_map}: metacooler reads a path that holds :: as file::group. Rename the map "
        "without ::.",
    ),
    "other lengths": (
        ["{shorter}"],
        "{shorter} is not a contact map of the contigs of {contig_data}: 1 contig(s) have another length in it: "
        "'{contig}' ({length} bp in the map, {size} bp in {contig_data}).",
    ),
    "both triangles": (
        ["{square}"],
        "{square} stores its contacts with storage-mode 'square': the contacts of contigs are read from a map that "
        "stores each contact once, in the upper triangle (storage-mode 'symmetric-upper'), as metacooler writes them.",
    ),
}


@pytest.mark.parametrize("case", list(QC_MAP_ERRORS))
def test_qc_command_from_a_map_refused(case, large_contigs, large_contigs_map, tmp_path, run_cli):
    """qc refuses a map with pairs files or another map, --threshold and
    --plot, which set and draw the events, a map that metacooler cannot
    open (a collection, which the error names the entries of, or a map
    whose name holds ::), and a map that cannot give the contacts of the
    contigs of the run (one of another length, or both triangles stored):
    the reason logged as an error in the log file, as MetacoolError, and
    no noise measured."""
    import h5py
    import metacooler

    community, files, pairs, directory = large_contigs
    _, metacool = large_contigs_map
    contig = community.contigs[0]
    names = {
        "map": metacool,
        "pairs": pairs,
        "collection": tmp_path / "maps.mmetacool",
        "colon_map": tmp_path / "run::1.metacool",
        "shorter": tmp_path / "shorter.metacool",
        "square": tmp_path / "square.metacool",
        "contig_data": directory / "contig_data_final.txt",
        "contig": contig.name,
        "length": contig.length - 1000,
        "size": contig.length,
    }
    metacooler.create_mmetacool(names["collection"], {"lib1": metacool})
    shutil.copyfile(metacool, names["colon_map"])
    shutil.copyfile(metacool, names["square"])
    with h5py.File(names["square"], "r+") as group:
        group.attrs["storage-mode"] = "square"
    truth = community.truth
    bins = pd.DataFrame({"chrom": truth.Name, "start": 0, "end": truth.Length.where(truth.Name != contig.name, names["length"])})
    metacooler.from_pairs(str(names["shorter"]), pairs=str(files["pairs"]), bins=bins, contig_genome_map=None, ps=False)
    inputs, message = QC_MAP_ERRORS[case]
    message = message.format(**names)
    args = ["-a", files["assembly"], "-e", "DpnII", "-O", directory, "-o", tmp_path / "out", "-T", tmp_path / "tmp"]
    pattern = f"^{re.escape(message)}" + (".*entries/lib1" if case == "a collection" else "")
    with pytest.raises(mtmc.MetacoolError, match=pattern):
        run_cli("qc", *args, *[value.format(**names) for value in inputs])
    (log,) = (tmp_path / "out").glob("metator_qc_*.log")
    assert f"ERROR :: {message}" in log.read_text() and "Noise" not in log.read_text()


@pytest.mark.parametrize(
    "case, message",
    [
        ("no MAG", "No high-quality MAG (weighted completeness above 70%, weighted redundancy under 1.15): no quality check."),
        ("no large contig", "The 2 high-quality MAGs have no contig longer than 100 kb: no quality check."),
        ("no pair", "No pair between the contigs longer than 100 kb of the high-quality MAGs: no quality check."),
    ],
)
def test_quality_check_from_a_map_without_pairs_says_why(case, message, large_contigs, tmp_path, logged):
    """From a map, as from pairs files, without pairs to check, the log says
    why. A map that lacks contigs of the high-quality MAGs (here, every one
    of them: no pair) gives them no pairs, with a warning that names them."""
    import metacooler

    _, files, _, directory = large_contigs
    contig_data = pd.read_csv(directory / "contig_data_final.txt", sep="\t")
    summary = pd.read_csv(directory / "cluster_summary.txt", sep="\t", index_col=0)
    contigs = list(contig_data.Name) if case != "no pair" else ["x1", "x2"]
    pairs = tmp_path / "other.pairs"
    pairs.write_text(f"## pairs format v1.0\nr1\t{contigs[0]}\t10\t{contigs[1]}\t20\t+\t-\n")
    metacool = tmp_path / "map.metacool"
    bins = pd.DataFrame({"chrom": contigs, "start": 0, "end": 200_000})
    metacooler.from_pairs(str(metacool), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    if case == "no MAG":
        summary["Weighted completeness"] = 0.4
    if case == "no large contig":
        contig_data["Size"] = 50_000
    contig_data.to_csv(tmp_path / "contig_data.txt", sep="\t", index=False)
    summary.to_csv(tmp_path / "cluster_summary.txt", sep="\t")
    mtq.quality_check(str(tmp_path / "contig_data.txt"), str(tmp_path / "cluster_summary.txt"), str(files["assembly"]), [],
                      str(tmp_path), str(tmp_path), "sample", False, ["DpnII"], None, metacool=str(metacool))  # fmt: skip
    found = logged()
    assert ("INFO", f"0 pairs have been extracted from {metacool}.") in found and ("INFO", message) in found
    absent = [entry for entry in found if entry[0] == "WARNING"]
    if case == "no pair":
        assert absent == [
            (
                "WARNING",
                f"3 of the contigs longer than 100 kb of the high-quality MAGs are not in {metacool}: they have no "
                f"contacts there, as contigs that no read maps to: {', '.join(repr(name) for name in contig_data.Name)}.",
            )
        ]
    else:
        assert absent == []


@pytest.fixture(scope="module")
def mixed_contigs(tmp_path_factory):
    """A metagenome whose contigs qc does not all count: two high-quality
    MAGs, of contigs of 20 to 250 kb, a low-quality MAG, and an MGE in no
    MAG, with pairs between all of them. Its files, its pairs sorted and
    indexed, the directory of its tables (contig data, bin summary), the
    .metacool map that metacooler builds from its pairs, one bin per
    contig, and the contigs that qc counts: those over 100 kb of the
    high-quality MAGs."""
    import metacooler

    genomes = (
        synthetic.GenomeSpec("hq_1", gc=0.40, abundance=1.0, length=450_000),
        synthetic.GenomeSpec("hq_2", gc=0.60, abundance=2.0, length=400_000),
        synthetic.GenomeSpec("lq", gc=0.50, abundance=1.5, length=300_000),
    )
    mges = (synthetic.MgeSpec("plasmid", host="hq_1", gc=0.45, abundance=3.0, contig_lengths=(30_000, 20_000)),)
    community = synthetic.make_metagenome(
        seed=12, genomes=genomes, mges=mges, n_pairs=40_000, contig_median=100_000, contig_bounds=(20_000, 250_000),
        p_noise=0.05,
    )  # fmt: skip
    directory = tmp_path_factory.mktemp("mixed_contigs")
    files = community.write_all(directory)
    pairs = mio.sort_and_index_pairs(str(files["pairs"]), tmp_dir=str(directory))
    truth = community.truth
    bins = {"hq_1": "metator_00001_00000", "hq_2": "metator_00002_00000", "lq": "metator_00003_00000"}
    contig_data = pd.DataFrame(
        {
            "ID": range(1, len(truth) + 1),
            "Name": truth.Name,
            "Size": truth.Length,
            "Final_cluster": truth.Genome.map(bins).fillna("ND"),
        }
    )
    contig_data.to_csv(directory / "contig_data_final.txt", sep="\t", index=False)
    summary = pd.DataFrame(
        {"Weighted completeness": [0.95, 0.92, 0.5], "Weighted redundancy": [1.01, 1.05, 1.0]}, index=list(bins.values())
    )
    summary.to_csv(directory / "cluster_summary.txt", sep="\t")
    metacool = directory / "map.metacool"
    map_bins = pd.DataFrame({"chrom": truth.Name, "start": 0, "end": truth.Length})
    metacooler.from_pairs(str(metacool), pairs=str(files["pairs"]), bins=map_bins, contig_genome_map=None, ps=False)
    counted = set(truth.Name[(truth.Length > 100_000) & truth.Genome.isin(["hq_1", "hq_2"])])
    return community, files, pairs, directory, metacool, counted


def test_noise_from_a_map_counts_the_large_contigs_of_the_high_quality_mags(mixed_contigs, tmp_path, run_cli):
    """qc from a map counts the pairs between two contigs over 100 kb of the
    high-quality MAGs, and only them, as from the pairs files: the pairs
    within contigs, within MAGs and between MAGs that it logs, and the noise
    contact ratio and score, are those of the pairs files, and those of
    the simulated pairs between two such contigs. A pair with one read on
    such a contig, and the other on a smaller contig of a high-quality MAG,
    on a contig of a low-quality MAG, or on a contig in no MAG, is left out:
    read with them, the noise contact ratio is 46.50%, where it is 1.04%.
    The tests above have only such contigs, and read them so too."""
    community, files, indexed, directory, metacool, counted = mixed_contigs
    truth = community.truth.set_index("Name")
    pairs = community.pairs
    genome = pairs.chr1.map(truth.Genome), pairs.chr2.map(truth.Genome)
    both = pairs.chr1.isin(counted) & pairs.chr2.isin(counted)
    one = pairs.chr1.isin(counted) ^ pairs.chr2.isin(counted)
    other = pairs.chr2.where(pairs.chr1.isin(counted), pairs.chr1)[one]
    # The pairs that are left out: with each kind of other contig.
    small = truth.Genome.isin(["hq_1", "hq_2"]) & (truth.Length <= 100_000)
    large_lq = (truth.Genome == "lq") & (truth.Length > 100_000)
    assert all(other.isin(truth.index[kind]).any() for kind in [small, large_lq, truth.Kind == "mge"])
    within_contigs = int((both & (pairs.chr1 == pairs.chr2)).sum())
    within_mags, between_mags = int((both & (genome[0] == genome[1])).sum()), int((both & (genome[0] != genome[1])).sum())
    assert between_mags > 0 and len(counted) < len(truth)
    common = ["-a", files["assembly"], "-e", "DpnII", "-O", directory, "-T", tmp_path / "tmp"]
    run_cli("qc", *common, "-o", tmp_path / "pairs", "-t", "3,3", indexed)
    run_cli("qc", *common, "-o", tmp_path / "map", metacool)
    (pairs_log,) = (tmp_path / "pairs").glob("metator_qc_*.log")
    (map_log,) = (tmp_path / "map").glob("metator_qc_*.log")

    def messages(log):
        return [line.split(" :: ", 2)[2] for line in log.read_text().splitlines() if " :: INFO :: " in line]

    expected = [
        "Number of MAGs: 2",
        f"Contacts within contigs: {within_contigs}.",
        f"Contacts within MAGs: {within_mags}.",
        f"Contacts between MAGs: {between_mags}.",
    ]
    assert _metrics(messages(map_log))[:4] == expected
    assert _metrics(messages(map_log)) == _metrics(messages(pairs_log))
    assert f"{int(both.sum())} pairs have been extracted from {metacool}." in messages(map_log)
