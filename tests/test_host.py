# Tests of metator host: the MAGs that host the MGE clusters of a map. A small
# community whose contacts are counted by hand (the network of the raw counts),
# and the mock metagenome of tests/synthetic.py: a phage assembled in two
# contigs and a plasmid in one, each bound to a different host genome.

import json
import re
import shutil

import metacooler
import numpy as np
import pandas as pd
import pytest

import metator.clusters as mtcl
import metator.host as mth
import metator.metacool as mtmc
import metator.mge as mtm
import metator.network as mtn
from metator.version import __version__
from tests.hdf5 import walk
from tests.logs import logs

# The community: two MAGs of six contigs (any names: those of another binner
# here), a contig in none, and four MGE clusters. P1, a contig of the first
# MGE cluster, is in the first MAG too, as a validation clusters a plasmid
# contig with its host.
LENGTHS = {
    **{f"G1{x}": 10_000 for x in "abcdef"},
    **{f"G2{x}": 10_000 for x in "abcdef"},
    "U1": 5_000,
    "P1": 6_000,
    "P2": 4_000,
    "V1": 8_000,
    "F1": 5_000,
    "N1": 3_000,
}
FINAL = {**{f"G1{x}": "semibin_1" for x in "abcdef"}, **{f"G2{x}": "semibin_2" for x in "abcdef"}, "P1": "semibin_1"}
MGES = {"MetaTOR_MGE_00001": ["P1", "P2"], "MetaTOR_MGE_00002": ["V1"], "MetaTOR_MGE_00003": ["F1"], "MetaTOR_MGE_00004": ["N1"]}
# The contacts of the community, each pair of contigs once.
CONTACTS = [
    # The first MGE cluster: with itself (which does not count), with the
    # six contigs of semibin_1 (20), one of semibin_2 (1) and the contig in
    # no MAG (4).
    ("P1", "P2", 10),
    *[("P1", f"G1{x}", 3) for x in "abcdef"],
    ("P2", "G1a", 2),
    ("P2", "G2a", 1),
    ("P2", "U1", 4),
    # The second, with five contigs of each MAG (10 each), and with P1, an MGE
    # contig that semibin_1 holds too, which counts for no MAG (7).
    *[("V1", f"G1{x}", 2) for x in "abcde"],
    *[("V1", f"G2{x}", 2) for x in "abcde"],
    ("V1", "P1", 7),
    # The third, mostly with the contig in no MAG (20), and with five contigs of
    # semibin_2 (5): a host from a fifth of its contacts.
    ("F1", "U1", 20),
    *[("F1", f"G2{x}", 1) for x in "abcde"],
    # Within the MAGs, and within the fourth MGE cluster, which has no other
    # contact.
    ("G1a", "G1b", 5),
    ("G2a", "G2b", 5),
    ("N1", "N1", 3),
]

# What metator host finds there, with its defaults: one row per MGE cluster and
# MAG in contact, the MGE clusters by name, the MAGs of each by decreasing
# signal (ties in the order of the genomes table).
PAIRS = pd.DataFrame(
    [
        ("MetaTOR_MGE_00001", "semibin_1", 20.0, 6, 100 * 20 / 21, 100 * 20 / 32, True),
        ("MetaTOR_MGE_00001", "semibin_2", 1.0, 1, 100 * 1 / 21, 100 * 1 / 32, False),
        ("MetaTOR_MGE_00002", "semibin_1", 10.0, 5, 50.0, 100 * 10 / 27, True),
        ("MetaTOR_MGE_00002", "semibin_2", 10.0, 5, 50.0, 100 * 10 / 27, True),
        ("MetaTOR_MGE_00003", "semibin_2", 5.0, 5, 100.0, 20.0, True),
    ],
    columns=mth.PAIRS_COLUMNS,
)
CLUSTERS = pd.DataFrame(
    {
        "name": list(MGES),
        "n_hosts": [1, 2, 1, 0],
        "hosts": ["semibin_1", "semibin_1,semibin_2", "semibin_2", None],
        "host_shares": ["95.24", "50.00,50.00", "100.00", None],
    }
)
GENOMES = pd.DataFrame(
    {
        "name": ["semibin_1", "semibin_2"],
        "n_mges": [2, 2],
        "mges": ["MetaTOR_MGE_00001,MetaTOR_MGE_00002", "MetaTOR_MGE_00002,MetaTOR_MGE_00003"],
    }
)


def _community(directory, final=FINAL, mges=MGES, network=True):
    """The map of the community, made from its pairs as metacooler create pairs
    makes it (one bin per contig), its contigs in the genomes of final, with
    the network of the raw counts that metator network stores (normalization
    None) and the MGE clusters mges as the round mge of metator mge. Returns
    the map and the assembly."""
    rng = np.random.default_rng(3)
    assembly = directory / "community.fa"
    with open(assembly, "w") as fasta:
        for name, length in LENGTHS.items():
            seq = "".join(rng.choice(list("ACGT"), size=length))
            fasta.write(f">{name}\n" + "".join(seq[i : i + 60] + "\n" for i in range(0, length, 60)))
    order = {name: i for i, name in enumerate(LENGTHS)}
    pairs = directory / "community.pairs"
    with open(pairs, "w") as handle:
        handle.write("## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        read = 0
        for first, second, count in CONTACTS:
            if order[first] > order[second]:
                first, second = second, first
            for _ in range(count):
                read += 1
                shift = 100 if first == second else 0
                handle.write(f"r{read}\t{first}\t{LENGTHS[first] // 2}\t{second}\t{LENGTHS[second] // 2 + shift}\t+\t-\n")
    path = directory / "community.metacool"
    bins = pd.DataFrame({"chrom": list(LENGTHS), "start": 0, "end": list(LENGTHS.values())})
    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=final, ps=False, keep_partial=False)
    if network:
        contacts = mtn.contig_contacts(str(path))
        names = metacooler.MetaCooler(str(path)).contigs["name"]
        table = mtn.contig_data_table(mtn.create_contig_data(str(assembly)), names, mtn.contig_hits(contacts, len(names)))
        mtmc.write_map(str(path), str(path), network=mtmc.Network(table, "None"), force=True)
    if mges:
        summary = pd.DataFrame(
            {
                "BinName": list(mges),
                "BinLength": [sum(LENGTHS[contig] for contig in contigs) for contigs in mges.values()],
                "GC": 50.0,
                "Hit": 1,
                "BinningScore": [0.9 if len(contigs) > 1 else np.nan for contigs in mges.values()],
                "ContigsNumber": [len(contigs) for contigs in mges.values()],
                "Contigs": [",".join(contigs) for contigs in mges.values()],
            }
        )
        entry = {"name": mtcl.MGE, "tool": "metaTOR", "step": "mge"}
        mtmc.write_map(str(path), str(path), mge=(mtcl.mge_round(summary), entry), force=True)
    return path, assembly


@pytest.fixture(scope="module")
def community(tmp_path_factory):
    """The map of the community (see _community) and its assembly."""
    return _community(tmp_path_factory.mktemp("community"))


def _listed(values):
    """The values of a column, a missing one (None, or NaN in a column of text of
    pandas 3) as None."""
    return [None if pd.isna(value) else value for value in values]


def _found(metacool, **thresholds):
    """The hosts that find_hosts finds on the map metacool, from the network of
    its raw counts."""
    return mth.find_hosts(str(metacool), metacooler.MetaCooler(str(metacool)), mtn.map_network(str(metacool), None), **thresholds)


def test_find_hosts(community):
    """A MAG hosts an MGE cluster with at least 10% of its signal with the MAGs
    and 5 of its contigs in contact with it: one host, two, or none. The
    contacts of an MGE cluster within itself, with the contigs in no MAG and
    with the other MGE clusters count for no MAG, and an MGE contig that a MAG
    holds too counts for its MGE cluster only (V1-P1 is no contact with
    semibin_1); but those with the contigs in no MAG count in share_of_all
    (F1 has a host from a fifth of its contacts)."""
    hosts = _found(community[0])
    pd.testing.assert_frame_equal(hosts.pairs, PAIRS)
    pd.testing.assert_frame_equal(hosts.clusters, CLUSTERS)
    pd.testing.assert_frame_equal(hosts.genomes, GENOMES)


@pytest.mark.parametrize(
    "thresholds, hosts",
    [
        ({"min_share": 60}, ["semibin_1", None, "semibin_2", None]),
        ({"min_contigs": 6}, ["semibin_1", None, None, None]),
        ({"min_share": 50, "min_contigs": 5}, ["semibin_1", "semibin_1,semibin_2", "semibin_2", None]),
        ({"min_share": 0, "min_contigs": 1}, ["semibin_1,semibin_2", "semibin_1,semibin_2", "semibin_2", None]),
    ],
)
def test_the_thresholds(thresholds, hosts, community):
    """A share and a number of contigs at the thresholds make a host."""
    found = _found(community[0], **thresholds)
    assert _listed(found.clusters["hosts"]) == hosts
    minimum = {"min_share": mth.MIN_SHARE, "min_contigs": mth.MIN_CONTIGS, **thresholds}
    expected = (PAIRS["share"] >= minimum["min_share"]) & (PAIRS["contigs"] >= minimum["min_contigs"])
    assert list(found.pairs["host"]) == list(expected)


def test_mge_clusters_without_contact_with_the_mags(tmp_path, run_cli):
    """MGE clusters in contact with no MAG have no host, and no row in
    mge_hosts.txt, which holds its header only; the MAGs host none."""
    metacool, _ = _community(tmp_path, mges={"MetaTOR_MGE_00001": ["N1"]})
    hosts = _found(metacool)
    assert hosts.pairs.empty and list(hosts.pairs.columns) == mth.PAIRS_COLUMNS
    assert list(hosts.clusters["n_hosts"]) == [0] and _listed(hosts.clusters["hosts"]) == [None]
    assert list(hosts.genomes["n_mges"]) == [0, 0] and list(hosts.genomes["mges"]) == ["", ""]
    run_cli("host", metacool, "-o", tmp_path / "out", "--inplace")
    assert (tmp_path / "out" / "tables" / mth.HOSTS_FILE).read_text() == "\t".join(mth.PAIRS_COLUMNS) + "\n"
    assert _stored(metacool)[0]["n_hosts"].tolist() == [0]


def test_genomes_named_after_an_mge_cluster_are_no_mags(tmp_path):
    """The MGE clusters made genomes of the map (named after them, as metacooler
    reassign --round mge makes them) are no MAGs: no MGE cluster gets one as a
    host, and they get no column of the MAGs."""
    final = {**FINAL, "V1": "MetaTOR_MGE_00002", "F1": "MetaTOR_MGE_00003"}
    metacool, _ = _community(tmp_path, final=final)
    assert "MetaTOR_MGE_00002" in set(metacooler.MetaCooler(str(metacool)).genomes["name"])
    hosts = _found(metacool)
    pd.testing.assert_frame_equal(hosts.pairs, PAIRS)
    pd.testing.assert_frame_equal(hosts.genomes, GENOMES)


@pytest.mark.parametrize("case", ["no round mge", "no MAG", "only MGE genomes"])
def test_a_map_without_mge_clusters_or_mags_is_refused(case, tmp_path):
    """A map without MGE clusters (no round mge), or without MAGs (its genomes
    __unbinned__ and MGE clusters made genomes), is refused, with what to run."""
    final = {"no round mge": FINAL, "no MAG": None, "only MGE genomes": {"V1": "MetaTOR_MGE_00002"}}[case]
    metacool, _ = _community(tmp_path, final=final, mges=None if case == "no round mge" else MGES)
    expected = "run metator mge on it first" if case == "no round mge" else "holds no MAG to host its MGE clusters"
    with pytest.raises(mtmc.MetacoolError, match=expected):
        _found(metacool)


def test_write_hosts(tmp_path):
    """mge_hosts.txt: a row per MGE cluster and MAG in contact, with a header,
    tab-separated, every number as it is; a symbolic link there is replaced,
    not its target."""
    target = tmp_path / "target.txt"
    target.write_text("kept")
    path = tmp_path / mth.HOSTS_FILE
    path.symlink_to(target)
    mth.write_hosts(PAIRS, str(path))
    assert not path.is_symlink() and target.read_text() == "kept"
    lines = path.read_text().splitlines()
    assert lines[0] == "mge\tmag\tsignal\tcontigs\tshare\tshare_of_all\thost"
    assert lines[2] == f"MetaTOR_MGE_00001\tsemibin_2\t1.0\t1\t{100 / 21!r}\t3.125\tFalse"
    read = pd.read_csv(path, sep="\t", **mtcl.EXACT_FLOATS)
    pd.testing.assert_frame_equal(read, PAIRS)


def test_plot_hosts(tmp_path, saved_figures):
    """mge_hosts.pdf: the signal of the hosts, on a log scale, and their share,
    from 0 to 100%: one count per host."""
    mth.plot_hosts(PAIRS, str(tmp_path / mth.PLOT_FILE))
    assert (tmp_path / mth.PLOT_FILE).stat().st_size > 0
    figure = saved_figures[mth.PLOT_FILE]
    signal, share = figure.axes
    assert figure._suptitle.get_text() == "The 4 hosts of 3 MGE clusters"
    assert signal.get_xscale() == "log" and share.get_xlim() == (0, 100)
    assert signal.get_xlabel() == "Signal of the MGE cluster with its host (log scale)"
    assert share.get_xlabel() == "Share of the signal of the MGE cluster with the MAGs (%)"
    assert [ax.get_ylabel() for ax in figure.axes] == ["Hosts", "Hosts"]
    assert sum(patch.get_height() for patch in signal.patches) == 4
    heights = {round(patch.get_x()): patch.get_height() for patch in share.patches if patch.get_height()}
    assert heights == {50: 2, 95: 2}


def test_plot_hosts_without_host(tmp_path, saved_figures):
    """Without a host, the axes say so."""
    mth.plot_hosts(PAIRS.assign(host=False), str(tmp_path / mth.PLOT_FILE))
    figure = saved_figures[mth.PLOT_FILE]
    assert figure._suptitle.get_text() == "The 0 hosts of 0 MGE clusters"
    assert [[text.get_text() for text in ax.texts] for ax in figure.axes] == [["No host"], ["No host"]]
    assert not any(ax.patches for ax in figure.axes)


# metator host: the hosts stored in the map of the run.


@pytest.fixture
def host_run(community, tmp_path):
    """A private copy of the map of the community, the output directory of a
    host run on it (tmp_path/out), and the assembly of the community."""
    metacool = tmp_path / "community.metacool"
    shutil.copyfile(community[0], metacool)
    return metacool, tmp_path / "out", community[1]


def _stored(metacool):
    """What metator host stored in the map: the host columns of the clusters of
    the round mge (in its order), those of the genomes that are MAGs (in the
    order of the genomes table), and the host entry."""
    mc = metacooler.MetaCooler(str(metacool))
    clusters = mc.clusters[(mc.clusters["round"] == mtcl.MGE).to_numpy()].reset_index(drop=True)
    genomes = mc.genomes[(mc.genomes["name"] != "__unbinned__").to_numpy()].reset_index(drop=True)
    return (
        clusters[["name", *mtmc.HOST_CLUSTER_COLUMNS]],
        genomes[["name", *mtmc.HOST_GENOME_COLUMNS]],
        mc.metadata["metator"]["host"],
    )


def _outside_hosts(path):
    """Every group, dataset and attribute of a map, as stored, but what metator
    host stores: the host columns of its clusters and genomes tables, the host
    entry of its metator metadata, and metacooler's record of the annotations
    that add or drop the host columns of the genomes table."""
    host = {f"/clusters/table/{column}" for column in mtmc.HOST_CLUSTER_COLUMNS}
    host |= {f"/genomes/{column}" for column in mtmc.HOST_GENOME_COLUMNS}
    annotated = {f"genomes.{column}" for column in mtmc.HOST_GENOME_COLUMNS}

    def of_hosts(entry):
        columns = {*entry.get("added", []), *entry.get("dropped", [])}
        return bool(columns) and columns <= annotated

    found = {name: value for name, value in walk(path).items() if name not in host}
    for name in ["/clusters/table", "/genomes"]:
        kind, keys, attributes = found[name]
        found[name] = (kind, [key for key in keys if f"{name}/{key}" not in host], attributes)
    kind, keys, attributes = found["/"]
    kept = []
    for key, dtype, value in attributes:
        if key == "metadata":
            metadata = json.loads(value)
            metadata["metator"].pop("host", None)
            metadata["annotate"] = [entry for entry in metadata.get("annotate", []) if not of_hosts(entry)]
            value = json.dumps(metadata, sort_keys=True)
        kept.append((key, dtype, value))
    found["/"] = (kind, keys, kept)
    return found


@pytest.mark.parametrize("output", ["the output directory", "--inplace", "-F the output directory", "an entry"])
def test_host_stores_its_hosts_in_the_map(output, host_run, tmp_path, run_cli):
    """host MAP stores the hosts of the MGE clusters in the map of its output
    directory (MAP itself, with --inplace; for an entry of a collection, a map
    of its own): the hosts of each MGE cluster in the round mge, the MGE
    clusters of each MAG in the genomes table, and the run under the host
    entry of the metator metadata of the map; everything else is kept as it
    was. It writes tables/mge_hosts.txt and plots/mge_hosts.pdf, and its log,
    in its output directory."""
    metacool, out, _ = host_run
    before = _outside_hosts(metacool)
    rounds = [entry for entry in metacooler.MetaCooler(str(metacool)).cluster_rounds]
    source = metacool
    if output == "an entry":
        source = f"{tmp_path / 'maps.mmetacool'}::entries/lib1"
        metacooler.create_mmetacool(tmp_path / "maps.mmetacool", {"lib1": metacool})
    options = ["--inplace"] if output == "--inplace" else []
    written = metacool if output == "--inplace" else out / ("lib1.metacool" if output == "an entry" else metacool.name)
    if output.startswith("-F"):
        out.mkdir()
        written.write_text("an older map")
        options.append("-F")
    unchanged = metacool.read_bytes()
    run_cli("host", source, "-o", out, *options)
    metacooler.validate(str(written), full=True)
    clusters, genomes, entry = _stored(written)
    expected = CLUSTERS.assign(n_hosts=CLUSTERS["n_hosts"].astype(float), hosts=["semibin_1", "semibin_1,semibin_2",
               "semibin_2", ""], host_shares=["95.24", "50.00,50.00", "100.00", ""])  # fmt: skip
    pd.testing.assert_frame_equal(clusters, expected)
    pd.testing.assert_frame_equal(genomes, GENOMES.assign(n_mges=GENOMES["n_mges"].astype(float)))
    assert {key: entry[key] for key in ["version", "normalization"]} == {"version": __version__, "normalization": "None"}
    assert entry["command"] == f"metator host {' '.join([str(source), '-o', str(out), *options])}"
    assert entry["parameters"]["min_share"] == 10.0 and entry["parameters"]["min_contigs"] == 5
    read = {"contigs": ["GC_content", "Hit"], "bins": [], "rounds": ["mge"], "genomes": ["name"]}
    assert entry["read"] == read
    assert entry["network_digest"] == mtn.stored_network(str(written)).digest
    # The rounds as they were, and everything outside what host stores.
    assert metacooler.MetaCooler(str(written)).cluster_rounds == rounds
    assert _outside_hosts(written) == before
    if output != "--inplace":
        assert metacool.read_bytes() == unchanged
    pairs = pd.read_csv(out / "tables" / mth.HOSTS_FILE, sep="\t", **mtcl.EXACT_FLOATS)
    pd.testing.assert_frame_equal(pairs, PAIRS)
    assert (out / "plots" / mth.PLOT_FILE).read_bytes().startswith(b"%PDF")
    (log,) = logs(out, "host")
    logged = log.read_text()
    assert "4 MGE clusters and 2 MAGs: 2 MGE clusters with one host, 1 with several, 1 with none" in logged
    ended = rf"INFO :: Contact map with the hosts of the MGE clusters written after 0h 0m [0-9]+s: {re.escape(str(written))}\n"
    assert re.search(ended, logged)
    assert not list(tmp_path.glob("*.lock")) and not list(out.glob("*.lock"))


def test_host_rerun_replaces_its_hosts(host_run, run_cli, logged):
    """A second run replaces the hosts of the first, with its thresholds: in
    place, with a warning before any work."""
    metacool, out, _ = host_run
    run_cli("host", metacool, "--inplace", "-o", out)
    run_cli("host", metacool, "--inplace", "-o", out, "--min-contigs", "6")
    warning = (
        f"{metacool} holds the hosts of its MGE clusters (metator host), which this run replaces, in place: leave out "
        f"--inplace to write the map of this run to another file, and leave {metacool} as it is."
    )
    assert ("WARNING", warning) in logged()
    clusters, genomes, entry = _stored(metacool)
    assert list(clusters["hosts"]) == ["semibin_1", "", "", ""] and list(clusters["n_hosts"]) == [1, 0, 0, 0]
    assert list(genomes["n_mges"]) == [1, 0] and list(genomes["mges"]) == ["MetaTOR_MGE_00001", ""]
    assert entry["parameters"]["min_contigs"] == 6
    pairs = pd.read_csv(out / "tables" / mth.HOSTS_FILE, sep="\t")
    assert list(pairs["host"]) == [True, False, False, False, False]


@pytest.mark.parametrize("case", ["no network", "no round mge", "no MAG"])
def test_host_refuses_a_map(case, tmp_path, run_cli, logged):
    """A map without the network of metator network, without MGE clusters or
    without MAGs is refused before anything is written, with what to run."""
    final = None if case == "no MAG" else FINAL
    metacool, _ = _community(tmp_path, final=final, mges=None if case == "no round mge" else MGES, network=case != "no network")
    before = metacool.read_bytes()
    with pytest.raises(mtmc.MetacoolError):
        run_cli("host", metacool, "-o", tmp_path / "out")
    expected = {
        "no network": "holds no network of its contigs: run metator network on it first",
        "no round mge": "run metator mge on it first",
        "no MAG": "holds no MAG to host its MGE clusters",
    }[case]
    assert any(level == "ERROR" and expected in message for level, message in logged())
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in (tmp_path / "out").iterdir())[0].startswith("metator_host_")
    assert [path.name for path in (tmp_path / "out").iterdir() if not path.name.endswith(".log")] == []


@pytest.mark.parametrize("case", ["mge_hosts.txt", "mge_hosts.png", "tables, a file"])
@pytest.mark.parametrize("force", [False, True], ids=["refused", "-F"])
def test_host_refuses_an_earlier_version_s_outputs(case, force, host_run, tmp_path, run_cli):
    """The table or the plot of the hosts that an earlier version wrote in the
    output directory (mge_hosts.txt, mge_hosts.png, which tables/mge_hosts.txt
    and plots/mge_hosts.pdf replace) would pass for this run's: refused with or
    without -F, before anything is written, and left as they are; and so is a
    tables that is no directory, where the run writes its table."""
    metacool, out, _ = host_run
    out.mkdir()
    blocking = out / ("tables" if case == "tables, a file" else case)
    blocking.write_text("an earlier run's\n")
    before = metacool.read_bytes()
    with pytest.raises(ValueError):
        run_cli("host", metacool, "-o", out, "--inplace", *(["-F"] if force else []))
    assert metacool.read_bytes() == before and blocking.read_text() == "an earlier run's\n"
    assert not (out / "plots").exists() and (case == "tables, a file" or not (out / "tables").exists())
    (log,) = logs(out, "host")
    message = (
        f"{blocking} is not a directory, where the run writes the table of the hosts"
        if case == "tables, a file"
        else f"{out} holds results of an earlier run, which this run does not write and does not replace: {case!r}"
    )
    assert f"ERROR :: {message}" in log.read_text()


def _dropped(where, why, outfile):
    """The warning about the hosts that a map no longer holds (where: the map, or
    the file written and the map), which described what the run replaced
    (why)."""
    return (
        "WARNING",
        f"{where} the hosts of the MGE clusters that metator host found (the host entry of its metator metadata, the "
        f"columns n_hosts, hosts, host_shares of its round mge and n_mges, mges of its genomes table), which described "
        f"{why}, which this run replaced. Run metator host again on {outfile} to find them.",
    )


def _hosts_held(metacool):
    """Whether the map holds anything that metator host stores."""
    mc = metacooler.MetaCooler(str(metacool))
    columns = [*mc.clusters.columns, *mc.genomes.columns]
    stored = [column for column in [*mtmc.HOST_CLUSTER_COLUMNS, *mtmc.HOST_GENOME_COLUMNS] if column in columns]
    return "host" in mc.metadata.get("metator", {}), stored


@pytest.mark.parametrize(
    "writer, why",
    [
        ("mge", "its round mge"),
        ("mge --random", None),
        ("mge --reassign", "its round mge"),
        ("partition", "its clustering rounds"),
        ("validation", "its genomes"),
        ("network", None),
        ("network -n length", "its network"),
    ],
)
def test_runs_that_replace_what_the_hosts_describe_drop_them(writer, why, host_run, run_cli, logged, debug_log):
    """The hosts describe the round mge, the genomes and the network of the map:
    a run that replaces one of them drops them all, once the map is written,
    with a warning (mge, but with --random; partition, which replaces the
    rounds; validation, the genomes; network, with another network). A run
    that leaves them as they were keeps them."""
    metacool, out, assembly = host_run
    run_cli("host", metacool, "--inplace", "-o", out)
    assert _hosts_held(metacool) == (True, [*mtmc.HOST_CLUSTER_COLUMNS, *mtmc.HOST_GENOME_COLUMNS])
    if writer.startswith("mge"):
        mges = metacool.parent / "mges.txt"
        mges.write_text("".join(f"{contig}\n" for contigs in MGES.values() for contig in contigs))
        run_cli("mge", metacool, "-a", assembly, "-m", mges, "--inplace", "-T", metacool.parent / "tmp",
                *writer.split()[1:])  # fmt: skip
    elif writer == "partition":
        mtmc.write_map(str(metacool), str(metacool), rounds=mtcl.Rounds().tables(), round_names=[mtcl.OVERLAPPING], force=True)
    elif writer == "validation":
        clusters = mtmc.FinalClusters({"G1a": "metator_00001_00000"}, None, None, {})
        mtmc.write_map(str(metacool), str(metacool), clusters=clusters, force=True)
    else:
        run_cli(
            "network", metacool, "-a", assembly, "--inplace", "-o", out, "-n", writer.split()[-1] if " " in writer else "None"
        )
    held = _hosts_held(metacool)
    warning = _dropped(f"{metacool} no longer holds", why, metacool)
    if why is None:
        assert held == (True, [*mtmc.HOST_CLUSTER_COLUMNS, *mtmc.HOST_GENOME_COLUMNS])
        assert not [message for level, message in logged() if level == "WARNING" and "hosts of the MGE clusters" in message]
    else:
        assert held == (False, [])
        assert warning in logged()
        messages = [message for _, message in logged()]
        assert messages.index(warning[1]) > max(i for i, message in enumerate(messages) if " is written: " in message)


# The mock metagenome: metator mge, then metator host, find each MGE its host
# genome.


@pytest.fixture(scope="module")
def mock_map(metagenome, metagenome_files, tmp_path_factory):
    """The map of the clusters of the mock metagenome, as metator validation
    writes it: one cluster per genome, the plasmid contig in the cluster of its
    host (a validation clusters it with its host), the phage contigs in none;
    with the network that metator network stores (empirical_hit). Returns the
    map, the MGE contigs and the cluster of each genome."""
    directory = tmp_path_factory.mktemp("mock_map")
    truth = metagenome.truth
    clusters = {genome.name: f"metator_{k:05d}_00000" for k, genome in enumerate(metagenome.genomes, start=1)}
    plasmid = metagenome.molecule("plasmid_1")
    final = {name: clusters.get(genome, "ND") for name, genome in zip(truth.Name, truth.Genome, strict=True)}
    final.update({contig.name: clusters[plasmid.host] for contig in plasmid.contigs})
    path = directory / "clusters.metacool"
    bins = pd.DataFrame({"chrom": truth.Name, "start": 0, "end": truth.Length})
    metacooler.from_pairs(str(path), pairs=str(metagenome_files["pairs"]), bins=bins, contig_genome_map=final, ps=False,
                          keep_partial=False)  # fmt: skip
    contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
    contacts = mtn.contig_contacts(str(path))
    names = metacooler.MetaCooler(str(path)).contigs["name"]
    table = mtn.contig_data_table(contig_data, names, mtn.contig_hits(contacts, len(names)))
    mtmc.write_map(str(path), str(path), network=mtmc.Network(table, "empirical_hit"), force=True)
    mges = directory / "mges.txt"
    mges.write_text("".join(f"{name}\n" for name in truth.Name[truth.Kind == "mge"]))
    return path, mges, clusters


def test_mge_then_host_find_each_mge_its_host(metagenome, mock_map, metagenome_inputs, tmp_path, run_cli):
    """metator mge clusters the phage and the plasmid, and metator host gives each
    the cluster of its host genome, and no other; the plasmid too, though the
    cluster of its host holds its contig."""
    source, mges, clusters = mock_map
    metacool = tmp_path / "clusters.metacool"
    shutil.copyfile(source, metacool)
    run_cli("mge", metacool, "-a", metagenome_inputs["assembly"], "-m", mges, "--inplace", "-T", tmp_path / "tmp")
    run_cli("host", metacool, "--inplace", "-o", tmp_path / "out")
    mc = metacooler.MetaCooler(str(metacool))
    members = mc.cluster_members(mtcl.MGE)
    cluster_of = dict(zip(members["contig"], members["cluster"], strict=True))
    clusters_found, genomes, _ = _stored(metacool)
    hosts = dict(zip(clusters_found["name"], clusters_found["hosts"], strict=True))
    for mge in metagenome.mges:
        (cluster,) = {cluster_of[contig.name] for contig in mge.contigs}
        assert hosts[cluster] == clusters[mge.host]
    pairs = pd.read_csv(tmp_path / "out" / "tables" / mth.HOSTS_FILE, sep="\t")
    assert (pairs.loc[pairs["host"], "share"] > 80).all()
    assert sorted(genomes.loc[genomes["n_mges"] > 0, "name"]) == sorted(clusters[mge.host] for mge in metagenome.mges)


def test_mge_genomes_change_no_host(mock_map, metagenome_inputs, tmp_path, run_cli):
    """metator mge --reassign makes the MGE clusters genomes of the map: metator
    host takes them for no MAG, gives them no column of the MAGs, and finds
    the hosts that it finds without --reassign (the contigs of an MGE cluster
    count for it only either way)."""
    source, mges, _ = mock_map
    found = {}
    for reassign in (False, True):
        directory = tmp_path / ("reassign" if reassign else "plain")
        directory.mkdir()
        metacool = directory / "clusters.metacool"
        shutil.copyfile(source, metacool)
        run_cli("mge", metacool, "-a", metagenome_inputs["assembly"], "-m", mges, "--inplace", "-T", directory / "tmp",
                "-o", directory / "out", *(["--reassign"] if reassign else []))  # fmt: skip
        run_cli("host", metacool, "--inplace", "-o", directory / "out")
        found[reassign] = pd.read_csv(directory / "out" / "tables" / mth.HOSTS_FILE, sep="\t")
        if reassign:
            mc = metacooler.MetaCooler(str(metacool))
            genomes = mc.genomes
            mge_genomes = (genomes["mge"] == "True").to_numpy()
            assert sorted(genomes.loc[mge_genomes, "name"]) == sorted(set(mc.cluster_members(mtcl.MGE)["cluster"]))
            assert genomes.loc[mge_genomes, "n_mges"].isna().all()
            assert not set(found[reassign]["mag"]) & set(genomes.loc[mge_genomes, "name"])
    assert len(found[True]) and found[True]["host"].any()
    pd.testing.assert_frame_equal(found[True], found[False], check_exact=False, rtol=1e-12)


def test_host_of_the_mge_contigs_of_a_validation(mock_map, metagenome_inputs, tmp_path):
    """The MGE cluster of the plasmid, whose contig the cluster of its host holds,
    has a host: its contacts with its own contig are no contacts with the
    cluster (the cluster is taken without it)."""
    source, mges, clusters = mock_map
    metacool = tmp_path / "clusters.metacool"
    shutil.copyfile(source, metacool)
    contigs = mges.read_text().split()
    summary, _ = mtm.mge_binning(str(metacool), contigs)
    mtmc.write_map(str(metacool), str(metacool), mge=(mtcl.mge_round(summary), {"name": mtcl.MGE}), force=True)
    mc = metacooler.MetaCooler(str(metacool))
    mag, members, _, mags = mth.map_bins(str(metacool), mc)
    positions = members["position"].to_numpy()
    assert (mag[positions] == "").all()
    assert mags == list(clusters.values())


@pytest.mark.parametrize(
    "parameters, expected",
    [
        ({"min_share": 20.0, "min_contigs": 3}, (20.0, 3)),
        ({"mge_min_share": 15.0, "mge_min_contigs": 2}, (15.0, 2)),
        ({}, (mth.MIN_SHARE, mth.MIN_CONTIGS)),
        (None, (mth.MIN_SHARE, mth.MIN_CONTIGS)),
    ],
    ids=["host", "pipeline --mges", "none recorded", "no host entry"],
)
def test_recorded_thresholds(parameters, expected):
    """The thresholds with which the hosts that a map holds were found, as the
    host entry of its metator metadata records the parameters of the run:
    those of metator host, or of pipeline --mges (--mge-min-share,
    --mge-min-contigs); the defaults of metator host for those it does not
    record, or for a map without hosts."""
    import types

    metator = {} if parameters is None else {"host": {"command": "metator host", "parameters": parameters}}
    mc = types.SimpleNamespace(metadata={"metator": metator})
    assert mth.recorded_thresholds(mc) == expected
