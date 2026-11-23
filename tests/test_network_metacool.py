# Tests of the network built from a .metacool contact map
# (metator.network.contig_contacts and map_network), against the network of
# the pairs the map is made of, counted straight from the pairs file
# (pairs_network), on the mock metagenome of tests/synthetic.py and on
# tests_data.

import collections
import gzip
import shutil
from pathlib import Path

import cooler
import h5py
import metacooler
import numpy as np
import pandas as pd
import pytest

import metator.metacool as mtmc
import metator.network as mtn
from tests.maps import store_assembly_columns

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"

# Pairs files of tests_data, with their assembly: alignment.pairs holds
# reads at position 0; subset.pairs.gz none, nor on the last base of a
# contig.
TESTS_DATA_PAIRS = {
    "alignment": (TESTS_DATA / "alignment.pairs", TESTS_DATA / "assembly.fa"),
    "subset": (TESTS_DATA / "subset.pairs.gz", TESTS_DATA / "subset.assembly.fa"),
}

# One bin per contig, or bins of 1 kb.
BIN_SIZES = [None, 1000]

PAIRS_HEADER = "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"


def make_metacool(path, pairs, assembly, binsize=None, filters=None):
    """A .metacool file of the pairs, its contigs in the order of the
    assembly, all of them unbinned: one bin per contig without binsize.
    With filters (metacooler's --filter), the contigs that they drop are not
    in the map, nor are the pairs with a read on them.

    metacooler indexes the assembly (.fai): it reads a private copy.
    """
    path = Path(path)
    assembly = shutil.copyfile(assembly, path.with_suffix(".fa"))
    if binsize:
        bins = {"binsize": binsize}
    else:
        contig_data = mtn.create_contig_data(str(assembly))
        lengths = [data["length"] for data in contig_data.values()]
        bins = {"bins": pd.DataFrame({"chrom": list(contig_data), "start": 0, "end": lengths})}
    metacooler.from_pairs(
        str(path),
        pairs=str(pairs),
        fasta=str(assembly),
        contig_genome_map=None,
        contig_ordering="fasta",
        ps=False,
        filters=filters,
        **bins,
    )
    return path


def filtered_metacool(path, pairs, assembly, binsize, shortest, how):
    """A map of the pairs that lacks the contigs shorter than shortest, and
    the pairs with a read on them, as metacooler makes it: with a contig
    filter at creation (create pairs --filter), or cut from the map of all
    the contigs (extract --filter, which writes a .cool file, made a map
    again by create cool)."""
    filters = [f"contigs.length >= {shortest}"]
    if how == "create":
        return make_metacool(path, pairs, assembly, binsize, filters)
    full = make_metacool(path.with_name("full.metacool"), pairs, assembly, binsize)
    cut = path.with_suffix(".cool")
    metacooler.extract(str(full), filters=filters, output=str(cut))
    metacooler.from_cool(str(cut), str(path), contig_genome_map=None, contig_ordering="cool")
    return path


def kept_pairs(pairs_file, contigs, path):
    """The pairs file (header included) of the pairs whose two reads lie on
    the given contigs."""
    with (gzip.open if str(pairs_file).endswith(".gz") else open)(pairs_file, "rt") as pairs, open(path, "w") as kept:
        for line in pairs:
            fields = line.split("\t")
            if line.startswith("#") or (fields[1] in contigs and fields[3] in contigs):
                kept.write(line)
    return path


def read_pairs(path):
    """The contig names and positions of the pairs of a pairs file."""
    with (gzip.open if str(path).endswith(".gz") else open)(path, "rt") as pairs:
        for line in pairs:
            if not line.startswith("#"):
                fields = line.split("\t")
                yield fields[1], int(fields[2]), fields[3], int(fields[4])


def outside_contacts(pairs_file, contig_data):
    """Inter-contig pairs with a read on the last base of its contig, and
    those with a read at position 0, by contig pair (ids): pairs on the
    edges of the contigs, which the network counts. Counted from the pairs
    file."""
    last_base, position_0 = collections.Counter(), collections.Counter()
    ids = {name: i for i, name in enumerate(contig_data, start=1)}
    for contig1, pos1, contig2, pos2 in read_pairs(pairs_file):
        id1, id2 = ids[contig1], ids[contig2]
        if id1 == id2:
            continue
        edge = (min(id1, id2), max(id1, id2))
        if pos1 == 0 or pos2 == 0:
            position_0[edge] += 1
        elif pos1 == contig_data[contig1]["length"] or pos2 == contig_data[contig2]["length"]:
            last_base[edge] += 1
    return last_base, position_0


def pairs_network(data):
    """The contacts and the hits of the pairs themselves, by contig name: the
    oracle of the network of a map of the pairs. Counted in Python, from the
    pairs file: each pair whose two reads lie on contigs of the assembly,
    within their length (at position 0 too, where metacooler puts the read
    on the first base), counts once, and once in the hits of each of its two
    contigs (a pair within a contig, twice in its hits). Returns the
    contacts, by pair of contig names (sorted, a contig with itself
    included), and the hits, by contig name."""
    contig_data = mtn.create_contig_data(str(data["assembly"]), data["depth"], data["enzyme"])
    contacts, hits = collections.Counter(), collections.Counter()
    for contig1, pos1, contig2, pos2 in read_pairs(data["pairs"]):
        if contig1 not in contig_data or contig2 not in contig_data:
            continue
        if pos1 > contig_data[contig1]["length"] or pos2 > contig_data[contig2]["length"]:
            continue
        contacts[tuple(sorted((contig1, contig2)))] += 1
        hits[contig1] += 1
        hits[contig2] += 1
    return contacts, hits


def map_contacts(metacool, **kwargs):
    """The contacts and the hits of a map (metator.network.contig_contacts
    and contig_hits), by contig names as pairs_network gives them, with the
    names of its contigs, in its order."""
    names = [str(name) for name in metacooler.MetaCooler(str(metacool)).contigs["name"]]
    table = mtn.contig_contacts(str(metacool), **kwargs)
    contacts = {
        tuple(sorted((names[i], names[j]))): count for i, j, count in zip(table["id1"], table["id2"], table["count"], strict=True)
    }
    hits = dict(zip(names, mtn.contig_hits(table, len(names)), strict=True))
    return contacts, hits, names


def edges_by_name(network, names):
    """The edges of a network (metator.network.map_network), by pair of contig
    names, sorted: their weights."""
    return {
        tuple(sorted((names[i - 1], names[j - 1]))): weight
        for i, j, weight in zip(network["contig1"], network["contig2"], network["weight"], strict=True)
    }


@pytest.fixture(params=["synthetic", "alignment"])
def dataset(request, metagenome_files):
    """Pairs, assembly, depth file and enzyme of a data set. The pairs of
    subset.pairs.gz, the largest, go through the tests by chunk."""
    if request.param == "synthetic":
        files = metagenome_files
        return {
            "name": "synthetic",
            "pairs": files["pairs"],
            "assembly": files["assembly"],
            "depth": str(files["depth"]),
            "enzyme": "DpnII",
        }
    pairs, assembly = TESTS_DATA_PAIRS[request.param]
    return {"name": request.param, "pairs": pairs, "assembly": assembly, "depth": None, "enzyme": None}


@pytest.mark.parametrize("binsize", BIN_SIZES)
def test_network_of_a_metacool(dataset, binsize, tmp_path):
    """The contacts and the hits of a metacool are those of its pairs: the
    pairs with a read on the last base of a contig (the mock metagenome) or
    at position 0 (alignment.pairs) count. A metacool sums its libraries:
    no network per library. Its network is its contacts between two
    contigs, raw without weights."""
    metacool = make_metacool(tmp_path / "map.metacool", dataset["pairs"], dataset["assembly"], binsize)
    expected_contacts, expected_hits = pairs_network(dataset)
    contacts, hits, names = map_contacts(metacool)
    assert contacts == expected_contacts
    assert hits == {name: expected_hits[name] for name in names}
    network = mtn.map_network(str(metacool), weights=None)
    assert edges_by_name(network, names) == {pair: n for pair, n in expected_contacts.items() if pair[0] != pair[1]}
    assert (network["contig1"] < network["contig2"]).all()
    # What each data set holds: the contacts on the edges of the contigs.
    contig_data = mtn.create_contig_data(str(dataset["assembly"]))
    last_base, position_0 = outside_contacts(dataset["pairs"], contig_data)
    assert (bool(last_base), bool(position_0)) == ((True, False) if dataset["name"] == "synthetic" else (False, True))


def test_normalized_network_of_a_metacool(dataset, tmp_path):
    """Normalized (empirical_hit, the default): each edge is its count times
    the weights of its two contigs, from their hits, whatever the bin size:
    the network of the pairs normalized so."""
    expected_contacts, expected_hits = pairs_network(dataset)
    networks = []
    for binsize in BIN_SIZES:
        metacool = make_metacool(tmp_path / f"map_{binsize}.metacool", dataset["pairs"], dataset["assembly"], binsize)
        _, hits, names = map_contacts(metacool)
        weights = mtn.normalization_weights(None, np.array([hits[name] for name in names]), None, None, "empirical_hit")
        network = edges_by_name(mtn.map_network(str(metacool), weights=weights), names)
        expected = {
            (first, second): n / np.sqrt(expected_hits[first]) / np.sqrt(expected_hits[second])
            for (first, second), n in expected_contacts.items()
            if first != second
        }
        assert network.keys() == expected.keys()
        for pair, weight in network.items():
            assert weight == pytest.approx(expected[pair], rel=1e-12, abs=0)
        networks.append(network)
    assert networks[0] == networks[1]


@pytest.mark.parametrize("binsize", BIN_SIZES)
@pytest.mark.parametrize(
    "dataset, how", [("synthetic", "create"), ("synthetic", "extract"), ("alignment", "create")], indirect=["dataset"]
)
def test_network_of_a_filtered_metacool(dataset, how, binsize, tmp_path):
    """A map that lacks contigs of the assembly, made with a contig filter,
    which leaves out the pairs with a read on the contigs it drops: its
    contacts and hits are those of the pairs it holds, whose two reads lie
    on the contigs it keeps. The contigs it lacks are no node of its
    network, and no row of the contig data of a run (see
    test_contig_data_table_of_a_map_lacking_contigs)."""
    contig_data = mtn.create_contig_data(str(dataset["assembly"]))
    lengths = {name: data["length"] for name, data in contig_data.items()}
    shortest = int(np.median(list(lengths.values())))
    kept = {name for name, length in lengths.items() if length >= shortest}
    metacool = filtered_metacool(tmp_path / "map.metacool", dataset["pairs"], dataset["assembly"], binsize, shortest, how)
    assert sorted(metacooler.MetaCooler(str(metacool)).contigs.name) == sorted(kept)
    assert 0 < len(kept) < len(lengths)
    pairs = {**dataset, "pairs": kept_pairs(dataset["pairs"], kept, tmp_path / "kept.pairs")}

    expected_contacts, expected_hits = pairs_network(pairs)
    contacts, hits, names = map_contacts(metacool)
    assert contacts == expected_contacts
    assert hits == {name: expected_hits[name] for name in names} and sum(hits.values()) > 0
    assert set(names) == kept


def test_contacts_by_chunk(tmp_path):
    """The contacts do not depend on the number of pixels read at a time (a
    pixel at a time with 1), and nor does the network."""
    pairs, assembly = TESTS_DATA_PAIRS["alignment"]
    metacool = make_metacool(tmp_path / "map.metacool", pairs, assembly)
    contacts = mtn.contig_contacts(str(metacool))
    for chunksize in [1, 997]:
        pd.testing.assert_frame_equal(mtn.contig_contacts(str(metacool), chunksize=chunksize), contacts)
    # Sorted, one row per contig pair, the upper triangle.
    assert (contacts.id1 <= contacts.id2).all()
    assert not contacts.duplicated(["id1", "id2"]).any()
    assert contacts.equals(contacts.sort_values(["id1", "id2"]))
    # Every pair counts twice in the hits.
    n_contigs = len(metacooler.MetaCooler(str(metacool)).contigs)
    hits = mtn.contig_hits(contacts, n_contigs)
    assert len(hits) == n_contigs and hits.sum() == 2 * contacts["count"].sum()


def test_network_by_chunk(tmp_path):
    """Nor does the network of a map of 1 kb bins, summed by contig pair in
    chunks of pixels. An equivalence guard."""
    pairs, assembly = TESTS_DATA_PAIRS["subset"]
    metacool = make_metacool(tmp_path / "map.metacool", pairs, assembly, 1000)
    tables = [mtn.contig_contacts(str(metacool), chunksize=chunksize) for chunksize in [mtn.METACOOL_CHUNKSIZE, 9973]]
    pd.testing.assert_frame_equal(tables[0], tables[1])
    networks = [mtn.map_network(str(metacool), weights=None, contacts=table) for table in tables]
    pd.testing.assert_frame_equal(networks[0], networks[1])
    assert len(networks[0]) > 2 * 9973


def test_network_of_an_mmetacool_entry(tmp_path):
    """An entry of a .mmetacool collection, given by its URI."""
    pairs, assembly = TESTS_DATA_PAIRS["subset"]
    metacool = make_metacool(tmp_path / "map.metacool", pairs, assembly, 1000)
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    networks = [mtn.map_network(uri, weights=None) for uri in [str(metacool), f"{collection}::entries/lib1"]]
    pd.testing.assert_frame_equal(networks[0], networks[1])
    assert len(networks[0])


# Maps that cannot give the network of an assembly.


def _write_fasta(path, contigs):
    path.write_text("".join(f">{name}\n{'ACGT' * (length // 4)}\n" for name, length in contigs.items()))
    return path


def _small_map(tmp_path, contigs, pairs, ordering="fasta"):
    """A one-bin-per-contig map of the given contigs and pairs, its contigs in
    their order (ordering fasta), or longest first (genome-then-length, the
    order of metacooler create pairs and hicstuff)."""
    fasta = _write_fasta(tmp_path / "contigs.fa", contigs)
    pairs_file = tmp_path / "map.pairs"
    pairs_file.write_text(
        PAIRS_HEADER + "".join(f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t-\n" for i, (c1, p1, c2, p2) in enumerate(pairs))
    )
    if ordering == "fasta":
        return make_metacool(tmp_path / "map.metacool", pairs_file, fasta)
    path = tmp_path / "map.metacool"
    bins = pd.DataFrame({"chrom": list(contigs), "start": 0, "end": list(contigs.values())})
    metacooler.from_pairs(str(path), pairs=str(pairs_file), bins=bins, contig_genome_map=None, ps=False)
    return path


def _contig_data(tmp_path, contigs):
    contig_data = mtn.create_contig_data(str(_write_fasta(tmp_path / "assembly.fa", contigs)))
    return contig_data


PAIRS = [("a", 10, "b", 20), ("b", 30, "c", 5), ("a", 4, "a", 40)]


@pytest.mark.parametrize(
    "assembly, message",
    [
        (
            {"a": 100, "b": 200, "d": 300, "e": 40},
            "is not a contact map of the contigs of the assembly: 1 of its contig(s) are not in the assembly: 'c'. "
            "Build it from pairs of reads aligned on this assembly.",
        ),
        ({"a": 100, "b": 200}, "1 of its contig(s) are not in the assembly: 'c'. Build it from pairs"),
        (
            {"a": 100, "b": 204, "c": 400},
            "is not a contact map of the contigs of the assembly: 2 contig(s) have another length in it: 'b' (200 bp "
            "in the map, 204 bp in the assembly), 'c' (300 bp in the map, 400 bp in the assembly).",
        ),
        (
            {"a": 100, "b": 204, "c": 300, "z": 8},
            "is not a contact map of the contigs of the assembly: 1 contig(s) have another length in it: 'b' (200 bp "
            "in the map, 204 bp in the assembly).",
        ),
    ],
)
def test_map_of_other_contigs(assembly, message, tmp_path):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    with pytest.raises(mtmc.MetacoolError) as error:
        mtn.contig_data_from_map(str(metacool), _contig_data(tmp_path, assembly))
    assert message in str(error.value)
    assert str(error.value).startswith(f"{metacool} is not a contact map")
    # The contigs of the assembly that the map lacks are not a reason.
    assert "of the assembly are not in" not in str(error.value)


def test_contig_data_table_of_a_map_lacking_contigs(tmp_path, logged):
    """A map may lack contigs of the assembly (a map made with a contig
    filter, for example): they are no row of the contig data of the run, nor
    nodes of its network, which an INFO message says. The IDs are the
    positions of the contigs in the map."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    contig_data = _contig_data(tmp_path, {"z": 8, "a": 100, "b": 200, "y": 12, "c": 300})
    contacts = mtn.contig_contacts(str(metacool))
    assert contacts.to_dict("list") == {"id1": [0, 0, 1], "id2": [0, 1, 2], "count": [1, 1, 1]}
    names = ["a", "b", "c"]
    table = mtn.contig_data_table(contig_data, names, mtn.contig_hits(contacts, 3), metacool=str(metacool))
    assert list(table.Name) == names and list(table.ID) == [1, 2, 3] and list(table.Hit) == [3, 2, 1]
    assert logged() == [
        (
            "INFO",
            f"2 contig(s) of the assembly are not in {metacool}: they get no contacts, as contigs that no read maps "
            "to, and are in no table of the run: 'z', 'y'.",
        )
    ]
    # The contact of a with itself is no edge.
    network = mtn.map_network(str(metacool), weights=None, contacts=contacts)
    assert network.to_dict("list") == {"contig1": [1, 2], "contig2": [2, 3], "weight": [1, 1]}


def test_complete_map_lacks_no_contig(tmp_path, logged):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    contacts = mtn.contig_contacts(str(metacool))
    mtn.contig_data_table(
        _contig_data(tmp_path, {"c": 300, "a": 100, "b": 200}), "abc", mtn.contig_hits(contacts, 3), metacool=str(metacool)
    )
    assert logged() == []


def test_many_contigs_are_cut_short(tmp_path, logged):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    assembly = {"a": 100, "b": 200, "c": 300, **{f"x{i}": 8 for i in range(7)}}
    contacts = mtn.contig_contacts(str(metacool))
    mtn.contig_data_table(_contig_data(tmp_path, assembly), "abc", mtn.contig_hits(contacts, 3), metacool=str(metacool))
    assert logged() == [
        (
            "INFO",
            f"7 contig(s) of the assembly are not in {metacool}: they get no contacts, as contigs that no read maps "
            "to, and are in no table of the run: 'x0', 'x1', 'x2', 'x3', 'x4', ... (+2 more).",
        )
    ]
    lengths = {"a": 104, "b": 204, "c": 304, **{f"x{i}": 8 for i in range(7)}}
    (tmp_path / "many").mkdir()
    small = _small_map(tmp_path / "many", {**{f"x{i}": 4 for i in range(7)}, "a": 100, "b": 200, "c": 300}, PAIRS)
    with pytest.raises(
        mtmc.MetacoolError, match=r"10 contig\(s\) have another length in it: .*bp in the assembly\), \.\.\. \(\+5 more\)\.$"
    ):
        mtn.contig_data_from_map(str(small), _contig_data(tmp_path, lengths))
    foreign = r"7 of its contig\(s\) are not in the assembly: 'x0', 'x1', 'x2', 'x3', 'x4', \.\.\. \(\+2 more\)\. Build it "
    with pytest.raises(mtmc.MetacoolError, match=foreign):
        mtn.contig_data_from_map(str(small), _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300}))


def test_small_map_network(tmp_path):
    """Hits and network of a map whose contigs are not in the order of the
    assembly, longest first, as metacooler create pairs and hicstuff order
    them: the IDs are the positions of the contigs in the map, whatever the
    order of the assembly."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, [*PAIRS, ("c", 300, "a", 100)], ordering="length")
    names = [str(name) for name in metacooler.MetaCooler(str(metacool)).contigs["name"]]
    assert names == ["c", "b", "a"]
    contacts = mtn.contig_contacts(str(metacool))
    # c 0, b 1, a 2.
    assert contacts.to_dict("list") == {"id1": [0, 0, 1, 2], "id2": [1, 2, 2, 2], "count": [1, 1, 1, 1]}
    assert list(mtn.contig_hits(contacts, 3)) == [2, 2, 4]
    table = mtn.contig_data_table(_contig_data(tmp_path, {"a": 100, "b": 200, "c": 300}), names, [2, 2, 4])
    assert list(table.Name) == ["c", "b", "a"] and list(table.ID) == [1, 2, 3] and list(table.Size) == [300, 200, 100]
    network = mtn.map_network(str(metacool), weights=None, contacts=contacts)
    assert network.to_dict("list") == {"contig1": [1, 1, 2], "contig2": [2, 3, 3], "weight": [1, 1, 1]}


def _stored_network(tmp_path, halves=False, normalization="length"):
    """A map of the contigs a, b and c, with the network (normalization
    length) that metator network stores in it (see
    metator.metacool.store_network): one bin per contig, or two (halves)."""
    contigs = {"a": 100, "b": 200, "c": 300}
    if halves:
        pairs_file = tmp_path / "map.pairs"
        pairs_file.write_text(
            PAIRS_HEADER + "".join(f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t-\n" for i, (c1, p1, c2, p2) in enumerate(PAIRS))
        )
        bins = pd.DataFrame(
            {
                "chrom": [name for name in contigs for _ in range(2)],
                "start": [start for length in contigs.values() for start in (0, length // 2)],
                "end": [end for length in contigs.values() for end in (length // 2, length)],
            }
        )
        metacool = tmp_path / "map.metacool"
        metacooler.from_pairs(str(metacool), pairs=str(pairs_file), bins=bins, contig_genome_map=None, ps=False)
    else:
        metacool = _small_map(tmp_path, contigs, PAIRS)
    names = [str(name) for name in metacooler.MetaCooler(str(metacool)).contigs["name"]]
    hits = mtn.contig_hits(mtn.contig_contacts(str(metacool)), len(names))
    table = mtn.contig_data_table(_contig_data(tmp_path, contigs), names, hits)
    mtmc.store_network(str(metacool), str(metacool), mtmc.Network(table, normalization), force=True)
    return str(metacool)


def test_network_stored_in_a_map(tmp_path):
    """A map that metator network stored the network in gives it to the
    steps after it: the network of its contacts, weighed with the stored
    weights, its contig data and its normalization."""
    metacool = _stored_network(tmp_path)
    stored = mtn.stored_network(metacool)
    assert stored.normalization == "length"
    weights = mtn.stored_weights(metacool)
    contacts = mtn.contig_contacts(metacool)
    pd.testing.assert_frame_equal(stored.network, mtn.map_network(metacool, weights, contacts))
    pd.testing.assert_frame_equal(stored.network, mtn.map_network(metacool))
    pd.testing.assert_frame_equal(stored.table, mtn.map_contig_table(metacool))
    # a: one contact with b, and one within a, counted twice.
    assert list(stored.table.Hit) == [3, 2, 1]


@pytest.mark.parametrize("weights", ["None", "length", "empirical_hit", ""])
def test_map_network_takes_no_normalization_name(weights, tmp_path):
    """map_network weighs the contacts with the stored weights ("stored"),
    weights given, or none (None): the name of a normalization is no
    weights, and is refused rather than read as "stored"."""
    metacool = _stored_network(tmp_path)
    with pytest.raises(ValueError, match=f"^map_network weighs the contacts with weights='stored'.*not {weights!r}\\.$"):
        mtn.map_network(metacool, weights)
    counts = mtn.map_network(metacool, None)
    assert counts.to_dict("list") == {"contig1": [1, 2], "contig2": [2, 3], "weight": [1, 1]}
    assert not counts.equals(mtn.map_network(metacool, "stored"))


@pytest.mark.parametrize("normalization", ["None", "length"])
def test_map_whose_contacts_changed_since_its_network(normalization, tmp_path):
    """The Hit column of a map must be the hits of its contacts: the maps of
    two libraries with their networks, summed by metacooler merge, keep the
    contigs table of the first, and its Hit, stale. The steps refuse the
    stored network, with the step to run again, whatever its normalization:
    under None, no weights were stored that the merge would drop."""
    from click.testing import CliRunner
    from metacooler.cli.main import cli as metacooler_cli

    metacool = _stored_network(tmp_path, normalization=normalization)
    library = shutil.copyfile(metacool, tmp_path / "library_2.metacool")
    merged = str(tmp_path / "merged.metacool")
    # The command, which carries the metadata of the maps (metator's
    # network entry among it), as the docs have the libraries merged.
    result = CliRunner().invoke(metacooler_cli, ["merge", metacool, str(library), "-o", merged, "-q"])
    assert result.exit_code == 0, result.output
    with pytest.raises(mtmc.MetacoolError) as error:
        mtn.stored_network(merged)
    assert str(error.value) == (
        f"The Hit column of {merged} is not the hits of its contacts for 3 contig(s): 'a', 'b', 'c'. Its contacts "
        "changed since metator network stored the network (a merge of maps keeps the contigs table of the first): "
        "run metator network on the map again."
    )


def test_map_whose_network_column_was_dropped(tmp_path):
    """A map whose Hit column an edit dropped since metator network stored
    it is refused, with the step to run again."""
    metacool = _stored_network(tmp_path)
    metacooler.annotate(metacool, drop=["contigs.Hit"], inplace=True)
    message = (
        f"{metacool} has no Hit column in its contigs table, which metator network stores (an edit of the map dropped "
        "it since): run metator network on the map again."
    )
    for read in (mtn.stored_network, mtn.map_contig_table):
        with pytest.raises(mtmc.MetacoolError) as error:
            read(metacool)
        assert str(error.value) == message


def test_map_whose_network_weights_were_dropped(tmp_path):
    """A map that lacks the weights of the normalization of its network is
    refused, with the step to run again."""
    metacool = _stored_network(tmp_path)
    with h5py.File(metacool, "r+") as group:
        del group["bins/metator_length"]
    message = (
        f"{metacool} lacks the weights metator_length of the length normalization of its network, which an edit of "
        "the map dropped since metator network stored them: run metator network on the map again."
    )
    for read in (mtn.stored_network, mtn.map_network):
        with pytest.raises(mtmc.MetacoolError) as error:
            read(metacool)
        assert str(error.value) == message


def test_weights_of_a_contig_that_differ_between_its_bins(tmp_path):
    """The weights of a map with several bins per contig are those of its
    contigs, the weight of each of their bins: weights that differ between
    the bins of a contig are not those that metator network stores."""
    metacool = _stored_network(tmp_path, halves=True)
    weights = mtn.stored_weights(metacool)
    assert len(weights) == 3
    with h5py.File(metacool, "r+") as group:
        group["bins/metator_length"][0] = 2 * group["bins/metator_length"][0]
    for read in (mtn.stored_weights, mtn.stored_network):
        with pytest.raises(mtmc.MetacoolError, match=r"metator_length of .* differ between the bins of a contig"):
            read(metacool)


def test_contigs_table_with_text_in_a_network_column(tmp_path):
    """A contigs table whose GC_content holds text, which metator network
    stores numbers in, is refused."""
    metacool = _stored_network(tmp_path)
    with h5py.File(metacool, "r+") as group:
        del group["contigs/GC_content"]
        group["contigs/GC_content"] = np.array([b"high", b"low", b"low"])
    for read in (mtn.map_contig_table, mtn.stored_network):
        with pytest.raises(mtmc.MetacoolError, match=r"holds a value that is not a number in a column that metator network"):
            read(metacool)


@pytest.mark.parametrize(
    "uri, message",
    [
        ("missing.metacool", "Cannot read the contact map missing.metacool: File not found: missing.metacool"),
        (str(TESTS_DATA / "subset.pairs.gz"), f"Cannot read the contact map {TESTS_DATA / 'subset.pairs.gz'}: "),
        ("maps.mmetacool", "Cannot read the contact map maps.mmetacool: No metacool entry at root of maps.mmetacool."),
        ("maps.mmetacool::entries/lib2", "Cannot read the contact map maps.mmetacool::entries/lib2: Group "),
        # A cooler that is not a map (what metacooler extract writes).
        ("plain.cool", "Cannot read the contact map plain.cool: "),
    ],
)
def test_unreadable_map(uri, message, tmp_path):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    metacooler.create_mmetacool(Path("maps.mmetacool"), {"lib1": metacool})
    cooler.create_cooler(
        "plain.cool",
        pd.DataFrame({"chrom": ["a"], "start": [0], "end": [100]}),
        pd.DataFrame({"bin1_id": [0], "bin2_id": [0], "count": [1]}),
    )
    with pytest.raises(mtmc.MetacoolError) as error:
        mtn.contig_contacts(uri)
    assert str(error.value).startswith(message)


def test_map_without_contacts(tmp_path):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, [])
    contacts = mtn.contig_contacts(str(metacool))
    assert contacts.empty and list(mtn.contig_hits(contacts, 3)) == [0, 0, 0]
    with pytest.raises(mtmc.MetacoolError, match=r"map\.metacool holds no contact: there is no network to build\."):
        mtn.log_contacts(str(metacool), contacts)


def test_contacts_are_logged(tmp_path, logged):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    mtn.log_contacts("map.metacool", mtn.contig_contacts(str(metacool)))
    assert [message for _, message in logged()] == [
        "Information of map.metacool:",
        "3 contacts in the library.",
        "2 contacts inter-contigs in the library.",
        f"3D ratio : {2 / 3}\n",
    ]


def test_map_of_a_square_matrix(tmp_path):
    """A map that stores both triangles of the matrix (which metacooler
    opens, but never writes) would count some contacts twice: refused."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    with h5py.File(metacool, "r+") as group:
        group.attrs["storage-mode"] = "square"
    assert metacooler.MetaCooler(str(metacool)).info["storage-mode"] == "square"
    with pytest.raises(
        mtmc.MetacoolError,
        match=r"map\.metacool stores its contacts with storage-mode 'square': a network is built from a map that stores "
        r"each contact once, in the upper triangle \(storage-mode 'symmetric-upper'\), as metacooler writes them\.$",
    ):
        mtn.contig_contacts(str(metacool))


def test_map_of_float_counts(tmp_path):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    with h5py.File(metacool, "r+") as group:
        counts = group["pixels"]["count"][:]
        del group["pixels"]["count"]
        group["pixels"].create_dataset("count", data=counts.astype(float))
    with pytest.raises(mtmc.MetacoolError, match="holds contact counts of type float64: a network is built from counts"):
        mtn.contig_contacts(str(metacool))


# The restriction sites and the coverage of the contigs, read from the
# Restriction_site and Shotgun_coverage columns of a map where --enzyme and
# --depth measure none (metator.network.contig_data_from_map).

NORMALIZATIONS = ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"]


def _synthetic(metagenome_files):
    return {
        "pairs": metagenome_files["pairs"],
        "assembly": metagenome_files["assembly"],
        "depth": str(metagenome_files["depth"]),
        "enzyme": "DpnII",
    }


def _map_with_columns(path, data):
    """The map of the pairs of a data set, with the Restriction_site and
    Shotgun_coverage columns that metaTOR stored in the maps it built (see
    tests.maps.store_assembly_columns): measured on the assembly with the
    enzyme and the depth file of the data set."""
    metacool = make_metacool(path, data["pairs"], data["assembly"])
    contig_data = mtn.create_contig_data(str(data["assembly"]), data["depth"], data["enzyme"])
    store_assembly_columns(metacool, contig_data)
    return metacool


def _annotated(metacool, tmp_path, columns):
    """The map, with columns (by column, the values of each contig by name;
    a contig left out has none) in its contigs table, as metacooler
    annotate adds them from a table."""
    names = list(dict.fromkeys(name for values in columns.values() for name in values))
    table = tmp_path / "columns.txt"
    with open(table, "w") as handle:
        handle.write("\t".join(["Name", *columns]) + "\n")
        for name in names:
            handle.write("\t".join([name, *(str(values.get(name, "")) for values in columns.values())]) + "\n")
    metacooler.annotate(str(metacool), contigs_table=str(table), contigs_table_format="Name", unmatched="ignore",
                        inplace=True)  # fmt: skip
    return metacool


def _sites_and_coverage(contig_data):
    return {name: (data["RS"], data["coverage"]) for name, data in contig_data.items()}


@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_network_from_the_columns_of_a_map(metagenome_files, normalization, tmp_path):
    """A map that holds the restriction sites and the coverage of its
    contigs gives, without --enzyme and --depth, the contig data and the
    network that they give."""
    data = _synthetic(metagenome_files)
    metacool = _map_with_columns(tmp_path / "columns.metacool", data)
    names = [str(name) for name in metacooler.MetaCooler(str(metacool)).contigs["name"]]
    contacts = mtn.contig_contacts(str(metacool))
    hits = mtn.contig_hits(contacts, len(names))
    tables = {}
    for how in ["options", "map"]:
        if how == "options":
            contig_data = mtn.create_contig_data(str(data["assembly"]), data["depth"], data["enzyme"])
        else:
            contig_data = mtn.create_contig_data(str(data["assembly"]))
            assert mtn.contig_data_from_map(str(metacool), contig_data) == ["Restriction_site", "Shotgun_coverage"]
        tables[how] = mtn.contig_data_table(contig_data, names, hits)
    pd.testing.assert_frame_equal(tables["map"], tables["options"])
    assert (tables["map"].Restriction_site >= 1).all() and (tables["map"].Shotgun_coverage > 0).all()
    networks = [mtn.map_network(str(metacool), mtn.network_weights(table, normalization), contacts) for table in tables.values()]
    pd.testing.assert_frame_equal(networks[0], networks[1])


def test_columns_of_a_map_read(tmp_path, logged):
    """The values of the map, as create_contig_data measures them: an int
    of restriction sites, a float of coverage (from a column of integers
    too), whatever the order of the contigs in the map."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    _annotated(metacool, tmp_path, {"Restriction_site": {"c": 7, "a": 1, "b": 3},
                                    "Shotgun_coverage": {"b": 2, "c": 0, "a": 12}})  # fmt: skip
    contig_data = _contig_data(tmp_path, {"c": 300, "a": 100, "b": 200})
    assert mtn.map_columns(str(metacool)) == ["Restriction_site", "Shotgun_coverage"]
    assert mtn.contig_data_from_map(str(metacool), contig_data) == ["Restriction_site", "Shotgun_coverage"]
    assert _sites_and_coverage(contig_data) == {"c": (7, 0.0), "a": (1, 12.0), "b": (3, 2.0)}
    assert [(type(data["RS"]), type(data["coverage"])) for data in contig_data.values()] == [(int, float)] * 3
    assert logged() == [
        ("INFO", f"The restriction sites of the contigs: the Restriction_site column of {metacool}."),
        ("INFO", f"The shotgun coverage of the contigs: the Shotgun_coverage column of {metacool}."),
    ]


def test_options_replace_the_columns_of_a_map(tmp_path, logged):
    """--enzyme and --depth replace the columns of the map: the contigs keep
    what they measure on the assembly, and a warning counts the contigs
    whose value in the map is another. None when they agree."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    _annotated(metacool, tmp_path, {"Restriction_site": {"a": 1, "b": 3, "c": 7},
                                    "Shotgun_coverage": {"a": 12.5, "b": 2.0, "c": 0.0}})  # fmt: skip
    contig_data = _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300})
    measured = {"a": (1, 12.5), "b": (2, 2.0), "c": (5, 0.0)}
    for name, (sites, coverage) in measured.items():
        contig_data[name].update(RS=sites, coverage=coverage)
    assert mtn.contig_data_from_map(str(metacool), contig_data, "DpnII", "depth.txt") == []
    assert _sites_and_coverage(contig_data) == measured
    assert logged() == [
        (
            "WARNING",
            f"--enzyme replaces the Restriction_site column of {metacool}: 2 contig(s) have other restriction sites "
            "in it: 'b', 'c'.",
        )
    ]


@pytest.mark.parametrize(
    "column, values, lacking",
    [
        ("Restriction_site", {"a": 1, "b": 3}, "'c'"),
        ("Restriction_site", {"a": 1, "b": 0, "c": 7}, "'b'"),
        ("Restriction_site", {"a": 1, "b": 2.5, "c": -1}, "'b', 'c'"),
        ("Restriction_site", {"a": "many", "b": 2, "c": 7}, "'a'"),
        ("Shotgun_coverage", {"a": 1.5, "b": -0.5, "c": 7}, "'b'"),
        ("Shotgun_coverage", {"a": 1.5, "b": "nan", "c": "inf"}, "'b', 'c'"),
        ("Shotgun_coverage", {"a": 1.5, "c": 2}, "'b'"),
    ],
    ids=["missing", "zero sites", "not an integer", "text", "negative", "not finite", "missing coverage"],
)
def test_columns_of_a_map_that_lack_contigs(column, values, lacking, tmp_path, logged):
    """A column that does not give every contig of the assembly a value (as
    create_contig_data measures them) is not read, as without the option,
    and a warning says which contigs it lacks."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    _annotated(metacool, tmp_path, {column: values})
    contig_data = _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300})
    assert mtn.contig_data_from_map(str(metacool), contig_data) == []
    assert set(_sites_and_coverage(contig_data).values()) == {("-", "-")}
    what, option, valid = {
        "Restriction_site": ("restriction sites", "--enzyme", "an integer of at least 1"),
        "Shotgun_coverage": ("shotgun coverage", "--depth", "a number of at least 0"),
    }[column]
    count = lacking.count("'") // 2
    assert logged() == [
        (
            "WARNING",
            f"The {column} column of {metacool} gives no {what} to {count} contig(s) of the assembly (not in the map, "
            f"or a value that is missing or not {valid}): {lacking}. The run does not read it, as without {option}.",
        )
    ]


def test_columns_of_a_map_lacking_contigs_of_the_assembly(tmp_path, logged):
    """A map that lacks contigs of the assembly (a map made with a contig
    filter) cannot give them a value: its columns are not read."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    _annotated(metacool, tmp_path, {"Restriction_site": {"a": 1, "b": 3, "c": 7}})
    contig_data = _contig_data(tmp_path, {"z": 8, "a": 100, "b": 200, "c": 300})
    assert mtn.contig_data_from_map(str(metacool), contig_data) == []
    assert contig_data["a"]["RS"] == "-"
    assert [level for level, message in logged()] == ["WARNING"]
    assert "gives no restriction sites to 1 contig(s) of the assembly" in logged()[0][1]


def test_columns_of_a_map_without_values(tmp_path, logged):
    """A column of ``-``, as metaTOR writes what it did not measure (a map
    made from contig_data_final.txt, metacooler create pairs
    --contigs-format metator), is as no column: no warning. A guard against
    reading too much: a version that reads nothing passes it."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    _annotated(metacool, tmp_path, {"Restriction_site": dict.fromkeys("abc", "-"),
                                    "Shotgun_coverage": dict.fromkeys("abc", "-")})  # fmt: skip
    assert "Restriction_site" in metacooler.MetaCooler(str(metacool)).contigs.columns
    contig_data = _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300})
    assert mtn.map_columns(str(metacool)) == []
    assert mtn.contig_data_from_map(str(metacool), contig_data, "DpnII") == []
    assert logged() == []


def test_columns_of_a_map_without_them(tmp_path, logged):
    """A guard against reading too much: a version that reads nothing
    passes it."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    contig_data = _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300})
    assert mtn.map_columns(str(metacool)) == []
    assert mtn.contig_data_from_map(str(metacool), contig_data) == []
    assert set(_sites_and_coverage(contig_data).values()) == {("-", "-")} and logged() == []


def test_columns_of_an_mmetacool_entry(tmp_path):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    _annotated(metacool, tmp_path, {"Restriction_site": {"a": 1, "b": 3, "c": 7}})
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    contig_data = _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300})
    uri = f"{collection}::entries/lib1"
    assert mtn.map_columns(uri) == ["Restriction_site"]
    assert mtn.contig_data_from_map(uri, contig_data) == ["Restriction_site"]
    assert [data["RS"] for data in contig_data.values()] == [1, 3, 7]


def test_columns_of_a_map_with_a_contig_of_no_length(tmp_path, logged):
    """A contig of no length, which no map holds (metacooler leaves it out),
    does not keep the columns of the map from being read: it gets what the
    options would give it: no site, so a Restriction_site of 1, and a
    coverage of 0."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    _annotated(metacool, tmp_path, {"Restriction_site": {"a": 1, "b": 3, "c": 7},
                                    "Shotgun_coverage": {"a": 1.5, "b": 2.0, "c": 0.5}})  # fmt: skip
    contig_data = _contig_data(tmp_path, {"a": 100, "e": 0, "b": 200, "c": 300})
    assert contig_data["e"]["length"] == 0
    assert mtn.contig_data_from_map(str(metacool), contig_data) == ["Restriction_site", "Shotgun_coverage"]
    assert _sites_and_coverage(contig_data) == {"a": (1, 1.5), "e": (1, 0.0), "b": (3, 2.0), "c": (7, 0.5)}
    assert [level for level, message in logged()] == ["INFO", "INFO"]


def test_columns_of_a_map_of_other_contigs(tmp_path, logged):
    """A map that holds contigs that the assembly lacks is refused for what
    it is, before its columns are read: no warning about them (it named
    the contigs that its columns lacked, and only the network step said
    why)."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    _annotated(metacool, tmp_path, {"Restriction_site": {"a": 1, "b": 3, "c": 7}})
    contig_data = _contig_data(tmp_path, {"a": 100, "b": 200, "x": 300})
    with pytest.raises(mtmc.MetacoolError, match="is not a contact map of the contigs of the assembly: 1 of its contig"):
        mtn.contig_data_from_map(str(metacool), contig_data)
    assert logged() == []


@pytest.mark.parametrize(
    "read", [mtn.map_columns, lambda map: mtn.contig_data_from_map(map, {})], ids=["map_columns", "contig_data_from_map"]
)
def test_columns_of_an_unreadable_map(read, tmp_path):
    with pytest.raises(mtmc.MetacoolError, match="Cannot read the contact map"):
        read(str(tmp_path / "missing.metacool"))
