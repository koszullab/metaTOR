# Tests of the network built from a .metacool contact map
# (metator.network.metacool_to_contacts), against the network of the pairs
# the map is made of, counted straight from the pairs file (pairs_network),
# on the mock metagenome of tests/synthetic.py and on tests_data.

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
    for contig1, pos1, contig2, pos2 in read_pairs(pairs_file):
        id1, id2 = contig_data[contig1]["id"], contig_data[contig2]["id"]
        if id1 == id2:
            continue
        edge = (min(id1, id2), max(id1, id2))
        if pos1 == 0 or pos2 == 0:
            position_0[edge] += 1
        elif pos1 == contig_data[contig1]["length"] or pos2 == contig_data[contig2]["length"]:
            last_base[edge] += 1
    return last_base, position_0


def pairs_network(data, out, normalization="None", self_contacts=False):
    """network.txt and contig_data_network.txt of the pairs themselves,
    written in out: the oracle of the network of a map of the pairs. The
    contacts and the hits are counted in Python, from the pairs file: each
    pair whose two reads lie on contigs of the assembly, within their
    length (at position 0 too, where metacooler puts the read on the first
    base), counts once, and once in the hits of each of its two contigs (a
    pair within a contig, twice in its hits). The edges are written by the
    writer of the network of a map (metator.network._write_map_network),
    whose order, normalization and format the test does not check: what it
    checks is that the pixels of the map give the contacts and the hits of
    its pairs."""
    contig_data = mtn.create_contig_data(str(data["assembly"]), data["depth"], data["enzyme"])
    contacts, hits = collections.Counter(), collections.Counter()
    for contig1, pos1, contig2, pos2 in read_pairs(data["pairs"]):
        if contig1 not in contig_data or contig2 not in contig_data:
            continue
        if pos1 > contig_data[contig1]["length"] or pos2 > contig_data[contig2]["length"]:
            continue
        id1, id2 = contig_data[contig1]["id"], contig_data[contig2]["id"]
        contacts[(min(id1, id2), max(id1, id2))] += 1
        hits[id1] += 1
        hits[id2] += 1
    names = np.empty(len(contig_data) + 1, dtype=object)
    for name, contig in contig_data.items():
        contig["hit"] = hits[contig["id"]]
        names[contig["id"]] = name
    edges = sorted(contacts.items())
    table = pd.DataFrame(
        {"id1": [edge[0] for edge, _ in edges], "id2": [edge[1] for edge, _ in edges], "count": [n for _, n in edges]}
    )
    out.mkdir(parents=True)
    network, contigs = out / "network.txt", out / "contig_data_network.txt"
    mtn._write_map_network(table, names, contig_data, str(network), str(out / "pre.txt"), str(out / "sorted.txt"),
                           str(out), 1, normalization, self_contacts, mtn.METACOOL_CHUNKSIZE)  # fmt: skip
    mtn.write_contig_data(contig_data, str(contigs))
    return {"network": network, "contig_data": contigs, "dir": out}


def both_networks(data, metacool, out, normalization="None", self_contacts=False):
    """network.txt and contig_data_network.txt of the pairs (see
    pairs_network) and of the metacool, written in out/pairs and
    out/metacool."""
    from_pairs = pairs_network(data, out / "pairs", normalization, self_contacts)
    directory = out / "metacool"
    directory.mkdir(parents=True)
    contig_data = mtn.create_contig_data(str(data["assembly"]), data["depth"], data["enzyme"])
    network, contigs = mtn.metacool_to_contacts(
        str(metacool), contig_data, str(directory), "network.txt", "contig_data_network.txt", str(directory), 1,
        normalization, self_contacts,
    )  # fmt: skip
    return from_pairs, {"network": Path(network), "contig_data": Path(contigs), "dir": directory}


@pytest.fixture(params=["synthetic", "alignment"])
def dataset(request, metagenome_files):
    """Pairs, assembly, depth file and enzyme of a data set. The pairs of
    subset.pairs.gz, the largest, go through test_commands.py
    (test_network_from_a_metacool) and the tests by chunk."""
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
    """The network of a metacool is that of its pairs, byte for byte, with
    the same contig data: the pairs with a read on the last base of a contig
    (the mock metagenome) or at position 0 (alignment.pairs) count. A
    metacool sums its libraries: no network per library."""
    metacool = make_metacool(tmp_path / "map.metacool", dataset["pairs"], dataset["assembly"], binsize)
    from_pairs, from_map = both_networks(dataset, metacool, tmp_path / "raw")
    assert from_map["network"].read_text() == from_pairs["network"].read_text()
    assert from_map["contig_data"].read_text() == from_pairs["contig_data"].read_text()
    # What each data set holds: the contacts on the edges of the contigs.
    contig_data = mtn.create_contig_data(str(dataset["assembly"]))
    last_base, position_0 = outside_contacts(dataset["pairs"], contig_data)
    assert (bool(last_base), bool(position_0)) == ((True, False) if dataset["name"] == "synthetic" else (False, True))
    assert sorted(path.name for path in from_map["dir"].glob("network*.txt")) == ["network.txt"]


def test_normalized_network_of_a_metacool(dataset, tmp_path):
    """Normalized (empirical_hit, the default): the same network, byte for
    byte."""
    metacool = make_metacool(tmp_path / "map.metacool", dataset["pairs"], dataset["assembly"], 1000)
    from_pairs, from_map = both_networks(dataset, metacool, tmp_path / "norm", normalization="empirical_hit")
    assert from_map["network"].read_text() == from_pairs["network"].read_text()


@pytest.mark.parametrize("binsize", BIN_SIZES)
def test_network_of_a_metacool_with_self_contacts(dataset, binsize, tmp_path):
    """With the intra-contig contacts: the same network, byte for byte."""
    metacool = make_metacool(tmp_path / "map.metacool", dataset["pairs"], dataset["assembly"], binsize)
    from_pairs, from_map = both_networks(dataset, metacool, tmp_path / "self", self_contacts=True)
    assert from_map["network"].read_text() == from_pairs["network"].read_text()
    assert from_map["contig_data"].read_text() == from_pairs["contig_data"].read_text()
    assert any(line.split("\t")[0] == line.split("\t")[1] for line in from_map["network"].read_text().splitlines())


@pytest.mark.parametrize("binsize", BIN_SIZES)
@pytest.mark.parametrize(
    "dataset, how", [("synthetic", "create"), ("synthetic", "extract"), ("alignment", "create")], indirect=["dataset"]
)
def test_network_of_a_filtered_metacool(dataset, how, binsize, tmp_path):
    """A map that lacks contigs of the assembly, made with a contig filter,
    which leaves out the pairs with a read on the contigs it drops: its
    network is that of the pairs it holds, whose two reads lie on the
    contigs it keeps, and the contigs it lacks get no hits and no edges, as
    contigs that no read maps to: the same contig data and network, byte for
    byte, with and without self contacts."""
    contig_data = mtn.create_contig_data(str(dataset["assembly"]))
    lengths = {name: data["length"] for name, data in contig_data.items()}
    shortest = int(np.median(list(lengths.values())))
    kept = {name for name, length in lengths.items() if length >= shortest}
    metacool = filtered_metacool(tmp_path / "map.metacool", dataset["pairs"], dataset["assembly"], binsize, shortest, how)
    assert sorted(metacooler.MetaCooler(str(metacool)).contigs.name) == sorted(kept)
    assert 0 < len(kept) < len(lengths)
    pairs = {**dataset, "pairs": kept_pairs(dataset["pairs"], kept, tmp_path / "kept.pairs")}

    for out, self_contacts in [("raw", False), ("self", True)]:
        from_pairs, from_map = both_networks(pairs, metacool, tmp_path / out, self_contacts=self_contacts)
        assert from_map["network"].read_text() == from_pairs["network"].read_text()
        assert from_map["contig_data"].read_text() == from_pairs["contig_data"].read_text()
    table = pd.read_csv(from_map["contig_data"], sep="\t")
    assert (table.Hit[~table.Name.isin(kept)] == 0).all()
    assert table.Hit[table.Name.isin(kept)].sum() > 0


def test_contacts_by_chunk(tmp_path):
    """The contacts do not depend on the number of pixels read at a time
    (a pixel at a time with 1)."""
    pairs, assembly = TESTS_DATA_PAIRS["alignment"]
    metacool = make_metacool(tmp_path / "map.metacool", pairs, assembly)
    contig_data = mtn.create_contig_data(str(assembly))
    hits, contacts = mtn.read_metacool_contacts(str(metacool), contig_data)
    for chunksize in [1, 997]:
        chunk_hits, chunk_contacts = mtn.read_metacool_contacts(str(metacool), contig_data, chunksize=chunksize)
        assert (chunk_hits == hits).all()
        pd.testing.assert_frame_equal(chunk_contacts, contacts)
    # Sorted, one row per contig pair, the upper triangle.
    assert (contacts.id1 <= contacts.id2).all()
    assert not contacts.duplicated(["id1", "id2"]).any()
    assert contacts.equals(contacts.sort_values(["id1", "id2"]))
    # Every pair counts twice in the hits.
    assert hits[0] == 0 and hits.sum() == 2 * contacts["count"].sum()


@pytest.mark.parametrize("self_contacts", [False, True])
def test_network_by_chunk(self_contacts, tmp_path):
    """Nor does the network, whose edges are written a block at a time.

    An equivalence guard: it passes without the change that writes the
    edges a block at a time, which keeps the network as it was."""
    pairs, assembly = TESTS_DATA_PAIRS["subset"]
    metacool = make_metacool(tmp_path / "map.metacool", pairs, assembly, 1000)
    networks = []
    for chunksize in [mtn.METACOOL_CHUNKSIZE, 9973]:
        contig_data = mtn.create_contig_data(str(assembly))
        out = tmp_path / str(chunksize)
        out.mkdir()
        network, contigs = mtn.metacool_to_contacts(
            str(metacool), contig_data, str(out), "network.txt", "contigs.txt", str(out), 1, "empirical_hit",
            self_contacts, chunksize=chunksize
        )  # fmt: skip
        networks.append((Path(network).read_text(), Path(contigs).read_text()))
    assert networks[0] == networks[1]
    assert len(networks[0][0].splitlines()) > 2 * 9973


def test_network_of_an_mmetacool_entry(tmp_path):
    """An entry of a .mmetacool collection, given by its URI."""
    pairs, assembly = TESTS_DATA_PAIRS["subset"]
    metacool = make_metacool(tmp_path / "map.metacool", pairs, assembly, 1000)
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    networks = []
    for i, uri in enumerate([str(metacool), f"{collection}::entries/lib1"]):
        contig_data = mtn.create_contig_data(str(assembly))
        out = tmp_path / f"out_{i}"
        out.mkdir()
        network, _ = mtn.metacool_to_contacts(uri, contig_data, str(out), "network.txt", "contig_data_network.txt",
                                              str(out), 1, "empirical_hit")  # fmt: skip
        networks.append(Path(network).read_text())
    assert networks[0] == networks[1] and networks[0]


# Maps that cannot give the network of an assembly.


def _write_fasta(path, contigs):
    path.write_text("".join(f">{name}\n{'ACGT' * (length // 4)}\n" for name, length in contigs.items()))
    return path


def _small_map(tmp_path, contigs, pairs):
    """A one-bin-per-contig map of the given contigs and pairs."""
    fasta = _write_fasta(tmp_path / "contigs.fa", contigs)
    pairs_file = tmp_path / "map.pairs"
    pairs_file.write_text(
        PAIRS_HEADER + "".join(f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t-\n" for i, (c1, p1, c2, p2) in enumerate(pairs))
    )
    return make_metacool(tmp_path / "map.metacool", pairs_file, fasta)


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
        mtn.read_metacool_contacts(str(metacool), _contig_data(tmp_path, assembly))
    assert message in str(error.value)
    assert str(error.value).startswith(f"{metacool} is not a contact map")
    # The contigs of the assembly that the map lacks are not a reason.
    assert "of the assembly are not in" not in str(error.value)


def test_map_lacking_contigs(tmp_path, logged):
    """A map may lack contigs of the assembly (a map made with a contig
    filter, for example): they get no contacts and no hits, as contigs that
    no read maps to, which an INFO message says. The ids are those of the
    assembly."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    # Ids: z 1, a 2, b 3, y 4, c 5.
    contig_data = _contig_data(tmp_path, {"z": 8, "a": 100, "b": 200, "y": 12, "c": 300})
    hits, contacts = mtn.read_metacool_contacts(str(metacool), contig_data)
    assert list(hits) == [0, 0, 3, 2, 0, 1]
    assert contacts.to_dict("list") == {"id1": [2, 2, 3], "id2": [2, 3, 5], "count": [1, 1, 1]}
    message = (
        f"2 contig(s) of the assembly are not in {metacool}: they get no contacts, as contigs that no read maps to: 'z', 'y'."
    )
    assert logged() == [("INFO", message)]
    network, contigs = mtn.metacool_to_contacts(str(metacool), contig_data, str(tmp_path), "network.txt",
                                                "contigs.txt", str(tmp_path), 1, "None", True)  # fmt: skip
    assert Path(network).read_text() == "2\t2\t1\n2\t3\t1\n3\t5\t1\n"
    table = pd.read_csv(contigs, sep="\t")
    assert list(table.Name) == ["z", "a", "b", "y", "c"] and list(table.Hit) == [0, 3, 2, 0, 1]
    assert [entry for entry in logged() if "of the assembly are not in" in entry[1]] == [("INFO", message)] * 2


def test_complete_map_lacks_no_contig(tmp_path, logged):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    mtn.read_metacool_contacts(str(metacool), _contig_data(tmp_path, {"c": 300, "a": 100, "b": 200}))
    assert logged() == []


def test_many_contigs_are_cut_short(tmp_path, logged):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    assembly = {"a": 100, "b": 200, "c": 300, **{f"x{i}": 8 for i in range(7)}}
    mtn.read_metacool_contacts(str(metacool), _contig_data(tmp_path, assembly))
    assert logged() == [
        (
            "INFO",
            f"7 contig(s) of the assembly are not in {metacool}: they get no contacts, as contigs that no read maps "
            "to: 'x0', 'x1', 'x2', 'x3', 'x4', ... (+2 more).",
        )
    ]
    lengths = {"a": 104, "b": 204, "c": 304, **{f"x{i}": 8 for i in range(7)}}
    (tmp_path / "many").mkdir()
    small = _small_map(tmp_path / "many", {**{f"x{i}": 4 for i in range(7)}, "a": 100, "b": 200, "c": 300}, PAIRS)
    with pytest.raises(
        mtmc.MetacoolError, match=r"10 contig\(s\) have another length in it: .*bp in the assembly\), \.\.\. \(\+5 more\)\.$"
    ):
        mtn.read_metacool_contacts(str(small), _contig_data(tmp_path, lengths))
    foreign = r"7 of its contig\(s\) are not in the assembly: 'x0', 'x1', 'x2', 'x3', 'x4', \.\.\. \(\+2 more\)\. Build it "
    with pytest.raises(mtmc.MetacoolError, match=foreign):
        mtn.read_metacool_contacts(str(small), _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300}))


def test_small_map_network(tmp_path):
    """Hits and network of a map whose contigs are not in the order of the
    assembly: ids are those of the assembly."""
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, [*PAIRS, ("c", 300, "a", 100)])
    contig_data = _contig_data(tmp_path, {"c": 300, "b": 200, "a": 100})
    hits, contacts = mtn.read_metacool_contacts(str(metacool), contig_data)
    # Ids: c 1, b 2, a 3.
    assert list(hits) == [0, 2, 2, 4]
    assert contacts.to_dict("list") == {"id1": [1, 1, 2, 3], "id2": [2, 3, 3, 3], "count": [1, 1, 1, 1]}
    network, _ = mtn.metacool_to_contacts(str(metacool), contig_data, str(tmp_path), "network.txt", "contigs.txt",
                                          str(tmp_path), 1, "None")  # fmt: skip
    # In the order of the contig names, as from pairs.
    assert Path(network).read_text() == "2\t3\t1\n1\t3\t1\n1\t2\t1\n"
    assert [contig_data[name]["hit"] for name in "abc"] == [4, 2, 2]


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
        mtn.read_metacool_contacts(uri, _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300}))
    assert str(error.value).startswith(message)


def test_map_without_contacts(tmp_path):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, [])
    contig_data = _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300})
    with pytest.raises(mtmc.MetacoolError, match=r"map\.metacool holds no contact: there is no network to build\."):
        mtn.metacool_to_contacts(str(metacool), contig_data, str(tmp_path), "network.txt", "contigs.txt", str(tmp_path),
                                 1, "None")  # fmt: skip


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
        mtn.read_metacool_contacts(str(metacool), _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300}))


def test_map_of_float_counts(tmp_path):
    metacool = _small_map(tmp_path, {"a": 100, "b": 200, "c": 300}, PAIRS)
    with h5py.File(metacool, "r+") as group:
        counts = group["pixels"]["count"][:]
        del group["pixels"]["count"]
        group["pixels"].create_dataset("count", data=counts.astype(float))
    with pytest.raises(mtmc.MetacoolError, match="holds contact counts of type float64: a network is built from counts"):
        mtn.read_metacool_contacts(str(metacool), _contig_data(tmp_path, {"a": 100, "b": 200, "c": 300}))


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
    contigs gives, without --enzyme and --depth, the network that they
    give, byte for byte, and the same contig data."""
    data = _synthetic(metagenome_files)
    metacool = _map_with_columns(tmp_path / "columns.metacool", data)
    outputs = {}
    for how in ["options", "map"]:
        directory = tmp_path / how
        directory.mkdir()
        if how == "options":
            contig_data = mtn.create_contig_data(str(data["assembly"]), data["depth"], data["enzyme"])
        else:
            contig_data = mtn.create_contig_data(str(data["assembly"]))
            assert mtn.contig_data_from_map(str(metacool), contig_data) == ["Restriction_site", "Shotgun_coverage"]
        outputs[how] = mtn.metacool_to_contacts(str(metacool), contig_data, str(directory), "network.txt",
                                                "contig_data_network.txt", str(directory), 1, normalization)  # fmt: skip
    for by_options, by_map in zip(outputs["options"], outputs["map"], strict=True):
        assert Path(by_map).read_bytes() == Path(by_options).read_bytes()
    table = pd.read_csv(outputs["map"][1], sep="\t")
    assert (table.Restriction_site >= 1).all() and (table.Shotgun_coverage > 0).all()


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
