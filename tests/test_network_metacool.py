# Tests of the network built from a .metacool contact map
# (metator.network.metacool_to_contacts), against the network of the pairs
# the map is made of (metator.network.alignment_to_contacts), on the mock
# metagenome of tests/synthetic.py and on tests_data.

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
        contig_data, _ = mtn.create_contig_data(str(assembly))
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
    those with a read at position 0, by contig pair (ids): the pairs that
    the network of the pairs leaves out (its edge test), and the metacool
    keeps. Counted from the pairs file."""
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


def both_networks(data, metacool, out, normalization="None", self_contacts=False):
    """network.txt and contig_data_network.txt of the pairs and of the
    metacool, written in out/pairs and out/metacool."""
    outputs = {}
    for start in ["pairs", "metacool"]:
        directory = out / start
        directory.mkdir(parents=True)
        contig_data, hit_data = mtn.create_contig_data(str(data["assembly"]), 1, data["depth"], data["enzyme"])
        outputs_args = (str(directory), "network.txt", "contig_data_network.txt", str(directory), 1, normalization)
        if start == "pairs":
            network, contigs = mtn.alignment_to_contacts(
                [str(data["pairs"])], contig_data, 0, hit_data, *outputs_args, self_contacts
            )
        else:
            network, contigs = mtn.metacool_to_contacts(str(metacool), contig_data, *outputs_args, self_contacts)
        outputs[start] = {"network": Path(network), "contig_data": Path(contigs), "dir": directory}
    return outputs["pairs"], outputs["metacool"]


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
    """The network of a metacool is that of its pairs at edge 0, byte for
    byte, with the same contig data: the pairs with a read on the last base
    of a contig (the mock metagenome) or at position 0 (alignment.pairs)
    count in both. A metacool sums its libraries: no network per library."""
    metacool = make_metacool(tmp_path / "map.metacool", dataset["pairs"], dataset["assembly"], binsize)
    from_pairs, from_map = both_networks(dataset, metacool, tmp_path / "raw")
    assert from_map["network"].read_text() == from_pairs["network"].read_text()
    assert from_map["contig_data"].read_text() == from_pairs["contig_data"].read_text()
    # What each data set holds: the contacts that the pairs path left out.
    contig_data, _ = mtn.create_contig_data(str(dataset["assembly"]))
    last_base, position_0 = outside_contacts(dataset["pairs"], contig_data)
    assert (bool(last_base), bool(position_0)) == ((True, False) if dataset["name"] == "synthetic" else (False, True))
    assert sorted(path.name for path in from_map["dir"].glob("network*.txt")) == ["network.txt"]
    assert sorted(path.name for path in from_pairs["dir"].glob("network*.txt")) == ["network.txt", "network_0.txt"]


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
    contig_data, _ = mtn.create_contig_data(str(dataset["assembly"]))
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
    (a row at a time with 1: one read per contig)."""
    pairs, assembly = TESTS_DATA_PAIRS["alignment"]
    metacool = make_metacool(tmp_path / "map.metacool", pairs, assembly)
    contig_data, _ = mtn.create_contig_data(str(assembly))
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
        contig_data, _ = mtn.create_contig_data(str(assembly))
        out = tmp_path / str(chunksize)
        out.mkdir()
        network, contigs = mtn.metacool_to_contacts(
            str(metacool), contig_data, str(out), "network.txt", "contigs.txt", str(out), 1, "empirical_hit",
            self_contacts, chunksize=chunksize
        )  # fmt: skip
        networks.append((Path(network).read_text(), Path(contigs).read_text()))
    assert networks[0] == networks[1]
    assert len(networks[0][0].splitlines()) > 2 * 9973


def test_pixel_row_blocks():
    # Pixels of the rows: 3, 0, 5, 1, 0.
    offsets = np.array([0, 3, 3, 8, 9, 9])
    assert list(mtn.pixel_row_blocks(offsets, 4)) == [(0, 2), (2, 3), (3, 5)]
    assert list(mtn.pixel_row_blocks(offsets, 100)) == [(0, 5)]
    assert list(mtn.pixel_row_blocks(offsets, 1)) == [(0, 1), (1, 2), (2, 3), (3, 5)]
    assert list(mtn.pixel_row_blocks(np.array([0]), 10)) == []


def test_network_of_an_mmetacool_entry(tmp_path):
    """An entry of a .mmetacool collection, given by its URI."""
    pairs, assembly = TESTS_DATA_PAIRS["subset"]
    metacool = make_metacool(tmp_path / "map.metacool", pairs, assembly, 1000)
    collection = tmp_path / "maps.mmetacool"
    metacooler.create_mmetacool(collection, {"lib1": metacool})
    networks = []
    for i, uri in enumerate([str(metacool), f"{collection}::entries/lib1"]):
        contig_data, _ = mtn.create_contig_data(str(assembly))
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
    contig_data, _ = mtn.create_contig_data(str(_write_fasta(tmp_path / "assembly.fa", contigs)))
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
