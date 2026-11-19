# Tests for the contact_map module and the contactmap command.

import re
import shutil
from pathlib import Path

import metator.contact_map as mtcm
import metator.io as mio
import metator.metacool as mtmc
import pandas as pd
import pytest
from tests import synthetic

CONTIG_DATA = "tests_data/contig_data_final.txt"


@pytest.mark.parametrize(
    "metator_object, name, column",
    [
        ("core_bin", "1", "Core_bin_ID"),
        ("core_bin", 1, "Core_bin_ID"),
        ("overlapping_bin", "1", "Overlapping_bin_ID"),
        ("overlapping_bin", "metator_00001_00000", "Overlapping_bin_ID"),
    ],
)
def test_bin_ids_match_zero_padded_ids(metator_object, name, column):
    """Bin ids are written zero-padded ("00001"): "1" must select that bin."""
    data = pd.read_csv(CONTIG_DATA, sep="\t")
    expected = data.loc[pd.to_numeric(data[column], errors="coerce") == 1, "Name"]
    assert len(expected) > 0
    obj = mtcm.MetatorObject(metator_object, name, "unused.fa", CONTIG_DATA, [], min_size=0)
    obj.set_contigs()
    assert obj.contigs == expected.tolist()


def test_contactmap_default_outdir(tmp_path, monkeypatch, run_cli):
    """Without --outdir the map goes to ./contact_map_<name>."""
    assembly = tmp_path / "assembly.fa"
    shutil.copy("tests_data/subset.assembly.fa", assembly)
    contig_data = Path(CONTIG_DATA).resolve()
    pairs = Path("tests_data/subset.pairs.gz").resolve()
    monkeypatch.chdir(tmp_path)
    args = f"-a {assembly} -c {contig_data} -e DpnII -n 1 -O overlapping_bin -s 0 {pairs}"
    run_cli("contactmap", *args.split(" "))
    outdir = tmp_path / "contact_map_1"
    assert (outdir / "1.fa").stat().st_size > 0


@pytest.mark.parametrize(
    "metator_object, name, column, expected",
    [
        ("contig", "NODE_1", "Name", "NODE_1"),
        ("core_bin", 3, "Core_bin_ID", "3"),
        ("overlapping_bin", "MetaTOR_7_0", "Overlapping_bin_ID", "7"),
        ("overlapping_bin", "MetaTOR_7", "Overlapping_bin_ID", "7"),
        ("overlapping_bin", "sample_metator_00007_00002", "Overlapping_bin_ID", "00007"),
        ("recursive_bin", "MetaTOR_7_2", "Recursive_bin_ID", "2"),
        ("recursive_bin", "sample_metator_00007_00002", "Recursive_bin_ID", "00002"),
        ("recursive_bin", "MetaTOR_7_0", "Recursive_bin_ID", "0"),  # logged as an error
        ("final_bin", "MetaTOR_7_2", "Final_bin", "MetaTOR_7_2"),
        ("Other", "12", "Other", 12),
        ("Other", "group_a", "Other", "group_a"),
    ],
)
def test_metator_object_names(metator_object, name, column, expected):
    obj = mtcm.MetatorObject(metator_object, name, "unused.fa", CONTIG_DATA, [])
    assert obj.object == column
    assert obj.name == expected


@pytest.mark.parametrize(
    "metator_object, name",
    [
        ("core_bin", "MetaTOR_1"),
        ("overlapping_bin", "MetaTOR_x_0"),
        ("recursive_bin", "MetaTOR_1_x"),
        # A recursive bin ID alone: recursive bins are numbered within their overlapping bin.
        ("recursive_bin", "4"),
        ("recursive_bin", "MetaTOR_4"),
    ],
)
def test_metator_object_invalid_names(metator_object, name):
    with pytest.raises(ValueError):
        mtcm.MetatorObject(metator_object, name, "unused.fa", CONTIG_DATA, [])


def test_recursive_bin_within_its_overlapping_bin(tmp_path):
    """Recursive bins are numbered within their overlapping bin: a recursive
    bin is the contigs of its recursive bin in its overlapping bin. The
    recursive bins of the same ID of every overlapping bin were selected,
    and a prefixed name (--prefix) failed."""
    table = tmp_path / "contig_data_final.txt"
    pd.DataFrame(
        {
            "ID": [1, 2, 3, 4],
            "Name": ["c1", "c2", "c3", "c4"],
            "Size": [10, 20, 30, 40],
            "Overlapping_bin_ID": ["00001", "00001", "00002", "00002"],
            "Recursive_bin_ID": ["00001", "00002", "00001", "00001"],
        }
    ).to_csv(table, sep="\t", index=False)
    obj = mtcm.MetatorObject("recursive_bin", "my_sample_metator_00002_00001", "unused.fa", str(table), [], min_size=0)
    obj.set_contigs()
    assert (obj.contigs, obj.contigs_size) == (["c3", "c4"], [30, 40])
    obj = mtcm.MetatorObject("overlapping_bin", "my_sample_metator_00001_00000", "unused.fa", str(table), [], min_size=0)
    obj.set_contigs()
    assert obj.contigs == ["c1", "c2"]


def test_metator_object_unknown_column():
    obj = mtcm.MetatorObject("Missing_column", "1", "unused.fa", CONTIG_DATA, [])
    with pytest.raises(ValueError):
        obj.set_contigs()


# Contact maps of a bin of the mock metagenome (tests/synthetic.py).


@pytest.fixture
def final_contig_data(metagenome, metagenome_network, tmp_path):
    """contig_data_final.txt with the true genomes as final bins."""
    data = pd.read_csv(metagenome_network["contig_data"], sep="\t")
    bins = {g.name: f"metator_{k:05d}_00000" for k, g in enumerate(metagenome.genomes, start=1)}
    data["Final_bin"] = data.Name.map(metagenome.truth.Host).map(bins)
    path = tmp_path / "contig_data_final.txt"
    data.to_csv(path, sep="\t", index=False)
    return path, bins


def _pairs_between(metagenome, contigs):
    pairs = metagenome.pairs
    return int((pairs.chr1.isin(contigs) & pairs.chr2.isin(contigs)).sum())


def _matrix_total(out_dir, name, mat_fmt):
    """Sum of the contacts of the matrix written by hicstuff."""
    if mat_fmt == "graal":
        return pd.read_csv(out_dir / f"{name}.mat.tsv", sep="\t").iloc[:, 2].sum()
    if mat_fmt == "bg2":
        return pd.read_csv(out_dir / f"{name}.bg2", sep="\t", header=None).iloc[:, 6].sum()
    import cooler

    return cooler.Cooler(str(out_dir / f"{name}.cool")).pixels()[:]["count"].sum()


@pytest.mark.parametrize("mat_fmt", ["graal", "bg2", "cool"])
def test_generate_contact_map(metagenome, metagenome_inputs, metagenome_indexed_pairs, final_contig_data, tmp_path, mat_fmt):
    contig_data, bins = final_contig_data
    name = bins["genome_3"]
    contigs = [c.name for c in metagenome.molecule("genome_3").contigs]
    out, tmp = tmp_path / "out", tmp_path / "tmp"
    out.mkdir()
    tmp.mkdir()
    n_pairs = mtcm.generate_contact_map(
        str(metagenome_inputs["assembly"]),
        str(contig_data),
        "DpnII",
        name,
        [str(metagenome_indexed_pairs)],
        str(out),
        str(tmp),
        mat_fmt=mat_fmt,
        min_size=0,
    )
    assert n_pairs == _pairs_between(metagenome, contigs)
    assert sorted(synthetic.read_fasta(out / f"{name}.fa")) == sorted(contigs)
    assert _matrix_total(out, name, mat_fmt) == n_pairs


def test_contact_map_min_size(metagenome, metagenome_inputs, metagenome_indexed_pairs, final_contig_data, tmp_path):
    """Contigs below --min-size are left out of the map."""
    contig_data, bins = final_contig_data
    large = [c.name for c in metagenome.molecule("genome_3").contigs if c.length >= 15_000]
    assert 0 < len(large) < len(metagenome.molecule("genome_3").contigs)
    out, tmp = tmp_path / "out", tmp_path / "tmp"
    out.mkdir()
    tmp.mkdir()
    n_pairs = mtcm.generate_contact_map(
        str(metagenome_inputs["assembly"]),
        str(contig_data),
        "DpnII",
        bins["genome_3"],
        [str(metagenome_indexed_pairs)],
        str(out),
        str(tmp),
        min_size=15_000,
    )
    assert n_pairs == _pairs_between(metagenome, large)
    assert sorted(synthetic.read_fasta(out / f"{bins['genome_3']}.fa")) == sorted(large)
    assert _matrix_total(out, bins["genome_3"], "graal") == n_pairs


def test_contact_map_without_pairs(metagenome, metagenome_inputs, metagenome_indexed_pairs, final_contig_data, tmp_path):
    """A contig smaller than --min-size: nothing to extract, no map."""
    contig_data, _ = final_contig_data
    small = min(metagenome.contigs, key=lambda c: c.length)
    n_pairs = mtcm.generate_contact_map(
        str(metagenome_inputs["assembly"]),
        str(contig_data),
        "DpnII",
        small.name,
        [str(metagenome_indexed_pairs)],
        str(tmp_path),
        str(tmp_path),
        metator_object="contig",
        min_size=small.length + 1,
    )
    assert n_pairs == 0
    assert not list(tmp_path.glob("*.mat.tsv"))


def test_metator_object_write_fasta(tmp_path, monkeypatch):
    """The fasta of an object, the records of its contigs, is written
    without the pyfastx command, which is not on the PATH here."""
    monkeypatch.setenv("PATH", str(tmp_path))
    # A private copy of the assembly: pyfastx indexes it (.fxi) in place.
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile("tests_data/assembly.fa", assembly)
    obj = mtcm.MetatorObject("contig", "NODE_522", str(assembly), CONTIG_DATA, [])
    obj.contigs = ["NODE_1404", "NODE_522"]
    obj.write_fasta(str(tmp_path))
    names = [line[1:] for line in Path(obj.fasta).read_text().splitlines() if line.startswith(">")]
    assert names == ["NODE_1404", "NODE_522"]


# The contact map of a bin from a .metacool contact map, in place of the pairs
# files: that of its contigs, as metacooler extract writes it.


@pytest.fixture(scope="module")
def edge_map(metagenome, metagenome_files, tmp_path_factory):
    """The pairs of the mock metagenome with four pairs on the edges of two
    contigs of genome_3, a and b (in the upper triangle, as metator.align
    writes them): a read at position 0 within a, a read on the last base
    of a with one on the first base of b, and a read past the end of a,
    within a and with b. Sorted and indexed for the pairs path, and the
    .metacool map that metacooler builds from them, one bin per contig,
    every contig in no genome yet. Returns both, and the names of a and b."""
    import metacooler

    directory = tmp_path_factory.mktemp("edge_map")
    rank = {contig.name: i for i, contig in enumerate(metagenome.contigs)}
    a, b = sorted((contig for contig in metagenome.molecule("genome_3").contigs), key=lambda contig: rank[contig.name])[:2]
    pairs = directory / "edges.pairs"
    shutil.copyfile(metagenome_files["pairs"], pairs)
    with open(pairs, "a") as handle:
        for read, pos1, contig2, pos2 in [
            ("edge_1", 0, a, 2000),
            ("edge_2", a.length, b, 1),
            ("edge_3", 100, a, a.length + 1),
            ("edge_4", a.length + 50, b, 7),
        ]:
            handle.write(f"{read}\t{a.name}\t{pos1}\t{contig2.name}\t{pos2}\t+\t-\n")
    indexed = mio.sort_and_index_pairs(str(pairs), tmp_dir=str(directory))
    truth = metagenome.truth
    metacool = directory / "map.metacool"
    bins = pd.DataFrame({"chrom": truth.Name, "start": 0, "end": truth.Length})
    metacooler.from_pairs(str(metacool), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    return indexed, metacool, (a.name, b.name)


def _by_contig_pair(contig1, contig2, count):
    """The contacts summed by pair of contigs, in either order."""
    table = pd.DataFrame({"contig1": list(contig1), "contig2": list(contig2), "count": list(count)})
    pairs = table[["contig1", "contig2"]].apply(sorted, axis=1, result_type="expand")
    return table.groupby([pairs[0], pairs[1]])["count"].sum().astype("int64").to_dict()


def _cool_by_contig_pair(path):
    """The contacts of a .cool file, summed by pair of contigs."""
    import cooler

    clr = cooler.Cooler(str(path))
    chrom = clr.bins()[:]["chrom"].astype(str).to_numpy()
    pixels = clr.pixels()[:]
    return _by_contig_pair(chrom[pixels["bin1_id"]], chrom[pixels["bin2_id"]], pixels["count"])


@pytest.mark.parametrize("min_size", [0, 15_000])
def test_contact_map_from_a_map_is_that_of_its_pairs(
    min_size, metagenome, metagenome_inputs, final_contig_data, edge_map, tmp_path
):
    """The contact map of a bin from a .metacool map (extract_contact_map:
    the .cool that metacooler extract writes, one bin per contig) is the
    map of its restriction fragments that hicstuff builds from the pairs
    files that the map was made of (generate_contact_map, written as a 2D
    bedgraph), summed by pair of contigs: the same contigs (all of them, or
    those of at least --min-size bp), the same contacts between each two of
    them and within each, and the same FASTA file of its contigs.

    What is compared: the pairs that the map holds, those whose two reads
    lie on the contigs, within their length. The pairs files hold four
    pairs on the edges of two contigs of the bin (see edge_map): a read at
    position 0, which metacooler puts on the first base and the pairs path
    queries from 0, and one on the last base count in both; a read past the
    end of its contig counts in neither, as metacooler leaves it out and
    the pairs path queries a contig up to its length. The pairs are in the
    upper triangle, as metator.align writes them: the pairs path reads one
    order of each two contigs, where the map holds both."""
    contig_data, bins = final_contig_data
    name = bins["genome_3"]
    contigs = [contig.name for contig in metagenome.molecule("genome_3").contigs if contig.length >= min_size]
    indexed, metacool, (a, b) = edge_map
    from_pairs, from_map, tmp = tmp_path / "pairs", tmp_path / "map", tmp_path / "tmp"
    for directory in (from_pairs, from_map, tmp):
        directory.mkdir()
    assembly = str(metagenome_inputs["assembly"])
    n_pairs = mtcm.generate_contact_map(assembly, str(contig_data), "DpnII", name, [str(indexed)], str(from_pairs), str(tmp),
                                        mat_fmt="bg2", min_size=min_size)  # fmt: skip
    written = mtcm.extract_contact_map(assembly, str(contig_data), name, str(metacool), str(from_map), min_size=min_size)
    fragments = pd.read_csv(from_pairs / f"{name}.bg2", sep="\t", header=None, usecols=[0, 3, 6])
    expected = _by_contig_pair(fragments[0], fragments[3], fragments[6])
    assert _cool_by_contig_pair(from_map / f"{name}.cool") == expected
    # The pairs on the edges of a and b that both count.
    edges = (a in contigs) + (a in contigs and b in contigs)
    assert n_pairs == sum(expected.values()) == _pairs_between(metagenome, contigs) + edges
    assert sorted(written) == sorted(contigs) and {contig for pair in expected for contig in pair} == set(contigs)
    assert (len(contigs) < len(metagenome.molecule("genome_3").contigs)) == (min_size > 0)
    assert synthetic.read_fasta(from_map / f"{name}.fa") == synthetic.read_fasta(from_pairs / f"{name}.fa")


@pytest.mark.parametrize("given", ["a map", "an entry of a collection"])
def test_contactmap_command_from_a_map(given, metagenome, metagenome_inputs, final_contig_data, edge_map, tmp_path, run_cli):
    """metator contactmap MAP, the map in place of the pairs files (a file,
    or an entry of a .mmetacool collection), writes the contact map of the
    bin as metacooler extract writes it, <name>.cool (--mat-fmt cool, the
    default then), at one bin per contig, and the FASTA file of its
    contigs; its log says what it wrote. --enzyme is not used. A second run
    is refused while <name>.cool is there, before any work (both files left
    as they were), and replaces it with -F."""
    import cooler
    import metacooler

    contig_data, bins = final_contig_data
    name = bins["genome_3"]
    _, metacool, _ = edge_map
    if given == "an entry of a collection":
        metacooler.create_mmetacool(tmp_path / "maps.mmetacool", {"lib1": metacool})
        metacool = f"{tmp_path / 'maps.mmetacool'}::entries/lib1"
    out = tmp_path / "out"
    args = ["-a", metagenome_inputs["assembly"], "-c", contig_data, "-e", "not_an_enzyme", "-n", name, "-s", "0", "-o", out]
    run_cli("contactmap", *args, metacool)
    cool = out / f"{name}.cool"
    contigs = sorted(contig.name for contig in metagenome.molecule("genome_3").contigs)
    clr = cooler.Cooler(str(cool))
    assert sorted(clr.chromnames) == contigs and len(clr.bins()[:]) == len(contigs)
    assert sorted(path.name for path in out.iterdir() if not path.name.endswith(".log")) == [f"{name}.cool", f"{name}.fa"]
    (log,) = out.glob("metator_contactmap_*.log")
    logged = log.read_text()
    assert f"INFO :: Generate HiC contact map for {name} from {metacool}" in logged
    assert f"INFO :: Created {cool}: bins={len(contigs)}, contigs={len(contigs)}. NNZ={clr.info['nnz']}" in logged
    first = _cool_by_contig_pair(cool)
    before = {path: path.stat().st_mtime_ns for path in out.iterdir() if not path.name.endswith(".log")}
    message = f"{cool} exists: give -F to replace it with the contact map of {name}."
    with pytest.raises(FileExistsError, match=f"^{re.escape(message)}$"):
        run_cli("contactmap", *args, metacool)
    assert {path: path.stat().st_mtime_ns for path in before} == before
    # In the log of the second run, which is another file in another second.
    assert f"ERROR :: {message}" in "".join(path.read_text() for path in out.glob("metator_contactmap_*.log"))
    run_cli("contactmap", *args, "-F", "-m", "cool", metacool)
    assert cool.stat().st_mtime_ns != before[cool] and _cool_by_contig_pair(cool) == first


@pytest.mark.parametrize("lacking", ["a contig", "every contig"])
def test_contact_map_from_a_map_lacking_contigs(
    lacking, metagenome, metagenome_files, metagenome_inputs, final_contig_data, tmp_path, logged
):
    """A contig of the bin that the map lacks (a map made with a contig
    filter, for example) is left out of its contact map, with a warning
    that names it; with none of them, no contact map is written."""
    import cooler
    import metacooler

    contig_data, bins = final_contig_data
    name = bins["genome_3"]
    contigs = [contig.name for contig in metagenome.molecule("genome_3").contigs]
    left_out = contigs[:1] if lacking == "a contig" else contigs
    truth = metagenome.truth[~metagenome.truth.Name.isin(left_out)]
    metacool = tmp_path / "map.metacool"
    map_bins = pd.DataFrame({"chrom": truth.Name, "start": 0, "end": truth.Length})
    metacooler.from_pairs(str(metacool), pairs=str(metagenome_files["pairs"]), bins=map_bins, contig_genome_map=None, ps=False)
    written = mtcm.extract_contact_map(
        str(metagenome_inputs["assembly"]), str(contig_data), name, str(metacool), str(tmp_path), min_size=0
    )
    # Named in the order of the contig data.
    left_out = [contig for contig in pd.read_csv(contig_data, sep="\t").Name if contig in left_out]
    shown = ", ".join(repr(contig) for contig in left_out[:5]) + (
        f", ... (+{len(left_out) - 5} more)" if len(left_out) > 5 else ""
    )
    warning = (
        "WARNING",
        f"{len(left_out)} of the contigs of {name} are not in {metacool}: they have no contacts there, as contigs that no "
        f"read maps to: {shown}.",
    )
    assert warning in logged()
    assert sorted(written) == sorted(set(contigs) - set(left_out))
    if written:
        assert sorted(cooler.Cooler(str(tmp_path / f"{name}.cool")).chromnames) == sorted(written)
    else:
        assert ("INFO", f"{metacool} holds no contig of {name} of at least 0 bp: no contact map is written.") in logged()
        assert not (tmp_path / f"{name}.cool").exists()


# What contactmap refuses of a map in place of the pairs files, the reason
# logged as an error, and no contact map written.
CONTACTMAP_MAP_ERRORS = {
    "--filter": (
        ["-f", "{map}"],
        "metator contactmap from the contact map {map} takes no -f/--filter: it writes the contact map of the contigs as "
        "metacooler extract writes it, a .cool file, from a map that holds no read positions, which the filters of "
        "hicstuff and its maps of restriction fragments need.",
    ),
    "--pcr-dup": (["-D", "{map}"], "metator contactmap from the contact map {map} takes no -D/--pcr-dup: "),
    "--mat-fmt graal": (["-m", "graal", "{map}"], "metator contactmap from the contact map {map} takes no -m/--mat-fmt graal: "),
    "--filter and --mat-fmt bg2": (
        ["-f", "-m", "bg2", "{map}"],
        "metator contactmap from the contact map {map} takes no -f/--filter or -m/--mat-fmt bg2: ",
    ),
    "a map and pairs files": (
        ["{map}", "{pairs}"],
        "metator contactmap reads pairs files, or one .metacool map: got 2 files, 1 of them .metacool maps.",
    ),
    "two maps": (
        ["{map}", "{map}"],
        "metator contactmap reads pairs files, or one .metacool map: got 2 files, 2 of them .metacool maps.",
    ),
    "other lengths": (
        ["{shorter}"],
        "{shorter} is not a contact map of the contigs of {contig_data}: 1 contig(s) have another length in it: "
        "'{contig}' ({length} bp in the map, {size} bp in {contig_data}).",
    ),
    "both triangles": (["{square}"], "{square} stores its contacts with storage-mode 'square': "),
}


@pytest.mark.parametrize("case", list(CONTACTMAP_MAP_ERRORS))
def test_contactmap_command_from_a_map_refused(case, metagenome, metagenome_files, metagenome_inputs, metagenome_indexed_pairs,
                                               final_contig_data, edge_map, tmp_path, run_cli):  # fmt: skip
    """contactmap refuses, with a map, the options that need the pairs files
    (--filter, --pcr-dup, a --mat-fmt other than cool), a map with pairs
    files or another map, and a map that cannot give the contacts of the
    contigs of the bin (one of another length, or both triangles stored):
    the reason logged as an error in the log file, as MetacoolError, and
    nothing written, neither the contact map nor the FASTA file of the
    contigs (written before the map was checked)."""
    import h5py
    import metacooler

    contig_data, bins = final_contig_data
    name = bins["genome_3"]
    _, metacool, _ = edge_map
    contig = metagenome.molecule("genome_3").contigs[0]
    names = {
        "map": metacool,
        "pairs": metagenome_indexed_pairs,
        "shorter": tmp_path / "shorter.metacool",
        "square": tmp_path / "square.metacool",
        "contig_data": contig_data,
        "contig": contig.name,
        "length": contig.length - 1000,
        "size": contig.length,
    }
    shutil.copyfile(metacool, names["square"])
    with h5py.File(names["square"], "r+") as group:
        group.attrs["storage-mode"] = "square"
    truth = metagenome.truth
    map_bins = pd.DataFrame(
        {"chrom": truth.Name, "start": 0, "end": truth.Length.where(truth.Name != contig.name, names["length"])}
    )
    metacooler.from_pairs(
        str(names["shorter"]), pairs=str(metagenome_files["pairs"]), bins=map_bins, contig_genome_map=None, ps=False
    )
    inputs, message = CONTACTMAP_MAP_ERRORS[case]
    message = message.format(**names)
    out = tmp_path / "out"
    args = ["-a", metagenome_inputs["assembly"], "-c", contig_data, "-e", "DpnII", "-n", name, "-s", "0", "-o", out]
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}"):
        run_cli("contactmap", *args, *[value.format(**names) for value in inputs])
    (log,) = out.glob("metator_contactmap_*.log")
    assert f"ERROR :: {message}" in log.read_text()
    assert [path.name for path in out.iterdir()] == [log.name]


@pytest.mark.parametrize("case", ["other lengths", "both triangles"])
def test_contactmap_command_from_a_refused_map_keeps_an_earlier_run(
    case, metagenome, metagenome_files, metagenome_inputs, final_contig_data, edge_map, tmp_path, run_cli
):
    """A map refused with -F (one of another length, or both triangles
    stored), after a run that wrote the contact map of the bin and the
    FASTA file of its contigs, leaves both as they were: the map is checked
    before an output is written or replaced. The contact map was removed,
    and the FASTA file written again, before the map was checked."""
    import h5py
    import metacooler

    contig_data, bins = final_contig_data
    name = bins["genome_3"]
    _, metacool, _ = edge_map
    out = tmp_path / "out"
    args = ["-a", metagenome_inputs["assembly"], "-c", contig_data, "-e", "DpnII", "-n", name, "-s", "0", "-o", out]
    run_cli("contactmap", *args, metacool)
    refused = tmp_path / "refused.metacool"
    if case == "both triangles":
        shutil.copyfile(metacool, refused)
        with h5py.File(refused, "r+") as group:
            group.attrs["storage-mode"] = "square"
    else:
        contig = metagenome.molecule("genome_3").contigs[0]
        truth = metagenome.truth
        map_bins = pd.DataFrame(
            {"chrom": truth.Name, "start": 0, "end": truth.Length.where(truth.Name != contig.name, contig.length - 1000)}
        )
        metacooler.from_pairs(str(refused), pairs=str(metagenome_files["pairs"]), bins=map_bins, contig_genome_map=None, ps=False)
    written = {path.name: path.read_bytes() for path in out.iterdir() if not path.name.endswith(".log")}
    assert sorted(written) == [f"{name}.cool", f"{name}.fa"]
    with pytest.raises(mtmc.MetacoolError):
        run_cli("contactmap", *args, "-F", refused)
    assert {path.name: path.read_bytes() for path in out.iterdir() if not path.name.endswith(".log")} == written
