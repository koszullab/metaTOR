# Test for CLI toos of metator
# Commands are simply run to test crashes.


from pathlib import Path
import pandas as pd
import metator.cli as mtc
import metator.metacool as mtmc
import metator.network as mtn
import metacooler
import pytest
from tests.maps import contig_map
import re
import shutil
import os


# Use global variables for temporary tmp/out folders
@pytest.fixture(scope="session")
def tmp_dir(tmp_path_factory):
    p = tmp_path_factory.mktemp("tmp")
    os.makedirs(p, exist_ok=True)
    return str(p)


global_args = {
    "DEPTH": "tests_data/depth.txt",
    "FASTA": "tests_data/assembly.fa",
    "PAIRS": "tests_data/alignment.pairs",
    "SUBSET_FASTA": "tests_data/subset.assembly.fa",
    "SUBSET_PAIRS": "tests_data/subset.pairs.gz",
}


@pytest.fixture(autouse=True)
def private_assemblies(private_link, monkeypatch):
    """Each test reads the assemblies of tests_data through links of its own (see
    private_link): two test sessions sharing tests_data removed each other's
    pyfastx index (FileNotFoundError), when metaTOR read them with pyfastx."""
    for key in ["FASTA", "SUBSET_FASTA"]:
        monkeypatch.setitem(global_args, key, private_link(global_args[key]))


@pytest.fixture(scope="session")
def tests_data_maps(tmp_path_factory):
    """The maps of the contigs of tests_data, in no cluster yet, as hicstuff or
    metacooler make them from the pairs (see tests.maps.contig_map): "full", of
    alignment.pairs on assembly.fa, and "subset", of subset.pairs.gz on
    subset.assembly.fa."""
    directory = tmp_path_factory.mktemp("tests_data_maps")
    maps = {}
    for name, assembly, pairs in [("full", "FASTA", "PAIRS"), ("subset", "SUBSET_FASTA", "SUBSET_PAIRS")]:
        contig_data = mtn.create_contig_data(global_args[assembly])
        maps[name] = contig_map(global_args[pairs], contig_data, directory / f"{name}.metacool", columns=False)
    return maps


@pytest.fixture(scope="session")
def tests_data_network_maps(tests_data_maps, tmp_path_factory):
    """The maps of tests_data (see tests_data_maps) with their network stored in
    them, as metator network -e HindIII,DpnII stores it (empirical_hit): what
    partition and validation read."""
    directory = tmp_path_factory.mktemp("tests_data_network_maps")
    maps = {}
    for name, assembly in [("full", "FASTA"), ("subset", "SUBSET_FASTA")]:
        path = Path(shutil.copyfile(tests_data_maps[name], directory / f"{name}.metacool"))
        contig_data = mtn.create_contig_data(global_args[assembly], None, "HindIII,DpnII")
        contacts = mtn.contig_contacts(str(path))
        names = metacooler.MetaCooler(str(path)).contigs["name"]
        table = mtn.contig_data_table(contig_data, names, mtn.contig_hits(contacts, len(names)))
        mtmc.write_map(str(path), str(path), network=mtmc.Network(table, "empirical_hit"), force=True)
        maps[name] = path
    return maps


@pytest.fixture
def a_network_map(tests_data_network_maps, tmp_path):
    """A copy of a map of tests_data with its network (see
    tests_data_network_maps), its own: runs update their map in place."""

    def copy(name="subset"):
        path = tmp_path / "network_maps" / f"{name}.metacool"
        path.parent.mkdir(exist_ok=True)
        return shutil.copyfile(tests_data_network_maps[name], path)

    return copy


@pytest.fixture
def a_map(tests_data_maps, tmp_path):
    """A copy of a map of tests_data (see tests_data_maps), its own: runs update
    their map in place."""

    def copy(name="subset"):
        path = tmp_path / "maps" / f"{name}.metacool"
        path.parent.mkdir(exist_ok=True)
        return shutil.copyfile(tests_data_maps[name], path)

    return copy


NORMALIZE = (
    "norm",
    ["None", "length", "abundance", "RS", "empirical_hit", "theoritical_hit"],
)
ALGORITHM = ("alg", ["louvain", "leiden"])


@pytest.mark.parametrize(*NORMALIZE)
def test_network_normalizations(norm, a_map, tmp_path, run_cli):
    args = "{MAP} -a {FASTA} -d {DEPTH} -o {OUT_TEST} -T {TMP} -n {0} -e HindIII,DpnII".format(
        norm,
        MAP=a_map("full"),
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    run_cli("network", *args.split(" "))


@pytest.mark.parametrize(*ALGORITHM)
def test_partition(alg, a_network_map, tmp_path, run_cli):
    args = ("{MAP} -i 5 -o {OUT_TEST} -t 8 -T {TMP} -A {0} -F --no-clean-up").format(
        alg,
        MAP=a_network_map(),
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator partition", *args.split(" "))
    run_cli("partition", *args.split(" "))
    # The partition is in the map, and in no table or FASTA file.
    assert sorted(path.name for path in Path(tmp_path, "out_test").iterdir() if path.suffix != ".log") == []


def _partitioned(metacool, tmp_path, run_cli):
    """The map metacool, in which metator partition stored its partition (two
    iterations), which validation reads."""
    run_cli("partition", metacool, "-i", "2", "-o", Path(tmp_path, "partition"), "-T", Path(tmp_path, "tmp"))
    return metacool


def test_validation(a_network_map, tmp_path, run_cli):
    """From the partition that partition MAP stores in the map of the subset: its
    overlapping clusters of 30 kb or more, written in the temporary directory
    of the run (kept with --no-clean-up), for miComplete."""
    args = ("{MAP} -a {SUBSET_FASTA} -i 5 -s 30000 -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up").format(
        MAP=_partitioned(a_network_map(), tmp_path, run_cli),
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator validation", *args.split(" "))
    run_cli("validation", *args.split(" "))
    assert list(Path(tmp_path, "tmp").glob("*/overlapping_clusters/*.fa"))
    assert list(Path(tmp_path, "tmp").glob("*/contig_data_partition.txt"))
    assert not Path(tmp_path, "out_test", "overlapping_clusters").exists()


def test_network_partition_validation(tmp_dir, a_map, run_cli):
    tmp = Path(tmp_dir, "tmp")
    out = Path(tmp_dir, "out_test")
    metacool = a_map()
    args = "{MAP} -a {SUBSET_FASTA} -o {OUT_TEST} -T {TMP} -e HindIII,DpnII --no-clean-up".format(
        MAP=metacool,
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator network", *args.split(" "))
    run_cli("network", *args.split(" "))

    args = ("{MAP} -i 5 -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up").format(
        MAP=metacool,
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator partition", *args.split(" "))
    run_cli("partition", *args.split(" "))

    args = ("{MAP} -a {SUBSET_FASTA} -s 30000 -i 5 -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up").format(
        MAP=metacool,
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator validation", *args.split(" "))
    run_cli("validation", *args.split(" "))
    # The network stored in the map by network, which the map of the
    # clusters keeps: no table of it.
    assert mtn.network_normalization(str(metacool)) == "empirical_hit"
    assert not list(out.glob("*network*.txt"))


def test_pipeline(a_map, tmp_path, run_cli):
    """Test the metator pipeline command with subset data."""
    metacool = a_map()
    args = ("{MAP} -a {SUBSET_FASTA} -F -o {OUT_TEST} -s 30000 -e HindIII,DpnII --no-clean-up -T {TMP}").format(
        MAP=metacool,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\n", "metator pipeline", *args.split(" "))
    run_cli("pipeline", *args.split(" "))
    # The contact map of the contigs, grouped by cluster: the map of the run,
    # updated in place.
    metacooler.validate(str(metacool))
    genomes = metacooler.MetaCooler(str(metacool)).genomes
    # The tables of the validation, in the temporary directory of the run,
    # which --no-clean-up keeps; none in the output directory.
    (summary,) = Path(tmp_path, "tmp").glob("*/cluster_summary.txt")
    assert set(genomes.name) - {"__unbinned__"} == {line.split("\t")[0] for line in summary.read_text().splitlines()[1:]}
    assert not list(Path(tmp_path, "out_test").glob("*.txt"))
    # The contig columns but Shotgun_coverage, which is - without a depth
    # file.
    contigs = metacooler.MetaCooler(str(metacool)).contigs
    assert [column for column in contigs.columns if not column.startswith("qc_")] == [
        "name", "length", "genome_id", "GC_content", "Hit", "Restriction_site"]  # fmt: skip


def test_network_from_a_metacool(a_map, tmp_path, run_cli):
    """metator network MAP, on a map of the pairs of tests_data at 1 kb: the
    network and contig data of the map with one bin per contig, without a
    network per library."""
    # metacooler indexes the assembly (.fai): a private copy.
    assembly = shutil.copyfile(global_args["SUBSET_FASTA"], tmp_path / "assembly.fa")
    metacool = tmp_path / "subset.metacool"
    metacooler.from_pairs(str(metacool), pairs=global_args["SUBSET_PAIRS"], fasta=str(assembly), binsize=1000,
                          contig_genome_map=None, contig_ordering="fasta", ps=False)  # fmt: skip
    networks, tables = {}, {}
    for binned, given in [("contig", a_map()), ("1kb", metacool)]:
        out = tmp_path / binned
        run_cli("network", given, "-a", assembly, "-e", "HindIII,DpnII", "-o", out, "-T", tmp_path / "tmp")
        names = [str(name) for name in metacooler.MetaCooler(str(given)).contigs["name"]]
        network = mtn.map_network(str(given))
        networks[binned] = {
            tuple(sorted((names[a - 1], names[b - 1]))): weight for a, b, weight in network.itertuples(index=False)
        }
        tables[binned] = mtn.map_contig_table(str(given)).drop(columns="ID").set_index("Name").sort_index()
        assert not list(out.glob("*network*.txt"))
    # The same network and contig data, by contig name: the maps hold their
    # contigs in other orders.
    assert networks["1kb"].keys() == networks["contig"].keys()
    for pair, weight in networks["1kb"].items():
        assert weight == pytest.approx(networks["contig"][pair], rel=1e-12)
    pd.testing.assert_frame_equal(tables["1kb"], tables["contig"])
    (log,) = (tmp_path / "1kb").glob("metator_network_*.log")
    assert f"INFO :: Information of {metacool}:" in log.read_text()


def test_network_from_the_columns_of_a_metacool(tmp_path, run_cli):
    """metator network MAP, from a map of the contigs with the columns that
    metaTOR measures with a depth file and an enzyme: without --depth and
    --enzyme, its Restriction_site and Shotgun_coverage columns give the
    network and the contig data that they give (theoritical_hit needs both),
    byte for byte."""
    contig_data = mtn.create_contig_data(global_args["FASTA"], global_args["DEPTH"], "HindIII,DpnII")
    metacool = contig_map(global_args["PAIRS"], contig_data, tmp_path / "contigs.metacool")
    maps = {how: shutil.copyfile(metacool, tmp_path / f"{how}.metacool") for how in ["options", "map"]}
    common = ["-a", global_args["FASTA"], "-n", "theoritical_hit", "-T", tmp_path / "tmp"]
    run_cli("network", maps["options"], *common, "-d", global_args["DEPTH"], "-e", "HindIII,DpnII", "-o", tmp_path / "options")
    run_cli("network", metacool, *common, "-o", tmp_path / "map", "--output-map", maps["map"], "-F")
    pd.testing.assert_frame_equal(mtn.map_network(str(maps["map"])), mtn.map_network(str(maps["options"])))
    pd.testing.assert_frame_equal(mtn.map_contig_table(str(maps["map"])), mtn.map_contig_table(str(maps["options"])))
    for path in maps.values():
        stored = metacooler.MetaCooler(str(path)).metadata["metator"]["network"]
        assert stored["read"] == {"contigs": ["Restriction_site", "Shotgun_coverage"] if path == maps["map"] else []}
    (log,) = (tmp_path / "map").glob("metator_network_*.log")
    logged = log.read_text()
    for what, column in [("restriction sites", "Restriction_site"), ("shotgun coverage", "Shotgun_coverage")]:
        assert f"INFO :: The {what} of the contigs: the {column} column of {metacool}." in logged


@pytest.mark.parametrize("data", ["full", "subset"])
def test_pipeline_network_is_that_of_network(data, a_map, network_from_a_map, tmp_path, run_cli):
    """A pipeline run reads its network from its map, as metator network does: the
    network and the contig data that its partition gets are those that network
    MAP stores and the partition of a run step by step reads (see
    metator.network.map_network and map_contig_table), on the map of the pairs
    of tests_data (alignment.pairs holds reads at position 0) and on the map of
    the subset."""
    assembly = global_args["FASTA" if data == "full" else "SUBSET_FASTA"]
    depth = ["-d", global_args["DEPTH"]] if data == "full" else []
    common = ["-a", assembly, *depth, "-e", "HindIII,DpnII", "-T", tmp_path / "tmp"]
    metacool = a_map(data)
    stored = tmp_path / "network.metacool"
    run_cli("network", metacool, *common, "-o", tmp_path / "network", "--output-map", stored)
    with pytest.raises(RuntimeError, match=r"^network built$"):
        run_cli("pipeline", metacool, *common, "-o", tmp_path / "pipeline")
    # Read by metator network, then by pipeline.
    assert network_from_a_map == [str(metacool), str(metacool)]
    contig_data, network = network_from_a_map.partitioned
    pd.testing.assert_frame_equal(network, mtn.map_network(str(stored)))
    pd.testing.assert_frame_equal(contig_data, mtn.map_contig_table(str(stored)))
    assert not list((tmp_path / "pipeline").glob("*network*.txt"))


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_metacool_of_other_contigs(command, tmp_path, run_cli):
    """A map that cannot give the network of the assembly stops the run with
    its MetacoolError, the reason logged as an error in the log file."""
    other = tmp_path / "other.fa"
    other.write_text(">a\n" + "ACGT" * 100 + "\n")
    pairs = tmp_path / "other.pairs"
    pairs.write_text("## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\nr1\ta\t10\ta\t300\t+\t-\n")
    metacool = tmp_path / "other.metacool"
    metacooler.from_pairs(str(metacool), pairs=str(pairs), fasta=str(other), binsize=100, contig_genome_map=None,
                          contig_ordering="fasta", ps=False)  # fmt: skip
    assembly = shutil.copyfile(global_args["SUBSET_FASTA"], tmp_path / "assembly.fa")
    out = tmp_path / "out"
    message = f"{metacool} is not a contact map of the contigs of the assembly: "
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}") as error:
        run_cli(command, metacool, "-a", assembly, "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: {error.value}\n" in log.read_text()
    assert not (out / "network.txt").exists()


def test_tmpdir_contents_survive(a_map, tmp_path, run_cli):
    """--tmpdir may already hold the user's files: a run removes only its own."""
    user_tmp = tmp_path / "scratch"
    user_tmp.mkdir()
    keep = user_tmp / "precious.txt"
    keep.write_text("keep me")
    args = [
        str(a_map()),
        "-a", str(Path(global_args["SUBSET_FASTA"]).resolve()),
        "-o", str(tmp_path / "out"),
        "-T", str(user_tmp),
    ]  # fmt: skip
    run_cli("network", *args)
    assert keep.read_text() == "keep me"
    assert list(user_tmp.iterdir()) == [keep]


def test_default_tmpdir_is_tidied(a_map, tmp_path, monkeypatch, run_cli):
    """Without -T, ./tmp is removed if the run created it, kept if it existed."""
    args = [
        str(a_map()),
        "-a", str(Path(global_args["SUBSET_FASTA"]).resolve()),
        "-o", "out",
    ]  # fmt: skip
    monkeypatch.chdir(tmp_path)
    run_cli("network", *args)
    assert not (tmp_path / "tmp").exists()

    (tmp_path / "tmp").mkdir()
    (tmp_path / "tmp" / "mine.txt").write_text("x")
    run_cli("network", *args)
    assert [p.name for p in (tmp_path / "tmp").iterdir()] == ["mine.txt"]


def test_network_default_outdir(a_map, tmp_path, monkeypatch, run_cli):
    """Without -o, the log of the run is written to the current directory: the
    network goes to the map, and no table of it anywhere."""
    args = [
        str(a_map()),
        "-a", str(Path(global_args["SUBSET_FASTA"]).resolve()),
    ]  # fmt: skip
    monkeypatch.chdir(tmp_path)
    run_cli("network", *args)
    assert mtn.network_normalization(args[0]) == "empirical_hit"
    assert not list(tmp_path.rglob("*network*.txt"))
    (log,) = tmp_path.glob("metator_network_*.log")
    assert f"## Command used: metator network {args[0]} -a {args[2]}\n" in log.read_text()


def test_invalid_prefix_is_rejected(a_network_map, tmp_path, run_cli):
    """miComplete drops dots from bin names, so a dotted --prefix is refused."""
    args = ("{MAP} -o {OUT} -T {TMP} -p my.sample").format(
        MAP=a_network_map(), OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError, match="Invalid --prefix"):
        run_cli("partition", *args.split(" "))


def test_spinglass_is_rejected_before_any_work(tmp_path, monkeypatch, run_cli):
    """spinglass is gone: refuse it up front, not after miComplete has run."""

    def not_reached(*args, **kwargs):
        raise AssertionError("recursive decontamination should not start")

    monkeypatch.setattr(mtc.mtv, "recursive_decontamination", not_reached)
    args = ("map.metacool -a {SUBSET_FASTA} -o {OUT} -T {TMP} -A spinglass").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    result = run_cli("validation", *args.split(" "), exit_code=2)
    assert "'spinglass' is not one of 'louvain', 'leiden'" in result.output
    assert not (tmp_path / "out").exists()
    assert not (tmp_path / "tmp").exists()


# Argument checks: each command refuses inconsistent options before working.
# The map of the subset has no Restriction_site or Shotgun_coverage column.

NETWORK_ERRORS = {
    "RS without enzyme": ["-n", "RS"],
    "theoritical_hit without enzyme": ["-n", "theoritical_hit", "-d", "depth.txt"],
    "abundance without depth": ["-n", "abundance"],
    "theoritical_hit without depth": ["-n", "theoritical_hit", "-e", "DpnII"],
    "assembly not a fasta": ["-a", "missing.fa"],
}


@pytest.mark.parametrize("case", list(NETWORK_ERRORS))
def test_network_argument_errors(case, a_map, tmp_path, run_cli):
    args = [a_map(), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), *NETWORK_ERRORS[case]]
    if "-a" not in args:
        args += ["-a", global_args["SUBSET_FASTA"]]
    with pytest.raises(ValueError):
        run_cli("network", *args)


# Unknown values are refused while the command line is read: nothing runs.
NETWORK_BAD_CHOICES = {
    "unknown normalization": ["-n", "cosine"],
}


@pytest.mark.parametrize("case", list(NETWORK_BAD_CHOICES))
def test_network_bad_choices(case, tmp_path, run_cli):
    args = ["map.metacool", "-a", global_args["FASTA"], "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"),
            *NETWORK_BAD_CHOICES[case]]  # fmt: skip
    result = run_cli("network", *args, exit_code=2)
    assert f"'{NETWORK_BAD_CHOICES[case][-1]}' is not one of" in result.output
    assert not (tmp_path / "out").exists()
    assert not (tmp_path / "tmp").exists()


PIPELINE_ERRORS = {
    "RS without enzyme": ["-m", "RS"],
    "abundance without depth": ["-m", "abundance"],
    "assembly not a fasta": ["-a", "missing.fa"],
}


@pytest.mark.parametrize("case", list(PIPELINE_ERRORS))
def test_pipeline_argument_errors(case, a_map, tmp_path, run_cli):
    args = [a_map(), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), *PIPELINE_ERRORS[case]]
    if "-a" not in args:
        args += ["-a", global_args["SUBSET_FASTA"]]
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)


PIPELINE_BAD_CHOICES = {
    "unknown algorithm": ["-A", "spinglass"],
    "unknown normalization": ["-m", "cosine"],
}


@pytest.mark.parametrize("case", list(PIPELINE_BAD_CHOICES))
def test_pipeline_bad_choices(case, tmp_path, run_cli):
    args = ["map.metacool", "-a", global_args["FASTA"], "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"),
            *PIPELINE_BAD_CHOICES[case]]  # fmt: skip
    result = run_cli("pipeline", *args, exit_code=2)
    assert f"'{PIPELINE_BAD_CHOICES[case][-1]}' is not one of" in result.output
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("existing", ["final_clusters"])
def test_pipeline_refuses_existing_outputs(existing, a_map, tmp_path, run_cli, logged):
    """A directory of clusters of an earlier run is refused without -F, named in
    the log."""
    (tmp_path / "out" / existing).mkdir(parents=True)
    args = [a_map(), "-a", global_args["SUBSET_FASTA"], "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp")]
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)
    message = f"{tmp_path / 'out' / existing} already existed. Remove directory or use -F argument to overwrite it."
    assert ("ERROR", message) in logged()


def test_pipeline_force_replaces_outputs(a_map, network_from_a_map, tmp_path, run_cli, logged):
    """With --force the output directories are emptied once the arguments and the
    assembly are checked, then the run goes on (here up to the partition, which
    network_from_a_map stops): final_clusters is replaced, and the recursive
    clusters that validation wrote there before they went to its temporary
    directory are removed (see metator.cli.EARLIER_OUTPUTS)."""
    for existing in ["final_clusters", "recursive_clusters"]:
        (tmp_path / "out" / existing).mkdir(parents=True)
        (tmp_path / "out" / existing / "stale.fa").write_text(">old\nA\n")
    others = ["overlapping_bin", "final_bin", "plot", "overlapping_clusters"]
    for existing in others:
        (tmp_path / "out" / existing).mkdir(parents=True)
        (tmp_path / "out" / existing / "earlier.fa").write_text(">old\nA\n")
    (tmp_path / "out" / "bin_summary.txt").write_text("earlier\n")
    args = [a_map(), "-a", global_args["SUBSET_FASTA"], "-F", "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp")]
    with pytest.raises(RuntimeError, match=r"^network built$"):
        run_cli("pipeline", *args)
    assert not list((tmp_path / "out").glob("*/stale.fa")) and not (tmp_path / "out" / "recursive_clusters").exists()
    assert all((tmp_path / "out" / name / "earlier.fa").exists() for name in others)
    assert (tmp_path / "out" / "bin_summary.txt").exists()
    assert not [message for _, message in logged() if message.startswith("Removed ")]


@pytest.mark.parametrize("existing", ["final_clusters"])
def test_validation_refuses_existing_outputs(existing, a_network_map, tmp_path, run_cli, logged):
    (tmp_path / "out" / existing).mkdir(parents=True)
    args = ("{MAP} -a {SUBSET_FASTA} -o {OUT} -T {TMP}").format(
        MAP=_partitioned(a_network_map(), tmp_path, run_cli), OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError):
        run_cli("validation", *args.split(" "))
    message = f"{tmp_path / 'out' / existing} already existed. Remove directory or use -F argument to overwrite it."
    assert ("ERROR", message) in logged()


def test_partition_prefix_names_the_clusters(metagenome_network, tmp_path, run_cli):
    """The --prefix names the overlapping clusters that partition MAP stores in
    the map; an unknown algorithm is refused before anything is written."""
    import metacooler

    metacool = shutil.copyfile(metagenome_network["map"], tmp_path / "network.metacool")
    before = Path(metacool).read_bytes()
    args = [str(metacool), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), "-i", "5", "-p", "sample-1"]
    run_cli("partition", *args, "-A", "spinglass", exit_code=2)
    assert Path(metacool).read_bytes() == before
    run_cli("partition", *args)
    clusters = metacooler.MetaCooler(str(metacool)).clusters
    overlapping = list(clusters.name[clusters["round"] == "overlapping"])
    assert overlapping and all(name.startswith("sample-1_metator_") and name.endswith("_00000") for name in overlapping)
    assert all(name.startswith("core_") for name in clusters.name[clusters["round"] == "core"])
