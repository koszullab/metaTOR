# Test for CLI toos of metator
# Commands are simply run to test crashes.


from pathlib import Path
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
    "NETWORK": "tests_data/network.txt",
    "CONTIGS": "tests_data/contig_data_network.txt",
    "CONTIGS_PARTITION": "tests_data/contig_data_partition.txt",
    "OVERLAPPING_BINS": "tests_data/overlapping_clusters/",
    "SUBSET_FASTA": "tests_data/subset.assembly.fa",
    "SUBSET_PAIRS": "tests_data/subset.pairs.gz",
    "SUBSET_MGES": "tests_data/subset.mges.txt",
}


@pytest.fixture(autouse=True)
def private_assemblies(private_link, monkeypatch):
    """Each test reads the assemblies of tests_data through links of its
    own (see private_link): two test sessions sharing tests_data removed
    each other's pyfastx index (FileNotFoundError)."""
    for key in ["FASTA", "SUBSET_FASTA"]:
        monkeypatch.setitem(global_args, key, private_link(global_args[key]))


@pytest.fixture(scope="session")
def tests_data_maps(tmp_path_factory):
    """The maps of the contigs of tests_data, in no cluster yet, as hicstuff
    or metacooler make them from the pairs (see tests.maps.contig_map):
    "full", of alignment.pairs on assembly.fa, and "subset", of
    subset.pairs.gz on subset.assembly.fa. Shared by the session: a test
    takes a copy of its own (see a_map)."""
    directory = tmp_path_factory.mktemp("tests_data_maps")
    maps = {}
    for name, assembly, pairs in [("full", "FASTA", "PAIRS"), ("subset", "SUBSET_FASTA", "SUBSET_PAIRS")]:
        contig_data = mtn.create_contig_data(global_args[assembly])
        maps[name] = contig_map(global_args[pairs], contig_data, directory / f"{name}.metacool", columns=False)
    return maps


@pytest.fixture
def a_map(tests_data_maps, tmp_path):
    """A copy of a map of tests_data (see tests_data_maps), its own: runs
    update their map in place."""

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
def test_partition(alg, tmp_path, run_cli):
    args = (
        "-a {SUBSET_FASTA} -c {CONTIGS} -i 5 -n {NETWORK} -o {OUT_TEST} -s 30000 -t 8 -T {TMP} -A {0} -F --no-clean-up"
    ).format(
        alg,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator partition", *args.split(" "))
    run_cli("partition", *args.split(" "))


def test_validation(tmp_path, run_cli):
    args = (
        "-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -i 5 -n {NETWORK} -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up"
    ).format(
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator validation", *args.split(" "))
    run_cli("validation", *args.split(" "))


def test_network_partition_validation(tmp_dir, a_map, run_cli):
    tmp = Path(tmp_dir, "tmp")
    out = Path(tmp_dir, "out_test")
    args = "{MAP} -a {SUBSET_FASTA} -o {OUT_TEST} -T {TMP} -e HindIII,DpnII --no-clean-up".format(
        MAP=a_map(),
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator network", *args.split(" "))
    run_cli("network", *args.split(" "))

    args = (
        "-a {SUBSET_FASTA} -c {OUT_TEST}/contig_data_network.txt -i 5 -n {OUT_TEST}/network.txt -o {OUT_TEST} -s 30000 -t 8 -T {TMP} -F --no-clean-up"
    ).format(
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator partition", *args.split(" "))
    run_cli("partition", *args.split(" "))

    args = (
        "-a {SUBSET_FASTA} -c {OUT_TEST}/contig_data_partition.txt -f {OUT_TEST}/overlapping_clusters/ -i 5 -n {OUT_TEST}/network.txt -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up"
    ).format(
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator validation", *args.split(" "))
    run_cli("validation", *args.split(" "))


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
    summary = Path(tmp_path, "out_test", "cluster_summary.txt").read_text().splitlines()
    assert set(genomes.name) - {"__unbinned__"} == {line.split("\t")[0] for line in summary[1:]}
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
    outputs = {}
    for binned, given in [("contig", a_map()), ("1kb", metacool)]:
        out = tmp_path / binned
        run_cli("network", given, "-a", assembly, "-e", "HindIII,DpnII", "-o", out, "-T", tmp_path / "tmp")
        outputs[binned] = out
    for name in ["network.txt", "contig_data_network.txt"]:
        assert (outputs["1kb"] / name).read_text() == (outputs["contig"] / name).read_text()
    assert sorted(path.name for path in outputs["1kb"].glob("network*.txt")) == ["network.txt"]
    (log,) = outputs["1kb"].glob("metator_network_*.log")
    assert f"INFO :: Information of {metacool}:" in log.read_text()


def test_network_from_the_columns_of_a_metacool(tmp_path, run_cli):
    """metator network MAP, from a map of the contigs with the columns that
    metaTOR measures with a depth file and an enzyme: without --depth and
    --enzyme, its Restriction_site and Shotgun_coverage columns give the
    network and the contig data that they give (theoritical_hit needs
    both), byte for byte."""
    contig_data = mtn.create_contig_data(global_args["FASTA"], global_args["DEPTH"], "HindIII,DpnII")
    metacool = contig_map(global_args["PAIRS"], contig_data, tmp_path / "contigs.metacool")
    common = [metacool, "-a", global_args["FASTA"], "-n", "theoritical_hit", "-T", tmp_path / "tmp"]
    run_cli("network", *common, "-d", global_args["DEPTH"], "-e", "HindIII,DpnII", "-o", tmp_path / "options")
    run_cli("network", *common, "-o", tmp_path / "map")
    for name in ["network.txt", "contig_data_network.txt"]:
        assert (tmp_path / "map" / name).read_bytes() == (tmp_path / "options" / name).read_bytes(), name
    (log,) = (tmp_path / "map").glob("metator_network_*.log")
    logged = log.read_text()
    for what, column in [("restriction sites", "Restriction_site"), ("shotgun coverage", "Shotgun_coverage")]:
        assert f"INFO :: The {what} of the contigs: the {column} column of {metacool}." in logged


@pytest.mark.parametrize("data", ["full", "subset"])
def test_pipeline_network_is_that_of_network(data, a_map, network_from_a_map, tmp_path, run_cli):
    """A pipeline run reads its network from its map, as metator network
    does: network.txt and contig_data_network.txt are those of network MAP,
    byte for byte, on the map of the pairs of tests_data (alignment.pairs
    holds reads at position 0) and on the map of the subset. The run stops
    at the partition (see network_from_a_map)."""
    assembly = global_args["FASTA" if data == "full" else "SUBSET_FASTA"]
    depth = ["-d", global_args["DEPTH"]] if data == "full" else []
    common = ["-a", assembly, *depth, "-e", "HindIII,DpnII", "-T", tmp_path / "tmp"]
    metacool = a_map(data)
    run_cli("network", metacool, *common, "-o", tmp_path / "network")
    with pytest.raises(RuntimeError, match=r"^network built$"):
        run_cli("pipeline", metacool, *common, "-o", tmp_path / "pipeline")
    # Read by metator network, then by pipeline.
    assert network_from_a_map == [str(metacool), str(metacool)]
    for name in ["network.txt", "contig_data_network.txt"]:
        assert (tmp_path / "pipeline" / name).read_bytes() == (tmp_path / "network" / name).read_bytes(), name
    assert sorted(path.name for path in (tmp_path / "pipeline").glob("network*.txt")) == ["network.txt"]


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
    """Without -o, the network is written to the current directory, with the
    log of the run."""
    args = [
        str(a_map()),
        "-a", str(Path(global_args["SUBSET_FASTA"]).resolve()),
    ]  # fmt: skip
    monkeypatch.chdir(tmp_path)
    run_cli("network", *args)
    assert (tmp_path / "network.txt").stat().st_size > 0
    assert (tmp_path / "contig_data_network.txt").stat().st_size > 0
    (log,) = tmp_path.glob("metator_network_*.log")
    assert f"## Command used: metator network {args[0]} -a {args[2]}\n" in log.read_text()


def test_invalid_prefix_is_rejected(tmp_path, run_cli):
    """miComplete drops dots from bin names, so a dotted --prefix is refused."""
    args = ("-a {SUBSET_FASTA} -c {CONTIGS} -n {NETWORK} -o {OUT} -T {TMP} -p my.sample").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError, match="Invalid --prefix"):
        run_cli("partition", *args.split(" "))


def test_spinglass_is_rejected_before_any_work(tmp_path, monkeypatch, run_cli):
    """spinglass is gone: refuse it up front, not after miComplete has run."""

    def not_reached(*args, **kwargs):
        raise AssertionError("recursive decontamination should not start")

    monkeypatch.setattr(mtc.mtv, "recursive_decontamination", not_reached)
    args = ("-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -n {NETWORK} -o {OUT} -T {TMP} -A spinglass").format(
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


@pytest.mark.parametrize("existing", ["overlapping_clusters", "recursive_clusters", "final_clusters"])
def test_pipeline_refuses_existing_outputs(existing, a_map, tmp_path, run_cli, logged):
    """A directory of clusters of an earlier run is refused without -F,
    named in the log."""
    (tmp_path / "out" / existing).mkdir(parents=True)
    args = [a_map(), "-a", global_args["SUBSET_FASTA"], "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp")]
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)
    message = f"{tmp_path / 'out' / existing} already existed. Remove directory or use -F argument to overwrite it."
    assert ("ERROR", message) in logged()


def test_pipeline_force_replaces_outputs(a_map, network_from_a_map, tmp_path, run_cli, logged):
    """With --force the output directories are emptied once the arguments
    and the assembly are checked, then the run goes on (here up to the
    partition, which network_from_a_map stops). The outputs that metaTOR
    1.x named after bins, and the plots of an earlier run, are no output
    of 2.0: they are neither refused nor removed."""
    for existing in ["overlapping_clusters", "recursive_clusters", "final_clusters"]:
        (tmp_path / "out" / existing).mkdir(parents=True)
        (tmp_path / "out" / existing / "stale.fa").write_text(">old\nA\n")
    others = ["overlapping_bin", "final_bin", "plot"]
    for existing in others:
        (tmp_path / "out" / existing).mkdir(parents=True)
        (tmp_path / "out" / existing / "earlier.fa").write_text(">old\nA\n")
    (tmp_path / "out" / "bin_summary.txt").write_text("earlier\n")
    args = [a_map(), "-a", global_args["SUBSET_FASTA"], "-F", "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp")]
    with pytest.raises(RuntimeError, match=r"^network built$"):
        run_cli("pipeline", *args)
    assert not list((tmp_path / "out").glob("*/stale.fa"))
    assert all((tmp_path / "out" / name / "earlier.fa").exists() for name in others)
    assert (tmp_path / "out" / "bin_summary.txt").exists()
    assert not [message for _, message in logged() if message.startswith("Removed ")]


@pytest.mark.parametrize("existing", ["recursive_clusters", "final_clusters"])
def test_validation_refuses_existing_outputs(existing, tmp_path, run_cli):
    (tmp_path / "out" / existing).mkdir(parents=True)
    args = ("-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -n {NETWORK} -o {OUT} -T {TMP}").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError):
        run_cli("validation", *args.split(" "))


def test_partition_existing_outputs(metagenome_inputs, metagenome_network, tmp_path, run_cli):
    """An existing overlapping_clusters is refused, or replaced with --force; the
    --prefix names the bins."""
    fasta_dir = tmp_path / "out" / "overlapping_clusters"
    fasta_dir.mkdir(parents=True)
    (fasta_dir / "stale.fa").write_text(">old\nA\n")
    args = ["-a", str(metagenome_inputs["assembly"]), "-c", str(metagenome_network["contig_data"]),
            "-n", str(metagenome_network["network"]), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"),
            "-i", "5", "-s", "50000", "-p", "sample-1"]  # fmt: skip
    with pytest.raises(ValueError):
        run_cli("partition", *args)
    # An unknown algorithm is refused before anything is replaced.
    run_cli("partition", *args, "-F", "-A", "spinglass", exit_code=2)
    assert (fasta_dir / "stale.fa").exists()
    run_cli("partition", *args, "-F")
    bins = sorted(p.name for p in fasta_dir.iterdir())
    assert "stale.fa" not in bins
    assert bins and all(name.startswith("sample-1_metator_") for name in bins)
