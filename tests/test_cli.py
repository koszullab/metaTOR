# Tests of the metator command line (metator.cli): parity with the docopt
# command line it replaced, help, version and usage errors, and the
# bookkeeping of the commands (log files, temporary directories, default
# outputs), with the heavy steps replaced by stubs.

import bz2
import gzip
import json
import os
import re
import shlex
import shutil
import sys
import zipfile
from pathlib import Path

import click
import pytest

import metator.cli as mtcli
import metator.io as mio
import metator.metacool as mtmc
from metator.version import __version__

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"

# The writers themselves, which the stub pipeline replaces.
mtmc_write_contig_map = mtmc.write_contig_map
mtmc_write_library_maps = mtmc.write_library_maps
mtmc_reassign_metacool = mtmc.reassign_metacool

# The docopt command line that rich-click replaced (that of 1.3.10, with the
# mge options renamed and host removed), as docopt read it: see the commit
# that added the file. Commands and options added since are added to it:
# metacool, the --no-metacool and --output-map options of pipeline, and the
# --output-map option of partition and validation. The 2.0 command line
# changed it on purpose (see docs/migration.md): the map that network and
# pipeline start from is their first argument (MAP), so their --forward is
# not required, and their --edge is gone; the map that partition and
# validation update is their first argument too.
INVENTORY = json.loads((TESTS_DATA / "cli_inventory.json").read_text())["metator"]

# Options whose value is one of a list, wherever they appear. --object takes
# any column name of the contig data file. The --start of network and
# pipeline took metacool before 2.0, which the map, their first argument,
# replaced.
CHOICES = {
    "--aligner": ["bowtie2", "bwa"],
    "--aligner-mode": ["normal", "iterative", "cutsite"],
    "--algorithm": ["louvain", "leiden"],
    "--normalization": ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"],
    "--start": ["fastq", "bam", "pair"],
    "--mat-fmt": ["graal", "bg2", "cool"],
}


def _info():
    ctx = mtcli.cli.make_context("metator", ["--help"], resilient_parsing=True)
    return mtcli.cli.to_info_dict(ctx)


def _options(params):
    """Options in the inventory's format, by long name (--help left out)."""
    options = {}
    for param in params:
        if param["param_type_name"] != "option" or param["name"] == "help":
            continue
        longs = [opt for opt in param["opts"] if opt.startswith("--")]
        shorts = [opt for opt in param["opts"] if not opt.startswith("--")]
        assert len(longs) == 1 and len(shorts) <= 1 and not param["secondary_opts"], param["opts"]
        default = param["default"]
        options[longs[0]] = {
            "short": shorts[0] if shorts else None,
            "long": longs[0],
            "takes_value": not param["is_flag"],
            "default": default if param["is_flag"] or default is None else str(default),
            "required": param["required"],
        }
    return options


def _positional(params):
    return [
        {"name": param["name"], "required": param["required"], "repeatable": param["nargs"] == -1}
        for param in params
        if param["param_type_name"] == "argument"
    ]


def _required_args(command):
    """A command line with every required option and argument of a command,
    with made-up values."""
    args = []
    for option in INVENTORY["commands"][command]["options"]:
        if option["required"]:
            args += [option["short"] or option["long"], "x"]
    args += ["p.pairs" for _ in INVENTORY["commands"][command]["positional"]]
    return args


def test_commands_are_the_docopt_commands():
    assert sorted(_info()["commands"]) == sorted(INVENTORY["commands"])


@pytest.mark.parametrize("command", sorted(INVENTORY["commands"]))
def test_options_are_the_docopt_options(command):
    """Same options, short aliases, defaults and required options as the
    docopt command line, and the same positional arguments."""
    params = _info()["commands"][command]["params"]
    frozen = {option["long"]: option for option in INVENTORY["commands"][command]["options"]}
    found = _options(params)
    assert sorted(found) == sorted(frozen)
    for long, option in frozen.items():
        expected = {key: option[key] for key in ["short", "long", "takes_value", "default", "required"]}
        assert found[long] == expected, long
    assert _positional(params) == INVENTORY["commands"][command]["positional"]


@pytest.mark.parametrize("command", sorted(INVENTORY["commands"]))
def test_value_names_and_choices(command):
    """Choice options list their values; the others show docopt's value
    name (FILE, INT...)."""
    frozen = {option["long"]: option for option in INVENTORY["commands"][command]["options"]}
    for param in mtcli.cli.commands[command].params:
        if not isinstance(param, click.Option):
            continue
        long = next(opt for opt in param.opts if opt.startswith("--"))
        if long in CHOICES:
            assert isinstance(param.type, click.Choice), long
            assert list(param.type.choices) == CHOICES[long], long
            assert frozen[long]["default"] in CHOICES[long]
        else:
            assert not isinstance(param.type, click.Choice), long
        if frozen[long]["takes_value"] and not isinstance(param.type, click.Choice):
            assert param.metavar == frozen[long]["metavar"], long


def test_global_options():
    found = {param["opts"][-1]: param["opts"] for param in _info()["params"]}
    frozen = {option["long"]: [option["short"], option["long"]] for option in INVENTORY["options"]}
    assert found == frozen


@pytest.mark.parametrize("command", sorted(INVENTORY["commands"]))
@pytest.mark.parametrize("help_option", ["-h", "--help"])
def test_command_help(command, help_option, run_cli):
    output = run_cli(command, help_option).output
    assert f"Usage: metator {command}" in output
    for option in INVENTORY["commands"][command]["options"]:
        assert option["long"] in output


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_help_of_a_run_from_a_map(command, run_cli):
    """The help shows both inputs: the map, the first argument, and the
    reads, alignments or pairs of -1; and that every pair counts."""
    output = " ".join(run_cli(command, "--help").output.split())
    assert f"Usage: metator {command} [OPTIONS] [MAP]" in output
    assert f"metator {command} MAP -a ASSEMBLY" in output and f"metator {command} -1 FORWARD -2 REVERSE" in output
    assert "every pair" in output


def test_group_help(run_cli):
    output = run_cli("-h").output
    for command in INVENTORY["commands"]:
        assert command in output
    assert "Clustering" in output


def test_version(run_cli):
    assert run_cli("-v").output == f"{__version__}\n"


@pytest.mark.parametrize(
    "command, missing",
    [
        (command, option["long"])
        for command, frozen in sorted(INVENTORY["commands"].items())
        for option in frozen["options"]
        if option["required"]
    ],
)
def test_missing_required_option(command, missing, run_cli):
    args = _required_args(command)
    short = next(option["short"] for option in INVENTORY["commands"][command]["options"] if option["long"] == missing)
    index = args.index(short)
    del args[index : index + 2]
    result = run_cli(command, *args, exit_code=2)
    assert f"Missing option '{short}' / '{missing}'" in result.output


@pytest.mark.parametrize(
    "command",
    [command for command, frozen in INVENTORY["commands"].items() if any(arg["required"] for arg in frozen["positional"])],
)
def test_missing_pairs_files(command, run_cli):
    args = _required_args(command)[:-1]
    result = run_cli(command, *args, exit_code=2)
    assert "Missing argument 'PAIRSFILE...'" in result.output


@pytest.mark.parametrize(
    "command, option",
    [
        (command, option["long"])
        for command, frozen in sorted(INVENTORY["commands"].items())
        for option in frozen["options"]
        if option["long"] in CHOICES
    ],
)
def test_bad_choice(command, option, run_cli, tmp_path):
    result = run_cli(command, *_required_args(command), option, "bogus", exit_code=2)
    assert f"'{option}': 'bogus' is not one of" in result.output
    # Refused before the command starts: no log, no temporary directory.
    assert sorted(os.listdir(tmp_path)) == ["tests_data"]


@pytest.mark.parametrize("value", ["8.5", "eight"])
def test_bad_number(value, run_cli):
    result = run_cli("pairs", "-t", value, "p.pairs", exit_code=2)
    assert f"'{value}' is not a valid integer" in result.output


def test_unknown_option(run_cli):
    """Long options cannot be abbreviated (docopt accepted --thread)."""
    result = run_cli("pairs", "--thread", "2", "p.pairs", exit_code=2)
    assert "No such option" in result.output and "--thread" in result.output


# Helpers


def test_private_tmp_dir(tmp_path):
    base = tmp_path / "tmp"
    with mtcli.private_tmp_dir(str(base)) as tmp_dir:
        assert Path(tmp_dir).parent == base
        (Path(tmp_dir) / "file").write_text("x")
    # The run created --tmpdir: both go.
    assert not base.exists()

    # --tmpdir already existed: only the run's directory goes.
    base.mkdir()
    with mtcli.private_tmp_dir(str(base)) as tmp_dir:
        pass
    assert base.exists() and not Path(tmp_dir).exists()

    # Something else was written to --tmpdir meanwhile: it stays.
    shutil.rmtree(base)
    with mtcli.private_tmp_dir(str(base)) as tmp_dir:
        (base / "other").write_text("x")
    assert sorted(os.listdir(base)) == ["other"]


def test_private_tmp_dir_kept(tmp_path, monkeypatch):
    """Kept with --no-clean-up, and when the run fails."""
    monkeypatch.chdir(tmp_path)
    with mtcli.private_tmp_dir(None, clean_up=False) as tmp_dir:
        assert Path(tmp_dir).parent == Path("tmp")
    assert Path(tmp_dir).is_dir()

    with pytest.raises(ValueError), mtcli.private_tmp_dir("failed") as tmp_dir:
        raise ValueError
    assert Path(tmp_dir).is_dir()


@pytest.mark.parametrize(
    "prefix, expected",
    [(None, "metator"), ("", "metator"), ("sample_1-a", "sample_1-a_metator")],
)
def test_bin_prefix(prefix, expected):
    assert mtcli.bin_prefix(prefix) == expected


@pytest.mark.parametrize("prefix", ["my.sample", "a/b", "é"])
def test_bin_prefix_refused(prefix):
    with pytest.raises(ValueError, match="Invalid --prefix"):
        mtcli.bin_prefix(prefix)


def test_log_header_records_the_command_line(tmp_path, monkeypatch):
    log = tmp_path / "run.log"
    monkeypatch.setattr(sys, "argv", ["/usr/bin/metator", "qc", "-p", "my sample", "a.pairs"])
    mtcli.generate_log_header(str(log))
    mtcli.generate_log_footer(str(log))
    lines = log.read_text().splitlines()
    assert lines[0] == f"## MetaTOR: v{__version__} log file"
    assert lines[1].startswith("## date: ")
    assert lines[2] == "## Command used: metator qc -p 'my sample' a.pairs"
    assert lines[-1].startswith("## date: ")


# The bookkeeping of the commands, with the heavy steps replaced by stubs.


@pytest.fixture
def stub_pipeline(monkeypatch, tmp_path):
    """Replace the steps of the pipeline by stubs that write their outputs
    and record their arguments."""
    calls = {}

    def get_contact_pairs(*args, sort_pairs=True, keep_temp=False):
        calls["get_contact_pairs"] = args
        calls["sort_pairs"] = sort_pairs
        calls["keep_temp"] = keep_temp
        libraries = range(len(args[0].split(",")))
        if sort_pairs:
            # Sorted, compressed and indexed by pypairix, in the output
            # directory, as the alignment step leaves them.
            pairs = [os.path.join(args[10], f"alignment_{i}_sorted.pairs.gz") for i in libraries]
        else:
            # As they are made, in the temporary directory.
            pairs = [os.path.join(args[11], f"alignment_{i}.pairs") for i in libraries]
        for path in pairs:
            Path(path).write_text("## pairs format v1.0\n")
            if sort_pairs:
                Path(f"{path}.px2").write_text("x")
        calls["pairs"] = pairs
        return pairs, {}, {}

    def create_contig_data(*args):
        calls["create_contig_data"] = args
        return {}, {}

    def contig_data_from_map(*args):
        # A run from a map reads no column of it (see
        # test_pipeline_reads_the_columns_of_its_map).
        calls["contig_data_from_map"] = args
        return []

    def alignment_to_contacts(alignment_files, contig_data, edge, hit_data, out_dir, network, contig_data_file, *args):
        calls["alignment_to_contacts"] = (alignment_files, edge, out_dir, *args)
        for name in [network, contig_data_file]:
            Path(out_dir, name).write_text("x")
        return os.path.join(out_dir, network), os.path.join(out_dir, contig_data_file)

    def metacool_to_contacts(metacool, contig_data, out_dir, network, contig_data_file, *args, libraries=None):
        calls["metacool_to_contacts"] = (metacool, out_dir, network, contig_data_file, *args)
        calls["libraries"] = libraries
        for name in [network, contig_data_file]:
            Path(out_dir, name).write_text("x")
        for i, _ in enumerate(libraries or []):
            Path(out_dir, f"network_{i}.txt").write_text("x")
        return os.path.join(out_dir, network), os.path.join(out_dir, contig_data_file)

    def partition(*args):
        calls["partition"] = args
        out_dir = args[6]
        # What the network step left in the output directory, and of the
        # pairs files that the run made.
        calls["partition_outdir"] = sorted(os.listdir(out_dir))
        if "pairs" in calls:
            calls["pairs_at_partition"] = [os.path.exists(path) for path in calls["pairs"]]
        Path(out_dir, "contig_data_partition.txt").write_text("x")
        return os.path.join(out_dir, "clustering_matrix_partition"), os.path.join(out_dir, "contig_data_partition.txt")

    def recursive_decontamination(*args, **kwargs):
        calls["recursive_decontamination"] = args
        # The map that the heatmap of the clusters is drawn from.
        calls["heatmap_map"] = kwargs.get("metacool")
        # pyfastx indexes the assembly while extracting the bins.
        Path(f"{args[1]}.fxi").write_text("x")
        return os.path.join(args[8], "clustering_matrix_validation")

    def generate_fasta_index(fasta, aligner, tmp_dir):
        calls["generate_fasta_index"] = (fasta, aligner, tmp_dir)
        return os.path.join(tmp_dir, "index")

    def write_metacool(pairs_files, contig_data, outfile, **kwargs):
        # metator metacool only: the pipeline writes its map with
        # write_contig_map and reassign_metacool.
        calls["write_metacool"] = (pairs_files, contig_data, outfile, kwargs)

    def write_contig_map(pairs_files, contig_data, outfile, **kwargs):
        # The pairs files must be there.
        calls["write_contig_map"] = (pairs_files, [os.path.exists(path) for path in pairs_files], contig_data, outfile, kwargs)
        Path(outfile).write_text("map of the contigs")

    def write_library_maps(pairs_files, contig_data, outfile, collection, **kwargs):
        # The pairs files must be there.
        existed = [os.path.exists(path) for path in pairs_files]
        calls["write_library_maps"] = (pairs_files, existed, contig_data, outfile, collection, kwargs)
        Path(outfile).write_text("map of the contigs")
        Path(collection).write_text("maps of the libraries")

    def reassign_metacool(metacool, contig_data, outfile, **kwargs):
        # The input map must still be there.
        calls["reassign_metacool"] = (metacool, os.path.exists(metacool), contig_data, outfile, kwargs)

    monkeypatch.setattr("metator.align.get_contact_pairs", get_contact_pairs)
    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    monkeypatch.setattr("metator.network.contig_data_from_map", contig_data_from_map)
    monkeypatch.setattr("metator.network.alignment_to_contacts", alignment_to_contacts)
    monkeypatch.setattr("metator.network.metacool_to_contacts", metacool_to_contacts)
    monkeypatch.setattr("metator.partition.partition", partition)
    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)
    monkeypatch.setattr("metator.io.generate_fasta_index", generate_fasta_index)
    monkeypatch.setattr("metator.metacool.write_metacool", write_metacool)
    monkeypatch.setattr("metator.metacool.write_contig_map", write_contig_map)
    monkeypatch.setattr("metator.metacool.write_library_maps", write_library_maps)
    monkeypatch.setattr("metator.metacool.reassign_metacool", reassign_metacool)
    return calls


@pytest.mark.parametrize("start", ["fastq", "bam"])
@pytest.mark.parametrize("clean_up", [True, False])
@pytest.mark.parametrize("libraries", [1, 2])
def test_pipeline_from_reads(start, clean_up, libraries, stub_pipeline, tmp_path, run_cli):
    """From reads or alignments: the pairs files and the contig data of the
    network step are intermediate files. The pairs files are written to the
    temporary directory as they are made, unsorted, and removed once the
    network is built, before the partition; with -N, they are sorted,
    compressed, indexed and kept in the output directory, as before (where
    the run sorted them always, then removed them at its end). The network
    is read from the map of the contigs. One library: copied to
    network_0.txt, the network of that library. Two: their sum, with a map
    of each library in metator.mmetacool, from which the network of each
    is read."""
    if start == "fastq":
        # A fasta: its index is built.
        assembly = tmp_path / "assembly.fa"
        shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    else:
        # A bowtie2 index: its fasta is retrieved.
        assembly = TESTS_DATA / "assembly"
    out = tmp_path / "out"
    forward, reverse = ",".join(["a1.fq", "b1.fq"][:libraries]), ",".join(["a2.fq", "b2.fq"][:libraries])
    args = ["-1", forward, "-2", reverse, "-a", assembly, "-S", start, "-B", "cutsite", "-q", "20", "-t", "2", "-o", out,
            "-T", tmp_path / "tmp"]  # fmt: skip
    if not clean_up:
        args.append("-N")
    run_cli("pipeline", *args)

    calls = stub_pipeline
    forward_arg, reverse_arg, index, fasta, aligner, mode, min_quality, start_arg, *_, threads = calls["get_contact_pairs"]
    assert (forward_arg, reverse_arg, aligner, start_arg) == (forward, reverse, "bowtie2", start)
    # cutsite needs an enzyme.
    assert mode == "iterative"
    assert min_quality == 20 and threads == 2
    if start == "fastq":
        assert calls["generate_fasta_index"][0] == str(assembly)
        assert (fasta, index) == (str(assembly), calls["generate_fasta_index"][2] + "/index")
    else:
        assert "generate_fasta_index" not in calls
        assert index == str(assembly) and fasta.endswith(".fa")
    assert calls["partition"][1] == fasta
    assert calls["recursive_decontamination"][6] == 10
    # The heatmap of the clusters is drawn from the map of the contigs.
    assert calls["heatmap_map"] == str(out / "metator.metacool")
    assert not Path(f"{fasta}.fxi").exists()
    assert not (out / "contig_data_network.txt").exists()
    assert not (out / "contig_data_partition.txt").exists()
    assert (out / "network.txt").exists()
    # The pairs files and their pypairix index go together.
    pairs = [Path(path) for path in calls["pairs"]]
    assert calls["sort_pairs"] == (not clean_up) and calls["keep_temp"] == (not clean_up)
    if clean_up:
        assert [path.name for path in pairs] == [f"alignment_{i}.pairs" for i in range(libraries)]
        assert {path.parent.parent for path in pairs} == {tmp_path / "tmp"}
    else:
        assert pairs == [out / f"alignment_{i}_sorted.pairs.gz" for i in range(libraries)]
    assert calls["pairs_at_partition"] == [not clean_up] * libraries
    assert [path.exists() for path in pairs] == [not clean_up] * libraries
    assert [Path(f"{path}.px2").exists() for path in pairs] == [not clean_up] * libraries
    assert (tmp_path / "tmp").exists() == (not clean_up)
    assert bool(list(out.glob("alignment_*"))) == (not clean_up)
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"## Command used: metator pipeline -1 {forward} -2 {reverse} -a {assembly}" in log.read_text()

    # The contact map of the contigs, built from the pairs files of the run
    # as soon as they are there, before the network, which is read from it.
    order = list(calls)
    builder = "write_contig_map" if libraries == 1 else "write_library_maps"
    assert order.index("get_contact_pairs") < order.index(builder) < order.index("metacool_to_contacts")
    assert {"write_contig_map", "write_library_maps"} & set(calls) == {builder}
    assert "alignment_to_contacts" not in calls
    source, out_dir, *files, tmp_dir, threads, normalization, self_contacts = calls["metacool_to_contacts"]
    assert (source, out_dir, files) == (str(out / "metator.metacool"), str(out), ["network.txt", "contig_data_network.txt"])
    assert (Path(tmp_dir).parent, threads, normalization, self_contacts) == (tmp_path / "tmp", 2, "empirical_hit", False)
    collection = str(out / "metator.mmetacool")
    if libraries == 1:
        assert calls["libraries"] is None
        assert (out / "network_0.txt").read_text() == (out / "network.txt").read_text()
        pairs_files, existed, contig_data, outfile, kwargs = calls["write_contig_map"]
    else:
        assert calls["libraries"] == [f"{collection}::entries/library_{i}" for i in range(libraries)]
        pairs_files, existed, contig_data, outfile, built, kwargs = calls["write_library_maps"]
        assert built == collection
    assert [(out / f"network_{i}.txt").exists() for i in range(libraries)] == [True] * libraries
    assert pairs_files == [str(path) for path in pairs] and existed == [True] * libraries
    assert (contig_data, outfile) == ({}, str(out / "metator.metacool"))
    assert (kwargs["assembly_name"], kwargs["threads"], kwargs["force"]) == ("assembly", 2, True)
    assert Path(kwargs["tmp_dir"]).parent == tmp_path / "tmp"
    metadata = kwargs["metadata"]
    assert metadata["command"].startswith(f"metator pipeline -1 {forward} -2 {reverse} -a")
    assert metadata["parameters"]["start"] == start and metadata["parameters"]["no_metacool"] is False
    assert "metator.metacool" in calls["partition_outdir"]
    # The map of the bins, written last over it, from the final tables: the
    # pairs files are not read again.
    assert order[-1] == "reassign_metacool" and "write_metacool" not in calls
    source, existed, contig_data, outfile, kwargs = calls["reassign_metacool"]
    assert (source, existed, outfile) == (str(out / "metator.metacool"), True, str(out / "metator.metacool"))
    assert contig_data == str(out / "contig_data_final.txt")
    assert (kwargs["bin_summary"], kwargs["force"], kwargs["metadata"]) == (str(out / "cluster_summary.txt"), True, metadata)
    # The bin summary goes to metacooler as it is: no temporary copy.
    assert "tmp_dir" not in kwargs


@pytest.mark.parametrize(
    "options, written",
    [([], "metator.metacool"), (["-p", "sample_1"], "sample_1.metacool"), (["--no-metacool"], "metator.metacool")],
)
def test_pipeline_metacool_name(options, written, stub_pipeline, tmp_path, run_cli):
    """<prefix>.metacool, metator.metacool without --prefix: the map of the
    contigs is built there, the network read from it, and the map of the
    bins written over it at the end. With --no-metacool, the map of the
    contigs is built and read all the same, and removed at the end."""
    out = tmp_path / "out"
    run_cli("pipeline", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp", *options)
    outfile = str(out / written)
    pairs_files, _, _, built, _ = stub_pipeline["write_contig_map"]
    assert pairs_files == ["a.pairs"] and built == outfile
    assert stub_pipeline["metacool_to_contacts"][0] == outfile and "alignment_to_contacts" not in stub_pipeline
    assert written in stub_pipeline["partition_outdir"]
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    if "--no-metacool" in options:
        assert "reassign_metacool" not in stub_pipeline and not (out / written).exists()
        assert f"INFO :: Removing {outfile}, the contact map of the contigs (--no-metacool)." in logged
    else:
        source, _, _, replaced, _ = stub_pipeline["reassign_metacool"]
        assert source == replaced == outfile


@pytest.mark.parametrize("no_metacool", [False, True])
def test_pipeline_maps_of_several_libraries(no_metacool, stub_pipeline, tmp_path, run_cli):
    """Two libraries: the maps of each, metator.mmetacool, next to the map of
    their sum, which an earlier run's are removed first (both: a run with
    one library left a sum, a run with two a collection too). The sum gets
    the map of the clusters at the end; with --no-metacool, both are
    removed."""
    out = tmp_path / "out"
    out.mkdir()
    earlier = [out / "metator.metacool", out / "metator.mmetacool"]
    for path in earlier:
        path.write_text("an earlier run's")
    pairs = [tmp_path / "a.pairs", tmp_path / "b.pairs"]
    for path in pairs:
        path.write_text("## pairs format v1.0\n")
    run_cli("pipeline", "-1", ",".join(map(str, pairs)), "-a", _assembly(tmp_path), "-S", "pair", "-o", out,
            "-T", tmp_path / "tmp", *(["--no-metacool"] if no_metacool else []))  # fmt: skip
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    for path in earlier:
        assert f"INFO :: Removing {path}, the contact map of an earlier run." in logged
    assert stub_pipeline["write_library_maps"][4] == str(out / "metator.mmetacool")
    if no_metacool:
        assert not any(path.exists() for path in earlier)
        assert f"INFO :: Removing {earlier[1]}, the contact maps of the libraries (--no-metacool)." in logged
    else:
        assert all(path.exists() for path in earlier)
        assert stub_pipeline["reassign_metacool"][3] == str(earlier[0])


def _assembly(tmp_path):
    """A copy of an assembly: the stub validation writes its pyfastx index."""
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    return assembly


def _final_tables(monkeypatch, contig, final_bin="metator_00001_00000", summary=None):
    """Have the stub validation also write the final contig table, where the
    contig ``contig`` is in the bin ``final_bin`` (ND: in none) and c2 is in
    none, and the bin summary (``summary``, if given)."""
    stub = mtcli.mtv.recursive_decontamination

    def recursive_decontamination(*args, **kwargs):
        out_dir = args[8]
        Path(out_dir, "contig_data_final.txt").write_text(
            "\tID\tName\tSize\tGC_content\tHit\tFinal_cluster\n"
            f"0\t1\t{contig}\t100\t40.5\t-\t{final_bin}\n1\t2\tc2\t50\t35.0\t-\tND\n"
        )
        rows = "" if final_bin == "ND" else f"{final_bin}\t0.9\n"
        Path(out_dir, "cluster_summary.txt").write_text(summary or f"\tCompleteness\n{rows}")
        return stub(*args, **kwargs)

    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)


def _contig_data(monkeypatch, contigs):
    """Have the stub create_contig_data give the contig data of these
    contigs: their length and GC content, by name, and nothing else
    measured (no depth file, no enzyme), as create_contig_data does."""

    def create_contig_data(*args):
        return {
            name: {"id": i, "length": length, "GC": gc, "hit": 0, "coverage": "-", "RS": "-"}
            for i, (name, (length, gc)) in enumerate(contigs.items(), start=1)
        }, None

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)


# The contigs of the stub validation's tables (see _final_tables).
CONTIGS = {"c1": (100, 40.5), "c2": (50, 35.0)}


def _pairs(tmp_path, contig="c1"):
    """A pairs file of two pairs, one between contig and c2, one on contig."""
    pairs = tmp_path / "a.pairs"
    pairs.write_text(f"## pairs format v1.0\nr1\t{contig}\t10\tc2\t20\t+\t-\nr2\t{contig}\t10\t{contig}\t90\t+\t-\n")
    return pairs


def _pair_run(tmp_path, out):
    """The arguments of a pipeline run from the pairs of _pairs."""
    return ["-1", _pairs(tmp_path), "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp"]


def _map_at_partition(monkeypatch, outfile):
    """Have the stub partition record the genomes of the map at outfile, as
    the network step leaves it, or None when there is none. Returns the
    record."""
    import metacooler

    seen = []
    stub = mtcli.mtp.partition

    def partition(*args):
        seen.append(list(metacooler.MetaCooler(str(outfile)).genomes.name) if os.path.exists(outfile) else None)
        return stub(*args)

    monkeypatch.setattr("metator.partition.partition", partition)
    return seen


@pytest.mark.parametrize("no_metacool", [False, True])
def test_pipeline_non_ascii_contigs_stop_the_run(no_metacool, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """Contig names that are not ASCII in the contig data (those of the
    assembly are refused before any alignment: see
    test_pipeline_refuses_non_ascii_contigs_before_any_alignment): the map
    of the contigs cannot be built, which stops the run before its network,
    --no-metacool or not, with the MetacoolError of write_contig_map, its
    reason logged as an error in the log file of the run (without a
    traceback: the error is expected, nor what metacooler needs of pairs
    files: they are not at fault). No .metacool is left."""
    monkeypatch.setattr("metator.metacool.write_contig_map", mtmc_write_contig_map)
    _contig_data(monkeypatch, {"contig_é": (100, 40.5), "c2": (50, 35.0)})
    out = tmp_path / "out"
    args = ["-1", _pairs(tmp_path, "contig_é"), "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp"]
    reason = (
        "A .metacool file needs ASCII contig names (cooler reads them as ASCII text), and 1 contig name(s) are not: 'contig_é'."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(reason)}$"):
        run_cli("pipeline", *args, *(["--no-metacool"] if no_metacool else []))
    assert not list(out.glob("*.metacool*"))
    assert not {"alignment_to_contacts", "metacool_to_contacts", "partition"} & set(stub_pipeline)
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"ERROR :: The contact map of the contigs cannot be built at {out / 'metator.metacool'}: {reason}\n" in logged
    assert "Traceback" not in logged and "WARNING" not in logged


@pytest.mark.parametrize("no_clean_up", [False, True])
def test_pipeline_scaffolds_from_sorted_pairs(no_clean_up, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """--scaffold reads the pairs files of the run, through their pypairix
    index: from reads or alignments, they are sorted, compressed and indexed
    in the output directory, kept until the scaffolding, then removed at
    the end of the run (kept with -N), as before. A guard: the run sorted
    them always."""
    _final_tables(monkeypatch, "c1")
    scaffolded = []

    def parallel_scaffold(bin_name, final_fasta_dir, alignment_files, scaffold_fasta_dir, junctions):
        scaffolded.append((bin_name, [(path, os.path.exists(path), os.path.exists(f"{path}.px2")) for path in alignment_files]))

    class Pool:
        """A pool of workers that runs its tasks in this process."""

        def __init__(self, processes, initializer, initargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *exc_info):
            return False

        def map(self, task, items):
            return [task(item) for item in items]

    monkeypatch.setattr("metator.scaffold.parallel_scaffold", parallel_scaffold)
    monkeypatch.setattr("metator.cli.mp.Pool", Pool)
    out = tmp_path / "out"
    args = ["-1", "a1.bam", "-2", "a2.bam", "-a", TESTS_DATA / "assembly", "-S", "bam", "--scaffold", "-o", out,
            "-T", tmp_path / "tmp"]  # fmt: skip
    run_cli("pipeline", *args, *(["-N"] if no_clean_up else []))
    pairs = out / "alignment_0_sorted.pairs.gz"
    assert stub_pipeline["sort_pairs"] and stub_pipeline["pairs_at_partition"] == [True]
    assert scaffolded == [("metator_00001_00000", [(str(pairs), True, True)])]
    assert pairs.exists() == no_clean_up and Path(f"{pairs}.px2").exists() == no_clean_up


@pytest.mark.parametrize("no_metacool", [False, True])
@pytest.mark.parametrize("libraries", [1, 2])
def test_pipeline_map_failure_stops_the_run(libraries, no_metacool, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the contigs that fails to build stops the run before its
    network, --no-metacool or not, whatever the number of libraries (with
    two, the maps of each library and their sum): an unexpected error (not
    a MetacoolError) is logged as an error, with its type and traceback,
    then raised as it is. The pairs files made from the alignments are
    left in the temporary directory of the run, as a run that stops leaves
    its files."""

    def write_contig_map(*args, **kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr("metator.metacool.write_contig_map", write_contig_map)
    monkeypatch.setattr("metator.metacool.write_library_maps", write_contig_map)
    out = tmp_path / "out"
    forward, reverse = ",".join(["a1.bam", "b1.bam"][:libraries]), ",".join(["a2.bam", "b2.bam"][:libraries])
    args = ["-1", forward, "-2", reverse, "-a", TESTS_DATA / "assembly", "-S", "bam", "-o", out, "-T", tmp_path / "tmp"]
    with pytest.raises(RuntimeError, match=r"^disk full$"):
        run_cli("pipeline", *args, *(["--no-metacool"] if no_metacool else []))
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    error = f"ERROR :: The contact map of the contigs cannot be built at {out / 'metator.metacool'}: RuntimeError: disk full."
    assert f"{error}\nTraceback" in logged
    assert logged.count("ERROR ::") == 1 and "WARNING" not in logged
    assert not {"alignment_to_contacts", "metacool_to_contacts", "partition"} & set(stub_pipeline)
    assert not (out / "metator.metacool").exists() and not (out / "network.txt").exists()
    pairs = [Path(path) for path in stub_pipeline["pairs"]]
    assert [path.name for path in pairs] == [f"alignment_{i}.pairs" for i in range(libraries)]
    assert all(path.exists() and path.parent.parent == tmp_path / "tmp" for path in pairs)


PAIRS_TEXT = "## pairs format v1.0\nr1\tc1\t10\tc2\t20\t+\t-\n"


def _write_pairs(path, codec):
    """A small pairs file at path: plain, compressed with gzip or bzip2, a
    zip archive holding it, empty, or the bytes given."""
    if codec == "plain":
        path.write_text(PAIRS_TEXT)
    elif codec == "empty":
        path.write_bytes(b"")
    elif codec == "gzip":
        path.write_bytes(gzip.compress(PAIRS_TEXT.encode()))
    elif codec == "bzip2":
        path.write_bytes(bz2.compress(PAIRS_TEXT.encode()))
    elif codec == "zip":
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("a.pairs", PAIRS_TEXT)
    else:
        path.write_bytes(codec)
    return path


# The first bytes of a zstd and an lz4 file, which metaTOR does not know.
ZSTD_START, LZ4_START = b"\x28\xb5\x2f\xfd", b"\x04\x22\x4d\x18"

# Pairs files that metacooler or metaTOR, which both find their compression
# from their first bytes, whatever their name, cannot read, given as the
# second library: their name, what they hold, and what the refusal says,
# after their path.
UNREADABLE_PAIRS = {
    "empty, named .gz": ("a.pairs.gz", "empty", "is empty: give the pairs file of its library, or leave it out of -1."),
    "empty": ("a.pairs", "empty", "is empty: give the pairs file of its library, or leave it out of -1."),
    "zip archive": (
        "a.pairs.zip",
        "zip",
        "is a zip archive, which metacooler, which builds the contact map of the contigs, cannot read: decompress "
        "it, and give the pairs file itself, plain or compressed with gzip.",
    ),
    "zip archive, named .gz": (
        "a.pairs.gz",
        "zip",
        "is a zip archive, which metacooler, which builds the contact map of the contigs, cannot read: decompress "
        "it, and give the pairs file itself, plain or compressed with gzip.",
    ),
}


@pytest.mark.parametrize("case", list(UNREADABLE_PAIRS))
def test_pipeline_refuses_pairs_files_metacooler_cannot_read(case, stub_pipeline, tmp_path, run_cli):
    """From pairs, a pairs file that metacooler, which builds the map of the
    contigs, cannot read is refused before any work, with what to do, in the
    log file of the run: the bins of an earlier run stay, even with -F. Each
    file is checked: the faulty one is the second library here."""
    name, codec, message = UNREADABLE_PAIRS[case]
    good = _write_pairs(tmp_path / "good.pairs.gz", "gzip")
    bad = _write_pairs(tmp_path / name, codec)
    out = tmp_path / "out"
    earlier = out / "overlapping_clusters" / "bin.fa"
    earlier.parent.mkdir(parents=True)
    earlier.write_text(">c1\nACGT\n")
    with pytest.raises(ValueError):
        run_cli("pipeline", "-1", f"{good},{bad}", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-F",
                "-T", tmp_path / "tmp")  # fmt: skip
    assert not stub_pipeline and earlier.exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"ERROR :: {bad} {message}\n" in logged
    assert logged.count("ERROR ::") == 1


@pytest.mark.parametrize("name, codec", [("b.pairs.zst", ZSTD_START), ("b.pairs.zstd", ZSTD_START), ("b.pairs.lz4", LZ4_START)])
def test_pipeline_zstd_and_lz4_pairs_of_several_libraries(name, codec, stub_pipeline, tmp_path, run_cli):
    """zstd and lz4 pairs files, which metacooler reads but not metaTOR, go
    on to the maps of the libraries, from several libraries too: their
    networks come from the maps (they were refused, as metaTOR read the
    network of several libraries from their pairs files)."""
    good = _write_pairs(tmp_path / "a.pairs.gz", "gzip")
    other = _write_pairs(tmp_path / name, codec)
    run_cli("pipeline", "-1", f"{good},{other}", "-a", _assembly(tmp_path), "-S", "pair", "-o", tmp_path / "out",
            "-T", tmp_path / "tmp")  # fmt: skip
    assert stub_pipeline["write_library_maps"][0] == [str(good), str(other)]


@pytest.mark.parametrize(
    "name, codec",
    [
        ("a.pairs", "plain"),
        ("a.txt", "plain"),
        ("a.pairs.gz", "gzip"),
        ("a.pairs.bgz", "gzip"),
        ("a.pairs.bz2", "bzip2"),
        ("a.pairs.GZ", "gzip"),
        ("a.pairs.Bz2", "bzip2"),
        ("a.pairs.zst", ZSTD_START),
        ("a.pairs.zstd", ZSTD_START),
        ("a.pairs.LZ4", LZ4_START),
        ("a.pairs", "gzip"),
        ("a.pairs", "bzip2"),
        ("a.pairs.zst", "gzip"),
        ("a.pairs.gz", "plain"),
        ("a.pairs.bz2", "plain"),
        ("a.pairs.GZ", "plain"),
        ("a.pairs", ZSTD_START),
        ("missing.pairs", None),
    ],
    ids=[
        "plain",
        "plain named .txt",
        "gzip",
        "bgzip named .bgz",
        "bzip2",
        "gzip named .GZ",
        "bzip2 named .Bz2",
        "zstd",
        "zstd named .zstd",
        "lz4 named .LZ4",
        "gzip without .gz",
        "bzip2 without .bz2",
        "gzip named .zst",
        "plain named .gz",
        "plain named .bz2",
        "plain named .GZ",
        "zstd without .zst",
        "missing",
    ],
)
def test_pipeline_pairs_files_metacooler_reads(name, codec, stub_pipeline, tmp_path, run_cli):
    """A pairs file that metacooler reads goes on to the map of the
    contigs, whatever its name: metacooler finds its compression from its
    first bytes, as metaTOR does (a gzip file without .gz, a plain file
    named .gz were refused, as metacooler went by the name). zstd and lz4,
    which metacooler reads but not metaTOR, go on too without --scaffold.
    A missing file is left to the map, which stops on it. zstd without .zst
    is a guard: it went on before too."""
    pairs = _write_pairs(tmp_path / name, codec) if codec else tmp_path / name
    run_cli("pipeline", "-1", pairs, "-a", _assembly(tmp_path), "-S", "pair", "-o", tmp_path / "out", "-T", tmp_path / "tmp")
    assert stub_pipeline["write_contig_map"][0] == [str(pairs)]


@pytest.mark.parametrize("libraries", [1, 2])
@pytest.mark.parametrize("name, codec, found", [("lib{i}.pairs.lz4", LZ4_START, "lz4"), ("lib{i}.pairs", ZSTD_START, "zstd")])
def test_pipeline_refuses_zstd_and_lz4_pairs_with_scaffold(libraries, name, codec, found, stub_pipeline, tmp_path, run_cli):
    """With --scaffold, which reads the pairs files, zstd and lz4 pairs
    files are refused before any work, whatever the number of libraries and
    their name: by their first bytes."""
    pairs = [_write_pairs(tmp_path / name.format(i=i), codec) for i in range(libraries)]
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli("pipeline", "-1", ",".join(map(str, pairs)), "-a", _assembly(tmp_path), "-S", "pair", "--scaffold",
                "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert not stub_pipeline
    message = (
        f"{pairs[0]} is compressed with {found}, which metacooler, which builds the contact map of the contigs, reads, "
        "but not metaTOR, which reads the pairs files for --scaffold: decompress it, or compress it with gzip "
        "(zstd and lz4 work without --scaffold)."
    )
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}\n" in log.read_text()


@pytest.mark.parametrize("codec", ["plain", "gzip", "bzip2", ZSTD_START, LZ4_START])
def test_metacooler_pairs_codecs(codec, tmp_path):
    """The compression of a pairs file is found from its first bytes, by
    metacooler, which builds the map from it, as by metaTOR: those that
    metaTOR cannot read (UNREAD_BY_METATOR) are those that metacooler finds
    as zstd and lz4, whatever the name. metacooler went by the name, which
    the previous version of this test pinned."""
    from metacooler.io.create import _pairs_codec

    path = _write_pairs(tmp_path / "a.pairs.gz", codec)
    found = _pairs_codec(path)
    unread = mtcli.UNREAD_BY_METATOR.get(path.read_bytes()[:4])
    assert unread == (found if found in ("zstd", "lz4") else None)
    assert {"gz": "gzip", "bz2": "bz2"}.get(mio.file_compression(str(path))) == (found if found in ("gzip", "bz2") else None)


def test_pipeline_reads_pairs_files_concatenated(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """Pairs files of one assembly concatenated with cat, whose second
    header sits amid the pairs, are one library: metacooler skips the
    header lines below the top, and the map of the contigs holds the pairs
    of both (metacooler 0.4.0; the run stopped on that header line)."""
    monkeypatch.setattr("metator.metacool.write_contig_map", mtmc_write_contig_map)
    _contig_data(monkeypatch, {"c1": (100, 40.5), "c2": (50, 35.0)})
    pairs = _pairs(tmp_path)
    pairs.write_text(pairs.read_text() * 2)
    out = tmp_path / "out"
    run_cli("pipeline", "-1", pairs, "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    import cooler

    # The two pairs of _pairs, twice.
    assert cooler.Cooler(str(out / "metator.metacool")).pixels()[:]["count"].sum() == 4


# What the run says that metacooler needs of the pairs files given to it
# (metator.cli.PAIRS_FILES_NEEDS), as it says it.
NEEDS_OF_PAIRS_FILES = (
    "metacooler, which builds it, needs the same number of fields on every row of each pairs file given to -1 (its "
    "header lines, which start with #, may come anywhere, as in pairs files of one assembly concatenated into one "
    "library), plain or compressed with gzip, bzip2 or, without --scaffold, zstd or lz4. Give the pairs files of "
    "several libraries comma-separated, one file per library."
)


@pytest.mark.parametrize("libraries", [1, 2])
def test_pipeline_says_what_metacooler_needs_of_pairs_files(libraries, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """From pairs, pairs files that metacooler cannot read stop the run on
    the map of the contigs, with the reason that metacooler gives (its CSV
    reader, pyarrow, stops on a row with one more field) and what it needs
    of them, in the log file of the run: with several libraries too, whose
    error also names the library. What it needs is written out here
    (NEEDS_OF_PAIRS_FILES): compared with the constant of metaTOR that says
    it, the log pinned no wording."""
    monkeypatch.setattr("metator.metacool.write_contig_map", mtmc_write_contig_map)
    monkeypatch.setattr("metator.metacool.write_library_maps", mtmc_write_library_maps)
    _contig_data(monkeypatch, {"c1": (100, 40.5), "c2": (50, 35.0)})
    pairs = _pairs(tmp_path)
    text = pairs.read_text()
    pairs.write_text(text.rstrip("\n") + "\tx\n")
    given = [pairs]
    if libraries == 2:
        good = tmp_path / "good.pairs"
        good.write_text(text)
        given = [good, pairs]
    out = tmp_path / "out"
    with pytest.raises(mtmc.MetacoolError) as raised:
        run_cli("pipeline", "-1", ",".join(map(str, given)), "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T",
                tmp_path / "tmp")  # fmt: skip
    if libraries == 2:
        assert str(raised.value).startswith(f"Library 1, {pairs}: ")
    (log,) = out.glob("metator_pipeline_*.log")
    reason = str(raised.value).rstrip(".")
    built = f"ERROR :: The contact map of the contigs cannot be built at {out / 'metator.metacool'}"
    assert f"{built}: {reason}. {NEEDS_OF_PAIRS_FILES}\n" in log.read_text()
    assert not list(out.glob("*.metacool*"))


@pytest.mark.parametrize("start", ["pair", "bam"])
@pytest.mark.parametrize("cause", ["pyarrow", "metaTOR"])
def test_pipeline_says_what_metacooler_needs_only_of_given_pairs(start, cause, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """What metacooler needs of pairs files is said when it cannot read
    those given to a run from pairs: on an error of pyarrow, which it reads
    them with. Not on another error (a guard: a missing pairs file here),
    nor from alignments, whose pairs files the run made itself (a guard)."""
    import pyarrow as pa

    def write_contig_map(*args, **kwargs):
        if cause == "pyarrow":
            raise mtmc.MetacoolError("CSV parse error: Expected 7 columns, got 1") from pa.ArrowInvalid("CSV parse error")
        raise mtmc.MetacoolError("Input file(s) not found: a.pairs.")

    monkeypatch.setattr("metator.metacool.write_contig_map", write_contig_map)
    inputs = ["-1", "a.pairs"] if start == "pair" else ["-1", "a.bam", "-2", "b.bam"]
    out = tmp_path / "out"
    with pytest.raises(mtmc.MetacoolError):
        run_cli("pipeline", *inputs, "-a", _assembly(tmp_path), "-S", start, "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob("metator_pipeline_*.log")
    (error,) = [line.split("ERROR :: ", 1)[1] for line in log.read_text().splitlines() if "ERROR :: " in line]
    assert error.endswith(NEEDS_OF_PAIRS_FILES) == (start == "pair" and cause == "pyarrow")
    assert error.startswith(f"The contact map of the contigs cannot be built at {out / 'metator.metacool'}: ")


NON_ASCII_ASSEMBLY = ">NODE_1 a contig\n" + "ACGT" * 30 + "\n>NODE_522_é\n" + "ACGT" * 30 + "\n"


@pytest.mark.parametrize(
    "start, inputs",
    [
        ("fastq", ["-1", "a.fq", "-2", "b.fq"]),
        ("fastq", ["-1", "a.fq", "-2", "b.fq"]),
        ("bam", ["-1", "a.bam", "-2", "b.bam"]),
        ("pair", ["-1", "a.pairs"]),
    ],
    ids=["fastq", "fastq, with a bowtie2 index", "bam", "pair"],
)
def test_pipeline_refuses_non_ascii_contigs_before_any_alignment(
    start, inputs, request, stub_pipeline, tmp_path, monkeypatch, run_cli
):
    """An assembly with a contig name that is not ASCII is refused as soon
    as its fasta is known, before its index is built or a read aligned,
    whatever the start but a map: the contact map of the contigs needs
    ASCII names, and could only be built once the reads are aligned. The
    message, in the log file of the run, says to rename them. The bins of
    an earlier run stay, even with -F."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(NON_ASCII_ASSEMBLY)
    assembly, named, rebuild = fasta, fasta, ""
    if "index" in request.node.callspec.id:
        # A bowtie2 index, whose fasta is retrieved (by bowtie2-inspect): the
        # message names the index, as given.
        assembly = tmp_path / "index"
        named, rebuild = f"the bowtie2 index {assembly}", ", build the index again"
        monkeypatch.setattr("metator.io.check_fasta_index", lambda ref, mode: str(ref))
        monkeypatch.setattr("metator.io.retrieve_fasta", lambda index, aligner, tmpdir: str(fasta))
    out = tmp_path / "out"
    earlier = out / "overlapping_clusters" / "bin.fa"
    earlier.parent.mkdir(parents=True)
    earlier.write_text(">c1\nACGT\n")
    with pytest.raises(ValueError):
        run_cli("pipeline", *inputs, "-a", assembly, "-S", start, "-o", out, "-F", "-T", tmp_path / "tmp")
    # No index built, no read aligned, no contig data, no map.
    assert not stub_pipeline and earlier.exists()
    (log,) = out.glob("metator_pipeline_*.log")
    message = (
        f"The contig names of {named} cannot be used: A .metacool file needs ASCII contig names (cooler reads them "
        "as ASCII text), and 1 contig name(s) are not: 'NODE_522_é'. Rename these contigs to ASCII names, in the "
        f"assembly and in the alignments or pairs files made from it{rebuild}, then run pipeline again."
    )
    assert f"ERROR :: {message}\n" in log.read_text()


def test_pipeline_from_a_map_keeps_the_contig_names(stub_pipeline, tmp_path, run_cli):
    """A guard: a run from a map does not check the contig names of the
    assembly, as it builds no map of the contigs."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(NON_ASCII_ASSEMBLY)
    metacool = _contig_map(tmp_path / "map.metacool")
    run_cli("pipeline", metacool, "-a", fasta, "-o", tmp_path / "out", "-T", tmp_path / "tmp")
    assert stub_pipeline["metacool_to_contacts"][0] == str(metacool)


@pytest.mark.parametrize(
    "case",
    [
        "pairs",
        "pairs, --no-metacool",
        "link to the map",
        "map a link",
        "reverse pairs",
        "bam",
        "reverse bam",
        "map elsewhere",
        "entry of a collection",
    ],
)
def test_pipeline_refuses_its_map_as_an_input(case, stub_pipeline, tmp_path, run_cli):
    """A run from reads, alignments or pairs refuses, before any work, a
    file given to it (-1, -2: also a -2 that a pair start does not read)
    that is its own contact map (the same path, or a symbolic link to it or
    from it), which it would remove before it builds the map of its contigs
    there, --no-metacool or not; and, from pairs, a contact map given as a
    pairs file, which a run takes as its first argument. The map is left as
    it was."""
    out = tmp_path / "out"
    out.mkdir()
    outfile = out / "metator.metacool"
    forward, reverse, start, options = outfile, None, "pair", []
    if case == "map a link":
        # The map of the run is a link to the map given as input.
        forward = _contig_map(tmp_path / "input.metacool")
        outfile.symlink_to(forward)
    else:
        _contig_map(outfile)
    if case == "pairs, --no-metacool":
        options = ["--no-metacool"]
    elif case == "link to the map":
        forward = tmp_path / "link.metacool"
        forward.symlink_to(outfile)
    elif case == "reverse pairs":
        forward, reverse = _pairs(tmp_path), outfile
    elif case == "bam":
        start, reverse = "bam", "b.bam"
    elif case == "reverse bam":
        start, forward, reverse = "bam", "a.bam", outfile
    elif case == "map elsewhere":
        forward = _contig_map(tmp_path / "elsewhere.metacool")
    elif case == "entry of a collection":
        import metacooler

        collection = tmp_path / "maps.mmetacool"
        metacooler.create_mmetacool(collection, {"lib1": _contig_map(tmp_path / "elsewhere.metacool")})
        forward = f"{collection}::entries/lib1"
    args = ["-1", forward, *(["-2", reverse] if reverse else []), "-a", _assembly(tmp_path), "-S", start, "-o", out,
            "-T", tmp_path / "tmp", *options]  # fmt: skip
    refused = reverse if case.startswith("reverse") else forward
    watched = [outfile, tmp_path / "maps.mmetacool" if case == "entry of a collection" else refused]
    before = {path: Path(path).read_bytes() for path in watched}
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)
    assert {path: Path(path).read_bytes() for path in watched} == before
    assert not stub_pipeline and not (out / "overlapping_clusters").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    if case in ["map elsewhere", "entry of a collection"]:
        message = f"{forward} is a contact map, not a pairs file: give it as the first argument to start from it."
    else:
        message = (
            f"{refused} is the contact map of this run ({outfile}), which the run removes before it builds the map of "
            "its contigs there: give it as the first argument to start from it, or write the run to another output "
            "directory or under another --prefix."
        )
    assert f"ERROR :: {message}" in log.read_text()


def test_pipeline_from_a_pairs_file_named_with_colons(stub_pipeline, tmp_path, run_cli):
    """A pairs file whose name holds :: is a pairs file, not the URI of a
    map inside a file: the run reads it, builds its map from it, and reads
    its network from that map."""
    pairs = tmp_path / "lib::L1.pairs"
    pairs.write_text("## pairs format v1.0\n")
    out = tmp_path / "out"
    run_cli("pipeline", "-1", pairs, "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    assert stub_pipeline["write_contig_map"][0] == [str(pairs)]
    assert stub_pipeline["metacool_to_contacts"][0] == str(out / "metator.metacool")
    assert "reassign_metacool" in stub_pipeline


def test_pipeline_interrupted_while_building_the_map(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """An interruption (Ctrl-C) while the map of the contigs is built stops
    the run, without a word about the map: it is not one of the failures of
    the map, whose reason the run logs as an error."""

    def write_contig_map(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr("metator.metacool.write_contig_map", write_contig_map)
    out = tmp_path / "out"
    result = run_cli("pipeline", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp",
                     exit_code=1)  # fmt: skip
    assert "Aborted" in result.output
    assert "alignment_to_contacts" not in stub_pipeline and "partition" not in stub_pipeline
    (log,) = out.glob("metator_pipeline_*.log")
    assert "cannot be built" not in log.read_text() and "ERROR ::" not in log.read_text()


@pytest.mark.parametrize("case", ["written", "failed", "no metacool", "link"])
def test_pipeline_earlier_metacool(case, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map left by an earlier run in the output directory is removed
    before the map of the contigs is built in its place, --no-metacool or
    not: it would pass for this run's map when the new one cannot be built
    (which stops the run). From then on, the map there is this run's, or
    there is none. A symbolic link there is removed, not its target."""
    out = tmp_path / "out"
    out.mkdir()
    earlier = out / "metator.metacool"
    target = tmp_path / "elsewhere.metacool"
    if case == "link":
        target.write_text("kept")
        earlier.symlink_to(target)
    else:
        earlier.write_text("earlier run")
    seen = []

    def write_contig_map(pairs_files, contig_data, outfile, **kwargs):
        seen.append(os.path.lexists(outfile))
        if case == "failed":
            raise RuntimeError("disk full")
        Path(outfile).write_text("this run")

    monkeypatch.setattr("metator.metacool.write_contig_map", write_contig_map)
    options = ["--no-metacool"] if case == "no metacool" else []
    args = ["-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp", *options]
    if case == "failed":
        with pytest.raises(RuntimeError, match=r"^disk full$"):
            run_cli("pipeline", *args)
    else:
        run_cli("pipeline", *args)
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert seen == [False]
    assert f"INFO :: Removing {earlier}, the contact map of an earlier run." in logged
    if case == "failed":
        assert not os.path.lexists(earlier) and "reassign_metacool" not in stub_pipeline
        assert f"ERROR :: The contact map of the contigs cannot be built at {earlier}: RuntimeError: disk full." in logged
    elif case == "no metacool":
        assert not os.path.lexists(earlier) and "reassign_metacool" not in stub_pipeline
    else:
        assert not earlier.is_symlink() and earlier.read_text() == "this run"
        assert stub_pipeline["reassign_metacool"][0] == str(earlier)
    if case == "link":
        assert target.read_text() == "kept"


def test_pipeline_from_a_map_leaves_another_map_of_its_output_directory(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """From a map stored elsewhere: the run updates that map in place, and a
    map that an earlier run left in the output directory is not this run's
    map: it is left as it is."""
    out = tmp_path / "out"
    out.mkdir()
    earlier = out / "metator.metacool"
    earlier.write_text("earlier run")
    metacool = _contig_map(tmp_path / "input.metacool")
    seen = []

    def reassign_metacool(metacool, contig_data, outfile, **kwargs):
        seen.append(outfile)

    monkeypatch.setattr("metator.metacool.reassign_metacool", reassign_metacool)
    run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert seen == [str(metacool)] and earlier.read_text() == "earlier run"
    assert "write_contig_map" not in stub_pipeline
    (log,) = out.glob("metator_pipeline_*.log")
    assert "Removing" not in log.read_text()


@pytest.mark.parametrize("no_metacool", [False, True])
def test_pipeline_metacool_written(no_metacool, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """With the real writers, from the pairs files and the tables of the
    (stub) validation: the map of the contigs, every contig in
    __unbinned__, is there once the network is built; at the end, the map
    of the bins replaces it or, with --no-metacool, it is removed."""
    monkeypatch.setattr("metator.metacool.write_contig_map", mtmc_write_contig_map)
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _contig_data(monkeypatch, CONTIGS)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    seen = _map_at_partition(monkeypatch, out / "metator.metacool")
    options = ["--no-metacool"] if no_metacool else []
    run_cli("pipeline", *_pair_run(tmp_path, out), *options)
    assert seen == [["__unbinned__"]]
    (log,) = out.glob("metator_pipeline_*.log")
    assert "WARNING" not in log.read_text()
    if no_metacool:
        assert not list(out.glob("*.metacool*"))
        return
    import cooler
    import metacooler

    metacooler.validate(str(out / "metator.metacool"), full=True)
    clr = cooler.Cooler(str(out / "metator.metacool"))
    assert clr.pixels()[:]["count"].sum() == 2
    assert clr.info["genome-assembly"] == "assembly"
    genomes = metacooler.MetaCooler(str(out / "metator.metacool")).genomes
    assert list(genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert list(genomes.Completeness.iloc[1:]) == [0.9]
    # The contig columns of the contig data: its GC_content (Hit is - on
    # every contig), not its ID or Final_cluster.
    contigs = metacooler.MetaCooler(str(out / "metator.metacool")).contigs.set_index("name")
    assert [column for column in contigs.columns if not column.startswith("qc_")] == ["length", "genome_id", "GC_content"]
    assert contigs.GC_content.to_dict() == {"c1": 40.5, "c2": 35.0}
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]


def test_pipeline_map_of_the_clusters_failure_keeps_the_map_of_the_contigs(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the clusters that cannot be written over the map of the
    contigs (here, a contig of another length in the final contig data)
    stops the run, with an error that says why and that the map of the
    contigs stays, as it was, in no cluster. The run's other outputs are
    kept."""
    monkeypatch.setattr("metator.metacool.write_contig_map", mtmc_write_contig_map)
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _contig_data(monkeypatch, {"c1": (120, 40.5), "c2": (50, 35.0)})
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    outfile = out / "metator.metacool"
    built = []
    stub = mtcli.mtp.partition

    def partition(*args):
        built.append(outfile.read_bytes())
        return stub(*args)

    monkeypatch.setattr("metator.partition.partition", partition)
    message = (
        f"{outfile} is not written: {outfile} is not a contact map of the contigs of "
        f"{out / 'contig_data_final.txt'}: 1 contig(s) have another length in it: 'c1' (120 bp in the map, 100 bp in the "
        "contig data). It is still the contact map of the contigs of the run, in no cluster. The other outputs of the "
        "run are kept."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_pair_run(tmp_path, out))
    assert outfile.read_bytes() == built[0]
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    assert (out / "network.txt").exists() and (out / "cluster_summary.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"ERROR :: {message}" in logged
    assert "Traceback" not in logged


@pytest.mark.parametrize("no_metacool", [False, True])
def test_pipeline_failure_leaves_the_map_of_the_contigs(no_metacool, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A run that stops with an error after it built the map of the contigs
    (in the validation, here) leaves it at the outfile, in no bin,
    --no-metacool or not, as it leaves its other files: a new run can start
    from it (the map as its first argument)."""
    monkeypatch.setattr("metator.metacool.write_contig_map", mtmc_write_contig_map)
    _contig_data(monkeypatch, CONTIGS)

    def recursive_decontamination(*args, **kwargs):
        raise RuntimeError("miComplete failed")

    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)
    out = tmp_path / "out"
    options = ["--no-metacool"] if no_metacool else []
    with pytest.raises(RuntimeError, match=r"^miComplete failed$"):
        run_cli("pipeline", *_pair_run(tmp_path, out), *options)
    import metacooler

    outfile = out / "metator.metacool"
    metacooler.validate(str(outfile), full=True)
    assert list(metacooler.MetaCooler(str(outfile)).genomes.name) == ["__unbinned__"]
    assert "reassign_metacool" not in stub_pipeline


@pytest.mark.parametrize("case", ["a contig of another assembly", "no contact"])
def test_pipeline_network_of_the_pairs_of_its_map(case, network_from_a_map, tmp_path, run_cli):
    """The network of a run with one library is that of the pairs that its
    map holds. A pair on a contig that the assembly lacks is left out, as
    metacooler warns (the pairs path stopped on it with a KeyError); a map
    that holds no contact (no pair within the contigs) stops the run with
    the MetacoolError of the network, its reason logged as an error (the
    pairs path divided by zero)."""
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">c1\n" + "ACGT" * 250 + "\n>c2\n" + "GGCC" * 100 + "\n")
    pairs = tmp_path / "a.pairs"
    rows = ["r1\tc1\t10\tc2\t20", "r2\tc1\t10\tzz\t90", "r3\tc1\t30\tc1\t500"]
    if case == "no contact":
        # A read past the end of c2.
        rows = ["r1\tc1\t10\tc2\t900"]
    pairs.write_text("## pairs format v1.0\n" + "".join(f"{row}\t+\t-\n" for row in rows))
    out = tmp_path / "out"
    args = ["-1", pairs, "-a", assembly, "-S", "pair", "-o", out, "-T", tmp_path / "tmp"]
    if case == "no contact":
        message = f"{out / 'metator.metacool'} holds no contact: there is no network to build."
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("pipeline", *args)
        (log,) = out.glob("metator_pipeline_*.log")
        assert f"ERROR :: {message}\n" in log.read_text()
        return
    with pytest.raises(RuntimeError, match=r"^network built$"):
        run_cli("pipeline", *args)
    assert network_from_a_map == [str(out / "metator.metacool")]
    # c1 and c2 (ids 1 and 2), one contact, normalized by their hits (3, 1).
    assert (out / "network.txt").read_text() == f"1\t2\t{1 / 3**0.5}\n"
    (log,) = out.glob("metator_pipeline_*.log")
    assert "WARNING :: Discarded 1 of 3 pair(s) (33.33%): 1 on unrecognised contigs (e.g. ['zz'])" in log.read_text()


def test_pipeline_from_pairs_keeps_them(stub_pipeline, tmp_path, run_cli):
    """From pairs files: they are inputs, kept with their pypairix index.
    Those of two libraries give the map of each library and of their sum,
    from which the networks are read."""
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    pairs = [tmp_path / "a_sorted.pairs.gz", tmp_path / "b_sorted.pairs.gz"]
    for path in pairs:
        # gzip, as their name says (see check_pairs_compression).
        path.write_bytes(gzip.compress(b"## pairs format v1.0\n"))
        Path(f"{path}.px2").write_text("x")
    out = tmp_path / "out"
    run_cli("pipeline", "-1", ",".join(map(str, pairs)), "-a", assembly, "-S", "pair", "-o", out, "-T", tmp_path / "tmp")

    calls = stub_pipeline
    assert "get_contact_pairs" not in calls and "alignment_to_contacts" not in calls
    assert calls["write_library_maps"][0] == [str(path) for path in pairs]
    assert calls["metacool_to_contacts"][0] == str(out / "metator.metacool")
    assert calls["libraries"] == [f"{out / 'metator.mmetacool'}::entries/library_{i}" for i in range(2)]
    assert all(path.exists() and Path(f"{path}.px2").exists() for path in pairs)
    assert not (tmp_path / "tmp").exists()


@pytest.mark.parametrize("no_metacool", [False, True])
def test_pipeline_from_a_metacool(no_metacool, stub_pipeline, tmp_path, run_cli):
    """From a contact map: the network is read from it, no pairs file is
    read or written, and the map of the bins is written back from it. When
    the input map is where the map of the bins goes, it is not removed
    first: it is replaced once the new map is validated."""
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool")
    before = metacool.read_bytes()
    options = ["--no-metacool"] if no_metacool else []
    run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-t", "3", "-m", "abundance",
            "-d", TESTS_DATA / "depth.txt", "-o", out, "-T", tmp_path / "tmp", *options)  # fmt: skip

    calls = stub_pipeline
    assert "get_contact_pairs" not in calls and "alignment_to_contacts" not in calls
    assert "write_metacool" not in calls
    fasta, nb_alignment, depth, enzyme = calls["create_contig_data"]
    assert (fasta, nb_alignment, depth, enzyme) == (str(tmp_path / "assembly.fa"), 1, str(TESTS_DATA / "depth.txt"), None)
    uri, out_dir, network, contig_data, tmp_dir, threads, normalization, self_contacts = calls["metacool_to_contacts"]
    assert (uri, out_dir, network, contig_data) == (str(metacool), str(out), "network.txt", "contig_data_network.txt")
    assert Path(tmp_dir).parent == tmp_path / "tmp"
    assert (threads, normalization, self_contacts) == (3, "abundance", False)
    assert calls["partition"][5] == str(out / "network.txt")
    # The heatmap of the clusters is drawn from the input map.
    assert calls["heatmap_map"] == str(metacool)
    assert metacool.read_bytes() == before
    # The contig data of the network step is kept, as from pairs files. A
    # map sums its libraries: no network_0.txt, even when the input map is
    # at the path of the run's own map.
    assert (out / "contig_data_network.txt").exists() and not (out / "network_0.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert "Removing" not in logged
    if no_metacool:
        assert "reassign_metacool" not in calls and "is the input map" not in logged
        return
    source, existed, contig_data, outfile, kwargs = calls["reassign_metacool"]
    assert (source, existed, contig_data, outfile) == (str(metacool), True, str(out / "contig_data_final.txt"), str(metacool))
    assert (kwargs["bin_summary"], kwargs["force"]) == (str(out / "cluster_summary.txt"), True)
    parameters = kwargs["metadata"]["parameters"]
    assert parameters["map"] == str(metacool) and parameters["start"] is None
    assert kwargs["metadata"]["command"].startswith(f"metator pipeline {metacool} -a")
    # No column was read from the map (the stub reads none).
    assert "columns_from_map" not in kwargs["metadata"]
    assert f"INFO :: {metacool} is the input map {metacool}: it is replaced once the map of the clusters is validated." in logged


def _contig_map(path, lengths=None):
    """A .metacool map of the contigs c1 and c2 of the stub validation's
    tables (see _final_tables), made before any binning: every contig in
    __unbinned__, one bin per contig."""
    import metacooler
    import pandas as pd

    lengths = lengths or {"c1": 100, "c2": 50}
    pairs = path.with_suffix(".pairs")
    pairs.write_text("## pairs format v1.0\nr1\tc1\t10\tc2\t20\t+\t-\nr2\tc1\t10\tc1\t90\t+\t-\n")
    bins = pd.DataFrame({"chrom": list(lengths), "start": 0, "end": list(lengths.values())})
    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    pairs.unlink()
    return path


# The restriction sites and the coverage of the contigs, read from the
# Restriction_site and Shotgun_coverage columns of the map that a run starts
# from, where --enzyme and --depth measure none.

mtn_contig_data_from_map = __import__("metator.network", fromlist=["contig_data_from_map"]).contig_data_from_map


def _map_with_columns(tmp_path, columns):
    """_contig_map, with columns (by column, the value of c1 and c2) in its
    contigs table."""
    import metacooler

    metacool = _contig_map(tmp_path / "map.metacool")
    if not columns:
        return metacool
    table = tmp_path / "columns.txt"
    rows = [f"{name}\t" + "\t".join(str(values[i]) for values in columns.values()) for i, name in enumerate(["c1", "c2"])]
    table.write_text("\t".join(["Name", *columns]) + "\n" + "".join(f"{row}\n" for row in rows))
    metacooler.annotate(str(metacool), contigs_table=str(table), contigs_table_format="Name", inplace=True)
    return metacool


SITES_AND_COVERAGE = {"Restriction_site": [3, 2], "Shotgun_coverage": [1.5, 0.5]}
RS_ERROR = 'For "RS" and "theoritical_hit" normalization, enzyme is required'
DEPTH_ERROR = 'For "abundance" and "theoritical_hit" normalization, depth is required'


@pytest.mark.parametrize(
    "normalization, options, columns, message",
    [
        ("RS", {}, {}, f"{RS_ERROR}, as {{map}} gives no restriction sites (no Restriction_site column with a valid value)."),
        ("abundance", {}, {"Restriction_site": [3, 2]},
         f"{DEPTH_ERROR}, as {{map}} gives no shotgun coverage (no Shotgun_coverage column with a valid value)."),
        ("RS", {}, {"Restriction_site": ["-", "-"]},
         f"{RS_ERROR}, as {{map}} gives no restriction sites (no Restriction_site column with a valid value)."),
        ("theoritical_hit", {}, SITES_AND_COVERAGE, None),
        ("RS", {}, {"Restriction_site": [3, 2]}, None),
        ("abundance", {}, {"Shotgun_coverage": [0, 2]}, None),
        ("theoritical_hit", {"depth": "depth.txt"}, {"Restriction_site": [3, 2]}, None),
        ("empirical_hit", {}, {}, None),
    ],
    ids=["map without columns", "map without coverage", "column of -", "map columns",
         "map sites", "map coverage", "map sites and depth", "no need"],
)  # fmt: skip
def test_check_normalization(normalization, options, columns, message, tmp_path, logged):
    """Before any work, from a map: --enzyme and --depth, or the columns of
    the map that give values."""
    metacool = _map_with_columns(tmp_path, columns)
    args = (normalization, options.get("enzyme"), options.get("depth"), str(metacool))
    if message is None:
        mtcli.check_normalization(*args)
        assert logged() == []
        return
    with pytest.raises(ValueError):
        mtcli.check_normalization(*args)
    assert logged() == [("ERROR", message.format(map=metacool))]


@pytest.mark.parametrize(
    "normalization, options, message",
    [("RS", [], f"{RS_ERROR}."), ("theoritical_hit", ["-e", "DpnII"], f"{DEPTH_ERROR}.")],
    ids=["RS", "depth"],
)
def test_normalization_refused_from_pairs(normalization, options, message, stub_pipeline, tmp_path, run_cli):
    """A guard: a run that does not start from a map is refused before any
    work, with the messages of before, which name the options only."""
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli("pipeline", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-m", normalization, *options,
                "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert "create_contig_data" not in stub_pipeline
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}\n" in log.read_text()


def test_check_normalization_once_the_map_is_read(tmp_path, logged):
    """Once the run has read the columns of its map (read): those that it
    could read, or a refusal for a column that lacks contigs."""
    metacool = str(tmp_path / "map.metacool")
    mtcli.check_normalization("theoritical_hit", None, None, metacool, ["Restriction_site", "Shotgun_coverage"])
    mtcli.check_normalization("RS", None, "depth.txt", metacool, ["Restriction_site"])
    with pytest.raises(ValueError):
        mtcli.check_normalization("theoritical_hit", None, None, metacool, ["Restriction_site"])
    assert logged() == [
        (
            "ERROR",
            f"{DEPTH_ERROR}, as the Shotgun_coverage column of {metacool} does not give every contig of the assembly its "
            "shotgun coverage.",
        )
    ]


@pytest.mark.parametrize("normalization", ["RS", "empirical_hit"])
def test_check_normalization_of_an_unreadable_map(normalization, tmp_path, logged):
    """A map that cannot be read is refused before any work when a
    normalization needs its columns, its reason logged; otherwise the run
    finds it at its network, as before."""
    metacool = str(tmp_path / "missing.metacool")
    if normalization == "empirical_hit":
        mtcli.check_normalization(normalization, None, None, metacool)
        assert logged() == []
        return
    with pytest.raises(mtmc.MetacoolError):
        mtcli.check_normalization(normalization, None, None, metacool)
    ((level, message),) = logged()
    assert level == "ERROR" and message.startswith(f"Cannot read the contact map {metacool}")


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_normalization_refused_before_any_work_from_a_map(command, stub_pipeline, tmp_path, run_cli):
    """A run from a map whose columns cannot give what its normalization
    needs is refused before the assembly is read."""
    metacool = _contig_map(tmp_path / "map.metacool")
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli(command, metacool, "-a", _assembly(tmp_path), "-n" if command == "network" else "-m",
                "abundance", "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert "create_contig_data" not in stub_pipeline
    (log,) = out.glob(f"metator_{command}_*.log")
    message = f"{DEPTH_ERROR}, as {metacool} gives no shotgun coverage (no Shotgun_coverage column with a valid value)."
    assert f"ERROR :: {message}\n" in log.read_text()


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_run_reads_the_columns_of_its_map(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A run from a map reads its columns once it has the contigs of the
    assembly (theoritical_hit needs both), and normalizes with them."""
    monkeypatch.setattr("metator.network.contig_data_from_map", mtn_contig_data_from_map)
    contig_data = {name: {"id": i, "length": length, "GC": 50.0, "hit": 0, "coverage": "-", "RS": "-"}
                   for i, (name, length) in enumerate([("c1", 100), ("c2", 50)], start=1)}  # fmt: skip
    monkeypatch.setattr("metator.network.create_contig_data", lambda *args: (contig_data, None))
    metacool = _map_with_columns(tmp_path, SITES_AND_COVERAGE)
    out = tmp_path / "out"
    run_cli(command, metacool, "-a", _assembly(tmp_path), "-n" if command == "network" else "-m",
            "theoritical_hit", "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert {name: (data["RS"], data["coverage"]) for name, data in contig_data.items()} == {"c1": (3, 1.5), "c2": (2, 0.5)}
    assert stub_pipeline["metacool_to_contacts"][-2:] == ("theoritical_hit", False)
    if command == "pipeline":
        # Where the values that the map of the bins stores came from.
        metadata = stub_pipeline["reassign_metacool"][4]["metadata"]
        assert metadata["columns_from_map"] == ["Restriction_site", "Shotgun_coverage"]


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_run_refused_for_a_column_that_lacks_contigs(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A column of the map that lacks a contig of the assembly cannot give
    the normalization what it needs: refused, with the warning that names
    the contig, before the network."""
    monkeypatch.setattr("metator.network.contig_data_from_map", mtn_contig_data_from_map)
    contig_data = {name: {"id": i, "length": length, "GC": 50.0, "hit": 0, "coverage": "-", "RS": "-"}
                   for i, (name, length) in enumerate([("c1", 100), ("c2", 50), ("c3", 30)], start=1)}  # fmt: skip
    monkeypatch.setattr("metator.network.create_contig_data", lambda *args: (contig_data, None))
    metacool = _map_with_columns(tmp_path, {"Restriction_site": [3, 2]})
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli(command, metacool, "-a", _assembly(tmp_path), "-n" if command == "network" else "-m",
                "RS", "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert "metacool_to_contacts" not in stub_pipeline
    (log,) = out.glob(f"metator_{command}_*.log")
    logged = log.read_text()
    assert "WARNING :: The Restriction_site column of" in logged and "contig(s) of the assembly" in logged
    message = f"{RS_ERROR}, as the Restriction_site column of {metacool} does not give every contig of the assembly its "
    assert f"ERROR :: {message}restriction sites.\n" in logged


def test_pipeline_refused_for_a_column_that_lacks_contigs_keeps_the_clusters(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A run from a map with -F (a rebinning in the output directory of an
    earlier run) that its columns refuse once the assembly is read leaves
    the clusters of the earlier run, under their names and under those of
    metaTOR 1.x, which -F removes too: the contig data is read before they
    are replaced (they were removed, then the run refused)."""
    monkeypatch.setattr("metator.network.contig_data_from_map", mtn_contig_data_from_map)
    contig_data = {name: {"id": i, "length": length, "GC": 50.0, "hit": 0, "coverage": "-", "RS": "-"}
                   for i, (name, length) in enumerate([("c1", 100), ("c2", 50), ("c3", 30)], start=1)}  # fmt: skip
    monkeypatch.setattr("metator.network.create_contig_data", lambda *args: (contig_data, None))
    metacool = _map_with_columns(tmp_path, {"Restriction_site": [3, 2]})
    out = tmp_path / "out"
    names = ["overlapping_clusters", "recursive_clusters", "final_clusters", "overlapping_bin", "recursive_bin",
             "final_bin_unscaffold"]  # fmt: skip
    clusters = [out / name / "cluster.fa" for name in names]
    for path in clusters:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    with pytest.raises(ValueError):
        run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-m", "RS", "-o", out, "-F",
                "-T", tmp_path / "tmp")  # fmt: skip
    assert all(path.exists() for path in clusters)


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_run_from_a_map_unreadable_for_its_columns(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map that cannot be read when the run reads its columns stops it,
    the reason logged."""

    def unreadable(metacool, *args):
        raise mtmc.MetacoolError(f"Cannot read the contact map {metacool}: gone.")

    monkeypatch.setattr("metator.network.contig_data_from_map", unreadable)
    metacool = _contig_map(tmp_path / "map.metacool")
    out = tmp_path / "out"
    with pytest.raises(mtmc.MetacoolError):
        run_cli(command, metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: Cannot read the contact map {metacool}: gone.\n" in log.read_text()


def _collection(path, *entries):
    """A .mmetacool collection at path, with a map of the contigs
    (_contig_map) as each entry of entries."""
    import metacooler

    maps = {name: _contig_map(path.with_name(f"{path.stem}_{name}.metacool")) for name in entries}
    metacooler.create_mmetacool(path, maps)
    return path


@pytest.mark.parametrize("where", ["elsewhere", "outdir", "output-map"])
def test_pipeline_metacool_written_from_a_map(where, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """With the real writer, from the tables of the (stub) validation: the
    map of the clusters, from the input map, which it replaces wherever it
    is; with --output-map, written to that file, the input map left as it
    was. No other map goes to the output directory."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map((out if where == "outdir" else tmp_path) / "metator.metacool")
    before = metacool.read_bytes()
    output_map = tmp_path / "copy.metacool"
    options = ["--output-map", output_map] if where == "output-map" else []
    run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp",
            *options)  # fmt: skip
    import cooler
    import metacooler

    outfile = output_map if where == "output-map" else metacool
    metacooler.validate(str(outfile), full=True)
    assert cooler.Cooler(str(outfile)).pixels()[:]["count"].sum() == 2
    mc = metacooler.MetaCooler(str(outfile))
    assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert list(mc.genomes.Completeness.iloc[1:]) == [0.9]
    assert dict(zip(mc.contigs.name, mc.contigs.GC_content, strict=True)) == {"c1": 40.5, "c2": 35.0}
    assert (metacool.read_bytes() == before) == (where == "output-map")
    # The heatmap is drawn from the input map, wherever the map of the
    # clusters goes: --output-map is not written yet.
    assert stub_pipeline["heatmap_map"] == str(metacool)
    expected = ["metator.metacool"] if where == "outdir" else []
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == expected
    assert not list(tmp_path.rglob("*.lock"))
    (log,) = out.glob("metator_pipeline_*.log")
    assert "is not written" not in log.read_text()


def test_pipeline_metacool_from_a_map_failure_keeps_it(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the clusters that cannot be written from the input map
    (here, a contig of another length in the contig data) stops the run,
    with an error that says why and that the input map is left as it was.
    The run's other outputs are kept."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool", {"c1": 120, "c2": 50})
    before = metacool.read_bytes()
    message = (
        f"{metacool} is not written: {metacool} is not a contact map of the contigs of "
        f"{out / 'contig_data_final.txt'}: 1 contig(s) have another length in it: 'c1' (120 bp in the map, 100 bp in the "
        "contig data). It is still the input map of the run, as it was. The other outputs of the run are kept."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    assert (out / "network.txt").exists() and (out / "cluster_summary.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"ERROR :: {message}" in logged
    assert "Traceback" not in logged


@pytest.mark.parametrize(
    "error", [ValueError("pixels out of order"), RuntimeError("unexpected failure")], ids=["ValueError", "RuntimeError"]
)
def test_pipeline_metacool_from_a_map_metacooler_failure(error, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """metacooler fails while it writes the map of the clusters over the
    input map (its validation, here): it removes its partial file and says
    so at INFO, and the run stops with an error that the map is not written,
    with the reason of one of the errors that metacooler is expected to
    raise (a ValueError here), or the traceback of another (a RuntimeError).
    The input map is left as it was, and no other file is left beside it."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool")
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise error

    monkeypatch.setattr("metacooler.io.reassign.validate", validate)
    with pytest.raises(mtmc.MetacoolError, match=r" is not written: "):
        run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    (staged,) = re.findall(r"INFO :: Creation failed; the partial file (\S+)\.partial was removed and \1 was not written", logged)
    assert Path(staged).parent == out and Path(staged).name.startswith(".metator.metacool.")
    assert "for inspection" not in logged
    kept = "It is still the input map of the run, as it was. The other outputs of the run are kept."
    if isinstance(error, ValueError):
        assert f"ERROR :: {metacool} is not written: pixels out of order. {kept}\n" in logged
        assert "Traceback" not in logged
    else:
        assert f"ERROR :: {metacool} is not written: RuntimeError: unexpected failure. {kept}\nTraceback" in logged
    assert logged.count("ERROR ::") == 1


def _map_run(metacool, tmp_path, *options):
    """The arguments of a pipeline run from the map metacool."""
    out = tmp_path / "out"
    return [metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp", *options]


def _lock(metacool, pid, host=None, started="2026-11-18T10:00:00+0000"):
    """A lock of the map metacool, as a run of pid on host holds it."""
    import socket

    holder = {"pid": pid, "host": host or socket.gethostname(), "started": started, "command": "metator pipeline"}
    lock = Path(f"{metacool}.metator.lock")
    lock.write_text(json.dumps(holder))
    return lock


def _rounds_tables(monkeypatch, members):
    """Have the stub validation also write the clustering rounds of the run
    (see metator.clusters): metator_00001_00000, an overlapping cluster of
    the contigs members, which is also the final cluster of c1."""
    import metator.clusters as mtcl

    stub = mtcli.mtv.recursive_decontamination

    def recursive_decontamination(*args, **kwargs):
        result = stub(*args, **kwargs)
        rounds = mtcl.Rounds()
        rounds.add(mtcl.OVERLAPPING, "metator_00001_00000", members, columns={"Completeness": 0.9, "Lineage": "x"})
        rounds.close_round(mtcl.OVERLAPPING)
        rounds.add(mtcl.FINAL, "metator_00001_00000", ["c1"], parent="metator_00001_00000", columns={"step": 0})
        rounds.write(args[8])
        return result

    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)


def test_pipeline_from_a_map_stores_its_rounds(stub_pipeline, tmp_path, monkeypatch, run_cli, caplog):
    """The map of the clusters holds the clustering rounds of the run in its
    /clusters group: its clusters, their parents and columns, and their
    contigs, those that the map lacks left out with a warning."""
    import metacooler

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    _rounds_tables(monkeypatch, ["c1", "c2", "c9"])
    metacool = _contig_map(tmp_path / "input.metacool")
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    metacooler.validate(str(metacool), full=True)
    mc = metacooler.MetaCooler(str(metacool))
    rounds = mc.cluster_rounds
    assert [entry["name"] for entry in rounds] == ["overlapping", "final"]
    assert {entry["tool"] for entry in rounds} == {"metaTOR"} and rounds[0]["step"] == "partition"
    assert rounds[1] == {"name": "final", "tool": "metaTOR", "version": __version__, "step": "validation"}
    clusters = mc.clusters
    assert list(clusters.parent_round.fillna("-")) == ["-", "overlapping"]
    assert list(clusters.Completeness.iloc[:1]) == [0.9] and list(clusters.Lineage) == ["x", ""]
    members = mc.cluster_members()
    assert sorted(zip(members["round"], members.contig, strict=True)) == [
        ("final", "c1"),
        ("overlapping", "c1"),
        ("overlapping", "c2"),
    ]
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert "1 contig(s) of the clusters of the run are not in" in log.read_text()


def test_pipeline_from_a_map_without_rounds_drops_those_it_held(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The rounds of an earlier run describe the clusters it replaced: a run
    that wrote none (here, the stub validation) leaves the map with none."""
    import metacooler
    import pandas as pd

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    metacooler.write_clusters(
        str(metacool),
        ["earlier"],
        pd.DataFrame({"round": ["earlier"], "name": ["x"]}),
        pd.DataFrame({"round": ["earlier"], "cluster": ["x"], "contig": ["c2"]}),
    )
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert metacooler.MetaCooler(str(metacool)).cluster_rounds == []


def _map_of_an_earlier_run(path):
    """A map of the contigs c1 and c2 of _contig_map as an earlier run left
    it: c1 in the cluster earlier_00001, c2 in __unbinned__, and the
    clustering rounds overlapping (made by partition), final (by
    validation) and mge (by mge)."""
    import metacooler
    import pandas as pd

    pairs = path.with_suffix(".pairs")
    pairs.write_text("## pairs format v1.0\nr1\tc1\t10\tc2\t20\t+\t-\nr2\tc1\t10\tc1\t90\t+\t-\n")
    bins = pd.DataFrame({"chrom": ["c1", "c2"], "start": 0, "end": [100, 50]})
    assignment = {"c1": "earlier_00001", "c2": "ND"}
    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=assignment, ps=False)
    pairs.unlink()
    rounds = [("overlapping", "partition", "c1"), ("final", "validation", "c1"), ("mge", "mge", "c2")]
    metacooler.write_clusters(
        str(path),
        [{"name": name, "tool": "metaTOR", "step": step} for name, step, _ in rounds],
        pd.DataFrame({"round": [name for name, _, _ in rounds], "name": ["o1", "earlier_00001", "m1"]}),
        pd.DataFrame(
            {
                "round": [name for name, _, _ in rounds],
                "cluster": ["o1", "earlier_00001", "m1"],
                "contig": [c for *_, c in rounds],
            }
        ),
    )
    return path


@pytest.mark.parametrize(
    "command, replaced",
    [
        ("pipeline", "1 cluster(s) and the clustering round(s) overlapping, final"),
        ("partition", "the clustering round(s) overlapping, final"),
        ("validation", "1 cluster(s) and the clustering round(s) final"),
        ("mge", "the clustering round(s) mge"),
    ],
)
def test_warn_of_what_a_run_replaces_in_a_map(command, replaced, tmp_path, logged):
    """What a run that updates a map in place replaces there, of an earlier
    run (or another binner): pipeline MAP its clusters and its rounds of
    the binning; partition MAP its rounds of the binning; validation MAP its
    clusters and its rounds of validation, as partition MAP, step by step,
    stores the round overlapping just before; mge --map its round mge. A map
    of no earlier run gives no warning."""
    metacool = _map_of_an_earlier_run(tmp_path / "input.metacool")
    mtcli.warn_of_what_is_replaced(command, str(metacool))
    mtcli.warn_of_what_is_replaced(command, str(_contig_map(tmp_path / "fresh.metacool")))
    message = (
        f"{metacool} holds {replaced}, which this run replaces, in place: give --output-map to write the map of this "
        f"run to another file, and leave {metacool} as it is."
    )
    assert [record for record in logged() if record[0] == "WARNING"] == [("WARNING", message)]


@pytest.mark.parametrize("output", ["in place", "--output-map"])
def test_pipeline_from_a_map_warns_of_what_it_replaces(output, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """pipeline MAP, in place, warns before any work that it replaces the
    clusters and the rounds that an earlier run left in the map: it
    replaced them without a word, and without -F (the map keeps no trace of
    them). With --output-map, which leaves the map as it was, it does not:
    a guard, as it gave no warning before either."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _map_of_an_earlier_run(tmp_path / "input.metacool")
    options = ["--output-map", tmp_path / "copy.metacool"] if output == "--output-map" else []
    run_cli("pipeline", *_map_run(metacool, tmp_path, *options))
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    replaced = (
        f"WARNING :: {metacool} holds 1 cluster(s) and the clustering round(s) overlapping, final, which this run "
        f"replaces, in place: give --output-map to write the map of this run to another file, and leave {metacool} "
        "as it is."
    )
    assert (replaced in log.read_text()) == (output == "in place")


@pytest.mark.parametrize("command", ["partition", "validation", "mge"])
def test_a_step_warns_of_what_it_replaces_in_place(command, tmp_path, logged):
    """partition MAP, validation MAP and mge --map MAP warn, as they take the
    lock of the map, of what they replace in it in place (see
    warn_of_what_is_replaced); with --output-map, they do not."""
    from contextlib import ExitStack

    metacool = _map_of_an_earlier_run(tmp_path / "input.metacool")
    option = "--map" if command == "mge" else None
    for output_map in [str(tmp_path / "copy.metacool"), None]:
        with ExitStack() as held:
            mtcli.step_map(command, str(metacool), output_map, False, held, option=option)
        warned = [message for level, message in logged() if level == "WARNING"]
        assert len(warned) == (output_map is None) and all(message.startswith(f"{metacool} holds ") for message in warned)


def test_pipeline_from_a_map_to_an_existing_output_map(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """--output-map replaces a file only with -F: refused before any work,
    else the map of the clusters replaces it, the input map left as it was."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    output_map = tmp_path / "copy.metacool"
    output_map.write_text("an older map")
    message = f"{output_map} exists: give -F to replace it with the map of this run, or another --output-map."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", output_map))
    assert "metacool_to_contacts" not in stub_pipeline and output_map.read_text() == "an older map"
    run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", output_map, "-F"))
    import metacooler

    metacooler.validate(str(output_map), full=True)
    assert metacool.read_bytes() == before


def test_pipeline_from_a_map_to_itself_as_output_map(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """--output-map naming the input map is the default: in place, no -F."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", metacool))
    assert metacool.read_bytes() != before


def test_pipeline_from_an_entry_needs_an_output_map(stub_pipeline, tmp_path, run_cli):
    """An entry of a .mmetacool collection is not updated in place: that
    would rewrite the collection. Refused before any work."""
    entry = f"{_collection(tmp_path / 'maps.mmetacool', 'rep1')}::entries/rep1"
    message = (
        f"{entry} is an entry of a .mmetacool collection, which a run does not update in place: write the map of its "
        "clusters to another file with --output-map."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(entry, tmp_path))
    assert "metacool_to_contacts" not in stub_pipeline


def test_pipeline_from_a_map_through_a_symbolic_link(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map given through a symbolic link: the link is replaced by the
    map of the clusters, and the map it pointed to is left as it was."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    target = _contig_map(tmp_path / "target.metacool")
    before = target.read_bytes()
    link = tmp_path / "link.metacool"
    link.symlink_to(target)
    run_cli("pipeline", *_map_run(link, tmp_path))
    import metacooler

    assert not link.is_symlink() and target.read_bytes() == before
    assert list(metacooler.MetaCooler(str(link)).genomes.name) == ["__unbinned__", "metator_00001_00000"]


@pytest.mark.parametrize("failure", ["partition", "Ctrl-C in the validation", "Ctrl-C in metacooler reassign"])
def test_pipeline_from_a_map_that_fails_leaves_it_as_it_was(failure, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A run that stops (an error mid-run, Ctrl-C, even while metacooler
    writes the new map) leaves the input map byte for byte as it was, with
    no file beside it: no staged map, no lock.

    A guard: a run from a map wrote a new map in its output directory, and
    left the input map alone, so this passes without the update in place
    too; it holds the update to replacing the map only once it is written.
    """
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    (tmp_path / "maps").mkdir()
    metacool = _contig_map(tmp_path / "maps" / "input.metacool")
    before = metacool.read_bytes()
    error = RuntimeError("partition failed") if failure == "partition" else KeyboardInterrupt()

    def stopped(*args, **kwargs):
        raise error

    target = {
        "partition": "metator.partition.partition",
        "Ctrl-C in the validation": "metator.validation.recursive_decontamination",
        "Ctrl-C in metacooler reassign": "metacooler.reassign",
    }[failure]
    monkeypatch.setattr(target, stopped)
    if failure == "partition":
        with pytest.raises(RuntimeError, match=r"^partition failed$"):
            run_cli("pipeline", *_map_run(metacool, tmp_path))
    else:
        assert "Aborted" in run_cli("pipeline", *_map_run(metacool, tmp_path), exit_code=1).output
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in (tmp_path / "maps").iterdir()) == ["input.metacool"]


def test_pipeline_from_a_map_refused_while_another_run_holds_it(stub_pipeline, tmp_path, run_cli):
    """A second run on a map that a running run holds (here, this process)
    is refused before any work, and leaves the lock to its holder."""
    metacool = _contig_map(tmp_path / "input.metacool")
    lock = _lock(metacool, os.getpid())
    message = (
        f"{metacool} is being written by another metaTOR run (pid {os.getpid()} on {json.loads(lock.read_text())['host']}, "
        f"since 2026-11-18T10:00:00+0000): wait for it to end, write this run's map elsewhere with --output-map, or, if "
        f"no such run is going on, remove {lock}."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert "metacool_to_contacts" not in stub_pipeline and lock.exists()
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("holder", ["another host", "unreadable", "not an object"])
def test_pipeline_from_a_map_refused_with_a_lock_it_cannot_check(holder, stub_pipeline, tmp_path, run_cli):
    """A lock of a run on another host, or one that says nothing readable
    (not JSON, or JSON but not an object, which raised AttributeError),
    cannot be told stale: the run is refused, and the lock left."""
    metacool = _contig_map(tmp_path / "input.metacool")
    if holder != "another host":
        lock = Path(f"{metacool}.metator.lock")
        lock.write_text("{not json" if holder == "unreadable" else "[1]")
        who = "a run that cannot be told"
    else:
        lock = _lock(metacool, 1, host="elsewhere")
        who = "pid 1 on elsewhere, since 2026-11-18T10:00:00+0000"
    with pytest.raises(mtmc.MetacoolError, match=re.escape(f"another metaTOR run ({who})")):
        run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert lock.exists()


def test_pipeline_from_a_map_whose_lock_cannot_be_created(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map whose directory is not writable cannot be locked, nor updated in
    place: the run is refused before any work, with what to do (it stopped
    on a bare PermissionError)."""
    metacool = _contig_map(tmp_path / "input.metacool")
    lock = f"{metacool}.metator.lock"
    opened = os.open

    def read_only(path, *args, **kwargs):
        if str(path) == lock:
            raise PermissionError(13, "Permission denied", path)
        return opened(path, *args, **kwargs)

    monkeypatch.setattr(mtmc.os, "open", read_only)
    message = (
        f"{lock} cannot be created (Permission denied): a run that writes {metacool} locks it there, so the directory "
        "of the map must exist and be writable; write this run's map elsewhere with --output-map."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert "metacool_to_contacts" not in stub_pipeline
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}" in log.read_text()


def test_the_lock_of_a_map_given_as_its_file_root(tmp_path):
    """A map given as file::/ (the file itself) is locked as the file: next to
    it, not at file::/.metator.lock."""
    metacool = _contig_map(tmp_path / "input.metacool")
    with mtmc.map_lock(f"{metacool}::/") as lock:
        assert lock == f"{metacool}.metator.lock" and Path(lock).exists()
    assert not Path(lock).exists()


def test_a_run_removes_its_own_lock_only(tmp_path):
    """A lock that another run took during this one (this run's removed by
    hand) is that run's: it stays when this run ends."""
    metacool = _contig_map(tmp_path / "input.metacool")
    with mtmc.map_lock(metacool) as lock:
        other = _lock(metacool, 1, host="elsewhere")
    assert other.exists() and Path(lock) == other


def test_a_stale_lock_is_taken_once(tmp_path):
    """Two runs that both found a lock stale: the first moves it aside and
    takes the lock; the second moves the first's lock aside, sees that it is
    not the stale one, puts it back, and is refused (both ran, each having
    removed the lock that the other had just taken)."""
    import socket

    metacool = _contig_map(tmp_path / "input.metacool")
    stale = json.loads(_lock(metacool, 999999999).read_text())
    first = {"pid": 1, "host": socket.gethostname(), "started": "now", "run": "first"}
    second = {"pid": 2, "host": socket.gethostname(), "started": "now", "run": "second"}
    lock = f"{metacool}.metator.lock"
    assert mtmc._take_stale_lock(lock, stale, first)
    assert not mtmc._take_stale_lock(lock, stale, second)
    assert json.loads(Path(lock).read_text()) == first
    assert [path.name for path in tmp_path.iterdir() if ".lock." in path.name] == []
    # Gone in the meantime (taken and released by another run): taken.
    Path(lock).unlink()
    assert mtmc._take_stale_lock(lock, stale, second)
    assert json.loads(Path(lock).read_text()) == second


def test_a_stale_lock_taken_by_a_third_run(tmp_path, monkeypatch):
    """A guard of the case that _take_stale_lock leaves open: a second run
    moves aside the lock that a first run took, and a third run creates the
    lock in the meantime. The second cannot put the first's back: it leaves the third's,
    and is refused (the first and the third then both run, the last to end
    replacing the map)."""
    import socket

    metacool = _contig_map(tmp_path / "input.metacool")
    stale = json.loads(_lock(metacool, 999999999).read_text())
    runs = {name: {"pid": i, "host": socket.gethostname(), "started": "now", "run": name}
            for i, name in enumerate(["first", "second", "third"], start=1)}  # fmt: skip
    lock = f"{metacool}.metator.lock"
    assert mtmc._take_stale_lock(lock, stale, runs["first"])
    rename = os.rename

    def then_a_third_run(source, target):
        rename(source, target)
        assert mtmc._create_lock(lock, runs["third"])

    monkeypatch.setattr(mtmc.os, "rename", then_a_third_run)
    assert not mtmc._take_stale_lock(lock, stale, runs["second"])
    assert json.loads(Path(lock).read_text()) == runs["third"]
    assert [path.name for path in tmp_path.iterdir() if ".lock." in path.name] == []


@pytest.mark.parametrize(
    "case",
    ["missing map", "no map written", "an entry", "a collection", "a directory", "the collection of the input"],
)
def test_pipeline_output_map_refusals(case, stub_pipeline, tmp_path, run_cli):
    """What --output-map cannot be, and an input map that is not there, are
    refused before any work (they failed at the end of the run, or on a
    bare error while locking). A map named as a .mmetacool collection is
    refused too: with several libraries, the collection of their maps
    would be the map itself (see
    test_pipeline_output_map_and_the_maps_of_its_libraries)."""
    metacool = _contig_map(tmp_path / "input.metacool")
    options, given = [], metacool
    if case == "missing map":
        given = tmp_path / "missing" / "map.metacool"
        message = f"Cannot read the contact map {given}: no such file."
    elif case == "no map written":
        output = tmp_path / "copy.metacool"
        options = ["--no-metacool", "--output-map", output]
        message = f"--no-metacool writes no map: drop --output-map {output}, or --no-metacool."
    elif case == "an entry":
        output = f"{tmp_path / 'maps.mmetacool'}::entries/rep2"
        options = ["--output-map", output]
        message = f"--output-map takes the path of a .metacool file, not an entry of a collection: {output}."
    elif case == "a collection":
        output = tmp_path / "copy.mmetacool"
        options = ["--output-map", output]
        message = f"--output-map takes the path of a .metacool file, not of a .mmetacool collection: {output}."
    elif case == "a directory":
        output = tmp_path / "maps"
        output.mkdir()
        options = ["--output-map", output, "-F"]
        message = f"{output} is a directory: --output-map takes the path of a .metacool file."
    else:
        collection = _collection(tmp_path / "maps.mmetacool", "rep1")
        given = f"{collection}::entries/rep1"
        options = ["--output-map", collection, "-F"]
        message = (
            f"{collection} is the collection of the input map {given}: a run does not replace a collection with the "
            "map of its clusters; give --output-map another file."
        )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(given, tmp_path, *options))
    assert not stub_pipeline


@pytest.mark.parametrize("case", ["one library", "two libraries", "named as a collection"])
def test_pipeline_output_map_and_the_maps_of_its_libraries(case, stub_pipeline, tmp_path, run_cli):
    """A run from reads, alignments or pairs writes the maps of its
    libraries next to --output-map, <name>.mmetacool, with several (it
    removes them with one, as an earlier run's). One there already is
    refused before any work, without -F: the run removed it, with one
    library as with several, and the user's own collection with it. With
    -F, it is replaced, or removed. An --output-map named as a collection,
    where the maps of the libraries would be the map itself, is refused
    before any work too: the run failed once it had aligned its reads, as
    the collection was written over the map ("No such file")."""
    results = tmp_path / "results"
    results.mkdir()
    pairs = [tmp_path / "a.pairs", tmp_path / "b.pairs"][: 1 if case == "one library" else 2]
    for path in pairs:
        path.write_text("## pairs format v1.0\n")
    output_map, collection = results / "sample.metacool", results / "sample.mmetacool"
    if case == "named as a collection":
        output_map = collection
        message = f"--output-map takes the path of a .metacool file, not of a .mmetacool collection: {collection}."
    elif case == "one library":
        message = (
            f"{collection} exists next to --output-map {output_map}, where a run with several libraries writes the "
            "maps of its libraries: this run removes it, as an earlier run's. Give -F to remove it, or another "
            "--output-map."
        )
    else:
        message = (
            f"{collection} exists, where this run writes the maps of its libraries, next to --output-map {output_map}: "
            "give -F to replace it, or another --output-map."
        )
    if case != "named as a collection":
        collection.write_text("the user's collection")
    args = ["-1", ",".join(map(str, pairs)), "-a", _assembly(tmp_path), "-S", "pair", "-o", tmp_path / "out",
            "-T", tmp_path / "tmp", "--output-map", output_map]  # fmt: skip
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *args)
    assert not stub_pipeline
    if case == "named as a collection":
        return
    assert collection.read_text() == "the user's collection"
    run_cli("pipeline", *args, "-F")
    assert stub_pipeline["reassign_metacool"][3] == str(output_map)
    if case == "one library":
        assert not collection.exists()
    else:
        assert collection.read_text() == "maps of the libraries"


def test_pipeline_output_map_in_a_new_directory(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """--output-map in a directory that is not there yet: the directory is
    made, as the map will be written there (the lock failed on it)."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    output = tmp_path / "new" / "copy.metacool"
    run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", output))
    assert output.exists() and not Path(f"{output}.metator.lock").exists()


@pytest.mark.parametrize("root", ["::/", "::"])
def test_pipeline_from_the_root_of_a_map(root, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map given as file::/ (or file::), the root of the file, is the file
    itself: it is updated in place as the file, next to which the run
    stages its map. The URI, as the path of the map to write, made
    metacooler write a group inside the file, truncating the input map."""
    import metacooler

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    run_cli("pipeline", *_map_run(f"{metacool}{root}", tmp_path))
    metacooler.validate(str(metacool), full=True)
    assert set(metacooler.MetaCooler(str(metacool)).genomes.name) == {"metator_00001_00000", "__unbinned__"}
    assert sorted(path.name for path in tmp_path.iterdir() if path.name.startswith(("input", ".input"))) == ["input.metacool"]


@pytest.mark.parametrize("writer", ["write_metacool", "reassign_metacool"])
def test_a_map_is_written_to_a_file_only(writer, tmp_path):
    """A map is not written to a URI (file::group), which metacooler reads
    as a group inside the file: an existing file would be truncated."""
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    contig_data = tmp_path / "contig_data_final.txt"
    contig_data.write_text("ID\tName\tSize\tFinal_cluster\n1\tc1\t100\tND\n")
    outfile = f"{metacool}::/"
    message = f"^{re.escape(outfile)} is not the path of a file: a map is written to a file of its own, not to a group"
    with pytest.raises(mtmc.MetacoolError, match=message):
        if writer == "write_metacool":
            pairs = tmp_path / "a.pairs"
            pairs.write_text("## pairs format v1.0\n")
            mtmc.write_metacool(str(pairs), str(contig_data), outfile, force=True)
        else:
            mtmc_reassign_metacool(str(metacool), str(contig_data), outfile, force=True)
    assert metacool.read_bytes() == before


def test_metator_metacool_to_the_root_of_a_map(tmp_path, run_cli):
    """metator metacool -o file::/ is refused, the file left as it was (it
    was truncated, -F or not)."""
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    contig_data = tmp_path / "contig_data_final.txt"
    contig_data.write_text("ID\tName\tSize\tFinal_cluster\n1\tc1\t100\tND\n")
    other = _contig_map(tmp_path / "other.metacool")
    with pytest.raises(mtmc.MetacoolError, match="is not the path of a file"):
        run_cli("metacool", "-c", contig_data, "-o", f"{metacool}::/", other)
    assert metacool.read_bytes() == before and not Path(f"{metacool}.metator.lock").exists()


def test_metator_metacool_to_a_uri_makes_no_directory(tmp_path, run_cli):
    """metator metacool -o dir::/sub/x.metacool is refused before the
    directory of the map is made (an empty dir::/sub was left)."""
    other = _contig_map(tmp_path / "other.metacool")
    outfile = tmp_path / "dir::" / "sub" / "x.metacool"
    with pytest.raises(mtmc.MetacoolError, match="is not the path of a file"):
        run_cli("metacool", "-c", tmp_path / "contig_data_final.txt", "-o", outfile, other)
    assert not (tmp_path / "dir::").exists()


def test_pipeline_map_in_an_output_directory_with_colons(stub_pipeline, tmp_path, run_cli):
    """A run from pairs whose map, <outdir>/<prefix>.metacool, would be a
    URI is refused before any work (it was refused once the map of the
    contigs was built, after the alignment of a run from reads)."""
    out = tmp_path / "out::x"
    with pytest.raises(mtmc.MetacoolError):
        run_cli("pipeline", "-1", _pairs(tmp_path), "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    assert not stub_pipeline
    (log,) = out.glob("metator_pipeline_*.log")
    message = (
        f"{out / 'metator.metacool'}, the map of this run, is not the path of a file: give -o/--outdir and --prefix without ::."
    )
    assert f"ERROR :: {message}\n" in log.read_text()


@pytest.mark.parametrize("command", ["pipeline", "metacool"])
def test_map_directory_that_cannot_be_made(command, stub_pipeline, tmp_path, run_cli):
    """The directory of the map to write cannot be made (a file is in its
    way): refused before any work, the reason logged."""
    blocking = tmp_path / "afile"
    blocking.write_text("x")
    output = blocking / "x.metacool"
    metacool = _contig_map(tmp_path / "input.metacool")
    option = "--output-map" if command == "pipeline" else "-o/--outfile"
    message = f"The directory of {option} {output} cannot be made (File exists): write the map elsewhere."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        if command == "pipeline":
            run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", output))
        else:
            run_cli("metacool", "-c", tmp_path / "contig_data_final.txt", "-o", output, metacool)
    assert not stub_pipeline


def test_sigterm_stops_a_run_as_ctrl_c_does(monkeypatch):
    """A job scheduler's SIGTERM stops the run with SystemExit(143), so that
    it cleans up after itself; the handler, and sys.unraisablehook, are put
    back after it. Under a handler that a host set from C, nothing is
    changed."""
    import signal
    import sys

    hook = sys.unraisablehook
    with pytest.raises(SystemExit) as stopped:
        with mtcli.terminated_as_interrupted():
            signal.raise_signal(signal.SIGTERM)
    assert stopped.value.code == 143
    assert signal.getsignal(signal.SIGTERM) is signal.SIG_DFL and sys.unraisablehook is hook
    monkeypatch.setattr(signal, "getsignal", lambda signum: None)
    monkeypatch.setattr(signal, "signal", lambda *args: pytest.fail("the handler must be left alone"))
    with mtcli.terminated_as_interrupted():
        assert sys.unraisablehook is hook


def test_a_process_forked_by_a_run_dies_of_sigterm():
    """A process that a run forks (a worker of a multiprocessing pool, which
    the pool ends with SIGTERM) dies of SIGTERM, as without the handler of
    the run: it raised SystemExit in whatever it ran, its exit handlers
    included, which printed a traceback ("Exception ignored in atexit
    callback"). Here the child would then exit with status 0."""
    import signal
    import time
    import warnings

    with mtcli.terminated_as_interrupted():
        # SIGTERM waits until the child is in its try block.
        signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM})
        with warnings.catch_warnings():
            # The test process may run threads: the child only sleeps.
            warnings.simplefilter("ignore", DeprecationWarning)
            pid = os.fork()
        if pid == 0:
            try:
                signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGTERM})
                time.sleep(30)
            finally:
                os._exit(0)
        signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGTERM})
        os.kill(pid, signal.SIGTERM)
        _, status = os.waitpid(pid, 0)
    assert os.WIFSIGNALED(status) and os.WTERMSIG(status) == signal.SIGTERM


@pytest.mark.parametrize("thread", ["main thread", "another thread"])
def test_a_sigterm_sent_as_a_process_is_forked_kills_it(thread):
    """A process that a run forks dies of a SIGTERM sent the instant it is
    forked, before it runs its after-fork code, from the main thread or
    another (the thread of a multiprocessing pool that replaces a worker).
    The handler of the run caught such a SIGTERM, then CPython dropped it:
    about one child in four lived on. A pool closed and dropped at once
    (metator.partition.get_hamming_distance) ends its workers so, and waited
    forever for the one that missed it. The handler and the signal mask of
    the run are left as they were."""
    import signal
    import threading
    import time
    import warnings

    statuses = []

    def fork_and_kill():
        for _ in range(100):
            pid = os.fork()
            if pid == 0:
                time.sleep(0.05)
                os._exit(0)
            os.kill(pid, signal.SIGTERM)
            statuses.append(os.waitpid(pid, 0)[1])

    with mtcli.terminated_as_interrupted():
        handler = signal.getsignal(signal.SIGTERM)
        with warnings.catch_warnings():
            # The test process may run threads: the children only sleep.
            warnings.simplefilter("ignore", DeprecationWarning)
            if thread == "main thread":
                fork_and_kill()
            else:
                forking = threading.Thread(target=fork_and_kill)
                forking.start()
                forking.join()
        assert signal.getsignal(signal.SIGTERM) is handler
        assert signal.SIGTERM not in signal.pthread_sigmask(signal.SIG_BLOCK, set())
    assert len(statuses) == 100
    assert [status for status in statuses if not (os.WIFSIGNALED(status) and os.WTERMSIG(status) == signal.SIGTERM)] == []


def _run_with_a_stop_handled_in(where, signame):
    """The exit status and outputs of a Python process that runs under the
    SIGTERM handler of a run (terminated_as_interrupted) and gets a SIGTERM,
    or a Ctrl-C (SIGINT, with Python's own handler, which raises
    KeyboardInterrupt), that is handled where CPython ignores what it
    raises: in a fork hook, or in a finalizer. It prints "the run went on"
    if it does."""
    import subprocess
    import sys

    stop = {
        # A fork hook that runs before metaTOR's sends it to the main thread.
        # SIGTERM is blocked there until metaTOR's unblocks it, in its fork
        # hook; SIGINT is not: it is handled at once, in this hook, as a
        # Ctrl-C is in that of logging when it comes as a pool forks.
        "a fork hook": (
            f"os.register_at_fork(after_in_parent=lambda: signal.pthread_kill(threading.get_ident(), signal.{signame}))",
            "    pid = os.fork()\n    if pid == 0:\n        os._exit(0)\n    os.waitpid(pid, 0)",
        ),
        # The finalizer of a dropped object, as a multiprocessing pool ends
        # its workers in its finalizer.
        "a finalizer": (
            "class Dropped:\n    pass",
            f"    dropped = Dropped()\n    weakref.finalize(dropped, signal.raise_signal, signal.{signame})\n    del dropped",
        ),
    }[where]
    script = (
        "import os, signal, threading, time, weakref\n"
        # Python's handler of Ctrl-C, even when the tests run with SIGINT
        # ignored (Python then leaves it ignored).
        "signal.signal(signal.SIGINT, signal.default_int_handler)\n"
        f"{stop[0]}\n"
        "from metator.cli import terminated_as_interrupted\n"
        "with terminated_as_interrupted():\n"
        f"{stop[1]}\n"
        "    time.sleep(2)\n"
        "print('the run went on')\n"
    )
    return subprocess.run([sys.executable, "-W", "ignore", "-c", script], capture_output=True, text=True, timeout=120)


@pytest.mark.parametrize("signame", ["SIGTERM", "SIGINT"])
@pytest.mark.parametrize("where", ["a fork hook", "a finalizer"])
def test_a_stop_handled_where_python_ignores_it_kills_the_run(where, signame):
    """A SIGTERM that the handler of the run handles in code whose exceptions
    CPython ignores kills the run at once, as without the handler: the run
    went on, its SIGTERM lost. That is a SIGTERM that comes as the run forks
    from its main thread, where SIGTERM is blocked across the fork and handled
    in metaTOR's fork hook once unblocked; or one that comes in a finalizer,
    as in that of a multiprocessing pool dropped at once, which ends its
    workers there (metator.partition.get_hamming_distance). A Ctrl-C raised
    there kills the run with SIGINT so: the run went on, and replaced the
    map that it was stopped from writing, with exit status 0."""
    import signal

    result = _run_with_a_stop_handled_in(where, signame)
    assert result.returncode == -getattr(signal, signame), result.stdout + result.stderr
    assert "the run went on" not in result.stdout


def test_metator_metacool_refused_while_a_run_holds_its_map(tmp_path, run_cli):
    """metator metacool locks the map it writes, as a pipeline run does: a
    map that a run holds is not written over (the last writer won)."""
    metacool = _contig_map(tmp_path / "input.metacool")
    lock = _lock(metacool, os.getpid())
    contig_data = tmp_path / "contig_data_final.txt"
    contig_data.write_text("ID\tName\tSize\tFinal_cluster\n1\tc1\t100\tND\n")
    message = "is being written by another metaTOR run .*: wait for it to end, write this run's map elsewhere with "
    with pytest.raises(mtmc.MetacoolError, match=message + "another -o/--outfile, or"):
        run_cli("metacool", "-c", contig_data, "-o", metacool, "-F", metacool)
    assert lock.exists()


def test_pipeline_from_a_map_takes_a_stale_lock(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The lock of a run that no longer runs on this host (killed with
    SIGKILL, say) is removed with a warning, and the run goes on; its own
    lock is removed at its end."""
    import subprocess

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    process = subprocess.Popen([sys.executable, "-c", "pass"])
    process.wait()
    lock = _lock(metacool, process.pid)
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert not lock.exists()
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert (
        f"WARNING :: Removing {lock}: the run that held it (pid {process.pid}, since 2026-11-18T10:00:00+0000) is gone."
    ) in log.read_text()


def test_pipeline_from_a_map_without_writing_it_needs_no_lock(stub_pipeline, tmp_path, run_cli):
    """--no-metacool from a map: the map is read, never written, so a run
    goes on while another holds its lock, and leaves the map as it was.

    A guard: without the lock, nothing refuses this run either; it holds
    the lock to the runs that write a map.
    """
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    lock = _lock(metacool, os.getpid())
    run_cli("pipeline", *_map_run(metacool, tmp_path, "--no-metacool"))
    assert metacool.read_bytes() == before and lock.exists()
    assert "reassign_metacool" not in stub_pipeline


def test_pid_is_running_when_it_cannot_be_signalled(monkeypatch):
    """A process of another user cannot be signalled, but runs."""

    def kill(pid, signal):
        raise PermissionError

    monkeypatch.setattr("os.kill", kill)
    assert mtmc._pid_is_running(1)


def test_network_from_a_metacool(monkeypatch, tmp_path, run_cli):
    """metator network MAP reads the network from the map, with the contig
    data of the assembly."""
    calls = {}

    def create_contig_data(*args):
        calls["create_contig_data"] = args
        return {"c1": {}}, None

    def contig_data_from_map(*args):
        calls["contig_data_from_map"] = args
        return []

    def metacool_to_contacts(*args, **kwargs):
        calls["metacool_to_contacts"] = args

    def not_called(*args):
        raise AssertionError("no alignment nor pairs from a metacool")

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    monkeypatch.setattr("metator.network.contig_data_from_map", contig_data_from_map)
    monkeypatch.setattr("metator.network.metacool_to_contacts", metacool_to_contacts)
    monkeypatch.setattr("metator.network.alignment_to_contacts", not_called)
    monkeypatch.setattr("metator.align.get_contact_pairs", not_called)
    out = tmp_path / "out"
    # The URI as typed, relative to the working directory.
    _collection(tmp_path / "maps.mmetacool", "lib1")
    monkeypatch.chdir(tmp_path)
    run_cli("network", "maps.mmetacool::entries/lib1", "-a", _assembly(tmp_path), "-e", "DpnII",
            "-s", "-t", "2", "-n", "RS", "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert calls["create_contig_data"] == (str(tmp_path / "assembly.fa"), 1, None, "DpnII")
    # The restriction sites and the coverage that the map gives where the
    # options give none.
    assert calls["contig_data_from_map"] == ("maps.mmetacool::entries/lib1", {"c1": {}}, "DpnII", None)
    uri, contig_data, *args = calls["metacool_to_contacts"]
    assert uri == "maps.mmetacool::entries/lib1" and contig_data == {"c1": {}}
    out_dir, network, contigs, tmp_dir, threads, normalization, self_contacts = args
    assert (out_dir, network, contigs, threads, normalization, self_contacts) == (
        str(out), "network.txt", "contig_data_network.txt", 2, "RS", True)  # fmt: skip
    assert Path(tmp_dir).parent == tmp_path / "tmp"


# What a network or pipeline run cannot take, refused before any work, the
# reason in its log file: a map with the options of a run from reads,
# alignments or pairs, a first argument that is no map, and no input.
RUN_FROM_A_MAP = "A run from the contact map {map} takes no {options}: these give the reads, alignments or pairs of a run that does not start from a map."
RUN_INPUT_ERRORS = {
    "map and -1": (["{map}", "-1", "a.pairs"], RUN_FROM_A_MAP.replace("{options}", "-1/--forward")),
    "map and -2": (["{map}", "-2", "b.fq"], RUN_FROM_A_MAP.replace("{options}", "-2/--reverse")),
    "map and --start": (["{map}", "-S", "pair"], RUN_FROM_A_MAP.replace("{options}", "--start")),
    "map and --start fastq": (["{map}", "--start=fastq"], RUN_FROM_A_MAP.replace("{options}", "--start")),
    "map, -1 and --start": (
        ["{map}", "-1", "a.pairs", "-S", "pair"],
        RUN_FROM_A_MAP.replace("{options}", "-1/--forward or --start"),
    ),
    "reads as the first argument": (
        ["{reads}"],
        "{reads} is not a .metacool contact map: give reads, alignments or pairs with -1/--forward (and --start).",
    ),
    "reads as the first argument, with -2": (
        ["{reads}", "-2", "b.fq"],
        "{reads} is not a .metacool contact map: give reads, alignments or pairs with -1/--forward (and --start).",
    ),
    "a map that is not there": (["{missing}"], "Cannot read the contact map {missing}: no such file."),
    "a URI of a directory": (["{directory}::/"], "Cannot read the contact map {directory}::/: {directory} is not a file."),
    "a pairs file named as a URI": (
        ["{reads}::/"],
        "{reads}::/ is not a .metacool contact map: give reads, alignments or pairs with -1/--forward (and --start).",
    ),
    "a pairs file whose name holds ::": (
        ["{colons}"],
        "{colons} is not a .metacool contact map: give reads, alignments or pairs with -1/--forward (and --start).",
    ),
    "a map whose name holds ::": (
        ["{colon_map}"],
        "Cannot read the contact map {colon_map}: metacooler reads a path that holds :: as file::group. Rename the map "
        "without ::.",
    ),
    "no input": (
        [],
        "Give the run a .metacool contact map, as the first argument, or its reads, alignments or pairs, with "
        "-1/--forward (and --start).",
    ),
}


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize("case", list(RUN_INPUT_ERRORS))
def test_run_input_refused(command, case, stub_pipeline, tmp_path, run_cli):
    """A run takes a .metacool map, as its first argument, or the reads,
    alignments or pairs of -1/--forward (and --start), not both nor neither:
    refused before any work, the reason in the log file. --start given with
    a map is refused even at its default value, fastq. A first argument
    that is no map, not there, or a URI whose file is no map nor a file, is
    refused too (it failed once the run read the map). A file whose name
    holds :: is that file: a pairs file is no map (it was refused as
    missing), and a map is refused with the reason, as metacooler cannot
    read it."""
    names = {
        "map": _contig_map(tmp_path / "input.metacool"),
        "reads": tmp_path / "R1.fq",
        "missing": tmp_path / "missing.metacool",
        "directory": tmp_path / "maps",
        "colons": tmp_path / "lib::1.pairs",
        "colon_map": tmp_path / "run::1.metacool",
    }
    names["reads"].write_text("@r1\nACGT\n+\nIIII\n")
    names["directory"].mkdir()
    names["colons"].write_text("## pairs format v1.0\n")
    shutil.copyfile(names["map"], names["colon_map"])
    args, message = RUN_INPUT_ERRORS[case]
    args = [arg.format(**names) for arg in args]
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli(command, *args, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert not stub_pipeline
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message.format(**names)}\n" in log.read_text()


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_a_collection_is_no_map(command, stub_pipeline, tmp_path, run_cli):
    """A .mmetacool collection given as the map, rather than one of its
    entries, is refused before any work, with metacooler's reason, which
    names its entries (it failed once the run read the map)."""
    collection = _collection(tmp_path / "maps.mmetacool", "lib1")
    out = tmp_path / "out"
    with pytest.raises(mtmc.MetacoolError, match=f"^Cannot read the contact map {re.escape(str(collection))}: .*entries/lib1"):
        run_cli(command, collection, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert not stub_pipeline
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: Cannot read the contact map {collection}: " in log.read_text()


# A first argument that is no map, given with -1/--forward: one argument too
# many, as click refused it before the map was an argument.
EXTRA_ARGUMENTS = {
    "the reverse reads without -2": (["-1", "R1.fq", "R2.fq"], "R2.fq"),
    "a value split by a missing space": (["-1", "a.pairs", "-S", "pair", "pairs"], "pairs"),
    "a map that is not there": (["-1", "a.pairs", "-S", "pair", "missing.metacool"], "missing.metacool"),
}


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize("case", list(EXTRA_ARGUMENTS))
def test_an_extra_argument_is_a_usage_error(command, case, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A first argument that is no map, given with -1/--forward, is one
    argument too many: a usage error, as click said before the map was an
    argument, so this passes against the previous commit too. Without
    Command's check, run_input refused it as a map. The command does not
    start."""
    monkeypatch.chdir(tmp_path)
    Path("R2.fq").write_text("@r1\nACGT\n+\nIIII\n")
    args, extra = EXTRA_ARGUMENTS[case]
    result = run_cli(command, *args, "-a", "a.fa", "-o", "out", exit_code=2)
    assert f"Got unexpected extra argument ({extra})" in result.output
    assert not stub_pipeline and not Path("out").exists()


@pytest.mark.parametrize("command", ["network", "pipeline", "partition", "validation"])
def test_a_directory_is_no_map(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A directory given as the map is a usage error, before the command
    starts (it failed once the run read the map, and pipeline -F had
    removed the clusters of the earlier run)."""
    monkeypatch.chdir(tmp_path)
    Path("maps").mkdir()
    bins = Path("out", "overlapping_clusters", "bin.fa")
    bins.parent.mkdir(parents=True)
    bins.write_text(">c1\nACGT\n")
    force = [] if command == "network" else ["-F"]
    result = run_cli(command, "maps", "-a", _assembly(tmp_path), "-o", "out", *force, exit_code=2)
    assert "File 'maps' is a directory." in " ".join(result.output.replace("│", " ").split())
    assert not stub_pipeline and bins.exists() and not list(Path("out").glob("*.log"))


def test_pipeline_from_a_map_cannot_scaffold(stub_pipeline, tmp_path, run_cli):
    """--scaffold reads the pairs files, which a run from a map has not:
    refused before any work."""
    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli("pipeline", metacool, "--scaffold", "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert not stub_pipeline
    (log,) = out.glob("metator_pipeline_*.log")
    assert "ERROR :: --scaffold reads the pairs files: it cannot run from a contact map.\n" in log.read_text()


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize("start", [["-S", "metacool"], ["--start=metacool"]], ids=["-S", "--start="])
def test_start_metacool_is_gone(command, start, stub_pipeline, tmp_path, run_cli):
    """--start metacool, how a run started from a map before metaTOR 2.0, is
    a usage error that says what replaced it: the map as the first argument.
    The command does not start: no log file."""
    result = run_cli(command, *start, "-1", "map.metacool", "-a", "a.fa", "-o", tmp_path / "out", exit_code=2)
    output = " ".join(result.output.replace("│", " ").split())
    message = (
        f"--start metacool is gone in metaTOR 2.0: give the map as the first argument, metator {command} MAP -a ASSEMBLY ..."
    )
    assert message in output
    assert not stub_pipeline and not (tmp_path / "out").exists()


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_help_says_from_a_map_not_start_metacool(command):
    """The help of the options that the columns of a map stand in for
    (--depth, --enzyme, the normalization) says what they do from a map
    (MAP): it still described --start metacool, a usage error since 2.0."""
    from metator.cli import cli

    helps = {param.name: param.help for param in cli.commands[command].params if getattr(param, "help", None)}
    assert [name for name, text in helps.items() if "--start metacool" in text] == []
    for name in ["depth", "enzyme", "normalization"]:
        assert "From a map (MAP), its " in helps[name], name


EDGE_GONE = (
    "--edge is gone in metaTOR 2.0: the network counts every pair, as its contact map does, which holds no read positions."
)


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize("edge", [["-E", "100"], ["--edge", "0"], ["--edge=-5"]], ids=["-E", "--edge 0", "--edge="])
def test_edge_is_gone(command, edge, stub_pipeline, tmp_path, run_cli):
    """--edge, which metaTOR 2.0 removed (the network counts every pair, as
    its map does), is a usage error that says so, at any value: not an
    unknown option. The command does not start."""
    result = run_cli(command, "-1", "a.pairs", "-S", "pair", *edge, "-a", "a.fa", "-o", tmp_path / "out", exit_code=2)
    assert EDGE_GONE in " ".join(result.output.replace("│", " ").split())
    assert not stub_pipeline and not (tmp_path / "out").exists()


# partition MAP and validation MAP: the steps run one by one on a map, with
# the steps replaced by stubs (see tests/test_binning_synthetic.py for the
# map that they give).

# The directories of clusters of each step, replaced with -F.
STEP_DIRS = {"partition": ["overlapping_clusters"], "validation": ["recursive_clusters", "final_clusters"]}


def _step_run(command, metacool, tmp_path, *options):
    """The arguments of a partition or validation run on the map metacool
    (none for None), from tables of the earlier steps that the stub steps
    do not read, but for the contig data, which the run checks the map
    against: the contigs c1 and c2 of _contig_map."""
    contigs = tmp_path / "contig_data.txt"
    if not contigs.exists():
        contigs.write_text("ID\tName\tSize\n1\tc1\t100\n2\tc2\t50\n")
    tables = ["-c", contigs, "-n", "network.txt", *(["-f", "overlapping_clusters"] if command == "validation" else [])]
    out = ["-o", tmp_path / "out", "-T", tmp_path / "tmp"]
    return [*([] if metacool is None else [metacool]), "-a", _assembly(tmp_path), *tables, *out, *options]


def _earlier_run(command, tmp_path):
    """The directories of clusters that an earlier run of command left in the
    output directory, each with a FASTA file, which the test checks."""
    files = [tmp_path / "out" / name / "earlier.fa" for name in STEP_DIRS[command]]
    for path in files:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    return files


# The outputs that an earlier run of each step left, named as metaTOR 1.x
# names them, with their names since 2.0.
STEP_OUTPUTS_OF_1X = {
    "partition": {"overlapping_bin": "overlapping_clusters"},
    "validation": {
        "recursive_bin": "recursive_clusters",
        "final_bin": "final_clusters",
        "final_bin_unscaffold": "final_clusters",
        "bin_summary.txt": "cluster_summary.txt",
    },
}


def _outputs_of_1x(out, names):
    """The outputs that an earlier run left in out, named as metaTOR 1.x
    names them: a directory with a FASTA file, or bin_summary.txt."""
    paths = []
    for name in names:
        path = out / name
        if name.endswith(".txt"):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("\tCompleteness\nmetator_00001_00000\t0.9\n")
        else:
            path.mkdir(parents=True)
            (path / "earlier.fa").write_text(">c1\nACGT\n")
        paths.append(path)
    return paths


def _partition_table(monkeypatch, size=150):
    """Have the stub partition write the contig data of a partition of the
    contigs c1 and c2 of _contig_map, as metator.partition writes it: one
    overlapping cluster of both, of size bp (00001: metator_00001_00000)."""
    stub = mtcli.mtp.partition

    def partition(*args):
        result = stub(*args)
        Path(args[6], "contig_data_partition.txt").write_text(
            "ID\tName\tSize\tOverlapping_cluster_ID\tOverlapping_cluster_contigs\tOverlapping_cluster_size\n"
            f"1\tc1\t100\t00001\t2\t{size}\n2\tc2\t50\t00001\t2\t{size}\n"
        )
        return result

    monkeypatch.setattr("metator.partition.partition", partition)


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_help_of_a_step_on_a_map(command, run_cli):
    """The help of partition and validation shows MAP, their optional first
    argument, and --output-map."""
    output = " ".join(run_cli(command, "--help").output.split())
    assert f"Usage: metator {command} [OPTIONS] [MAP]" in output
    assert f"metator {command} MAP" in output and "--output-map" in output


# What partition MAP and validation MAP refuse before any work, with -F (an
# earlier run's directories of clusters are kept), the reason in the log file
# of the run: the map that run_input refuses, and what run_map_path refuses.
STEP_MAP_ERRORS = {
    "a map that is not there": ("{missing}", [], "Cannot read the contact map {missing}: no such file."),
    "a URI of a directory": ("{directory}::/", [], "Cannot read the contact map {directory}::/: {directory} is not a file."),
    "no map": (
        "{pairs}",
        [],
        "{pairs} is not a .metacool contact map: give the map of the contigs of the run as the first argument of "
        "{command}, or leave it out.",
    ),
    "a collection": ("{collection}", [], "Cannot read the contact map {collection}: "),
    "a map of other contigs": (
        "{other_map}",
        [],
        "{other_map} is not a contact map of the contigs of {contigs}: 1 of its contig(s) are not in it: 'c3'.",
    ),
    "a map of other lengths": (
        "{long_map}",
        [],
        "{long_map} is not a contact map of the contigs of {contigs}: 1 contig(s) have another length in it: 'c1' (120 "
        "bp in the map, 100 bp in {contigs}).",
    ),
    "an input as --output-map": (
        "{map}",
        ["--output-map", "{contigs}"],
        "--output-map {contigs} is {contigs}, an input of the run: write the map elsewhere.",
    ),
    "a map whose name holds ::": (
        "{colon_map}",
        [],
        "Cannot read the contact map {colon_map}: metacooler reads a path that holds :: as file::group. Rename the map "
        "without ::.",
    ),
    "an entry without --output-map": (
        "{entry}",
        [],
        "{entry} is an entry of a .mmetacool collection, which a run does not update in place: write the map of its "
        "clusters to another file with --output-map.",
    ),
    "--output-map without a map": (
        None,
        ["--output-map", "{copy}"],
        "--output-map writes the map that {command} MAP updates: give the map as the first argument, or drop --output-map.",
    ),
    "a directory as --output-map": (
        "{map}",
        ["--output-map", "{directory}"],
        "{directory} is a directory: --output-map takes the path of a .metacool file.",
    ),
    "a lock held by another run": (
        "{map}",
        [],
        "{map} is being written by another metaTOR run (pid {pid} on {host}, since 2026-11-18T10:00:00+0000): wait "
        "for it to end, write this run's map elsewhere with --output-map, or, if no such run is going on, remove "
        "{map}.metator.lock.",
    ),
}


@pytest.mark.parametrize("command", ["partition", "validation"])
@pytest.mark.parametrize("case", list(STEP_MAP_ERRORS))
def test_step_on_a_map_refused(command, case, stub_pipeline, tmp_path, run_cli):
    """partition MAP and validation MAP check the map as network MAP and
    pipeline MAP do, then lock the map that they write, before any work:
    refused, the reason logged as an error in the log file of the run and
    raised as a MetacoolError, before the step and before -F replaces the
    directories of clusters of an earlier run. A lock that a running run
    holds (here, this process) is left to it. A collection is refused with
    metacooler's reason, which names its entries."""
    import socket

    names = {
        "command": command,
        "map": _contig_map(tmp_path / "input.metacool"),
        "missing": tmp_path / "missing.metacool",
        "directory": tmp_path / "maps",
        "pairs": _pairs(tmp_path),
        "collection": _collection(tmp_path / "maps.mmetacool", "lib1"),
        "copy": tmp_path / "copy.metacool",
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "colon_map": tmp_path / "run::1.metacool",
        "other_map": _contig_map(tmp_path / "other.metacool", {"c1": 100, "c3": 50}),
        "long_map": _contig_map(tmp_path / "long.metacool", {"c1": 120, "c2": 50}),
        "contigs": tmp_path / "contig_data.txt",
    }
    names["contigs"].write_text("ID\tName\tSize\n1\tc1\t100\n2\tc2\t50\n")
    shutil.copyfile(names["map"], names["colon_map"])
    names["entry"] = f"{names['collection']}::entries/lib1"
    names["directory"].mkdir()
    metacool, options, message = STEP_MAP_ERRORS[case]
    if case == "a lock held by another run":
        lock = _lock(names["map"], os.getpid())
    earlier = _earlier_run(command, tmp_path)
    metacool = None if metacool is None else metacool.format(**names)
    message = message.format(**names)
    pattern = f"^{re.escape(message)}" + (".*entries/lib1" if case == "a collection" else "$")
    with pytest.raises(mtmc.MetacoolError, match=pattern):
        run_cli(command, *_step_run(command, metacool, tmp_path, *[option.format(**names) for option in options], "-F"))
    assert not stub_pipeline and all(path.exists() for path in earlier)
    (log,) = (tmp_path / "out").glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}" in log.read_text()
    if case == "a lock held by another run":
        assert lock.exists()
    assert not (tmp_path / "copy.metacool").exists()


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_step_on_a_map_refuses_an_existing_output_map(command, stub_pipeline, tmp_path, run_cli):
    """--output-map replaces a file only with -F, as in pipeline: refused
    before any work, the file left as it was."""
    metacool = _contig_map(tmp_path / "input.metacool")
    output_map = tmp_path / "copy.metacool"
    output_map.write_text("an older map")
    message = f"{output_map} exists: give -F to replace it with the map of this run, or another --output-map."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *_step_run(command, metacool, tmp_path, "--output-map", output_map))
    assert not stub_pipeline and output_map.read_text() == "an older map"
    (log,) = (tmp_path / "out").glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("output", ["in place", "--output-map", "-F --output-map"])
def test_partition_on_a_map_stores_its_overlapping_clusters(output, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """partition MAP stores the overlapping clusters that reach --size in the
    clustering rounds of the map, the round overlapping, with their contigs
    and without the columns of miComplete, in place of the rounds the map
    held; its genomes are left as they were. In place, the map keeps its
    mode; with --output-map, the map is left as it was, and the copy is
    written to the file (an existing one, with -F, is replaced: removed
    first, as pipeline does). The map is locked while the run writes it."""
    import metacooler
    import pandas as pd

    _partition_table(monkeypatch)
    metacool = _contig_map(tmp_path / "input.metacool")
    metacooler.write_clusters(
        str(metacool),
        ["earlier"],
        pd.DataFrame({"round": ["earlier"], "name": ["x"]}),
        pd.DataFrame({"round": ["earlier"], "cluster": ["x"], "contig": ["c2"]}),
    )
    metacool.chmod(0o640)
    before = metacool.read_bytes()
    output_map = tmp_path / "copy.metacool"
    options = [] if output == "in place" else ["--output-map", output_map]
    if output.startswith("-F"):
        output_map.write_text("an older map")
        options.append("-F")
    locks = []
    write_rounds = mtmc.write_rounds

    def locked(*args, **kwargs):
        locks.append(sorted(path.name for path in tmp_path.glob("*.lock")))
        return write_rounds(*args, **kwargs)

    monkeypatch.setattr("metator.metacool.write_rounds", locked)
    run_cli("partition", *_step_run("partition", metacool, tmp_path, "-s", "100", *options))
    written = metacool if output == "in place" else output_map
    assert locks == [[f"{written.name}.metator.lock"]] and not list(tmp_path.glob("*.lock"))
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    assert mc.cluster_rounds == [{"name": "overlapping", "tool": "metaTOR", "version": __version__, "step": "partition"}]
    assert list(mc.clusters.name) == ["metator_00001_00000"] and list(mc.clusters.n_contigs) == [2]
    # genome, which metacooler computes: no cluster is a genome of the map yet.
    assert list(mc.clusters.columns) == ["round", "name", "parent_round", "parent", "n_contigs", "length", "genome"]
    assert mc.clusters.genome.isna().all()
    assert sorted(mc.cluster_members().contig) == ["c1", "c2"]
    assert list(mc.genomes.name) == ["__unbinned__"]
    assert (metacool.read_bytes() == before) == (output != "in place")
    if output == "in place":
        assert metacool.stat().st_mode & 0o777 == 0o640
    assert not [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]
    (log,) = (tmp_path / "out").glob("metator_partition_*.log")
    logged = log.read_text()
    assert f"INFO :: {written} is written: the contact map, with the clustering rounds of the run." in logged
    assert ("Removing" in logged) == output.startswith("-F")
    # What replaces the map, in place: not the map of the clusters.
    replaced = f"INFO :: {metacool} is the input map {metacool}: it is replaced once the map with its overlapping clusters"
    assert (replaced in logged) == (output == "in place")


def test_partition_on_a_map_without_a_cluster_of_its_size(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """No overlapping cluster reaches --size (the default, 500 kb, here): the
    map holds the round overlapping, with no cluster, in place of the rounds
    it held."""
    import metacooler

    _partition_table(monkeypatch)
    metacool = _contig_map(tmp_path / "input.metacool")
    run_cli("partition", *_step_run("partition", metacool, tmp_path))
    mc = metacooler.MetaCooler(str(metacool))
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping"]
    assert mc.clusters.empty and mc.cluster_members().empty


def test_partition_on_an_entry_of_a_collection(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """From an entry of a .mmetacool collection, --output-map is a map of its
    own, with the overlapping clusters, and the collection is left as it
    was."""
    import metacooler

    _partition_table(monkeypatch)
    collection = _collection(tmp_path / "maps.mmetacool", "lib1", "lib2")
    before = collection.read_bytes()
    output_map = tmp_path / "lib1.metacool"
    run_cli(
        "partition", *_step_run("partition", f"{collection}::entries/lib1", tmp_path, "-s", "100", "--output-map", output_map)
    )
    metacooler.validate(str(output_map), full=True)
    mc = metacooler.MetaCooler(str(output_map))
    assert list(mc.clusters.name) == ["metator_00001_00000"] and list(mc.genomes.name) == ["__unbinned__"]
    assert collection.read_bytes() == before


def test_partition_on_a_map_through_a_symbolic_link(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map given through a symbolic link: the link is replaced by the
    map with the overlapping clusters, and the map it pointed to is left as
    it was, as pipeline replaces it."""
    import metacooler

    _partition_table(monkeypatch)
    target = _contig_map(tmp_path / "target.metacool")
    before = target.read_bytes()
    link = tmp_path / "link.metacool"
    link.symlink_to(target)
    run_cli("partition", *_step_run("partition", link, tmp_path, "-s", "100"))
    assert not link.is_symlink() and target.read_bytes() == before
    assert list(metacooler.MetaCooler(str(link)).clusters.name) == ["metator_00001_00000"]


@pytest.mark.parametrize("error", [ValueError("not valid"), KeyboardInterrupt()], ids=["ValueError", "Ctrl-C"])
def test_partition_on_a_map_that_cannot_be_written(error, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map with the overlapping clusters that cannot be written (here, its
    validation fails) stops the run with an error that says why and that the
    map is left as it was, as pipeline says it; Ctrl-C stops it too. The map
    is left byte for byte, with no file beside it (no staged map, no lock);
    the other outputs of the run are kept."""
    _partition_table(monkeypatch)
    (tmp_path / "maps").mkdir()
    metacool = _contig_map(tmp_path / "maps" / "input.metacool")
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise error

    monkeypatch.setattr("metacooler.validate", validate)
    args = _step_run("partition", metacool, tmp_path, "-s", "100")
    message = (
        f"{metacool} is not written: not valid. It is still the input map of the run, as it was. The other outputs of "
        "the run are kept."
    )
    if isinstance(error, ValueError):
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("partition", *args)
        (log,) = (tmp_path / "out").glob("metator_partition_*.log")
        assert f"ERROR :: {message}\n" in log.read_text()
    else:
        assert "Aborted" in run_cli("partition", *args, exit_code=1).output
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in (tmp_path / "maps").iterdir()) == ["input.metacool"]
    assert (tmp_path / "out" / "contig_data_partition.txt").exists()


@pytest.mark.parametrize("output", ["in place", "--output-map"])
def test_validation_on_a_map_writes_the_map_of_the_clusters(output, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """validation MAP writes the map of the clusters as pipeline MAP writes it
    at its end, from the tables of the run, with its clustering rounds: in
    place, keeping the mode of the map, or to --output-map, the map left as
    it was. Its metadata are those of the validation: the map, and no
    --start, which validation does not have."""
    import metacooler

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    _rounds_tables(monkeypatch, ["c1", "c2"])
    metacool = _contig_map(tmp_path / "input.metacool")
    metacool.chmod(0o640)
    before = metacool.read_bytes()
    output_map = tmp_path / "copy.metacool"
    options = [] if output == "in place" else ["--output-map", output_map]
    run_cli("validation", *_step_run("validation", metacool, tmp_path, *options))
    written = metacool if output == "in place" else output_map
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping", "final"]
    metator = mc.metadata["metator"]
    assert metator["command"].startswith(f"metator validation {metacool} -a")
    assert metator["parameters"]["map"] == str(metacool) and "start" not in metator["parameters"]
    assert metator["parameters"]["output_map"] == (None if output == "in place" else str(output_map))
    assert (metacool.read_bytes() == before) == (output != "in place")
    # The heatmap of the clusters is drawn from the map, as pipeline MAP
    # draws it.
    assert stub_pipeline["heatmap_map"] == str(metacool)
    if output == "in place":
        assert metacool.stat().st_mode & 0o777 == 0o640
    assert not list(tmp_path.glob("*.lock")) and not [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]


def test_validation_on_a_map_that_cannot_be_written(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the clusters that cannot be written (here, its validation
    fails) stops validation MAP as it stops pipeline MAP: the map is left
    as it was."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise ValueError("not valid")

    monkeypatch.setattr("metacooler.validate", validate)
    with pytest.raises(mtmc.MetacoolError, match=r"is not written: .* It is still the input map of the run, as it was\."):
        run_cli("validation", *_step_run("validation", metacool, tmp_path))
    assert metacool.read_bytes() == before


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_a_step_on_a_map_runs_under_the_sigterm_handler(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """partition MAP and validation MAP, which write a map, run their step
    under the SIGTERM handler of a run (see terminated_as_interrupted), so
    that SIGTERM stops them as Ctrl-C does and they clean up after
    themselves, their lock and staged map. Without it, SIGTERM killed them
    and left both. It passes against this commit's first version, which
    had the handler and no test of it; it fails without the handler."""
    import signal

    seen = {}

    def step(*args, **kwargs):
        seen["handler"] = signal.getsignal(signal.SIGTERM)
        raise RuntimeError("stopped")

    target = {"partition": "metator.partition.partition", "validation": "metator.validation.recursive_decontamination"}
    monkeypatch.setattr(target[command], step)
    metacool = _contig_map(tmp_path / "input.metacool")
    with pytest.raises(RuntimeError, match="stopped"):
        run_cli(command, *_step_run(command, metacool, tmp_path))
    assert callable(seen["handler"]) and seen["handler"] not in (signal.SIG_DFL, signal.SIG_IGN)
    assert signal.getsignal(signal.SIGTERM) is signal.SIG_DFL and not list(tmp_path.glob("*.lock"))


@pytest.mark.parametrize("with_map", [False, True], ids=["without a map", "with a map"])
@pytest.mark.parametrize("command", ["partition", "validation"])
def test_a_step_refuses_a_bad_prefix_before_replacing_the_clusters(command, with_map, stub_pipeline, tmp_path, run_cli):
    """An invalid --prefix is refused before -F replaces the directories of
    clusters of an earlier run, as pipeline refuses it: it emptied them
    first, then was refused."""
    earlier = _earlier_run(command, tmp_path) + _outputs_of_1x(tmp_path / "out", STEP_OUTPUTS_OF_1X[command])
    metacool = _contig_map(tmp_path / "input.metacool") if with_map else None
    with pytest.raises(ValueError, match="Invalid --prefix"):
        run_cli(command, *_step_run(command, metacool, tmp_path, "-F", "-p", "bad.prefix"))
    assert all(path.exists() for path in earlier) and not stub_pipeline


def test_pipeline_refuses_a_file_of_the_index_as_its_output_map(stub_pipeline, tmp_path, run_cli):
    """pipeline --output-map that names a file of the index of the assembly,
    given as its basename (-a), is refused before any work: with -F, the map
    replaced it at the end of the run."""
    metacool = _contig_map(tmp_path / "input.metacool")
    index = tmp_path / "idx.1.bt2"
    index.write_text("index")
    message = f"--output-map {index} is {index}, an input of the run: write the map elsewhere."
    args = [metacool, "-a", tmp_path / "idx", "-o", tmp_path / "out", "-T", tmp_path / "tmp", "--output-map", index, "-F"]
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *args)
    assert index.read_text() == "index" and not stub_pipeline


def test_pipeline_replaces_an_earlier_map_named_after_its_assembly(stub_pipeline, tmp_path, run_cli):
    """pipeline --output-map -F replaces the map of an earlier run named
    after the assembly (assembly.fa.metacool, for -a assembly.fa): a map is
    no file of the index of the assembly, which --output-map must not name.
    It was refused as one: every file named assembly.fa.<anything> was."""
    metacool = _contig_map(tmp_path / "input.metacool")
    output_map = _contig_map(tmp_path / "assembly.fa.metacool")
    run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", output_map, "-F"))
    assert stub_pipeline["reassign_metacool"][3] == str(output_map)


def test_pipeline_refuses_an_input_as_its_output_map(stub_pipeline, tmp_path, run_cli):
    """pipeline --output-map that names the assembly of the run is refused
    before any work: with -F, the map of the clusters replaced it at the
    end of the run."""
    metacool = _contig_map(tmp_path / "input.metacool")
    assembly = Path(_assembly(tmp_path))
    before = assembly.read_bytes()
    message = f"--output-map {assembly} is {assembly}, an input of the run: write the map elsewhere."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", assembly, "-F"))
    assert assembly.read_bytes() == before and not stub_pipeline


@pytest.mark.parametrize("command, name", [(command, name) for command, names in STEP_OUTPUTS_OF_1X.items() for name in names])
def test_a_step_refuses_the_outputs_of_1x(command, name, stub_pipeline, tmp_path, run_cli, logged):
    """An output of an earlier run that metaTOR 1.x named after bins
    (overlapping_bin, bin_summary.txt...) is refused without -F, as those of
    2.0 are: the run would leave it next to its own. -F removes it."""
    (path,) = _outputs_of_1x(tmp_path / "out", [name])
    with pytest.raises(ValueError):
        run_cli(command, *_step_run(command, None, tmp_path))
    assert path.exists() and not stub_pipeline
    # Refused before the run makes its own directories.
    assert not [name for name in STEP_DIRS[command] if (tmp_path / "out" / name).exists()]
    new = STEP_OUTPUTS_OF_1X[command][name]
    message = (
        f"{path} already existed: an output of an earlier run, named as metaTOR 1.x names it ({new} now). Remove it "
        "or use -F argument to remove it."
    )
    assert ("ERROR", message) in logged()
    run_cli(command, *_step_run(command, None, tmp_path, "-F"))
    assert not path.exists() and stub_pipeline
    assert ("INFO", f"Removed {path}, an output of an earlier run named as metaTOR 1.x names it ({new} now).") in logged()


@pytest.mark.parametrize("dangling", [False, True], ids=["link", "dangling link"])
def test_force_removes_a_link_of_1x_not_its_target(dangling, stub_pipeline, tmp_path, run_cli):
    """An output of 1.x that is a symbolic link, to a directory elsewhere or
    to nothing, is refused without -F, and removed with it, its target
    kept."""
    target = tmp_path / "elsewhere"
    if not dangling:
        target.mkdir()
        (target / "earlier.fa").write_text(">c1\nACGT\n")
    link = tmp_path / "out" / "overlapping_bin"
    link.parent.mkdir()
    link.symlink_to(target)
    with pytest.raises(ValueError):
        run_cli("partition", *_step_run("partition", None, tmp_path))
    assert os.path.lexists(link) and not stub_pipeline
    run_cli("partition", *_step_run("partition", None, tmp_path, "-F"))
    assert not os.path.lexists(link) and stub_pipeline
    assert dangling or (target / "earlier.fa").exists()


@pytest.mark.parametrize("name", ["final_bin", "final_bin_unscaffold", "final_clusters", "recursive_bin/step_1"])
def test_validation_refuses_its_input_in_what_force_removes(name, stub_pipeline, tmp_path, run_cli, logged):
    """The clusters to validate (-f) in a directory that -F removes, of 1.x
    or of 2.0: refused before anything is removed. -F removed them, then
    the run failed, or found no cluster."""
    fasta = tmp_path / "out" / name
    fasta.mkdir(parents=True)
    (fasta / "metator_00001_00000.fa").write_text(">c1\nACGT\n")
    args = _step_run("validation", None, tmp_path, "-F")
    args[args.index("overlapping_clusters")] = fasta
    with pytest.raises(ValueError):
        run_cli("validation", *args)
    removed = tmp_path / "out" / name.split("/")[0]
    message = f"{fasta}, given as -f, is in {removed}, which -F removes: write the run to another directory (-o)."
    assert ("ERROR", message) in logged()
    assert (fasta / "metator_00001_00000.fa").exists() and not stub_pipeline


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_steps_without_a_map_are_unchanged(command, stub_pipeline, tmp_path, run_cli):
    """A guard: without MAP, partition and validation run as before MAP was
    their argument, so this passes against the previous commit too. The
    step gets the same arguments, the run writes the same files (no map, no
    lock, no rounds) and the same log file, line for line but for its
    dates: no line about a map. -F replaces the directories of clusters of
    an earlier run before the log file starts, as before; without -F, they
    are refused then, with no log file."""
    earlier = _earlier_run(command, tmp_path)
    with pytest.raises(ValueError):
        run_cli(command, *_step_run(command, None, tmp_path))
    assert all(path.exists() for path in earlier) and not list((tmp_path / "out").glob("*.log"))

    args = _step_run(command, None, tmp_path, "-F", "-s", "1000", "-p", "sample")
    run_cli(command, *args)
    assert not any(path.exists() for path in earlier)
    step = {"partition": "partition", "validation": "recursive_decontamination"}[command]
    out, assembly, contigs = str(tmp_path / "out"), str(tmp_path / "assembly.fa"), str(tmp_path / "contig_data.txt")
    expected = {
        "partition": (
            "louvain", assembly, False, contigs, 100, "network.txt", out, f"{out}/overlapping_clusters", 0.8, 1.0,
            1000, 1, "sample_metator",
        ),
        "validation": (
            "louvain", assembly, False, contigs, f"{out}/final_clusters", "overlapping_clusters", 10, "network.txt", out,
            0.9, f"{out}/recursive_clusters", 1.0, 1000, "tmp_dir", 1, "sample_metator",
        ),
    }[command]  # fmt: skip
    called = list(stub_pipeline[step])
    if command == "validation":
        # The private temporary directory of the run, inside --tmpdir.
        assert Path(called[13]).parent == tmp_path / "tmp"
        called[13] = "tmp_dir"
    assert tuple(called) == expected
    # Without a map, the heatmap of validation is drawn from the pairs files
    # of its output directory (no map), as before.
    assert set(stub_pipeline) == ({step, "partition_outdir"} if command == "partition" else {step, "heatmap_map"})
    assert stub_pipeline.get("heatmap_map") is None
    files = {
        "partition": ["contig_data_partition.txt", "overlapping_clusters"],
        "validation": ["final_clusters", "recursive_clusters"],
    }
    (log,) = (tmp_path / "out").glob(f"metator_{command}_*.log")
    assert sorted(path.name for path in (tmp_path / "out").iterdir()) == sorted([log.name, *files[command]])
    lines = [re.sub(r"\d{4}-\d\d-\d\d[ ,]\d\d:\d\d:\d\d", "DATE", line) for line in log.read_text().splitlines()]
    command_line = shlex.join(["metator", command, *(str(arg) for arg in args)])
    assert lines == [
        f"## MetaTOR: v{__version__} log file",
        "## date: DATE",
        f"## Command used: {command_line}",
        "",
        "---",
        "",
        "",
        "---",
        "## date: DATE",
    ]
    assert not Path(f"{assembly}.fxi").exists()


# The commands of the tutorial of metaTOR 1.x (4253551,
# docs/example/metator_tutorial.md) that its command line accepted: its -S
# pairs and -S network were no values of --start.
COMMANDS_OF_4253551 = [
    "pipeline -1 hic_reads_for.fastq.gz -2 hic_reads_rev.fastq.gz -a assembly.fa -o metator_folder --threads 16",
    "pipeline -1 hic_reads_for.fastq.gz -2 hic_reads_rev.fastq.gz -a assembly_bowtie2_index -S fastq",
    "pipeline -1 for.bam -2 rev.bam -a assembly.fa -S bam",
    "pipeline -1 hic_reads_for.fastq.gz -2 hic_reads_rev.fastq.gz -a assembly_bowtie2_index -v",
    "network -1 for.bam -2 rev.bam -a assembly.fa",
    "partition -c contig_data_network.txt -n network.txt -a assembly.fa",
    "validation -c contig_data_partition.txt -n network.txt -a assembly.fa -f overlapping_bin",
    "network -1 for.bam -2 rev.bam -a assembly.fa -n None",
]


@pytest.mark.parametrize("command", COMMANDS_OF_4253551)
def test_the_commands_of_metator_1_are_unchanged(command):
    """A guard: the command lines of metaTOR 1.x that its tutorial gives are
    those of 2.0, read as they were (parsed here, not run): the input of a
    run from reads, alignments or pairs is -1 (-2, --start), as before."""
    name, *args = command.split()
    ctx = click.Context(mtcli.cli, info_name="metator")
    params = mtcli.cli.get_command(ctx, name).make_context(name, list(args), parent=ctx).params
    if name in ["network", "pipeline"]:
        assert params.get("map") is None and params["forward"] == args[args.index("-1") + 1]


def test_an_unknown_option_is_still_unknown(run_cli):
    """A guard: an option that metaTOR never had is refused as unknown, as
    before."""
    result = run_cli("pipeline", "-1", "a.pairs", "-S", "pair", "--bogus", "-a", "a.fa", exit_code=2)
    assert "No such option" in result.output and "--bogus" in result.output and "gone in metaTOR 2.0" not in result.output


@pytest.mark.parametrize("no_clean_up", [False, True])
def test_network_from_alignments_sorts_its_pairs(no_clean_up, stub_pipeline, tmp_path, run_cli):
    """metator network from alignments writes the pairs files sorted and
    indexed in its output directory, for the steps that follow it, as
    before; the tables of the alignments of the mates are kept with -N
    only, as every temporary file. Without -N, a guard: the stub reads the
    call of the previous version, which passed no keep_temp, as False."""
    out = tmp_path / "out"
    args = ["-1", "a1.bam", "-2", "a2.bam", "-a", TESTS_DATA / "assembly", "-S", "bam", "-o", out, "-T", tmp_path / "tmp"]
    run_cli("network", *args, *(["-N"] if no_clean_up else []))
    assert stub_pipeline["sort_pairs"] and stub_pipeline["keep_temp"] == no_clean_up
    assert stub_pipeline["pairs"] == [str(out / "alignment_0_sorted.pairs.gz")]
    assert (out / "alignment_0_sorted.pairs.gz").exists() and (out / "alignment_0_sorted.pairs.gz.px2").exists()


def test_network_from_pairs_counts_every_pair(stub_pipeline, tmp_path, run_cli):
    """A guard: metator network gives the pairs network no edge region,
    every pair counts (alignment_to_contacts at edge 0), as in the map. So
    did its default, --edge 0, before --edge was removed."""
    run_cli("network", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", tmp_path / "out", "-T", tmp_path / "tmp")
    alignment_files, edge_arg, *_ = stub_pipeline["alignment_to_contacts"]
    assert (alignment_files, edge_arg) == (["a.pairs"], 0)


@pytest.mark.parametrize(
    "args",
    [
        ["{map}", "--scaffold"],
        ["{map}", "-1", "a.pairs"],
        ["{missing}"],
        ["{collection}"],
        ["-S", "pair", "-1", "a.pairs", "-m", "abundance"],
        ["-S", "fastq", "-1", "a.fq"],
        ["-S", "pair", "-1", "a.pairs", "-p", "bad.prefix"],
        ["-S", "pair", "-1", "a.pairs", "-a", "missing.fa"],
    ],
    ids=[
        "scaffold",
        "map and -1",
        "map not there",
        "collection",
        "normalization without depth",
        "no reverse reads",
        "prefix",
        "no assembly",
    ],
)
def test_pipeline_refused_keeps_the_bins_of_an_earlier_run(args, tmp_path, run_cli):
    """A refused run with --force (a rebinning into the same directory)
    stops before replacing the bins of the earlier run."""
    names = {
        "map": _contig_map(tmp_path / "input.metacool"),
        "missing": tmp_path / "missing.metacool",
        "collection": _collection(tmp_path / "maps.mmetacool", "lib1"),
    }
    args = [arg.format(**names) for arg in args]
    out = tmp_path / "out"
    # The outputs of a run of 2.0, and of 1.x, named after bins.
    names = ["overlapping_clusters", "recursive_clusters", "final_clusters", "overlapping_bin", "final_bin_unscaffold"]
    bins = [out / name / "bin.fa" for name in names]
    for path in bins:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    (out / "bin_summary.txt").write_text("earlier\n")
    with pytest.raises(ValueError):
        run_cli("pipeline", "-a", TESTS_DATA / "assembly.fa", *args, "-o", out, "-F", "-T", tmp_path / "tmp")
    assert all(path.exists() for path in bins) and (out / "bin_summary.txt").exists()


@pytest.mark.parametrize(
    "existing, message",
    [
        ("recursive_clusters", "already existed."),
        ("recursive_bin", "already existed: an output of an earlier run, named as metaTOR 1.x names it"),
    ],
)
def test_pipeline_refuses_earlier_bins_before_the_index(existing, message, tmp_path, monkeypatch, run_cli):
    """Without --force, the bins of an earlier run are refused before the
    index of the assembly is built (an hour on a large metagenome)."""

    def generate_fasta_index(*args):
        raise AssertionError("index built before the refusal")

    monkeypatch.setattr("metator.io.generate_fasta_index", generate_fasta_index)
    out = tmp_path / "out"
    (out / existing).mkdir(parents=True)
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    with pytest.raises(ValueError):
        run_cli("pipeline", "-1", "a.fq", "-2", "b.fq", "-a", assembly, "-S", "fastq", "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"ERROR :: {out / existing} {message}" in log.read_text()


def test_pipeline_without_pyfastx_index(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """No bin passed the size threshold, so pyfastx never indexed the
    assembly: the run ends without removing an index that is not there (it
    used to fail at its very end, with FileNotFoundError)."""
    monkeypatch.setattr(
        "metator.validation.recursive_decontamination", lambda *args, **kwargs: os.path.join(args[8], "clustering_matrix")
    )
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    pairs = tmp_path / "a_sorted.pairs.gz"
    pairs.write_bytes(gzip.compress(b"## pairs format v1.0\n"))
    out = tmp_path / "out"
    run_cli("pipeline", "-1", pairs, "-a", assembly, "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    assert (out / "network.txt").exists() and not Path(f"{assembly}.fxi").exists()


def test_mge_without_pyfastx_index(tmp_path, monkeypatch, run_cli):
    """No MGE bin, so pyfastx never indexed the assembly: mge ends without
    removing an index that is not there."""
    monkeypatch.setattr("metator.mge.mge_binning", lambda **kwargs: None)
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    run_cli("mge", "-a", assembly, "-c", "tests_data/contig_data_final.txt", "-b", "tests_data/binning.txt",
            "-m", "tests_data/subset.mges.txt", "a.pairs")  # fmt: skip
    assert not Path(f"{assembly}.fxi").exists()


def test_qc_default_outdir(tmp_path, monkeypatch, run_cli):
    """Without --outdir, qc writes to ./metator_qc_<prefix>."""
    calls = []
    monkeypatch.setattr("metator.quality_check.quality_check", lambda *args, **kwargs: calls.append((args, kwargs)))
    for name in ["cluster_summary.txt", "contig_data_final.txt"]:
        (tmp_path / name).write_text("x")
    run_cli("qc", "-a", "assembly.fa", "-e", "DpnII,HinfI", "-p", "sample", "-t", "3,4", "a.pairs", "b.pairs")
    ((args, kwargs),) = calls
    assert args[:5] == (
        "./contig_data_final.txt",
        "./cluster_summary.txt",
        "assembly.fa",
        ["a.pairs", "b.pairs"],
        "./metator_qc_sample",
    )
    assert args[6:] == ("sample", False)
    assert kwargs == {"enzyme": ["DpnII", "HinfI"], "threshold": ["3", "4"], "metacool": None}
    assert len(list((tmp_path / "metator_qc_sample").glob("metator_qc_*.log"))) == 1


def test_scaffold_default_outputs(tmp_path, monkeypatch, run_cli):
    calls = []
    monkeypatch.setattr("metator.scaffold.get_scaffolds", lambda *args, **kwargs: calls.append((args, kwargs)))
    run_cli("scaffold", "-b", "bin_1", "-i", "bin_1.fa", "a.pairs")
    ((args, kwargs),) = calls
    assert args == ("bin_1", "bin_1.fa", ["a.pairs"], "./bin_1_scaffolded.fa", "./bin_1_info_frags.txt")
    assert kwargs == {"threshold": 0.05, "threads": 1, "window_size": 5000, "junctions": "NNNNN"}
    assert len(list(tmp_path.glob("metator_scaffold_*.log"))) == 1


def test_mge_default_outdir(tmp_path, monkeypatch, run_cli):
    """Without --outdir, mge writes to the current directory; the pyfastx
    index of the assembly is removed at the end."""
    calls = []

    def mge_binning(**kwargs):
        calls.append(kwargs)
        Path(f"{kwargs['fasta_mges_contigs']}.fxi").write_text("x")

    monkeypatch.setattr("metator.mge.mge_binning", mge_binning)
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    run_cli("mge", "-a", assembly, "-c", "tests_data/contig_data_final.txt", "-b", "tests_data/binning.txt",
            "-m", "tests_data/subset.mges.txt", "a.pairs")  # fmt: skip
    (kwargs,) = calls
    assert kwargs["out_dir"] == "." and kwargs["pairs_files"] == ["a.pairs"]
    assert kwargs["threshold_bin"] == 0.8 and kwargs["random"] is False
    assert not Path(f"{assembly}.fxi").exists()
    assert not (tmp_path / "tmp").exists()
