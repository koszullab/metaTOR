# Tests of the metator command line (metator.cli): parity with the docopt
# command line it replaced, help, version and usage errors, and the
# bookkeeping of the commands (log files, temporary directories, default
# outputs), with the heavy steps replaced by stubs.

import json
import os
import re
import shutil
import sys
from pathlib import Path

import click
import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import metator.cli as mtcli
import metator.io as mio
import metator.metacool as mtmc
import metator.report as mtr
from metator.version import __version__
from tests.locks import described, held_lock, killed_lock, lock_file, refusal

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"

# The writer itself, which the stub pipeline replaces for the map of the
# clusters.
mtmc_write_map = mtmc.write_map

# The docopt command line that rich-click replaced (that of 1.3.10, with the
# mge options renamed and host removed), as docopt read it: see the commit
# that added the file. Commands and options added since are added to it:
# metacool, the --output-map option of pipeline (and --no-metacool, gone
# since), and the --output-map option of partition and validation, which
# --inplace replaced: a run writes its map to its output directory. The 2.0
# command line
# changed it on purpose (see docs/migration.md): the map that network and
# pipeline start from is their first argument (MAP), required, and their
# options of reads, alignments and pairs (--forward, --reverse, --start,
# --aligner, --aligner-mode, --min-quality) and --edge are gone; the map
# that partition and validation update is their first argument too; the
# commands qc, contactmap, metacool, scaffold and pairs are gone, and so are
# the --scaffold and --junctions options of pipeline, --self-contacts of
# network and --cluster-matrix of partition, validation and pipeline. 2.0
# adds dump, which writes the tables of a run from its map, and host again,
# which finds the hosts of the MGE clusters of a map (1.x's host read tables).
# network and pipeline then took pairs and BAM files again, as their
# arguments (INPUT..., the map or the files that the run makes it from), with
# -q/--min-quality for the BAM files, and network -t/--threads to read them;
# and pipeline the --mges, --mge-threshold-bin, --mge-min-share and
# --mge-min-contigs options, with which it runs mge and host after the
# validation; then partition, validation, pipeline and mge took --seed,
# the seed of their random numbers. export writes the outputs of a run
# again from its map, and report its report.
INVENTORY = json.loads((TESTS_DATA / "cli_inventory.json").read_text())["metator"]

# Options whose value is one of a list, wherever they appear.
CHOICES = {
    "--algorithm": ["louvain", "leiden"],
    "--normalization": ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"],
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
    """A command line with every required option and argument of a command, with
    made-up values."""
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
    """The help shows the inputs, the arguments: the map, which hicstuff or
    metacooler make, or the pairs or BAM files that the run makes it from;
    and that every pair counts."""
    output = " ".join(run_cli(command, "--help").output.split())
    assert f"Usage: metator {command} [OPTIONS] INPUT..." in output
    assert f"metator {command} INPUT... -a ASSEMBLY" in output and "-1 FORWARD" not in output
    assert "hicstuff pipeline -M metacool" in output and "metacooler create pairs" in output
    assert "the pairs files or the BAM files" in output
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
    """A required argument left out: the map, which the arguments come
    after."""
    args = _required_args(command)
    del args[args.index("p.pairs")]
    # network and pipeline take one input or more (INPUT...), the others a map.
    missing = "INPUT..." if INVENTORY["commands"][command]["positional"][0]["repeatable"] else "MAP"
    result = run_cli(command, *args, exit_code=2)
    assert f"Missing argument '{missing}'" in result.output


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
    result = run_cli("partition", "-t", value, "p.metacool", exit_code=2)
    assert f"'{value}' is not a valid integer" in result.output


def test_unknown_option(run_cli):
    """Long options cannot be abbreviated (docopt accepted --thread)."""
    result = run_cli("network", "--thread", "2", "p.metacool", exit_code=2)
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
    with mtcli.private_tmp_dir("", clean_up=False) as tmp_dir:
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
    """The log starts as hicstuff's does: the version, the date and the command
    line (cmd), then a line of dashes, with no empty line; it ends with what
    the run wrote and how long it took, as hicstuff's does too."""
    log = tmp_path / "run.log"
    monkeypatch.setattr(sys, "argv", ["/usr/bin/metator", "partition", "-p", "my sample", "a.metacool"])
    started = mtcli.generate_log_header(str(log))
    mtcli.log_end(started, "Contact map with the partition written", "out/a.metacool")
    lines = log.read_text().splitlines()
    assert lines[0] == f"## MetaTOR: v{__version__} log file"
    assert lines[1].startswith("## date: ")
    assert lines[2] == "## cmd: metator partition -p 'my sample' a.metacool"
    assert lines[3:-1] == ["---"]
    assert re.fullmatch(r"\S+ :: INFO :: Contact map with the partition written after 0h 0m [0-9]s: out/a\.metacool", lines[-1])


def test_log_header_lists_the_parameters(stub_pipeline, tmp_path, run_cli):
    """The header of the log of a run lists the value of each parameter of the
    command, by its long option, in the order of the help, as hicstuff lists
    its own: as given, by default, or as the run resolved them (its output
    directory)."""
    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=True)
    run_cli("partition", metacool, "-i", "4", "-T", tmp_path / "tmp")
    (log,) = tmp_path.glob("metator_partition_*.log")
    lines = log.read_text().splitlines()
    assert lines[2] == f"## cmd: metator partition {metacool} -i 4 -T {tmp_path / 'tmp'}"
    assert lines[3 : lines.index("---")] == [
        f"## map: {metacool}",
        "## algorithm: louvain",
        "## iterations: 4",
        "## overlap: 80",
        "## res-param: 1.0",
        "## outdir: ./",
        "## prefix: None",
        "## force: False",
        f"## tmpdir: {tmp_path / 'tmp'}",
        "## threads: 1",
        "## no-clean-up: False",
        "## inplace: False",
        "## seed: None",
    ]
    assert "" not in lines
    assert re.search(r" :: INFO :: Contact map with the partition written after 0h 0m [0-9]s: input\.metacool$", lines[-1])


# The bookkeeping of the commands, with the heavy steps replaced by stubs.


# The contacts of a map as metator.network reads them, which stub_pipeline
# records: the maps that the tests make read them unrecorded.
CONTIG_CONTACTS = mtcli.mtn.contig_contacts

# What metator.network measures on the assembly, which stub_pipeline replaces:
# the tests of what it refuses run it.
CREATE_CONTIG_DATA = mtcli.mtn.create_contig_data


def _no_step(calls):
    """No step of the stub pipeline ran (see stub_pipeline): a partition or
    validation run reads the contacts of its map as it reads the network stored
    in it (see metator.network.stored_network), before anything else."""
    return not set(calls) - {"contig_contacts"}


class _AnyContig(dict):
    """The contig data of the stub create_contig_data: every contig of the map, of
    1000 bp, measured on nothing but its GC content (no depth file, no enzyme),
    as create_contig_data gives a contig."""

    def __missing__(self, name):
        return {"length": 1000, "GC": 50.0, "coverage": "-", "RS": "-"}


@pytest.fixture
def stub_pipeline(monkeypatch, tmp_path):
    """Replace the steps of the pipeline by stubs that write their outputs and
    record their arguments."""
    calls = {}
    contig_contacts = mtcli.mtn.contig_contacts

    def create_contig_data(*args):
        calls["create_contig_data"] = args
        return _AnyContig()

    def contig_data_from_map(*args):
        # A run from a map reads no column of it (see
        # test_pipeline_reads_the_columns_of_its_map).
        calls["contig_data_from_map"] = args
        return []

    def recorded_contig_contacts(metacool, *args, **kwargs):
        calls["contig_contacts"] = (metacool, *args)
        return contig_contacts(metacool, *args, **kwargs)

    def partition(*args):
        calls["partition"] = args
        return _one_cluster(args[1])

    def recursive_decontamination(*args, **kwargs):
        calls["recursive_decontamination"] = args
        calls["recursive_decontamination_kwargs"] = kwargs
        # The contig data table of the partition, in the temporary directory
        # of the run, which is removed at its end.
        calls["partition_table"] = pd.read_csv(args[2], sep="\t", dtype={"Core_cluster_ID": str})
        # No final cluster: every contig of the partition in none, and no
        # column of a summary.
        return pd.DataFrame({"Name": calls["partition_table"]["Name"], "Final_cluster": "ND"}), mio.SummaryLayout([], [], {})

    def write_map(metacool, outfile, **stores):
        # The map of the clusters of a validation or pipeline run, recorded;
        # the others (network, partition, mge) written for real.
        if stores.get("clusters") is None:
            return mtmc_write_map(metacool, outfile, **stores)
        # The input map must still be there.
        calls["write_map"] = (metacool, os.path.exists(metacool), outfile, stores)

    def run_export(*args, **kwargs):
        # The step export of a pipeline run, its outputs written from the map
        # of its clusters, which the stub validation does not write: recorded.
        calls["run_export"] = (args, kwargs)

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    monkeypatch.setattr("metator.cli.run_export", run_export)
    monkeypatch.setattr("metator.network.contig_data_from_map", contig_data_from_map)
    monkeypatch.setattr("metator.network.contig_contacts", recorded_contig_contacts)
    monkeypatch.setattr("metator.partition.partition", partition)
    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)
    monkeypatch.setattr("metator.metacool.write_map", write_map)
    return calls


def _one_cluster(network):
    """The partition of the stub partition: one core and overlapping cluster
    (metator_00001_00000) of the contigs of the network, in the order of their
    IDs (none for a network without an edge)."""
    import metator.clusters as mtcl

    contigs = sorted({int(contig) for contig in network.contig1} | {int(contig) for contig in network.contig2})
    clusters = {1: contigs} if contigs else {}
    return mtcl.Partition(core=dict(clusters), overlapping=dict(clusters), names={1: "metator_00001_00000"} if contigs else {})


def _assembly(tmp_path):
    """A copy of an assembly of tests_data, of the test's own."""
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    return assembly


def _final_tables(monkeypatch, contig, final_bin="metator_00001_00000", summary=None):
    """Have the stub validation give the final clusters: the contig ``contig`` in
    the bin ``final_bin`` (ND: in none) and c2 in none, write the FASTA file of
    the bin in the directory of the final clusters, and the bin summary
    (``summary``, if given) in the directory of its tables, the temporary
    directory of the run."""
    stub = mtcli.mtv.recursive_decontamination

    def recursive_decontamination(*args, **kwargs):
        stub(*args, **kwargs)
        rows = "" if final_bin == "ND" else f"{final_bin}\t0.9\n"
        if final_bin != "ND":
            Path(args[3], f"{final_bin}.fa").write_text(f">{contig}\nACGT\n")
        # The summary as the run writes it, and at full precision, which the
        # map stores (the same here).
        for name in ("cluster_summary.txt", mio.SUMMARY_FULL_FILE):
            Path(args[7], name).write_text(summary or f"\tCompleteness\n{rows}")
        layout = mio.SummaryLayout([] if final_bin == "ND" else [final_bin], ["Completeness"], {"Completeness": "%.2f"})
        return pd.DataFrame({"Name": [contig, "c2"], "Final_cluster": [final_bin, "ND"]}), layout

    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)


def _contig_data(monkeypatch, contigs):
    """Have the stub create_contig_data give the contig data of these contigs:
    their length and GC content, by name, and nothing else measured (no depth
    file, no enzyme), as create_contig_data does."""

    def create_contig_data(*args):
        return {name: {"length": length, "GC": gc, "coverage": "-", "RS": "-"} for name, (length, gc) in contigs.items()}

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)


# The contigs of the stub validation's tables (see _final_tables).
CONTIGS = {"c1": (100, 40.5), "c2": (50, 35.0)}


def _pairs(tmp_path, contig="c1"):
    """A pairs file of two pairs, one between contig and c2, one on contig."""
    pairs = tmp_path / "a.pairs"
    pairs.write_text(f"## pairs format v1.0\nr1\t{contig}\t10\tc2\t20\t+\t-\nr2\t{contig}\t10\t{contig}\t90\t+\t-\n")
    return pairs


NON_ASCII_ASSEMBLY = ">NODE_1 a contig\n" + "ACGT" * 30 + "\n>NODE_522_é\n" + "ACGT" * 30 + "\n"


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize("compression", ["bzip2", "zip"])
def test_a_run_refuses_an_assembly_that_metacooler_cannot_read(compression, command, stub_pipeline, tmp_path, run_cli, logged):
    """An assembly compressed with bzip2 or zip, which metaTOR reads but
    metacooler, which writes the FASTA files of the clusters, does not, is
    refused before any work: the clusters of an earlier run stay, even with
    -F (pipeline)."""
    import bz2
    import zipfile

    fasta = tmp_path / "assembly.fa"
    if compression == "bzip2":
        fasta.write_bytes(bz2.compress(b">c1\nACGT\n"))
    else:
        with zipfile.ZipFile(fasta, "w") as archive:
            archive.writestr("assembly.fa", ">c1\nACGT\n")
    earlier = tmp_path / "out" / "fasta" / "mags_with-mges" / "bin.fa"
    earlier.parent.mkdir(parents=True)
    earlier.write_text(">old\nA\n")
    metacool = _contig_map(tmp_path / "map.metacool")
    with pytest.raises(ValueError):
        force = ["-F"] if command == "pipeline" else []
        run_cli(command, metacool, "-a", fasta, *force, "-o", tmp_path / "out", "-T", tmp_path / "tmp")
    assert earlier.exists() and not stub_pipeline
    message = f"{fasta} is compressed with {compression}: give the assembly as a FASTA file, plain or compressed with gzip."
    assert ("ERROR", message) in logged()


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize("inputs", ["map", "pairs"])
def test_a_run_refuses_an_assembly_that_names_a_contig_by_several_records(
    inputs, command, stub_pipeline, tmp_path, monkeypatch, run_cli, logged
):
    """An assembly that names a contig by several records, the last of which gave
    the contig its length and GC content, is refused before any work: no map
    is made from the pairs, the map of the output directory is not written,
    and the clusters of an earlier run stay, even with -F (pipeline)."""
    monkeypatch.setattr("metator.network.create_contig_data", CREATE_CONTIG_DATA)
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">c1\n" + "A" * 100 + "\n>c2\n" + "C" * 50 + "\n>c1 again\n" + "G" * 100 + "\n")
    source = _contig_map(tmp_path / "map.metacool") if inputs == "map" else _pairs(tmp_path)
    out = tmp_path / "out"
    earlier = out / "fasta" / "mags_with-mges" / "bin.fa"
    earlier.parent.mkdir(parents=True)
    earlier.write_text(">old\nA\n")
    what = f"1 contig(s) are named by several records of the assembly {fasta}"
    with pytest.raises(ValueError, match=re.escape(f"{what}.")):
        force = ["-F"] if command == "pipeline" else []
        run_cli(command, source, "-a", fasta, *force, "-o", out, "-T", tmp_path / "tmp")
    assert _no_step(stub_pipeline) and not list(out.glob("*.metacool")) and earlier.exists()
    message = f"{what}: 'c1'. A contig is read by its name, the first word of its header: give each record a name of its own."
    assert ("ERROR", message) in logged()


def test_pipeline_from_a_map_keeps_the_contig_names(stub_pipeline, tmp_path, run_cli):
    """A guard: a run from a map does not check the contig names of the assembly,
    as it builds no map of the contigs (nor does the map of its clusters: see
    test_a_non_ascii_contig_that_the_map_lacks_is_left_out in
    test_metacool.py)."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(NON_ASCII_ASSEMBLY)
    metacool = _contig_map(tmp_path / "map.metacool")
    run_cli("pipeline", metacool, "-a", fasta, "-o", tmp_path / "out", "-T", tmp_path / "tmp")
    assert stub_pipeline["contig_contacts"][0] == str(metacool)


def test_pipeline_from_a_map_leaves_another_map_of_its_output_directory(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """From a map stored elsewhere: the run writes its map to the output
    directory, under the name of that map, and a map of another name that an
    earlier run left there (that of 1.4.0) is not this run's map: it is left
    as it is."""
    out = tmp_path / "out"
    out.mkdir()
    earlier = out / "metator.metacool"
    earlier.write_text("earlier run")
    metacool = _contig_map(tmp_path / "input.metacool")
    run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert stub_pipeline["write_map"][2] == str(out / "input.metacool") and earlier.read_text() == "earlier run"
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert "Removing" not in logged
    message = (
        f"INFO :: The run writes its map to {out / 'input.metacool'}, and leaves {metacool} as it is (--inplace updates it)."
    )
    assert message in logged


def test_pipeline_map_stores_the_normalization_of_the_network(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """With the real writers: the map of the clusters stores the normalization of
    the network of the run, empirical_hit by default, as the weights
    metator_empirical_hit of its bins, one bin per contig: 1/sqrt(hits), the
    hits counted on the contacts of the map (c1: 3, c2: 1), which its Hit
    column stores."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _contig_data(monkeypatch, CONTIGS)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "map.metacool")
    out = tmp_path / "out"
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    (log,) = out.glob("metator_pipeline_*.log")
    assert "WARNING" not in log.read_text()
    import metacooler

    path = str(out / "map.metacool")
    names = list(metacooler.MetaCooler(path).contigs.name)
    weights = metacooler.read_weights(path, "metator_empirical_hit")
    assert dict(zip(names, weights, strict=True)) == pytest.approx({"c1": 1 / 3**0.5, "c2": 1.0})
    hits = metacooler.MetaCooler(path).contigs.set_index("name")["Hit"]
    assert hits.to_dict() == {"c1": 3, "c2": 1}


@pytest.mark.parametrize("case", ["a contig of another assembly", "no contact"])
def test_pipeline_network_of_the_pairs_of_its_map(case, network_from_a_map, tmp_path, run_cli):
    """The network of a run is that of the pairs that its map holds."""
    import metacooler
    import pandas as pd

    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">c1\n" + "ACGT" * 250 + "\n>c2\n" + "GGCC" * 100 + "\n")
    pairs = tmp_path / "a.pairs"
    rows = ["r1\tc1\t10\tc2\t20", "r2\tc1\t10\tzz\t90", "r3\tc1\t30\tc1\t500"]
    if case == "no contact":
        # A read past the end of c2.
        rows = ["r1\tc1\t10\tc2\t900"]
    pairs.write_text("## pairs format v1.0\n" + "".join(f"{row}\t+\t-\n" for row in rows))
    metacool = tmp_path / "map.metacool"
    bins = pd.DataFrame({"chrom": ["c1", "c2"], "start": 0, "end": [1000, 400]})
    metacooler.from_pairs(str(metacool), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    out = tmp_path / "out"
    args = [metacool, "-a", assembly, "-o", out, "-T", tmp_path / "tmp"]
    if case == "no contact":
        message = f"{metacool} holds no contact: there is no network to build."
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("pipeline", *args)
        (log,) = out.glob("metator_pipeline_*.log")
        assert f"ERROR :: {message}\n" in log.read_text()
        return
    with pytest.raises(RuntimeError, match=r"^network built$"):
        run_cli("pipeline", *args)
    assert network_from_a_map == [str(metacool)]
    # c1 and c2 (IDs 1 and 2: their positions in the map), one contact,
    # normalized by their hits (3, 1).
    contig_data, network = network_from_a_map.partitioned
    assert list(contig_data.Name) == ["c1", "c2"] and list(contig_data.Hit) == [3, 1]
    assert network.to_dict("list") == {"contig1": [1], "contig2": [2], "weight": [pytest.approx(1 / 3**0.5, rel=1e-15)]}
    # No file of the network.
    assert not list(out.glob("*network*"))


def test_pipeline_from_a_metacool(stub_pipeline, tmp_path, run_cli):
    """From a contact map: the network is read from it, and the map of the
    clusters is written back from it, once: what the three steps store, in one
    edit session."""
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool")
    before = metacool.read_bytes()
    run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-t", "3", "-m", "length",
            "-d", TESTS_DATA / "depth.txt", "-o", out, "-T", tmp_path / "tmp")  # fmt: skip

    calls = stub_pipeline
    fasta, depth, enzyme = calls["create_contig_data"]
    assert (fasta, depth, enzyme) == (str(tmp_path / "assembly.fa"), str(TESTS_DATA / "depth.txt"), None)
    assert calls["contig_contacts"] == (str(metacool),)
    # The network of the map, in memory: c1 and c2 (the stub contig data
    # gives them 1000 bp), one contact; the contig data of its partition,
    # which the validation reads, in the temporary directory of the run.
    network = calls["partition"][1]
    assert network.to_dict("list") == {
        "contig1": [1], "contig2": [2], "weight": [pytest.approx((1000 / 3) ** 0.5 * 1000**0.5, rel=1e-15)]
    }  # fmt: skip
    assert calls["partition"][5] == 3
    table = calls["partition_table"]
    assert list(table.Name) == ["c1", "c2"] and list(table.Hit) == [3, 1]
    assert list(table.Core_cluster_ID) == ["00001", "00001"] and list(table.Overlapping_cluster_size) == [2000, 2000]
    # The validation works in the temporary directory of the run: the
    # contig data of the partition, the tables of the validation and its
    # recursive clusters; the final clusters go to fasta/mags_with-mges of
    # the output directory.
    tables = Path(calls["recursive_decontamination"][7])
    assert Path(calls["recursive_decontamination"][2]).parent == tables and tables.parent == tmp_path / "tmp"
    assert calls["recursive_decontamination"][3] == str(out / "fasta" / "mags_with-mges")
    assert Path(calls["recursive_decontamination"][9]) == tables / "recursive_clusters"
    assert calls["recursive_decontamination_kwargs"]["metacool"] == str(metacool)
    # The plots of the final clusters go to plots/ of the output directory,
    # which the validation makes (the stub validation draws none).
    assert calls["recursive_decontamination_kwargs"]["plot_outdir"] == str(out) and not (out / "plots").exists()
    # The outputs of the run, exported from its map in its output directory
    # (in place of those of its steps), from the assembly (plain: as it is),
    # in its temporary directory.
    (metacool_file, assembly, outdir, tmp_dir), options = calls["run_export"]
    assert (metacool_file, assembly, outdir, options) == (str(metacool), str(tmp_path / "assembly.fa"), str(out), {"force": True})
    assert Path(tmp_dir).parent == tmp_path / "tmp"
    assert not (out / "overlapping_clusters").exists() and not (out / "contig_data_partition.txt").exists()
    assert metacool.read_bytes() == before
    # No table of the network: it is in the map.
    assert not list(out.glob("*network*"))
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert "Removing" not in logged
    source, existed, outfile, stores = calls["write_map"]
    assert (source, existed, outfile) == (str(metacool), True, str(metacool))
    # The map as the run read it under its lock.
    assert stores["force"] and os.path.samefile(stores["planned"].filename, metacool)
    # The network, stored as metator network stores it, and its
    # normalization.
    network = stores["network"]
    assert network.normalization == "length"
    pd.testing.assert_frame_equal(network.table, table[network.table.columns], check_dtype=False)
    assert network.entry["command"].startswith(f"metator pipeline {metacool} -a") and network.entry["version"] == __version__
    parameters = network.entry["parameters"]
    assert parameters["map"] == str(metacool) and "start" not in parameters
    assert Path(network.tmp_dir).parent == tmp_path / "tmp"
    # What its partition was, as metator partition stores it, and its
    # validation, as metator validation stores it.
    assert stores["partition"]["normalization"] == "length"
    assert stores["partition"]["command"] == network.entry["command"] and stores["partition"]["version"] == __version__
    clusters = stores["clusters"]
    assert clusters.entry["command"] == network.entry["command"] and clusters.entry["version"] == __version__
    assert clusters.entry["normalization"] == "length" and clusters.entry["summary_columns"] == []
    # The coverage of the depth file is not measured by the stub contig data.
    assert clusters.entry["read"] == {"contigs": ["GC_content", "Hit"], "bins": ["metator_length"],
                                      "rounds": ["core", "overlapping"]}  # fmt: skip
    # The stub validation clusters no contig, and writes no table.
    assert clusters.assignment == {} and clusters.rounds is None and clusters.summary == str(tables / mio.SUMMARY_FULL_FILE)
    # No column was read from the map (the stub reads none).
    assert network.entry["read"] == {"contigs": []}
    assert stores["partition"]["read"] == {"contigs": ["GC_content", "Hit"], "bins": ["metator_length"]}
    # The map of the output directory, updated in place, as the log says.
    in_place = f"INFO :: The run updates {metacool} in place: it is the map of its output directory (give another one"
    assert in_place in logged and "is the input map" not in logged


# The steps of a pipeline run, as its log marks them, without --mges.
STEPS = [("1", "4", "network"), ("2", "4", "partition"), ("3", "4", "validation"), ("4", "4", "export")]


def test_pipeline_without_mges_has_no_steps_mge_and_host(stub_pipeline, tmp_path, run_cli):
    """Without --mges, a pipeline run from a map has four steps ([1/4]
    network... [4/4] export), and says, once it has written the map of its
    clusters, that it has no steps mge and host, and how to run them: with
    --mges. It writes no MGE fasta and no hosts, and then exports its outputs,
    from the map of its clusters."""
    metacool = _contig_map(tmp_path / "input.metacool")
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    out = tmp_path / "out"
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert re.findall(r" :: INFO :: \[([0-9])/([0-9])\] ([a-z]+): ", logged) == STEPS
    stop = (
        "INFO :: No list of MGE contigs (--mges): the run has no steps mge and host, which cluster the MGE contigs and "
        "find their hosts: give --mges to run them.\n"
    )
    assert logged.index("[3/4] validation") < logged.index(stop) < logged.index("[4/4] export")
    assert stub_pipeline["write_map"] and stub_pipeline["run_export"]
    assert not (out / "fasta" / "mges").exists() and not (out / "tables").exists()


# What the steps mge and host of a pipeline run with --mges would refuse
# once its clusters are made, refused before any work, the reason logged as an
# error: a list of MGE contigs that cannot be read, or that is the map of the
# output directory; a map without the contact statistics that mge reads, or
# that holds none of the MGE contigs; and the summary, the FASTA file of the MGE
# clusters and the hosts of an earlier version in the output directory, which
# mge and host refuse.
PIPELINE_MGES_ERRORS = {
    "a list that is not there": (
        "{missing}",
        "Cannot read the list of MGE contigs {missing} (No such file or directory): give a text file, one contig name per line.",
    ),
    "a directory as the list": (
        "{directory}",
        "Cannot read the list of MGE contigs {directory} (Is a directory): give a text file, one contig name per line.",
    ),
    "the list as the map of the output directory": (
        "{out_map}",
        "{out_map}, where the run writes its map, is {out_map}, an input of the run: give the run another output "
        "directory with -o.",
    ),
    "no MGE contig in the map": ("{other_mges}", "None of the MGE contigs is in {map}: there is no MGE contig to cluster."),
    "a map without contact statistics": (
        "{mges}",
        "{map} holds no contact statistics of the reads of its contigs (metacooler's QC block, with "
        "qc_intra_contig_long in its contigs table): metator mge counts the contacts within each MGE contig between "
        "reads more than 1 kb apart, which the pixels of a map do not tell. ",
    ),
    "contact statistics of the pixels": (
        "{mges}",
        "{map} has the contact statistics of its contigs measured from its pixels (QC source 'pixels'), not from the "
        "pairs: metator mge counts",
    ),
    "the summary of an earlier version": (
        "{mges}",
        "{out} holds results of an earlier run, which this run does not write and does not replace: "
        "'mges_bin_summary.tsv'. Move them away, or give the run another output directory with -o.",
    ),
    "the MGE FASTA file and the hosts of an earlier version": (
        "{mges}",
        "{out} holds results of an earlier run, which this run does not write and does not replace: "
        "'mge_hosts.txt', 'mges_binned.fa'. Move them away, or give the run another output directory with -o.",
    ),
}


@pytest.mark.parametrize("case", list(PIPELINE_MGES_ERRORS))
def test_pipeline_with_mges_refused_before_any_work(case, stub_pipeline, tmp_path, run_cli):
    """What the steps mge and host of a pipeline run with --mges would refuse,
    after hours of work, is refused before any (see PIPELINE_MGES_ERRORS): the
    map and the output directory are left as they were."""
    import h5py

    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "out"
    out.mkdir()
    if case == "a map without contact statistics":
        with h5py.File(metacool, "r+") as handle:
            metadata = json.loads(handle.attrs["metadata"])
            del metadata["qc"]
            handle.attrs["metadata"] = json.dumps(metadata)
    elif case == "contact statistics of the pixels":
        with h5py.File(metacool, "r+") as handle:
            metadata = json.loads(handle.attrs["metadata"])
            metadata["qc"]["source"] = "pixels"
            handle.attrs["metadata"] = json.dumps(metadata)
    elif case == "the summary of an earlier version":
        (out / "mges_bin_summary.tsv").write_text("an earlier summary")
    elif case == "the MGE FASTA file and the hosts of an earlier version":
        (out / "mges_binned.fa").write_text(">MetaTOR_MGE_00001\nACGT\n")
        (out / "mge_hosts.txt").write_text("mge\tmag\n")
    elif case == "the list as the map of the output directory":
        shutil.copyfile(metacool, out / "input.metacool")
    mges = tmp_path / "mges.txt"
    mges.write_text("c2\n")
    other_mges = tmp_path / "other_mges.txt"
    other_mges.write_text("c3\nc4\n")
    paths = {
        "map": metacool,
        "out": out,
        "out_map": out / "input.metacool",
        "mges": mges,
        "other_mges": other_mges,
        "missing": tmp_path / "missing.txt",
        "directory": tmp_path,
    }
    given, message = (text.format(**paths) for text in PIPELINE_MGES_ERRORS[case])
    before = metacool.read_bytes()
    with pytest.raises(ValueError):
        run_cli("pipeline", *_map_run(metacool, tmp_path, "--mges", given, "-F"))
    assert _no_step(stub_pipeline) and metacool.read_bytes() == before
    assert not (out / "fasta").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}" in log.read_text()


def test_pipeline_with_mges_warns_that_it_replaces_the_round_mge(tmp_path, logged):
    """pipeline MAP --mges, which updates MAP in place, warns before any work that
    it replaces the round mge of the map too, which its step mge stores, with
    the clusters and the rounds of the binning that it replaces."""
    metacool = _map_of_an_earlier_run(tmp_path / "input.metacool")
    mtcli.warn_of_what_is_replaced("pipeline", str(metacool), rounds=["mge"], inplace=True)
    message = (
        f"{metacool} holds 1 cluster(s) and the clustering round(s) overlapping, final, mge, which this run replaces, in "
        f"place: leave out --inplace to write the map of this run to another file, and leave {metacool} as it is."
    )
    assert [record for record in logged() if record[0] == "WARNING"] == [("WARNING", message)]


def _contig_map(path, lengths=None, network=False, normalization="empirical_hit", partition=False):
    """A .metacool map of the contigs c1 and c2 of the stub validation's tables
    (see _final_tables), made before any binning: every contig in __unbinned__,
    one bin per contig."""
    import metacooler
    import pandas as pd

    lengths = lengths or {"c1": 100, "c2": 50}
    pairs = path.with_suffix(".pairs")
    pairs.write_text("## pairs format v1.0\nr1\tc1\t10\tc2\t20\t+\t-\nr2\tc1\t10\tc1\t90\t+\t-\n")
    bins = pd.DataFrame({"chrom": list(lengths), "start": 0, "end": list(lengths.values())})
    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    pairs.unlink()
    if network:
        _store_contig_network(path, normalization)
    if partition:
        import metator.clusters as mtcl

        names = [str(name) for name in metacooler.MetaCooler(str(path)).contigs["name"]]
        clusters = _one_cluster(pd.DataFrame({"contig1": [1], "contig2": [2]}))
        entry = {"version": __version__, "normalization": normalization}
        if network:
            entry[mtmc.NETWORK_DIGEST] = mtcli.mtn.stored_network(str(path)).digest
        mtmc.write_map(
            str(path),
            str(path),
            partition=entry,
            rounds=mtcl.partition_rounds(clusters, names),
            round_names=list(mtmc.PARTITION_ROUNDS),
            force=True,
        )
    return path


def _store_contig_network(path, normalization):
    """Store the network of the contigs of the map at path (see _contig_map) in
    it, as metator network stores it (see metator.metacool.write_map), with
    normalization: the network of an earlier run is replaced, its partition
    kept."""
    import metacooler

    names = metacooler.MetaCooler(str(path)).contigs["name"]
    contacts = CONTIG_CONTACTS(str(path))
    table = mtcli.mtn.contig_data_table(_AnyContig(), names, mtcli.mtn.contig_hits(contacts, len(names)))
    mtmc.write_map(str(path), str(path), network=mtmc.Network(table, normalization), force=True)


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
    """Before any work, from a map: --enzyme and --depth, or the columns of the
    map that give values."""
    metacool = _map_with_columns(tmp_path, columns)
    args = (normalization, options.get("enzyme"), options.get("depth"), str(metacool))
    if message is None:
        mtcli.check_normalization(*args)
        assert logged() == []
        return
    with pytest.raises(ValueError):
        mtcli.check_normalization(*args)
    assert logged() == [("ERROR", message.format(map=metacool))]


def test_check_normalization_once_the_map_is_read(tmp_path, logged):
    """Once the run has read the columns of its map (read): those that it could
    read, or a refusal for a column that lacks contigs."""
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
    """A map that cannot be read is refused before any work when a normalization
    needs its columns, its reason logged; otherwise the run finds it at its
    network, as before."""
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
    contig_data = {name: {"length": length, "GC": 50.0, "coverage": "-", "RS": "-"}
                   for name, length in [("c1", 100), ("c2", 50)]}  # fmt: skip
    monkeypatch.setattr("metator.network.create_contig_data", lambda *args: contig_data)
    metacool = _map_with_columns(tmp_path, SITES_AND_COVERAGE)
    out = tmp_path / "out"
    run_cli(command, metacool, "-a", _assembly(tmp_path), "-n" if command == "network" else "-m",
            "theoritical_hit", "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert {name: (data["RS"], data["coverage"]) for name, data in contig_data.items()} == {"c1": (3, 1.5), "c2": (2, 0.5)}
    # Where the values that the map stores came from: the entry network of
    # its metadata.
    if command == "pipeline":
        stores = stub_pipeline["write_map"][3]
        assert stores["network"].normalization == stores["partition"]["normalization"] == "theoritical_hit"
        entry = stores["network"].entry
    else:
        import metacooler

        entry = metacooler.MetaCooler(str(out / metacool.name)).metadata["metator"]["network"]
        assert entry["normalization"] == "theoritical_hit"
    assert entry["read"] == {"contigs": ["Restriction_site", "Shotgun_coverage"]}


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_run_refused_for_a_column_that_lacks_contigs(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A column of the map that lacks a contig of the assembly cannot give
    the normalization what it needs: refused, with the warning that names
    the contig, before the network."""
    monkeypatch.setattr("metator.network.contig_data_from_map", mtn_contig_data_from_map)
    contig_data = {name: {"length": length, "GC": 50.0, "coverage": "-", "RS": "-"}
                   for name, length in [("c1", 100), ("c2", 50), ("c3", 30)]}  # fmt: skip
    monkeypatch.setattr("metator.network.create_contig_data", lambda *args: contig_data)
    metacool = _map_with_columns(tmp_path, {"Restriction_site": [3, 2]})
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli(command, metacool, "-a", _assembly(tmp_path), "-n" if command == "network" else "-m",
                "RS", "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert "contig_contacts" not in stub_pipeline
    (log,) = out.glob(f"metator_{command}_*.log")
    logged = log.read_text()
    assert "WARNING :: The Restriction_site column of" in logged and "contig(s) of the assembly" in logged
    message = f"{RS_ERROR}, as the Restriction_site column of {metacool} does not give every contig of the assembly its "
    assert f"ERROR :: {message}restriction sites.\n" in logged


def test_pipeline_refused_for_a_column_that_lacks_contigs_keeps_the_clusters(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A run from a map with -F (a rebinning in the output directory of an earlier
    run) that its columns refuse once the assembly is read leaves the clusters
    of the earlier run: the contig data is read before they are replaced (they
    were removed, then the run refused)."""
    monkeypatch.setattr("metator.network.contig_data_from_map", mtn_contig_data_from_map)
    contig_data = {name: {"length": length, "GC": 50.0, "coverage": "-", "RS": "-"}
                   for name, length in [("c1", 100), ("c2", 50), ("c3", 30)]}  # fmt: skip
    monkeypatch.setattr("metator.network.create_contig_data", lambda *args: contig_data)
    metacool = _map_with_columns(tmp_path, {"Restriction_site": [3, 2]})
    out = tmp_path / "out"
    # The directories that the run makes anew, which -F replaces.
    names = [os.path.join("fasta", "mags_with-mges"), os.path.join("fasta", "mags"), "cool"]
    clusters = [out / name / "cluster.fa" for name in names]
    for path in clusters:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    with pytest.raises(ValueError):
        run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-m", "RS", "-o", out, "-F",
                "-T", tmp_path / "tmp")  # fmt: skip
    assert all(path.exists() for path in clusters)
    (log,) = out.glob("metator_pipeline_*.log")
    assert "does not give every contig of the assembly its restriction sites." in log.read_text()


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


def _merged_map(tmp_path, normalization):
    """The map of two libraries (two copies of _contig_map with the network of
    their contigs, normalized by normalization) summed by metacooler merge,
    which keeps the contigs table of the first, and the metadata of both: the
    network that metator network stored in the first, whose Hit is no longer
    that of the contacts."""
    from metacooler.cli.main import cli as metacooler_cli

    library = _contig_map(tmp_path / "library_1.metacool", network=True, normalization=normalization)
    copy = shutil.copyfile(library, tmp_path / "library_2.metacool")
    merged = tmp_path / "merged.metacool"
    result = CliRunner().invoke(metacooler_cli, ["merge", str(library), str(copy), "-o", str(merged), "-q"])
    assert result.exit_code == 0, result.output
    return merged


def _collection(path, *entries, network=False, partition=False):
    """A .mmetacool collection at path, with a map of the contigs (_contig_map,
    with its network with network, and its partition with partition) as each
    entry of entries."""
    import metacooler

    maps = {
        name: _contig_map(path.with_name(f"{path.stem}_{name}.metacool"), network=network, partition=partition)
        for name in entries
    }
    metacooler.create_mmetacool(path, maps)
    return path


@pytest.mark.parametrize("where", ["elsewhere", "outdir", "--inplace"])
def test_pipeline_metacool_written_from_a_map(where, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """With the real writer, from the final clusters of the (stub) validation: the
    map of the clusters, from the input map, with the columns of its network
    (the GC content of the stub contig data), written to the output directory
    under the name of the input map, which is left as it was; but the input
    map is replaced when it is that map already, and with --inplace."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _contig_data(monkeypatch, CONTIGS)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map((out if where == "outdir" else tmp_path) / "metator.metacool")
    before = metacool.read_bytes()
    options = ["--inplace"] if where == "--inplace" else []
    run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp",
            *options)  # fmt: skip
    import cooler
    import metacooler

    outfile = out / "metator.metacool" if where != "--inplace" else metacool
    metacooler.validate(str(outfile), full=True)
    assert cooler.Cooler(str(outfile)).pixels()[:]["count"].sum() == 2
    mc = metacooler.MetaCooler(str(outfile))
    assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert list(mc.genomes.Completeness.iloc[1:]) == [0.9]
    assert dict(zip(mc.contigs.name, mc.contigs.GC_content, strict=True)) == {"c1": 40.5, "c2": 35.0}
    assert (metacool.read_bytes() == before) == (where == "elsewhere")
    expected = [] if where == "--inplace" else ["metator.metacool"]
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == expected
    assert not list(tmp_path.rglob("*.lock"))
    (log,) = out.glob("metator_pipeline_*.log")
    assert "is not written" not in log.read_text()


def test_pipeline_metacool_from_a_map_failure_keeps_it(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the clusters that cannot be written from the input map (here, a
    final cluster named as metacooler marks a contig in no cluster) stops the
    run, with an error that says why and that the input map is left as it was."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1", final_bin="NA")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool")
    before = metacool.read_bytes()
    message = (
        f"{metacool} is not written: A .metacool file cannot have clusters named as metacooler marks a contig in no "
        "cluster (NA, NaN, N/A, N/D, -1 or __unbinned__ in any case, and nd in any case other than ND, which marks a "
        "contig in no cluster): it would put their contigs in __unbinned__. Rename the 1 cluster(s): 'NA'. It is still "
        "the input map of the run, as it was. The other outputs of the run are kept."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    assert len(list((tmp_path / "tmp").rglob("cluster_summary.txt"))) == 1
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"ERROR :: {message}" in logged
    assert "Traceback" not in logged


@pytest.mark.parametrize(
    "error", [ValueError("pixels out of order"), RuntimeError("unexpected failure")], ids=["ValueError", "RuntimeError"]
)
def test_pipeline_metacool_from_a_map_metacooler_failure(error, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """metacooler fails while it writes the map of the clusters over the input map
    (the validation of its edit session, here): it removes its staged copy, and
    the run stops with an error that the map is not written, with the reason of
    one of the errors that metacooler is expected to raise (a ValueError here),
    or the traceback of another (a RuntimeError)."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool")
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise error

    monkeypatch.setattr("metacooler.schema.validate", validate)
    with pytest.raises(mtmc.MetacoolError, match=r" is not written: "):
        run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert metacool.read_bytes() == before
    # Staged by metacooler in a private directory next to the map, removed.
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
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


def _lock(metacool, pid, host=None, start="2026-11-18T10:00:00+0000"):
    """A lock file of the map metacool that names a process pid of host (this host
    by default), but not the machine and the PID namespace where it runs, which
    metacooler's locks name: metacooler cannot check it (as a lock of another
    host, or one that metaTOR wrote before its map writes went through
    metacooler), and never takes it over."""
    import socket

    holder = {"pid": pid, "host": host or socket.gethostname(), "start": start, "command": "metator pipeline",
              "uuid": "0123456789abcdef0123456789abcdef"}  # fmt: skip
    lock = Path(f"{metacool}.lock")
    lock.write_text(json.dumps(holder))
    return lock


def _rounds_tables(monkeypatch, members):
    """Have the stub validation also write the clustering rounds of the run (see
    metator.clusters): metator_00001_00000, an overlapping cluster of the
    contigs members, which is also the final cluster of c1."""
    import metator.clusters as mtcl

    stub = mtcli.mtv.recursive_decontamination

    def recursive_decontamination(*args, **kwargs):
        result = stub(*args, **kwargs)
        rounds = mtcl.Rounds()
        rounds.add(mtcl.OVERLAPPING, "metator_00001_00000", members, columns={"Completeness": 0.9, "Lineage": "x"})
        rounds.close_round(mtcl.OVERLAPPING)
        rounds.add(mtcl.FINAL, "metator_00001_00000", ["c1"], parent="metator_00001_00000", columns={"step": 0})
        rounds.write(args[7])
        return result

    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)


def test_pipeline_from_a_map_stores_its_rounds(stub_pipeline, tmp_path, monkeypatch, run_cli, caplog):
    """The map of the clusters holds the clustering rounds of the run in its
    /clusters group: its clusters, their parents and columns, and their
    contigs, those that the map lacks left out with a warning."""
    import metacooler

    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    _rounds_tables(monkeypatch, ["c1", "c2", "c9"])
    metacool = _contig_map(tmp_path / "input.metacool")
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    written = tmp_path / "out" / "input.metacool"
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
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
    """The rounds of an earlier run describe the clusters it replaced: a run that
    wrote none (here, the stub validation) leaves the map with none."""
    import metacooler
    import pandas as pd

    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    metacooler.write_clusters(
        str(metacool),
        ["earlier"],
        pd.DataFrame({"round": ["earlier"], "name": ["x"]}),
        pd.DataFrame({"round": ["earlier"], "cluster": ["x"], "contig": ["c2"]}),
    )
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert metacooler.MetaCooler(str(tmp_path / "out" / "input.metacool")).cluster_rounds == []
    # The input map, left as it was, keeps them.
    assert [entry["name"] for entry in metacooler.MetaCooler(str(metacool)).cluster_rounds] == ["earlier"]


def _map_of_an_earlier_run(path):
    """A map of the contigs c1 and c2 of _contig_map as an earlier run left it: c1
    in the cluster earlier_00001, c2 in __unbinned__, and the clustering rounds
    overlapping (made by partition), final (by validation) and mge (by mge)."""
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
@pytest.mark.parametrize("inplace", [False, True], ids=["the map of the output directory", "--inplace"])
def test_warn_of_what_a_run_replaces_in_a_map(command, replaced, inplace, tmp_path, logged):
    """What a run that updates a map in place replaces there, of an earlier
    run (or another binner): pipeline MAP its clusters and its rounds of
    the binning; partition MAP its rounds of the binning; validation MAP its
    clusters and its rounds of validation, as partition MAP, step by step,
    stores the round overlapping just before; mge MAP its round mge. A map
    of no earlier run gives no warning. The warning says how to leave the map
    as it is: without --inplace, or, for the map of the output directory of
    the run, with another one."""
    metacool = _map_of_an_earlier_run(tmp_path / "input.metacool")
    mtcli.warn_of_what_is_replaced(command, str(metacool), inplace=inplace)
    mtcli.warn_of_what_is_replaced(command, str(_contig_map(tmp_path / "fresh.metacool")), inplace=inplace)
    elsewhere = "leave out --inplace" if inplace else "give the run another output directory with -o"
    message = (
        f"{metacool} holds {replaced}, which this run replaces, in place: {elsewhere} to write the map of this run to "
        f"another file, and leave {metacool} as it is."
    )
    assert [record for record in logged() if record[0] == "WARNING"] == [("WARNING", message)]


@pytest.mark.parametrize("output", ["--inplace", "the output directory"])
def test_pipeline_from_a_map_warns_of_what_it_replaces(output, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """pipeline MAP --inplace warns before any work that it replaces the clusters
    and the rounds that an earlier run left in the map: it replaced them
    without a word, and without -F (the map keeps no trace of them). A run
    that writes its map to its output directory leaves them in MAP, and does
    not warn."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _map_of_an_earlier_run(tmp_path / "input.metacool")
    options = ["--inplace"] if output == "--inplace" else []
    run_cli("pipeline", *_map_run(metacool, tmp_path, *options))
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    replaced = (
        f"WARNING :: {metacool} holds 1 cluster(s) and the clustering round(s) overlapping, final, which this run "
        f"replaces, in place: leave out --inplace to write the map of this run to another file, and leave {metacool} "
        "as it is."
    )
    assert (replaced in log.read_text()) == (output == "--inplace")


@pytest.mark.parametrize("command", ["partition", "validation", "mge"])
def test_a_step_warns_of_what_it_replaces_in_place(command, tmp_path, logged):
    """partition MAP, validation MAP and mge MAP warn, as they take the lock
    of the map, of what they replace in it in place (see
    warn_of_what_is_replaced): with --inplace, or when MAP is the map of their
    output directory; when they write their map to another, they do not."""
    from contextlib import ExitStack

    metacool = _map_of_an_earlier_run(tmp_path / "input.metacool")
    # The warnings so far, after a run to another output directory, then with
    # --inplace, then in the directory of the map.
    for outdir, inplace, warnings in [(tmp_path / "out", False, 0), (tmp_path / "out", True, 1), (tmp_path, False, 2)]:
        with ExitStack() as held:
            mtcli.step_map(command, str(metacool), str(outdir), inplace, False, held)
        warned = [message for level, message in logged() if level == "WARNING"]
        assert len(warned) == warnings and all(message.startswith(f"{metacool} holds ") for message in warned)


def test_pipeline_to_an_existing_map_of_its_output_directory(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map of the output directory of the run, which an earlier run left
    there, is replaced only with -F: refused before any work, else the map of
    the clusters replaces it, the input map left as it was."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    output_map = tmp_path / "out" / "input.metacool"
    output_map.parent.mkdir()
    output_map.write_text("an older map")
    message = f"{output_map} exists: give -F to replace it with the map of this run, or another output directory with -o."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert "contig_contacts" not in stub_pipeline and output_map.read_text() == "an older map"
    run_cli("pipeline", *_map_run(metacool, tmp_path, "-F"))
    import metacooler

    metacooler.validate(str(output_map), full=True)
    assert metacool.read_bytes() == before


@pytest.mark.parametrize("how", ["the map of the output directory", "--inplace"])
def test_pipeline_updates_its_map_in_place(how, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map that is the map of the output directory of the run already (as the
    steps of a run, one by one, in one output directory, leave it), and any
    map with --inplace: updated in place, without -F."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    directory = tmp_path / "out" if how == "the map of the output directory" else tmp_path
    directory.mkdir(exist_ok=True)
    metacool = _contig_map(directory / "input.metacool")
    before = metacool.read_bytes()
    run_cli("pipeline", *_map_run(metacool, tmp_path, *(["--inplace"] if how == "--inplace" else [])))
    assert metacool.read_bytes() != before
    assert sorted(path.name for path in directory.iterdir() if "metacool" in path.name) == ["input.metacool"]


def test_pipeline_from_an_entry(stub_pipeline, tmp_path, run_cli):
    """An entry of a .mmetacool collection is not updated in place (--inplace),
    which would rewrite the collection: refused before any work. Without
    --inplace, the map of the run is a map of its own, named after the entry,
    in the output directory."""
    entry = f"{_collection(tmp_path / 'maps.mmetacool', 'rep1')}::entries/rep1"
    message = (
        f"{entry} is an entry of a .mmetacool collection, which a run does not update in place: leave out --inplace, "
        "and the run writes its map, a map of its own, in its output directory."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(entry, tmp_path, "--inplace"))
    assert "contig_contacts" not in stub_pipeline
    run_cli("pipeline", *_map_run(entry, tmp_path))
    source, _, outfile, _ = stub_pipeline["write_map"]
    assert (source, outfile) == (entry, str(tmp_path / "out" / "rep1.metacool"))


@pytest.mark.parametrize("how", ["--inplace", "the output directory"])
def test_pipeline_from_a_map_through_a_symbolic_link(how, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map given through a symbolic link: with --inplace, the map it points
    to is updated in place, as metacooler edits the target of a link, and the
    link keeps pointing at it; without, the map of the run is written to the
    output directory, under the name of the link, the map it points to left
    as it was."""
    import metacooler

    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    (tmp_path / "maps").mkdir()
    target = _contig_map(tmp_path / "maps" / "target.metacool")
    before = target.read_bytes()
    link = tmp_path / "link.metacool"
    link.symlink_to(target)
    run_cli("pipeline", *_map_run(link, tmp_path, *(["--inplace"] if how == "--inplace" else [])))
    written = target if how == "--inplace" else tmp_path / "out" / "link.metacool"
    assert link.is_symlink() and link.resolve() == target and (target.read_bytes() == before) == (how != "--inplace")
    assert list(metacooler.MetaCooler(str(written)).genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert sorted(path.name for path in (tmp_path / "maps").iterdir()) == ["target.metacool"]


@pytest.mark.parametrize("failure", ["partition", "Ctrl-C in the validation", "Ctrl-C in metacooler reassign"])
def test_pipeline_from_a_map_that_fails_leaves_it_as_it_was(failure, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A run that stops (an error mid-run, Ctrl-C, even while metacooler writes
    the new map) leaves the input map byte for byte as it was, with no file
    beside it: no staged map, no lock."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
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
    """A second run on a map that another run holds (here, another process that
    holds its lock) is refused before any work, with metacooler's reason (who
    holds it, and the files to delete if it is gone) and where else to write
    the map, and leaves the lock to its holder."""
    metacool = _contig_map(tmp_path / "input.metacool")
    with held_lock(metacool) as holder:
        message = f"{refusal(metacool, holder)} Or write this run's map elsewhere: leave out --inplace."
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("pipeline", *_map_run(metacool, tmp_path, "--inplace"))
        assert "contig_contacts" not in stub_pipeline
        assert json.loads(lock_file(metacool).read_text()) == holder
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("holder", ["another host", "no machine named", "unreadable", "not an object"])
def test_pipeline_from_a_map_refused_with_a_lock_it_cannot_check(holder, stub_pipeline, tmp_path, run_cli):
    """A lock that metacooler cannot tell stale is never taken over: that of a
    process of another host; one that does not name the machine and the PID
    namespace of its process (as metaTOR's own locks did, before its map writes
    went through metacooler), even of a process that is gone; or one that says
    nothing readable (not JSON, or JSON but not an object)."""
    import subprocess

    metacool = _contig_map(tmp_path / "input.metacool")
    lock = Path(f"{metacool}.lock")
    if holder in ("unreadable", "not an object"):
        lock.write_text("{not json" if holder == "unreadable" else "[1]")
        reason = refusal(metacool)
        assert "does not say which process holds" in reason
    else:
        gone = subprocess.Popen([sys.executable, "-c", "pass"])
        gone.wait()
        _lock(metacool, 1, host="elsewhere") if holder == "another host" else _lock(metacool, gone.pid)
        reason = refusal(metacool, json.loads(lock.read_text()))
        assert "a process this host cannot check" in reason
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(reason)} Or write this run's map elsewhere"):
        run_cli("pipeline", *_map_run(metacool, tmp_path, "--inplace"))
    assert lock.exists()


def test_pipeline_from_a_map_whose_lock_cannot_be_created(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map whose directory is not writable cannot be locked, nor updated in
    place: the run is refused before any work, with metacooler's reason and
    where else to write the map."""
    metacool = _contig_map(tmp_path / "input.metacool")
    lock = f"{metacool}.lock"
    opened = os.open

    def read_only(path, *args, **kwargs):
        if str(path) == lock:
            raise PermissionError(13, "Permission denied", path)
        return opened(path, *args, **kwargs)

    monkeypatch.setattr(os, "open", read_only)
    message = (
        f"{lock} cannot be created (Permission denied): an in-place write of {metacool} locks it there, so the "
        "directory of the map must exist and be writable. Write a copy of the map instead. Or write this run's map "
        "elsewhere: leave out --inplace."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path, "--inplace"))
    assert "contig_contacts" not in stub_pipeline
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("how", ["--inplace", "the map of the output directory"])
@pytest.mark.parametrize("command", ["pipeline", "partition", "validation"])
def test_a_read_only_map_is_refused_before_any_work(command, how, stub_pipeline, tmp_path, run_cli):
    """A map without write permission that the run would update in place (with
    --inplace, or as the map of its output directory) is refused before any
    work, with how to write a copy instead: metacooler refuses to replace it,
    which it would do at the end of the run. The map and its permission are
    left as they were, and so is the lock."""
    directory = tmp_path / "out" if how == "the map of the output directory" else tmp_path
    directory.mkdir(exist_ok=True)
    metacool = _contig_map(directory / "input.metacool")
    metacool.chmod(0o444)
    before = metacool.read_bytes()
    elsewhere = "leave out --inplace" if how == "--inplace" else "give the run another output directory with -o"
    message = (
        f"{metacool} is read-only (mode -r--r--r--): the run would replace it with the map of its clusters. Make it "
        f"writable (chmod u+w), or write the map of the run to another file: {elsewhere}."
    )
    options = ["--inplace"] if how == "--inplace" else []
    arguments = (
        _map_run(metacool, tmp_path, *options) if command == "pipeline" else _step_run(command, metacool, tmp_path, *options)
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *arguments)
    assert not stub_pipeline
    assert metacool.read_bytes() == before and metacool.stat().st_mode & 0o777 == 0o444
    assert not Path(f"{metacool}.lock").exists()


@pytest.mark.parametrize(
    "case",
    [
        "missing map",
        "an entry with --inplace",
        "a directory",
        "the collection of the input",
        "an output directory that holds ::",
    ],
)
def test_pipeline_map_refusals(case, stub_pipeline, tmp_path, run_cli):
    """What the map of the output directory of the run cannot be, an entry of a
    collection with --inplace, and an input map that is not there, are refused
    before any work (they failed at the end of the run, or on a bare error
    while locking)."""
    metacool = _contig_map(tmp_path / "input.metacool")
    options, given = [], metacool
    out = tmp_path / "out"
    if case == "missing map":
        given = tmp_path / "missing" / "map.metacool"
        message = f"Cannot read the input {given}: no such file."
    elif case == "an entry with --inplace":
        given = f"{_collection(tmp_path / 'maps.mmetacool', 'rep1')}::entries/rep1"
        options = ["--inplace"]
        message = (
            f"{given} is an entry of a .mmetacool collection, which a run does not update in place: leave out "
            "--inplace, and the run writes its map, a map of its own, in its output directory."
        )
    elif case == "a directory":
        (out / "input.metacool").mkdir(parents=True)
        options = ["-F"]
        message = (
            f"{out / 'input.metacool'} is a directory, where the run writes its map: move it away, or give the run "
            "another output directory with -o."
        )
    elif case == "the collection of the input":
        # A collection named as the map of the output directory of a run from
        # its entry.
        out.mkdir()
        collection = _collection(out / "rep1.metacool", "rep1")
        given = f"{collection}::entries/rep1"
        options = ["-F"]
        message = (
            f"{collection} is the collection of the input map {given}: a run does not replace a collection with its "
            "map. Give the run another output directory with -o."
        )
    else:
        options = ["-o", tmp_path / "a::b"]
        message = (
            f"{tmp_path / 'a::b' / 'input.metacool'}, the map of the run in its output directory, holds ::, which "
            "metacooler reads as file::group: give the run another output directory with -o."
        )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(given, tmp_path, *options))
    assert not stub_pipeline


def test_pipeline_in_a_new_output_directory(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """An output directory that is not there yet: it is made, with its parents,
    and the map of the run is written there, its lock removed."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    output = tmp_path / "new" / "deeper" / "input.metacool"
    run_cli("pipeline", *_map_run(metacool, tmp_path, "-o", output.parent))
    assert output.exists() and not Path(f"{output}.lock").exists()


@pytest.mark.parametrize("how", ["--inplace", "the output directory"])
@pytest.mark.parametrize("root", ["::/", "::"])
def test_pipeline_from_the_root_of_a_map(root, how, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map given as file::/ (or file::), the root of the file, is the file
    itself: with --inplace, it is updated in place as the file, next to which
    the run stages its map; without, the map of the output directory has the
    name of the file."""
    import metacooler

    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    run_cli("pipeline", *_map_run(f"{metacool}{root}", tmp_path, *(["--inplace"] if how == "--inplace" else [])))
    written = metacool if how == "--inplace" else tmp_path / "out" / "input.metacool"
    metacooler.validate(str(written), full=True)
    assert set(metacooler.MetaCooler(str(written)).genomes.name) == {"metator_00001_00000", "__unbinned__"}
    assert sorted(path.name for path in tmp_path.iterdir() if path.name.startswith(("input", ".input"))) == ["input.metacool"]
    assert (metacool.read_bytes() == before) == (how != "--inplace")


def test_a_map_is_written_to_a_file_only(tmp_path):
    """A map is not written to a URI (file::group), which metacooler reads as a
    group inside the file: an existing file would be truncated."""
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    clusters = mtmc.FinalClusters({}, str(tmp_path / "cluster_summary.txt"), None, {})
    outfile = f"{metacool}::/"
    message = f"^{re.escape(outfile)} is not the path of a file: a map is written to a file of its own, not to a group"
    with pytest.raises(mtmc.MetacoolError, match=message):
        mtmc_write_map(str(metacool), outfile, clusters=clusters, force=True)
    assert metacool.read_bytes() == before


@pytest.mark.parametrize("command", ["network", "pipeline", "partition", "validation", "mge", "dump"])
def test_an_output_directory_that_cannot_be_made(command, stub_pipeline, tmp_path, run_cli):
    """The output directory of a run, where it writes its map (and its log), cannot
    be made (a file is in its way): refused before any work, the reason
    logged."""
    blocking = tmp_path / "afile"
    blocking.write_text("x")
    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=True)
    message = f"The output directory {blocking} cannot be made (File exists): give the run another one with -o."
    arguments = {
        "network": [metacool, "-a", _assembly(tmp_path)],
        "pipeline": [metacool, "-a", _assembly(tmp_path)],
        "partition": [metacool],
        "validation": [metacool, "-a", _assembly(tmp_path)],
        "mge": [metacool, "-a", _assembly(tmp_path), "-m", _assembly(tmp_path)],
        "dump": [metacool],
    }[command]
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *arguments, "-o", blocking, *([] if command == "dump" else ["-T", tmp_path / "tmp"]))
    assert _no_step(stub_pipeline) and blocking.read_text() == "x"


@pytest.mark.parametrize(
    ("command", "output"),
    [("pipeline", False), ("pipeline", True), ("partition", False)],
    ids=["pipeline --inplace", "pipeline to its output directory", "partition --inplace"],
)
def test_a_map_changed_during_the_run_is_refused_at_its_end(command, output, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map changed while a run works from it (here, rewritten by a program
    that takes no lock of metacooler, once the partition is done) is refused
    when the run writes its map: the run read the network from the map as it
    was at its start, and would write its clusters into another map. The
    map of the run is not written, with the reason, in place or in its output
    directory, and the map is left as the other program left it."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool", network=command == "partition")
    partition = mtcli.mtp.partition

    def partition_then_rewrite(*args, **kwargs):
        written = partition(*args, **kwargs)
        rewritten = tmp_path / "rewritten.tmp"
        shutil.copyfile(metacool, rewritten)
        os.replace(rewritten, metacool)
        return written

    monkeypatch.setattr(mtcli.mtp, "partition", partition_then_rewrite)
    options = [] if output else ["--inplace"]
    arguments = (
        _map_run(metacool, tmp_path, *options) if command == "pipeline" else _step_run(command, metacool, tmp_path, *options)
    )
    reason = (
        f"{metacool} was replaced or written since it was read for this run (by another run, or a program that "
        "takes no lock of metacooler): the clusters of the run come from the map as it was. Run it again on the map "
        "as it is now."
    )
    with pytest.raises(mtmc.MetacoolError, match=re.escape(reason)):
        run_cli(command, *arguments)
    rewritten = metacool.read_bytes()
    assert not (tmp_path / "out" / "input.metacool").exists()
    import metacooler

    assert list(metacooler.MetaCooler(str(metacool)).genomes.name) == ["__unbinned__"]
    assert metacool.read_bytes() == rewritten and not lock_file(metacool).exists()


def _output_map_run(command, metacool, tmp_path, *options):
    """The arguments of a pipeline or partition run from the map metacool, with
    options (--inplace, -F...)."""
    if command == "pipeline":
        return _map_run(metacool, tmp_path, *options)
    return _step_run(command, metacool, tmp_path, *options)


@pytest.mark.parametrize("command", ["pipeline", "partition"])
@pytest.mark.parametrize("kind", ["symbolic link", "hard link"])
def test_the_map_of_the_output_directory_a_link_to_the_input_map(kind, command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map of the output directory, a link to the input map. A symbolic link
    there is a map of its own: refused without -F, as any map there, and with
    -F replaced by the map of the run, the input map left as it was. A hard
    link is the input map under another name: refused before any work, with
    -F too."""
    import metacooler

    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    (tmp_path / "maps").mkdir()
    metacool = _contig_map(tmp_path / "maps" / "input.metacool", network=command == "partition")
    before = metacool.read_bytes()
    out = tmp_path / "out" / "input.metacool"
    out.parent.mkdir()
    if kind == "symbolic link":
        out.symlink_to(metacool)
    else:
        os.link(metacool, out)
    if kind == "hard link":
        message = (
            f"{out} is a hard link of the input map {metacool}: the same file under another name, where the run writes "
            "its map. Give the run another output directory with -o, or --inplace to update the map in place."
        )
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli(command, *_output_map_run(command, metacool, tmp_path, "-F"))
        assert "partition" not in stub_pipeline
        assert metacool.read_bytes() == before and os.path.samefile(metacool, out)
        return
    message = (
        f"{out} exists (a symbolic link to the input map, which -F replaces, not the map): give -F to replace it with "
        "the map of this run, or another output directory with -o."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *_output_map_run(command, metacool, tmp_path))
    run_cli(command, *_output_map_run(command, metacool, tmp_path, "-F"))
    assert metacool.read_bytes() == before and not out.is_symlink()
    mc = metacooler.MetaCooler(str(out))
    if command == "pipeline":
        assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    else:
        assert list(mc.clusters.name) == ["core_00001", "metator_00001_00000"]
    assert sorted(path.name for path in (tmp_path / "maps").iterdir()) == ["input.metacool"]
    assert not [path.name for path in out.parent.iterdir() if path.name.startswith(".") or path.suffix == ".lock"]


@pytest.mark.parametrize("case", ["a map that another process holds", "a link to nothing"])
def test_the_map_of_the_output_directory_a_symbolic_link_locks_its_own_name(case, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map of the output directory, a symbolic link: the run locks the
    link's own name, as it replaces the link (with -F), never what it points
    to."""
    from contextlib import ExitStack

    import metacooler

    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "out" / "input.metacool"
    out.parent.mkdir()
    with ExitStack() as held:
        if case == "a link to nothing":
            out.symlink_to(tmp_path / "gone" / "old.metacool")
        else:
            (tmp_path / "elsewhere").mkdir()
            target = _contig_map(tmp_path / "elsewhere" / "other.metacool")
            out.symlink_to(target)
            held.enter_context(held_lock(target))
            before = target.read_bytes()
        run_cli("pipeline", *_map_run(metacool, tmp_path, "-F"))
        if case == "a map that another process holds":
            assert target.read_bytes() == before
    assert not out.is_symlink()
    assert list(metacooler.MetaCooler(str(out)).genomes.name) == ["__unbinned__", "metator_00001_00000"]


def test_hold_lock_follows_the_links_of_the_map_in_place_only(monkeypatch):
    """hold_lock locks the map that a run writes: in place, the file that a
    symbolic link points to, as metacooler edits it (metacooler.lock follows
    links); otherwise the map of the output directory, its own name
    (follow_symlinks=False), as a link there is replaced, not its target."""
    import contextlib

    import metacooler

    seen = []

    @contextlib.contextmanager
    def lock(uri, *, command=None, follow_symlinks=True):
        seen.append((uri, follow_symlinks))
        yield uri

    monkeypatch.setattr(metacooler, "lock", lock)
    with contextlib.ExitStack() as held:
        mtcli.hold_lock(held, "map.metacool", True)
        mtcli.hold_lock(held, "copy.metacool", False)
    assert seen == [("map.metacool", True), ("copy.metacool", False)]


@pytest.mark.parametrize("command", ["pipeline", "partition"])
def test_a_writer_before_the_lock_is_seen_by_the_run(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """Another writer that takes metacooler's lock (here metacooler annotate
    --inplace) and ends after the run has checked its map, before the run
    takes the map's lock: the run reads its map once it holds the lock, so
    it works from the map that writer left, and writes its own, which keeps
    that writer's column. It read the map before taking the lock, did all
    its work, and was refused at its end."""
    import metacooler

    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool", network=command == "partition")
    table = tmp_path / "extra.txt"
    table.write_text("Name\tnote\nc1\tx\nc2\ty\n")
    run_map_path = mtcli.run_map_path

    def another_writer_then(*args, **kwargs):
        metacooler.annotate(str(metacool), contigs_table=str(table), contigs_table_format="Name", inplace=True)
        return run_map_path(*args, **kwargs)

    monkeypatch.setattr(mtcli, "run_map_path", another_writer_then)
    run_cli(command, *_output_map_run(command, metacool, tmp_path, "--inplace"))
    mc = metacooler.MetaCooler(str(metacool))
    assert "note" in mc.contigs.columns
    if command == "pipeline":
        assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    else:
        assert list(mc.clusters.name) == ["core_00001", "metator_00001_00000"]


def _replaced_once_opened_by(monkeypatch, metacool, reader):
    """Have the map metacool replaced by another file (a copy of itself), as a run
    that updates it in place replaces it, once the function reader opens it,
    before it reads it."""
    opened = mtmc.open_map

    def open_map(path):
        mc = opened(path)
        if sys._getframe(1).f_code.co_name == reader:
            copy = Path(f"{metacool}.copy")
            shutil.copyfile(metacool, copy)
            os.replace(copy, metacool)
        return mc

    monkeypatch.setattr(mtmc, "open_map", open_map)


REPLACED_WHILE_READ = (
    "{map} was replaced or written while it was read, by a run that updated it: what was read of it would mix two "
    "maps. Read it again once that run has ended."
)


@pytest.mark.parametrize(
    ("reader", "command", "options"),
    [
        ("mge_binning", "mge", []),
        ("map_columns", "pipeline", ["--normalization", "RS"]),
        ("warn_of_what_is_replaced", "partition", ["--inplace"]),
    ],
)
def test_a_map_replaced_as_the_run_checks_it(reader, command, options, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map replaced as the run reads it to check it (what mge reads of it,
    its columns against the normalization, what it holds that the run
    replaces) stops the run with metaTOR's error, logged as an ERROR: it
    stopped on metacooler's MapChangedError, raised as it is, with nothing
    in the log."""
    metacool = _contig_map(tmp_path / "input.metacool")
    _replaced_once_opened_by(monkeypatch, metacool, reader)
    message = REPLACED_WHILE_READ.format(map=metacool)
    options = [tmp_path / option if option.endswith(".metacool") else option for option in options]
    if command == "mge":
        mges = tmp_path / "mges.txt"
        mges.write_text("c2\n")
        arguments = [metacool, "-a", _assembly(tmp_path), "-m", mges, "-o", tmp_path / "out", *options]
    else:
        arguments = _output_map_run(command, metacool, tmp_path, *options)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *arguments)
    assert "partition" not in stub_pipeline and "recursive_decontamination" not in stub_pipeline
    # mge writes its log in the current directory.
    (log,) = (Path.cwd() if command == "mge" else tmp_path / "out").glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}" in log.read_text()


def test_a_map_replaced_as_validation_reads_its_partition(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map replaced as a validation run reads the partition that partition MAP
    stored in it (see metator.metacool.read_partition) stops the run with
    metaTOR's error, logged as an ERROR, before any work."""
    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=True)
    read_partition = mtmc.read_partition

    def replaced_then_read(*args, **kwargs):
        copy = Path(f"{metacool}.copy")
        shutil.copyfile(metacool, copy)
        os.replace(copy, metacool)
        return read_partition(*args, **kwargs)

    monkeypatch.setattr(mtmc, "read_partition", replaced_then_read)
    message = REPLACED_WHILE_READ.format(map=metacool)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("validation", *_step_run("validation", metacool, tmp_path))
    assert "recursive_decontamination" not in stub_pipeline
    (log,) = (tmp_path / "out").glob("metator_validation_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("reader", ["map_columns", "contig_data_from_map"])
def test_a_map_replaced_as_its_columns_are_read(reader, tmp_path, monkeypatch):
    """The columns of a map, read for the network (map_columns,
    contig_data_from_map), from a map replaced as they are read: metaTOR's
    error, from metacooler's MapChangedError."""
    import metacooler

    import metator.network as mtn

    metacool = _contig_map(tmp_path / "input.metacool")
    _replaced_once_opened_by(monkeypatch, metacool, reader)
    contig_data = {
        name: {"length": length, "GC": 0.5, "coverage": "-", "RS": "-"} for name, length in {"c1": 100, "c2": 50}.items()
    }
    read = {"map_columns": lambda: mtn.map_columns(str(metacool)), "contig_data_from_map": lambda: mtn_contig_data_from_map(str(metacool), contig_data)}  # fmt: skip
    message = REPLACED_WHILE_READ.format(map=metacool)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$") as raised:
        read[reader]()
    assert isinstance(raised.value.__cause__, metacooler.MapChangedError)


def test_the_map_of_the_output_directory_held_by_another_run_is_refused_without_the_hint(stub_pipeline, tmp_path, run_cli):
    """A run whose map of the output directory another process holds (another
    run to the same output directory) is refused before any work with
    metacooler's reason only: the hint to write the map elsewhere, which a run
    in place gets, is no help to a run that does."""
    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "out" / "input.metacool"
    out.parent.mkdir()
    with held_lock(out) as holder:
        message = refusal(out, holder)
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert "contig_contacts" not in stub_pipeline
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}\n" in log.read_text()


def test_pipeline_from_a_map_takes_a_stale_lock(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The lock of a run that is gone (killed with SIGKILL, or SIGTERM: neither
    metaTOR nor metacooler handles signals), and the copy that it staged next
    to the map, are taken over and removed by the next run on this host, each
    with a warning in its log; the run goes on, and removes its own lock at its
    end."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    holder = killed_lock(metacool)
    staged = tmp_path / f".input.metacool.{holder['uuid']}.1.d"
    staged.mkdir()
    (staged / "input.metacool").write_text("the copy that the killed run staged")
    run_cli("pipeline", *_map_run(metacool, tmp_path, "--inplace"))
    assert not lock_file(metacool).exists() and not staged.exists()
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"WARNING :: Took over {lock_file(metacool)}: the process that held it ({described(holder)}) is gone." in logged
    assert f"WARNING :: Removed {staged.name}, left beside {metacool} by a writer that no longer holds its lock." in logged


def test_network_from_a_metacool(monkeypatch, tmp_path, run_cli):
    """metator network MAP reads the network from the map, with the contig data of
    the assembly, and stores it in the map of the output directory: for an
    entry of a collection, which is left as it was, a map of its own, named
    after the entry."""
    import metacooler

    calls = {}
    contig_data = {name: {"length": length, "GC": 50.0, "coverage": "-", "RS": 2}
                   for name, length in [("c1", 100), ("c2", 50)]}  # fmt: skip

    def create_contig_data(*args):
        calls["create_contig_data"] = args
        return contig_data

    def contig_data_from_map(*args):
        calls["contig_data_from_map"] = args
        return []

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    monkeypatch.setattr("metator.network.contig_data_from_map", contig_data_from_map)
    out = tmp_path / "out"
    collection = _collection(tmp_path / "maps.mmetacool", "lib1")
    before = collection.read_bytes()
    # The URI as typed, relative to the working directory.
    monkeypatch.chdir(tmp_path)
    run_cli("network", "maps.mmetacool::entries/lib1", "-a", _assembly(tmp_path), "-e", "DpnII", "-n", "RS", "-o", out,
            "-T", tmp_path / "tmp")  # fmt: skip
    assert calls["create_contig_data"] == (str(tmp_path / "assembly.fa"), None, "DpnII")
    # The restriction sites and the coverage that the map gives where the
    # options give none.
    assert calls["contig_data_from_map"] == ("maps.mmetacool::entries/lib1", contig_data, "DpnII", None)
    assert collection.read_bytes() == before
    written = out / "lib1.metacool"
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    # The hits of c1 and c2 (a contact between them, one within c1), their
    # restriction sites, and the RS normalization: 10 * sqrt(sites / hits).
    assert list(mc.contigs.name) == ["c1", "c2"] and list(mc.contigs.Hit) == [3, 1]
    assert list(mc.contigs.Restriction_site) == [2, 2] and list(mc.contigs.GC_content) == [50.0, 50.0]
    assert list(metacooler.read_weights(str(written), "metator_RS")) == pytest.approx([10 * (2 / 3) ** 0.5, 10 * 2**0.5])
    entry = mc.metadata["metator"]["network"]
    assert entry["normalization"] == "RS" and entry["version"] == __version__
    assert entry["command"].startswith("metator network maps.mmetacool::entries/lib1 -a")
    assert entry["parameters"]["inplace"] is False
    assert mc.cluster_rounds == [] and list(mc.genomes.name) == ["__unbinned__"]
    (log,) = out.glob("metator_network_*.log")
    assert sorted(path.name for path in out.iterdir()) == sorted([log.name, "lib1.metacool"])
    assert not list(tmp_path.glob("*.lock")) and not list(out.glob("*.lock")) and not (tmp_path / "tmp").exists()


def test_network_on_an_entry_with_inplace(stub_pipeline, tmp_path, run_cli):
    """metator network writes its map: an entry of a collection, which a run does
    not update in place, is refused with --inplace, before any work."""
    collection = _collection(tmp_path / "maps.mmetacool", "lib1")
    before = collection.read_bytes()
    entry = f"{collection}::entries/lib1"
    message = (
        f"{entry} is an entry of a .mmetacool collection, which a run does not update in place: leave out --inplace, "
        "and the run writes its map, a map of its own, in its output directory."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("network", entry, "-a", _assembly(tmp_path), "-o", tmp_path / "out", "-T", tmp_path / "tmp", "--inplace")
    assert collection.read_bytes() == before and not stub_pipeline


def test_network_updates_its_map_in_place(stub_pipeline, tmp_path, run_cli):
    """metator network MAP --inplace stores the network in MAP, in place, under
    its lock, and replaces the network of an earlier run (its columns and its
    weights): run again with another normalization, the map holds that one
    only."""
    import metacooler

    metacool = _contig_map(tmp_path / "input.metacool")
    metacool.chmod(0o640)
    for normalization in ["empirical_hit", "length", "None"]:
        run_cli("network", metacool, "-a", _assembly(tmp_path), "-n", normalization, "-o", tmp_path / "out",
                "-T", tmp_path / "tmp", "--inplace")  # fmt: skip
        metacooler.validate(str(metacool), full=True)
        mc = metacooler.MetaCooler(str(metacool))
        weights = [column for column in mc.bins().columns if column.startswith("metator_")]
        assert weights == ([] if normalization == "None" else [f"metator_{normalization}"])
        assert mc.metadata["metator"]["network"]["normalization"] == normalization
        assert list(mc.contigs.Hit) == [3, 1]
    assert metacool.stat().st_mode & 0o777 == 0o640
    assert not list(tmp_path.glob("*.lock")) and not [path for path in tmp_path.iterdir() if path.name.startswith(".")]
    # The map is the network that partition reads: the counts, without
    # weights.
    import metator.network as mtn

    assert mtn.map_network(str(metacool)).to_dict("list") == {"contig1": [1], "contig2": [2], "weight": [1]}


def test_network_keeps_the_partition_of_the_map(stub_pipeline, tmp_path, run_cli):
    """metator network MAP keeps the clusters of the map, its partition among
    them, with what partition recorded of it."""
    import metacooler

    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=True)
    mc = metacooler.MetaCooler(str(metacool))
    rounds, members, entry = mc.cluster_rounds, mc.cluster_members(), mc.metadata["metator"]["partition"]
    warning = (
        f"WARNING :: {metacool} holds a partition of another network than the one stored in it now, which the run "
        f"kept: metator validation refuses it until metator partition {metacool} runs again, which partitions this "
        "network."
    )
    for normalization, stale in [("empirical_hit", False), ("None", True)]:
        run_cli("network", metacool, "-a", _assembly(tmp_path), "-n", normalization, "-o", tmp_path / normalization,
                "-T", tmp_path / "tmp", "--inplace")  # fmt: skip
        mc = metacooler.MetaCooler(str(metacool))
        assert mc.cluster_rounds == rounds and mc.cluster_members().equals(members)
        assert mc.metadata["metator"]["partition"] == entry
        (log,) = (tmp_path / normalization).glob("metator_network_*.log")
        assert (warning in log.read_text()) == stale
        digest = mtcli.mtn.stored_network(str(metacool)).digest
        assert (digest != entry[mtmc.NETWORK_DIGEST]) == stale
    with pytest.raises(mtmc.MetacoolError, match=re.escape(f"{metacool} holds a partition of another network")):
        run_cli("validation", *_step_run("validation", metacool, tmp_path, "--inplace"))
    run_cli("partition", *_step_run("partition", metacool, tmp_path, "--inplace"))
    assert metacooler.MetaCooler(str(metacool)).metadata["metator"]["partition"][mtmc.NETWORK_DIGEST] == digest
    run_cli("validation", *_step_run("validation", metacool, tmp_path, "--inplace"))
    assert [entry["name"] for entry in metacooler.MetaCooler(str(metacool)).cluster_rounds][:2] == ["core", "overlapping"]


# What the inputs of a network or pipeline run cannot be, refused before any
# work, the reason in its log file.
NOT_A_MAP = (
    "{name} is not a .metacool contact map: give the map of the contigs (made by hicstuff pipeline -M metacool or "
    "metacooler create pairs), or the pairs files or the BAM files of the run, which the run makes it from."
)
READS = (
    "{reads} holds reads (FASTQ): metaTOR aligns no reads. Give the BAM files or the pairs files of the run, or the map "
    "of its contigs, which hicstuff makes from the reads (hicstuff pipeline -M metacool)."
)
RUN_INPUT_ERRORS = {
    "reads as the first argument": (["{reads}"], READS),
    "a map that is not there": (["{missing}"], "Cannot read the input {missing}: no such file."),
    "a URI of a directory": (["{directory}::/"], "Cannot read the contact map {directory}::/: {directory} is not a file."),
    "a file named as a URI": (["{reads}::/"], NOT_A_MAP.replace("{name}", "{reads}::/")),
    "a map whose name holds ::": (
        ["{colon_map}"],
        "Cannot read the contact map {colon_map}: metacooler reads a path that holds :: as file::group. Rename the map "
        "without ::.",
    ),
}


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize("case", list(RUN_INPUT_ERRORS))
def test_run_input_refused(command, case, stub_pipeline, tmp_path, run_cli):
    """A run takes a .metacool map, or pairs or BAM files, as its arguments:
    reads, a file that is not there, or a URI whose file is no map nor a
    file, are refused before any work, the reason in the log file (it failed
    once the run read the map). A map whose name holds :: is refused with the
    reason, as metacooler cannot read it (a pairs file whose name holds ::
    is that file: see test_a_pairs_file_whose_name_holds_colons)."""
    names = {
        "map": _contig_map(tmp_path / "input.metacool"),
        "reads": tmp_path / "R1.fq",
        "missing": tmp_path / "missing.metacool",
        "directory": tmp_path / "maps",
        "colon_map": tmp_path / "run::1.metacool",
    }
    names["reads"].write_text("@r1\nACGT\n+\nIIII\n")
    names["directory"].mkdir()
    shutil.copyfile(names["map"], names["colon_map"])
    args, message = RUN_INPUT_ERRORS[case]
    args = [arg.format(**names) for arg in args]
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli(command, *args, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert not stub_pipeline
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message.format(**names)}\n" in log.read_text()


def test_a_pairs_file_whose_name_holds_colons(tmp_path, run_cli):
    """A pairs file whose name holds :: is that file, which the run makes its
    map from (it was refused as no map): metacooler reads it as a path."""
    import metacooler

    pairs = shutil.copyfile(TESTS_DATA / "alignment.pairs", tmp_path / "lib::1.pairs")
    out = tmp_path / "out"
    run_cli("network", pairs, "-a", TESTS_DATA / "assembly.fa", "-o", out, "-T", tmp_path / "tmp")
    assert metacooler.MetaCooler(str(out / "assembly.metacool")).info["sum"] == 41563


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


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_a_map_is_the_one_input(command, stub_pipeline, tmp_path, run_cli):
    """A map given with another input (a pairs file after it, which the run
    would make a map from) is refused before any work, the reason in the log
    file: it was a usage error, one argument too many."""
    metacool = _contig_map(tmp_path / "input.metacool")
    pairs = shutil.copyfile(TESTS_DATA / "alignment.pairs", tmp_path / "lib.pairs")
    out = tmp_path / "out"
    with pytest.raises(mtmc.MetacoolError):
        run_cli(command, metacool, pairs, "-a", _assembly(tmp_path), "-o", out, "-T", tmp_path / "tmp")
    assert not stub_pipeline
    (log,) = out.glob(f"metator_{command}_*.log")
    assert (
        f"ERROR :: A .metacool map is the one input of a run, given alone, and this run has 2: {metacool}, {pairs}.\n"
        in log.read_text()
    )
    assert sorted(path.name for path in out.iterdir()) == [log.name]


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_a_run_needs_its_map(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A network or pipeline run needs its inputs (the map, or pairs or BAM
    files): without them, a usage error, and the command does not start."""
    monkeypatch.chdir(tmp_path)
    result = run_cli(command, "-a", "a.fa", "-o", "out", exit_code=2)
    assert "Missing argument 'INPUT...'" in result.output
    assert not stub_pipeline and not Path("out").exists()


@pytest.mark.parametrize("command", ["network", "pipeline", "partition", "validation"])
def test_a_directory_is_no_map(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A directory given as the map is a usage error, before the command
    starts (it failed once the run read the map, and pipeline -F had
    removed the clusters of the earlier run)."""
    monkeypatch.chdir(tmp_path)
    Path("maps").mkdir()
    bins = Path("out", "fasta", "mags_with-mges", "bin.fa")
    bins.parent.mkdir(parents=True)
    bins.write_text(">c1\nACGT\n")
    force = [] if command == "network" else ["-F"]
    assembly = [] if command == "partition" else ["-a", _assembly(tmp_path)]
    result = run_cli(command, "maps", *assembly, "-o", "out", *force, exit_code=2)
    assert "File 'maps' is a directory." in " ".join(result.output.replace("│", " ").split())
    assert not stub_pipeline and bins.exists() and not list(Path("out").glob("*.log"))


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_help_says_what_the_columns_of_the_map_stand_for(command):
    """The help of the options that the columns of a map stand in for
    (--depth, --enzyme, the normalization) names those columns of the
    map."""
    from metator.cli import cli

    helps = {param.name: param.help for param in cli.commands[command].params if getattr(param, "help", None)}
    assert [name for name, text in helps.items() if "--start" in text] == []
    said = {
        "depth": "The Shotgun_coverage column of the map",
        "enzyme": "The Restriction_site column of the map",
        "normalization": "The Restriction_site and Shotgun_coverage columns of the map",
    }
    for name, words in said.items():
        assert words in helps[name], name


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize(
    "option",
    [["-E", "100"], ["--edge", "0"], ["-1", "R1.fq"], ["-S", "pair"], ["-b", "bwa"]],
    ids=["-E", "--edge", "-1", "-S", "-b"],
)
def test_options_of_reads_and_pairs_are_gone(command, option, stub_pipeline, tmp_path, run_cli):
    """The options of a run from reads, alignments or pairs, and --edge,
    which metaTOR 2.0 removed, are unknown options: click's usage error,
    with no hint. The command does not start. (-q, --min-quality is back:
    the mapping quality of the alignments of BAM files, which the arguments
    of a run may be.)"""
    result = run_cli(command, "map.metacool", *option, "-a", "a.fa", "-o", tmp_path / "out", exit_code=2)
    assert f"No such option '{option[0]}'." in " ".join(result.output.replace("│", " ").split())
    assert not stub_pipeline and not (tmp_path / "out").exists()


# partition MAP and validation MAP: the steps run one by one on a map, with
# the steps replaced by stubs (see tests/test_binning_synthetic.py for the
# map that they give).

# The directories of clusters of each step, replaced with -F: partition writes
# none.
# The directories of clusters of an earlier run: those of an earlier version
# (recursive_clusters), refused, and those of the run (fasta/mags_with-mges),
# which -F replaces.
MAGS_WITH_MGES = os.path.join("fasta", "mags_with-mges")
STEP_DIRS = {"partition": [], "validation": ["recursive_clusters", MAGS_WITH_MGES]}


def _step_run(command, metacool, tmp_path, *options):
    """The arguments of a partition or validation run on the map metacool (for
    None, a map of the contigs c1 and c2 with their network, which metator
    network stores, and their partition, which metator partition stores)."""
    assembly = ["-a", _assembly(tmp_path)] if command == "validation" else []
    out = ["-o", tmp_path / "out", "-T", tmp_path / "tmp"]
    if metacool is None:
        metacool = tmp_path / "network.metacool"
        if not metacool.exists():
            _contig_map(metacool, network=True, partition=True)
    return [metacool, *assembly, *out, *options]


def _earlier_run(command, tmp_path):
    """The directories of clusters that an earlier run of command left in the
    output directory, each with a FASTA file, which the test checks."""
    files = [tmp_path / "out" / name / "earlier.fa" for name in STEP_DIRS[command]]
    for path in files:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    return files


# The results of a run of an earlier version in the output directory of a
# validation or pipeline run (see metator.cli.EARLIER_OUTPUTS): the tables,
# final clusters and plots of 2.0's development builds, and the outputs of 1.x,
# files and directories, the plots of 1.x that 2.0 no longer draws among them
# (the run draws its own in plots/: see test_the_plots_of_an_earlier_run).
EARLIER_FILES = [
    "bin_summary.txt",
    "binning.txt",
    "cluster_members.txt",
    "cluster_summary.txt",
    "clustering_matrix.npz",
    "clusters.txt",
    "contig_data_final.txt",
    "contig_data_network.txt",
    "contig_data_partition.txt",
    "network.txt",
    "overlapping_micomplete_results.txt",
    "plot/MAGs_GC_distribution.png",
    "plot/bins_distribution.png",
    "plot/bins_size_distribution.png",
    "plot/MAGs_HiC_cov_distribution.png",
    "plot/MAGs_SG_cov_distribution.png",
    "plot/network_heatmap.png",
    "recursive_micomplete_step_1.txt",
]
EARLIER_DIRECTORIES = [
    "final_bin",
    "final_bin_unscaffold",
    "final_clusters",
    "overlapping_bin",
    "recursive_bin",
    "recursive_clusters",
    "scaffold_bin",
]
EARLIER_OUTPUTS = sorted([*EARLIER_FILES, *EARLIER_DIRECTORIES])


@pytest.mark.parametrize("command", ["validation", "pipeline"])
@pytest.mark.parametrize("force", [False, True], ids=["refused", "-F"])
def test_the_outputs_of_an_earlier_version(command, force, stub_pipeline, tmp_path, run_cli):
    """An earlier run's results in the output directory, or a user's own files
    of those names, are refused with or without -F, before any work, and
    nothing is removed, a link's target included."""
    metacool = _contig_map(tmp_path / "network.metacool", network=True, partition=True)
    out = tmp_path / "out"
    for directory in EARLIER_DIRECTORIES:
        (out / directory).mkdir(parents=True)
        (out / directory / "earlier.fa").write_text(">c1\nACGT\n")
    for name in EARLIER_FILES:
        (out / name).parent.mkdir(parents=True, exist_ok=True)
        (out / name).write_text("an earlier run's\n")
    (out / "binning_other.txt").write_text("another tool's\n")
    (tmp_path / "kept.txt").write_text("the target of a link\n")
    (out / "cluster_summary.txt").unlink()
    (out / "cluster_summary.txt").symlink_to(tmp_path / "kept.txt")
    before = sorted(path.relative_to(out) for path in out.rglob("*"))
    if command == "validation":
        args = _step_run("validation", metacool, tmp_path, *(["-F"] if force else []))
    else:
        args = _map_run(metacool, tmp_path, *(["-F"] if force else []))
    with pytest.raises(ValueError):
        run_cli(command, *args)
    assert sorted(path.relative_to(out) for path in out.rglob("*") if not path.name.endswith(".log")) == before
    assert (tmp_path / "kept.txt").read_text() == "the target of a link\n"
    assert _no_step(stub_pipeline)
    (log,) = out.glob(f"metator_{command}_*.log")
    shown = ", ".join(repr(name) for name in EARLIER_OUTPUTS[:5])
    assert (
        f"ERROR :: {out} holds results of an earlier run, which this run does not write and does not replace: "
        f"{shown}, ... (+{len(EARLIER_OUTPUTS) - 5} more). Move them away, or give the run another output directory "
        "with -o."
    ) in log.read_text()


@pytest.mark.parametrize("command", ["validation", "pipeline"])
@pytest.mark.parametrize("force", [False, True], ids=["refused", "-F"])
def test_the_plots_of_an_earlier_run(command, force, stub_pipeline, tmp_path, run_cli):
    """The plots of an earlier run in plots/ are refused without -F, before any
    work, and removed with -F once the run is checked, for the validation to
    draw them again: a symbolic link, not its target. The other files of plots/
    are left as they are, with or without -F: a file of the user's, and the map
    of the run itself."""
    out = tmp_path / "out"
    plot = out / "plots"
    plot.mkdir(parents=True)
    metacool = _contig_map(plot / "network.metacool", network=True, partition=True)
    before = metacool.read_bytes()
    (plot / "notes.txt").write_text("the user's\n")
    (plot / "bins_distribution.pdf").write_text("an earlier run's\n")
    (tmp_path / "kept.pdf").write_text("the target of a link\n")
    (plot / "bins_size_distribution.pdf").symlink_to(tmp_path / "kept.pdf")
    options = ["-F"] if force else []
    args = (
        _step_run("validation", metacool, tmp_path, *options)
        if command == "validation"
        else _map_run(metacool, tmp_path, *options)
    )
    if not force:
        with pytest.raises(ValueError):
            run_cli(command, *args)
        assert all(os.path.lexists(plot / name) for name in ["bins_distribution.pdf", "bins_size_distribution.pdf"])
        assert _no_step(stub_pipeline) and not (out / "fasta").exists()
        (log,) = out.glob(f"metator_{command}_*.log")
        assert (
            f"ERROR :: {plot} holds the plots of an earlier run ('bins_distribution.pdf', 'bins_size_distribution.pdf'), "
            "which this run draws again: give -F to replace them, or give the run another output directory with -o."
        ) in log.read_text()
    else:
        run_cli(command, *args)
        # Removed, for the validation to draw them (the stub draws none).
        assert not any(os.path.lexists(plot / name) for name in ["bins_distribution.pdf", "bins_size_distribution.pdf"])
        assert stub_pipeline["recursive_decontamination_kwargs"]["plot_outdir"] == str(out)
    assert (plot / "notes.txt").read_text() == "the user's\n" and metacool.read_bytes() == before
    assert (tmp_path / "kept.pdf").read_text() == "the target of a link\n"


@pytest.mark.parametrize("force", [False, True], ids=["refused", "-F"])
def test_the_tables_of_an_earlier_run(force, stub_pipeline, tmp_path, run_cli):
    """The tables of an earlier run in tables/, which a pipeline run writes again
    from the map of its clusters (metator dump), are refused without -F, before
    any work, and removed with -F once the run is checked: they would describe
    another map, should the run stop before it writes its own. The other files
    of tables/ are left as they are, with or without -F."""
    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "out"
    tables = out / "tables"
    tables.mkdir(parents=True)
    (tables / "binning.txt").write_text("an earlier run's\n")
    (tables / "notes.txt").write_text("the user's\n")
    if not force:
        with pytest.raises(ValueError):
            run_cli("pipeline", *_map_run(metacool, tmp_path))
        assert (tables / "binning.txt").read_text() == "an earlier run's\n"
        assert _no_step(stub_pipeline) and not (out / "fasta").exists()
        (log,) = out.glob("metator_pipeline_*.log")
        assert (
            f"ERROR :: {tables} holds the tables of an earlier run ('binning.txt'), which this run writes again: give -F "
            "to replace them, or give the run another output directory with -o."
        ) in log.read_text()
    else:
        run_cli("pipeline", *_map_run(metacool, tmp_path, "-F"))
        # Removed, for the run to write them again (its stub outputs write none).
        assert not os.path.lexists(tables / "binning.txt") and stub_pipeline["run_export"]
    assert (tables / "notes.txt").read_text() == "the user's\n"


@pytest.mark.parametrize("force", [False, True], ids=["refused", "-F"])
def test_the_report_of_an_earlier_run(force, stub_pipeline, tmp_path, run_cli):
    """The report of an earlier run in the output directory, which a pipeline run
    writes again in its step export, is refused without -F, before any work,
    and removed with -F once the run is checked, as its tables are (see
    test_the_tables_of_an_earlier_run): a symbolic link too, not its target.
    A directory where a file of the report goes is refused with or without
    -F."""
    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "out"
    out.mkdir()
    (out / mtr.REPORT_TXT).write_text("an earlier run's\n")
    (tmp_path / "kept.html").write_text("the target of a link\n")
    (out / mtr.REPORT_HTML).symlink_to(tmp_path / "kept.html")
    if not force:
        with pytest.raises(ValueError):
            run_cli("pipeline", *_map_run(metacool, tmp_path))
        assert (out / mtr.REPORT_TXT).read_text() == "an earlier run's\n" and (out / mtr.REPORT_HTML).is_symlink()
        assert _no_step(stub_pipeline) and not (out / "fasta").exists()
        (log,) = out.glob("metator_pipeline_*.log")
        assert (
            f"ERROR :: {out} holds the report of an earlier run ('report.txt', 'report.html'), which this run writes "
            "again: give -F to replace it, or give the run another output directory with -o."
        ) in log.read_text()
    else:
        run_cli("pipeline", *_map_run(metacool, tmp_path, "-F"))
        # Removed, for the run to write it again (its stub export writes none).
        assert not any(os.path.lexists(out / name) for name in mtr.REPORT_FILES) and stub_pipeline["run_export"]
    assert (tmp_path / "kept.html").read_text() == "the target of a link\n"
    # A directory where a file of the report goes, with or without -F.
    (out / mtr.REPORT_MQC).mkdir()
    with pytest.raises(ValueError):
        run_cli("pipeline", *_map_run(metacool, tmp_path, "-F"))
    (log,) = sorted(out.glob("metator_pipeline_*.log"))[-1:]
    assert (
        f"ERROR :: {out / mtr.REPORT_MQC} is a directory, where the run writes its report: move it away, or give the run "
        "another output directory with -o."
    ) in log.read_text()


@pytest.mark.parametrize("command", ["validation", "pipeline"])
@pytest.mark.parametrize("what", ["plots, a file", "a plot, a directory"])
def test_where_the_run_draws_no_plot(command, what, stub_pipeline, tmp_path, run_cli):
    """A plots/ that is no directory, or a directory where a plot goes, is refused
    with or without -F, before any work, and left as it is: the run would draw
    no plot there."""
    metacool = _contig_map(tmp_path / "network.metacool", network=True, partition=True)
    out = tmp_path / "out"
    out.mkdir()
    if what == "plots, a file":
        where = out / "plots"
        where.write_text("a file\n")
        message = f"{where} is not a directory, where the run draws its plots"
    else:
        where = out / "plots" / "bins_size_distribution.pdf"
        where.mkdir(parents=True)
        message = f"{where} is a directory, where the run draws a plot"
    args = _step_run("validation", metacool, tmp_path, "-F") if command == "validation" else _map_run(metacool, tmp_path, "-F")
    with pytest.raises(ValueError):
        run_cli(command, *args)
    assert where.is_file() if what == "plots, a file" else where.is_dir()
    assert _no_step(stub_pipeline) and not (out / "fasta").exists()
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}: move it away, or give the run another output directory with -o." in log.read_text()


@pytest.mark.parametrize("command", ["validation", "pipeline"])
def test_the_map_of_the_run_is_no_earlier_output(command, stub_pipeline, tmp_path, run_cli):
    """The map of the run, in a directory named as a 1.x output, is no earlier result."""
    out = tmp_path / "out"
    (out / "final_bin").mkdir(parents=True)
    metacool = _contig_map(out / "final_bin" / "network.metacool", network=True, partition=True)
    if command == "validation":
        args = _step_run("validation", metacool, tmp_path)
    else:
        args = _map_run(metacool, tmp_path)
    run_cli(command, *args)
    assert (out / "fasta" / "mags_with-mges").is_dir()


@pytest.mark.parametrize("command", ["validation", "pipeline"])
@pytest.mark.parametrize(
    "where",
    [
        "in fasta/mags_with-mges",
        "recursive_clusters",
        "recursive_bin",
        "through a link of -o",
        "through a link of fasta/mags_with-mges",
    ],
)
def test_tmpdir_in_the_clusters_refused(command, where, stub_pipeline, tmp_path, run_cli):
    """A --tmpdir in fasta/mags_with-mges, which the run makes anew, or in a
    directory of an earlier run's results, as written or through links, is
    refused as a usage error before anything is created or removed."""
    metacool = _contig_map(tmp_path / "network.metacool", network=True, partition=True)
    earlier = _earlier_run("validation", tmp_path)
    out = tmp_path / "out"
    if where == "through a link of -o":
        (tmp_path / "alias").symlink_to(out)
    if where == "through a link of fasta/mags_with-mges":
        shutil.move(out / MAGS_WITH_MGES, tmp_path / "elsewhere")
        (out / MAGS_WITH_MGES).symlink_to(tmp_path / "elsewhere")
        earlier = [tmp_path / "elsewhere" / "earlier.fa", earlier[0]]
    outdir = tmp_path / "alias" if where == "through a link of -o" else out
    tmpdir = {
        "in fasta/mags_with-mges": out / MAGS_WITH_MGES / "tmp",
        "recursive_clusters": out / "recursive_clusters",
        "recursive_bin": out / "recursive_bin" / "tmp",
        "through a link of -o": out / MAGS_WITH_MGES / "tmp",
        "through a link of fasta/mags_with-mges": out / MAGS_WITH_MGES / "tmp",
    }[where]
    clusters = outdir / (where if where in ("recursive_clusters", "recursive_bin") else MAGS_WITH_MGES)
    if command == "validation":
        args = _step_run("validation", metacool, tmp_path, "-o", outdir, "-T", tmpdir, "-F")
    else:
        args = [metacool, "-a", _assembly(tmp_path), "-o", outdir, "-T", tmpdir, "-F"]
    before = sorted(path.relative_to(tmp_path) for path in tmp_path.rglob("*"))
    # The map made, nothing read yet.
    stub_pipeline.clear()
    result = run_cli(command, *args, exit_code=2)
    where = (
        "where an earlier version wrote its results, which the run refuses"
        if where in ("recursive_clusters", "recursive_bin")
        else "which the run makes anew (-F replaces an earlier one)"
    )
    message = f"Invalid value for '-T' / '--tmpdir': {tmpdir} is in {clusters}, {where}: give a temporary directory outside it."
    assert message in " ".join(result.output.replace("│", " ").split())
    assert all(path.exists() for path in earlier)
    assert sorted(path.relative_to(tmp_path) for path in tmp_path.rglob("*")) == before
    assert not stub_pipeline


@pytest.mark.parametrize(
    ("directory", "mges"),
    [(os.path.join("fasta", "mags"), False), ("cool", False), (os.path.join("fasta", "mges"), True)],
    ids=["fasta/mags", "cool", "fasta/mges with --mges"],
)
def test_tmpdir_in_an_output_of_the_pipeline_refused(directory, mges, stub_pipeline, tmp_path, run_cli):
    """A --tmpdir in a directory that pipeline makes anew for its outputs, which
    -F replaces (fasta/mags, cool, and fasta/mges with --mges), is refused as a
    usage error before anything is created or removed; fasta/mges is no
    output of a run without --mges."""
    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "out"
    tmpdir = out / directory / "tmp"
    options = ["--mges", tmp_path / "mges.txt"] if mges else []
    (tmp_path / "mges.txt").write_text("c2\n")
    result = run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-o", out, "-T", tmpdir, "-F", *options, exit_code=2)
    message = (
        f"Invalid value for '-T' / '--tmpdir': {tmpdir} is in {out / directory}, which the run makes anew (-F replaces "
        "an earlier one): give a temporary directory outside it."
    )
    assert message in " ".join(result.output.replace("│", " ").split())
    assert not out.exists() and not stub_pipeline


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_help_of_a_step_on_a_map(command, run_cli):
    """The help of partition and validation shows MAP, their first argument,
    which metator network stores the network in, and --inplace."""
    output = " ".join(run_cli(command, "--help").output.split())
    assert f"Usage: metator {command} [OPTIONS] MAP" in output and "metator network" in output
    assert f"metator {command} MAP" in output and "--inplace" in output


# What partition MAP and validation MAP refuse before any work, with -F (an
# earlier run's directories of clusters are kept), the reason in the log file
# of the run: the map that run_input refuses, and what run_map_path refuses.
STEP_MAP_ERRORS = {
    "a map that is not there": ("{missing}", [], "Cannot read the contact map {missing}: no such file."),
    "a URI of a directory": ("{directory}::/", [], "Cannot read the contact map {directory}::/: {directory} is not a file."),
    "no map": (
        "{pairs}",
        [],
        "{pairs} is not a .metacool contact map: give the map of the contigs of the run as the first argument of {command}.",
    ),
    "a map without a network": (
        "{contig_map}",
        [],
        "{contig_map} holds no network of its contigs: run metator network on it first, which stores the network in the map.",
    ),
    # A network that the run reads, then refuses, before -F replaces anything.
    "a network with text in a column": (
        "{text_map}",
        [],
        "The contigs table of {text_map} holds a value that is not a number in a column that metator network stores "
        "(could not convert string to float: 'high'): run metator network on the map again.",
    ),
    "a network of other contacts": (
        "{merged_map}",
        [],
        "The Hit column of {merged_map} is not the hits of its contacts for 2 contig(s): 'c1', 'c2'. Its contacts changed "
        "since metator network stored the network (a merge of maps keeps the contigs table of the first): run metator "
        "network on the map again.",
    ),
    "a collection": ("{collection}", [], "Cannot read the contact map {collection}: "),
    "a map without a partition": (
        "{map}",
        [],
        "{map} holds no partition of its network (the clustering round(s) core and overlapping): run metator partition "
        "on it first, which stores its core and overlapping clusters in the map.",
    ),
    # metator network stored another network over the partition, which it
    # keeps (see test_network_keeps_the_partition_of_the_map).
    "a partition of another network": (
        "{stale_map}",
        [],
        "{stale_map} holds a partition of another network than the one stored in it now: metator network stored "
        "another network in the map since metator partition ran on it. Run metator partition {stale_map} again, which "
        "partitions the network stored in the map.",
    ),
    "a partition that records no network": (
        "{unrecorded_map}",
        [],
        "{unrecorded_map} does not record the network that its partition was made on, which metator partition "
        "records. Run metator partition {unrecorded_map} again, which partitions the network stored in the map.",
    ),
    # The map of the output directory, a link to the assembly.
    "an input as the map of the output directory": (
        "{map}",
        [],
        "{out_map}, where the run writes its map, is {assembly}, an input of the run: give the run another output "
        "directory with -o.",
    ),
    "a map whose name holds ::": (
        "{colon_map}",
        [],
        "Cannot read the contact map {colon_map}: metacooler reads a path that holds :: as file::group. Rename the map "
        "without ::.",
    ),
    "an entry with --inplace": (
        "{entry}",
        ["--inplace"],
        "{entry} is an entry of a .mmetacool collection, which a run does not update in place: leave out --inplace, "
        "and the run writes its map, a map of its own, in its output directory.",
    ),
    "a directory as the map of the output directory": (
        "{map}",
        [],
        "{out_map} is a directory, where the run writes its map: move it away, or give the run another output directory with -o.",
    ),
    "a lock held by another run": (
        "{map}",
        ["--inplace"],
        # metacooler's reason (see tests.locks.refusal).
        "{refusal} Or write this run's map elsewhere: leave out --inplace.",
    ),
}


# The refusals of a validation alone: it reads the partition that partition
# stores, and the assembly, which partition does not read.
VALIDATION_ONLY = {
    "a map without a partition",
    "a partition of another network",
    "a partition that records no network",
    "an input as the map of the output directory",
}


STEP_CASES = [
    (command, case)
    for case in STEP_MAP_ERRORS
    for command in ["partition", "validation"]
    if command == "validation" or case not in VALIDATION_ONLY
]


@pytest.mark.parametrize("command, case", STEP_CASES)
def test_step_on_a_map_refused(command, case, stub_pipeline, tmp_path, run_cli):
    """partition MAP and validation MAP check the map as network MAP and
    pipeline MAP do, then lock the map that they write, before any work:
    refused, the reason logged as an error in the log file of the run and
    raised as a MetacoolError, before the step and before -F replaces the
    directories of clusters of an earlier run. A lock that a running run
    holds (here, another process) is left to it. A collection is refused
    with metacooler's reason, which names its entries. validation MAP
    refuses a map that partition MAP has not stored its partition in, or
    whose partition is not that of the network stored in it."""
    from contextlib import ExitStack

    import metacooler

    names = {
        "command": command,
        "map": _contig_map(tmp_path / "input.metacool", network=True),
        "contig_map": _contig_map(tmp_path / "contigs.metacool"),
        "missing": tmp_path / "missing.metacool",
        "directory": tmp_path / "maps",
        "pairs": _pairs(tmp_path),
        "collection": _collection(tmp_path / "maps.mmetacool", "lib1"),
        "colon_map": tmp_path / "run::1.metacool",
        # The copy of the assembly that a validation run reads (see _step_run).
        "assembly": tmp_path / "assembly.fa",
        # The map that a run from {map} writes in its output directory.
        "out_map": tmp_path / "out" / "input.metacool",
    }
    if case == "an input as the map of the output directory":
        names["out_map"].parent.mkdir()
        names["out_map"].symlink_to(_assembly(tmp_path))
    if case == "a directory as the map of the output directory":
        names["out_map"].mkdir(parents=True)
    if case == "a network with text in a column":
        names["text_map"] = _contig_map(tmp_path / "text.metacool", network=True)
        with h5py.File(names["text_map"], "r+") as group:
            del group["contigs/GC_content"]
            group["contigs/GC_content"] = np.array([b"high", b"low"])
    if case == "a network of other contacts":
        names["merged_map"] = _merged_map(tmp_path, "None")
    if case == "a partition of another network":
        names["stale_map"] = _contig_map(tmp_path / "stale.metacool", network=True, partition=True)
        _store_contig_network(names["stale_map"], "None")
    if case == "a partition that records no network":
        names["unrecorded_map"] = _contig_map(tmp_path / "unrecorded.metacool", network=True, partition=True)
        mc = metacooler.MetaCooler(str(names["unrecorded_map"]))
        metadata = dict(mc.metadata)
        del metadata["metator"]["partition"][mtmc.NETWORK_DIGEST]
        with mc.open("r+") as group:
            group.attrs["metadata"] = json.dumps(metadata)
    shutil.copyfile(names["map"], names["colon_map"])
    names["entry"] = f"{names['collection']}::entries/lib1"
    names["directory"].mkdir()
    metacool, options, message = STEP_MAP_ERRORS[case]
    with ExitStack() as held:
        if case == "a lock held by another run":
            holder = held.enter_context(held_lock(names["map"]))
            names["refusal"] = refusal(names["map"], holder)
        earlier = _earlier_run(command, tmp_path)
        metacool = metacool.format(**names)
        message = message.format(**names)
        pattern = f"^{re.escape(message)}" + (".*entries/lib1" if case == "a collection" else "$")
        with pytest.raises(mtmc.MetacoolError, match=pattern):
            run_cli(command, *_step_run(command, metacool, tmp_path, *[option.format(**names) for option in options], "-F"))
        if case == "a lock held by another run":
            assert json.loads(lock_file(names["map"]).read_text()) == holder
    assert _no_step(stub_pipeline) and all(path.exists() for path in earlier)
    (log,) = (tmp_path / "out").glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}" in log.read_text()
    # No map is written.
    if "the map of the output directory" not in case:
        assert not names["out_map"].exists()


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_step_on_a_map_refuses_an_existing_map_of_its_output_directory(command, stub_pipeline, tmp_path, run_cli):
    """The map of the output directory is replaced only with -F, as in pipeline:
    refused before any work, the file left as it was."""
    metacool = _contig_map(tmp_path / "input.metacool", network=True)
    output_map = tmp_path / "out" / "input.metacool"
    output_map.parent.mkdir()
    output_map.write_text("an older map")
    message = f"{output_map} exists: give -F to replace it with the map of this run, or another output directory with -o."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *_step_run(command, metacool, tmp_path))
    assert not stub_pipeline and output_map.read_text() == "an older map"
    (log,) = (tmp_path / "out").glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("output", ["--inplace", "the output directory", "-F the output directory"])
def test_partition_on_a_map_stores_its_partition(output, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """partition MAP stores its partition in the clustering rounds of the map:
    every core cluster (core_00001...) and every overlapping cluster
    (metator_00001_00000...), whatever its size, as the rounds core and
    overlapping, with their contigs and without the columns of miComplete, in
    place of the rounds the map held, and what the run was, as the entry
    partition of the metator metadata of the map; its genomes are left as they
    were, and it writes no table and no FASTA file."""
    import metacooler
    import pandas as pd

    metacool = _contig_map(tmp_path / "input.metacool", network=True)
    metacooler.write_clusters(
        str(metacool),
        ["earlier"],
        pd.DataFrame({"round": ["earlier"], "name": ["x"]}),
        pd.DataFrame({"round": ["earlier"], "cluster": ["x"], "contig": ["c2"]}),
    )
    metacool.chmod(0o640)
    before = metacool.read_bytes()
    out = tmp_path / "out"
    output_map = out / "input.metacool"
    options = ["--inplace"] if output == "--inplace" else []
    if output.startswith("-F"):
        out.mkdir()
        output_map.write_text("an older map")
        options.append("-F")
    locks = []
    write_map = mtmc.write_map

    def locked(*args, **kwargs):
        locks.append(sorted(str(path.relative_to(tmp_path)) for path in [*tmp_path.glob("*.lock"), *out.glob("*.lock")]))
        return write_map(*args, **kwargs)

    monkeypatch.setattr("metator.metacool.write_map", locked)
    run_cli("partition", *_step_run("partition", metacool, tmp_path, "-i", "4", *options))
    written = metacool if output == "--inplace" else output_map
    assert locks == [[f"{written.relative_to(tmp_path)}.lock"]]
    assert not list(tmp_path.glob("*.lock")) and not list(out.glob("*.lock"))
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    entries = [{"name": name, "tool": "metaTOR", "version": __version__, "step": "partition"} for name in ["core", "overlapping"]]
    assert mc.cluster_rounds == entries
    assert list(mc.clusters.name) == ["core_00001", "metator_00001_00000"] and list(mc.clusters.n_contigs) == [2, 2]
    # genome, which metacooler computes: no cluster is a genome of the map yet.
    assert list(mc.clusters.columns) == ["round", "name", "parent_round", "parent", "n_contigs", "length", "genome"]
    assert mc.clusters.genome.isna().all()
    assert sorted(mc.cluster_members().contig) == ["c1", "c1", "c2", "c2"]
    assert list(mc.genomes.name) == ["__unbinned__"]
    # What the run was: the network that it partitioned, its command line.
    partition = mc.metadata["metator"]["partition"]
    assert partition["normalization"] == "empirical_hit" and partition["version"] == __version__
    assert partition["command"].startswith(f"metator partition {metacool}") and partition["parameters"]["iterations"] == 4
    assert mc.metadata["metator"]["network"]["normalization"] == "empirical_hit"
    # The stub partition got the network that metator network stored.
    algorithm, network, iterations, overlap, res_param, threads, prefix = stub_pipeline["partition"]
    assert (algorithm, iterations, overlap, res_param, threads, prefix) == ("louvain", 4, 0.8, 1.0, 1, "metator")
    assert network.to_dict("list") == {"contig1": [1], "contig2": [2], "weight": [pytest.approx(1 / 3**0.5, rel=1e-15)]}
    assert (metacool.read_bytes() == before) == (output != "--inplace")
    if output == "--inplace":
        assert metacool.stat().st_mode & 0o777 == 0o640
    assert not [path.name for path in [*tmp_path.iterdir(), *out.iterdir()] if path.name.startswith(".")]
    # No table, no FASTA file: the partition is in the map.
    maps = [] if output == "--inplace" else ["input.metacool"]
    assert sorted(path.name for path in out.iterdir() if not path.name.endswith(".log")) == maps
    (log,) = out.glob("metator_partition_*.log")
    logged = log.read_text()
    assert re.search(rf"INFO :: Contact map with the partition written after 0h 0m [0-9]+s: {re.escape(str(written))}\n", logged)
    assert ("Removing" in logged) == output.startswith("-F")
    # Which map the run writes, from its start.
    if output == "--inplace":
        assert f"INFO :: The run updates {metacool} in place (--inplace).\n" in logged
    else:
        assert f"INFO :: The run writes its map to {written}, and leaves {metacool} as it is (--inplace updates it).\n" in logged


def test_partition_on_an_entry_of_a_collection(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """From an entry of a .mmetacool collection, the map of the output directory
    is a map of its own, named after the entry, with the partition, and the
    collection is left as it was."""
    import metacooler

    collection = _collection(tmp_path / "maps.mmetacool", "lib1", "lib2", network=True)
    before = collection.read_bytes()
    output_map = tmp_path / "out" / "lib1.metacool"
    run_cli("partition", *_step_run("partition", f"{collection}::entries/lib1", tmp_path))
    metacooler.validate(str(output_map), full=True)
    mc = metacooler.MetaCooler(str(output_map))
    assert list(mc.clusters.name) == ["core_00001", "metator_00001_00000"] and list(mc.genomes.name) == ["__unbinned__"]
    assert collection.read_bytes() == before


def test_partition_on_a_map_through_a_symbolic_link(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map given through a symbolic link, with --inplace: the map it points to
    stores the partition, as metacooler edits the target of a link, and the
    link keeps pointing at it, as pipeline updates it."""
    import metacooler

    target = _contig_map(tmp_path / "target.metacool", network=True)
    link = tmp_path / "link.metacool"
    link.symlink_to(target)
    run_cli("partition", *_step_run("partition", link, tmp_path, "--inplace"))
    assert link.is_symlink() and link.resolve() == target
    assert list(metacooler.MetaCooler(str(target)).clusters.name) == ["core_00001", "metator_00001_00000"]


@pytest.mark.parametrize("error", [ValueError("not valid"), KeyboardInterrupt()], ids=["ValueError", "Ctrl-C"])
def test_partition_on_a_map_that_cannot_be_written(error, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map with the partition that cannot be written (here, its validation
    fails) stops the run with an error that says why and that the map is left
    as it was, as pipeline says it; Ctrl-C stops it too."""
    (tmp_path / "maps").mkdir()
    metacool = _contig_map(tmp_path / "maps" / "input.metacool", network=True)
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise error

    # The full validation of the staged map, before metacooler renames it.
    monkeypatch.setattr("metacooler.schema.validate", validate)
    args = _step_run("partition", metacool, tmp_path, "--inplace")
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


@pytest.mark.parametrize("output", ["--inplace", "the output directory"])
def test_validation_on_a_map_writes_the_map_of_the_clusters(output, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """validation MAP writes the map of the clusters as pipeline MAP writes it at
    its end, from the final clusters of the run, with its clustering rounds: to
    the output directory, the map left as it was, or in place (--inplace),
    keeping the mode of the map."""
    import metacooler

    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    _rounds_tables(monkeypatch, ["c1", "c2"])
    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=True)
    metacool.chmod(0o640)
    before = metacool.read_bytes()
    output_map = tmp_path / "out" / "input.metacool"
    options = ["--inplace"] if output == "--inplace" else []
    run_cli("validation", *_step_run("validation", metacool, tmp_path, *options))
    written = metacool if output == "--inplace" else output_map
    metacooler.validate(str(written), full=True)
    mc = metacooler.MetaCooler(str(written))
    assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert [entry["name"] for entry in mc.cluster_rounds] == ["overlapping", "final"]
    metator = mc.metadata["metator"]
    assert list(metator) == ["network", "partition", "validation"]
    validation = metator["validation"]
    assert validation["command"].startswith(f"metator validation {metacool} -a")
    assert validation["parameters"]["map"] == str(metacool) and "start" not in validation["parameters"]
    assert validation["parameters"]["inplace"] is (output == "--inplace")
    assert validation["summary_columns"] == ["Completeness"] and validation["normalization"] == "empirical_hit"
    assert (metacool.read_bytes() == before) == (output != "--inplace")
    if output == "--inplace":
        assert metacool.stat().st_mode & 0o777 == 0o640
    directories = [tmp_path, tmp_path / "out"]
    assert not [path.name for directory in directories for path in directory.iterdir() if path.name.startswith(".")]
    assert not [path.name for directory in directories for path in directory.glob("*.lock")]


def test_validation_on_a_map_that_cannot_be_written(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the clusters that cannot be written (here, its validation fails)
    stops validation MAP as it stops pipeline MAP: the map is left as it was."""
    monkeypatch.setattr("metator.metacool.write_map", mtmc_write_map)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=True)
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise ValueError("not valid")

    # The full validation of the staged map, before metacooler renames it.
    monkeypatch.setattr("metacooler.schema.validate", validate)
    with pytest.raises(mtmc.MetacoolError, match=r"is not written: .* It is still the input map of the run, as it was\."):
        run_cli("validation", *_step_run("validation", metacool, tmp_path, "--inplace"))
    assert metacool.read_bytes() == before


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_a_step_on_a_map_holds_its_lock(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """partition MAP and validation MAP, which write a map (the map of their
    output directory), hold its lock (metacooler's) through their step, and
    set no signal handler: a run killed by SIGTERM dies of it, and leaves its
    lock to the next writer (see hold_lock). A step that fails releases the
    lock."""
    import signal

    seen = {}
    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=command == "validation")

    def step(*args, **kwargs):
        lock = lock_file(tmp_path / "out" / "input.metacool")
        seen["held"] = lock.exists() and json.loads(lock.read_text())["pid"] == os.getpid()
        seen["handler"] = signal.getsignal(signal.SIGTERM)
        raise RuntimeError("stopped")

    target = {"partition": "metator.partition.partition", "validation": "metator.validation.recursive_decontamination"}
    monkeypatch.setattr(target[command], step)
    with pytest.raises(RuntimeError, match="stopped"):
        run_cli(command, *_step_run(command, metacool, tmp_path))
    assert seen == {"held": True, "handler": signal.SIG_DFL}
    assert not list(tmp_path.glob("*.lock")) and not list((tmp_path / "out").glob("*.lock"))


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_a_step_refuses_a_bad_prefix_before_replacing_the_clusters(command, stub_pipeline, tmp_path, run_cli):
    """An invalid --prefix is refused before -F replaces the directories of
    clusters of an earlier run, as pipeline refuses it: it emptied them
    first, then was refused."""
    earlier = _earlier_run(command, tmp_path)
    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=command == "validation")
    with pytest.raises(ValueError, match="Invalid --prefix"):
        run_cli(command, *_step_run(command, metacool, tmp_path, "-F", "-p", "bad.prefix"))
    assert all(path.exists() for path in earlier) and _no_step(stub_pipeline)


def test_pipeline_replaces_an_earlier_map_named_after_its_assembly(stub_pipeline, tmp_path, run_cli):
    """pipeline -F replaces the map of an earlier run in its output directory,
    named after the assembly (assembly.fa.metacool, for -a assembly.fa, from a
    map of that name): the assembly itself is the input that the map must not
    be."""
    metacool = _contig_map(tmp_path / "assembly.fa.metacool")
    (tmp_path / "out").mkdir()
    output_map = _contig_map(tmp_path / "out" / "assembly.fa.metacool")
    run_cli("pipeline", *_map_run(metacool, tmp_path, "-F"))
    assert stub_pipeline["write_map"][2] == str(output_map)


def test_pipeline_refuses_an_input_as_its_map(stub_pipeline, tmp_path, run_cli):
    """A map of the output directory that is the assembly of the run (an
    assembly named as the map) is refused before any work: with -F, the map of
    the clusters would replace it at the end of the run."""
    metacool = _contig_map(tmp_path / "input.metacool")
    assembly = tmp_path / "out" / "input.metacool"
    assembly.parent.mkdir()
    shutil.copyfile(_assembly(tmp_path), assembly)
    before = assembly.read_bytes()
    message = (
        f"{assembly}, where the run writes its map, is {assembly}, an input of the run: give the run another output "
        "directory with -o."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path, "-a", assembly, "-F"))
    assert assembly.read_bytes() == before and not stub_pipeline


@pytest.mark.parametrize("dangling", [False, True], ids=["link", "dangling link"])
def test_force_replaces_a_link_of_clusters_not_its_target(dangling, stub_pipeline, tmp_path, run_cli):
    """A directory of clusters of an earlier run that is a symbolic link, to a
    directory elsewhere or to nothing, is refused without -F, and replaced with
    -F: the link is removed, its target kept. -F failed with OSError (rmtree on
    a link), after the other directory was replaced."""
    target = tmp_path / "elsewhere"
    if not dangling:
        target.mkdir()
        (target / "earlier.fa").write_text(">c1\nACGT\n")
    link = tmp_path / "out" / MAGS_WITH_MGES
    link.parent.mkdir(parents=True)
    link.symlink_to(target)
    with pytest.raises(ValueError):
        run_cli("validation", *_step_run("validation", None, tmp_path))
    assert os.path.lexists(link) and _no_step(stub_pipeline)
    run_cli("validation", *_step_run("validation", None, tmp_path, "-F"))
    assert link.is_dir() and not link.is_symlink() and stub_pipeline
    assert dangling or (target / "earlier.fa").exists()


# The commands of the tutorial of metaTOR 1.x (4253551,
# docs/example/metator_tutorial.md) that its command line accepted: its -S
# pairs and -S network were no values of --start. Its pipeline ... -v, which
# the command line read as --scaffold, is gone with --scaffold.
COMMANDS_OF_4253551 = [
    "pipeline -1 hic_reads_for.fastq.gz -2 hic_reads_rev.fastq.gz -a assembly.fa -o metator_folder --threads 16",
    "pipeline -1 hic_reads_for.fastq.gz -2 hic_reads_rev.fastq.gz -a assembly_bowtie2_index -S fastq",
    "pipeline -1 for.bam -2 rev.bam -a assembly.fa -S bam",
    "network -1 for.bam -2 rev.bam -a assembly.fa",
    "partition -c contig_data_network.txt -n network.txt -a assembly.fa",
    "validation -c contig_data_partition.txt -n network.txt -a assembly.fa -f overlapping_bin",
    "network -1 for.bam -2 rev.bam -a assembly.fa -n None",
]


@pytest.mark.parametrize("command", COMMANDS_OF_4253551)
def test_the_commands_of_metator_1(command):
    """The command lines of metaTOR 1.x that its tutorial gives (parsed here,
    not run) are usage errors, with no hint: a network or pipeline run from
    reads, alignments or pairs (-1, -2, --start), -1 an unknown option, as
    metaTOR 2.0 starts from a map, which hicstuff makes; partition and
    validation from the tables of the network and of the partition (-c, -n
    and validation's -f), -c an unknown option of both, as they read the
    network that metator network stores in the map, and validation the
    partition that metator partition stores there."""
    name, *args = command.split()
    ctx = click.Context(mtcli.cli, info_name="metator")
    option = {"network": "-1", "pipeline": "-1", "partition": "-c", "validation": "-c"}[name]
    with pytest.raises(click.NoSuchOption, match=rf"^No such option '{option}'\.$"):
        mtcli.cli.get_command(ctx, name).make_context(name, list(args), parent=ctx)


def test_an_unknown_option_is_still_unknown(run_cli):
    """A guard: an option that metaTOR never had is refused as unknown, as before,
    by click's own message."""
    result = run_cli("pipeline", "map.metacool", "--bogus", "-a", "a.fa", exit_code=2)
    assert "No such option" in result.output and "--bogus" in result.output


@pytest.mark.parametrize(
    "args",
    [
        ["{missing}"],
        ["{collection}"],
        ["{map}", "-m", "abundance"],
        ["{map}", "-p", "bad.prefix"],
        ["{map}", "-a", "missing.fa"],
    ],
    ids=[
        "map not there",
        "collection",
        "normalization without depth",
        "prefix",
        "no assembly",
    ],
)
def test_pipeline_refused_keeps_the_bins_of_an_earlier_run(args, tmp_path, run_cli):
    """A refused run with --force (a rebinning into the same directory) stops
    before replacing the bins of the earlier run."""
    names = {
        "map": _contig_map(tmp_path / "input.metacool"),
        "missing": tmp_path / "missing.metacool",
        "collection": _collection(tmp_path / "maps.mmetacool", "lib1"),
    }
    args = [arg.format(**names) for arg in args]
    out = tmp_path / "out"
    # The outputs of an earlier run: its clusters, and its plots.
    names = ["overlapping_clusters", "recursive_clusters", "final_clusters"]
    bins = [out / name / "bin.fa" for name in names] + [out / "plot" / "bins_distribution.png"]
    for path in bins:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    with pytest.raises(ValueError):
        run_cli("pipeline", "-a", TESTS_DATA / "assembly.fa", *args, "-o", out, "-F", "-T", tmp_path / "tmp")
    assert all(path.exists() for path in bins)


def test_pipeline_refuses_earlier_bins_before_reading_the_assembly(tmp_path, monkeypatch, run_cli):
    """Without --force, the clusters of an earlier run are refused before the
    assembly is read (minutes on a large metagenome)."""

    def create_contig_data(*args):
        raise AssertionError("assembly read before the refusal")

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    out = tmp_path / "out"
    (out / MAGS_WITH_MGES).mkdir(parents=True)
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    with pytest.raises(ValueError):
        run_cli("pipeline", _contig_map(tmp_path / "map.metacool"), "-a", assembly, "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"ERROR :: {out / MAGS_WITH_MGES} already existed." in log.read_text()


def test_final_clusters_are_those_of_their_fasta_files(tmp_path, caplog):
    """The final clusters of a validation are those of the FASTA files it wrote: a
    file of a cluster without a contig, or a cluster without a file, is
    refused, named, and the map is not written."""
    directory = tmp_path / "final_clusters"
    directory.mkdir()
    for name in ["m_00001_00000", "m_00002_00000"]:
        (directory / f"{name}.fa").write_text(">c\nACGT\n")
    mtcli.check_final_clusters({"c1": "m_00001_00000", "c2": "m_00002_00000"}, str(directory))
    with pytest.raises(ValueError) as refused:
        mtcli.check_final_clusters({"c1": "m_00001_00000", "c3": "m_00003_00000"}, str(directory))
    assert str(refused.value) == (
        f"The final clusters of the run are not those of their FASTA files in {directory}: no contig in the final "
        "cluster(s) 'm_00002_00000', and no FASTA file for the final cluster(s) 'm_00003_00000'. The map of the "
        "clusters is not written."
    )
    assert caplog.records[-1].levelname == "ERROR" and caplog.records[-1].getMessage() == str(refused.value)


def test_pipeline_leaves_the_directory_of_the_assembly_alone(stub_pipeline, tmp_path, run_cli):
    """The run writes nothing beside the assembly, and removes nothing there: the
    pyfastx index of another program (<assembly>.fxi), which the run removed at
    its end, is left as it was."""
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    index = Path(f"{assembly}.fxi")
    index.write_text("another program's index")
    out = tmp_path / "out"
    run_cli("pipeline", _contig_map(tmp_path / "map.metacool"), "-a", assembly, "-o", out, "-T", tmp_path / "tmp")
    assert stub_pipeline["contig_contacts"] and index.read_text() == "another program's index"
    # tests_data: the link to the data of the tests (see tests/conftest.py);
    # the temporary directory of the run is removed at its end.
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "assembly.fa", "assembly.fa.fxi", "map.metacool", "out", "tests_data"
    ]  # fmt: skip


@pytest.mark.parametrize("flag", [None, "--random", "--reassign"], ids=["mge", "--random", "--reassign"])
def test_mge_default_outdir(flag, tmp_path, monkeypatch, run_cli):
    """Without --outdir, mge writes the FASTA files of its MGE clusters to
    fasta/mges of the current directory (fasta/mges_random with --random),
    which it makes, with the assembly that the map was
    made from, and stores the round mge (mge_random) in the map, and its MGE
    clusters as genomes with --reassign; its temporary directory is removed at
    its end, and nothing is written next to the assembly."""
    import metacooler
    import pandas as pd

    import metator.clusters as mtcl

    random = flag == "--random"
    metacool = _contig_map(tmp_path / "map.metacool")
    calls = {}
    summary = pd.DataFrame({"BinName": ["MetaTOR_MGE_00001"], "Contigs": ["c2"]})

    def mge_binning(metacool, mges_list, threshold_bin, random):
        calls["binning"] = (metacool, mges_list, threshold_bin, random)
        return summary, {1: {"Contigs": ["c2"], "Score": float("nan")}}

    def write_mges_fasta(metacool, assembly, mge_bins, out_dir, temp_dir):
        calls["fasta"] = (metacool, assembly, out_dir)
        assert os.path.isdir(temp_dir) and os.path.isdir(out_dir)

    def write_map(metacool, outfile, planned, what, mge, reassign):
        (clusters, members), entry = mge
        calls["map"] = (metacool, outfile, list(clusters["round"]), list(members["contig"]), entry["name"], reassign)
        assert isinstance(planned, metacooler.MetaCooler)

    monkeypatch.setattr("metator.mge.mge_binning", mge_binning)
    monkeypatch.setattr("metator.mge.write_mges_fasta", write_mges_fasta)
    monkeypatch.setattr("metator.cli.write_map", write_map)
    assembly = tmp_path / "inputs" / "assembly.fa"
    assembly.parent.mkdir()
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    mges = tmp_path / "mges.txt"
    mges.write_text("c2\n")
    run_cli("mge", metacool, "-a", assembly, "-m", mges, *([flag] if flag else []))
    assert calls["binning"] == (str(metacool), ["c2"], 0.8, random)
    fasta = "mges_random" if random else "mges"
    assert calls["fasta"] == (str(metacool), str(assembly), os.path.join(".", "fasta", fasta))
    name = mtcl.MGE_RANDOM if random else mtcl.MGE
    # The map of the output directory, the current one: the map itself, in
    # place.
    assert calls["map"] == (str(metacool), "map.metacool", [name], ["c2"], name, flag == "--reassign")
    assert sorted(path.name for path in assembly.parent.iterdir()) == ["assembly.fa"]
    assert not (tmp_path / "tmp").exists()


@pytest.mark.parametrize(
    ("command", "reader"),
    [("partition", "metator.network.map_network"), ("validation", "metator.network.map_contig_table")],
)
def test_a_step_that_cannot_read_what_the_map_stores(command, reader, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The network (partition), or the contig data (validation), that the map
    cannot give as the run reads the network stored in it stops the run with
    its reason, logged as an error, and leaves the map as it was."""
    metacool = _contig_map(tmp_path / "input.metacool", network=True, partition=command == "validation")
    before = metacool.read_bytes()

    def unreadable(*args, **kwargs):
        raise mtmc.MetacoolError("unreadable")

    monkeypatch.setattr(reader, unreadable)
    with pytest.raises(mtmc.MetacoolError, match=r"^unreadable$"):
        run_cli(command, *_step_run(command, metacool, tmp_path))
    (log,) = (tmp_path / "out").glob(f"metator_{command}_*.log")
    assert "ERROR :: unreadable\n" in log.read_text()
    assert metacool.read_bytes() == before and "partition" not in stub_pipeline
