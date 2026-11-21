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
import metator.metacool as mtmc
from metator.version import __version__
from tests.locks import described, held_lock, killed_lock, lock_file, refusal

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"

# The writer itself, which the stub pipeline replaces.
mtmc_reassign_metacool = mtmc.reassign_metacool

# The docopt command line that rich-click replaced (that of 1.3.10, with the
# mge options renamed and host removed), as docopt read it: see the commit
# that added the file. Commands and options added since are added to it:
# metacool, the --no-metacool and --output-map options of pipeline, and the
# --output-map option of partition and validation. The 2.0 command line
# changed it on purpose (see docs/migration.md): the map that network and
# pipeline start from is their first argument (MAP), required, and their
# options of reads, alignments and pairs (--forward, --reverse, --start,
# --aligner, --aligner-mode, --min-quality) and --edge are gone; the map
# that partition and validation update is their first argument too; the
# commands qc, contactmap, metacool, scaffold and pairs are gone, and so are
# the --scaffold and --junctions options of pipeline, --self-contacts of
# network and --cluster-matrix of partition, validation and pipeline.
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
    """The help shows the one input, the map, the first argument, which
    hicstuff or metacooler make; and that every pair counts."""
    output = " ".join(run_cli(command, "--help").output.split())
    assert f"Usage: metator {command} [OPTIONS] MAP" in output
    assert f"metator {command} MAP -a ASSEMBLY" in output and "-1 FORWARD" not in output
    assert "hicstuff pipeline -M metacool" in output and "metacooler create pairs" in output
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
    missing = "MAP"
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
    """Kept with --no-clean-up, and when the run fails. An empty --tmpdir
    is ./tmp."""
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
    log = tmp_path / "run.log"
    monkeypatch.setattr(sys, "argv", ["/usr/bin/metator", "partition", "-p", "my sample", "a.metacool"])
    mtcli.generate_log_header(str(log))
    mtcli.generate_log_footer(str(log))
    lines = log.read_text().splitlines()
    assert lines[0] == f"## MetaTOR: v{__version__} log file"
    assert lines[1].startswith("## date: ")
    assert lines[2] == "## Command used: metator partition -p 'my sample' a.metacool"
    assert lines[-1].startswith("## date: ")


# The bookkeeping of the commands, with the heavy steps replaced by stubs.


# The contacts of a map as metator.network reads them, which stub_pipeline
# records: the maps that the tests make read them unrecorded.
CONTIG_CONTACTS = mtcli.mtn.contig_contacts


def _no_step(calls):
    """No step of the stub pipeline ran (see stub_pipeline): a partition or
    validation run reads the contacts of its map as it reads the network
    stored in it (see metator.network.stored_network), before anything
    else."""
    return not set(calls) - {"contig_contacts"}


class _AnyContig(dict):
    """The contig data of the stub create_contig_data: every contig of the
    map, of 1000 bp, measured on nothing but its GC content (no depth file,
    no enzyme), as create_contig_data gives a contig."""

    def __missing__(self, name):
        return {"length": 1000, "GC": 50.0, "coverage": "-", "RS": "-"}


@pytest.fixture
def stub_pipeline(monkeypatch, tmp_path):
    """Replace the steps of the pipeline by stubs that write their outputs
    and record their arguments. The network of the map is read and stored
    for real (see metator.network), from the contig data of the stub
    create_contig_data."""
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
        out_dir = args[5]
        # What the network step left in the output directory.
        calls["partition_outdir"] = sorted(os.listdir(out_dir))
        Path(out_dir, "contig_data_partition.txt").write_text("x")
        return os.path.join(out_dir, "contig_data_partition.txt")

    def recursive_decontamination(*args):
        calls["recursive_decontamination"] = args
        # pyfastx indexes the assembly while extracting the bins.
        Path(f"{args[1]}.fxi").write_text("x")

    def reassign_metacool(metacool, contig_data, outfile, **kwargs):
        # The input map must still be there.
        calls["reassign_metacool"] = (metacool, os.path.exists(metacool), contig_data, outfile, kwargs)

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    monkeypatch.setattr("metator.network.contig_data_from_map", contig_data_from_map)
    monkeypatch.setattr("metator.network.contig_contacts", recorded_contig_contacts)
    monkeypatch.setattr("metator.partition.partition", partition)
    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)
    monkeypatch.setattr("metator.metacool.reassign_metacool", reassign_metacool)
    return calls


def _assembly(tmp_path):
    """A copy of an assembly: the stub validation writes its pyfastx index."""
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    return assembly


def _final_tables(monkeypatch, contig, final_bin="metator_00001_00000", summary=None, hits=("-", "-")):
    """Have the stub validation also write the final contig table, where the
    contig ``contig`` is in the bin ``final_bin`` (ND: in none) and c2 is in
    none, with these hits (none counted, by default), and the bin summary
    (``summary``, if given)."""
    stub = mtcli.mtv.recursive_decontamination

    def recursive_decontamination(*args, **kwargs):
        out_dir = args[7]
        Path(out_dir, "contig_data_final.txt").write_text(
            "\tID\tName\tSize\tGC_content\tHit\tFinal_cluster\n"
            f"0\t1\t{contig}\t100\t40.5\t{hits[0]}\t{final_bin}\n1\t2\tc2\t50\t35.0\t{hits[1]}\tND\n"
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
def test_a_run_refuses_an_assembly_that_pyfastx_cannot_read(compression, command, stub_pipeline, tmp_path, run_cli, logged):
    """An assembly compressed with bzip2 or zip, which metaTOR reads but
    pyfastx, in the partition and the validation, does not, is refused
    before any work: the clusters of an earlier run stay, even with -F
    (pipeline)."""
    import bz2
    import zipfile

    fasta = tmp_path / "assembly.fa"
    if compression == "bzip2":
        fasta.write_bytes(bz2.compress(b">c1\nACGT\n"))
    else:
        with zipfile.ZipFile(fasta, "w") as archive:
            archive.writestr("assembly.fa", ">c1\nACGT\n")
    earlier = tmp_path / "out" / "final_clusters" / "bin.fa"
    earlier.parent.mkdir(parents=True)
    earlier.write_text(">old\nA\n")
    metacool = _contig_map(tmp_path / "map.metacool")
    with pytest.raises(ValueError):
        force = ["-F"] if command == "pipeline" else []
        run_cli(command, metacool, "-a", fasta, *force, "-o", tmp_path / "out", "-T", tmp_path / "tmp")
    assert earlier.exists() and not stub_pipeline
    message = f"{fasta} is compressed with {compression}: give the assembly as a FASTA file, plain or compressed with gzip."
    assert ("ERROR", message) in logged()


def test_pipeline_from_a_map_keeps_the_contig_names(stub_pipeline, tmp_path, run_cli):
    """A guard: a run from a map does not check the contig names of the
    assembly, as it builds no map of the contigs (nor does the map of its
    clusters: see test_a_non_ascii_contig_that_the_map_lacks_is_left_out in
    test_metacool.py)."""
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(NON_ASCII_ASSEMBLY)
    metacool = _contig_map(tmp_path / "map.metacool")
    run_cli("pipeline", metacool, "-a", fasta, "-o", tmp_path / "out", "-T", tmp_path / "tmp")
    assert stub_pipeline["contig_contacts"][0] == str(metacool)


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
    (log,) = out.glob("metator_pipeline_*.log")
    assert "Removing" not in log.read_text()


def test_pipeline_map_stores_the_normalization_of_the_network(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """With the real writers: the map of the clusters stores the normalization
    of the network of the run, empirical_hit by default, as the weights
    metator_empirical_hit of its bins, one bin per contig: 1/sqrt(hits)."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _contig_data(monkeypatch, CONTIGS)
    _final_tables(monkeypatch, "c1", hits=(4, 1))
    metacool = _contig_map(tmp_path / "map.metacool")
    out = tmp_path / "out"
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    (log,) = out.glob("metator_pipeline_*.log")
    assert "WARNING" not in log.read_text()
    import metacooler

    path = str(metacool)
    names = list(metacooler.MetaCooler(path).contigs.name)
    weights = metacooler.read_weights(path, "metator_empirical_hit")
    assert dict(zip(names, weights, strict=True)) == pytest.approx({"c1": 1 / 2, "c2": 1.0})


@pytest.mark.parametrize("case", ["a contig of another assembly", "no contact"])
def test_pipeline_network_of_the_pairs_of_its_map(case, network_from_a_map, tmp_path, run_cli):
    """The network of a run is that of the pairs that its map holds. A pair
    on a contig that the assembly lacks is not in the map (metacooler left
    it out as it made the map); a map that holds no contact (no pair within
    the contigs) stops the run with the MetacoolError of the network, its
    reason logged as an error."""
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


@pytest.mark.parametrize("no_metacool", [False, True])
def test_pipeline_from_a_metacool(no_metacool, stub_pipeline, tmp_path, run_cli):
    """From a contact map: the network is read from it, and the map of the
    clusters is written back from it. The input map is not removed first:
    it is replaced once the new map is validated."""
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool")
    before = metacool.read_bytes()
    options = ["--no-metacool"] if no_metacool else []
    run_cli("pipeline", metacool, "-a", _assembly(tmp_path), "-t", "3", "-m", "length",
            "-d", TESTS_DATA / "depth.txt", "-o", out, "-T", tmp_path / "tmp", *options)  # fmt: skip

    calls = stub_pipeline
    fasta, depth, enzyme = calls["create_contig_data"]
    assert (fasta, depth, enzyme) == (str(tmp_path / "assembly.fa"), str(TESTS_DATA / "depth.txt"), None)
    assert calls["contig_contacts"] == (str(metacool),)
    # The contig data and the network of the map, in memory: c1 and c2 (the
    # stub contig data gives them 1000 bp), one contact.
    table, network = calls["partition"][2], calls["partition"][4]
    assert list(table.Name) == ["c1", "c2"] and list(table.Hit) == [3, 1]
    assert network.to_dict("list") == {
        "contig1": [1], "contig2": [2], "weight": [pytest.approx((1000 / 3) ** 0.5 * 1000**0.5, rel=1e-15)]
    }  # fmt: skip
    assert calls["partition"][10] == 3
    assert metacool.read_bytes() == before
    # No table of the network: it is in the map.
    assert not list(out.glob("*network*"))
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
    assert parameters["map"] == str(metacool) and "start" not in parameters
    assert kwargs["metadata"]["command"].startswith(f"metator pipeline {metacool} -a")
    # The network, stored as metator network stores it, and its
    # normalization.
    network = kwargs["network"]
    assert (kwargs["normalization"], network.normalization) == ("length", "length")
    pd.testing.assert_frame_equal(network.table, table)
    assert network.entry["command"] == kwargs["metadata"]["command"] and network.entry["version"] == __version__
    assert Path(network.tmp_dir).parent == tmp_path / "tmp"
    # No column was read from the map (the stub reads none).
    assert "columns_from_map" not in kwargs["metadata"] and "columns_from_map" not in network.entry
    assert f"INFO :: {metacool} is the input map {metacool}: it is replaced once the map of the clusters is validated." in logged


def _contig_map(path, lengths=None, network=False, normalization="empirical_hit"):
    """A .metacool map of the contigs c1 and c2 of the stub validation's
    tables (see _final_tables), made before any binning: every contig in
    __unbinned__, one bin per contig. With network, the network of its
    contigs stored in it, as metator network stores it (see
    metator.metacool.store_network), which partition and validation read:
    the contacts of c1 and c2 normalized by their hits (empirical_hit, by
    default)."""
    import metacooler
    import pandas as pd

    lengths = lengths or {"c1": 100, "c2": 50}
    pairs = path.with_suffix(".pairs")
    pairs.write_text("## pairs format v1.0\nr1\tc1\t10\tc2\t20\t+\t-\nr2\tc1\t10\tc1\t90\t+\t-\n")
    bins = pd.DataFrame({"chrom": list(lengths), "start": 0, "end": list(lengths.values())})
    metacooler.from_pairs(str(path), pairs=str(pairs), bins=bins, contig_genome_map=None, ps=False)
    pairs.unlink()
    if network:
        names = metacooler.MetaCooler(str(path)).contigs["name"]
        contacts = CONTIG_CONTACTS(str(path))
        table = mtcli.mtn.contig_data_table(_AnyContig(), names, mtcli.mtn.contig_hits(contacts, len(names)))
        mtmc.store_network(str(path), str(path), mtmc.Network(table, normalization), force=True)
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
        kwargs = stub_pipeline["reassign_metacool"][4]
        assert kwargs["normalization"] == kwargs["network"].normalization == "theoritical_hit"
        entry = kwargs["network"].entry
    else:
        import metacooler

        entry = metacooler.MetaCooler(str(metacool)).metadata["metator"]["network"]
        assert entry["normalization"] == "theoritical_hit"
    assert entry["columns_from_map"] == ["Restriction_site", "Shotgun_coverage"]


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
    """A run from a map with -F (a rebinning in the output directory of an
    earlier run) that its columns refuse once the assembly is read leaves
    the clusters of the earlier run: the contig data is read before they
    are replaced (they were removed, then the run refused)."""
    monkeypatch.setattr("metator.network.contig_data_from_map", mtn_contig_data_from_map)
    contig_data = {name: {"length": length, "GC": 50.0, "coverage": "-", "RS": "-"}
                   for name, length in [("c1", 100), ("c2", 50), ("c3", 30)]}  # fmt: skip
    monkeypatch.setattr("metator.network.create_contig_data", lambda *args: contig_data)
    metacool = _map_with_columns(tmp_path, {"Restriction_site": [3, 2]})
    out = tmp_path / "out"
    names = ["overlapping_clusters", "recursive_clusters", "final_clusters"]
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


def _merged_map(tmp_path, normalization):
    """The map of two libraries (two copies of _contig_map with the network
    of their contigs, normalized by normalization) summed by metacooler
    merge, which keeps the contigs table of the first, and the metadata of
    both: the network that metator network stored in the first, whose Hit
    is no longer that of the contacts. A normalization other than None
    loses its weights in the merge; None stores none to lose."""
    from metacooler.cli.main import cli as metacooler_cli

    library = _contig_map(tmp_path / "library_1.metacool", network=True, normalization=normalization)
    copy = shutil.copyfile(library, tmp_path / "library_2.metacool")
    merged = tmp_path / "merged.metacool"
    result = CliRunner().invoke(metacooler_cli, ["merge", str(library), str(copy), "-o", str(merged), "-q"])
    assert result.exit_code == 0, result.output
    return merged


def _collection(path, *entries, network=False):
    """A .mmetacool collection at path, with a map of the contigs
    (_contig_map, with its network with network) as each entry of
    entries."""
    import metacooler

    maps = {name: _contig_map(path.with_name(f"{path.stem}_{name}.metacool"), network=network) for name in entries}
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
    assert (out / "cluster_summary.txt").exists()
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
    # Staged by metacooler in a private directory next to the map, removed.
    assert Path(staged).parent.parent == out and Path(staged).parent.name.startswith(".metator.metacool.")
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
    """A lock file of the map metacool that names a process pid of host (this
    host by default), but not the machine and the PID namespace where it
    runs, which metacooler's locks name: metacooler cannot check it (as a
    lock of another host, or one that metaTOR wrote before its map writes
    went through metacooler), and never takes it over."""
    import socket

    holder = {"pid": pid, "host": host or socket.gethostname(), "start": start, "command": "metator pipeline",
              "uuid": "0123456789abcdef0123456789abcdef"}  # fmt: skip
    lock = Path(f"{metacool}.lock")
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
        rounds.write(args[7])
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
    stores the round overlapping just before; mge MAP its round mge. A map
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
    """partition MAP, validation MAP and mge MAP warn, as they take the lock
    of the map, of what they replace in it in place (see
    warn_of_what_is_replaced); with --output-map, they do not."""
    from contextlib import ExitStack

    metacool = _map_of_an_earlier_run(tmp_path / "input.metacool")
    for output_map in [str(tmp_path / "copy.metacool"), None]:
        with ExitStack() as held:
            mtcli.step_map(command, str(metacool), output_map, False, held)
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
    assert "contig_contacts" not in stub_pipeline and output_map.read_text() == "an older map"
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
    assert "contig_contacts" not in stub_pipeline


def test_pipeline_from_a_map_through_a_symbolic_link(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map given through a symbolic link: the map it points to is
    updated in place, as metacooler edits the target of a link, and the
    link keeps pointing at it. Its lock and its staged copy are next to the
    target, and gone once the run has ended."""
    import metacooler

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    (tmp_path / "maps").mkdir()
    target = _contig_map(tmp_path / "maps" / "target.metacool")
    before = target.read_bytes()
    link = tmp_path / "link.metacool"
    link.symlink_to(target)
    run_cli("pipeline", *_map_run(link, tmp_path))
    assert link.is_symlink() and link.resolve() == target and target.read_bytes() != before
    assert list(metacooler.MetaCooler(str(target)).genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert sorted(path.name for path in (tmp_path / "maps").iterdir()) == ["target.metacool"]


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
    """A second run on a map that another run holds (here, another process
    that holds its lock) is refused before any work, with metacooler's
    reason (who holds it, and the files to delete if it is gone) and where
    else to write the map, and leaves the lock to its holder."""
    metacool = _contig_map(tmp_path / "input.metacool")
    with held_lock(metacool) as holder:
        message = f"{refusal(metacool, holder)} Or write this run's map elsewhere with --output-map."
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("pipeline", *_map_run(metacool, tmp_path))
        assert "contig_contacts" not in stub_pipeline
        assert json.loads(lock_file(metacool).read_text()) == holder
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("holder", ["another host", "no machine named", "unreadable", "not an object"])
def test_pipeline_from_a_map_refused_with_a_lock_it_cannot_check(holder, stub_pipeline, tmp_path, run_cli):
    """A lock that metacooler cannot tell stale is never taken over: that of
    a process of another host; one that does not name the machine and the
    PID namespace of its process (as metaTOR's own locks did, before its
    map writes went through metacooler), even of a process that is gone;
    or one that says nothing readable (not JSON, or JSON but not an
    object). The run is refused, with the files to delete, and the lock
    left."""
    import subprocess

    metacool = _contig_map(tmp_path / "input.metacool")
    lock = Path(f"{metacool}.lock")
    if holder in ("unreadable", "not an object"):
        lock.write_text("{not json" if holder == "unreadable" else "[1]")
        reason = f"{lock} does not say which process holds {metacool}: if none is writing it, delete {lock}."
    else:
        gone = subprocess.Popen([sys.executable, "-c", "pass"])
        gone.wait()
        _lock(metacool, 1, host="elsewhere") if holder == "another host" else _lock(metacool, gone.pid)
        record = json.loads(lock.read_text())
        reason = (
            f"{metacool} is locked by a process this host cannot check ({described(record)}: another host, machine "
            f"or PID namespace): if it is no longer running, delete {lock} and "
            f"{tmp_path}/.input.metacool.{record['uuid']}.*.d/ (its staged copies)."
        )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(reason)} Or write this run's map elsewhere"):
        run_cli("pipeline", *_map_run(metacool, tmp_path))
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
        "elsewhere with --output-map."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert "contig_contacts" not in stub_pipeline
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}" in log.read_text()


@pytest.mark.parametrize("command", ["pipeline", "partition", "validation"])
def test_a_read_only_map_is_refused_before_any_work(command, stub_pipeline, tmp_path, run_cli):
    """A map without write permission is refused before any work, with how
    to write a copy instead: metacooler refuses to replace it, which it
    would do at the end of the run. The map and its permission are left as
    they were, and so is the lock."""
    metacool = _contig_map(tmp_path / "input.metacool")
    metacool.chmod(0o444)
    before = metacool.read_bytes()
    message = (
        f"{metacool} is read-only (mode -r--r--r--): the run would replace it with the map of its clusters. Make it "
        "writable (chmod u+w), or write the map of the run to another file with --output-map."
    )
    arguments = _map_run(metacool, tmp_path) if command == "pipeline" else _step_run(command, metacool, tmp_path)
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *arguments)
    assert not stub_pipeline
    assert metacool.read_bytes() == before and metacool.stat().st_mode & 0o777 == 0o444
    assert not Path(f"{metacool}.lock").exists()


@pytest.mark.parametrize(
    "case",
    [
        "missing map",
        "no map written",
        "an entry",
        "a directory",
        "the collection of the input",
        "the pyfastx index of the assembly",
    ],
)
def test_pipeline_output_map_refusals(case, stub_pipeline, tmp_path, run_cli):
    """What --output-map cannot be, and an input map that is not there, are
    refused before any work (they failed at the end of the run, or on a
    bare error while locking)."""
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
    elif case == "a directory":
        output = tmp_path / "maps"
        output.mkdir()
        options = ["--output-map", output, "-F"]
        message = f"{output} is a directory: --output-map takes the path of a .metacool file."
    elif case == "the pyfastx index of the assembly":
        # Not there yet: the run writes it, then removes it at its end.
        output = f"{tmp_path / 'assembly.fa'}.fxi"
        options = ["--output-map", output]
        message = (
            f"--output-map {output} is the pyfastx index of the assembly {tmp_path / 'assembly.fa'}, which the run "
            "removes at its end: write the map elsewhere."
        )
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


def test_pipeline_output_map_in_a_new_directory(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """--output-map in a directory that is not there yet: the directory is
    made, as the map will be written there (the lock failed on it)."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    output = tmp_path / "new" / "copy.metacool"
    run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", output))
    assert output.exists() and not Path(f"{output}.lock").exists()


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


def test_a_map_is_written_to_a_file_only(tmp_path):
    """A map is not written to a URI (file::group), which metacooler reads
    as a group inside the file: an existing file would be truncated."""
    metacool = _contig_map(tmp_path / "input.metacool")
    before = metacool.read_bytes()
    contig_data = tmp_path / "contig_data_final.txt"
    contig_data.write_text("ID\tName\tSize\tFinal_cluster\n1\tc1\t100\tND\n")
    outfile = f"{metacool}::/"
    message = f"^{re.escape(outfile)} is not the path of a file: a map is written to a file of its own, not to a group"
    with pytest.raises(mtmc.MetacoolError, match=message):
        mtmc_reassign_metacool(str(metacool), str(contig_data), outfile, force=True)
    assert metacool.read_bytes() == before


def test_map_directory_that_cannot_be_made(stub_pipeline, tmp_path, run_cli):
    """The directory of the map to write cannot be made (a file is in its
    way): refused before any work, the reason logged."""
    blocking = tmp_path / "afile"
    blocking.write_text("x")
    output = blocking / "x.metacool"
    metacool = _contig_map(tmp_path / "input.metacool")
    message = f"The directory of --output-map {output} cannot be made (File exists): write the map elsewhere."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", output))
    assert not stub_pipeline


@pytest.mark.parametrize(
    ("command", "output"),
    [("pipeline", False), ("pipeline", True), ("partition", False)],
    ids=["pipeline in place", "pipeline --output-map", "partition in place"],
)
def test_a_map_changed_during_the_run_is_refused_at_its_end(command, output, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map changed while a run works from it (here, rewritten by a program
    that takes no lock of metacooler, once the partition is done) is refused
    when the run writes its map: the run read the network from the map as it
    was at its start, and would write its clusters into another map. The
    map of the run is not written, with the reason; --output-map is not
    written either, and the map is left as the other program left it."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    if command == "partition":
        _partition_table(monkeypatch)
    metacool = _contig_map(tmp_path / "input.metacool", network=command == "partition")
    partition = mtcli.mtp.partition

    def partition_then_rewrite(*args, **kwargs):
        written = partition(*args, **kwargs)
        rewritten = tmp_path / "rewritten.tmp"
        shutil.copyfile(metacool, rewritten)
        os.replace(rewritten, metacool)
        return written

    monkeypatch.setattr(mtcli.mtp, "partition", partition_then_rewrite)
    options = ["--output-map", tmp_path / "copy.metacool"] if output else []
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
    assert not (tmp_path / "copy.metacool").exists()
    import metacooler

    assert list(metacooler.MetaCooler(str(metacool)).genomes.name) == ["__unbinned__"]
    assert metacool.read_bytes() == rewritten and not lock_file(metacool).exists()


def _output_map_run(command, metacool, tmp_path, *options):
    """The arguments of a pipeline or partition run from the map metacool,
    with options (--output-map, -F...)."""
    if command == "pipeline":
        return _map_run(metacool, tmp_path, *options)
    return _step_run(command, metacool, tmp_path, "-s", "100", *options)


@pytest.mark.parametrize("command", ["pipeline", "partition"])
@pytest.mark.parametrize("kind", ["symbolic link", "hard link"])
def test_output_map_a_link_to_the_input_map(kind, command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """--output-map FILE, a link to the input map. A symbolic link at FILE is
    a map of its own: refused without -F, as any FILE that exists, and with
    -F replaced by the map of the run, the input map left as it was. A hard
    link is the input map under another name: refused before any work,
    with -F too. Both updated the input map in place, and the hard link
    kept the old map under FILE, though the log said FILE was written."""
    import metacooler

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    if command == "partition":
        _partition_table(monkeypatch)
    (tmp_path / "maps").mkdir()
    metacool = _contig_map(tmp_path / "maps" / "input.metacool", network=command == "partition")
    before = metacool.read_bytes()
    out = tmp_path / "run.metacool"
    if kind == "symbolic link":
        out.symlink_to(metacool)
    else:
        os.link(metacool, out)
    if kind == "hard link":
        message = (
            f"--output-map {out} is a hard link of the input map {metacool}: the same file under another name. Give "
            "--output-map another file, or leave it out to update the map in place."
        )
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli(command, *_output_map_run(command, metacool, tmp_path, "--output-map", out, "-F"))
        assert "partition" not in stub_pipeline
        assert metacool.read_bytes() == before and os.path.samefile(metacool, out)
        return
    message = f"{out} exists: give -F to replace it with the map of this run, or another --output-map."
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *_output_map_run(command, metacool, tmp_path, "--output-map", out))
    run_cli(command, *_output_map_run(command, metacool, tmp_path, "--output-map", out, "-F"))
    assert metacool.read_bytes() == before and not out.is_symlink()
    mc = metacooler.MetaCooler(str(out))
    if command == "pipeline":
        assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    else:
        assert list(mc.clusters.name) == ["metator_00001_00000"]
    assert sorted(path.name for path in (tmp_path / "maps").iterdir()) == ["input.metacool"]
    assert not [path.name for path in tmp_path.iterdir() if path.name.startswith(".") or path.suffix == ".lock"]


@pytest.mark.parametrize("case", ["a map that another process holds", "a link to nothing"])
def test_output_map_a_symbolic_link_locks_its_own_name(case, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """--output-map FILE, a symbolic link: the run locks FILE's own name, as
    it replaces the link, never what it points to. A link to a map that
    another process holds (a run that updates it in place) is no refusal,
    and a link to nothing is replaced (-F). The run locked what the link
    points to: it was refused on the lock of that map, or on a lock that
    cannot be created beside nothing."""
    from contextlib import ExitStack

    import metacooler

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "run.metacool"
    with ExitStack() as held:
        if case == "a link to nothing":
            out.symlink_to(tmp_path / "gone" / "old.metacool")
        else:
            (tmp_path / "elsewhere").mkdir()
            target = _contig_map(tmp_path / "elsewhere" / "other.metacool")
            out.symlink_to(target)
            held.enter_context(held_lock(target))
            before = target.read_bytes()
        run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", out, "-F"))
        if case == "a map that another process holds":
            assert target.read_bytes() == before
    assert not out.is_symlink()
    assert list(metacooler.MetaCooler(str(out)).genomes.name) == ["__unbinned__", "metator_00001_00000"]


def test_hold_lock_follows_the_links_of_the_map_in_place_only(monkeypatch):
    """hold_lock locks the map that a run writes: in place, the file that a
    symbolic link points to, as metacooler edits it (metacooler.lock follows
    links); with --output-map, FILE's own name (follow_symlinks=False), as
    a link there is replaced, not its target."""
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

    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    if command == "partition":
        _partition_table(monkeypatch)
    metacool = _contig_map(tmp_path / "input.metacool", network=command == "partition")
    table = tmp_path / "extra.txt"
    table.write_text("Name\tnote\nc1\tx\nc2\ty\n")
    run_map_path = mtcli.run_map_path

    def another_writer_then(*args, **kwargs):
        metacooler.annotate(str(metacool), contigs_table=str(table), contigs_table_format="Name", inplace=True)
        return run_map_path(*args, **kwargs)

    monkeypatch.setattr(mtcli, "run_map_path", another_writer_then)
    run_cli(command, *_output_map_run(command, metacool, tmp_path))
    mc = metacooler.MetaCooler(str(metacool))
    assert "note" in mc.contigs.columns
    if command == "pipeline":
        assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    else:
        assert list(mc.clusters.name) == ["metator_00001_00000"]


def _replaced_once_opened_by(monkeypatch, metacool, reader):
    """Have the map metacool replaced by another file (a copy of itself), as
    a run that updates it in place replaces it, once the function reader
    opens it, before it reads it."""
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
        ("check_map_of_contigs", "validation", ["--output-map", "out.metacool"]),
        ("map_columns", "pipeline", ["--normalization", "RS"]),
        ("warn_of_what_is_replaced", "partition", []),
    ],
)
def test_a_map_replaced_as_the_run_checks_it(reader, command, options, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map replaced as the run reads it to check it (its contigs against
    the contig data of a validation, its columns against the normalization,
    what it holds that the run replaces) stops the run with metaTOR's error,
    logged as an ERROR: it stopped on metacooler's MapChangedError, raised as
    it is, with nothing in the log."""
    metacool = _contig_map(tmp_path / "input.metacool")
    _replaced_once_opened_by(monkeypatch, metacool, reader)
    message = REPLACED_WHILE_READ.format(map=metacool)
    options = [tmp_path / option if option.endswith(".metacool") else option for option in options]
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli(command, *_output_map_run(command, metacool, tmp_path, *options))
    assert "partition" not in stub_pipeline and "recursive_decontamination" not in stub_pipeline
    (log,) = (tmp_path / "out").glob(f"metator_{command}_*.log")
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


def test_output_map_held_by_another_run_is_refused_without_the_hint(stub_pipeline, tmp_path, run_cli):
    """A run with --output-map FILE, which another process holds, is refused
    before any work with metacooler's reason only: the hint to write the map
    elsewhere with --output-map, which a run in place gets, is no help to a
    run that gives one."""
    metacool = _contig_map(tmp_path / "input.metacool")
    out = tmp_path / "copy.metacool"
    with held_lock(out) as holder:
        message = refusal(out, holder)
        with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
            run_cli("pipeline", *_map_run(metacool, tmp_path, "--output-map", out))
    assert "contig_contacts" not in stub_pipeline
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    assert f"ERROR :: {message}\n" in log.read_text()


def test_pipeline_from_a_map_takes_a_stale_lock(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The lock of a run that is gone (killed with SIGKILL, or SIGTERM:
    neither metaTOR nor metacooler handles signals), and the copy that it
    staged next to the map, are taken over and removed by the next run on
    this host, each with a warning in its log; the run goes on, and removes
    its own lock at its end."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool")
    holder = killed_lock(metacool)
    staged = tmp_path / f".input.metacool.{holder['uuid']}.1.d"
    staged.mkdir()
    (staged / "input.metacool").write_text("the copy that the killed run staged")
    run_cli("pipeline", *_map_run(metacool, tmp_path))
    assert not lock_file(metacool).exists() and not staged.exists()
    (log,) = (tmp_path / "out").glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"WARNING :: Took over {lock_file(metacool)}: the process that held it ({described(holder)}) is gone." in logged
    assert f"WARNING :: Removed {staged.name}, left beside {metacool} by a writer that no longer holds its lock." in logged


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


def test_network_from_a_metacool(monkeypatch, tmp_path, run_cli):
    """metator network MAP reads the network from the map, with the contig
    data of the assembly, and stores it in the map: in --output-map for an
    entry of a collection, which is left as it was. The output directory
    holds the log of the run, and nothing else."""
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
            "-T", tmp_path / "tmp", "--output-map", "network.metacool")  # fmt: skip
    assert calls["create_contig_data"] == (str(tmp_path / "assembly.fa"), None, "DpnII")
    # The restriction sites and the coverage that the map gives where the
    # options give none.
    assert calls["contig_data_from_map"] == ("maps.mmetacool::entries/lib1", contig_data, "DpnII", None)
    assert collection.read_bytes() == before
    written = tmp_path / "network.metacool"
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
    assert entry["parameters"]["output_map"] == "network.metacool"
    assert mc.cluster_rounds == [] and list(mc.genomes.name) == ["__unbinned__"]
    (log,) = out.glob("metator_network_*.log")
    assert [path.name for path in out.iterdir()] == [log.name]
    assert not list(tmp_path.glob("*.lock")) and not (tmp_path / "tmp").exists()


def test_network_on_an_entry_without_output_map(stub_pipeline, tmp_path, run_cli):
    """metator network writes its map: an entry of a collection, which a run
    does not update in place, needs --output-map, refused before any
    work."""
    collection = _collection(tmp_path / "maps.mmetacool", "lib1")
    before = collection.read_bytes()
    entry = f"{collection}::entries/lib1"
    message = (
        f"{entry} is an entry of a .mmetacool collection, which a run does not update in place: write the map of its "
        "clusters to another file with --output-map."
    )
    with pytest.raises(mtmc.MetacoolError, match=f"^{re.escape(message)}$"):
        run_cli("network", entry, "-a", _assembly(tmp_path), "-o", tmp_path / "out", "-T", tmp_path / "tmp")
    assert collection.read_bytes() == before and not stub_pipeline


def test_network_updates_its_map_in_place(stub_pipeline, tmp_path, run_cli):
    """metator network MAP stores the network in MAP, in place, under its
    lock, and replaces the network of an earlier run (its columns and its
    weights): run again with another normalization, the map holds that
    one only."""
    import metacooler

    metacool = _contig_map(tmp_path / "input.metacool")
    metacool.chmod(0o640)
    for normalization in ["empirical_hit", "length", "None"]:
        run_cli("network", metacool, "-a", _assembly(tmp_path), "-n", normalization, "-o", tmp_path / "out",
                "-T", tmp_path / "tmp")  # fmt: skip
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


# What the first argument of a network or pipeline run cannot be, refused
# before any work, the reason in its log file.
NOT_A_MAP = (
    "{name} is not a .metacool contact map: make the map of the contigs with hicstuff pipeline -M metacool, from "
    "the reads of the run, or metacooler create pairs, from its pairs files."
)
RUN_INPUT_ERRORS = {
    "reads as the first argument": (["{reads}"], NOT_A_MAP.replace("{name}", "{reads}")),
    "a map that is not there": (["{missing}"], "Cannot read the contact map {missing}: no such file."),
    "a URI of a directory": (["{directory}::/"], "Cannot read the contact map {directory}::/: {directory} is not a file."),
    "a pairs file named as a URI": (["{reads}::/"], NOT_A_MAP.replace("{name}", "{reads}::/")),
    "a pairs file whose name holds ::": (["{colons}"], NOT_A_MAP.replace("{name}", "{colons}")),
    "a map whose name holds ::": (
        ["{colon_map}"],
        "Cannot read the contact map {colon_map}: metacooler reads a path that holds :: as file::group. Rename the map "
        "without ::.",
    ),
}


@pytest.mark.parametrize("command", ["network", "pipeline"])
@pytest.mark.parametrize("case", list(RUN_INPUT_ERRORS))
def test_run_input_refused(command, case, stub_pipeline, tmp_path, run_cli):
    """A run takes a .metacool map, as its first argument: a first argument
    that is no map, not there, or a URI whose file is no map nor a file, is
    refused before any work, the reason in the log file (it failed once the
    run read the map). A file whose name holds :: is that file: a pairs
    file is no map (it was refused as missing), and a map is refused with
    the reason, as metacooler cannot read it."""
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


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_an_extra_argument_is_a_usage_error(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """An argument after the map is one argument too many: a usage error,
    and the command does not start."""
    monkeypatch.chdir(tmp_path)
    result = run_cli(command, "map.metacool", "R2.fq", "-a", "a.fa", "-o", "out", exit_code=2)
    assert "Got unexpected extra argument (R2.fq)" in result.output
    assert not stub_pipeline and not Path("out").exists()


@pytest.mark.parametrize("command", ["network", "pipeline"])
def test_a_run_needs_its_map(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """The map is the one input of a network or pipeline run: without it, a
    usage error, and the command does not start."""
    monkeypatch.chdir(tmp_path)
    result = run_cli(command, "-a", "a.fa", "-o", "out", exit_code=2)
    assert "Missing argument 'MAP'" in result.output
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
    [["-E", "100"], ["--edge", "0"], ["-1", "R1.fq"], ["-S", "pair"], ["-b", "bwa"], ["-q", "30"]],
    ids=["-E", "--edge", "-1", "-S", "-b", "-q"],
)
def test_options_of_reads_and_pairs_are_gone(command, option, stub_pipeline, tmp_path, run_cli):
    """The options of a run from reads, alignments or pairs, and --edge,
    which metaTOR 2.0 removed, are unknown options: click's usage error,
    with no hint. The command does not start."""
    result = run_cli(command, "map.metacool", *option, "-a", "a.fa", "-o", tmp_path / "out", exit_code=2)
    assert f"No such option '{option[0]}'." in " ".join(result.output.replace("│", " ").split())
    assert not stub_pipeline and not (tmp_path / "out").exists()


# partition MAP and validation MAP: the steps run one by one on a map, with
# the steps replaced by stubs (see tests/test_binning_synthetic.py for the
# map that they give).

# The directories of clusters of each step, replaced with -F.
STEP_DIRS = {"partition": ["overlapping_clusters"], "validation": ["recursive_clusters", "final_clusters"]}


def _step_run(command, metacool, tmp_path, *options):
    """The arguments of a partition or validation run on the map metacool
    (for None, a map of the contigs c1 and c2 with their network, which
    metator network stores). A validation reads the tables of the partition, which
    the stub steps do not read, but for the contig data, which the run
    checks the map against: the contigs c1 and c2 of _contig_map, at their
    positions in the map (c1, the longest, first)."""
    tables = []
    if command == "validation":
        contigs = tmp_path / "contig_data.txt"
        if not contigs.exists():
            contigs.write_text(
                "ID\tName\tSize\tHit\tShotgun_coverage\tOverlapping_cluster_ID\n1\tc1\t100\t3\t-\t00001\n2\tc2\t50\t1\t-\t00001\n"
            )
        tables = ["-c", contigs, "-f", "overlapping_clusters"]
    out = ["-o", tmp_path / "out", "-T", tmp_path / "tmp"]
    if metacool is None:
        metacool = tmp_path / "network.metacool"
        if not metacool.exists():
            _contig_map(metacool, network=True)
    return [metacool, "-a", _assembly(tmp_path), *tables, *out, *options]


def _earlier_run(command, tmp_path):
    """The directories of clusters that an earlier run of command left in the
    output directory, each with a FASTA file, which the test checks."""
    files = [tmp_path / "out" / name / "earlier.fa" for name in STEP_DIRS[command]]
    for path in files:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    return files


def _partition_table(monkeypatch, size=150):
    """Have the stub partition write the contig data of a partition of the
    contigs c1 and c2 of _contig_map, as metator.partition writes it: one
    overlapping cluster of both, of size bp (00001: metator_00001_00000)."""
    stub = mtcli.mtp.partition

    def partition(*args):
        result = stub(*args)
        Path(args[5], "contig_data_partition.txt").write_text(
            "ID\tName\tSize\tOverlapping_cluster_ID\tOverlapping_cluster_contigs\tOverlapping_cluster_size\n"
            f"1\tc1\t100\t00001\t2\t{size}\n2\tc2\t50\t00001\t2\t{size}\n"
        )
        return result

    monkeypatch.setattr("metator.partition.partition", partition)


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_help_of_a_step_on_a_map(command, run_cli):
    """The help of partition and validation shows MAP, their first argument,
    which metator network stores the network in, and --output-map."""
    output = " ".join(run_cli(command, "--help").output.split())
    assert f"Usage: metator {command} [OPTIONS] MAP" in output and "metator network" in output
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
    "a directory as --output-map": (
        "{map}",
        ["--output-map", "{directory}"],
        "{directory} is a directory: --output-map takes the path of a .metacool file.",
    ),
    "a lock held by another run": (
        "{map}",
        [],
        # metacooler's reason (see tests.locks.refusal).
        "{refusal} Or write this run's map elsewhere with --output-map.",
    ),
}


# The refusals of a validation alone: partition reads the contig data from the
# map, and takes no table that the map could not be the map of.
VALIDATION_ONLY = {"a map of other contigs", "a map of other lengths", "an input as --output-map"}


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
    with metacooler's reason, which names its entries."""
    from contextlib import ExitStack

    names = {
        "command": command,
        "map": _contig_map(tmp_path / "input.metacool", network=True),
        "contig_map": _contig_map(tmp_path / "contigs.metacool"),
        "missing": tmp_path / "missing.metacool",
        "directory": tmp_path / "maps",
        "pairs": _pairs(tmp_path),
        "collection": _collection(tmp_path / "maps.mmetacool", "lib1"),
        "copy": tmp_path / "copy.metacool",
        "colon_map": tmp_path / "run::1.metacool",
        "other_map": _contig_map(tmp_path / "other.metacool", {"c1": 100, "c3": 50}),
        "long_map": _contig_map(tmp_path / "long.metacool", {"c1": 120, "c2": 50}),
        "contigs": tmp_path / "contig_data.txt",
    }
    names["contigs"].write_text("ID\tName\tSize\tHit\tShotgun_coverage\tOverlapping_cluster_ID\n1\tc1\t100\t3\t-\t00001\n"
                                "2\tc2\t50\t1\t-\t00001\n")  # fmt: skip
    if case == "a network with text in a column":
        names["text_map"] = _contig_map(tmp_path / "text.metacool", network=True)
        with h5py.File(names["text_map"], "r+") as group:
            del group["contigs/GC_content"]
            group["contigs/GC_content"] = np.array([b"high", b"low"])
    if case == "a network of other contacts":
        names["merged_map"] = _merged_map(tmp_path, "None")
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
    assert not (tmp_path / "copy.metacool").exists()


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_step_on_a_map_refuses_an_existing_output_map(command, stub_pipeline, tmp_path, run_cli):
    """--output-map replaces a file only with -F, as in pipeline: refused
    before any work, the file left as it was."""
    metacool = _contig_map(tmp_path / "input.metacool", network=True)
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
    metacool = _contig_map(tmp_path / "input.metacool", network=True)
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
    assert locks == [[f"{written.name}.lock"]] and not list(tmp_path.glob("*.lock"))
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
    metacool = _contig_map(tmp_path / "input.metacool", network=True)
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
    collection = _collection(tmp_path / "maps.mmetacool", "lib1", "lib2", network=True)
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
    """The map given through a symbolic link: the map it points to stores
    the overlapping clusters, as metacooler edits the target of a link, and
    the link keeps pointing at it, as pipeline updates it."""
    import metacooler

    _partition_table(monkeypatch)
    target = _contig_map(tmp_path / "target.metacool", network=True)
    link = tmp_path / "link.metacool"
    link.symlink_to(target)
    run_cli("partition", *_step_run("partition", link, tmp_path, "-s", "100"))
    assert link.is_symlink() and link.resolve() == target
    assert list(metacooler.MetaCooler(str(target)).clusters.name) == ["metator_00001_00000"]


@pytest.mark.parametrize("error", [ValueError("not valid"), KeyboardInterrupt()], ids=["ValueError", "Ctrl-C"])
def test_partition_on_a_map_that_cannot_be_written(error, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map with the overlapping clusters that cannot be written (here, its
    validation fails) stops the run with an error that says why and that the
    map is left as it was, as pipeline says it; Ctrl-C stops it too. The map
    is left byte for byte, with no file beside it (no staged map, no lock);
    the other outputs of the run are kept."""
    _partition_table(monkeypatch)
    (tmp_path / "maps").mkdir()
    metacool = _contig_map(tmp_path / "maps" / "input.metacool", network=True)
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise error

    # The full validation of the staged map, before metacooler renames it.
    monkeypatch.setattr("metacooler.schema.validate", validate)
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
    metacool = _contig_map(tmp_path / "input.metacool", network=True)
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
    if output == "in place":
        assert metacool.stat().st_mode & 0o777 == 0o640
    assert not list(tmp_path.glob("*.lock")) and not [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]


def test_validation_on_a_map_that_cannot_be_written(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the clusters that cannot be written (here, its validation
    fails) stops validation MAP as it stops pipeline MAP: the map is left
    as it was."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    metacool = _contig_map(tmp_path / "input.metacool", network=True)
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise ValueError("not valid")

    # The full validation of the staged map, before metacooler renames it.
    monkeypatch.setattr("metacooler.schema.validate", validate)
    with pytest.raises(mtmc.MetacoolError, match=r"is not written: .* It is still the input map of the run, as it was\."):
        run_cli("validation", *_step_run("validation", metacool, tmp_path))
    assert metacool.read_bytes() == before


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_a_step_on_a_map_holds_its_lock(command, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """partition MAP and validation MAP, which write a map, hold its lock
    (metacooler's) through their step, and set no signal handler: a run
    killed by SIGTERM dies of it, and leaves its lock to the next writer
    (see hold_lock). A step that fails releases the lock."""
    import signal

    seen = {}
    metacool = _contig_map(tmp_path / "input.metacool", network=True)

    def step(*args, **kwargs):
        lock = lock_file(metacool)
        seen["held"] = lock.exists() and json.loads(lock.read_text())["pid"] == os.getpid()
        seen["handler"] = signal.getsignal(signal.SIGTERM)
        raise RuntimeError("stopped")

    target = {"partition": "metator.partition.partition", "validation": "metator.validation.recursive_decontamination"}
    monkeypatch.setattr(target[command], step)
    with pytest.raises(RuntimeError, match="stopped"):
        run_cli(command, *_step_run(command, metacool, tmp_path))
    assert seen == {"held": True, "handler": signal.SIG_DFL}
    assert not list(tmp_path.glob("*.lock"))


@pytest.mark.parametrize("command", ["partition", "validation"])
def test_a_step_refuses_a_bad_prefix_before_replacing_the_clusters(command, stub_pipeline, tmp_path, run_cli):
    """An invalid --prefix is refused before -F replaces the directories of
    clusters of an earlier run, as pipeline refuses it: it emptied them
    first, then was refused."""
    earlier = _earlier_run(command, tmp_path)
    metacool = _contig_map(tmp_path / "input.metacool", network=True)
    with pytest.raises(ValueError, match="Invalid --prefix"):
        run_cli(command, *_step_run(command, metacool, tmp_path, "-F", "-p", "bad.prefix"))
    assert all(path.exists() for path in earlier) and _no_step(stub_pipeline)


def test_pipeline_replaces_an_earlier_map_named_after_its_assembly(stub_pipeline, tmp_path, run_cli):
    """pipeline --output-map -F replaces the map of an earlier run named
    after the assembly (assembly.fa.metacool, for -a assembly.fa): the
    assembly itself is the input that --output-map must not name."""
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


@pytest.mark.parametrize("name", ["final_clusters", "recursive_clusters/step_1"])
def test_validation_refuses_its_input_in_what_force_removes(name, stub_pipeline, tmp_path, run_cli, logged):
    """The clusters to validate (-f) in a directory that -F removes: refused
    before anything is removed. -F removed them, then the run failed, or
    found no cluster."""
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


@pytest.mark.parametrize("dangling", [False, True], ids=["link", "dangling link"])
def test_force_replaces_a_link_of_clusters_not_its_target(dangling, stub_pipeline, tmp_path, run_cli):
    """A directory of clusters of an earlier run that is a symbolic link, to
    a directory elsewhere or to nothing, is refused without -F, and replaced
    with -F: the link is removed, its target kept. -F failed with OSError
    (rmtree on a link), after the other directory was replaced."""
    target = tmp_path / "elsewhere"
    if not dangling:
        target.mkdir()
        (target / "earlier.fa").write_text(">c1\nACGT\n")
    link = tmp_path / "out" / "final_clusters"
    link.parent.mkdir()
    link.symlink_to(target)
    with pytest.raises(ValueError):
        run_cli("validation", *_step_run("validation", None, tmp_path))
    assert os.path.lexists(link) and _no_step(stub_pipeline)
    run_cli("validation", *_step_run("validation", None, tmp_path, "-F"))
    assert link.is_dir() and not link.is_symlink() and stub_pipeline
    assert dangling or (target / "earlier.fa").exists()


def test_validation_refuses_its_input_through_a_link_that_force_replaces(stub_pipeline, tmp_path, run_cli, logged):
    """The clusters to validate (-f) in the target of a directory of
    clusters that is a symbolic link: refused with -F, before anything is
    replaced, as in the directory itself. The run let it through, as the
    link alone is removed, then replaced the directory it was in."""
    target = tmp_path / "elsewhere"
    fasta = target / "clusters"
    fasta.mkdir(parents=True)
    (fasta / "metator_00001_00000.fa").write_text(">c1\nACGT\n")
    link = tmp_path / "out" / "final_clusters"
    link.parent.mkdir()
    link.symlink_to(target)
    args = _step_run("validation", None, tmp_path, "-F")
    args[args.index("overlapping_clusters")] = fasta
    with pytest.raises(ValueError):
        run_cli("validation", *args)
    message = f"{fasta}, given as -f, is in {link}, which -F removes: write the run to another directory (-o)."
    assert ("ERROR", message) in logged()
    assert link.is_symlink() and (fasta / "metator_00001_00000.fa").exists() and not stub_pipeline


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
    validation from the tables of the network (-n, and partition's -c), -n
    an unknown option of validation, -c of partition, as both read the
    network that metator network stores in the map."""
    name, *args = command.split()
    ctx = click.Context(mtcli.cli, info_name="metator")
    option = {"network": "-1", "pipeline": "-1", "partition": "-c", "validation": "-n"}[name]
    with pytest.raises(click.NoSuchOption, match=rf"^No such option '{option}'\.$"):
        mtcli.cli.get_command(ctx, name).make_context(name, list(args), parent=ctx)


def test_an_unknown_option_is_still_unknown(run_cli):
    """A guard: an option that metaTOR never had is refused as unknown, as
    before, by click's own message."""
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
    """A refused run with --force (a rebinning into the same directory)
    stops before replacing the bins of the earlier run."""
    names = {
        "map": _contig_map(tmp_path / "input.metacool"),
        "missing": tmp_path / "missing.metacool",
        "collection": _collection(tmp_path / "maps.mmetacool", "lib1"),
    }
    args = [arg.format(**names) for arg in args]
    out = tmp_path / "out"
    # The outputs of an earlier run.
    names = ["overlapping_clusters", "recursive_clusters", "final_clusters"]
    bins = [out / name / "bin.fa" for name in names]
    for path in bins:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    with pytest.raises(ValueError):
        run_cli("pipeline", "-a", TESTS_DATA / "assembly.fa", *args, "-o", out, "-F", "-T", tmp_path / "tmp")
    assert all(path.exists() for path in bins)


def test_pipeline_refuses_earlier_bins_before_reading_the_assembly(tmp_path, monkeypatch, run_cli):
    """Without --force, the clusters of an earlier run are refused before
    the assembly is read (minutes on a large metagenome)."""

    def create_contig_data(*args):
        raise AssertionError("assembly read before the refusal")

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    out = tmp_path / "out"
    (out / "recursive_clusters").mkdir(parents=True)
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    with pytest.raises(ValueError):
        run_cli("pipeline", _contig_map(tmp_path / "map.metacool"), "-a", assembly, "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"ERROR :: {out / 'recursive_clusters'} already existed." in log.read_text()


def test_pipeline_without_pyfastx_index(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """No bin passed the size threshold, so pyfastx never indexed the
    assembly: the run ends without removing an index that is not there (it
    used to fail at its very end, with FileNotFoundError)."""
    monkeypatch.setattr("metator.validation.recursive_decontamination", lambda *args, **kwargs: None)
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    out = tmp_path / "out"
    run_cli("pipeline", _contig_map(tmp_path / "map.metacool"), "-a", assembly, "-o", out, "-T", tmp_path / "tmp")
    assert stub_pipeline["contig_contacts"] and not Path(f"{assembly}.fxi").exists()


@pytest.mark.parametrize("random", [False, True], ids=["mge", "--random"])
def test_mge_default_outdir(random, tmp_path, monkeypatch, run_cli):
    """Without --outdir, mge writes its MGE fasta to the current directory
    (mges_binned_random.fa with --random), with the assembly that the map
    was made from, and stores the round mge (mge_random) in the map; its
    temporary directory is removed at its end, and nothing is written next
    to the assembly."""
    import metacooler
    import pandas as pd

    import metator.clusters as mtcl

    metacool = _contig_map(tmp_path / "map.metacool")
    calls = {}
    summary = pd.DataFrame({"BinName": ["MetaTOR_MGE_00001"], "Contigs": ["c2"]})

    def mge_binning(metacool, mges_list, threshold_bin, random):
        calls["binning"] = (metacool, mges_list, threshold_bin, random)
        return summary, {1: {"Contigs": ["c2"], "Score": float("nan")}}

    def write_mges_fasta(metacool, assembly, mge_bins, out_file, temp_dir):
        calls["fasta"] = (metacool, assembly, out_file)
        assert os.path.isdir(temp_dir)

    def mge_metacool(metacool, outfile, stored, name, planned):
        calls["map"] = (metacool, outfile, stored, name)
        assert isinstance(planned, metacooler.MetaCooler)

    monkeypatch.setattr("metator.mge.mge_binning", mge_binning)
    monkeypatch.setattr("metator.mge.write_mges_fasta", write_mges_fasta)
    monkeypatch.setattr("metator.cli.mge_metacool", mge_metacool)
    assembly = tmp_path / "inputs" / "assembly.fa"
    assembly.parent.mkdir()
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    mges = tmp_path / "mges.txt"
    mges.write_text("c2\n")
    run_cli("mge", metacool, "-a", assembly, "-m", mges, *(["--random"] if random else []))
    assert calls["binning"] == (str(metacool), ["c2"], 0.8, random)
    fasta = "mges_binned_random.fa" if random else "mges_binned.fa"
    assert calls["fasta"] == (str(metacool), str(assembly), os.path.join(".", fasta))
    assert calls["map"] == (str(metacool), str(metacool), summary, mtcl.MGE_RANDOM if random else mtcl.MGE)
    assert sorted(path.name for path in assembly.parent.iterdir()) == ["assembly.fa"]
    assert not (tmp_path / "tmp").exists()
