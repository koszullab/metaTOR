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
import pytest

import metator.cli as mtcli
import metator.metacool as mtmc
from metator.version import __version__

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"

# The writers themselves, which the stub pipeline replaces.
mtmc_write_contig_map = mtmc.write_contig_map
mtmc_reassign_metacool = mtmc.reassign_metacool

# The docopt command line that rich-click replaced (that of 1.3.10, with the
# mge options renamed and host removed), as docopt read it: see the commit
# that added the file. Commands and options added since are added to it:
# metacool, and the --no-metacool option of pipeline; and the --forward of
# pipeline is required, as every start needs it.
INVENTORY = json.loads((TESTS_DATA / "cli_inventory.json").read_text())["metator"]

# Options whose value is one of a list, wherever they appear. --object takes
# any column name of the contig data file. Values added since docopt:
# metacool (--start of network and pipeline).
CHOICES = {
    "--aligner": ["bowtie2", "bwa"],
    "--aligner-mode": ["normal", "iterative", "cutsite"],
    "--algorithm": ["louvain", "leiden"],
    "--normalization": ["None", "abundance", "length", "RS", "empirical_hit", "theoritical_hit"],
    "--start": ["fastq", "bam", "pair", "metacool"],
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


def test_network_help_of_a_metacool_start(run_cli):
    """--edge must be 0 with --start metacool, as a map holds no read
    positions: the help says so right after what a map gives."""
    output = " ".join(run_cli("network", "--help").output.split())
    assert "that of the pairs it was made of. --edge must then be 0." in output


def test_group_help(run_cli):
    output = run_cli("-h").output
    for command in INVENTORY["commands"]:
        assert command in output
    assert "Binning" in output


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


@pytest.mark.parametrize("command", [command for command, frozen in INVENTORY["commands"].items() if frozen["positional"]])
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

    def get_contact_pairs(*args):
        calls["get_contact_pairs"] = args
        out_dir = args[10]
        # Sorted, compressed and indexed by pypairix, as the alignment step
        # leaves them.
        pairs = [os.path.join(out_dir, f"alignment_{i}_sorted.pairs.gz") for i in range(len(args[0].split(",")))]
        for path in pairs:
            Path(path).write_text("## pairs format v1.0\n")
            Path(f"{path}.px2").write_text("x")
        return pairs, {}, {}

    def create_contig_data(*args):
        calls["create_contig_data"] = args
        return {}, {}

    def alignment_to_contacts(alignment_files, contig_data, edge, hit_data, out_dir, network, contig_data_file, *args):
        calls["alignment_to_contacts"] = (alignment_files, edge, out_dir, *args)
        for name in [network, contig_data_file]:
            Path(out_dir, name).write_text("x")
        return os.path.join(out_dir, network), os.path.join(out_dir, contig_data_file)

    def metacool_to_contacts(metacool, contig_data, out_dir, network, contig_data_file, *args):
        calls["metacool_to_contacts"] = (metacool, out_dir, network, contig_data_file, *args)
        for name in [network, contig_data_file]:
            Path(out_dir, name).write_text("x")
        return os.path.join(out_dir, network), os.path.join(out_dir, contig_data_file)

    def partition(*args):
        calls["partition"] = args
        out_dir = args[6]
        # What the network step left in the output directory.
        calls["partition_outdir"] = sorted(os.listdir(out_dir))
        Path(out_dir, "contig_data_partition.txt").write_text("x")
        return os.path.join(out_dir, "clustering_matrix_partition"), os.path.join(out_dir, "contig_data_partition.txt")

    def recursive_decontamination(*args):
        calls["recursive_decontamination"] = args
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

    def reassign_metacool(metacool, contig_data, outfile, **kwargs):
        # The input map must still be there.
        calls["reassign_metacool"] = (metacool, os.path.exists(metacool), contig_data, outfile, kwargs)

    monkeypatch.setattr("metator.align.get_contact_pairs", get_contact_pairs)
    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    monkeypatch.setattr("metator.network.alignment_to_contacts", alignment_to_contacts)
    monkeypatch.setattr("metator.network.metacool_to_contacts", metacool_to_contacts)
    monkeypatch.setattr("metator.partition.partition", partition)
    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)
    monkeypatch.setattr("metator.io.generate_fasta_index", generate_fasta_index)
    monkeypatch.setattr("metator.metacool.write_metacool", write_metacool)
    monkeypatch.setattr("metator.metacool.write_contig_map", write_contig_map)
    monkeypatch.setattr("metator.metacool.reassign_metacool", reassign_metacool)
    return calls


@pytest.mark.parametrize("start", ["fastq", "bam"])
@pytest.mark.parametrize("clean_up", [True, False])
def test_pipeline_from_reads(start, clean_up, stub_pipeline, tmp_path, run_cli):
    """From reads or alignments: the pairs files and the contig data of the
    network step are intermediate files."""
    if start == "fastq":
        # A fasta: its index is built.
        assembly = tmp_path / "assembly.fa"
        shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    else:
        # A bowtie2 index: its fasta is retrieved.
        assembly = TESTS_DATA / "assembly"
    out = tmp_path / "out"
    args = ["-1", "a1.fq,b1.fq", "-2", "a2.fq,b2.fq", "-a", assembly, "-S", start, "-B", "cutsite", "-q", "20",
            "-t", "2", "-o", out, "-T", tmp_path / "tmp"]  # fmt: skip
    if not clean_up:
        args.append("-N")
    run_cli("pipeline", *args)

    calls = stub_pipeline
    forward, reverse, index, fasta, aligner, mode, min_quality, start_arg, *_, threads = calls["get_contact_pairs"]
    assert (forward, reverse, aligner, start_arg) == ("a1.fq,b1.fq", "a2.fq,b2.fq", "bowtie2", start)
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
    assert not Path(f"{fasta}.fxi").exists()
    assert not (out / "contig_data_network.txt").exists()
    assert not (out / "contig_data_partition.txt").exists()
    assert (out / "network.txt").exists()
    # The pairs files and their pypairix index go together.
    pairs = [out / "alignment_0_sorted.pairs.gz", out / "alignment_1_sorted.pairs.gz"]
    assert [path.exists() for path in pairs] == [not clean_up] * 2
    assert [Path(f"{path}.px2").exists() for path in pairs] == [not clean_up] * 2
    assert (tmp_path / "tmp").exists() == (not clean_up)
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"## Command used: metator pipeline -1 a1.fq,b1.fq -2 a2.fq,b2.fq -a {assembly}" in log.read_text()

    # The contact map of the contigs, built from the pairs files of the run
    # as soon as they are there, before the network.
    order = list(calls)
    assert order.index("get_contact_pairs") < order.index("write_contig_map") < order.index("alignment_to_contacts")
    pairs_files, existed, contig_data, outfile, kwargs = calls["write_contig_map"]
    assert pairs_files == [str(path) for path in pairs] and existed == [True, True]
    assert (contig_data, outfile) == ({}, str(out / "metator.metacool"))
    assert (kwargs["assembly_name"], kwargs["threads"], kwargs["force"]) == ("assembly", 2, True)
    assert Path(kwargs["tmp_dir"]).parent == tmp_path / "tmp"
    metadata = kwargs["metadata"]
    assert metadata["command"].startswith("metator pipeline -1 a1.fq,b1.fq -2 a2.fq,b2.fq -a")
    assert metadata["parameters"]["start"] == start and metadata["parameters"]["no_metacool"] is False
    assert "metator.metacool" in calls["partition_outdir"]
    # The map of the bins, written last over it, from the final tables: the
    # pairs files are not read again.
    assert order[-1] == "reassign_metacool" and "write_metacool" not in calls
    source, existed, contig_data, outfile, kwargs = calls["reassign_metacool"]
    assert (source, existed, outfile) == (str(out / "metator.metacool"), True, str(out / "metator.metacool"))
    assert contig_data == str(out / "contig_data_final.txt")
    assert (kwargs["bin_summary"], kwargs["force"], kwargs["metadata"]) == (str(out / "bin_summary.txt"), True, metadata)
    # The bin summary goes to metacooler as it is: no temporary copy.
    assert "tmp_dir" not in kwargs


@pytest.mark.parametrize(
    "options, written",
    [([], "metator.metacool"), (["-p", "sample_1"], "sample_1.metacool"), (["--no-metacool"], "metator.metacool")],
)
def test_pipeline_metacool_name(options, written, stub_pipeline, tmp_path, run_cli):
    """<prefix>.metacool, metator.metacool without --prefix: the map of the
    contigs is built there, and the map of the bins written over it at the
    end. With --no-metacool, the map of the contigs is built all the same,
    and removed at the end."""
    out = tmp_path / "out"
    run_cli("pipeline", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp", *options)
    outfile = str(out / written)
    pairs_files, _, _, built, _ = stub_pipeline["write_contig_map"]
    assert pairs_files == ["a.pairs"] and built == outfile
    assert written in stub_pipeline["partition_outdir"]
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    if "--no-metacool" in options:
        assert "reassign_metacool" not in stub_pipeline and not (out / written).exists()
        assert f"INFO :: Removing {outfile}, the contact map of the contigs (--no-metacool)." in logged
    else:
        source, _, _, replaced, _ = stub_pipeline["reassign_metacool"]
        assert source == replaced == outfile


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

    def recursive_decontamination(*args):
        out_dir = args[8]
        Path(out_dir, "contig_data_final.txt").write_text(
            "\tID\tName\tSize\tGC_content\tHit\tFinal_bin\n"
            f"0\t1\t{contig}\t100\t40.5\t-\t{final_bin}\n1\t2\tc2\t50\t35.0\t-\tND\n"
        )
        rows = "" if final_bin == "ND" else f"{final_bin}\t0.9\n"
        Path(out_dir, "bin_summary.txt").write_text(summary or f"\tCompleteness\n{rows}")
        return stub(*args)

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
def test_pipeline_non_ascii_contigs_skip_the_metacool(no_metacool, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """Contig names that are not ASCII: the map of the contigs cannot be
    built. A warning, or, with --no-metacool, which asks for no map, a line
    at INFO; no .metacool, no map of the bins at the end, and the run
    succeeds with its other outputs."""
    monkeypatch.setattr("metator.metacool.write_contig_map", mtmc_write_contig_map)
    _contig_data(monkeypatch, {"contig_é": (100, 40.5), "c2": (50, 35.0)})
    _final_tables(monkeypatch, "contig_é")
    out = tmp_path / "out"
    args = ["-1", _pairs(tmp_path, "contig_é"), "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp"]
    run_cli("pipeline", *args, *(["--no-metacool"] if no_metacool else []))
    assert not list(out.glob("*.metacool*"))
    assert "reassign_metacool" not in stub_pipeline
    assert (out / "network.txt").exists() and (out / "contig_data_final.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    reason = (
        f"{out / 'metator.metacool'} is not written: A .metacool file needs ASCII contig names (cooler reads them as "
        "ASCII text), and 1 contig name(s) are not: 'contig_é'. The run goes on without it"
    )
    if no_metacool:
        assert f"INFO :: {reason}.\n" in logged and "WARNING" not in logged
    else:
        assert f"WARNING :: {reason}, and writes no map of the bins.\n" in logged
    assert "Traceback" not in logged


@pytest.mark.parametrize("no_metacool", [False, True])
def test_pipeline_metacool_failure_is_a_warning(no_metacool, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the contigs that fails to build costs a warning, with the
    traceback in the log, --no-metacool or not when the error is
    unexpected, not the run: the network is built, no map of the bins is
    written, the exit status is 0 and the pairs files go as usual."""

    def write_contig_map(*args, **kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr("metator.metacool.write_contig_map", write_contig_map)
    out = tmp_path / "out"
    args = ["-1", "a.fq", "-2", "b.fq", "-a", TESTS_DATA / "assembly", "-S", "bam", "-o", out, "-T", tmp_path / "tmp"]
    run_cli("pipeline", *args, *(["--no-metacool"] if no_metacool else []))
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    consequence = "The run goes on without it." if no_metacool else "The run goes on without it, and writes no map of the bins."
    warning = f"WARNING :: {out / 'metator.metacool'} is not written: RuntimeError: disk full. {consequence}"
    assert f"{warning}\nTraceback" in logged
    assert logged.count("WARNING ::") == 1
    assert "reassign_metacool" not in stub_pipeline
    assert (out / "network.txt").exists() and not (out / "metator.metacool").exists()
    assert not (out / "alignment_0_sorted.pairs.gz").exists()
    assert not (out / "alignment_0_sorted.pairs.gz.px2").exists()


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
    pairs file, which --start metacool reads. The map is left as it was."""
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
    assert not stub_pipeline and not (out / "overlapping_bin").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    if case in ["map elsewhere", "entry of a collection"]:
        message = f"{forward} is a contact map, not a pairs file: start from it with --start metacool."
    else:
        message = (
            f"{refused} is the contact map of this run ({outfile}), which the run removes before it builds the map of "
            "its contigs there: start from it with --start metacool, or write the run to another output directory or "
            "under another --prefix."
        )
    assert f"ERROR :: {message}" in log.read_text()


def test_pipeline_from_a_pairs_file_named_with_colons(stub_pipeline, tmp_path, run_cli):
    """A pairs file whose name holds :: is a pairs file, not the URI of a
    map inside a file: the run reads it, and builds its map from it."""
    pairs = tmp_path / "lib::L1.pairs"
    pairs.write_text("## pairs format v1.0\n")
    out = tmp_path / "out"
    run_cli("pipeline", "-1", pairs, "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    assert stub_pipeline["write_contig_map"][0] == [str(pairs)]
    assert stub_pipeline["alignment_to_contacts"][0] == [str(pairs)]
    assert "reassign_metacool" in stub_pipeline


def test_pipeline_interrupted_while_building_the_map(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """An interruption (Ctrl-C) while the map of the contigs is built stops
    the run: it is not one of the failures of the map that the run goes on
    without."""

    def write_contig_map(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr("metator.metacool.write_contig_map", write_contig_map)
    out = tmp_path / "out"
    result = run_cli("pipeline", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp",
                     exit_code=1)  # fmt: skip
    assert "Aborted" in result.output
    assert "alignment_to_contacts" not in stub_pipeline and "partition" not in stub_pipeline
    (log,) = out.glob("metator_pipeline_*.log")
    assert "not written" not in log.read_text()


@pytest.mark.parametrize("case", ["written", "failed", "no metacool", "link"])
def test_pipeline_earlier_metacool(case, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map left by an earlier run in the output directory is removed
    before the map of the contigs is built in its place, --no-metacool or
    not: it would pass for this run's map when the new one cannot be built.
    From then on, the map there is this run's, or there is none. A symbolic
    link there is removed, not its target."""
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
    run_cli("pipeline", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp", *options)
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert seen == [False]
    assert f"INFO :: Removing {earlier}, the contact map of an earlier run." in logged
    if case == "failed":
        assert not os.path.lexists(earlier) and "reassign_metacool" not in stub_pipeline
        assert f"WARNING :: {earlier} is not written: RuntimeError: disk full." in logged
    elif case == "no metacool":
        assert not os.path.lexists(earlier) and "reassign_metacool" not in stub_pipeline
    else:
        assert not earlier.is_symlink() and earlier.read_text() == "this run"
        assert stub_pipeline["reassign_metacool"][0] == str(earlier)
    if case == "link":
        assert target.read_text() == "kept"


def test_pipeline_from_a_map_earlier_metacool(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """From a map stored elsewhere: a map left by an earlier run where the
    map of the bins goes is removed before it is written, as it no longer
    matches the bins."""
    out = tmp_path / "out"
    out.mkdir()
    earlier = out / "metator.metacool"
    earlier.write_text("earlier run")
    metacool = tmp_path / "input.metacool"
    metacool.write_text("input map")
    seen = []

    def reassign_metacool(metacool, contig_data, outfile, **kwargs):
        seen.append(os.path.lexists(outfile))

    monkeypatch.setattr("metator.metacool.reassign_metacool", reassign_metacool)
    run_cli("pipeline", "-1", metacool, "-a", _assembly(tmp_path), "-S", "metacool", "-o", out, "-T", tmp_path / "tmp")
    assert seen == [False] and not earlier.exists() and metacool.read_text() == "input map"
    assert "write_contig_map" not in stub_pipeline
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"INFO :: Removing {earlier}, the contact map of an earlier run." in log.read_text()


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
    # every contig), not its ID or Final_bin.
    contigs = metacooler.MetaCooler(str(out / "metator.metacool")).contigs.set_index("name")
    assert [column for column in contigs.columns if not column.startswith("qc_")] == ["length", "genome_id", "GC_content"]
    assert contigs.GC_content.to_dict() == {"c1": 40.5, "c2": 35.0}
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]


def test_pipeline_map_of_the_bins_failure_keeps_the_map_of_the_contigs(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the bins that cannot be written over the map of the contigs
    (here, a contig of another length in the final contig data): a warning
    that says why and that the map of the contigs stays, as it was, in no
    bin, and the run keeps its other outputs."""
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
    run_cli("pipeline", *_pair_run(tmp_path, out))
    assert outfile.read_bytes() == built[0]
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    assert (out / "network.txt").exists() and (out / "bin_summary.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert (
        f"WARNING :: {outfile} is not written: {outfile} is not a contact map of the contigs of "
        f"{out / 'contig_data_final.txt'}: 1 contig(s) have another length in it: 'c1' (120 bp in the map, 100 bp in the "
        "contig data). It is still the contact map of the contigs of the run, in no bin. The other outputs of the run "
        "are kept."
    ) in logged
    assert "Traceback" not in logged


@pytest.mark.parametrize("no_metacool", [False, True])
def test_pipeline_failure_leaves_the_map_of_the_contigs(no_metacool, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A run that stops with an error after it built the map of the contigs
    (in the validation, here) leaves it at the outfile, in no bin,
    --no-metacool or not, as it leaves its other files: a new run can start
    from it (--start metacool)."""
    monkeypatch.setattr("metator.metacool.write_contig_map", mtmc_write_contig_map)
    _contig_data(monkeypatch, CONTIGS)

    def recursive_decontamination(*args):
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


def test_pipeline_from_pairs_keeps_them(stub_pipeline, tmp_path, run_cli):
    """From pairs files: they are inputs, kept with their pypairix index."""
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    pairs = [tmp_path / "a_sorted.pairs.gz", tmp_path / "b_sorted.pairs.gz"]
    for path in pairs:
        path.write_text("## pairs format v1.0\n")
        Path(f"{path}.px2").write_text("x")
    out = tmp_path / "out"
    run_cli("pipeline", "-1", ",".join(map(str, pairs)), "-a", assembly, "-S", "pair", "-o", out, "-T", tmp_path / "tmp")

    calls = stub_pipeline
    assert "get_contact_pairs" not in calls
    assert calls["alignment_to_contacts"][0] == [str(path) for path in pairs]
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
    metacool = out / "metator.metacool"
    metacool.write_text("input map")
    options = ["--no-metacool"] if no_metacool else []
    run_cli("pipeline", "-1", metacool, "-a", _assembly(tmp_path), "-S", "metacool", "-t", "3", "-m", "abundance",
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
    assert metacool.read_text() == "input map"
    # The contig data of the network step is kept, as from pairs files.
    assert (out / "contig_data_network.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert "Removing" not in logged
    if no_metacool:
        assert "reassign_metacool" not in calls and "is the input map" not in logged
        return
    source, existed, contig_data, outfile, kwargs = calls["reassign_metacool"]
    assert (source, existed, contig_data, outfile) == (str(metacool), True, str(out / "contig_data_final.txt"), str(metacool))
    assert (kwargs["bin_summary"], kwargs["force"]) == (str(out / "bin_summary.txt"), True)
    assert kwargs["metadata"]["parameters"]["start"] == "metacool"
    assert kwargs["metadata"]["command"].startswith(f"metator pipeline -1 {metacool} -a")
    assert f"INFO :: {metacool} is the input map {metacool}: it is replaced once the map of the bins is validated." in logged


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


@pytest.mark.parametrize("where", ["elsewhere", "outfile"])
def test_pipeline_metacool_written_from_a_map(where, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """With the real writer, from the tables of the (stub) validation: the
    map of the bins, from the input map, which it replaces when it is the
    outfile."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map((out if where == "outfile" else tmp_path) / "metator.metacool")
    before = metacool.read_bytes()
    run_cli("pipeline", "-1", metacool, "-a", _assembly(tmp_path), "-S", "metacool", "-o", out, "-T", tmp_path / "tmp")
    import cooler
    import metacooler

    outfile = out / "metator.metacool"
    metacooler.validate(str(outfile), full=True)
    assert cooler.Cooler(str(outfile)).pixels()[:]["count"].sum() == 2
    mc = metacooler.MetaCooler(str(outfile))
    assert list(mc.genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert list(mc.genomes.Completeness.iloc[1:]) == [0.9]
    assert dict(zip(mc.contigs.name, mc.contigs.GC_content, strict=True)) == {"c1": 40.5, "c2": 35.0}
    assert (metacool.read_bytes() == before) == (where == "elsewhere")
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    (log,) = out.glob("metator_pipeline_*.log")
    assert "is not written" not in log.read_text()


def test_pipeline_metacool_from_a_map_failure_keeps_it(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map of the bins that cannot be written from the input map, which is
    the outfile (here, a contig of another length in the contig data): a
    warning that says why and that it is still the input map, left as it
    was, and the run keeps its other outputs."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool", {"c1": 120, "c2": 50})
    before = metacool.read_bytes()
    run_cli("pipeline", "-1", metacool, "-a", _assembly(tmp_path), "-S", "metacool", "-o", out, "-T", tmp_path / "tmp")
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    assert (out / "network.txt").exists() and (out / "bin_summary.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert (
        f"WARNING :: {metacool} is not written: {metacool} is not a contact map of the contigs of "
        f"{out / 'contig_data_final.txt'}: 1 contig(s) have another length in it: 'c1' (120 bp in the map, 100 bp in the "
        "contig data). It is still the input map of the run. The other outputs of the run are kept."
    ) in logged
    assert "Traceback" not in logged


@pytest.mark.parametrize(
    "error", [ValueError("pixels out of order"), RuntimeError("unexpected failure")], ids=["ValueError", "RuntimeError"]
)
def test_pipeline_metacool_from_a_map_metacooler_failure(error, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """metacooler fails while it writes the map of the bins over the input
    map (its validation, here): it removes its partial file and says so at
    INFO, and the run warns that the map is not written, with the reason of
    one of the errors that metacooler is expected to raise (a ValueError
    here), or the traceback of another (a RuntimeError). The input map is
    left as it was, and no other file is left beside it."""
    monkeypatch.setattr("metator.metacool.reassign_metacool", mtmc_reassign_metacool)
    _final_tables(monkeypatch, "c1")
    out = tmp_path / "out"
    out.mkdir()
    metacool = _contig_map(out / "metator.metacool")
    before = metacool.read_bytes()

    def validate(*args, **kwargs):
        raise error

    monkeypatch.setattr("metacooler.io.reassign.validate", validate)
    run_cli("pipeline", "-1", metacool, "-a", _assembly(tmp_path), "-S", "metacool", "-o", out, "-T", tmp_path / "tmp")
    assert metacool.read_bytes() == before
    assert sorted(path.name for path in out.iterdir() if "metacool" in path.name) == ["metator.metacool"]
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    (staged,) = re.findall(r"INFO :: Creation failed; the partial file (\S+)\.partial was removed and \1 was not written", logged)
    assert Path(staged).parent == out and Path(staged).name.startswith(".metator.metacool.")
    assert "for inspection" not in logged and "ERROR ::" not in logged
    kept = "It is still the input map of the run. The other outputs of the run are kept."
    if isinstance(error, ValueError):
        assert f"WARNING :: {metacool} is not written: pixels out of order. {kept}\n" in logged
        assert "Traceback" not in logged
    else:
        assert f"WARNING :: {metacool} is not written: RuntimeError: unexpected failure. {kept}\nTraceback" in logged
    assert logged.count("WARNING ::") == 1


def test_network_from_a_metacool(monkeypatch, tmp_path, run_cli):
    """metator network --start metacool reads the network from the map, with
    the contig data of the assembly."""
    calls = {}

    def create_contig_data(*args):
        calls["create_contig_data"] = args
        return {"c1": {}}, None

    def metacool_to_contacts(*args):
        calls["metacool_to_contacts"] = args

    def not_called(*args):
        raise AssertionError("no alignment nor pairs from a metacool")

    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    monkeypatch.setattr("metator.network.metacool_to_contacts", metacool_to_contacts)
    monkeypatch.setattr("metator.network.alignment_to_contacts", not_called)
    monkeypatch.setattr("metator.align.get_contact_pairs", not_called)
    out = tmp_path / "out"
    run_cli("network", "-1", "maps.mmetacool::entries/lib1", "-a", _assembly(tmp_path), "-S", "metacool", "-e", "DpnII",
            "-s", "-t", "2", "-n", "RS", "-o", out, "-T", tmp_path / "tmp")  # fmt: skip
    assert calls["create_contig_data"] == (str(tmp_path / "assembly.fa"), 1, None, "DpnII")
    uri, contig_data, *args = calls["metacool_to_contacts"]
    assert uri == "maps.mmetacool::entries/lib1" and contig_data == {"c1": {}}
    out_dir, network, contigs, tmp_dir, threads, normalization, self_contacts = args
    assert (out_dir, network, contigs, threads, normalization, self_contacts) == (
        str(out), "network.txt", "contig_data_network.txt", 2, "RS", True)  # fmt: skip
    assert Path(tmp_dir).parent == tmp_path / "tmp"


METACOOL_START_ERRORS = {
    "several maps": (
        ["-1", "a.metacool,b.metacool,c.metacool"],
        "--start metacool takes one .metacool file (its libraries are already summed), given with -1/--forward: got 3 files.",
    ),
}
# pipeline refuses an edge region whatever the start: see
# test_pipeline_refuses_an_edge_region.
NETWORK_METACOOL_START_ERRORS = {
    "edge region": (
        ["-1", "map.metacool", "-E", "100"],
        "--edge 100 cannot be applied with --start metacool: a contact map holds the contacts of bins, not the "
        "positions of the reads. Leave --edge at 0, or start from the pairs files.",
    ),
}
# A run without a map is refused by click: --forward is required.
PIPELINE_METACOOL_START_ERRORS = {
    "scaffold": (["-1", "map.metacool", "--scaffold"], "--scaffold reads the pairs files: it cannot run with --start metacool."),
}


@pytest.mark.parametrize(
    "command, case",
    [(command, case) for command in ["network", "pipeline"] for case in METACOOL_START_ERRORS]
    + [("network", case) for case in NETWORK_METACOOL_START_ERRORS]
    + [("pipeline", case) for case in PIPELINE_METACOOL_START_ERRORS],
)
def test_metacool_start_refused(command, case, monkeypatch, tmp_path, run_cli):
    """Options a contact map cannot serve are refused before any work, with
    a message in the log file."""

    def not_called(*args):
        raise AssertionError("refused before any work")

    monkeypatch.setattr("metator.network.create_contig_data", not_called)
    args, message = {**METACOOL_START_ERRORS, **NETWORK_METACOOL_START_ERRORS, **PIPELINE_METACOOL_START_ERRORS}[case]
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli(command, *args, "-a", TESTS_DATA / "assembly.fa", "-S", "metacool", "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}" in log.read_text()


PIPELINE_EDGE_ERROR = (
    "--edge {} cannot be applied by pipeline, which counts every pair, as its contact map does: a map holds the "
    "contacts of bins, not the positions of the reads. Leave --edge at 0, or build the network with metator network, "
    "then run metator partition and metator validation."
)


@pytest.mark.parametrize(
    "start, inputs",
    [
        ("fastq", ["-1", "a.fq", "-2", "b.fq"]),
        ("bam", ["-1", "a.bam", "-2", "b.bam"]),
        ("pair", ["-1", "a.pairs"]),
        ("pair", ["-1", "a.pairs,b.pairs"]),
        ("metacool", ["-1", "map.metacool"]),
    ],
    ids=["fastq", "bam", "pair", "pairs of two libraries", "metacool"],
)
@pytest.mark.parametrize("edge", ["100", "-5"], ids=["above 0", "negative"])
def test_pipeline_refuses_an_edge_region(start, inputs, edge, stub_pipeline, tmp_path, run_cli):
    """pipeline refuses an --edge other than 0 before any work, whatever
    the start and the number of libraries: it counts every pair, as its
    contact map does. A negative --edge, which counted every pair as 0
    does, is refused too, as the help says "Only 0". --start metacool gets
    this message too, not the one of metator network, which would send it
    to the pairs files."""
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli("pipeline", *inputs, "-a", _assembly(tmp_path), "-S", start, f"--edge={edge}", "-o", out, "-T", tmp_path / "tmp")
    assert not stub_pipeline and not (out / "overlapping_bin").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"ERROR :: {PIPELINE_EDGE_ERROR.format(edge)}\n" in logged
    assert logged.count("ERROR ::") == 1


@pytest.mark.parametrize("edge", [100, -5], ids=["above 0", "negative"])
def test_network_keeps_the_edge_region(edge, stub_pipeline, tmp_path, run_cli):
    """A guard: metator network still applies --edge to the pairs files
    (alignment_to_contacts), a negative one included."""
    run_cli("network", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", f"--edge={edge}", "-o", tmp_path / "out",
            "-T", tmp_path / "tmp")  # fmt: skip
    alignment_files, edge_arg, *_ = stub_pipeline["alignment_to_contacts"]
    assert (alignment_files, edge_arg) == (["a.pairs"], edge)


@pytest.mark.parametrize(
    "args",
    [
        ["-S", "metacool", "-1", "map.metacool", "-E", "100"],
        ["-S", "pair", "-1", "a.pairs", "-E", "100"],
        ["-S", "pair", "-1", "a.pairs", "--edge=-5"],
        ["-S", "metacool", "-1", "map.metacool", "--scaffold"],
        ["-S", "pair", "-1", "a.pairs", "-m", "abundance"],
        ["-S", "fastq", "-1", "a.fq"],
        ["-S", "pair", "-1", "a.pairs", "-p", "bad.prefix"],
        ["-S", "pair", "-1", "a.pairs", "-a", "missing.fa"],
    ],
    ids=[
        "edge region",
        "edge region from pairs",
        "negative edge region",
        "scaffold",
        "normalization without depth",
        "no reverse reads",
        "prefix",
        "no assembly",
    ],
)
def test_pipeline_refused_keeps_the_bins_of_an_earlier_run(args, tmp_path, run_cli):
    """A refused run with --force (a rebinning into the same directory)
    stops before replacing the bins of the earlier run."""
    out = tmp_path / "out"
    bins = [out / name / "bin.fa" for name in ["overlapping_bin", "recursive_bin", "final_bin_unscaffold"]]
    for path in bins:
        path.parent.mkdir(parents=True)
        path.write_text(">c1\nACGT\n")
    with pytest.raises(ValueError):
        run_cli("pipeline", "-a", TESTS_DATA / "assembly.fa", *args, "-o", out, "-F", "-T", tmp_path / "tmp")
    assert all(path.exists() for path in bins)


def test_pipeline_refuses_earlier_bins_before_the_index(tmp_path, monkeypatch, run_cli):
    """Without --force, the bins of an earlier run are refused before the
    index of the assembly is built (an hour on a large metagenome)."""

    def generate_fasta_index(*args):
        raise AssertionError("index built before the refusal")

    monkeypatch.setattr("metator.io.generate_fasta_index", generate_fasta_index)
    out = tmp_path / "out"
    (out / "recursive_bin").mkdir(parents=True)
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    with pytest.raises(ValueError):
        run_cli("pipeline", "-1", "a.fq", "-2", "b.fq", "-a", assembly, "-S", "fastq", "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob("metator_pipeline_*.log")
    assert f"ERROR :: {out / 'recursive_bin'} already existed." in log.read_text()


def test_pipeline_without_pyfastx_index(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """No bin passed the size threshold, so pyfastx never indexed the
    assembly: the run ends without removing an index that is not there (it
    used to fail at its very end, with FileNotFoundError)."""
    monkeypatch.setattr("metator.validation.recursive_decontamination", lambda *args: os.path.join(args[8], "clustering_matrix"))
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile(TESTS_DATA / "subset.assembly.fa", assembly)
    pairs = tmp_path / "a_sorted.pairs.gz"
    pairs.write_text("## pairs format v1.0\n")
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
    for name in ["bin_summary.txt", "contig_data_final.txt"]:
        (tmp_path / name).write_text("x")
    run_cli("qc", "-a", "assembly.fa", "-e", "DpnII,HinfI", "-p", "sample", "-t", "3,4", "a.pairs", "b.pairs")
    ((args, kwargs),) = calls
    assert args[:5] == (
        "./contig_data_final.txt",
        "./bin_summary.txt",
        "assembly.fa",
        ["a.pairs", "b.pairs"],
        "./metator_qc_sample",
    )
    assert args[6:] == ("sample", False)
    assert kwargs == {"enzyme": ["DpnII", "HinfI"], "threshold": ["3", "4"]}
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
