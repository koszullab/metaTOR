# Tests of the metator command line (metator.cli): parity with the docopt
# command line it replaced, help, version and usage errors, and the
# bookkeeping of the commands (log files, temporary directories, default
# outputs), with the heavy steps replaced by stubs.

import json
import os
import shutil
import sys
from pathlib import Path

import click
import pytest

import metator.cli as mtcli
import metator.metacool as mtmc
from metator.version import __version__

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"

# The writer itself, which the stub pipeline replaces.
mtmc_write_metacool = mtmc.write_metacool

# The docopt command line that rich-click replaced (that of 1.3.10, with the
# mge options renamed and host removed), as docopt read it: see the commit
# that added the file. Commands and options added since are added to it:
# metacool, and the --no-metacool option of pipeline.
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
    positions: the help says so right after what a map gives, not after
    the contigs that a map may lack, which do not call for it."""
    output = " ".join(run_cli("network", "--help").output.split())
    assert "which it keeps. --edge must then be 0. The map may lack contigs of the assembly" in output


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
        # The pairs files must still be there.
        calls["write_metacool"] = (pairs_files, [os.path.exists(path) for path in pairs_files], contig_data, outfile, kwargs)

    monkeypatch.setattr("metator.align.get_contact_pairs", get_contact_pairs)
    monkeypatch.setattr("metator.network.create_contig_data", create_contig_data)
    monkeypatch.setattr("metator.network.alignment_to_contacts", alignment_to_contacts)
    monkeypatch.setattr("metator.network.metacool_to_contacts", metacool_to_contacts)
    monkeypatch.setattr("metator.partition.partition", partition)
    monkeypatch.setattr("metator.validation.recursive_decontamination", recursive_decontamination)
    monkeypatch.setattr("metator.io.generate_fasta_index", generate_fasta_index)
    monkeypatch.setattr("metator.metacool.write_metacool", write_metacool)
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

    # The contact map is written last, from the pairs files of the run.
    pairs_files, existed, contig_data, outfile, kwargs = calls["write_metacool"]
    assert pairs_files == [str(path) for path in pairs] and existed == [True, True]
    assert (contig_data, outfile) == (str(out / "contig_data_final.txt"), str(out / "metator.metacool"))
    assert kwargs["bin_summary"] == str(out / "bin_summary.txt")
    assert (kwargs["assembly_name"], kwargs["threads"], kwargs["force"]) == ("assembly", 2, True)
    # The bin summary goes to metacooler as it is: no temporary copy.
    assert "tmp_dir" not in kwargs
    metadata = kwargs["metadata"]
    assert metadata["command"].startswith("metator pipeline -1 a1.fq,b1.fq -2 a2.fq,b2.fq -a")
    assert metadata["parameters"]["start"] == start and metadata["parameters"]["no_metacool"] is False


@pytest.mark.parametrize(
    "options, written",
    [([], "metator.metacool"), (["-p", "sample_1"], "sample_1.metacool"), (["--no-metacool"], None)],
)
def test_pipeline_metacool_name(options, written, stub_pipeline, tmp_path, run_cli):
    """<prefix>.metacool, metator.metacool without --prefix, and none with
    --no-metacool."""
    out = tmp_path / "out"
    run_cli("pipeline", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp", *options)
    if written:
        pairs_files, _, _, outfile, _ = stub_pipeline["write_metacool"]
        assert pairs_files == ["a.pairs"] and outfile == str(out / written)
    else:
        assert "write_metacool" not in stub_pipeline
    assert not list(out.glob("*.metacool"))


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


def test_pipeline_non_ascii_contigs_skip_the_metacool(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """Contig names that are not ASCII: a warning, no .metacool, and the run
    succeeds with its other outputs."""
    monkeypatch.setattr("metator.metacool.write_metacool", mtmc_write_metacool)
    _final_tables(monkeypatch, "contig_é")
    pairs = tmp_path / "a.pairs"
    pairs.write_text("## pairs format v1.0\nr1\tcontig_é\t10\tc2\t20\t+\t-\n")
    out = tmp_path / "out"
    run_cli("pipeline", "-1", pairs, "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    assert not list(out.glob("*.metacool*"))
    assert (out / "network.txt").exists() and (out / "contig_data_final.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert (
        f"WARNING :: {out / 'metator.metacool'} is not written: A .metacool file needs ASCII contig names (cooler reads "
        "them as ASCII text), and 1 contig name(s) are not: 'contig_é'. The other outputs of the run are kept."
    ) in logged
    assert "Traceback" not in logged


def test_pipeline_misaligned_bin_summary_is_a_warning(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """metacooler refuses a bin summary whose header has more fields than its
    rows: a warning with its reason, without traceback, no map, and the
    other outputs of the run."""
    monkeypatch.setattr("metator.metacool.write_metacool", mtmc_write_metacool)
    _final_tables(monkeypatch, "c1", summary="\tCompleteness\tContigs\nmetator_00001_00000\t0.9\n")
    pairs = tmp_path / "a.pairs"
    pairs.write_text("## pairs format v1.0\nr1\tc1\t10\tc2\t20\t+\t-\n")
    out = tmp_path / "out"
    run_cli("pipeline", "-1", pairs, "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    assert not list(out.glob("*.metacool*"))
    assert (out / "network.txt").exists() and (out / "bin_summary.txt").exists()
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert (
        f"WARNING :: {out / 'metator.metacool'} is not written: Genomes table is misaligned: {out / 'bin_summary.txt'}\n"
        in logged
    )
    assert "The other outputs of the run are kept." in logged
    assert "Traceback" not in logged


def test_pipeline_metacool_failure_is_a_warning(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A .metacool that fails to write costs a warning, with the traceback in
    the log, not the run: the exit status is 0 and the pairs files go as
    usual."""

    def write_metacool(*args, **kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr("metator.metacool.write_metacool", write_metacool)
    out = tmp_path / "out"
    run_cli("pipeline", "-1", "a.fq", "-2", "b.fq", "-a", TESTS_DATA / "assembly", "-S", "bam", "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    assert f"WARNING :: {out / 'metator.metacool'} is not written: RuntimeError: disk full. The other outputs" in logged
    assert "Traceback" in logged
    assert (out / "network.txt").exists()
    assert not (out / "alignment_0_sorted.pairs.gz").exists()
    assert not (out / "alignment_0_sorted.pairs.gz.px2").exists()


@pytest.mark.parametrize("case", ["written", "failed", "skipped"])
def test_pipeline_earlier_metacool(case, stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A map left by an earlier run in the output directory is removed
    before the new one is written: it would pass for this run's map when
    the new one fails. --no-metacool leaves it alone."""
    out = tmp_path / "out"
    out.mkdir()
    earlier = out / "metator.metacool"
    earlier.write_text("earlier run")
    seen = []

    def write_metacool(pairs_files, contig_data, outfile, **kwargs):
        seen.append(os.path.exists(outfile))
        if case == "failed":
            raise RuntimeError("disk full")
        Path(outfile).write_text("this run")

    monkeypatch.setattr("metator.metacool.write_metacool", write_metacool)
    options = ["--no-metacool"] if case == "skipped" else []
    run_cli("pipeline", "-1", "a.pairs", "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp", *options)
    (log,) = out.glob("metator_pipeline_*.log")
    logged = log.read_text()
    if case == "skipped":
        assert seen == [] and earlier.read_text() == "earlier run"
        assert "Removing" not in logged
        return
    assert seen == [False]
    assert f"INFO :: Removing {earlier}, the contact map of an earlier run." in logged
    if case == "failed":
        assert not earlier.exists()
        assert f"WARNING :: {earlier} is not written: RuntimeError: disk full." in logged
    else:
        assert earlier.read_text() == "this run"


def test_pipeline_metacool_written(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """With the real writer, from the tables of the (stub) validation."""
    monkeypatch.setattr("metator.metacool.write_metacool", mtmc_write_metacool)
    _final_tables(monkeypatch, "c1")
    pairs = tmp_path / "a.pairs"
    pairs.write_text("## pairs format v1.0\nr1\tc1\t10\tc2\t20\t+\t-\nr2\tc1\t10\tc1\t90\t+\t-\n")
    out = tmp_path / "out"
    run_cli("pipeline", "-1", pairs, "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    import cooler
    import metacooler

    metacooler.validate(str(out / "metator.metacool"))
    clr = cooler.Cooler(str(out / "metator.metacool"))
    assert clr.pixels()[:]["count"].sum() == 2
    assert clr.info["genome-assembly"] == "assembly"
    genomes = metacooler.MetaCooler(str(out / "metator.metacool")).genomes
    assert list(genomes.name) == ["__unbinned__", "metator_00001_00000"]
    assert list(genomes.Completeness.iloc[1:]) == [0.9]
    # The contig columns of the contig data: its GC_content (Hit is - on
    # every contig), not its ID or Final_bin.
    contigs = metacooler.MetaCooler(str(out / "metator.metacool")).contigs.set_index("name")
    assert list(contigs.columns) == ["length", "genome_id", "GC_content"]
    assert contigs.GC_content.to_dict() == {"c1": 40.5, "c2": 35.0}


def test_pipeline_metacool_without_final_bin(stub_pipeline, tmp_path, monkeypatch, run_cli):
    """A run with no final bin keeps its map: every contig in __unbinned__,
    with its contacts, and no column of the bin summary, which has no row."""
    monkeypatch.setattr("metator.metacool.write_metacool", mtmc_write_metacool)
    _final_tables(monkeypatch, "c1", final_bin="ND")
    pairs = tmp_path / "a.pairs"
    pairs.write_text("## pairs format v1.0\nr1\tc1\t10\tc2\t20\t+\t-\nr2\tc1\t10\tc1\t90\t+\t-\n")
    out = tmp_path / "out"
    run_cli("pipeline", "-1", pairs, "-a", _assembly(tmp_path), "-S", "pair", "-o", out, "-T", tmp_path / "tmp")
    import cooler
    import metacooler

    metacooler.validate(str(out / "metator.metacool"), full=True)
    assert cooler.Cooler(str(out / "metator.metacool")).pixels()[:]["count"].sum() == 2
    mc = metacooler.MetaCooler(str(out / "metator.metacool"))
    assert list(mc.genomes.name) == ["__unbinned__"] and "Completeness" not in mc.genomes.columns
    assert sorted(mc.contigs.name) == ["c1", "c2"]
    assert list(mc.contigs.columns) == ["name", "length", "genome_id", "GC_content"]
    (log,) = out.glob("metator_pipeline_*.log")
    assert "is not written" not in log.read_text()


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
    read or written, and the map of the bins is not written (it is made
    from pairs files). The input map stays, even where the map of the bins
    would go."""
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
    message = f"The contact map of the bins is not written: it is made from the pairs files, and this run started from the contact map {metacool}."
    assert (message in log.read_text()) == (not no_metacool)


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
    "edge region": (
        ["-1", "map.metacool", "-E", "100"],
        "--edge 100 cannot be applied with --start metacool: a contact map holds the contacts of bins, not the "
        "positions of the reads. Leave --edge at 0, or start from the pairs files.",
    ),
    "several maps": (
        ["-1", "a.metacool,b.metacool,c.metacool"],
        "--start metacool takes one .metacool file (its libraries are already summed), given with -1/--forward: got 3 files.",
    ),
}
PIPELINE_METACOOL_START_ERRORS = {
    "no map": ([], "--start metacool needs the .metacool file, given with -1/--forward."),
    "scaffold": (["-1", "map.metacool", "--scaffold"], "--scaffold reads the pairs files: it cannot run with --start metacool."),
}


@pytest.mark.parametrize(
    "command, case",
    [(command, case) for command in ["network", "pipeline"] for case in METACOOL_START_ERRORS]
    + [("pipeline", case) for case in PIPELINE_METACOOL_START_ERRORS],
)
def test_metacool_start_refused(command, case, monkeypatch, tmp_path, run_cli):
    """Options a contact map cannot serve are refused before any work, with
    a message in the log file."""

    def not_called(*args):
        raise AssertionError("refused before any work")

    monkeypatch.setattr("metator.network.create_contig_data", not_called)
    args, message = {**METACOOL_START_ERRORS, **PIPELINE_METACOOL_START_ERRORS}[case]
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run_cli(command, *args, "-a", TESTS_DATA / "assembly.fa", "-S", "metacool", "-o", out, "-T", tmp_path / "tmp")
    (log,) = out.glob(f"metator_{command}_*.log")
    assert f"ERROR :: {message}" in log.read_text()


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
