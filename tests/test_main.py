# Tests of the metator command line entry points.
#
# They check what a user sees (the commands listed, the options of a command,
# a failure on a missing or unknown command), not the layout of the help
# text, which belongs to the command line parsing library. The `metator`
# command runs metator.cli.cli, and so does `python -m metator`.

import re
import runpy
import sys

import pytest

import metator.cli as mtcli
from metator.version import __version__
from tests.ansi import plain

COMMANDS = ["network", "partition", "validation", "pipeline", "mge"]


def _run(monkeypatch, *argv):
    """Run metator.cli.cli() on the command line `metator *argv`."""
    monkeypatch.setattr(sys, "argv", ["metator", *argv])
    mtcli.cli()


def _shown(exit, capsys):
    """Everything shown to the user: printed output and exit message."""
    captured = capsys.readouterr()
    message = "" if isinstance(exit.value.code, int) or exit.value.code is None else str(exit.value.code)
    return plain(captured.out + captured.err + message)


def test_help(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "--help")
    assert not exit.value.code
    out = plain(capsys.readouterr().out)
    for command in COMMANDS:
        # Listed at the start of a line, after indentation or table borders.
        assert re.search(rf"^\W*{command}\b", out, flags=re.MULTILINE), command
    assert not re.search(r"^\W*host\b", out, flags=re.MULTILINE)
    assert "--version" in out


@pytest.mark.parametrize("option", ["--version", "-v"])
def test_version(monkeypatch, capsys, option):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, option)
    assert not exit.value.code
    assert capsys.readouterr().out == f"{__version__}\n"


def test_no_command(monkeypatch, capsys):
    """Without a command, the usage is shown and nothing runs."""
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch)
    assert "usage" in _shown(exit, capsys).lower()


# host: removed in 1.4.0, it always failed; the others in 2.0.0.
@pytest.mark.parametrize("command", ["binning", "host", "qc", "contactmap", "metacool", "scaffold", "pairs"])
def test_unknown_command(monkeypatch, capsys, command):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, command)
    # A failure: a non-zero status, or an error message as exit code.
    assert exit.value.code not in (0, None)
    assert "command" in _shown(exit, capsys).lower()


def test_command_help(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "partition", "--help")
    assert not exit.value.code
    out = plain(capsys.readouterr().out)
    for option in ["--force", "--inplace", "--threads"]:
        assert option in out
    assert "map" in out.lower() and "--network" not in out


def test_command_runs(monkeypatch, tmp_path, metagenome_network):
    """partition MAP, on a map whose network metator network stored (see
    metagenome_network): it stores its partition in the map of its output
    directory."""
    import shutil

    import metacooler

    out = tmp_path / "out"
    metacool = shutil.copyfile(metagenome_network["map"], tmp_path / "network.metacool")
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "partition", str(metacool), "-i", "2", "-o", str(out), "-T", str(tmp_path / "tmp"))
    assert not exit.value.code
    assert [entry["name"] for entry in metacooler.MetaCooler(str(out / "network.metacool")).cluster_rounds] == [
        "core",
        "overlapping",
    ]


def test_python_m(monkeypatch, capsys):
    """`python -m metator` runs the command line."""
    monkeypatch.setattr(sys, "argv", ["metator", "--version"])
    # Run from a fresh copy, as python does.
    monkeypatch.delitem(sys.modules, "metator.__main__", raising=False)
    with pytest.raises(SystemExit) as exit:
        runpy.run_module("metator", run_name="__main__", alter_sys=True)
    assert not exit.value.code
    assert capsys.readouterr().out == f"{__version__}\n"


def test_python_m_metator_main_is_gone():
    """`python -m metator.main`, the former entry point, is gone in 2.0."""
    with pytest.raises(ImportError):
        runpy.run_module("metator.main", run_name="__main__", alter_sys=True)
