# Tests of the metator command line entry points.
#
# They check what a user sees (the commands listed, the options of a command,
# a failure on a missing or unknown command), not the layout of the help
# text, which belongs to the command line parsing library. The `metator`
# command runs metator.cli.cli; metator.main.main, the former entry point,
# and `python -m metator` run the same command line.

import re
import runpy
import sys

import pytest

import metator.cli as mtcli
import metator.main as mtm
from metator.version import __version__
from tests.ansi import plain

COMMANDS = ["network", "partition", "validation", "pipeline", "mge"]


def _run(monkeypatch, *argv):
    """Run metator.main.main() on the command line `metator *argv`."""
    monkeypatch.setattr(sys, "argv", ["metator", *argv])
    mtm.main()


def _shown(exit, capsys):
    """Everything shown to the user: printed output and exit message."""
    captured = capsys.readouterr()
    message = "" if isinstance(exit.value.code, int) or exit.value.code is None else str(exit.value.code)
    return plain(captured.out + captured.err + message)


def test_main_is_the_cli():
    assert mtm.main is mtcli.cli


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
    for option in ["--force", "--network", "--threads"]:
        assert option in out
    assert "map" in out.lower()


def test_command_runs(monkeypatch, tmp_path, private_link):
    out = tmp_path / "out"
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "partition", "-a", private_link("tests_data/subset.assembly.fa"), "-c",
             "tests_data/contig_data_network.txt", "-n", "tests_data/network.txt", "-i", "2", "-s", "30000", "-o",
             str(out), "-T", str(tmp_path / "tmp"))  # fmt: skip
    assert not exit.value.code
    assert (out / "contig_data_partition.txt").exists()
    assert any((out / "overlapping_clusters").iterdir())


@pytest.mark.parametrize("module", ["metator", "metator.main"])
def test_python_m(monkeypatch, capsys, module):
    """`python -m metator` and `python -m metator.main` run the command line."""
    monkeypatch.setattr(sys, "argv", [module, "--version"])
    # Run from a fresh copy, as python does.
    monkeypatch.delitem(sys.modules, "metator.__main__" if module == "metator" else module, raising=False)
    with pytest.raises(SystemExit) as exit:
        runpy.run_module(module, run_name="__main__", alter_sys=True)
    assert not exit.value.code
    assert capsys.readouterr().out == f"{__version__}\n"
