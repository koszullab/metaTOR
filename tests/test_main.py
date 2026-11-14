# Tests of the metator command line entry point.
#
# They check what a user sees (the commands listed, the options of a command,
# a failure on a missing or unknown command), not the layout of the help
# text, which belongs to the command line parsing library.

import re
import shutil
import sys

import pytest

import metator.main as mtm
from metator.version import __version__

COMMANDS = ["network", "partition", "validation", "pipeline", "qc", "contactmap", "scaffold", "pairs", "host", "mge"]


def _run(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["metator", *argv])
    mtm.main()


def _shown(exit, capsys):
    """Everything shown to the user: printed output and exit message."""
    captured = capsys.readouterr()
    message = "" if isinstance(exit.value.code, int) or exit.value.code is None else str(exit.value.code)
    return captured.out + captured.err + message


def test_help(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "--help")
    assert not exit.value.code
    out = capsys.readouterr().out
    for command in COMMANDS:
        # Listed at the start of a line, after indentation or table borders.
        assert re.search(rf"^\W*{command}\b", out, flags=re.MULTILINE), command
    assert "--version" in out


def test_version(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "--version")
    assert not exit.value.code
    assert __version__ in capsys.readouterr().out


def test_no_command(monkeypatch, capsys):
    """Without a command, the usage is shown and nothing runs."""
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch)
    assert "usage" in _shown(exit, capsys).lower()


def test_unknown_command(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "binning")
    # A failure: a non-zero status, or an error message as exit code.
    assert exit.value.code not in (0, None)
    assert "command" in _shown(exit, capsys).lower()


def test_command_help(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "pairs", "--help")
    assert not exit.value.code
    out = capsys.readouterr().out
    for option in ["--force", "--remove", "--threads"]:
        assert option in out
    assert "pairsfile" in out.lower()


@pytest.mark.skipif(shutil.which("pairtools") is None, reason="pairtools is not installed")
def test_command_runs(monkeypatch, tmp_path):
    pairs = tmp_path / "sample.pairs"
    shutil.copyfile("tests_data/alignment.pairs", pairs)
    _run(monkeypatch, "pairs", "--remove", str(pairs))
    assert (tmp_path / "sample_sorted.pairs.gz").exists()
    assert (tmp_path / "sample_sorted.pairs.gz.px2").exists()
    assert not pairs.exists()
