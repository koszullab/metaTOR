# Tests of the metator command line entry point.

import shutil
import sys

import pytest

import metator.main as mtm
from metator.version import __version__


def _run(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["metator", *argv])
    mtm.main()


def test_help(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "--help")
    assert not exit.value.code
    out = capsys.readouterr().out
    assert "metator [-hv] <command> [<args>...]" in out
    for command in ["network", "partition", "validation", "pipeline", "qc", "contactmap", "scaffold", "pairs", "host", "mge"]:
        assert f"    {command} " in out


def test_version(monkeypatch, capsys):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch, "--version")
    assert not exit.value.code
    assert capsys.readouterr().out.strip() == __version__


def test_no_command(monkeypatch):
    with pytest.raises(SystemExit) as exit:
        _run(monkeypatch)
    assert "usage:" in str(exit.value.code)


def test_unknown_command(monkeypatch, capsys):
    with pytest.raises(SystemExit):
        _run(monkeypatch, "binning")
    assert "Unknown command." in capsys.readouterr().out


def test_command_help(monkeypatch, capsys):
    with pytest.raises(SystemExit):
        _run(monkeypatch, "pairs", "--help")
    assert "pairs [--force] [--remove] [--threads=1] <pairsfile>..." in capsys.readouterr().out


@pytest.mark.skipif(shutil.which("pairtools") is None, reason="pairtools is not installed")
def test_command_runs(monkeypatch, tmp_path):
    pairs = tmp_path / "sample.pairs"
    shutil.copyfile("tests_data/alignment.pairs", pairs)
    _run(monkeypatch, "pairs", "--remove", str(pairs))
    assert (tmp_path / "sample_sorted.pairs.gz").exists()
    assert (tmp_path / "sample_sorted.pairs.gz.px2").exists()
    assert not pairs.exists()
