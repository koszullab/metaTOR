# Tests of metator._micomplete, the shim that runs miComplete with recent
# Biopython and Python, exercised in-process by simulating what is missing.

import shutil
import sys
import types

import Bio.SeqUtils
import pytest


@pytest.fixture
def shim(monkeypatch):
    """The shim module. Importing it applies the patches: whatever they add
    to this process is removed afterwards."""
    had_gc, gc = hasattr(Bio.SeqUtils, "GC"), getattr(Bio.SeqUtils, "GC", None)
    modules = {name: sys.modules.get(name) for name in ("distutils", "distutils.spawn")}
    from metator import _micomplete

    yield _micomplete
    if had_gc:
        Bio.SeqUtils.GC = gc
    elif hasattr(Bio.SeqUtils, "GC"):
        del Bio.SeqUtils.GC
    for name, module in modules.items():
        if module is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = module


def test_gc_restored_when_missing(shim, monkeypatch):
    """Biopython >= 1.82 has no Bio.SeqUtils.GC: the shim adds the old one."""
    monkeypatch.delattr(Bio.SeqUtils, "GC", raising=False)
    shim.patch_biopython_gc()
    assert Bio.SeqUtils.GC("GGCCAT") == pytest.approx(100 * 4 / 6)
    # As in Biopython <= 1.81: S counts as G/C, ambiguous bases count in the
    # length, lower case is accepted.
    assert Bio.SeqUtils.GC("GSNN") == pytest.approx(50)
    assert Bio.SeqUtils.GC("gcat") == pytest.approx(50)


def test_gc_kept_when_present(shim, monkeypatch):
    def existing(seq):
        return 42.0

    monkeypatch.setattr(Bio.SeqUtils, "GC", existing, raising=False)
    shim.patch_biopython_gc()
    assert Bio.SeqUtils.GC is existing


def test_distutils_spawn_added_when_missing(shim, monkeypatch):
    """Python >= 3.12 has no distutils: the shim provides
    distutils.spawn.find_executable."""
    monkeypatch.setitem(sys.modules, "distutils", None)
    monkeypatch.setitem(sys.modules, "distutils.spawn", None)
    with pytest.raises(ImportError):
        from distutils import spawn
    shim.patch_distutils_spawn()
    from distutils import spawn

    assert spawn is sys.modules["distutils.spawn"]
    assert spawn.find_executable("sh") == shutil.which("sh")
    assert spawn.find_executable("no-such-program-anywhere") is None


def test_distutils_kept_when_present(shim, monkeypatch):
    existing = types.ModuleType("distutils.spawn")
    package = types.ModuleType("distutils")
    package.spawn = existing
    monkeypatch.setitem(sys.modules, "distutils", package)
    monkeypatch.setitem(sys.modules, "distutils.spawn", existing)
    shim.patch_distutils_spawn()
    assert sys.modules["distutils.spawn"] is existing
    assert not hasattr(existing, "find_executable")


def test_main_runs_micomplete_cli(shim, monkeypatch):
    """main() hands the arguments to miComplete's CLI, under its own name."""
    calls = []
    package = types.ModuleType("micomplete")
    cli = types.ModuleType("micomplete.micomplete")
    cli.main = lambda: calls.append(list(sys.argv)) or 3
    package.micomplete = cli
    monkeypatch.setitem(sys.modules, "micomplete", package)
    monkeypatch.setitem(sys.modules, "micomplete.micomplete", cli)
    monkeypatch.setattr(sys, "argv", ["/somewhere/_micomplete.py", "seqs.tab", "--hmms", "Bact105"])
    assert shim.main() == 3
    assert calls == [["miComplete", "seqs.tab", "--hmms", "Bact105"]]
