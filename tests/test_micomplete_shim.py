# Tests for the miComplete compatibility shim (metator._micomplete).
#
# metator runs the shim as `python -m metator._micomplete`, in a subprocess
# that pytest-cov does not measure, so these tests import it in-process.
# Importing it patches Bio.SeqUtils and may put a stand-in distutils into
# sys.modules: the `shim_state` fixture undoes all of it after each test, and
# `_nothing_leaks` checks that the module as a whole left that state as it
# found it.

import importlib
import runpy
import shutil
import sys
import types

import Bio.SeqUtils
import metator
import pytest
from Bio.Seq import Seq

SHIM = "metator._micomplete"
# sys.modules entries that importing the shim, or miComplete through main(),
# may add or replace.
WATCHED_MODULES = ("distutils", "micomplete", SHIM)
_ABSENT = object()

# Python < 3.12 warns on every import of distutils, and the fixtures make the
# shim import it afresh in each test.
pytestmark = pytest.mark.filterwarnings("ignore:The distutils package is deprecated:DeprecationWarning")


def _watched_modules():
    return {
        name: module
        for name, module in sys.modules.items()
        if any(name == root or name.startswith(root + ".") for root in WATCHED_MODULES)
    }


def _snapshot():
    """The global state that importing or running the shim can change."""
    return {
        "Bio.SeqUtils attributes": dict(vars(Bio.SeqUtils)),
        "sys.modules entries": _watched_modules(),
        "metator package attribute": {"_micomplete": getattr(metator, "_micomplete", _ABSENT)},
    }


def _differences(before, after):
    """Per part of two snapshots, the names not bound to the very same object."""
    changed = {}
    for part in before:
        names = before[part].keys() | after[part].keys()
        moved = sorted(name for name in names if before[part].get(name, _ABSENT) is not after[part].get(name, _ABSENT))
        if moved:
            changed[part] = moved
    return changed


def _restore(snapshot):
    sequtils = snapshot["Bio.SeqUtils attributes"]
    for name in set(vars(Bio.SeqUtils)) - set(sequtils):
        delattr(Bio.SeqUtils, name)
    for name, value in sequtils.items():
        if vars(Bio.SeqUtils).get(name, _ABSENT) is not value:
            setattr(Bio.SeqUtils, name, value)

    modules = snapshot["sys.modules entries"]
    for name in set(_watched_modules()) - set(modules):
        del sys.modules[name]
    sys.modules.update(modules)

    package_attribute = snapshot["metator package attribute"]["_micomplete"]
    if package_attribute is _ABSENT:
        if hasattr(metator, "_micomplete"):
            delattr(metator, "_micomplete")
    else:
        metator._micomplete = package_attribute


@pytest.fixture(scope="module", autouse=True)
def _nothing_leaks():
    """Fail if the tests below, fixtures included, leave any shim patch behind."""
    before = _snapshot()
    yield
    assert _differences(before, _snapshot()) == {}


@pytest.fixture
def shim_state():
    """Put back Bio.SeqUtils, sys.modules and the metator package after a test."""
    before = _snapshot()
    yield
    _restore(before)
    assert _differences(before, _snapshot()) == {}


def _import_shim():
    """Import the shim afresh, so that its import-time patches run again."""
    sys.modules.pop(SHIM, None)
    return importlib.import_module(SHIM)


def _remove_gc():
    if hasattr(Bio.SeqUtils, "GC"):
        delattr(Bio.SeqUtils, "GC")


def _make_distutils_unimportable():
    # What `from distutils import spawn` meets on Python >= 3.12 without
    # setuptools: an ImportError.
    sys.modules["distutils"] = None
    sys.modules.pop("distutils.spawn", None)


def test_import_applies_both_patches(shim_state):
    """The patches run at import, which is how the worker processes of
    miComplete's multiprocessing pool get them under spawn and forkserver."""
    _remove_gc()
    _make_distutils_unimportable()
    _import_shim()
    from distutils import spawn

    assert Bio.SeqUtils.GC("GCAT") == 50.0
    assert spawn.find_executable("sh") == shutil.which("sh")


@pytest.mark.parametrize(
    ("seq", "expected"),
    [
        # 100 000 nt, 30 015 of them G or C: 100 * (30015 / 100000), what the
        # gc_fraction-based shim returned, is 30.014999999999997.
        ("G" * 15_000 + "C" * 15_015 + "A" * 69_985, 30015 * 100.0 / 100000),
        # Both cases count, S counts as G or C, N counts in the length only.
        ("ggccSsNNAT", 60.0),
        ("NNNN", 0.0),
        ("", 0.0),
        # 14.375 rounds to 14.38, as miComplete does with Biopython 1.80;
        # 100 * (23 / 160) is 14.374999999999998, which rounds to 14.37.
        ("G" * 23 + "A" * 137, 14.375),
    ],
)
def test_gc_is_computed_as_biopython_1_80(shim_state, seq, expected):
    shim = _import_shim()
    _remove_gc()
    shim.patch_biopython_gc()
    gc = Bio.SeqUtils.GC
    assert gc(seq) == expected
    assert gc(Seq(seq)) == expected
    assert round(gc(seq), 2) == round(expected, 2)


@pytest.mark.parametrize("seq", ["G" * 15_000 + "C" * 15_015 + "A" * 69_985, "G" * 23 + "A" * 137])
def test_gc_is_not_100_gc_fraction(shim_state, seq):
    """The two differ in the last bit on these sequences: dividing first, as
    gc_fraction does, is what moved miComplete's rounded GC-content by 0.01."""
    shim = _import_shim()
    _remove_gc()
    shim.patch_biopython_gc()
    assert Bio.SeqUtils.GC(seq) != 100.0 * Bio.SeqUtils.gc_fraction(seq, ambiguous="ignore")


def test_gc_left_alone_when_present(shim_state):
    shim = _import_shim()

    def native_gc(seq):
        return 42.0

    Bio.SeqUtils.GC = native_gc
    shim.patch_biopython_gc()
    assert Bio.SeqUtils.GC is native_gc


def test_distutils_left_alone_when_importable(shim_state):
    shim = _import_shim()
    distutils = types.ModuleType("distutils")
    spawn = types.ModuleType("distutils.spawn")
    distutils.spawn = spawn
    sys.modules["distutils"] = distutils
    sys.modules["distutils.spawn"] = spawn
    shim.patch_distutils_spawn()
    assert sys.modules["distutils"] is distutils
    assert sys.modules["distutils.spawn"] is spawn
    assert not hasattr(spawn, "find_executable")


def test_distutils_spawn_provided_when_missing(shim_state, tmp_path):
    shim = _import_shim()
    _make_distutils_unimportable()
    shim.patch_distutils_spawn()
    from distutils import spawn

    assert sys.modules["distutils.spawn"] is spawn
    assert sys.modules["distutils"].spawn is spawn

    tool = tmp_path / "metator-shim-test-tool"
    tool.write_text("#!/bin/sh\n")
    tool.chmod(0o755)
    path = str(tmp_path)
    assert spawn.find_executable(tool.name, path=path) == shutil.which(tool.name, path=path) == str(tool)
    assert spawn.find_executable(tool.name, path) == str(tool)
    # Not on PATH.
    assert spawn.find_executable(tool.name) == shutil.which(tool.name)
    assert spawn.find_executable("sh") == shutil.which("sh") is not None


def test_main_runs_micomplete_cli(shim_state, monkeypatch, capsys):
    shim = _import_shim()
    argv = ["x", "-h"]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit) as exit_info:
        shim.main()
    assert exit_info.value.code == 0
    assert argv[0] == "miComplete"
    assert capsys.readouterr().out.startswith("usage: miComplete")


def test_run_as_module(shim_state, monkeypatch, capsys):
    """`python -m metator._micomplete -h`, as validation runs it."""
    # runpy warns when the module it runs as __main__ is already imported.
    sys.modules.pop(SHIM, None)
    monkeypatch.setattr(sys, "argv", ["x", "-h"])
    with pytest.raises(SystemExit) as exit_info:
        runpy.run_module(SHIM, run_name="__main__")
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.startswith("usage: miComplete")
