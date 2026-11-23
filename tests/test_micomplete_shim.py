# Tests for the miComplete compatibility shim (metator._micomplete).
#
# metator runs the shim as `python -m metator._micomplete`, in a subprocess:
# the last tests run it so, the others import it in-process. Importing it
# patches Bio.SeqUtils, may put a stand-in distutils into sys.modules, and
# imports miComplete to replace its prodigal and hmmsearch calls: the
# `shim_state` fixture undoes all of it after each test, and `_nothing_leaks`
# checks that the module as a whole left that state as it found it.
#
# The pyrodigal and pyhmmer stand-ins are checked against what prodigal 2.6.3
# and HMMER 3.4 (hmmsearch) gave on the bins of tests_data.

import importlib
import logging
import os
import runpy
import shutil
import subprocess
import sys
import types
from pathlib import Path

import Bio.SeqUtils
import metator
import pytest
from Bio import SeqIO
from Bio.Seq import Seq

SHIM = "metator._micomplete"
# sys.modules entries that importing the shim, or miComplete through main(),
# may add or replace.
WATCHED_MODULES = ("distutils", "micomplete", SHIM)
# The names of miComplete that the shim replaces: (module, class, name).
REPLACED_NAMES = (
    ("micomplete.micomplete", None, "create_proteome"),
    ("micomplete.micomplete", None, "shutil"),
    ("micomplete.completeness", "calcCompleteness", "hmm_search"),
)
_ABSENT = object()

BINS = Path(__file__).resolve().parent.parent / "tests_data" / "overlapping_bin"

# Python < 3.12 warns on every import of distutils, and the fixtures make the
# shim import it afresh in each test.
pytestmark = pytest.mark.filterwarnings("ignore:The distutils package is deprecated:DeprecationWarning")


def _watched_modules():
    return {
        name: module
        for name, module in sys.modules.items()
        if any(name == root or name.startswith(root + ".") for root in WATCHED_MODULES)
    }


def _replaced_owner(module_name, class_name):
    """The module, or its class, that holds a name the shim replaces, if the
    module is imported."""
    module = sys.modules.get(module_name)
    if module is None or class_name is None:
        return module
    return vars(module).get(class_name)


def _replaced_names():
    """What the names the shim replaces are bound to, in the miComplete
    modules imported."""
    bound = {}
    for module_name, class_name, name in REPLACED_NAMES:
        owner = _replaced_owner(module_name, class_name)
        if owner is not None:
            bound[".".join(filter(None, (module_name, class_name, name)))] = vars(owner).get(name, _ABSENT)
    return bound


def _snapshot():
    """The global state that importing or running the shim can change."""
    return {
        "Bio.SeqUtils attributes": dict(vars(Bio.SeqUtils)),
        "sys.modules entries": _watched_modules(),
        "metator package attribute": {"_micomplete": getattr(metator, "_micomplete", _ABSENT)},
        "miComplete names": _replaced_names(),
        # The shim gives miComplete a which() of its own, not to everyone.
        "shutil attributes": {"which": shutil.which},
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

    # miComplete modules imported before: put back what the shim replaced in
    # them (those imported since are gone from sys.modules).
    names = snapshot["miComplete names"]
    for module_name, class_name, name in REPLACED_NAMES:
        key = ".".join(filter(None, (module_name, class_name, name)))
        if key not in names:
            continue
        owner = _replaced_owner(module_name, class_name)
        if names[key] is _ABSENT:
            if name in vars(owner):
                delattr(owner, name)
        else:
            setattr(owner, name, names[key])


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


def test_import_applies_the_patches(shim_state):
    """The patches run at import, which is how the worker processes of
    miComplete's multiprocessing pool get them under spawn and forkserver."""
    _remove_gc()
    _make_distutils_unimportable()
    shim = _import_shim()
    from distutils import spawn

    import micomplete.micomplete
    from micomplete.completeness import calcCompleteness

    assert Bio.SeqUtils.GC("GCAT") == 50.0
    assert spawn.find_executable("sh") == shutil.which("sh")
    assert micomplete.micomplete.create_proteome is shim.create_proteome
    assert calcCompleteness.hmm_search is shim.hmm_search
    assert micomplete.micomplete.shutil.which is shim.which


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


def test_main_passes_arguments_and_exit_status(shim_state, monkeypatch):
    """main() hands the arguments to miComplete's CLI, under its own name, and
    returns its exit status (miComplete stubbed)."""
    shim = _import_shim()
    calls = []
    package = types.ModuleType("micomplete")
    cli = types.ModuleType("micomplete.micomplete")
    cli.main = lambda: calls.append(list(sys.argv)) or 3
    package.micomplete = cli
    # Not monkeypatch.setitem: its undo, which comes after shim_state's, would
    # put back the miComplete that importing the shim imported. shim_state
    # removes the stubs.
    sys.modules["micomplete"] = package
    sys.modules["micomplete.micomplete"] = cli
    monkeypatch.setattr(sys, "argv", ["/somewhere/_micomplete.py", "seqs.tab", "--hmms", "Bact105"])
    assert shim.main() == 3
    assert calls == [["miComplete", "seqs.tab", "--hmms", "Bact105"]]


def test_run_as_module(shim_state, monkeypatch, capsys):
    """`python -m metator._micomplete -h`, as validation runs it."""
    # runpy warns when the module it runs as __main__ is already imported.
    sys.modules.pop(SHIM, None)
    monkeypatch.setattr(sys, "argv", ["x", "-h"])
    with pytest.raises(SystemExit) as exit_info:
        runpy.run_module(SHIM, run_name="__main__")
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.startswith("usage: miComplete")


# -- pyrodigal and pyhmmer, in place of prodigal and hmmsearch ----------------


def _marker_hmms(markers):
    """Path to a marker set of miComplete (Bact105 or Arch131)."""
    import micomplete

    return str(Path(micomplete.__file__).parent / "share" / f"{markers}.hmm")


def _first_contigs(out, length):
    """Write to out the first contigs of bin metator_00009_00000, as long as
    they total less than length bp."""
    contigs, total = [], 0
    for contig in SeqIO.parse(BINS / "metator_00009_00000.fa", "fasta"):
        if total + len(contig) >= length:
            break
        contigs.append(contig)
        total += len(contig)
    SeqIO.write(contigs, out, "fasta")
    return out


def test_create_proteome_finds_the_genes_of_prodigal(shim_state, recwarn):
    """The 46 genes that prodigal finds in this bin, from the first to the
    last, written as prodigal writes them, to <fasta name>_prodigal.faa."""
    shim = _import_shim()
    proteome = shim.create_proteome(str(BINS / "metator_00082_00000.fa"))
    assert proteome == "metator_00082_00000_prodigal.faa"
    headers = [line for line in Path(proteome).read_text().splitlines() if line.startswith(">")]
    assert len(headers) == 46
    # Gene, start, end and strand.
    assert headers[0].split(" # ")[:4] == [">Contig_200577_1", "1", "975", "-1"]
    assert headers[-1].split(" # ")[:4] == [">Contig_486059_3", "1188", "1664", "-1"]
    # pyrodigal warns on this bin, under 100 kb, as prodigal does, but on a
    # stderr that miComplete discards: the warning is not shown.
    assert not [warning for warning in recwarn if issubclass(warning.category, UserWarning)]


@pytest.mark.parametrize("length", [0, 15_000])
def test_create_proteome_empty_under_20kb(shim_state, tmp_path, length):
    """prodigal needs 20 kb to train on: with less, or no contig at all, it
    leaves the proteome file empty."""
    shim = _import_shim()
    fasta = _first_contigs(tmp_path / "bin.fa", length)
    assert shim.create_proteome(str(fasta), "small") == "small_prodigal.faa"
    assert Path("small_prodigal.faa").read_text() == ""


# hmmsearch -E 4e-10 --tblout (HMMER 3.4) on the proteins of bin
# metator_00063_00000: target, marker, E-value, score, bias and E-value of the
# best domain, the columns of the table that miComplete reads.
HMMSEARCH_HITS = {
    "Bact105": [
        ("Contig_314115_2", "Ribosomal_L14", "5.3e-54", "173.3", "0.4", "5.8e-54"),
        ("Contig_140928_1", "Ribosomal_L16", "8.4e-53", "169.8", "3.7", "9.3e-53"),
        ("Contig_140928_2", "Ribosomal_L29", "8.6e-28", "88.0", "6.9", "1.1e-27"),
        ("Contig_314115_4", "Ribosomal_L5_C", "1.3e-41", "132.4", "0.0", "2e-41"),
        ("Contig_314115_1", "Ribosomal_S17", "1.5e-31", "100.4", "3.1", "1.7e-31"),
        ("Contig_140928_3", "Ribosomal_S17", "2.3e-12", "38.9", "3.6", "2.7e-12"),
        ("Contig_314115_5", "Ribosomal_S8", "8.5e-38", "121.2", "0.0", "9.4e-38"),
        ("Contig_419692_1", "tRNA-synt_1d", "8.3e-85", "276.6", "0.1", "1.2e-84"),
        ("Contig_522952_1", "TruB_N", "4.1e-49", "158.6", "0.0", "4.9e-49"),
    ],
    "Arch131": [
        ("Contig_314115_2", "Ribosomal_L14", "5.3e-54", "173.3", "0.4", "5.8e-54"),
        ("Contig_140928_1", "Ribosomal_L16", "8.4e-53", "169.8", "3.7", "9.3e-53"),
        ("Contig_314115_4", "Ribosomal_L5", "1.2e-26", "84.4", "0.1", "3.4e-26"),
        ("Contig_314115_1", "Ribosomal_S17", "1.5e-31", "100.4", "3.1", "1.7e-31"),
        ("Contig_140928_3", "Ribosomal_S17", "2.3e-12", "38.9", "3.6", "2.7e-12"),
        ("Contig_314115_5", "Ribosomal_S8", "8.5e-38", "121.2", "0.0", "9.4e-38"),
        ("Contig_419692_1", "tRNA-synt_1d", "8.3e-85", "276.6", "0.1", "1.2e-84"),
        ("Contig_140928_2", "Ribosomal_L29", "8.6e-28", "88.0", "6.9", "1.1e-27"),
        ("Contig_522952_1", "TruB_N", "4.1e-49", "158.6", "0.0", "4.9e-49"),
    ],
}


@pytest.mark.parametrize("markers", ["Bact105", "Arch131"])
def test_hmm_search_finds_the_hits_of_hmmsearch(shim_state, markers):
    shim = _import_shim()
    from micomplete.completeness import calcCompleteness

    proteome = shim.create_proteome(str(BINS / "metator_00063_00000.fa"))
    # The e-value of miComplete's command line.
    search = calcCompleteness(proteome, "bin", _marker_hmms(markers), evalue=4e-10)
    assert search.hmm_search() == ("bin.tblout", 0)
    hits = [line.split() for line in Path("bin.tblout").read_text().splitlines() if not line.startswith("#")]
    assert sorted(tuple(hit[i] for i in (0, 2, 4, 5, 6, 7)) for hit in hits) == sorted(HMMSEARCH_HITS[markers])


def test_hmm_search_refuses_an_empty_proteome(shim_state, caplog, capsys):
    """As hmmsearch, which fails on a proteome without any protein and writes
    no table: miComplete warns and counts no marker."""
    _import_shim()
    from micomplete.completeness import calcCompleteness

    Path("bin.faa").touch()
    logger = logging.getLogger("metator-shim-test")
    search = calcCompleteness("bin.faa", "bin", _marker_hmms("Bact105"), logger=logger)
    with caplog.at_level(logging.WARNING, logger=logger.name):
        assert search.hmm_search() == ("bin.tblout", 1)
        assert search.get_completeness() == (0, 0, 0)
    assert "is bin.faa empty?" in caplog.text
    assert not Path("bin.tblout").exists()
    # Without a logger, the warning goes to stderr.
    assert calcCompleteness("bin.faa", "bin", _marker_hmms("Bact105")).hmm_search() == ("bin.tblout", 1)
    assert "Warning: Error thrown by HMMER, is bin.faa empty?" in capsys.readouterr().err


def test_micomplete_finds_prodigal_and_hmmsearch(shim_state, tmp_path, monkeypatch):
    """miComplete stops if prodigal or hmmsearch is not on the PATH: its
    which() finds them, as their stand-ins, and other programs as
    shutil.which does. shutil itself is left alone."""
    shim = _import_shim()
    import micomplete.micomplete

    monkeypatch.setenv("PATH", str(tmp_path))
    cli_shutil = micomplete.micomplete.shutil
    assert shutil.which("prodigal") is None and shutil.which("hmmsearch") is None
    assert cli_shutil.which("prodigal") and cli_shutil.which("hmmsearch")
    assert cli_shutil.which("sh") is None
    assert cli_shutil.which("sh", path=os.defpath) == shutil.which("sh", path=os.defpath) is not None
    assert cli_shutil.rmtree is shutil.rmtree
    assert shutil.which is not shim.which


# The rows of miComplete for three bins of tests_data/overlapping_bin, and two
# bins of the first contigs of metator_00009_00000, as prodigal 2.6.3 and
# HMMER 3.4 gave them. too_short holds 13 351 bp, too little for prodigal to
# train on; just_enough holds 19 914 bp, which is enough once prodigal joins
# its 9 contigs with linkers. Columns: Length, GC-content, Present Markers,
# Completeness, Redundancy, Weighted completeness, Weighted redundancy,
# Contigs, N50, L50, N90, L90, CDs.
WITH_THE_PROGRAMS = {
    "Bact105": {
        "metator_00002_00000": "3361245 38.39 104 0.9905 1.8365 0.9907 1.5089 687 13185 72 1950 325 4026",
        "metator_00063_00000": "38265 41.83 8 0.0762 1.1250 0.0284 1.0009 26 2232 7 588 20 55",
        "metator_00082_00000": "36558 42.09 2 0.0190 1.5000 0.0224 1.5847 13 5301 2 1039 8 46",
        "too_short": "13351 42.49 0 0 0 0 0 5 4609 2 1294 3 0",
        "just_enough": "19914 42.62 0 0.0000 0 0 0 9 4609 2 1018 7 29",
    },
    "Arch131": {
        "metator_00002_00000": "3361245 38.39 56 0.4275 2.1429 0.3657 1.6738 687 13185 72 1950 325 4026",
        "metator_00063_00000": "38265 41.83 8 0.0611 1.1250 0.0457 1.0025 26 2232 7 588 20 55",
        "metator_00082_00000": "36558 42.09 1 0.0076 2.0000 0.0073 2.0 13 5301 2 1039 8 46",
        "too_short": "13351 42.49 0 0 0 0 0 5 4609 2 1294 3 0",
        "just_enough": "19914 42.62 1 0.0076 1.0000 0.0024 1.0 9 4609 2 1018 7 29",
    },
}


@pytest.mark.parametrize("markers", ["Bact105", "Arch131"])
def test_bins_validated_as_with_the_programs(tmp_path, markers):
    """`python -m metator._micomplete`, as validation runs it, gives the rows
    that prodigal and hmmsearch gave, with stand-ins of both first on the
    PATH, which fail, and leave a mark, if anything runs them."""
    stubs = tmp_path / "stubs"
    stubs.mkdir()
    for program in ("prodigal", "hmmsearch"):
        (stubs / program).write_text(f"#!/bin/sh\ntouch '{tmp_path}/{program} was run'\nexit 1\n")
        (stubs / program).chmod(0o755)
    fastas = [BINS / f"{name}.fa" for name in ("metator_00002_00000", "metator_00063_00000", "metator_00082_00000")]
    fastas += [_first_contigs(tmp_path / "too_short.fa", 15_000), _first_contigs(tmp_path / "just_enough.fa", 20_000)]
    seq_tab = tmp_path / "seq.tab"
    seq_tab.write_text("".join(f"{fasta}\tfna\n" for fasta in fastas))
    outfile = tmp_path / f"{markers}.tsv"
    cmd = [sys.executable, "-m", SHIM, str(seq_tab), "--hmms", markers, "--weights", markers, "--threads", "2",
           "--outfile", str(outfile)]  # fmt: skip
    env = {**os.environ, "PATH": f"{stubs}{os.pathsep}{os.environ.get('PATH', '')}"}
    subprocess.run(cmd, cwd=tmp_path, env=env, check=True)

    rows = {}
    for line in outfile.read_text().splitlines():
        if not line.startswith(("#", "Name\t")):
            name, *fields = line.rstrip("\t").split("\t")
            rows[name] = fields
    assert rows == {name: row.split() for name, row in WITH_THE_PROGRAMS[markers].items()}
    assert not list(tmp_path.glob("* was run"))
