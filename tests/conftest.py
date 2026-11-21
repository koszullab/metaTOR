import logging
import os
import shutil
from pathlib import Path

import matplotlib

# Figures are only written to files: never open a window.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from tests import synthetic
from tests.ansi import plain

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"

# Seed of the mock metagenome shared by the synthetic-data tests. Set
# METATOR_TEST_SEED to check that they do not depend on one lucky draw.
SYNTHETIC_SEED = int(os.environ.get("METATOR_TEST_SEED", "20251114"))


@pytest.fixture(autouse=True)
def _run_in_tmp_path(tmp_path, monkeypatch):
    """Run each test from its own temporary directory.

    Tests read their inputs as "tests_data/...": a symlink keeps those paths
    valid, while everything a test or a command writes to the working
    directory (logs, temporary directories, outputs) stays out of the
    repository.
    """
    os.symlink(TESTS_DATA, tmp_path / "tests_data")
    monkeypatch.chdir(tmp_path)


@pytest.fixture
def private_link(tmp_path):
    """Make links to files of tests_data, each test its own, and check at
    teardown that no file was written to through them.

    pyfastx writes its index (.fxi) next to the path it is given, and the
    commands remove it at the end, so two test sessions reading one
    assembly of tests_data removed it from under each other. The links are
    in a directory of their own: in tmp_path itself, a link would take the
    name of a file the test writes there (a copy named assembly.fa), which
    would then be written through the link, into tests_data."""
    directory = tmp_path / "private_links"
    targets = {}

    def link(path):
        target = Path(path).resolve()
        directory.mkdir(exist_ok=True)
        new = directory / Path(path).name
        new.symlink_to(target)
        targets[target] = target.stat()
        return str(new)

    yield link
    for target, before in targets.items():
        after = target.stat()
        assert (after.st_size, after.st_mtime_ns) == (before.st_size, before.st_mtime_ns), (
            f"{target} was written to through its link"
        )


@pytest.fixture(autouse=True)
def _close_figures():
    """Close the figures a test leaves open.

    Some plotting functions draw on pyplot's current figure (and some never
    close theirs): without this, a figure left open by one test ends up in
    the plots of the next one.
    """
    yield
    plt.close("all")


@pytest.fixture
def saved_figures(monkeypatch):
    """The figures saved during a test, by file name.

    The plotting functions close the figures they save: this keeps them,
    so that a test can check what they show.
    """
    from matplotlib.figure import Figure

    saved = {}
    savefig = Figure.savefig

    def record(figure, fname, *args, **kwargs):
        saved[Path(fname).name] = figure
        return savefig(figure, fname, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", record)
    return saved


class _PlainResult:
    """A click Result whose output has no ANSI escape codes (see tests.ansi)."""

    def __init__(self, result):
        self._result = result
        self.output = plain(result.output)

    def __getattr__(self, name):
        return getattr(self._result, name)


@pytest.fixture
def run_cli(monkeypatch):
    """Run a metator command line in-process, as typed after ``metator``.

    ``sys.argv`` is set as the shell would set it (the log files record the
    command line from it). An exception raised by the command propagates;
    a usage error (exit code 2) is returned as the click result, whose exit
    code is checked against ``exit_code``.
    """
    import sys

    from click.testing import CliRunner

    from metator.cli import cli

    def run(*args, exit_code=0):
        args = [str(arg) for arg in args]
        monkeypatch.setattr(sys, "argv", ["metator", *args])
        # As wide as the help gets (100 columns): short error messages do not wrap.
        result = _PlainResult(CliRunner(env={"COLUMNS": "200"}).invoke(cli, args, catch_exceptions=False))
        assert result.exit_code == exit_code, result.output
        return result

    return run


@pytest.fixture
def network_from_a_map(monkeypatch):
    """Have a pipeline run stop once its network is built from its .metacool
    map: the partition raises RuntimeError("network built"). Returns the
    maps that metacool_to_contacts read, in order."""
    import metator.network as mtn

    read = []
    metacool_to_contacts = mtn.metacool_to_contacts

    def from_the_map(metacool, *args, **kwargs):
        read.append(metacool)
        return metacool_to_contacts(metacool, *args, **kwargs)

    def partition(*args):
        raise RuntimeError("network built")

    monkeypatch.setattr("metator.network.metacool_to_contacts", from_the_map)
    monkeypatch.setattr("metator.partition.partition", partition)
    return read


@pytest.fixture(scope="session")
def metagenome():
    """The mock community (see tests/synthetic.py) and its ground truth."""
    return synthetic.make_metagenome(seed=SYNTHETIC_SEED)


@pytest.fixture(scope="session")
def metagenome_files(metagenome, tmp_path_factory):
    """Assembly, depth.txt, pairs and truth of the mock community.

    Shared by the whole session: tests must not write next to these files
    (use ``metagenome_inputs`` for commands that do), which is checked at
    the end of the session.
    """
    directory = tmp_path_factory.mktemp("metagenome")
    files = metagenome.write_all(directory)
    yield files
    written = sorted(set(os.listdir(directory)) - {path.name for path in files.values()})
    assert not written, f"tests wrote next to the shared mock metagenome: {written}"


@pytest.fixture(scope="session")
def metagenome_indexed_pairs(metagenome_files, tmp_path_factory):
    """The pairs sorted, compressed and indexed for pypairix, as metator mge
    reads them: by chr1, chr2, pos1 and pos2 (the contig names byte by
    byte, as pairtools sorts them, so that each forms one block; the
    positions as numbers), with "#sorted: chr1-chr2-pos1-pos2" after the
    "## pairs format" line, compressed with BGZF (Biopython) and indexed
    with pypairix."""
    import pypairix
    from Bio import bgzf

    header, rows = [], []
    with open(metagenome_files["pairs"]) as pairs:
        for line in pairs:
            (header if line.startswith("#") else rows).append(line)

    def key(line):
        fields = line.split("\t")
        return fields[1].encode(), fields[3].encode(), int(fields[2]), int(fields[4])

    rows.sort(key=key)
    header = [line for line in header if not line.startswith("#sorted")]
    header.insert(1 if header and header[0].startswith("##") else 0, "#sorted: chr1-chr2-pos1-pos2\n")
    path = tmp_path_factory.mktemp("metagenome_pairs") / "synthetic_sorted.pairs.gz"
    with bgzf.BgzfWriter(str(path), "wb") as out:
        out.write("".join(header + rows).encode())
    pypairix.build_index(str(path), force=1)
    return path


@pytest.fixture(scope="session")
def metagenome_network(metagenome_files, tmp_path_factory):
    """Network (empirical_hit) and contig data of the mock community, read
    from the map of its pairs (see tests.maps.contig_map)."""
    import metator.network as mtn
    from tests.maps import contig_map

    directory = tmp_path_factory.mktemp("metagenome_network")
    tmp_dir = directory / "tmp"
    tmp_dir.mkdir()
    contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
    metacool = contig_map(metagenome_files["pairs"], contig_data, directory / "contigs.metacool", columns=False)
    network, contig_data_file = mtn.metacool_to_contacts(
        str(metacool),
        contig_data,
        str(directory),
        "network.txt",
        "contig_data_network.txt",
        str(tmp_dir),
        1,
        "empirical_hit",
    )
    return {"network": Path(network), "contig_data": Path(contig_data_file)}


@pytest.fixture
def metagenome_inputs(metagenome_files, tmp_path):
    """Private copies of the mock community files, in ``tmp_path/inputs``."""
    directory = tmp_path / "inputs"
    directory.mkdir()
    copies = {}
    for key, path in metagenome_files.items():
        copies[key] = directory / path.name
        shutil.copyfile(path, copies[key])
    return copies


@pytest.fixture
def logged():
    """The messages that metaTOR logs during the test, with their level."""
    from metator.log import logger

    records = []
    handler = logging.Handler()
    handler.emit = records.append
    logger.addHandler(handler)
    yield lambda: [(record.levelname, record.getMessage()) for record in records]
    logger.removeHandler(handler)
