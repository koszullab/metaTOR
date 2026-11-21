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

    pyfastx, which metaTOR read the assemblies with, wrote its index (.fxi)
    next to the path it was given, and the commands removed it at the end,
    so two test sessions reading one assembly of tests_data removed it from
    under each other; no command writes beside an assembly now, which the
    check at teardown keeps so. The links are
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


class _Read(list):
    """The maps that a run read the contacts of, in order, the contig data of
    the run, and what its partition was given (see network_from_a_map)."""

    table = None
    partitioned = None


@pytest.fixture
def network_from_a_map(monkeypatch):
    """Have a pipeline run stop once its network is built from its .metacool
    map: the partition raises RuntimeError("network built"). Returns the
    maps whose contacts were read (metator.network.contig_contacts), in
    order, with the contig data of the run (the last that
    metator.network.contig_data_table made) and the network that the
    partition was given (their attribute partitioned)."""
    import metator.network as mtn

    read = _Read()
    contig_contacts = mtn.contig_contacts
    contig_data_table = mtn.contig_data_table

    def from_the_map(metacool, *args, **kwargs):
        read.append(metacool)
        return contig_contacts(metacool, *args, **kwargs)

    def table(*args, **kwargs):
        read.table = contig_data_table(*args, **kwargs)
        return read.table

    def partition(algorithm, network, *args):
        read.partitioned = (read.table, network)
        raise RuntimeError("network built")

    monkeypatch.setattr("metator.network.contig_contacts", from_the_map)
    monkeypatch.setattr("metator.network.contig_data_table", table)
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
def metagenome_network(metagenome_files, tmp_path_factory):
    """The network (empirical_hit) of the mock community, and the data of its
    contigs, as metator network stores them in the map of its pairs (see
    tests.maps.contig_map and metator.metacool.write_map): the map
    (map, which a test copies before it writes it), the network and the
    contig data as the steps read them from the map (network_table, see
    metator.network.map_network, and contig_table, see
    metator.network.map_contig_table), and both written as tables
    (network.txt and contig_data_network.txt: network and contig_data), for
    the tests of the steps that read tables. The IDs are the positions of
    the contigs in the map."""
    import metacooler

    import metator.metacool as mtmc
    import metator.network as mtn
    from tests.maps import contig_map

    directory = tmp_path_factory.mktemp("metagenome_network")
    contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
    metacool = contig_map(metagenome_files["pairs"], contig_data, directory / "contigs.metacool", columns=False)
    contacts = mtn.contig_contacts(str(metacool))
    names = metacooler.MetaCooler(str(metacool)).contigs["name"]
    table = mtn.contig_data_table(contig_data, names, mtn.contig_hits(contacts, len(names)))
    mtmc.write_map(str(metacool), str(metacool), network=mtmc.Network(table, "empirical_hit"), force=True)
    network, contig_table = mtn.map_network(str(metacool)), mtn.map_contig_table(str(metacool))
    network.to_csv(directory / "network.txt", sep="\t", header=False, index=False)
    contig_table.to_csv(directory / "contig_data_network.txt", sep="\t", index=False)
    return {
        "map": metacool,
        "network": directory / "network.txt",
        "contig_data": directory / "contig_data_network.txt",
        "network_table": network,
        "contig_table": contig_table,
    }


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
