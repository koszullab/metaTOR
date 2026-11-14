import os
import shutil
from pathlib import Path

import matplotlib

# Figures are only written to files: never open a window.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from tests import synthetic

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


@pytest.fixture(autouse=True)
def _close_figures():
    """Close the figures a test leaves open.

    Some plotting functions draw on pyplot's current figure (and some never
    close theirs): without this, a figure left open by one test ends up in
    the plots of the next one.
    """
    yield
    plt.close("all")


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
    """The pairs sorted, bgzipped and indexed for pypairix."""
    if shutil.which("pairtools") is None:
        pytest.skip("pairtools is not installed")
    import metator.io as mio

    directory = tmp_path_factory.mktemp("metagenome_pairs")
    pairs = directory / "synthetic.pairs"
    shutil.copyfile(metagenome_files["pairs"], pairs)
    return Path(mio.sort_pairs_pairtools(str(pairs), tmp_dir=str(directory)))


@pytest.fixture(scope="session")
def metagenome_network(metagenome_files, tmp_path_factory):
    """Network (empirical_hit) and contig data of the mock community."""
    import metator.network as mtn

    directory = tmp_path_factory.mktemp("metagenome_network")
    tmp_dir = directory / "tmp"
    tmp_dir.mkdir()
    contig_data, hit_data = mtn.create_contig_data(str(metagenome_files["assembly"]), 1, str(metagenome_files["depth"]), "DpnII")
    network, contig_data_file = mtn.alignment_to_contacts(
        [str(metagenome_files["pairs"])],
        contig_data,
        0,
        hit_data,
        str(directory),
        "network.txt",
        "contig_data_network.txt",
        str(tmp_dir),
        1,
        "empirical_hit",
        False,
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
