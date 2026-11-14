import os
from pathlib import Path

import pytest

TESTS_DATA = Path(__file__).resolve().parents[1] / "tests_data"


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
