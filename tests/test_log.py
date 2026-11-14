# Tests of the logging helpers.

import logging

import pytest

import metator.log as mtl


@pytest.fixture
def restore_handlers():
    """Put the metator logger handlers back as they were."""
    saved = mtl.logger.handlers[:]
    yield
    for handler in mtl.logger.handlers:
        if handler not in saved:
            handler.close()
    mtl.logger.handlers[:] = saved


def test_set_file_handler(tmp_path, restore_handlers):
    """A new log file replaces the previous one."""
    mtl.set_file_handler(str(tmp_path / "first.log"))
    mtl.set_file_handler(str(tmp_path / "second.log"))
    files = [h.baseFilename for h in mtl.logger.handlers if isinstance(h, logging.FileHandler)]
    assert files == [str(tmp_path / "second.log")]
    mtl.logger.info("hello")
    assert (tmp_path / "second.log").read_text().endswith(" :: INFO :: hello\n")
    assert (tmp_path / "first.log").read_text() == ""
