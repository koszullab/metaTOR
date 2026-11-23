# Tests of the logging helpers.

import logging
import subprocess
import sys

import pytest

import metator.log as mtl


@pytest.fixture
def restore_handlers():
    """Put the handlers of the metator logger, and of Python's warnings, back as
    they were."""
    saved = {log: log.handlers[:] for log in (mtl.logger, mtl.warnings_logger)}
    yield
    for log, handlers in saved.items():
        for handler in log.handlers:
            if handler not in handlers:
                handler.close()
        log.handlers[:] = handlers


def test_set_file_handler(tmp_path, restore_handlers):
    """A new log file replaces the previous one."""
    mtl.set_file_handler(str(tmp_path / "first.log"))
    mtl.set_file_handler(str(tmp_path / "second.log"))
    files = [h.baseFilename for h in mtl.logger.handlers if isinstance(h, logging.FileHandler)]
    assert files == [str(tmp_path / "second.log")]
    mtl.logger.info("hello")
    assert (tmp_path / "second.log").read_text().endswith(" :: INFO :: hello\n")
    assert (tmp_path / "first.log").read_text() == ""
    assert mtl.warnings_logger.handlers == [mtl.stream_handler, *[h for h in mtl.logger.handlers if h is not mtl.stream_handler]]


def test_python_warnings_are_shown_and_logged(tmp_path):
    """Python's warnings (warnings.warn, which cooler gives the problems of a map
    with) are shown, and written to the log file of the run, as metaTOR's own
    messages are: without a handler, logging dropped them."""
    log = tmp_path / "run.log"
    code = (
        "import warnings\n"
        "import metator.cli\n"
        "import metator.log as mtl\n"
        f"mtl.set_file_handler({str(log)!r})\n"
        "warnings.warn('a problem of the map')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stderr.startswith("WARNING :: ") and "UserWarning: a problem of the map" in result.stderr
    logged = log.read_text()
    assert " :: WARNING :: " in logged and "UserWarning: a problem of the map" in logged


@pytest.mark.parametrize(
    "seconds, written",
    [(0, "0h 0m 0s"), (0.6, "0h 0m 1s"), (59.4, "0h 0m 59s"), (1051, "0h 17m 31s"), (3 * 3600 + 61, "3h 1m 1s")],
)
def test_duration(seconds, written):
    """A duration as hicstuff writes it at the end of its log."""
    assert mtl.duration(seconds) == written
