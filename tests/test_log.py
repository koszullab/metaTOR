# Tests of the logging helpers.

import json
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


def test_text_logging(tmp_path, monkeypatch, restore_handlers):
    """Errors, and only errors, are sent to the text notification API."""
    credentials = tmp_path / "text_credentials.json"
    credentials.write_text(json.dumps({"api_service": "https://sms.invalid/send", "user": "me"}))
    sent = []
    monkeypatch.setattr(mtl.requests, "get", lambda url, params: sent.append((url, dict(params))))
    mtl.setup_text_logging(str(credentials))
    mtl.logger.warning("just a warning")
    mtl.logger.error("disk full")
    assert sent == [("https://sms.invalid/send", {"user": "me", "msg": "disk full"})]


def test_text_logging_needs_credentials(tmp_path):
    with pytest.raises(FileNotFoundError):
        mtl.setup_text_logging(str(tmp_path / "missing.json"))
