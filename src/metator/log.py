#!/usr/bin/env python3
# coding: utf-8

"""Basic logging setup for metator

Basic logging setup with two handlers (stdout and log file). By default the
log file is disabled, but can be enabled or changed using set_file_handler.
Python's warnings (warnings.warn, which cooler and pandas give theirs with)
go to the same handlers, as metaTOR's own messages do.
"""

import logging
import logging.handlers


CURRENT_LOG_LEVEL = logging.INFO


logging.captureWarnings(True)

logger = logging.getLogger("metator_logger")
logger.setLevel(CURRENT_LOG_LEVEL)

# The logger that captureWarnings sends Python's warnings to: without a
# handler, logging gives it one that drops them (hicstuff, which metaTOR no
# longer imports, gave the root logger one that printed them).
warnings_logger = logging.getLogger("py.warnings")

logfile_formatter = logging.Formatter("%(asctime)s :: %(levelname)s :: %(message)s", datefmt="%Y-%m-%d,%H:%M:%S")
stdout_formatter = logging.Formatter("%(levelname)s :: %(message)s")

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(stdout_formatter)
for _logger in (logger, warnings_logger):
    _logger.addHandler(stream_handler)
    _logger.propagate = False


def set_file_handler(log_path, formatter=logfile_formatter):
    """Change the file handler for custom log file location: that of metaTOR's
    messages and of Python's warnings."""

    filehandler = logging.FileHandler(log_path, "a")
    filehandler.setLevel(logging.INFO)
    filehandler.setFormatter(formatter)
    for _logger in (logger, warnings_logger):
        for hdlr in _logger.handlers[:]:  # remove the existing file handlers
            if isinstance(hdlr, logging.FileHandler):
                _logger.removeHandler(hdlr)
        _logger.addHandler(filehandler)  # set the new handler
