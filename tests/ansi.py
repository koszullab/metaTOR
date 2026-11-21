"""Plain text of what rich-click prints."""

import re

# CSI sequences (colours, bold...) and OSC sequences (hyperlinks).
_ANSI = re.compile(r"\x1b(?:\[[0-?]*[ -/]*[@-~]|\][^\x07\x1b]*(?:\x07|\x1b\\))")


def plain(text):
    """``text`` without its ANSI escape codes."""
    return _ANSI.sub("", text)
