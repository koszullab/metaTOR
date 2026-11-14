"""Plain text of what rich-click prints.

rich-click styles its output with ANSI escape codes whenever it believes it
writes to a terminal, which it does as soon as GITHUB_ACTIONS, FORCE_COLOR or
PY_COLORS is set, even under click's CliRunner. The tests read the text, so
they drop the codes first.
"""

import re

# CSI sequences (colours, bold...) and OSC sequences (hyperlinks).
_ANSI = re.compile(r"\x1b(?:\[[0-?]*[ -/]*[@-~]|\][^\x07\x1b]*(?:\x07|\x1b\\))")


def plain(text):
    """``text`` without its ANSI escape codes."""
    return _ANSI.sub("", text)
