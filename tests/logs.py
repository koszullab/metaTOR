"""The log files of the runs: <name>.metator_<date>.log, named after the map of
the run, in its output directory (the current directory for mge)."""

from pathlib import Path


def logs(directory, command=None):
    """The logs in directory, oldest first; with command, those of its runs only,
    by the command line that their header records (## cmd: metator <command>)."""
    found = sorted(Path(directory).glob("*.metator_*.log"), key=lambda log: (log.stat().st_mtime_ns, log.name))
    if command is None:
        return found
    return [log for log in found if f"## cmd: metator {command} " in log.read_text()]
