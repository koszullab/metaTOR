"""Locks of maps as metacooler writes them, held by another process, as by another
run of metaTOR: alive, or killed while it held it."""

import json
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

# A process that holds the lock of a map until its standard input closes.
HOLD = (
    "import sys, metacooler\n"
    "with metacooler.lock(sys.argv[1], command=sys.argv[2]):\n"
    "    print('held', flush=True)\n"
    "    sys.stdin.read()\n"
)
# A process killed while it held the lock of a map: it ends without leaving
# the with block that would remove it (the context manager kept, as its
# generator, collected, would release the lock).
KILLED = "import os, sys, metacooler\nheld = metacooler.lock(sys.argv[1], command=sys.argv[2])\nheld.__enter__()\nos._exit(0)\n"


def lock_file(metacool):
    """The lock file of the map metacool: next to its file (next to the target of
    a symbolic link)."""
    return Path(f"{Path(str(metacool).split('::', 1)[0]).resolve()}.lock")


@contextmanager
def held_lock(metacool, command="metator pipeline"):
    """Hold the lock of the map metacool in another process for the with block:
    yields what its lock file says of the process (pid, host, start, command,
    uuid...)."""
    process = subprocess.Popen(
        [sys.executable, "-c", HOLD, str(metacool), command], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
    )
    stdin, stdout = process.stdin, process.stdout
    assert stdin is not None and stdout is not None
    try:
        assert stdout.readline().strip() == "held"
        yield json.loads(lock_file(metacool).read_text())
    finally:
        stdin.close()
        process.wait(timeout=60)


def killed_lock(metacool, command="metator pipeline"):
    """The lock of the map metacool left by a process that held it, and is gone:
    what its lock file says of the process."""
    subprocess.run([sys.executable, "-c", KILLED, str(metacool), command], check=True, timeout=60)
    return json.loads(lock_file(metacool).read_text())


def described(holder):
    """A holder of a lock, as metacooler's messages name it."""
    return f"pid {holder['pid']} on {holder['host']}, since {holder['start']}: {holder['command']}"


def refusal(metacool, holder):
    """metacooler's refusal of the lock of the map metacool, which a live process
    of this host holds (holder: what its lock file says)."""
    path = lock_file(metacool).with_suffix("")
    return (
        f"{path} is being written by another process ({described(holder)}): wait for it to end, or, if no such "
        f"process is running, delete {lock_file(metacool)} and {path.parent}/.{path.name}.{holder['uuid']}.*.d/ "
        "(its staged copies)."
    )
