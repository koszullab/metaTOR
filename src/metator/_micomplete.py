"""Run miComplete with compatibility patches for recent Biopython and Python.

miComplete (<= 1.1.1, no longer updated upstream) relies on two names that
have since been removed:

    - ``Bio.SeqUtils.GC``, removed in Biopython 1.82 (replaced by
      ``Bio.SeqUtils.gc_fraction``);
    - ``distutils.spawn``, removed from the standard library in Python 3.12.

This module puts both names back before miComplete is imported, so that
metator does not have to pin an old Biopython. The patches are applied at
import time (not only under ``__main__``) so that the worker processes of
miComplete's multiprocessing pool also get them whatever the start method.

The patches only add missing names and are only applied in the miComplete
subprocess: the Python process running metator is left untouched.

Usage::

    python -m metator._micomplete <miComplete arguments>
"""

import shutil
import sys
import types


def patch_biopython_gc():
    """Restore ``Bio.SeqUtils.GC`` if the installed Biopython lacks it."""
    import Bio.SeqUtils

    if hasattr(Bio.SeqUtils, "GC"):
        return

    def GC(seq):
        """GC content in percent, as computed by ``Bio.SeqUtils.GC`` in
        Biopython <= 1.81: G, C and S are counted, and every character
        (including ambiguous ones) counts in the length.
        """
        return 100.0 * Bio.SeqUtils.gc_fraction(seq, ambiguous="ignore")

    Bio.SeqUtils.GC = GC


def patch_distutils_spawn():
    """Provide ``distutils.spawn.find_executable`` if distutils is missing."""
    try:
        from distutils import spawn
    except ImportError:

        def find_executable(executable, path=None):
            return shutil.which(executable, path=path)

        spawn = types.ModuleType("distutils.spawn")
        spawn.find_executable = find_executable
        distutils = types.ModuleType("distutils")
        distutils.spawn = spawn
        sys.modules["distutils"] = distutils
        sys.modules["distutils.spawn"] = spawn


patch_biopython_gc()
patch_distutils_spawn()


def main():
    """Run the miComplete command line interface."""
    from micomplete.micomplete import main as micomplete_main

    sys.argv[0] = "miComplete"
    return micomplete_main()


if __name__ == "__main__":
    sys.exit(main())
