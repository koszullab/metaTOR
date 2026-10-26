#!/usr/bin/env python3
from pathlib import Path
import sys
from .version import __version__ as version
from . import *

__author__ = "Amaury Bignaud, Jacques Serizay, Lyam Baudry, Théo Foutel-Rodier, Martial Marbouty"
__copyright__ = "Copyright © 2017-2025, Institut Pasteur, Paris, France"
__credits__ = [
    "Amaury Bignaud",
    "Jacques Serizay",
    "Lyam Baudry",
    "Théo Foutel-Rodier",
    "Martial Marbouty",
    "Axel Cournac",
    "Vittore Scolari",
    "Romain Koszul",
]
__license__ = "GPLv3"
__maintainer__ = "Amaury Bignaud"
__email__ = "amaury.bignaud@pasteur.fr"
__status__ = "Alpha"
__version__ = version


def find_bin_dir():
    """Find the directory with the louvain binaries and the leiden jar.

    Regular install: next to this file. Editable install: at the repository
    root, or in site-packages/metator/bin.
    """
    here = Path(__file__).parent
    candidates = [here / "bin", here.parents[1] / "bin"]
    candidates += [Path(p, "metator", "bin") for p in sys.path if p]
    for bin_dir in candidates:
        if (bin_dir / "networkanalysis-1.3.0.jar").is_file():
            return bin_dir
    return here / "bin"


__bin_dir__ = find_bin_dir()
LEIDEN_PATH = str(__bin_dir__ / "networkanalysis-1.3.0.jar")
LOUVAIN_PATH = str(__bin_dir__)
