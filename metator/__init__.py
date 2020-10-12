#!/usr/bin/env python3

from .version import __version__ as version
from .scripts import bins, fasta_utils, figures, hicstuff, log, network, hamming

__author__ = "Lyam Baudry, Théo Foutel-Rodier, Martial Marbouty"
__copyright__ = "Copyright © 2017-2018, Institut Pasteur, Paris, France"
__credits__ = [
    "Lyam Baudry",
    "Théo Foutel-Rodier",
    "Martial Marbouty",
    "Axel Cournac",
    "Vittore Scolari",
    "Romain Koszul",
]
__license__ = "GPLv3"
__maintainer__ = "Lyam Baudry"
__email__ = "lyam.baudry@pasteur.fr"
__status__ = "Alpha"
__version__ = version
