import importlib.util
from pathlib import Path
import os

__metator_source__ = os.path.dirname(importlib.util.find_spec("metator").origin)  # type: ignore
__metator_root__ = os.path.abspath(os.path.join(__metator_source__, "../../"))
__leiden_dir__ = Path(__metator_root__, "external", "artifacts", "networkanalysis", "build", "libs")
LEIDEN_PATH = str(next(__leiden_dir__.glob("networkanalysis-1.3.0*.jar")))
LOUVAIN_PATH = str(Path(__metator_root__, "external", "artifacts", "gen-louvain"))
PAIRIX_PATH = str(Path(__metator_root__, "external", "artifacts", "pairix", "bin", "pairix"))
