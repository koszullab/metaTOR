import importlib.util
from pathlib import Path
import os
import site


def is_editable_install():
    """Check if the metator package was installed in editable mode."""
    site_packages = site.getsitepackages()
    for site_package in site_packages:
        pth_file = os.path.join(site_package, "_metator.pth")
        if os.path.isfile(pth_file):
            return True
    return False


__metator_source__ = os.path.dirname(importlib.util.find_spec("metator").origin)  # type: ignore
__metator_root__ = __metator_source__
if is_editable_install():
    __metator_root__ = os.path.abspath(os.path.join(__metator_source__, "../../"))
__bin_dir__ = Path(__metator_root__, "bin")
LEIDEN_PATH = str(next(__bin_dir__.glob("networkanalysis-1.3.0*.jar")))
LOUVAIN_PATH = str(__bin_dir__)
