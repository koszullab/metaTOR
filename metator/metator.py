#!/usr/bin/env python3

"""metaTOR - a set of scripts that streamlines the processing and binning of
metagenomic 3C datasets.
"""


import subprocess
import sys
import pkg_resources
import pathlib
import requests
import tarfile

BASE_PRODIGAL = (
    "https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal"
)
HMM_URL = "http://dl.pasteur.fr/fop/LItxiFe9/hmm_databases.tgz"


def download_and_install_dependencies():
    """Setup URLS and download dependencies.
    """

    dependencies = {"hmm_databases": HMM_URL}

    if sys.platform.startswith("linux") or "bsd" in sys.platform:

        dependencies["prodigal"] = "{}.linux".format(BASE_PRODIGAL)
        dependencies["louvain"] = (
            "https://lip6.github.io/Louvain-BinaryBuild/"
            "louvain_linux.tar.gz"
        )
    elif sys.platform == "darwin":

        dependencies["prodigal"] = "{}.osx.10.9.5".format(BASE_PRODIGAL)
        dependencies["louvain"] = (
            "https://github.com/lip6/Louvain-BinaryBuilds/raw/osx/"
            "louvain_osx.tar.gz"
        )
    elif sys.platform.startswith("win") or sys.platform == "cygwin":

        dependencies["prodigal"] = "{}.windows.exe"
        dependencies["louvain"] = (
            "https://ci.appveyor.com/api/projects/yanntm/"
            "Louvain-BinaryBuild/artifacts/website/"
            "louvain_windows.tar.gz"
        )

    else:
        raise NotImplementedError(
            "Your platform is not supported: {}".format(sys.platform)
        )

    cache_dir = pathlib.Path.cwd() / pathlib.Path("cache")

    try:
        print("Downloading dependencies...")
        cache_dir.mkdir()
        for dependency_name, url in dependencies.items():
            print("Downloading {} at {}".format(dependency_name, url))
            request = requests.get(url)
            basename = url.split("/")[-1]

            with open(cache_dir / basename, "wb") as handle:
                print(dependency_name, basename, cache_dir / basename)
                handle.write(request.content)
    except FileExistsError:
        print("Using cached dependencies...")

    share_dir = pathlib.Path.cwd()
    tools_dir = share_dir / "tools"
    louvain_dir = tools_dir / "louvain"

    louvain_dir.mkdir(parents=True, exist_ok=True)

    louvain_basename = dependencies["louvain"].split("/")[-1]
    louvain_path = louvain_dir / louvain_basename
    (cache_dir / louvain_basename).replace(louvain_path)

    with tarfile.open(louvain_path, "r:gz") as tar:
        tar.extractall()

    hmm_basename = dependencies["hmm_databases"].split("/")[-1]
    hmm_path = share_dir / hmm_basename
    (cache_dir / hmm_basename).replace(hmm_path)

    prodigal_basename = dependencies["prodigal"].split("/")[-1]
    prodigal_path = tools_dir / "prodigal"
    (cache_dir / prodigal_basename).replace(prodigal_path)


def main():
    """This module just acts as an entry point to the bulk of the pipeline.
    All argument parsing is delegated to metator.sh
    """

    metator_args = sys.argv[1:]
    entry_point = pkg_resources.resource_filename("metator", "bin/metator.sh")

    metator_process = subprocess.Popen((entry_point, *metator_args))
    metator_process.wait()


if __name__ == "__main__":
    main()
