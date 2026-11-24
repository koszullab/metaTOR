# Installation

metaTOR needs Python 3.10 to 3.14, and nothing outside PyPI.

## With conda (recommended)

```bash
conda create -n metator bioconda::metator
conda activate metator
```

!!! tip "Check it worked"

    ```bash
    metator --version
    ```

## With pip

```bash
pip install metator
```

pip installs metaTOR's Python dependencies, among them
[metacooler](https://github.com/js2264/metacooler),
[cooler](https://github.com/open2c/cooler),
[miComplete](https://pypi.org/project/micomplete/),
[pyrodigal](https://pypi.org/project/pyrodigal/) and
[pyhmmer](https://pypi.org/project/pyhmmer/). Install it in its own
environment, so that it does not change the versions other tools rely on.

To make the contact map from the reads, install
[hicstuff](https://koszullab.github.io/hicstuff/) and its aligner (bowtie2 or
bwa) separately. metaTOR can also make the map from pairs or BAM files: see
[Quickstart](quickstart.md).

## With Docker

Images are published to the GitHub Container Registry for each release:

```bash
docker pull ghcr.io/koszullab/metator:latest
docker pull ghcr.io/koszullab/metator:v<version>
```

## Development version

The development version is managed with [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/koszullab/metator.git
cd metator
uv sync
uv run metator --help
```

`uv sync --extra docs` adds the tools that build this documentation; see
[Contributing](contributing.md).

## Licence

metaTOR is distributed under the GNU General Public License v3 (GPLv3).
