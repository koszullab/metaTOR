# Installation

metaTOR needs Python 3.10 to 3.14, and two programs that are not on PyPI:
the prodigal and hmmer that miComplete runs to validate the clusters. It
starts from a contact map of the contigs, which
[hicstuff](https://koszullab.github.io/hicstuff/) makes from the reads
(`hicstuff pipeline -M metacool`), or
[metacooler](https://github.com/js2264/metacooler) from pairs files
(`metacooler create pairs`): see [Quickstart](quickstart.md).

## With conda (recommended)

bioconda installs metaTOR with the programs it runs:

```bash
conda create -n metator bioconda::metator
conda activate metator
```

!!! tip "Check it worked"

    ```bash
    metator --version
    ```

## With pip

Install the programs first, then metaTOR:

```bash
# Debian or Ubuntu; on macOS: brew install hmmer prodigal
sudo apt install hmmer prodigal

pip install metator
```

| Program | Used for |
|---|---|
| `prodigal`, `hmmer` | Finding the marker genes of the clusters, with miComplete |

pip installs metaTOR's Python dependencies, among them
[metacooler](https://github.com/js2264/metacooler),
[cooler](https://github.com/open2c/cooler) and
[miComplete](https://pypi.org/project/micomplete/); install metaTOR in an
environment of its own so that they do not change the versions other tools
rely on. hicstuff, which makes the map from the reads, is not one of them:
install it with its aligner, bowtie2 or bwa, where the reads are aligned (see
[hicstuff's documentation](https://koszullab.github.io/hicstuff/)).

## With Docker

Docker images, with metaTOR and the programs it runs, are published to the
GitHub Container Registry when a release is tagged:

```bash
docker pull ghcr.io/koszullab/metator:latest
docker pull ghcr.io/koszullab/metator:v<version>
```

## Development version

The development version is managed with [uv](https://docs.astral.sh/uv/).
With prodigal and hmmer installed:

```bash
git clone https://github.com/koszullab/metator.git
cd metator
uv sync
uv run metator --help
```

`uv sync --extra docs` adds what builds this documentation; see
[Contributing](contributing.md).

## Licence

metaTOR is distributed under the GNU General Public License v3 (GPLv3).
