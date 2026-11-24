# Installation

metaTOR needs Python 3.10 to 3.14, and no program that is not on PyPI. It
starts from a contact map of the contigs, which
[hicstuff](https://koszullab.github.io/hicstuff/) makes from the reads
(`hicstuff pipeline -M metacool`), or
[metacooler](https://github.com/js2264/metacooler) from pairs files
(`metacooler create pairs`), or which `metator network` and `metator pipeline`
make themselves from pairs files or BAM files (read with pysam, without a
samtools program): see [Quickstart](quickstart.md).

## With conda (recommended)

bioconda installs metaTOR:

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
[miComplete](https://pypi.org/project/micomplete/), and
[pyrodigal](https://pypi.org/project/pyrodigal/) and
[pyhmmer](https://pypi.org/project/pyhmmer/), with which miComplete runs
Prodigal and HMMER to validate the clusters; install metaTOR in an
environment of its own so that they do not change the versions other tools
rely on. hicstuff, which makes the map from the reads, is not one of them:
install it with its aligner, bowtie2 or bwa, where the reads are aligned (see
[hicstuff's documentation](https://koszullab.github.io/hicstuff/)).

## With Docker

Docker images of metaTOR are published to the GitHub Container Registry
when a release is tagged:

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

`uv sync --extra docs` adds what builds this documentation; see
[Contributing](contributing.md).

## Licence

metaTOR is distributed under the GNU General Public License v3 (GPLv3).
