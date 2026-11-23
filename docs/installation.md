# Installation

metaTOR needs Python 3.10 to 3.14, and two programs that are not on PyPI:
an aligner (bwa, recommended, or bowtie2) and samtools.

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
# Debian or Ubuntu; on macOS: brew install bowtie2 bwa samtools
sudo apt install bowtie2 bwa samtools

pip install metator
```

| Program | Used for |
|---|---|
| `bwa` | Aligning the reads with `--aligner bwa`, the recommended aligner |
| `bowtie2` | Aligning the reads with `--aligner bowtie2`, the default |
| `samtools` | Handling the alignments |

Only the aligner a run uses is needed: `--aligner bwa` needs bwa, not
bowtie2. A run that starts from BAM files, a pairs file or a contact map, with
the assembly given as a FASTA file, needs no aligner. pip installs metaTOR's Python dependencies, among them
[hicstuff](https://koszullab.github.io/hicstuff/),
[metacooler](https://github.com/js2264/metacooler),
[cooler](https://github.com/open2c/cooler),
[miComplete](https://pypi.org/project/micomplete/), and
[pyrodigal](https://pypi.org/project/pyrodigal/) and
[pyhmmer](https://pypi.org/project/pyhmmer/), with which miComplete runs
Prodigal and HMMER to validate the bins; install metaTOR in an
environment of its own so that they do not change the versions other tools
rely on.

## With Docker

Docker images, with metaTOR and the programs it runs, are published to the
GitHub Container Registry when a release is tagged:

```bash
docker pull ghcr.io/koszullab/metator:latest
docker pull ghcr.io/koszullab/metator:v<version>
```

## Development version

The development version is managed with [uv](https://docs.astral.sh/uv/).
With the programs above installed:

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
