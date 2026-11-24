# MetaTOR

[![Bioconda](https://img.shields.io/conda/vn/bioconda/metator?color=blue)](https://anaconda.org/bioconda/metator)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/metator.svg)
[![Build Status](https://github.com/koszullab/metator/actions/workflows/ci.yml/badge.svg)](https://github.com/koszullab/metator/actions)
[![codecov](https://codecov.io/gh/koszullab/metator/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/metator)
[![License: GPLv3](https://img.shields.io/badge/License-GPL%203-0298c3.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Metagenomic Tridimensional Organisation-based Reassembly: metaTOR clusters the
contigs of a metagenome assembly into metagenome-assembled genomes (MAGs) from
the Hi-C contacts between them (metaHiC). It also clusters the contigs of
mobile genetic elements (MGEs) and finds their host MAGs.

<p align="center">
  <img src="https://raw.githubusercontent.com/koszullab/metaTOR/master/docs/images/metator_logo.png" width="200">
</p>

**Documentation: <https://koszullab.github.io/metator/>**:
[quickstart](https://koszullab.github.io/metator/quickstart/),
[the end-to-end pipeline](https://koszullab.github.io/metator/end-to-end/),
[concepts](https://koszullab.github.io/metator/concepts/),
[CLI reference](https://koszullab.github.io/metator/cli/) and
[file formats](https://koszullab.github.io/metator/formats/).

## Installation

```sh
conda create -n metator bioconda::metator   # recommended
pip install metator                          # Python 3.10 to 3.14
```

Docker images are at `ghcr.io/koszullab/metator`. See
[Installation](https://koszullab.github.io/metator/installation/) for the
development version.

## Usage

metaTOR starts from a `.metacool` contact map of the contigs, which
[hicstuff](https://github.com/koszullab/hicstuff) (4.0.0 or later) makes from
the Hi-C reads:

```sh
# The contact map of the contigs: hicstuff/sample.metacool
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz

# The clusters: binning/sample.metacool, the contact map of the clusters
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -d depth.txt -o binning
```

`metator pipeline` also takes pairs or BAM files in place of the map, and with
`--mges mges.txt` it clusters the MGE contigs and finds their hosts. It writes,
in `binning/`:

- `sample.metacool`: the contact map of the clusters, with their quality;
- `fasta/`, `cool/`, `tables/` and `plots/`: the MAGs, their contact maps, the
  cluster summary and other tables, and quality plots;
- `sample.report.html`, `sample.report.txt` and `sample.report_mqc.yaml`: the
  report of the run, for a browser, as text, and for MultiQC;
- `sample.metator_<date>.log`: the log.

| Command | What it does |
|---|---|
| `pipeline` | Runs `network`, `partition` and `validation`; with `--mges`, then `mge` and `host`. |
| `network` | Builds the network of the contigs and stores it in the map. |
| `partition` | Partitions the network many times with Louvain or Leiden. |
| `validation` | Checks the clusters with miComplete, and partitions the contaminated ones again. |
| `mge` | Clusters the contigs of mobile genetic elements. |
| `host` | Finds the MAGs that host the MGE clusters. |
| `dump` | Writes the tables of a run from the map of its clusters. |
| `export` | Writes every output of a run again from the map of its clusters. |
| `report` | Writes the report of a run from the map of its clusters. |

`metator --help` lists the commands, and `metator <command> --help` their
options. The [quickstart](https://koszullab.github.io/metator/quickstart/)
runs metaTOR on the test data of this repository.

## Migrating from metaTOR 1.x

metaTOR 2.0 no longer aligns reads, and removed the commands that hicstuff,
metacooler, instaGRAAL or pairix replace (`qc`, `contactmap`, `metacool`,
`scaffold`, `pairs`). [Migrating to 2.0](https://koszullab.github.io/metator/migration/)
lists every change, and [Other tools](https://koszullab.github.io/metator/other-tools/)
what replaces each command.

## Contributing

Contributions are welcome: see
[CONTRIBUTING.md](https://github.com/koszullab/metaTOR/blob/master/CONTRIBUTING.md).

## Citation

> Lyam Baudry, Théo Foutel-Rodier, Agnès Thierry, Romain Koszul and Martial
> Marbouty. 2019. MetaTOR: A Computational Pipeline to Recover High-Quality
> Metagenomic Bins From Mammalian Gut Proximity-Ligation (meta3C) Libraries.
> *Frontiers in Genetics* 10: 753.
> <https://doi.org/10.3389/fgene.2019.00753>

Related work ([BibTeX and the tools metaTOR runs](https://koszullab.github.io/metator/citation/)):

* [Metagenomic chromosome conformation capture (meta3C) unveils the diversity of chromosome organization in microorganisms](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4381813/), Martial Marbouty, Axel Cournac, Jean-François Flot, Hervé Marie-Nelly, Julien Mozziconacci, and Romain Koszul. eLife, 2014.
* [Scaffolding bacterial genomes and probing host-virus interactions in gut microbiome by proximity ligation (chromosome capture) assay](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5315449/), Martial Marbouty, Lyam Baudry, Axel Cournac, and Romain Koszul. Science Advances, 2017.
* [MetaHiC phage-bacteria infection network reveals active cycling phages of the healthy human gut](https://elifesciences.org/articles/60608), Martial Marbouty, Agnès Thierry, Gaël A Millot, Romain Koszul. eLife, 2021.
* [Phages with a broad host range are common across ecosystems](https://hal.science/pasteur-05269627), Amaury Bignaud, Devon Conti, Agnès Thierry, Jacques Serizay, Karine Labadie, Julie Poulain, Olivia Cheny, Maritrini Colón-González, Laurent Debarbieux, Marianna Guerrero-Osornio, Sophie Helaine, Peter Hill, Gwenaelle Le Tinier, Gael A Millot, Lucia Morales, Andrés Parada, Nadia Riera, Gregorio Iraola, Romain Koszul & Martial Marbouty. Nature Microbiology, 2025.

## License

metaTOR is distributed under the GNU General Public License v3 (GPLv3): see
[LICENSE](https://github.com/koszullab/metaTOR/blob/master/LICENSE).

## Contact

* martial.marbouty@pasteur.fr
* jacques.serizay@pasteur.fr
* romain.koszul@pasteur.fr
* amaury.bignaud@pasteur.fr
* lyam.baudry@pasteur.fr
* thfoutel@pasteur.fr

[Spatial Regulation of Genomes](https://research.pasteur.fr/en/team/spatial-regulation-of-genomes/) (Institut Pasteur, Paris)
