# MetaTOR

[![Bioconda](https://img.shields.io/conda/vn/bioconda/metator?color=blue)](https://anaconda.org/bioconda/metator)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/metator.svg)
[![Build Status](https://github.com/koszullab/metator/actions/workflows/ci.yml/badge.svg)](https://github.com/koszullab/metator/actions)
[![codecov](https://codecov.io/gh/koszullab/metator/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/metator)
[![License: GPLv3](https://img.shields.io/badge/License-GPL%203-0298c3.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Metagenomic Tridimensional Organisation-based Reassembly: metaTOR clusters the
contigs of a metagenome assembly into metagenome-assembled genomes (MAGs) from
the Hi-C contacts between them (metaHiC).

<p align="center">
  <img src="https://raw.githubusercontent.com/koszullab/metaTOR/master/docs/images/metator_logo.png" width="200">
</p>

The documentation is at <https://koszullab.github.io/metator/>: a
[quickstart](https://koszullab.github.io/metator/quickstart/), the [concepts](https://koszullab.github.io/metator/concepts/), a guide per
step, the [CLI reference](https://koszullab.github.io/metator/cli/) and the [file formats](https://koszullab.github.io/metator/formats/).

1. [Installation](#installation)
2. [Usage](#usage)
   1. [From the reads to the clusters](#from-the-reads-to-the-clusters)
   2. [Step by step](#step-by-step)
   3. [Outputs](#outputs)
   4. [Mobile genetic elements](#mobile-genetic-elements)
3. [Migrating from metaTOR 1.x](#migrating-from-metator-1x)
4. [Contributions](#contributions)
5. [References](#references)
6. [Contact](#contact)

## Installation

### With conda

`conda` is the recommended way to install the latest `metator` release:

```sh
conda create -n metator bioconda::metator
```

### With pip

`metator` is on PyPI, with everything it runs: miComplete finds the marker
genes of the clusters with Prodigal and HMMER through their Python bindings,
pyrodigal and pyhmmer. metaTOR needs Python 3.10 to 3.14. Install it in an
environment of its own, as pip may change the versions of packages that
other tools rely on.

```sh
pip3 install metator
```

### Development version

The development version is managed with [uv](https://docs.astral.sh/uv/):

```sh
git clone https://github.com/koszullab/metaTOR
cd metaTOR
uv sync
uv run metator --help
```

See [CONTRIBUTING.md](https://github.com/koszullab/metaTOR/blob/master/CONTRIBUTING.md) to run the tests and checks.

## Usage

metaTOR starts from a `.metacool` contact map of the contigs: a
[cooler](https://github.com/open2c/cooler) file whose chromosomes are the
contigs, which [hicstuff](https://github.com/koszullab/hicstuff) makes from the
reads (`hicstuff pipeline -M metacool`) and
[metacooler](https://github.com/js2264/metacooler) from pairs files
(`metacooler create pairs`). Every command takes the map as its first argument:
the clustering steps and `mge` store what they find in it, so that the map is
what the steps share, and it becomes the contact map of the clusters; `dump`
only reads it.

```sh
metator {network|partition|validation|pipeline|mge|dump} MAP [OPTIONS]
```

| Command | What it does |
|---|---|
| `network` | Builds the network of the contigs from the map, normalizes it, and stores it in the map. |
| `partition` | Partitions the network many times with Louvain or Leiden, and stores the core and overlapping clusters in the map. |
| `validation` | Checks the overlapping clusters with miComplete, partitions the contaminated ones again, writes the FASTA files of the final clusters, and makes the map the contact map of the clusters. |
| `pipeline` | Runs `network`, `partition` and `validation` in one command. |
| `mge` | Clusters the contigs of mobile genetic elements (MGEs) from their contacts, and stores the MGE clusters in the map. |
| `dump` | Writes the tables of a run (contig data, cluster summary, binning, clustering rounds, network) from the map of its clusters. |

`metator --help` lists the commands, `metator <command> --help` their options,
and `metator --version` prints the version. A run writes its map to its output
directory (`-o`), under the name of `MAP`, which it leaves as it was:
metacooler copies `MAP` there, the run edits the copy, under its lock, and
renames it into place once it is validated. `--inplace` updates `MAP` itself
instead, and so does a run whose output directory holds `MAP` already (the
steps of a run, one by one, in one output directory). [The contact map of the clusters](https://koszullab.github.io/metator/guides/metacool/)
says what the map holds, how it is written, and what a run refuses.

### From the reads to the clusters

```sh
# 1. The map of each library: one bin per restriction fragment, every contig
#    in the genome __unbinned__.
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff_L1 -P L1 L1_R1.fq.gz L1_R2.fq.gz
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff_L2 -P L2 L2_R1.fq.gz L2_R2.fq.gz

# 2. One map of both libraries, their contacts summed.
metacooler merge hicstuff_L1/L1.metacool hicstuff_L2/L2.metacool -o sample.metacool

# 3. The clusters: binning/sample.metacool, the contact map of the clusters.
#    sample.metacool is left as it was.
metator pipeline sample.metacool -a assembly.fa -e DpnII -d depth.txt -o binning

# 4. The tables of the run, from the map.
metator dump binning/sample.metacool -o tables
```

`hicstuff pipeline -M metacool` needs hicstuff 4.0.0 or later; with an
earlier hicstuff, make the map from its pairs files,
`<prefix>.valid_idx.pairs.gz`, with `metacooler create pairs`.
metacooler makes the map from the pairs files of any aligner:

```sh
metacooler create pairs --pairs L1.pairs.gz --pairs L2.pairs.gz --bins 1000 --no-ps \
    --assembly-name assembly -o sample.metacool
```

The map can have any bin size, from one bin per contig to fixed bins.
[Other tools](https://koszullab.github.io/metator/other-tools/#from-the-reads-to-the-map) gives the details.

### Step by step

The three steps of `pipeline` can run one by one, in one output directory,
which gives the map that `pipeline` writes: the first writes the map there,
and the next ones update it in place.

```sh
# 1. The network, stored in binning/sample.metacool.
metator network sample.metacool -a assembly.fa -e DpnII -d depth.txt -o binning
# 2. The partition: the map stores its core and overlapping clusters.
metator partition binning/sample.metacool -o binning
# 3. The validation: the final clusters. The map becomes the map of the clusters.
metator validation binning/sample.metacool -a assembly.fa -o binning
```


### Outputs

**In the output directory** (`-o`): the map of the run, named as `MAP`;
`final_clusters/`, one FASTA file per final cluster; and the log of the run,
`metator_<command>_<date>.log`, which starts with the version, the command line
and the value of each option, and ends with what the run wrote and how long it
took.

**In the map**, which `validation` and `pipeline` make the contact map of the
clusters:

* one genome per final cluster, and `__unbinned__` for the contigs in no
  cluster: the contigs grouped genome by genome, the longest first;
* the genomes table holds the summary of the clusters (size, GC content,
  miComplete's completeness and redundancy, N50, Hi-C and shotgun abundances);
* the contigs table holds what the network measured on each contig
  (`GC_content`, `Hit`, `Shotgun_coverage`, `Restriction_site`), and the
  contact statistics of its reads (`qc_intra_contig_short`,
  `qc_intra_contig_long`);
* the bins table holds the normalization of the network, as the weights
  `metator_<normalization>`, the way cooler stores balancing weights;
* the `/clusters` group holds every clustering round of the run (`core`,
  `overlapping`, `recursive_<k>`, `final`, and `mge` after `metator mge`);
* the metadata hold, under the `metator` key, one entry per step: its version,
  command line, parameters, and what it read of the map.

cooler and the tools built on it read the map as any `.cool` file, and
metacooler adds the genomes:

```python
import cooler
import metacooler

clr = cooler.Cooler("binning/sample.metacool")
clr.matrix(balance="metator_empirical_hit")  # the normalized contacts
mc = metacooler.MetaCooler("binning/sample.metacool")
mc.genomes  # one row per cluster: the summary of the clusters
mc.contigs  # one row per contig: its GC content, hits... and its cluster
```

**The tables of a run**, which `metator dump MAP -o DIR` writes from the map
of the clusters: `contig_data_final.txt`, `cluster_summary.txt`,
`binning.txt` (which anvi'o imports), `clusters.txt`, `cluster_members.txt`
and `network.txt`. A run keeps the first five in its temporary directory,
which `--no-clean-up` (`-N`) keeps, with miComplete's results and the FASTA
files of the clusters it checked; `network.txt` comes from the map only. [File formats](https://koszullab.github.io/metator/formats/) describes each one.

### Mobile genetic elements

![mge](https://raw.githubusercontent.com/koszullab/metaTOR/master/docs/images/metator_mge.png)

`metator mge MAP -a assembly.fa -m mges.txt` clusters the contigs that the list
`mges.txt` names as MGEs, from the contacts of the map between them, and those
within each, between reads more than 1 kb apart, which the map counts when
hicstuff or metacooler made it from the pairs (metacooler 0.4.0 or later). It
reads the map after `metator network`, `pipeline` or `validation`, writes
`mges_binned.fa`, one record per MGE cluster (`MetaTOR_MGE_00001`...), its
contigs joined by 200 `N`, and stores the MGE clusters in the map of its output
directory, as the round `mge` of its clustering rounds (`-o binning` updates
`binning/sample.metacool` in place). See
[Clustering mobile genetic elements](https://koszullab.github.io/metator/guides/mge/).

## Migrating from metaTOR 1.x

metaTOR 2.0 reads no reads, alignments or pairs files: hicstuff and metacooler
make the map. It removed the commands that did what hicstuff, metacooler,
instaGRAAL or pairix do (`qc`, `contactmap`, `metacool`, `scaffold`, `pairs`)
and the plots: [Other tools](https://koszullab.github.io/metator/other-tools/) says what replaces each one.
Its outputs are named after clusters, where 1.x named them after bins
(`cluster_summary.txt` was `bin_summary.txt`). [Migrating to 2.0](https://koszullab.github.io/metator/migration/)
lists every change to the command line and the outputs, and the
[changelog](https://koszullab.github.io/metator/changelog/) has the details.

## Contributions

We welcome contributions to the project. Please refer to the [CONTRIBUTING.md](https://github.com/koszullab/metaTOR/blob/master/CONTRIBUTING.md) file for more information.

## References

* [Metagenomic chromosome conformation capture (meta3C) unveils the diversity of chromosome organization in microorganisms](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4381813/), Martial Marbouty, Axel Cournac, Jean-François Flot, Hervé Marie-Nelly, Julien Mozziconacci, and Romain Koszul. eLife, 2014.
* [Scaffolding bacterial genomes and probing host-virus interactions in gut microbiome by proximity ligation (chromosome capture) assay](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5315449/), Martial Marbouty, Lyam Baudry, Axel Cournac, and Romain Koszul. Science Advances, 2017.
* [MetaHiC phage-bacteria infection network reveals active cycling phages of the healthy human gut](https://elifesciences.org/articles/60608), Martial Marbouty, Agnès Thierry, Gaël A Millot, Romain Koszul. eLife, 2021.
* [Phages with a broad host range are common across ecosystems](https://hal.science/pasteur-05269627), Amaury Bignaud, Devon Conti, Agnès Thierry, Jacques Serizay, Karine Labadie, Julie Poulain, Olivia Cheny, Maritrini Colón-González, Laurent Debarbieux, Marianna Guerrero-Osornio, Sophie Helaine, Peter Hill, Gwenaelle Le Tinier, Gael A Millot, Lucia Morales, Andrés Parada, Nadia Riera, Gregorio Iraola, Romain Koszul & Martial Marbouty.
Nature Microbiology, 2025

## Contact

### Authors

* martial.marbouty@pasteur.fr
* jacques.serizay@pasteur.fr
* romain.koszul@pasteur.fr
* amaury.bignaud@pasteur.fr
* lyam.baudry@pasteur.fr
* thfoutel@pasteur.fr

### Research lab

[Spatial Regulation of Genomes](https://research.pasteur.fr/en/team/spatial-regulation-of-genomes/) (Institut Pasteur, Paris)
