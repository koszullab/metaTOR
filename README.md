# metaTOR

[![PyPI version](https://badge.fury.io/py/metator.svg)](https://badge.fury.io/py/metator)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/metator.svg)
[![Build Status](https://github.com/koszullab/metator/actions/workflows/python-package.yml/badge.svg)](https://github.com/koszullab/metaTOR/actions)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/koszullab/metator)](https://hub.docker.com/r/koszullab/metator)
[![Read the docs](https://readthedocs.org/projects/metator/badge)](https://metator.readthedocs.io)
[![License: GPLv3](https://img.shields.io/badge/License-GPL%203-0298c3.svg)](https://opensource.org/licenses/bo-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Metagenomic Tridimensional Organisation-based Reassembly - A set of scripts that 
streamlines the processing and binning of metagenomic metaHiC datasets.

## Installation

### Requirements:

* Python 3.6 or later is required.
* The following librairies are required but will be automatically installed with
 the pip installation: ```numpy```, ```scipy```, ```sklearn```, ```pandas```, 
 ```docopt```, ```networkx``` ```biopython``` ```pyfastx``` and ```pysam```.
* The following software should be installed separetely if you used the pip 
installation:
    * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    * [samtools](http://www.htslib.org/)
    * [louvain](https://sourceforge.net/projects/louvain/) (original
        implementation).
    * [networkanalysis](https://github.com/vtraag/networkanalysis) (not 
    necessary only if you want to use Leiden algorithm to partition the network)
    * [checkm](https://github.com/Ecogenomics/CheckM)

### Using pip:

```sh
   pip3 install metator
```

or, to use the latest version:

```sh
   pip3 install -e git+https://github.com/koszullab/metator.git@master#egg=metator
```

In order to use Louvain or Leiden it's necessary to set a global variable 
```LOUVAIN_PATH``` and ```LEIDEN_PATH``` depending on which algorithm you wan to 
use with the absolute path where the executable are.

For Louvain algorithm in the directory where you have the archive file 
(available in the external directory of this repository):

```sh
YOUR_DIRECTORY=$(pwd)
tar -xvzf louvain-generic.tar.gz
cd gen-louvain
make
export LOUVAIN_PATH=$YOUR_DIRECTORY/gen-louvain/
```

For Leiden algorithm, clone the networkanalysis repository from github and build
the Java script. Then you can export the Leiden path:

```sh
export LEIDEN_PATH=/networkanalysis_repository_path/build/libs/networkanalysis-1.2.0.jar
```
### Using docker container:

A dockerfile is also available if that is of interest. You may fetch the image by running the following:

```sh
    docker pull koszullab/metator
```

## Usage

```sh
    metator {network|partition|validation|pipeline} [parameters]
```

A metaTOR command takes the form ```metator action --param1 arg1 --param2
arg2 #etc.```

There are three actions/steps in the metaTOR pipeline, which must be run in the 
following order:

* ```network``` : Generate metaHiC contigs network from fastq reads or bam files
 and normalize it.
* ```partition``` : Perform the Louvain or Leiden community detection algorithm 
many times to bin contigs together according to the metaHiC signal between 
contigs.

* ```validation``` : Use CheckM to validate the bins, then do a recursive decontamination step to remove contamination.

After the last step is completed there should be a set of bins and a table with
various descriptors of the bins.

There are a number of other, optional, miscellaneous actions:

* ```pipeline``` : Run all three of the above actions sequentially or only some 
of them depending on the arguments given. This can take a while.

* ```version``` : display current version number.

* ```help``` : display help message.

## References

* [Metagenomic chromosome conformation capture (meta3C) unveils the diversity of chromosome organization in microorganisms](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4381813/), Martial Marbouty, Axel Cournac, Jean-François Flot, Hervé Marie-Nelly, Julien Mozziconacci, and Romain Koszul, eLife, 2014
* [Meta3C analysis of a mouse gut microbiome](https://www.biorxiv.org/content/early/2015/12/17/034793), Martial Marbouty, Lyam Baudry, Axel Cournac, Romain Koszul, 2015
* [Scaffolding bacterial genomes and probing host-virus interactions in gut microbiome by proximity ligation (chromosome capture) assay](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5315449/), Martial Marbouty, Lyam Baudry, Axel Cournac, and Romain Koszul, Science Advances, 2017

## Contact

### Authors

* amaury.bignaud@pasteur.fr
* lyam.baudry@pasteur.fr
* thfoutel@pasteur.fr
* martial.marbouty@pasteur.fr
* romain.koszul@pasteur.fr

### Research lab

[Spatial Regulation of Genomes](https://research.pasteur.fr/en/team/spatial-regulation-of-genomes/) (Institut Pasteur, Paris)
