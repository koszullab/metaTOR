# MetaTOR

[![PyPI version](https://badge.fury.io/py/metator.svg)](https://badge.fury.io/py/metator)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/metator.svg)
[![Build Status](https://github.com/koszullab/metator/actions/workflows/python-package.yml/badge.svg)](https://github.com/koszullab/metaTOR/actions)
[![codecov](https://codecov.io/gh/koszullab/metator/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/metator)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/koszullab/metator)](https://hub.docker.com/r/koszullab/metator)
[![Read the docs](https://readthedocs.org/projects/metator/badge)](https://metator.readthedocs.io)
[![License: GPLv3](https://img.shields.io/badge/License-GPL%203-0298c3.svg)](https://opensource.org/licenses/bo-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Metagenomic Tridimensional Organisation-based Reassembly - A set of scripts that streamlines the processing and binning of metagenomic metaHiC datasets.

## Table of contents

* [Installation](#Installation)
  * [Requirements](#Requirements)
  * [Using pip](#Using-pip)
  * [Using docker container](#Using-docker-container)
* [Usage](#Usage)
* [Output files](#Output-files)
* [References](#References)
* [Contact](#Contact)

## Installation

### Requirements

* Python 3.6 or later is required.
* The following librairies are required but will be automatically installed with the pip installation: `numpy`, `scipy`, `sklearn`, `pandas`, `docopt`, `networkx` `biopython` `pyfastx` and `pysam`.
* The following software should be installed separately if you used the pip installation:
  * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [samtools](http://www.htslib.org/)
  * [louvain](https://sourceforge.net/projects/louvain/) (original
        implementation).
  * [networkanalysis](https://github.com/vtraag/networkanalysis) (not
    necessary only if you want to use Leiden algorithm to partition the network)
  * [checkm](https://github.com/Ecogenomics/CheckM)
  * [pairix](https://github.com/4dn-dcic/pairix) (not necessary if you do not use metator contact map)

### Using pip

```sh
   pip3 install metator
```

or, to use the latest version:

```sh
   pip3 install -e git+https://github.com/koszullab/metator.git@master#egg=metator
```

In order to use Louvain or Leiden it's necessary to set a global variable `LOUVAIN_PATH` and `LEIDEN_PATH` depending on which algorithm you wan to use with the absolute path where the executable are.

For Louvain algorithm in the directory where you have the archive file (available in the external directory of this repository):

```sh
YOUR_DIRECTORY=$(pwd)
tar -xvzf louvain-generic.tar.gz
cd gen-louvain
make
export LOUVAIN_PATH=$YOUR_DIRECTORY/gen-louvain/
```

For Leiden algorithm, clone the networkanalysis repository from github and build the Java script. Then you can export the Leiden path:

```sh
export LEIDEN_PATH=/networkanalysis_repository_path/build/libs/networkanalysis-1.2.0.jar
```

### Using docker container

A dockerfile is also available if that is of interest. You may fetch the image by running the following:

```sh
    docker pull koszullab/metator
```

## Usage

```sh
    metator {network|partition|validation|pipeline} [parameters]
```

A metaTOR command takes the form `metator action --param1 arg1 --param2
arg2 #etc.`

There are three actions/steps in the metaTOR pipeline, which must be run in the
following order:

* `network` : Generate metaHiC contigs network from fastq reads or bam files and normalize it.
* `partition` : Perform the Louvain or Leiden community detection algorithm many times to bin contigs according to the metaHiC signal between contigs.

* `validation` : Use CheckM to validate the bins, then do a recursive decontamination step to remove contamination.

There are a number of other, optional, miscellaneous actions:

* `pipeline` : Run all three of the above actions sequentially or only some of them depending on the arguments given. This can take a while.
* `contactmap` : Generates a contact map from one bin from the final ouptut of metaTOR.

* `version` : display current version number.

* `help` : display help message.

A tutorial is available [here](example/metator_tutorial.md) to explain how to use metaTOR. More advanced tutorials to analyze the output files are also available:

* [Anvio](https://merenlab.org/software/anvio/) manual curation of the contaminated bins. Available [here](example/manual_curation_of_metator_MAGs.md).
* Visualization and scaffolding of the MAGs with the contactmap modules of MetaTOR. Available [here](example/MAG_visualization_and_scaffolding.md).

Principle of MetaTOR pipeline:

![metator_pipeline](example/images/metator_figure.png)

## Output files

The output files will be in the ouput directory given as parmater or in working directory if no paramater were given. Depending on the command used, different files will be in the ouptut:

| Files/Commands | description | network | partition | validation | pipeline |
| - | :-: | :-: | :-: | :-: | :-: |
| alignment_N_for.bam |Bam file of the forward alignment |X|||X|
| alignment_N_rev.bam |Bam file of the reverse alignment|X|||X|
| alignment_N.pairs |Pairs of the merge alignment|X|||X|
| network.txt |Normalized network of the metaHiC library|X|||X|
| contig_data_network.txt |Information on contigs after network step|X||||
| clustering_matrix_partition.txt |Matrix of clustering from the partition iterations||X|||
| contig_data_partition.txt |Information on contigs after partition step||X|||
| overlapping_checkm_results.txt |CheckM results summary file from the partition step|||X|X|
| overlapping_checkm_taxonomy.txt |CheckM taxonomy file from the partition step|||X|X|
| recursive_checkm_results.txt |CheckM results summary file from the recursive step|||X|X|
| recursive_checkm_taxonomy.txt |CheckM taxonomy file from the recursive step|||X|X|
| clustering_matrix_validation.txt |Matrix of clustering from the recursive iterations|||X||
| clustering_matrix.txt |Matrix of clustering from the partition and recursive iterations||||X|
| contig_data_final.txt |Information on contigs after whole pipeline|||X|X|
| **bin_summary.txt** |**Information on the final bins**|||**X**|**X**|
| binning.txt |File with contigs names and their final clustering|||X|X|
| overlapping_bin |Directory with the fasta of the partition bins||X|X|X|
| recursive_bin |Directory with the fasta of the recursive bins|||X|X|
| **final_bin** |**Directory with the fasta of the final bins**|||**X**|**X**|

**Bam alignment files**
For the bam alignments files, only the aligned reads are kept and the bam are sorted by name. The N value correspond to the id (order of the given fastq started at 0)

**Pairs aligment files**
This format is used to store the relevant information of mapping of the merged alignment. It's a s The N value correspond to the id (order of the given fastq started at 0). It is a tab-separated format holding informations about Hi-C pairs. It has an [official specification](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) defined by the 4D Nucleome data coordination and integration center. Here we kept 7 columns readID-chr1-pos1-chr2-pos2-strand1-strand2.

**Network file**
This is a tsv file of the network with edgelist form: Id of the first contig, id of the second contig and the weigth of edge normalized.

**Contig data files**
These are the files with all the informations from the contigs:
|ID|Name|Size|GC_content|Hit|Shotgun_coverage|Restriction_site|Core_bin_ID|Core_bin_contigs|Core_bin_size|Overlapping_bin_ID|Overlapping_bin_contigs|Overlapping_bin_size|Recursive_bin_ID|Recursive_bin_contigs|Recursive_bin_size|Final_bin|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|1|NODE_1|642311|38.6876450815882|3837|41.1565|2006|1|65|2175226|1|396|6322353|1|52|2158803|MetaTOR_1_1|
|2|NODE_2|576356|30.235826468363303|1724|24.509|1256|2|40|1735419|2|401|735419|0|-|-|MetaTOR_2_0|
|3|NODE_3|540571|42.305266098255366|2188|14.5855|3405|3|127|6409484|3|431|13615480|1|112|6385126|MetaTOR_3_1|

They have to have the header when they are use as input but the order of the columns nd if they are others columns doesn't matter when they are used as input files.
Depending on which step of the pipeline have been launch they just have some of these columns:

* `contig_data_network.txt`: columns: ID, Name, Size, GC content, Hit, Shotgun_coverage and Restriction Site only.
* `contig_data_partition.txt`: The same as the previous with the information of core bins and overlapping bins.
* `contig_data_final.txt`: All the columns.

The shotgun coverage will be filled only if the depth.txt file is given, otherwise it will be filled with `-`. This column is only necessarry for the *abundance* and the *theoritical_hit* normalization. The restriction will also be filled with `-` if no enzyme are given. This column is only necessary for the *RS* and the *theoritical_hit* normalization.
Moreover, if the contig is not binned (no HiC reads mapped on it) all the columns with binning information will be filled with `-`, and if a bin is not recursively decontamined because it's smaller than the size threshold or it doesn't have any contamination the recusive bin information will be filled `0`, `-`, `-`. Finally, if the bin is not in a final bin, it will be annotated `ND` in the last column (for not determined).

**clustering matrix files**
The clustering matrix files are at the `.npz` format which is a compresed foramt for sparsed matrix. This sparsed matrix contains the ratio of time each pairs of contigs are clusterize together by the algorithm of clustering (either Louvain or Leiden). The partition matrix contains the information for the partition step, the recursive one for the recursive step and the general one is the mean of both matrices. Be careful the index of the contigs are **zero-based** and not one-based as in the contig data file.

It's possible to read them in python using the `scipy.sparse.load_npz()` function. If the users wants a tsv file instead, he or she could load the matrix in python using load_npz and make sure to transform the matrix in the `scipy.sparse.coo_matrix` function and used the function from metator `metator.io.save_sparse_matrix` to save it as tsv file.

**CheckM results**
Files from [checkM](https://github.com/Ecogenomics/CheckM) output. Two types of files one with the main results of checkM `checkm_results.txt` and one with the taxonomy  `checkm_taxonomy.txt` for both the partition and the recurisve bins.

**binning.txt file**
This is a tsv file with two columns: the contig name and the final were the contig is. It only contains contigs which are binned. It could be use a an input to import a binning results in anvio.

**Bin summary file**
This is the summary of the data of the final bins build with all the step of metaTOR. The HiC coverage is the number of contacts (intra and inter contigs) per kilobase in the whole bin. The Shotgun coverage is the mean coverage normalized by the size of the shotgun reads from the depth file.

||lineage|completness|contamination|size|contigs|N50|longest_contig|GC|coding_density|taxonomy|Coverage|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaTOR_8_1|o__Clostridiales|68.29|2.46|1431612|15|116129|291620|26.36|87.97|k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales|146.46719755483332|
|MetaTOR_8_2|o__Clostridiales|58.42|2.01|1396934|58|41290|174682|28.89|83.70|k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales|22.252416224710686|
|MetaTOR_8_3|o__Clostridiales|49.37|0.94|1420821|82|33095|89964|30.29|83.24|k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae_3;g__Clostridium_3|44.27369196532141|

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
