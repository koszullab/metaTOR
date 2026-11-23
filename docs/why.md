# Why metaTOR?

A metagenome assembly is a pile of contigs from many genomes. Binning sorts
them back into genomes: the metagenome-assembled genomes, or MAGs. Coverage
and composition, which most binners use, separate genomes poorly when they
are as abundant as each other or close in sequence.

A metaHiC (or meta3C) library adds another signal. Cross-linking captures
DNA molecules that lie close in space, and the DNA of one cell, one genome,
lies close together: the read pairs of the library join contigs of the same
genome far more often than contigs of different genomes. metaTOR clusters on
that signal. It was described in [Baudry *et al.*, *Frontiers in Genetics*,
2019](https://doi.org/10.3389/fgene.2019.00753), and is developed by the
[Koszul laboratory](https://github.com/koszullab) at Institut Pasteur.

![The metaTOR logo](images/metator_logo.png){ width="200" }

## What it does

- **One command, from a contact map to clusters.** `metator pipeline` builds
  the network of contigs from the `.metacool` contact map of a metaHiC
  library, which hicstuff makes from the reads, partitions it, validates the
  clusters, and writes one FASTA file per cluster with a table of their
  quality.
- **Robust partitions.** The network is partitioned many times with the
  Louvain (or Leiden) algorithm: the contigs that always fall together form
  core clusters, and the core clusters that fall together in most runs form a
  cluster, so that a single run of a randomised algorithm decides nothing.
- **Validation and recursive decontamination.** Each cluster is checked with
  [miComplete](https://pypi.org/project/micomplete/): the clusters it finds
  contaminated are partitioned again, on their own, and replaced by their
  parts when the most complete part keeps enough of their completeness.
- **The contact map of the clusters.** A run ends with a
  [metacooler](https://github.com/js2264/metacooler) file, a
  [cooler](https://github.com/open2c/cooler) of all the contigs grouped by
  cluster, with the clusters' quality in its genomes table, which cooler and
  the tools built on it read.
- **What comes next.** Clustering of mobile genetic elements and their
  association with their hosts. metacooler, hicstuff and instaGRAAL take the
  map and the clusters further: contact maps of single clusters, scaffolds,
  and a check of the library (see [Other tools](other-tools.md)).

## A quick look

```bash
# The contact map of the contigs, from the reads, with hicstuff
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz

# The map and the assembly to clusters
metator pipeline hicstuff/sample.metacool --enzyme DpnII --depth depth.txt --threads 16 \
    --assembly assembly.fa --outdir binning

# The clusters, their completeness and contamination
head binning/cluster_summary.txt

# The contact map of the clusters, with metacooler
metacooler dump -t genomes -H hicstuff/sample.metacool
```

`-M metacool` comes with hicstuff's next release: see
[Other tools](other-tools.md#from-the-reads-to-the-map). The
[Quickstart](quickstart.md) runs the same steps on the test data of the
repository, and [Concepts](concepts.md) explains the network, the partition
and the validation.
