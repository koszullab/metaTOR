# Why metaTOR?

Binning sorts the contigs of a metagenome assembly back into genomes, the
metagenome-assembled genomes (MAGs). Most binners use coverage and sequence
composition, which separate genomes poorly when they are equally abundant or
close in sequence.

A metaHiC (or meta3C) library adds another signal: cross-linking captures DNA
that lies close in space, and the DNA of one cell lies close together. Its
read pairs join contigs of the same genome far more often than contigs of
different genomes. metaTOR clusters on that signal. It was described in
[Baudry *et al.*, *Frontiers in Genetics*, 2019](https://doi.org/10.3389/fgene.2019.00753),
and is developed by the [Koszul laboratory](https://github.com/koszullab) at
Institut Pasteur.

![The metaTOR logo](images/metator_logo.png){ width="200" }

## What it does

- **One command, from contact map to MAGs.** `metator pipeline` builds the
  network of contigs from a `.metacool` contact map (or from pairs or BAM
  files), partitions it, validates the clusters, and writes their FASTA files
  and a quality table.
- **Robust partitions.** The network is partitioned many times with Louvain
  (or Leiden), and only contigs that fall together consistently are
  clustered, so no single random run decides.
- **Recursive decontamination.** [miComplete](https://pypi.org/project/micomplete/)
  checks each cluster; contaminated ones are partitioned again on their own.
- **The contact map of the clusters.** A run ends with a
  [metacooler](https://github.com/js2264/metacooler) file, a
  [cooler](https://github.com/open2c/cooler) of all the contigs grouped by
  cluster, which cooler-based tools read.
- **Mobile genetic elements.** metaTOR clusters MGE contigs and finds their
  host MAGs. hicstuff, metacooler and instaGRAAL take the results further
  (see [Other tools](other-tools.md)).

## A quick look

```bash
# The contact map of the contigs, from the reads, with hicstuff
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz

# The map and the assembly to clusters: binning/sample.metacool
metator pipeline hicstuff/sample.metacool --enzyme DpnII --depth depth.txt --threads 16 \
    --assembly assembly.fa --outdir binning

# The clusters, with their completeness and contamination
metacooler dump -t genomes -H binning/sample.metacool
```

`-M metacool` needs hicstuff 4.0.0 or later (see
[Other tools](other-tools.md#from-the-reads-to-the-map)). The
[Quickstart](quickstart.md) runs these steps on test data,
[The end-to-end pipeline](end-to-end.md) goes from the reads to MAGs and MGE
hosts, and [Concepts](concepts.md) explains how it works.
