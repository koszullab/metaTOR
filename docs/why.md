# Why metaTOR?

Binning sorts the contigs of a metagenome assembly into genomes, the
metagenome-assembled genomes (MAGs). Coverage and sequence composition, which
most binners use, separate genomes poorly when they are equally abundant or
similar in sequence.

metaHiC (or meta3C) adds another signal: its read pairs join contigs of the
same cell, and so of the same genome, far more often than contigs of
different genomes. metaTOR clusters on that signal. It was described in
[Baudry *et al.*, *Frontiers in Genetics*, 2019](https://doi.org/10.3389/fgene.2019.00753),
and is developed by the [Koszul laboratory](https://github.com/koszullab) at
Institut Pasteur.

![The metaTOR logo](images/metator_logo.png){ width="200" }

## What it does

- **One command, from contact map to MAGs.** `metator pipeline` builds the
  network of contigs from a `.metacool` contact map (or from pairs or BAM
  files), partitions it, validates the clusters, and writes their FASTA files
  and quality.
- **Robust partitions.** Louvain (or Leiden) runs many times, and only
  contigs that consistently fall together are clustered.
- **Recursive decontamination.** [miComplete](https://pypi.org/project/micomplete/)
  checks each cluster; contaminated ones are partitioned again.
- **A contact map of the clusters.** A
  [metacooler](https://github.com/js2264/metacooler) file, a
  [cooler](https://github.com/open2c/cooler) of all the contigs grouped by
  cluster.
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

`-M metacool` needs hicstuff 4.0.0 or later. The [Quickstart](quickstart.md)
runs these steps on test data, [The end-to-end pipeline](end-to-end.md) goes
from the reads to MAGs and MGE hosts, and [Concepts](concepts.md) explains
how it works.
