# Clustering mobile genetic elements

Plasmids and phages, the mobile genetic elements (MGEs) of a community, are
clustered poorly by the partition: their contigs are small, an MGE may be in
contact with several hosts, and they carry none of the marker genes that
validation counts. `metator mge` clusters them on their own: given the contigs
of the assembly that are MGEs, found beforehand with an annotation tool, it
groups those that belong to the same element from the Hi-C contacts between
them, which the contact map of the run holds, and stores the MGE clusters in
it.

![The clustering of the MGE contigs](../images/metator_mge.png)

## How it works

1. **The contacts.** The contacts between each two MGE contigs are those of
   the contact map of the run, summed by pair of contigs, whatever its bin
   size. Within a contig, only the pairs whose mates lie more than 1 kb apart
   count: closer ones say nothing about the element. A map does not tell
   them by its pixels, but by its contact statistics: those that metacooler
   counts per contig as it makes the map from the pairs files (its QC block,
   `qc_intra_contig_long`: the pairs of reads 1001 bp apart or more).
2. **The scores.** The score of two contigs *i* and *j* is their contacts
   divided by the geometric mean of the contacts within each:
   *A<sub>i,j</sub>* = *x<sub>i,j</sub>* / √(*x<sub>i,i</sub>* × *x<sub>j,j</sub>*).
3. **The merges.** The two contigs with the highest score are merged, when
   it is above `--threshold-bin` (0.8 by default): the merged cluster takes the
   contacts of both, and its scores with the other contigs are computed again.
   This repeats until no score is above the threshold. Two equal scores are
   merged in the order of the MGE contigs: the longest first, ties by name,
   whatever the order of the map.

Each contig that was never merged is an MGE cluster of its own. The clusters
are numbered in the order they are made: the merged ones first, then the
others, in the order of the MGE contigs. So the MGE clusters, and their
names, are the same on the map of the contigs and on the map of the clusters
written from it, grouped by cluster.

A pair with a read at position 0 of its contig, which a pairs file, 1-based,
should not hold (but some do), counts within its contig as if that read were
on the first base: metacooler moves it there as it makes the map. Its mate
1001 bp away is then counted as 1000 bp away, below the threshold, where
metaTOR 1.x counted it from the pairs files.

## Running it

```bash
# The map of the run, once metator network stored its network in it (by
# network, pipeline or validation)
metator mge hicstuff/sample.metacool -a assembly.fa -m mges.txt -o mge
```

| Input | What it is |
|---|---|
| `MAP` | The contact map of the run, which `metator network` stored the network of its contigs in (that of `metator network MAP`, `pipeline MAP` or `validation MAP`, or their `--output-map`): it gives the contacts of the MGE contigs, and the size, GC content and hits of each, which `metator network` stored in it. Its clusters, if it has any, are not read. It must hold the contact statistics of the reads of its contigs, at 1001 bp: `hicstuff pipeline -M metacool` and `metacooler create pairs` count them as they make the map from the pairs files, with metacooler 0.4.0 or later (`metacooler create pairs --qc-threshold 1001` with an earlier one) |
| `--assembly` (`-a`) | The assembly that the map was made from, as a FASTA file, plain or compressed with gzip |
| `--mges` (`-m`) | The MGE contigs, one name per line. Those that the map lacks (a map made with a contig filter) are left out, with a warning that names them |

| Option | Default | Effect |
|---|---|---|
| `--threshold-bin` (`-s`) | 0.8 | The score above which two contigs, or clusters, are merged. A lower threshold merges more |
| `--random` (`-r`) | off | Shuffle the contigs between the clusters, keeping the number of contigs of each: a random clustering of the same shape, to compare the real one with. Stored as the round `mge_random`, and written to `mges_binned_random.fa` |
| `--outdir` (`-o`) | current directory | Where `mges_binned.fa` goes |
| `--output-map` | none | Write the map with the MGE clusters to this file, and leave `MAP` as it was (see [In the contact map of the run](#in-the-contact-map-of-the-run)) |

A map that `mge` cannot cluster the MGE contigs of is refused before anything
is written, with the reason and how to make one: a map without the contact
statistics of its reads, with them measured from its pixels (`metacooler qc
--store`) or at another distance than 1001 bp, or without the columns of
`metator network` (after a run with `--output-map`, give `mge` the map that
the run wrote); and a map that holds none of the MGE contigs. So is an
assembly that is not a FASTA file, plain or compressed with gzip, as
`network`, `validation` and `pipeline` refuse it.

## Outputs

| Output | What it holds |
|---|---|
| `mges_binned.fa` | One record per MGE cluster, `MetaTOR_MGE_00001` and on: its contigs, the longest first (ties by name), joined by 200 `N`, in lines of 80 bases |
| The round `mge` of the map | One cluster per MGE cluster, named as its record: its contigs, length, GC content, hits and clustering score (`BinningScore`) (see [In the contact map of the run](#in-the-contact-map-of-the-run)) |
| `metator_mge_<date>.log` | The log of the run, in the current directory |

The clustering score of an MGE cluster is the lowest score at which its
contigs were merged (missing for a single contig). The record of a cluster
holds a single sequence, so that the tools that check the quality of an MGE
genome, which expect one sequence per genome, read each cluster as one.
metacooler writes it (`metacooler.write_fasta`, stitched).

`mge` writes no summary table: the round `mge` holds what
`mges_bin_summary.tsv` held (`mges_cluster_summary.tsv` in a development
build). An earlier run's, in the output directory, would pass for this run's:
it is refused without `-F`, before any work, and `-F` removes it.

## In the contact map of the run

`mge` stores the MGE clusters in the contact map, `MAP` (see
[The contact map of the clusters](metacool.md)): they become the round `mge`
of its clustering rounds (see
[Clustering rounds](metacool.md#clustering-rounds)), after the rounds that
it holds, those of a clustering run (`core`, `overlapping`,
`recursive_<k>`, `final`), which are kept: their clusters, contigs and
values. Each MGE
cluster is named as its record of `mges_binned.fa` (`MetaTOR_MGE_00001`),
with its contigs and its columns: its length and number of contigs, which
metacooler computes from its contigs, and its GC content (`GC`), hits
(`Hit`) and clustering score (`BinningScore`). The entry `mge` of the
`metator` key of the metadata of the map records the run: the version of
metaTOR, its command line and parameters, and what it read of the map, under
`read` (the `GC_content`, `Hit` and `qc_intra_contig_long` columns of its
contigs table, and the `qc` key of its metadata). The clusters of every round
are one table, where a column that a round
lacks is missing in its rows: so a column of whole numbers, or of booleans,
of the rounds of the clustering run, which the round `mge` lacks, is stored
as decimals, or as text (`True`, `False`), from then on. The genomes and
contigs tables of the map, the final MAGs if it has any, are left as they
were, and so is the rest of the map: its contacts, contact statistics and
metadata.

```python
import metacooler

mc = metacooler.MetaCooler("hicstuff/sample.metacool")
mc.clusters[mc.clusters["round"] == "mge"]  # one row per MGE cluster
mc.cluster_members("mge")                   # their contigs, in map order
mc.metadata["metator"]["mge"]               # the run
```

`metator dump` writes them out with the other rounds of the map of the
clusters of a run, in `clusters.txt` and `cluster_members.txt` (see
[`dump`](../cli.md#dump)), and the round `mge_random` of a run with
`--random` too.

A second run replaces the round `mge`, and only it; `--random` stores the
round `mge_random`, which replaces only an earlier `mge_random`. The map is
updated as `metator partition MAP` updates it (see
[Step by step](metacool.md#step-by-step)): in an edit session of metacooler,
validated in full, then renamed over it, keeping its mode (through a symbolic
link, the file it points to is updated), and locked while the run writes it;
`--output-map FILE` writes the map to FILE instead (replaced only with `-F`),
and leaves `MAP` as it was, which a run from an entry of a `.mmetacool`
collection needs. What `partition MAP` refuses is refused before any work,
with the reason in the log: a map that is not there or that metacooler cannot
open, an entry without `--output-map`, an `--output-map` that is one of the
inputs of the run (even with `-F`) or that exists without `-F`, and a map
that another run holds.

`mge` reads no cluster of the map: it runs as well on the map of the contigs,
once `metator network` stored its network there, as on the map of the
clusters, and gives the same MGE clusters on both. But a run that replaces
the clusters of the map drops the rounds `mge` and `mge_random`, and their
entries, with a warning that says so, as it drops every round that it does
not make: `partition MAP`, `validation MAP` and `pipeline MAP`. Run `mge`
again on the map that they write to store the MGE clusters there.

## The hosts

`metator mge` does not look for the host MAGs of the MGE clusters: the
round `mge` has no column for them (`mges_bin_summary.tsv` had three,
always empty). The Python module `metator.host` (see [the API](../api/mge.md))
holds the functions, still in development, that associate MGE clusters with
their hosts from the network; no command runs them yet. They join the
network to `contig_data_final.txt` by the IDs of the contigs. The IDs of the
network of a map (`metator.network.map_network`) are the positions of its
contigs in the map as it is read, and a run that regroups the map
(`validation`, `pipeline`) moves them: the network of the map of the
clusters has other IDs than the tables of the run. Give the network the IDs
of the table first, by contig name, with `metator.network.relabel_network`.

The phages of a community and their hosts are the subject of two of the
studies listed in [Citing metaTOR](../citation.md): Marbouty et al., 2021,
and Bignaud et al., 2025.
