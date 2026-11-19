# Binning mobile genetic elements

Plasmids and phages, the mobile genetic elements (MGEs) of a community, are
binned poorly by the partition: their contigs are small, an MGE may be in
contact with several hosts, and they carry none of the marker genes that
validation counts. `metator mge` bins them on their own: given the contigs of
the assembly that are MGEs, found beforehand with an annotation tool, it
groups those that belong to the same element from the Hi-C contacts between
them.

![The binning of the MGE contigs](../images/metator_mge.png)

## How it works

1. **The contacts.** The contacts between each two MGE contigs are counted
   from the pairs files. Within a contig, only the pairs whose mates lie more
   than 1 kb apart are counted: closer ones say nothing about the element.
2. **The scores.** The score of two contigs *i* and *j* is their contacts
   divided by the geometric mean of the contacts within each:
   *A<sub>i,j</sub>* = *x<sub>i,j</sub>* / √(*x<sub>i,i</sub>* × *x<sub>j,j</sub>*).
3. **The merges.** The two contigs with the highest score are merged, when
   it is above `--threshold-bin` (0.8 by default): the merged bin takes the
   contacts of both, and its scores with the other contigs are computed again.
   This repeats until no score is above the threshold.

Each contig that was never merged is an MGE bin of its own.

## Running it

```bash
metator mge -a assembly.fa -c binning/contig_data_final.txt -b binning/binning.txt \
    -m mges.txt -o mge binning/alignment_0_sorted.pairs.gz
```

| Input | What it is |
|---|---|
| The pairs files | The pairs files of the run, sorted and indexed with pairix as `metator network` writes them (`metator pipeline` deletes them, unless `--no-clean-up`). An uncompressed pairs file without index is read whole, which is slower: `metator pairs` indexes it |
| `--assembly` (`-a`) | The assembly, as a FASTA file |
| `--contigs` (`-c`) | The contig data table of the run: the size, GC content and hits of each contig |
| `--binning` (`-b`) | The bins of the contigs, two columns: the contig, and its bin. `binning.txt`, or the collection exported by anvi'o, whose split names (`NODE_1_split_00001`) are read as their contig (`NODE_1`) |
| `--mges` (`-m`) | The MGE contigs, one name per line |

| Option | Default | Effect |
|---|---|---|
| `--threshold-bin` (`-s`) | 0.8 | The score above which two contigs, or bins, are merged. A lower threshold merges more |
| `--random` (`-r`) | off | Shuffle the contigs between the bins, keeping the number of contigs of each: a random binning of the same shape, to compare the real one with |
| `--outdir` (`-o`) | current directory | Where the two outputs go |
| `--map` | none | Also store the MGE clusters in this contact map of the run (see [In the contact map of the run](#in-the-contact-map-of-the-run)) |

## Outputs

| File | What it holds |
|---|---|
| `mges_binned.fa` | One record per MGE bin, `MetaTOR_MGE_00001` and on: its contigs, joined by 200 `N` |
| `mges_bin_summary.tsv` | One row per MGE bin: its length, GC content, hits, binning score and contigs (see [File formats](../formats.md#mge-bin-summary)) |
| `metator_mge_<date>.log` | The log of the run, in the current directory |

The binning score of an MGE bin is the lowest score at which its contigs were
merged (`NA` for a single contig). The record of a bin holds a single
sequence, so that the tools that check the quality of an MGE genome, which
expect one sequence per genome, read each bin as one.

## In the contact map of the run

`--map MAP` also stores the MGE clusters in the `.metacool` contact map of
the run, that of `metator pipeline` (see
[The contact map of the bins](metacool.md)):

```bash
metator mge -a assembly.fa -c binning/contig_data_final.txt -b binning/binning.txt \
    -m mges.txt -o mge --map binning/metator.metacool binning/alignment_0_sorted.pairs.gz
```

They become the round `mge` of the clustering rounds of the map (see
[Clustering rounds](metacool.md#clustering-rounds)), after the rounds of the
binning run (`overlapping`, `recursive_<k>`, `final`), which are kept: their
clusters, contigs and values. Each row of `mges_bin_summary.tsv` is a
cluster, named as its record (`MetaTOR_MGE_00001`), with its contigs and the
columns of the summary: its length and number of contigs, which metacooler
computes from its contigs, and the GC content, hits, binning score and host
columns, as written (`NA` is a missing value: no number, or an empty text).
The round records the version of metaTOR and the command line and
parameters of the run. The clusters of every round are one table, where a
column that a round lacks is missing in its rows: so a column of whole
numbers, or of booleans, of the rounds of the binning run, which the round
`mge` lacks, is stored as decimals, or as text (`True`, `False`), from then
on. The genomes and contigs tables of the map, the final MAGs, are left as
they were, and so is the rest of the map: its contacts, contact statistics
and metadata.

```python
import metacooler

mc = metacooler.MetaCooler("binning/metator.metacool")
mc.clusters[mc.clusters["round"] == "mge"]  # one row per MGE cluster
mc.cluster_members("mge")                   # their contigs
```

A second run replaces the round `mge`, and only it. The map is updated as
`metator partition MAP` updates it (see
[Step by step](metacool.md#step-by-step)): written next to it, validated in
full, then renamed over it, keeping its mode (a symbolic link is replaced,
not its target), and locked while the run writes it; `--output-map FILE`
writes the map to FILE instead (replaced only with `-F`), and leaves `MAP` as
it was, which a run from an entry of a `.mmetacool` collection needs. What
`partition MAP` refuses is refused before any work, with the reason in the
log: a map that is not there or that metacooler cannot open, a map of other
contigs or lengths than the contig data (`-c`), an entry without
`--output-map`, an `--output-map` that is one of the inputs of the run (even
with `-F`) or that exists without `-F`, a map that another run holds, and
`--output-map` without `--map`. So are `--random`, a random binning to
compare the real one with, which is not stored in a map, and a map given in
place of the pairs files: the contacts within an MGE contig are counted
between reads more than 1 kb apart, which a map does not tell, so `mge`
reads the pairs files, with or without `--map`.

The round `mge` belongs to the binning that it was made with (`-b`, the
final clusters of the run). A run that replaces the clusters of the map
drops it, with a warning that says so: `partition MAP`, `validation MAP`,
`pipeline MAP` and `metator metacool MAP`. Run `mge --map` again with the
new clusters to store their MGE clusters.

## The hosts

The last three columns of the summary, `MetagenomicBin`, `AssociationScore`
and `AssociatedBins`, are for the host MAGs of each MGE bin: `metator mge`
leaves them empty (`None`, `NA`, `NA`). The Python module `metator.host`
(see [the API](../api/mge.md)) holds the functions, still in development, that
associate MGE bins with their hosts from the network; no command runs them
yet.

The phages of a community and their hosts are the subject of two of the
studies listed in [Citing metaTOR](../citation.md): Marbouty et al., 2021,
and Bignaud et al., 2025.
