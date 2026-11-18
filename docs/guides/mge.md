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
