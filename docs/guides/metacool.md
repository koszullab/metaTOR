# The contact map of the clusters

A run starts from a `.metacool` contact map of the contigs and ends with the
contact map of its clusters. `.metacool` is the format of
[metacooler](https://github.com/js2264/metacooler): a
[cooler](https://github.com/open2c/cooler) whose chromosomes are the contigs,
grouped into genomes. `cooler` reads it as any `.cool` file; metacooler adds
the genomes.

- **The map of the contigs**, made before any clustering, has every contig
  in `__unbinned__`. hicstuff, metacooler or metaTOR itself make it (see
  [From the reads to the map](#from-the-reads-to-the-map) and
  [From pairs or BAM files](#from-pairs-or-bam-files)).
- **The map of the clusters** is what `pipeline` (or `validation`) writes at
  its end: the same contacts, the contigs grouped by final cluster, with
  their summary, every clustering round and the network. It goes to the
  output directory, under the name of the input map, which is left as it
  was. `--inplace` updates the input map instead.

## What it holds

| Part | What it holds |
|---|---|
| Bins and contacts | Those of the input map, at its bin size: the raw contacts, without the filters and the normalization of the network |
| Genomes | One per final cluster, and `__unbinned__` for the other contigs (all of them, for a run with no final cluster). The contigs come genome by genome, the longest first |
| Genomes table | The columns of `cluster_summary.txt` (completeness, redundancy, abundances...), the length, contig and bin counts of each genome, and its contact statistics (`qc_*`) |
| Contigs table | The network columns `GC_content`, `Hit`, `Shotgun_coverage` and `Restriction_site` (the last two when the run had them), the contact statistics within each contig (`qc_intra_contig_short`, `qc_intra_contig_long`), and miComplete's marker genes on each contig (`markers_Bact105`, `markers_Arch131`: see [Marker genes of the contigs](../formats.md#marker-genes-of-the-contigs)). The genome of a contig is its final cluster |
| Bins table | The weights of the normalization, `metator_<normalization>` (`--normalization`, `empirical_hit` by default; none for `None`), stored as `cooler` stores balancing weights |
| Clustering rounds | Every cluster of every round: see [Clustering rounds](#clustering-rounds) |
| Metadata | Under the `metator` key, one entry per step (see below), and the assembly name and command history of the input map |

**The normalization.** The normalized contact of two bins is their raw count
times both weights. Summed by pair of contigs, these are the edges of the
network that `partition` and `validation` read (edges of weight 0 are
dropped). The weights are computed from the `Hit`, `Shotgun_coverage` and
`Restriction_site` columns: metacooler drops them when an edit changes one
of these, and `partition` and `validation` then refuse the map until
`metator network` stores them again.

**The metadata.** Each entry of the `metator` key holds the version of
metaTOR, the command line and parameters of the run that made it (the step,
or `pipeline`), and what it read of the map (`read`).

- `network`: the normalization, and the map columns read in place of
  `--enzyme` and `--depth`.
- `partition`: the normalization and `network_digest`, the SHA-256 of what
  the network is computed from, which `validation` checks against the
  network stored in the map.
- `validation`: how `cluster_summary.txt` is written from the genomes table
  (`summary_rows`, `summary_columns`, `summary_formats`: see
  [Cluster summary](../formats.md#cluster-summary)).
- `markers`: miComplete's version and thresholds, how the genes were found,
  and the markers and weights of both marker sets (see
  [Marker genes of the contigs](../formats.md#marker-genes-of-the-contigs)).

**Contact statistics.** Those of the input map (the `qc` block of the
metadata, the `qc_*` columns) are carried over: the counts within each
contig are kept, and the contacts between contigs are counted again for the
new genomes.

```python
import cooler
import metacooler

# The map of the clusters, which metator pipeline sample.metacool -o binning wrote
clr = cooler.Cooler("binning/sample.metacool")
mc = metacooler.MetaCooler("binning/sample.metacool")
mc.genomes  # one row per cluster, and __unbinned__, with the columns of cluster_summary.txt
mc.contigs  # one row per contig, with its GC content and Hi-C hits
clr.matrix(balance="metator_empirical_hit")[:]  # the contacts, normalized as the network
```

```bash
metacooler info binning/sample.metacool
metacooler dump -t genomes -H binning/sample.metacool
# The tables of the run (contig_data_final.txt, cluster_summary.txt,
# binning.txt, clusters.txt, cluster_members.txt and network.txt)
metator dump binning/sample.metacool -o tables
# The network of the contigs, normalized as metaTOR's (without the contacts
# within each contig, which metacooler writes too and are no edge)
metacooler network binning/sample.metacool --level contig --balance metator_empirical_hit --header \
    | awk -F '\t' '$1 != $2 && $3 != 0'
```

## Writing the map

`metacooler reassign` writes the map of the clusters anew from the input
map, in an edit session of metacooler. It is renamed into place only once
it is validated in full, so a run that stops, for any reason, leaves the
input map as it was; and a map updated run after run does not grow. For
locks, staged copies and links, see
[Writing a map](../cli.md#conventions).

- From an entry of a `.mmetacool` collection, the map of the clusters is
  `NAME.metacool`. `--inplace` is refused for an entry.
- A map that an earlier run left in the output directory is replaced only
  with `-F`.
- When the output directory holds the input map, the run updates it in
  place: this is how the steps run one by one (see
  [Step by step](#step-by-step)).
- A run that updates a map in place first warns of what it replaces there:
  clusters, clustering rounds, and hosts of MGE clusters.
- When the map cannot be written, the run logs the reason and stops. Its
  other outputs are kept, and the input map is left as it was: run again
  from it, with `-F`.

**Names.** Contig and cluster names must be ASCII, as `cooler` reads them as
ASCII. `ND`, `NA`, `NaN`, `N/A`, `N/D`, `-1` and `__unbinned__`, in any
case, mark a contig in no cluster: no cluster may be named so. metaTOR's
own names never are (`--prefix` takes letters, digits, `_` and `-`).

## The map of a 1.x run

A run whose `partition` and `validation` were not given the map (metaTOR 1.x)
has no map of its clusters. `metacooler reassign` writes it from the map of
the contigs, `contig_data_final.txt` and the summary of the clusters
(`bin_summary.txt` for 1.x), and `metacooler annotate` adds the columns that
metaTOR measured on the contigs (see [Other tools](../other-tools.md)).
`--contigs-format metator` reads the cluster of each contig in
`Final_cluster`, then in 1.x's `Final_bin`.

```bash
metacooler reassign sample.metacool --contigs binning/contig_data_final.txt \
    --contigs-format metator --genomes binning/bin_summary.txt \
    --genomes-format 1 -o binning/clusters.metacool
metacooler annotate binning/clusters.metacool --contigs binning/contig_data_final.txt \
    --contigs-format Name --inplace \
    --columns contigs.GC_content,contigs.Hit,contigs.Shotgun_coverage,contigs.Restriction_site
```

## Starting from a contact map

`metator network` and `metator pipeline` take the map of the contigs as
their first argument: a file, or an entry of a `.mmetacool` collection
(`file.mmetacool::entries/<name>`). `partition` and `validation` take it too,
to run the steps one by one (see [Step by step](#step-by-step)).
`network` and `pipeline` can also make it from pairs or BAM files (see
[From pairs or BAM files](#from-pairs-or-bam-files)).

```bash
metator pipeline sample.metacool -a assembly.fa -d depth.txt -e DpnII -o binning
```

A map that is missing, or that metacooler cannot open as one map (a
collection rather than one of its entries), is refused before any work.

- **Any bins.** One bin per restriction fragment (as hicstuff makes it),
  fixed bins (`metacooler create pairs --bins 1000`) or one bin per contig.
  The network and the hits come from the pairs that the map holds, whatever
  its bins.
- **Contigs of the assembly.** The contigs of the map must be contigs of the
  assembly, with their lengths; GC contents come from the assembly. The
  contigs of the run are those of the map, in its order, and their IDs are
  their positions in it. A map of other contigs or lengths, or that stores
  both triangles of the matrix, stops the run.
- **Missing contigs.** A map may lack contigs of the assembly, as with a
  contig filter (`metacooler create pairs --filter 'contigs.length >= 1kb'`).
  Those are in no table of the run, the log counts them, and the shotgun
  abundance of a cluster is its share of the sequenced bases of the contigs
  of the map.
- **Restriction sites and coverage** can come from the map: its
  `Restriction_site` and `Shotgun_coverage` columns, which metaTOR stores
  with `--enzyme` and `--depth`. A run from such a map needs neither option,
  and gives the same network and tables.
    - `--enzyme` and `--depth` fill in a missing column, or replace one, with
      a warning that counts the contigs whose value differs. The `network`
      metadata entry records which columns were read (`read`).
    - A column is read only if it gives every contig a valid value: at least
      1 restriction fragment, a coverage of at least 0. A column of `-` (what
      metaTOR writes for what it did not measure) counts as none; one that
      lacks some contigs is skipped, with a warning that names them.
    - A normalization that needs a column that neither the map nor the
      options give (`RS`, `abundance`, `theoritical_hit`) is refused before
      the run replaces anything.

## Step by step

The steps of `pipeline MAP` are commands of their own, which take the map as
their first argument: the map is what they share. Run in one output
directory, `network` writes the map there, and `partition` and `validation`
update it in place. The result is the map that `pipeline MAP` writes from
the same clusters, but for the dates and command lines it stores.

```bash
# 1. The network, stored in binning/sample.metacool; sample.metacool is left
#    as it was.
metator network sample.metacool -a assembly.fa -d depth.txt -e DpnII -o binning

# 2. The partition of the network of the map: the map stores its core and
#    overlapping clusters, in place.
metator partition binning/sample.metacool -o binning

# 3. The validation: the final clusters. The map becomes the map of the clusters.
metator validation binning/sample.metacool -a assembly.fa -o binning

# 4. The marker genes of the contigs, found in their clusters, in place.
metator markers binning/sample.metacool -a assembly.fa -o binning
```

| Step | What it stores in the map |
|---|---|
| `network MAP` | The network: the `GC_content`, `Hit`, `Shotgun_coverage` and `Restriction_site` columns, the weights of its normalization, and the `network` metadata entry, replacing those of an earlier network. Clusters and genomes are left as they were |
| `partition MAP` | The rounds `core` and `overlapping`, replacing the rounds the map held, and the `partition` metadata entry. No miComplete columns, no table, no FASTA file. Genomes and contigs are left as they were until `validation` |
| `validation MAP` | The map of the clusters, as `pipeline MAP` writes it: the final clusters as genomes, every round, the miComplete columns of the overlapping clusters it checked, and the network again, with the contigs grouped by cluster. The FASTA files of the overlapping clusters that reach `--size` go to its temporary directory, for miComplete |
| `markers MAP` | The marker genes of each contig, found in its final cluster: the columns `markers_Bact105` and `markers_Arch131` and the `markers` metadata entry, which every later run keeps but `markers` and `pipeline`, which find them again |
| `mge MAP --reassign`, then `host MAP` | The MGE clusters, as genomes, and their hosts, as `pipeline MAP --mges` does (see [mge on the map](#mge-on-the-map)) |

- **Seeds.** Louvain and Leiden draw other random numbers at each run, so
  steps run separately differ from `pipeline` as two `pipeline` runs do.
  `partition --seed S`, then `validation --seed S`, give the clusters of
  `pipeline --seed S`.
- **Checks.** `partition` and `validation` refuse a map without the network
  that `network` stores (or whose weights or columns an edit dropped).
  `validation` refuses a map without the partition, or whose partition was
  made on another network: `network MAP` run again keeps the partition, with
  a warning, until `partition MAP` partitions the new network. Both refuse a
  map of the output directory that is one of their inputs, and an invalid
  `--prefix`.
- **Contig IDs** are their positions in the map as the validation reads it:
  a validation run again on the map of the clusters numbers them by
  cluster.
- **Writing.** Each step writes the map as `pipeline MAP` does (see
  [Writing the map](#writing-the-map)), and refuses what `pipeline MAP`
  refuses before any work.
- **No tables.** metaTOR 2.0 removed the runs of `partition` and
  `validation` from the tables of the network and of the partition (`-n`,
  `-c`, `-f`): both need the map.

## Clustering rounds

A run reaches its final clusters in rounds:

| Round | Clusters |
|---|---|
| `core` | The core clusters of the partition, `core_00001`...: contigs that fall together at every iteration of Louvain or Leiden. No parent |
| `overlapping` | The overlapping clusters, `metator_00001_00000`..., each merging several core clusters. Those that reach `--size` hold miComplete's columns |
| `recursive_1`, `recursive_2`... | The clusters that each step of the validation makes by partitioning a contaminated cluster, their parent, again |
| `final` | The final clusters, each the cluster of the same name of the round it comes from |
| `mge` | The MGE clusters of `metator mge` (`mge_random` with `--random`): see [mge on the map](#mge-on-the-map) |

A contig is in at most one cluster of a round. The map stores the rounds in
its `/clusters` group (metacooler 0.4.0): `metacooler dump -t clusters` and
`-t members` write them out, and `MetaCooler.clusters` and
`MetaCooler.cluster_members()` read them. `metator dump MAP -o DIR` writes
them as `clusters.txt` and `cluster_members.txt` (see
[Clustering rounds](../formats.md#clustering-rounds)), as the validation does
in its temporary directory (`-N` keeps it).

A map holds the rounds of the run that wrote it, and those only:
`partition MAP` stores `core` and `overlapping`, then `validation MAP` every
round.

## mge on the map

`metator mge MAP` clusters the MGE contigs of a map once `metator network`
stored its network there: the map of the contigs, or of the clusters. It
reads no cluster. It reads the contacts between contigs from the pixels, and
those within a contig from its contact statistics (reads more than 1 kb
apart, which hicstuff and `metacooler create pairs` store with metacooler
0.4.0 or later).

- It stores its MGE clusters as the round `mge` (`mge_random` with
  `--random`), after the others, and the run as the `mge` metadata entry.
  Genomes and contigs are left as they were, but with `--reassign`, which
  makes each MGE cluster a genome (see [MGE genomes](mge.md#mge-genomes)).
- `metator host MAP` then stores the hosts: the columns `n_hosts`, `hosts`
  and `host_shares` of the round `mge`, and `n_mges` and `mges` of the
  genomes table (see [The hosts](mge.md#the-hosts)).
- A run that replaces the clusters of the map (`partition`, `validation`,
  `pipeline`) drops the rounds `mge` and `mge_random`, with a warning: run
  `mge` again on its map.

The map is written as `partition MAP` writes it (see
[In the contact map of the run](mge.md#in-the-contact-map-of-the-run)).

## From the reads to the map

hicstuff makes the map of a library from its reads: one bin per restriction
fragment, every contig in `__unbinned__`. metacooler sums the maps of
several libraries:

```bash
# 1. The map of each library, without -f or -D: they leave out pairs (3C events
#    that hicstuff filters, PCR duplicates) that metaTOR 1.x kept.
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff_L1 -P L1 \
    L1_R1.fq.gz L1_R2.fq.gz
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff_L2 -P L2 \
    L2_R1.fq.gz L2_R2.fq.gz

# 2. One map of both libraries, their contacts summed.
metacooler merge hicstuff_L1/L1.metacool hicstuff_L2/L2.metacool -o sample.metacool

# 3. Cluster the contigs: binning/sample.metacool, the map of the clusters.
metator pipeline sample.metacool -a assembly.fa -e DpnII -o binning
```

!!! note "`-M metacool` and mapping mode `5spm` need hicstuff 4.0.0"

    hicstuff 3.3.0 and earlier write no `.metacool` map and have no mapping
    mode `5spm`: with them, make the map from their pairs files with
    `metacooler create pairs`.

hicstuff aligns the reads with bowtie2 by default. `--mapping 5spm` aligns
them with `bwa mem -5SPM`, close to metaTOR 1.x's `--aligner bwa`
(`bwa mem -5SP`): on the test reads of 1.x, it finds the same 28,510
contacts between contigs as 1.x did.

`metacooler create pairs` makes the map from the pairs files of hicstuff
(`<prefix>.valid_idx.pairs.gz`) or of another aligner:

```bash
# One map of both libraries, one bin per contig, every contig in __unbinned__.
metacooler create pairs --pairs hicstuff_L1/L1.valid_idx.pairs.gz \
    --pairs hicstuff_L2/L2.valid_idx.pairs.gz --no-ps \
    --assembly-name assembly -o sample.metacool
```

| Option | Effect |
|---|---|
| no `--bins` | One bin per contig, the lengths from the `#chromsize` lines of the pairs files or from `--fasta`: the map that `network` and `pipeline` make from pairs files. Needs metacooler 0.5.0; with 0.4.0, give a BED file of the contigs, `--bins contigs.bed` (below) |
| `--bins 1000` | Bins of 1 kb |
| `--no-ps` | Do not write `sample.ps.tsv`, a P(s) table of `__unbinned__` |
| `--assembly-name` | The name of the assembly of the map (`unknown` otherwise) |
| `--contig-order fasta` | Keep the contigs in the order of the assembly (by default, the longest first). Needs `--fasta` or `#chromsize` lines; `--fasta` writes the index of the assembly (`.fai`) next to it |

With metacooler 0.4.0, the BED file of the contigs, from the `#chromsize`
lines of a pairs file of hicstuff:

```bash
zcat hicstuff_L1/L1.valid_idx.pairs.gz \
    | awk -v OFS='\t' '!/^#/ {exit} /^#chromsize:/ {print $2, 0, $3}' > contigs.bed
```

## From pairs or BAM files

`metator network` and `metator pipeline` also make the map of the contigs
themselves, from the pairs files or the BAM files of the libraries given in
place of the map: one bin per contig, every contig in `__unbinned__`, the
libraries summed, with the contact statistics that `metator mge` reads. See
[The inputs of a run](../cli.md#inputs) for what each argument may be.

```bash
# binning/sample.metacool, made from the pairs of both libraries, then the map of the clusters
metator pipeline hicstuff_L1/L1.valid_idx.pairs.gz,hicstuff_L2/L2.valid_idx.pairs.gz -p sample \
    -a assembly.fa -e DpnII -o binning

# binning/assembly.metacool, from BAM files of both mates (bwa mem -5SP), one per library
metator pipeline L1.bam L2.bam -a assembly.fa -e DpnII -o binning
```

- **The map** is `<prefix>.metacool` (`pipeline --prefix`), or is named
  after the assembly (`assembly.metacool` from `assembly.fa.gz`), in the
  output directory. The run then updates it in place, as the steps one by one
  do (`metator network L1.bam L2.bam -o binning`, then
  `metator partition binning/assembly.metacool -o binning`...). A run that
  stops later leaves it there, to start again without reading the pairs.
- **A map already there** is replaced only with `-F`. Inputs that give no
  contact on the contigs of the assembly leave no map.
- **From BAM files,** a pair is the primary alignments of both mates of a
  read, each with a mapping quality of at least `--min-quality` (30), on
  contigs of the assembly, flagged neither as a duplicate nor as failing QC,
  each at the 5' end of its alignment (see [Pairs](../formats.md#pairs)).
- **The sequences** of the BAM files (`@SQ`) and of the pairs files
  (`#chromsize`) must be contigs of the assembly, with their lengths (see
  [The inputs of a run](../cli.md#inputs)).
