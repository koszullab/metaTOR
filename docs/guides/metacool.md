# The contact map of the clusters

A run starts from a `.metacool` contact map of the contigs, made before any
clustering: every contig in `__unbinned__`, with its contacts. hicstuff makes
it from the reads of a library (`hicstuff pipeline -M metacool`), and
metacooler from pairs files (`metacooler create pairs`), or from the maps of
several libraries, which `metacooler merge` sums: see
[From the reads to the map](#from-the-reads-to-the-map). `metator network`
and `metator pipeline` read the network of the contigs from it (see
[Starting from a contact map](#starting-from-a-contact-map)). At its end,
`pipeline` writes the contact map of all the contigs, grouped by cluster, with
its final clusters, their summary, every clustering round and the network of
the contigs, with `metacooler reassign` from the map as the run read it: to
its output directory, under the name of the map, which it leaves as it was
(`--inplace` writes it over the map instead). It is written through
metacooler, anew from the map, so a map updated in place holds at most the
space that the edits of its last run left behind: it does not grow with each
run. A `.metacool` file is the contact map format of
[metacooler](https://github.com/js2264/metacooler): a
[cooler](https://github.com/open2c/cooler) whose chromosomes are the contigs,
grouped by genome. `cooler` and the tools built on it read it as any `.cool`
file, and metacooler adds the genomes.

## What it holds

- **The bins of the map the run started from, and one genome per final
  cluster.** The contigs that are in no final cluster are kept, with all
  their contacts, in the
  `__unbinned__` genome: all the contigs, for a run with no final cluster. The
  contigs come genome by genome, the longest first.
- **The contacts** of the map the run started from: the raw contacts, without
  the filters and the normalisation of the network.
- **The genomes table** holds the columns of `cluster_summary.txt` (completeness,
  redundancy, Hi-C and shotgun abundances…), next to the length, contig and bin
  counts of each genome and the contact statistics that metacooler computed
  while making the map.
- **The contigs table** (the chromosomes table of `cooler`) holds what
  metaTOR measured on each contig, the network of the contigs, last: the
  `GC_content`, `Hit`, `Shotgun_coverage` and `Restriction_site` columns,
  but those that the run did not measure (`Shotgun_coverage` of a run without
  `--depth`, `Restriction_site` of a run without `--enzyme`, when the map
  that a run starts from does not give them either), and the contact
  statistics that metacooler computes within each contig
  (`qc_intra_contig_short`, `qc_intra_contig_long`). The genome of a contig
  is its final cluster.
- **The normalization of the network**, which `metator network` stores in
  the map, as `metator pipeline` does: the weights `metator_<normalization>`
  of its bins (`--normalization`, `empirical_hit` by default; none for
  `None`), stored the way `cooler` stores balancing weights. The normalized
  contact of two bins is their raw count times both weights, so the contacts
  normalized and summed by pair of contigs are the edges of the network that
  `partition` and `validation` read, but those of weight 0 (a contig without
  hits, or without coverage, has no weight). They are computed from the
  `Hit`, `Shotgun_coverage` and `Restriction_site` columns of the contigs
  table, which the column records, so that metacooler drops it when an edit
  changes one of them: `partition` and `validation` then refuse the map, and
  `metator network` stores them again. A normalization that metacooler fails
  to store fails the map, which is then not written, with the reason.
- **The metadata** hold, under the `metator` key, one entry per step, each
  with the version of metaTOR, the command line and its parameters of the
  run that made it (the step, or `pipeline`), and what it read of the map
  (`read`: columns of the contigs table, weights of the bins, clustering
  rounds). The `network` entry is the run that stored the network, with its
  normalization (and, under `read`, the columns of the map that it read in
  place of `--enzyme` and `--depth`); the `partition` entry, the run that
  stored the partition, with the normalization and the identity of the
  network that it partitioned (`network_digest`: the SHA-256 of what the
  network is computed from, by contig name), which `validation` checks
  against the network stored in the map; the `validation` entry, the run
  that wrote the map of the clusters, with how the summary of the clusters,
  `cluster_summary.txt`, is written from the genomes table: its rows
  (`summary_rows`, the clusters in order), its columns (`summary_columns`)
  and the format of each (`summary_formats`, as printf writes it: `%.2f`,
  `%.4f`, `%d` or `%s`; a missing value is an empty field).

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

## The map of a run made without it

A run of an earlier version whose `partition` and `validation` were not given
the map has no map of its clusters: `metacooler reassign` writes it from the map of its contigs,
`contig_data_final.txt` and the summary of the clusters, and `metacooler
annotate` adds the columns that metaTOR measured on the contigs (see
[Other tools](../other-tools.md)). For the tables of a run of 1.x, the summary
is `bin_summary.txt`, and `--contigs-format metator` reads the cluster of each
contig in `Final_cluster`, then in 1.x's `Final_bin` (the summary of a 2.0
development build is `cluster_summary.txt`):

```bash
metacooler reassign sample.metacool --contigs binning/contig_data_final.txt \
    --contigs-format metator --genomes binning/bin_summary.txt \
    --genomes-format 1 -o binning/clusters.metacool
metacooler annotate binning/clusters.metacool --contigs binning/contig_data_final.txt \
    --contigs-format Name --inplace \
    --columns contigs.GC_content,contigs.Hit,contigs.Shotgun_coverage,contigs.Restriction_site
```

## Names and failures

The contig and cluster names must be ASCII, as `cooler` reads them as ASCII
text. No cluster may be named as metacooler marks a contig in no cluster:
`NA`, `NaN`, `N/A`, `N/D`, `-1` or `__unbinned__` in any case, and `nd` in any
case other than `ND`, which marks a contig in no cluster. metaTOR's own cluster
names never are (`--prefix` takes letters, digits, `_` and `-`). The clusters
of a run are clusters of the contigs of its map, by their names in it: a
contig of a final FASTA file that the run cannot find back by its name
stops the run, as a bug of metaTOR would.

When the map of the clusters cannot be written, `metator pipeline` logs the
reason as an error and stops with it, and keeps its other outputs. The map
that the run started from is then left as it was, as it is when the run stops
with an error before (in the partition or the validation, for example): a new
run can start from it again, with `-F`, which replaces the final clusters
that the stopped run left in the output directory.

## Starting from a contact map

`metator network` and `metator pipeline` build the network from a
`.metacool` contact map of the contigs of the assembly, their first argument,
which they require (one file, or an entry of a `.mmetacool` collection,
`file.mmetacool::entries/<name>`). `metator partition` and `metator
validation` take it as their first argument too, to run the steps one by one
on the map (see [Step by step](#step-by-step)):

```bash
metator pipeline sample.metacool -a assembly.fa -d depth.txt -e DpnII -o binning
```

A first argument that is not there, or that metacooler cannot open as one map
(a collection rather than one of its entries, which the error names), is
refused before any work.

- **Any bins**: one bin per restriction fragment (as hicstuff makes the
  map), fixed bins (`metacooler create pairs --bins 1000`), or one bin per
  contig. The contacts between two contigs are those between their bins, and
  the hits of a contig those of its bins: the network and the hits are those
  of the pairs the map holds, whatever its bins. On the test reads of
  metaTOR 1.x, the map of hicstuff (bowtie2, one bin per restriction
  fragment) and the map that metacooler makes from
  `tests_data/alignment.pairs` (1 kb bins) give the same network.
- **Its contigs must be contigs of the assembly**, with their lengths: the
  GC contents still come from the assembly. The contigs of the run are those
  of the map, in its order (longest first, as hicstuff and metacooler make
  it), and their IDs in the tables of the run their positions in it. It may
  lack contigs of the assembly, as a map made with a contig filter does
  (`metacooler create pairs --filter 'contigs.length >= 1kb'`, which also
  leaves out the pairs with a read on the contigs that it drops): those have
  no contacts, and are in no table of the run, and the log says how many they
  are; the shotgun abundance of a cluster is its share of the sequenced bases
  of the contigs of the map. A map of other contigs or lengths, or that stores
  both triangles of the matrix, stops the run with an error.
- **Its restriction sites and coverage**: the `Restriction_site` and
  `Shotgun_coverage` columns of its contigs table, which metaTOR stores in
  the maps it writes with `--enzyme` and `--depth`, give the restriction
  sites and the shotgun coverage of the contigs. So a run from such a map
  needs neither option: it gives the same network and tables, byte for byte
  (with pandas 2, for coverages of at most 15 significant digits and 15
  digits after the decimal point, as jgi's depth files have: metacooler can
  store a longer one one unit off in its last place, which changes the last
  digits of the contig data files and of a network normalized by it, and
  has `--depth` warn that it replaces the map's value).
  - `--enzyme` and `--depth` fill in a column that the map lacks, and replace
    one that it has, with a warning that counts the contigs whose value
    differs. The map of the clusters stores what the run used, and its metadata
    which columns it read (`read`, under the `network` entry).
  - A column is read only if it gives every contig of the assembly a value: a
    number of restriction sites (plus one) of at least 1, a coverage of at
    least 0. A contig of no length, which no map holds, gets 1 and 0. A
    column with no valid value, such as one of `-`, as metaTOR writes what it
    did not measure, is as no column. One that lacks some contigs is not
    read, with a warning that names them, and the map of the clusters does not
    carry it.
  - A normalization that needs a column that neither the map nor the options
    give (`RS`, `abundance`, `theoritical_hit`) is refused before the run
    replaces anything: before any work when the map has no such column, once
    the assembly is read when the column lacks contigs.
At its end, `pipeline` writes the map of its clusters to its output directory,
under the name of the input map, which `metacooler reassign` writes anew from
the input map (`--inplace` updates the input map in place with it instead):

- the bins and contacts of the input map, at its bin size; one genome per
  final cluster and `__unbinned__` for the other contigs, grouped genome by
  genome, the longest first;
- the columns of `cluster_summary.txt`, in the genomes table, and the
  network of the contigs, last in the contigs and bins tables: they replace
  those the input map had (the map of an earlier clustering, for example),
  and so do every clustering round of the run;
- the contact statistics (the `qc` block of the metadata and the `qc_*`
  columns of the genomes and contigs tables) of the input map, when it has
  them, carried over to the clusters: metacooler keeps the counts within each
  contig (measured on the reads, for a map made from pairs files) and counts
  the contacts between contigs again. A map written by metacooler 0.2.0 holds
  no count per contig: its statistics are computed again from the pixels,
  which hold no read positions (at one bin per contig, every contact within a
  contig then counts as short-range, and the numbers of pairs read and left
  out are lost, which metacooler warns about for a map made from pairs files);
- the assembly name and the history of the commands of the input map; what
  metacooler's regroup does not carry (an extra group, or an extra
  attribute of the file or of its tables) is left out, with a warning that
  names it.

The map of the clusters is renamed into place (over the input map, in place)
only once it is written and validated in full. metaTOR writes it in an edit
session of metacooler (`metacooler.staged_edit`): every edit of the step goes
into one copy of the input map, staged in a private directory next to the map
that the run writes (`.<map>.<uuid>.<n>.d/`), which metacooler validates in
full, then renames onto that map. The copy starts fresh, so the map holds at
most the space that the edits of the last run left behind (HDF5 does not
reclaim the space of what an edit replaces): a map updated run after run does
not grow.
When it cannot be written, the run stops with the reason, and the input map is
left as it was, as it is when the run stops for any other reason (Ctrl-C
included). A map updated in place through a symbolic link is updated where
the link points: the link keeps pointing at the new map, and the lock and the
staged copy are next to the file it points to. A map updated in place keeps
its mode, its owner and group, and its extended attributes (POSIX ACLs
included), where the run may set them: a warning names those it cannot keep
(the owner of a map of another user, written through its group, for
example). The directory of the map that the run writes must be writable, with
room for the copy of the edit session and for the new map that `metacooler
reassign` writes from the input map, plus cooler's temporary files; a map
updated in place without write permission is refused before any work.
The map of the output directory is a map of its own, which leaves the input
map as it was: from an entry of a `.mmetacool` collection, a map named after
the entry (`NAME.metacool`), as a run does not update a collection in place
(`--inplace` is refused for it). One that an earlier run left there is
replaced only with `-F`; a symbolic link there is replaced, not its target
(even one that points to the input map), and a hard link of the input map
(the same file under another name) is refused before any work. When the
output directory holds the input map itself, the run updates it in place: the
steps of a run, one by one, in one output directory. A run that updates a map
in place (as the map of its output directory, or with `--inplace`) warns,
before any work, of what it replaces there that an earlier run left (or
another binner): the clusters of the map (`pipeline MAP`, `validation MAP`)
and its clustering rounds (the rounds of the binning for `pipeline MAP` and
`partition MAP`, those of the validation for `validation MAP`, the round
`mge` for `mge MAP`, `mge_random` for `mge MAP --random`), and the hosts of
its MGE clusters (`host MAP`, and the runs that replace what they describe:
`mge MAP`, `partition MAP`, `validation MAP`, `pipeline MAP`). The map keeps no
record of them: a run that writes its map to another output directory, without
`--inplace`, leaves them in the input map.

A run holds metacooler's lock of the map it writes (`<map>.lock`, next to it)
from its checks to its end, and reads the map under it: a second run, or a
metacooler command that writes the map in place (`annotate --inplace`,
`reassign --inplace`, `normalize`, the `metadata` edits), stops before any
work, naming the process that holds it, and the files to delete if it is
gone. A run that writes the map of its output directory locks that map, its
own name: runs from one map, each in its own output directory, run side by
side, and the input map is not locked. A map that changed since the run read it (another run
replaced it, or a program that writes it without metacooler's lock, `cooler
balance` or an h5py script, wrote it) is refused when the run writes its map,
rather than mixed with the clusters of the run, read from the map as it was:
the run stops with the reason, its other outputs kept, and the map keeps that
change; run it again. A command that reads the map to check it (its contigs,
its columns) stops so too when the map is replaced as it reads it. A command that reads a map
(`network MAP`, the input map of a run that writes its map to its output
directory) stops with an error when another run replaces it, or writes it, as
it reads it. metacooler
reads a map one chunk of pixels at a time, and a map updated in place has its
bins in another order: such a read would mix the two maps. metacooler checks
the file of the map before and after each read. Run it again once the other
run has ended.

Neither metaTOR nor metacooler handles signals. A run stopped by Ctrl-C
removes its staged copy and its lock; its temporary directory, with the
working directory of miComplete in it, stays, as after any failure. miComplete runs in the process group of the run: a Ctrl-C in the
terminal, or the signal that a job scheduler sends to the job, reaches it too.
A run killed by SIGTERM (a job scheduler's time limit, or a cancelled job) or
SIGKILL (`kill -9`), or on a machine that goes down, leaves its lock, and its
staging directory next to the map. The next run, or metacooler command, that
writes the map on the same host takes the lock over, with a warning, and
removes every staging directory of the map that another writer left, with a
warning each. A lock of another host, or of another machine or PID namespace
(a container) under the same host name, is never taken over: once that run is
over, delete the lock file and the staging directories that the error names.
The map itself is left as it was.

### Step by step

The three steps of `pipeline MAP` are commands of their own, which take the
map as their first argument too: the map is the network that they share.
Run one after the other in one output directory, the first writes the map
there, and the next ones read it and update it in place, as `pipeline MAP`
does, and write the map that it writes from the same clusters,
but for the date the map was written and the command lines and parameters
that it stores (those of `network`, `partition` and `validation`). Louvain draws
random numbers: two runs give the same clusters only from the same random
stream, as the tests check in one process, so steps run as separate commands
differ from `pipeline` as two `pipeline` runs do. `partition MAP` on a map
that holds the clusters of an earlier run replaces its clustering rounds with
its own, and the `partition` entry of the `metator` key of its metadata with
what the run was (the network it partitioned among it), and leaves its
genomes, contigs and the rest of its metadata as the earlier run wrote them,
until `validation MAP` writes the new ones. `validation MAP` refuses a
partition made on another network than the one stored in the map:
`network MAP` run again keeps the partition, with a warning, until
`partition MAP` partitions the new network:

```bash
# 1. The network, stored in binning/sample.metacool; sample.metacool is left
#    as it was.
metator network sample.metacool -a assembly.fa -d depth.txt -e DpnII -o binning

# 2. The partition of the network of the map: the map stores its core and
#    overlapping clusters, in place.
metator partition binning/sample.metacool -o binning

# 3. The validation: the final clusters. The map becomes the map of the clusters.
metator validation binning/sample.metacool -a assembly.fa -o binning
```

- `network MAP` stores the network of the contigs in the map: the
  `GC_content`, `Hit`, `Shotgun_coverage` and `Restriction_site` columns of
  its contigs table, the weights of its normalization, and the `network`
  entry of the `metator` key of its metadata, in place of those of an earlier
  network. Its clusters and genomes are left as they were.
- `partition MAP` stores its partition: every core and overlapping cluster,
  with its contigs, as the rounds `core` and `overlapping` of the clustering
  rounds of the map (see [Clustering rounds](#clustering-rounds)), in place
  of the rounds that it held, without the columns of miComplete, which the
  validation runs, and what the run was, as the `partition` entry of the
  `metator` key of its metadata. It writes no table and no FASTA file. The
  genomes and contigs tables of the map are left as they were.
- `validation MAP` reads the partition from the map, writes the FASTA files
  of the overlapping clusters that reach its `--size` from the assembly, in
  its temporary directory, for miComplete, and writes the map of the
  clusters, as `pipeline MAP` does at its end: the final clusters are its
  genomes, and it holds every round (`core`, `overlapping`, `recursive_<k>`,
  `final`), the overlapping clusters that miComplete checked with its
  columns, and the network of the map (its columns and the weights of its
  normalization), which it stores again once it has grouped the contigs by
  cluster, last in the contigs and bins tables, as `pipeline MAP` stores it.
  The IDs of the contigs are
  their positions in the map as the validation reads it: a validation run
  again on the map of the clusters, whose contigs are grouped by cluster,
  numbers them in that order.

`network MAP`, `partition MAP` and `validation MAP` write the map as
`pipeline MAP` does: to their output directory, through a map written next to
it and validated in full, which leaves `MAP` as it was (give the next step
that map), or in place, when `MAP` is the map of their output directory or
with `--inplace`. They lock
the map while they run, and refuse what `pipeline MAP` refuses, before any
work and before `-F` replaces the directories of clusters of an earlier run:
a map that is not there or that metacooler cannot open, a lock that another
run holds, an entry of a collection with `--inplace`, and a map of an earlier
run in the output directory, without `-F`. `partition` and `validation` refuse a map without the
network that `metator network` stores (or one whose weights or columns an
edit dropped since), and `validation` a map without the partition that
`metator partition` stores (the rounds `core` and `overlapping`, their
clusters named as `partition` names them); they refuse a map of the output
directory that is one of their inputs, and an invalid `--prefix`. Both need the map:
metaTOR 2.0 removed their runs from the tables of the network and of the
partition (`-n`, `-c`, `-f`).

## Clustering rounds

A run reaches its final clusters in rounds:

- `core`: the core clusters of the partition, `core_00001`..., the contigs
  that fall together at every iteration of Louvain or Leiden, with no parent:
  an overlapping cluster merges several of them;
- `overlapping`: every overlapping cluster of the partition,
  `metator_00001_00000`...; miComplete checks those that reach `--size`, which
  hold its columns;
- `recursive_1`, `recursive_2`...: the clusters that each step of the
  validation makes by partitioning a contaminated cluster again, its parent;
- `final`: the final clusters, each the cluster of the same name of the round
  it comes from.

`metator mge` adds one more after them, `mge`: the MGE clusters (or
`mge_random`, a random clustering of the same shape, with `--random`; see
[mge on the map](#mge-on-the-map)).

The validation writes them all next to its other tables, in its temporary
directory, which `-N` keeps, and `metator dump MAP -o DIR` writes them from
the map of the clusters: `clusters.txt`, one row per cluster of every round
(its round and name, the round and name of its parent, `NA` for none, its
number of contigs and their length, `n_contigs` and `length`, as the map
computes them, and the columns of miComplete), and `cluster_members.txt`, one
row per contig of a cluster (its round, cluster and contig). A contig is in at
most one cluster of a round.

The map of the clusters stores both in its `/clusters` group (metacooler
0.4.0), next to the final clusters, which its genomes and contigs tables hold:
`metacooler dump -t clusters` and `-t members` write them out, and
`MetaCooler.clusters` and `MetaCooler.cluster_members()` read them. A map
holds the rounds of the run that wrote it, and those only: `metator
partition MAP` stores the rounds `core` and `overlapping` alone, which
`metator validation MAP` reads, until it stores every round (see
[Step by step](#step-by-step)).
`metator mge` then adds the round `mge`, or replaces the one of an earlier
`mge` run, and keeps the others; and `metator host` adds the hosts of its MGE
clusters as columns of the round `mge` (`n_hosts`, `hosts`, `host_shares`),
and the MGE clusters of each final cluster as columns of the genomes table
(`n_mges`, `mges`).

## mge on the map

`metator mge MAP` clusters the MGE contigs of a map that `metator network`
stored its network in (the map of the contigs after `network`, or the map of
the clusters) from its contacts: between two contigs, its pixels, summed by
pair of contigs; within a contig, between reads more than 1 kb apart, its
contact statistics (metacooler's QC block, counted from the pairs at
1001 bp, which `hicstuff pipeline -M metacool` and `metacooler create pairs`
store with metacooler 0.4.0 or later). It reads no cluster of the map. It
stores its MGE clusters in the map, as the round `mge` of its clustering
rounds, after the others, and the run in the `mge` entry of its metadata
(`mge_random` for both with `--random`); the genomes and contigs tables are
left as they were. A run that replaces the clusters of the map drops the
rounds `mge` and `mge_random`, and their entries, with a warning: run `mge`
again on its map. The map is written as `partition MAP` writes it: to the
output directory, or in place (see
[Clustering mobile genetic elements](mge.md#in-the-contact-map-of-the-run)).
`metator host MAP` then finds the MAGs that host the MGE clusters, and
stores them in the map, as `mge` stores its MGE clusters (see
[The hosts](mge.md#the-hosts)).

## From the reads to the map

hicstuff makes the map to start from, with no clustering yet, from the reads
of a library: one bin per restriction fragment, and every contig in
`__unbinned__`. metacooler sums the maps of several libraries:

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

hicstuff also writes the pairs of each library, `<prefix>.valid_idx.pairs.gz`,
from which `metacooler create pairs` makes a map too, as it does from the
pairs files of another aligner:

```bash
# One map of both libraries, one bin per contig, every contig in __unbinned__.
python -c 'import sys; from Bio import SeqIO; [print(r.id, 0, len(r), sep="\t") for r in SeqIO.parse(sys.argv[1], "fasta")]' \
    assembly.fa > contigs.bed
metacooler create pairs --pairs hicstuff_L1/L1.valid_idx.pairs.gz \
    --pairs hicstuff_L2/L2.valid_idx.pairs.gz --bins contigs.bed --no-ps \
    --assembly-name assembly -o sample.metacool
```

!!! note "`-M metacool` and mapping mode `5spm` need hicstuff 4.0.0"

    hicstuff 3.3.0 and earlier write no `.metacool` map and have no mapping
    mode `5spm`: with them, make the map from their pairs files with
    `metacooler create pairs`.

- `--bins 1000` makes bins of 1 kb instead, taking the lengths of the contigs
  from the `#chromsize` lines of hicstuff's pairs files, or from `--fasta`.
- `--no-ps`: otherwise metacooler also writes `sample.ps.tsv`, a P(s) table of
  the single `__unbinned__` genome.
- `--assembly-name`: otherwise the assembly of the map, and of the map of the
  clusters, is `unknown`.
- `--contig-order fasta` keeps the contigs in the order of the assembly (by
  default, the longest first). It needs `--fasta` or the `#chromsize` lines of
  the pairs files: metacooler does not keep the order of a BED file.
  `--fasta` writes the index of the assembly (`.fai`) next to it, in a
  directory that must be writable.

hicstuff aligns the reads with bowtie2 by default; `--mapping 5spm` aligns the
mates with `bwa mem -5SPM`, close to metaTOR 1.x's `--aligner bwa`, which ran
`bwa mem -5SP`. On the test reads of metaTOR 1.x, `--mapping 5spm` finds the
same 28,510 contacts between contigs as metaTOR 1.x with `--aligner bwa`, and
bowtie2 the 41,563 contacts of `tests_data/alignment.pairs`, 18,632 of them
between contigs.
