# The contact map of the clusters

`metator pipeline` writes the contact map of all the contigs, grouped by
cluster, as `metator.metacool` in the output directory (`<prefix>.metacool`
with `--prefix`). It reads the pairs files into the map once, as soon as it
has them: the map of the contigs, in no cluster yet (every contig in
`__unbinned__`, the longest first), with their contacts, the contact statistics
that metacooler measures on the reads and what metaTOR measures on the assembly
(`GC_content`, and `Shotgun_coverage` and `Restriction_site` with `--depth`
and `--enzyme`). A run with one library reads its network from this map, as a
run from a map does (see
[Starting from a contact map](#starting-from-a-contact-map)):
`network.txt` and `contig_data_network.txt` are those of its pairs files, byte
for byte, as `pipeline` counts every pair, whatever the input. It still
writes `network_0.txt`, the network of its library, as a copy
of `network.txt`: its pairs files give the same bytes. A run with several
libraries builds the map of each library, the entries `library_0`,
`library_1`... of `metator.mmetacool` (`<prefix>.mmetacool`), and the map of
their sum, `metator.metacool`, merged from them (their contacts and contact
statistics summed; its counts are stored as 64-bit integers): it reads
`network.txt` from the sum, and the network of each library (`network_N.txt`)
and `hit_data_alignment.txt` from its map, the files that its pairs files give,
byte for byte. So a pair on a contig that the assembly lacks is left out, as
metacooler warns (it counts them, with the pairs that have a read past the end
of their contig), and the run goes on, whatever the number of libraries. At the
end of the run, `pipeline` writes the map of the clusters over the map of
the contigs, from the final tables, with `metacooler reassign`; `--no-metacool`
removes the maps instead (the `.mmetacool` collection too). The maps of the
libraries stay maps of the contigs, in no cluster, without the columns measured on
the assembly, which the sum holds. It is written through a fresh file, which holds none of the space that the edits making it leave behind, so a map updated in place does not grow with each run. A `.metacool` file is the contact map format of
[metacooler](https://github.com/js2264/metacooler): a
[cooler](https://github.com/open2c/cooler) whose chromosomes are the contigs,
grouped by genome. `cooler` and the tools built on it read it as any `.cool`
file, and metacooler adds the genomes.

## What it holds

- **One bin per contig, and one genome per final cluster.** The contigs that
  are in no final cluster (`ND` in `contig_data_final.txt`) are kept, with
  all their contacts, in the `__unbinned__` genome: all the contigs, for a run
  with no final cluster. The contigs come genome by genome, the longest first.
- **The contacts.** Each pair whose two ends lie on contigs of the assembly is
  one contact between these two contigs (an intra-contig contact when both
  ends are on the same contig). These are the raw contacts of the pairs
  files, without the filters and the normalisation of the network, and those
  of several libraries are summed.
- **The genomes table** holds the columns of `cluster_summary.txt` (completeness,
  redundancy, Hi-C and shotgun abundances…), next to the length, contig and bin
  counts of each genome and the contact statistics that metacooler computes
  while reading the pairs.
- **The contigs table** (the chromosomes table of `cooler`) holds what
  metaTOR measured on each contig: the `GC_content`, `Hit`,
  `Shotgun_coverage` and `Restriction_site` columns of
  `contig_data_final.txt`, except those that are `-` on every contig
  (`Shotgun_coverage` of a run without `--depth`, `Restriction_site` of a run
  without `--enzyme`, when the map that a run starts from does not give them
  either), and the contact statistics that metacooler computes
  within each contig (`qc_intra_contig_short`, `qc_intra_contig_long`). The
  genome of a contig is its final cluster.
- **The normalization of the network**, in the map of a `metator pipeline`
  run: the weights `metator_<normalization>` of its bins (`--normalization`,
  `empirical_hit` by default), stored the way `cooler` stores balancing
  weights. The normalized contact of two bins is their raw count times both
  weights, so the contacts normalized and summed by pair of contigs are the
  edges of `network.txt`, but those of weight 0 (a contig without hits, or
  without coverage, has no weight). They are computed from the `Hit`,
  `Shotgun_coverage` and `Restriction_site` columns of the contigs table,
  which the column records, so that metacooler drops it when an edit changes
  one of them. A run whose map lacks a column that the normalization reads
  (`Restriction_site` without `--enzyme`) stores none, with a warning, and so
  does `--normalization None`, without one. It needs
  `metacooler.store_weights`, from metacooler's next release: with metacooler
  0.3.0, the map stores none, which the log says.
- **The metadata** hold the version of metaTOR, the command line and its
  parameters, under the `metator` key.

```python
import cooler
import metacooler

clr = cooler.Cooler("binning/metator.metacool")
mc = metacooler.MetaCooler("binning/metator.metacool")
mc.genomes  # one row per cluster, and __unbinned__, with the columns of cluster_summary.txt
mc.contigs  # one row per contig, with its GC content and Hi-C hits
clr.matrix(balance="metator_empirical_hit")[:]  # the contacts, normalized as the network
```

```bash
metacooler info binning/metator.metacool
metacooler dump -t genomes -H binning/metator.metacool
# The network of the contigs, normalized as metaTOR's
metacooler network binning/metator.metacool --level contig --balance metator_empirical_hit
```

## The map of a run made without it

A run whose `partition` and `validation` were not given the map has no map of
its clusters: `metacooler reassign` writes it from the map of its contigs,
`contig_data_final.txt` and `cluster_summary.txt`, and `metacooler annotate`
adds the columns that metaTOR measured on the contigs (see
[Other tools](../other-tools.md)):

```bash
metacooler reassign sample.metacool --contigs binning/contig_data_final.txt \
    --contigs-format Name,Final_cluster --genomes binning/cluster_summary.txt \
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
names never are, but those of a hand-edited `contig_data_final.txt` can be.
Each row of `cluster_summary.txt` must have as many fields as its header, as
metaTOR writes it.

When the map of the clusters cannot be written, `metator pipeline` logs the
reason as an error and stops with it, and keeps its other outputs. `pipeline` removes a map that
an earlier run left in the output directory before it builds the map of the
contigs in its place (`--no-metacool` or not), as it no longer matches the
clusters: from then on, the map there is this run's, or there is none. So a run
from reads, alignments or pairs refuses, before any work, an input that is
that map (or a symbolic link to it or from it), and, from pairs, a contact
map: start from a map by giving it as the first argument. A map of the contigs
that cannot be built (for example, from pairs files that metacooler cannot
read, see below) stops the run before its network, `--no-metacool` or not,
whatever the number of libraries: the reason is logged as an error, in the log
file of the run, and no map is left. When the map of the clusters cannot be
written over it (for example, from a misaligned `cluster_summary.txt`), the
map of the contigs stays, in no cluster. So it does when the run stops with an
error once it is built (in the partition or the validation, for example),
`--no-metacool` or not: a new run can start from it, without aligning the
reads again, with `metator pipeline <outdir>/metator.metacool`
(`<outdir>/<prefix>.metacool` with `--prefix`) and `-F`, which replaces the
cluster directories that the stopped run left in the output directory
(see [Starting from a contact map](#starting-from-a-contact-map)).

`pipeline` checks what the map of the contigs needs before any work. Its
contig names must be ASCII: a run that is not from a map refuses an assembly
with other contig names, which it gives, before it builds the index of the
assembly or aligns a read; rename these contigs to ASCII names, in the
assembly and in the alignments or pairs files made from it (and build its
index again, for an assembly given as an index). metacooler, which builds the
map, finds the compression of a pairs file from its first bytes, whatever its
name: it reads plain, gzip, bzip2, zstd and LZ4 files. So a run from pairs
refuses, before any work, an empty pairs file, and a zip archive, which
metacooler cannot read: decompress the file, or compress it with gzip.
metacooler also needs the same number of
fields on every row of a pairs file; its header lines (those that start with
`#`) may come anywhere, so pairs files of one assembly concatenated with `cat`
are one library. Give the pairs files of several libraries comma-separated,
one file per library. A run from pairs that metacooler cannot read stops with
the reason and these needs.

## Starting from a contact map

`metator network` and `metator pipeline` can build the network from a
`.metacool` contact map of the contigs of the assembly, instead of reads,
alignments or pairs: the map is their first argument (one file, or an entry
of a `.mmetacool` collection, `file.mmetacool::entries/<name>`), in place of
`-1`, `-2` and `--start`. `metator partition` and `metator validation` take
it as their first argument too, to run the steps one by one on the map (see
[Step by step](#step-by-step)):

```bash
metator pipeline sample.metacool -a assembly.fa -d depth.txt -e DpnII -o binning
```

A first argument that is not there, or that metacooler cannot open as one map
(a collection rather than one of its entries, which the error names), is
refused before any work.

- **Any bin size**, from one bin per contig (as the maps written by metaTOR)
  to fixed bins. The contacts between two contigs are those between their
  bins, and the hits of a contig those of its bins: the network and
  `contig_data_network.txt` are those of the pairs the map holds, byte for
  byte, as both count every pair.
- **Its contigs must be contigs of the assembly**, with their lengths: the
  contig IDs and GC contents still come from the assembly. It may lack some
  of them, as a map made with a contig filter does (`metacooler create pairs
  --filter 'contigs.length >= 1kb'`, which also leaves out the pairs with a
  read on the contigs that it drops): those get no contacts and no hits, and
  the log says how many they are. A map of other contigs or lengths, or that
  stores both triangles of the matrix, stops the run with an error.
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
    that `pipeline` writes which columns it read (`columns_from_map`).
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
- **Its libraries are already summed**: no `network_<i>.txt` is written.

At the end of such a run, `pipeline` updates the input map in place with the
map of its clusters, with `metacooler reassign` (`--output-map FILE` writes it
to FILE instead):

- the bins and contacts of the input map, at its bin size; one genome per
  final cluster and `__unbinned__` for the other contigs, grouped genome by
  genome, the longest first, as in the map of a run from pairs;
- the columns of `cluster_summary.txt` and `contig_data_final.txt`, last, as in
  the map of a run from pairs: they replace those the input map had (the map
  of an earlier clustering, for example);
- the contact statistics (the `qc` block of the metadata and the `qc_*`
  columns of the genomes and contigs tables) of the input map, when it has
  them, carried over to the clusters: metacooler keeps the counts within each
  contig (measured on the reads, for a map made from pairs files) and counts
  the contacts between contigs again. A map written by metacooler 0.2.0 holds
  no count per contig: its statistics are computed again from the pixels,
  which hold no read positions (at one bin per contig, every contact within a
  contig then counts as short-range, and the numbers of pairs read and left
  out are lost, which metacooler warns about for a map made from pairs files);
- the assembly name and the history of the commands of the input map.

The input map is replaced only once the map of the clusters is written next to
it and validated in full. When it cannot be, the run stops with the reason, and
the input map is left as it was, as it is when the run stops for any other
reason (Ctrl-C included). A symbolic link is replaced, not its target.
`--output-map FILE` leaves the input map as it was: a run from an entry of a
`.mmetacool` collection needs it, as a run does not update a collection in
place, and a FILE that exists is replaced only with `-F`. A run that updates a
map in place warns, before any work, of what it replaces there that an earlier
run left (or another binner): the clusters of the map (`pipeline MAP`,
`validation MAP`) and its clustering rounds (the rounds of the binning for
`pipeline MAP` and `partition MAP`, those of the validation for `validation
MAP`, the round `mge` for `mge --map`). The map keeps no record of them:
`--output-map` writes the map of the run to another file instead.

While a run writes a map, a lock file next to it (`<map>.lock`) holds it: a
second run that would write the same map stops before any work, naming the run
that holds it. A lock that a run left behind on the same host when it died is
removed, with a warning; one from another host is not: remove it once that
run is over. Only the runs that write a map take its lock: a command that reads
a map (`network MAP`, the input map of a run with `--output-map`) stops with an error when another run replaces it,
or writes it, as it reads it. metacooler opens the map for each chunk of
pixels that it reads, and a map updated in place has its bins in another
order: such a read mixed the two maps, and gave wrong contacts. Run it again
once the other run has ended.

A run stopped by Ctrl-C or SIGTERM (a job scheduler's time limit, or a
cancelled job) removes its lock and the files of the map it was writing next
to the map, and ends miComplete, with the processes that it started, if it
runs; its temporary directory, and the working directory of miComplete in the
output directory, stay, as after any failure. A run killed outright (`kill
-9`, or a machine that goes down) cleans nothing up, nor does one stopped
where Python cannot stop it cleanly (as a pool of processes starts or ends: it
then dies of the signal at once): besides its lock, it can leave the files of
the map it was writing next to the map, which can be deleted:
`.<map>.<random>.tmp` or `.<map>.<random>.tmp.<random>.fresh`, and those that
metacooler and cooler write while it is written
(`.<map>.<random>.tmp.partial`, `.<map>.<random>.tmp.partial.__writing__`,
`tmp<random>.multi.cool`). The map itself is left as it was.

### Step by step

The three steps of `pipeline MAP` are commands of their own, which take the
map as their first argument too. Run one after the other on one map, they
read and update it as `pipeline MAP` does, and write the map that it writes
from the same clusters, but for the date the map was written and the command
line and parameters that it stores (those of `validation`). Louvain draws
random numbers: two runs give the same clusters only from the same random
stream, as the tests check in one process, so steps run as separate commands
differ from `pipeline` as two `pipeline` runs do. `partition MAP` on a map
that holds the clusters of an earlier run replaces its clustering rounds with
its own, and leaves its genomes, contigs and the `metator` key of its
metadata as the earlier run wrote them, until `validation MAP` writes the new
ones:

```bash
# 1. The network: network.txt and contig_data_network.txt. The map is only read.
metator network sample.metacool -a assembly.fa -d depth.txt -e DpnII -o binning

# 2. The partition: overlapping_clusters/ and contig_data_partition.txt. The map
#    stores the overlapping clusters.
metator partition sample.metacool -a assembly.fa -c binning/contig_data_network.txt \
    -n binning/network.txt -o binning

# 3. The validation: the final clusters. The map becomes the map of the clusters.
metator validation sample.metacool -a assembly.fa -c binning/contig_data_partition.txt \
    -n binning/network.txt -f binning/overlapping_clusters -o binning
```

- `network MAP` reads the map, and leaves it as it was.
- `partition MAP` stores its overlapping clusters that reach `--size`, those
  of `overlapping_clusters/`, with their contigs, as the round `overlapping` of the
  clustering rounds of the map (see [Clustering rounds](#clustering-rounds)),
  in place of the rounds that it held, without the columns of miComplete,
  which the validation runs. The genomes and contigs tables of the map are
  left as they were.
- `validation MAP` writes the map of the clusters, as `pipeline MAP` does at
  its end: the final clusters are its genomes, and it holds every round
  (`overlapping`, `recursive_<k>`, `final`).

`partition MAP` and `validation MAP` write the map as `pipeline MAP` does: in
place, through a map written next to it and validated in full, or to
`--output-map FILE`, which leaves the map as it was (give the next step that
FILE; a run from an entry of a `.mmetacool` collection needs it). They lock
the map while they run, and refuse what `pipeline MAP` refuses, before any
work and before `-F` replaces the directories of clusters of an earlier run:
a map that is not there or that metacooler cannot open, a lock that another
run holds, an entry of a collection without `--output-map`, and a FILE that
exists, without `-F`. They also refuse a map that is no map of the contigs
of their contig data (`-c`: other contigs or lengths), an `--output-map`
that is one of their inputs, and an invalid `--prefix`. Without the map,
`partition` and `validation` run as before, from the tables of the earlier
steps, and refuse `--output-map`; an invalid `--prefix` is now refused
before `-F` replaces their directories of clusters, without a log file.

## Clustering rounds

A run reaches its final clusters in rounds:

- `overlapping`: the overlapping clusters of the partition that reach
  `--size`, which miComplete checks;
- `recursive_1`, `recursive_2`...: the clusters that each step of the
  validation makes by partitioning a contaminated cluster again, its parent;
- `final`: the final clusters, each the cluster of the same name of the round
  it comes from.

`metator mge --map` adds one more after them, `mge`: the MGE clusters (see
[mge on the map](#mge-on-the-map)).

The validation writes them all next to its other tables, in the output
directory: `clusters.txt`, one row per cluster of every round (its round and
name, the round and name of its parent, `NA` for none, and the columns of
miComplete), and `cluster_members.txt`, one row per contig of a cluster (its
round, cluster and contig). A contig is in at most one cluster of a round.

The map of the clusters stores both in its `/clusters` group (metacooler
0.4.0), next to the final clusters, which its genomes and contigs tables hold:
`metacooler dump -t clusters` and `-t members` write them out, and
`MetaCooler.clusters` and `MetaCooler.cluster_members()` read them. A map
holds the rounds of the run that wrote it, and those only: `metator
partition MAP` stores the round `overlapping` alone, until
`metator validation MAP` stores every round (see [Step by step](#step-by-step)).
`metator mge --map` then adds the round `mge`, or replaces the one of an
earlier `mge` run, and keeps the others.

## mge on the map

`metator mge --map MAP` counts the contacts within each MGE contig between
reads more than 1 kb apart, which a map does not tell: it reads the pairs
files, and refuses a map in their place. It stores its MGE clusters in the
map, as the round `mge` of its clustering rounds, after those of the run; the
genomes and contigs tables, the final clusters, are left as they were. A run
that replaces the clusters of the map drops the round `mge`, with a warning.
The map is written as `partition MAP` writes it: in place, or to
`--output-map FILE` (see
[Clustering mobile genetic elements](mge.md#in-the-contact-map-of-the-run)).

## From hicstuff to metaTOR through metacooler

A map to start from can be made from the pairs files of
[hicstuff](https://koszullab.github.io/hicstuff/), aligned in mapping mode
`5spm`, with no clustering yet:

```bash
# 1. Align each library, without -f or -D: they leave out pairs (3C events that
#    hicstuff filters, PCR duplicates) that metaTOR's own alignment keeps.
hicstuff pipeline --mapping 5spm -g assembly.fa -e DpnII -o hicstuff_L1 -P L1 \
    L1_R1.fq.gz L1_R2.fq.gz
hicstuff pipeline --mapping 5spm -g assembly.fa -e DpnII -o hicstuff_L2 -P L2 \
    L2_R1.fq.gz L2_R2.fq.gz

# 2. One map of both libraries, one bin per contig, every contig in __unbinned__.
samtools faidx assembly.fa
awk -v OFS='\t' '{print $1, 0, $2}' assembly.fa.fai > contigs.bed
metacooler create pairs --pairs hicstuff_L1/L1.valid_idx.pairs.gz \
    --pairs hicstuff_L2/L2.valid_idx.pairs.gz --bins contigs.bed --no-ps \
    --assembly-name assembly -o sample.metacool

# 3. Cluster the contigs: sample.metacool becomes the map of the clusters, in place.
metator pipeline sample.metacool -a assembly.fa -e DpnII -o binning
```

!!! note "Mapping mode `5spm` needs a hicstuff release after 3.3.0"

    hicstuff 3.3.0, the latest release so far, has no mapping mode `5spm`:
    with it, leave `--mapping 5spm` out.

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

On the test library of the [Quickstart](../quickstart.md), this route finds
the same 28,510 contacts between contigs as metaTOR's own alignment with
`--aligner bwa`.
