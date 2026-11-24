# File formats

The files metaTOR reads and writes. Every table is tab-separated.

## Inputs

### Assembly

A FASTA file of the contigs (`--assembly`), gzip-compressed or not: the
assembly the contact map was made from.

- A contig is read by its name, the first word of its header. An assembly
  that names a contig in several records is refused.
- `--enzyme` counts restriction sites in upper and lower case, so the
  repeats of a soft-masked assembly count too.

### Depth

The output of MetaBAT2's `jgi_summarize_bam_contig_depths` (`--depth`): a
header line, then one line per contig.

```text
contigName	contigLen	totalAvgDepth	sample.bam	sample.bam-var
NODE_522	22786	4.76595	4.76595	15.5294
```

- metaTOR reads the first three columns: name, length and coverage.
- It must hold every contig of the assembly, with its length. The depth file
  of another assembly stops the run.
- Lengths of a million bases or more are written to 6 significant digits
  (`1.21754e+06` for 1,217,543 bases): they match the contig's length
  rounded the same way.

### Pairs

The 4D Nucleome
[pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md),
with 7 columns: read name, contig and position of each mate, and the strand
of each. hicstuff writes them with the map (`<prefix>.valid_idx.pairs.gz`),
sorted, bgzip-compressed, with two more columns (the restriction fragment of
each mate).

```text
## pairs format v1.0
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#chromsize: NODE_522 22786
...
H9:1:HGJYYBBXY:8:1102:3792:13693	NODE_100063	161	NODE_100063	402	-	-
```

`metacooler create pairs`, `metator network` and `metator pipeline` make a
map from them (see [The inputs of a run](cli.md#inputs)). They read columns
2 to 5, in any order, and the `#chromsize` lines, which must give the
contigs of the assembly their lengths.

From BAM files, a run writes the pairs of each library in its temporary
directory, `library_<i>.pairs`: read name `.`, each mate at the 5' end of
its alignment (1-based), the contig first in the assembly first, and no
`#chromsize` line.

```text
## pairs format v1.0
#shape: upper triangle
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2
.	NODE_100063	161	NODE_100063	402	-	-
```

### Lists of contigs

`metator mge --mges` reads one contig name per line (the first word of each
line).

## The contig data tables

`contig_data_final.txt` has one row per contig of the map, in its order,
with a header. [`metator dump`](cli.md#dump) writes it from the map of the
clusters. The validation writes its own in its temporary directory (`-N`
keeps it), with `contig_data_partition.txt`, the table of the partition.

| Column | Written by | What it holds |
|---|---|---|
| `ID` | network | The ID of the contig, from 1: its position in the map as the run read it, and its node in the network. `validation` and `pipeline` regroup the map, so join a table to another map's network by name (`metator.network.relabel_network`), not by ID. Validation tables only: `metator dump` leaves it out |
| `Name` | network | Its name |
| `Size` | network | Its length, in bp |
| `GC_content` | network | Its GC content, in percent |
| `Hit` | network | Its hits: the ends of the contacts of the map that lie on it (a contact within the contig counts twice) |
| `Shotgun_coverage` | network | Its coverage, from the depth file or the map's `Shotgun_coverage` column (when it gives every contig a value); `-` otherwise |
| `Restriction_site` | network | Its number of restriction sites plus one, from `--enzyme` or the map's `Restriction_site` column (when it gives every contig a value); `-` otherwise |
| `Core_cluster_ID`, `Core_cluster_contigs`, `Core_cluster_size` | partition | Its core cluster: ID, number of contigs and size in bp |
| `Overlapping_cluster_ID`, `Overlapping_cluster_contigs`, `Overlapping_cluster_size` | partition | Its overlapping cluster, likewise |
| `Recursive_cluster_ID`, `Recursive_cluster_contigs`, `Recursive_cluster_size` | validation | Its recursive cluster, likewise; `00000`, `-`, `-` when its cluster was not partitioned again |
| `Final_cluster` | validation | Its final cluster, or `ND` (not determined) |

```text
	ID	Name	Size	GC_content	Hit	Shotgun_coverage	Restriction_site	Core_cluster_ID	Core_cluster_contigs	Core_cluster_size	Overlapping_cluster_ID	Overlapping_cluster_contigs	Overlapping_cluster_size	Recursive_cluster_ID	Recursive_cluster_contigs	Recursive_cluster_size	Final_cluster
0	1	NODE_522	22786	62.0775915035548	1683	4.76595	162	00001	2	23301	00001	15	49021	00000	-	-	metator_00001_00000
1	2	NODE_1404	15015	64.2823842823842	952	4.00585	118	00002	1	15015	00002	5	19711	00000	-	-	ND
```

- The validation's table starts with an unnamed column (the row number, from
  0), then `ID`. That of `metator dump` starts with `Name`, its rows in the
  order of the map of the clusters, grouped by cluster.
- A contig with no edge in the network (no hits, or no contact with another
  contig) has `-` for its core and overlapping clusters.
- IDs have five digits (`00001`), but `Core_cluster_ID` is a plain number
  when every contig is in a core cluster.
- As inputs, the tables are read by column header: other columns (such as
  the bins of another binner) may be added, in any order.

## The network

`metator network` stores the network in the map (see
[Building the network](guides/network.md#building-the-network)).
[`metator dump`](cli.md#dump) writes it as `network.txt`, without a header,
one line per pair of contigs in contact, in map order (each contig, then the
contigs after it):

| Column | What it holds |
|---|---|
| `contig1`, `contig2` | The names of the two contigs (1.x wrote their IDs) |
| weight | The weight of the edge, with the normalization that `metator network` stored (the contact count for `None`) |

Contacts within a contig and edges of weight 0 are not edges. metacooler
exports the network of any map that holds it, but with both: drop them.

```console
metacooler network sample.metacool --level contig --balance metator_empirical_hit --header \
    | awk -F '\t' '$1 != $2 && $3 != 0' > network.tsv
```

`network.tsv` has the columns `node1`, `node2` and `value` (raw contact
counts without `--balance`). Without `--header`, there is no header line.

## The clusters

### Cluster summary

`cluster_summary.txt` has one row per final cluster, named in its first
column, with miComplete's columns first:

| Column | What it holds |
|---|---|
| `Length` | The size of the cluster, in bp |
| `GC-content` | Its GC content, in percent |
| `Present Markers` | The single-copy marker genes found in it |
| `Completeness`, `Redundancy` | The share of the markers found, and their mean number of copies (1 for no duplicate) |
| `Weighted completeness`, `Weighted redundancy` | The same, with miComplete's marker weights |
| `Contigs`, `N50`, `L50`, `N90`, `L90` | Its number of contigs, and the N50, L50, N90 and L90 of their sizes |
| `CDs` | Its number of coding sequences, found by Prodigal (pyrodigal) |
| `recursive` | Whether the cluster comes from the recursive partition |
| `step` | The recursive step that made it, 0 for an overlapping cluster |
| `parent` | The cluster it was split from (an overlapping cluster at step 1, a recursive cluster after), empty otherwise |
| `HiC_abundance` | The share of the hits of all the contigs that are on the cluster, in percent |
| `SG_abundance` | The share of the shotgun reads (size × coverage) of the contigs of the map that are on the cluster, in percent; with `--depth` only. With a contig filter, the clusters add up to 100% of the map, where 1.x counted the whole assembly |

```text
	Length	GC-content	Present Markers	Completeness	Redundancy	Weighted completeness	Weighted redundancy	Contigs	N50	L50	N90	L90	CDs	recursive	step	parent	HiC_abundance	SG_abundance
metator_00001_00000	49021.00	62.03	1.00	0.01	1.00	0.02	1.00	15.00	9019.00	2.00	1390.00	9.00	70.00	False	0		2.4736	2.6314
```

- Both marker sets, bacterial (105 markers) and archaeal (131), are tried;
  the one with the higher weighted completeness is kept.
- Contamination is the weighted redundancy minus one (see
  [Cluster quality](concepts.md#cluster-quality)).
- miComplete's numbers have two decimals (whole ones too), the abundances
  four, `step` none.
- The map's genomes table holds the same values, and the `validation` entry
  of its `metator` metadata how to write them: `summary_rows` (the clusters
  in order), `summary_columns` and `summary_formats` (printf formats: `%.2f`,
  `%.4f`, `%d`, `%s`; a missing value is an empty field).
  [`metator dump`](cli.md#dump) writes the table again from them, byte for
  byte.
- After `metator host`, `n_mges` and `mges` come last (see
  [Hosts of the MGE clusters](#hosts-of-the-mge-clusters)).

After [`metator mge --reassign`](guides/mge.md#mge-genomes), `metator dump`
writes the final clusters (without those whose contigs all went to MGE
genomes), then the MGE genomes, with four more columns. The final clusters
keep the values of their validation, which count the MGE contigs that they
lost.

| Column | What it holds |
|---|---|
| `GC` | The GC content of an MGE genome, in percent; empty for a final cluster |
| `Hit` | Its hits, the sum of those of its contigs |
| `BinningScore` | Its clustering score: the lowest score at which its contigs were merged; empty for a single contig |
| `mge` | `True` for an MGE genome, `False` for a final cluster |

### binning.txt

Each contig of a genome of the map and its genome, without a header: the
collection that anvi'o imports (see
[Curating MAGs with anvi'o](guides/curation.md)). The genomes are the final
clusters, plus the MGE genomes after `metator mge --reassign`.

```text
NODE_522	metator_00001_00000
NODE_2487	metator_00007_00000
```

### miComplete results

`overlapping_micomplete_results.txt` and `recursive_micomplete_step_<n>.txt`,
in the validation's temporary directory, hold miComplete's results on the
overlapping clusters and on the recursive clusters of each step: three
comment lines (`##`), the columns of the cluster summary up to `CDs`, then
the marker set kept, `Bacteria` or `Archaea`.

### Clustering rounds

`clusters.txt` has one row per cluster of every round, in the order they
were made: `core`, `overlapping`, `recursive_<k>` (step `k` of the
validation) and `final`; from `metator dump`, also `mge` and `mge_random`
when `metator mge` stored them. Columns, under a header:

- `round`, `name`;
- `parent_round`, `parent`: the cluster it was split from, or, for a final
  cluster, the cluster it is; `NA` for none;
- `n_contigs`, `length`, as the map computes them from its contigs;
- miComplete's columns, `NA` where the cluster has none, and those of the
  MGE clusters (see [MGE clusters](#mge-clusters)).

```text
round	name	parent_round	parent	n_contigs	length	Length	GC-content	...	Completeness	...
overlapping	metator_00001_00000	NA	NA	31	2817353	2817353.0	52.30	...	0.9655	...
recursive_1	metator_00002_00001	overlapping	metator_00002_00000	18	1504321	1504321.0	41.80	...	0.7103	...
final	metator_00001_00000	overlapping	metator_00001_00000	31	2817353	2817353.0	52.30	...	0.9655	...
```

`cluster_members.txt` has one row per contig of a cluster: its round,
cluster and contig, under a header. The map stores both (see
[Clustering rounds](guides/metacool.md#clustering-rounds)).

### FASTA files of the clusters

| Directory | Written by | What it holds |
|---|---|---|
| `fasta/mags_with-mges/` | `validation`, `pipeline` | One file per final cluster, `<cluster>.fa`: the records of its contigs as in the assembly, MGE contigs included. A copy of the file that miComplete checked |
| `fasta/mags/` | `pipeline` | Each MAG (a genome of the map, but the MGE genomes) twice: `<MAG>.fa`, the records of its contigs in map order, and `<MAG>.stitched.fa`, one record named after it, its contigs longest first, joined by 200 `N`, in lines of 80 bases |
| `fasta/mges/` | `mge`, `pipeline --mges` | Each MGE cluster the same two ways, its contigs longest first, ties by name |

In its temporary directory (`-N` keeps it), the validation writes the
clusters that miComplete checks (`metacooler.write_fasta`): the overlapping
clusters in `overlapping_clusters/` (their contigs by core cluster, then by
ID), and the recursive clusters of step `n` in
`recursive_clusters/step_<n>/`. A gzipped
assembly is decompressed there once per run. No pyfastx index is made.

### Cool files of the MAGs

`cool/` (`pipeline`) has one [cooler](https://github.com/open2c/cooler)
file per MAG, `<MAG>.cool`, extracted from the map by `metacooler extract`:
the contacts of its contigs (with `--mges`, also of the MGE clusters it
hosts), at the bins of the map, each contig a chromosome. `cooler balance`
balances it with its defaults (the `weight` column,
`clr.matrix(balance=True)`); a file too small to balance keeps `NaN`
weights, with a warning.

## The contact map of the clusters

A `.metacool` file of [metacooler](https://github.com/js2264/metacooler): a
cooler whose contigs are grouped into genomes, one per final cluster. See
[The contact map of the clusters](guides/metacool.md).

## MGE clusters

`metator mge` (and `pipeline --mges`) stores its MGE clusters in the map as
the clustering round `mge` (`mge_random` for `--random`), one cluster per MGE
cluster (see
[Clustering mobile genetic elements](guides/mge.md#in-the-contact-map-of-the-run)):

| Column | What it holds |
|---|---|
| `name` | `MetaTOR_MGE_00001` and on, as its FASTA files in `fasta/mges/` |
| `length`, `n_contigs` | Its size and number of contigs, computed by metacooler |
| `GC`, `Hit` | Its GC content and hits |
| `BinningScore` | The lowest score at which its contigs were merged; missing for a single contig |
| `n_hosts`, `hosts`, `host_shares` | After `metator host`: its number of hosts, their names by decreasing share, and their shares in percent with two decimals, each list comma-separated (see [Hosts of the MGE clusters](#hosts-of-the-mge-clusters)) |

Its contigs are the members of the round (`mc.cluster_members("mge")`), in
map order. `metator dump` writes the round in `clusters.txt` and
`cluster_members.txt`.

## Hosts of the MGE clusters

`metator host` (and `pipeline --mges`) writes `tables/mge_hosts.txt`: one row
per MGE cluster and MAG in contact, the MGE clusters in the order of the
round `mge`, the MAGs of each by decreasing signal (see
[Clustering mobile genetic elements](guides/mge.md#the-hosts)). Numbers are
written in full.

| Column | What it holds |
|---|---|
| `mge` | The MGE cluster (`MetaTOR_MGE_00001`...) |
| `mag` | The MAG: a genome of the map |
| `signal` | The sum of the weights of the network edges between their contigs, with the normalization that `metator network` stored (the contact count for `None`) |
| `contigs` | The number of contigs of the MAG in contact with the MGE cluster |
| `share` | The signal, in percent of the signal of the MGE cluster with all the MAGs |
| `share_of_all` | The signal, in percent of the signal of the MGE cluster with every other contig (in a MAG, in none, or in another MGE cluster) |
| `host` | Whether the MAG hosts the MGE cluster (`True`, `False`): a `share` of at least `--min-share` and at least `--min-contigs` contigs |

```text
mge	mag	signal	contigs	share	share_of_all	host
MetaTOR_MGE_00013	metator_00005_00000	0.7090148667776763	37	93.8246167981684	76.06906868895449	True
MetaTOR_MGE_00013	metator_00020_00000	0.005037475316584434	4	0.6666139362589284	0.5404626529450921	False
```

`plots/mge_hosts.pdf` plots the signal of the hosts (log scale) and their
share. The map holds the hosts in the round `mge` (see
[MGE clusters](#mge-clusters)), and the MGE clusters of each MAG in its
genomes table, which `metator dump` writes last in `cluster_summary.txt`:

| Column | What it holds |
|---|---|
| `n_mges` | The number of MGE clusters that the MAG hosts; missing for a genome that is not a MAG (empty for the MGE genomes of `--reassign`) |
| `mges` | Their names, in the order of the round `mge`, comma-separated |

## Marker genes of the contigs

[`metator markers`](cli.md#markers), and the step `markers` of `pipeline`,
store in the map the marker genes of miComplete that each contig carries, as
miComplete finds them in its genome (see
[Concepts](concepts.md#marker-genes-of-the-contigs)): two columns of the
contigs table, and an entry of the metadata. Summed over the contigs of a
genome of the map, as the run found their genes, they give the four numbers
of miComplete that the [cluster summary](#cluster-summary) holds for the
final clusters: completeness, redundancy, weighted completeness and weighted
redundancy. Summed over any other group of contigs, they give an estimate of
them (see *Recomputing miComplete's numbers*, below).

| Column | What it holds |
|---|---|
| `markers_Bact105` | The markers of the contig in miComplete's bacterial set (105 markers) |
| `markers_Arch131` | Those in its archaeal set (131 markers, 35 of them in both sets) |

Each holds ASCII text: the markers of the contig sorted by name, in byte
order (upper case first), comma-separated, each followed by `*` and its
count when that is not 1. The count is the number of proteins of the contig
that miComplete keeps as hits of the marker. An empty text is no marker.

```text
name	markers_Bact105	markers_Arch131
NODE_522	Ribosomal_L14,Ribosomal_L22,SecY*2	Ribosomal_L14,Ribosomal_L22,SecY*2
NODE_1404
NODE_2487	Cons_hypoth95*0,Ribosomal_S9	Ribosomal_S9
```

- A marker name matches `[A-Za-z0-9_-]+`, of at most 15 characters.
- `*0`: miComplete sees a hit of the marker on the contig (a line of the
  table of hmmsearch, `--tblout`) but drops every one, its bias too high for
  its score or its best domain too weak. The marker counts in the weighted
  completeness of a bin, as miComplete counts every marker it sees, and
  nowhere else.
- Both columns hold text, a fixed-length ASCII string dataset, even when no
  contig of the map has a marker of a set: an empty text on every contig.
- Each column is as wide as its longest text, that of the contig with the
  most markers (852 and 547 bytes on a map of 20,566 contigs; 105 KB and
  73 KB on disk, compressed). cooler, when it opens a map, and metacooler,
  when it reads its contigs table, decode every column of the table at 4
  bytes a character: on that map, the two columns make each take about 0.3 s
  and 150 MB more, and `metator export`, which reads the map once per MAG,
  70 s in place of 47 s. The cost grows with the number of contigs (some 3 s
  and 1.5 GB a read for 200,000 contigs).
- `metacooler dump -t contigs -H MAP` writes both columns with the rest of
  the contigs table; the run writes them in its temporary directory too
  (`contig_markers.txt`, which `-N` keeps).

The entry `metator` → `markers` of the metadata records the run, as the
entries of the other steps do (`version`, `command`, `parameters`, `read`),
and:

| Key | What it holds |
|---|---|
| `micomplete` | The version of miComplete whose marker sets and rules the run applied |
| `thresholds` | `evalue`, `bias` and `domain_cutoff`: miComplete's `--evalue`, `--bias` and `--domain-cutoff` |
| `genes` | How the genes were found: the versions of `pyrodigal` and `pyhmmer`; `single` and `metagenome`, the rule of each of pyrodigal's two modes, in words (see [`metator markers`](cli.md#markers)); `trained_genomes`, the genomes pyrodigal was trained on, by name, in the order of the genomes table (single mode: those of `min_training` bp, 20000, or more, with the linkers between their contigs); `metagenome_contigs`, the number of contigs in metagenome mode (in no genome, or in a shorter one); `metagenome_z`, the number of proteins their E-values were counted over: the median number of proteins of the genomes trained in single mode, every one of them, MGE genomes included (2,111 on a map of the mock community MM133 whose 66 trained genomes are 50 MAGs and 16 MGE genomes, where the median of the MAGs is 2,974.5; without a trained genome, the proteins of these contigs, as one bin; `null` without such contigs) |
| `sets` | By set (`Bact105`, `Arch131`): its `column`, its `markers` in the order of its HMMs, and their `weights`, in the same order, as miComplete reads them |
| `kept` | The set that metaTOR keeps for a bin: `by` the higher `weighted_completeness`, `ties` to `Bact105` |
| `assignment_digest` | The binning that the genes were found in: the SHA-256, in hexadecimal, of the name of each contig of the map and of its genome then (`__unbinned__` for none), the contigs sorted by name, each name followed by a NUL byte (`metator.metacool.assignment_digest`). The same whatever the order of the contigs; another once a contig is in another genome, or a genome is renamed |

Every other run keeps the markers, those that regroup the contigs too
(`validation`, `mge --reassign`, `metacooler reassign`). A metaTOR run that
regroups them into other genomes than those of `assignment_digest` says once,
in its log, that the markers predate its new genomes, and `metator report`
says so too (`markers_predate_genomes`, in the section `run` of
[the report](#the-report-of-a-run); `metator.metacool.markers_predate`): each
contig keeps its markers, but summed over a genome whose contigs changed,
they give an estimate. `metator markers`, run again, finds them in the
genomes as they are.

### Recomputing miComplete's numbers

For a group of contigs (a genome of the map, or any curated bin, whose numbers
are an estimate: see below), in each set (`metator.markers.marker_quality`):

1. Sum the counts of each marker over the contigs. A marker written without
   `*` counts 1; one at `*0` counts 0 but is seen.
2. P is the number of markers of a sum of 1 or more, N the number of markers
   of the set. Without such a marker (P = 0), all four numbers are 0.
3. Completeness = P / N, and redundancy = the sum of the counts / P, each
   rounded to 4 decimals as Python's `round` does: to the nearest, an exact
   tie to the even digit (33 / 32 = 1.03125 gives 1.0312, where JavaScript's
   `toFixed(4)` gives 1.0313).
4. WC is the sum of the weights of the markers seen (those of a sum of 0
   too). The weighted redundancy is (the sum of the weights of the markers of
   a sum above 1, each once, + WC) / WC, rounded to 4 decimals; the weighted
   completeness is WC rounded to 4 decimals, 0.0001 if that is 0.

metaTOR keeps the numbers of Bact105, unless those of Arch131 have the higher
weighted completeness (`Bacteria` or `Archaea` in the column `Markers` of
miComplete's results; `metator.markers.kept_quality`).

The markers, their weights and the columns are those that the entry records
(`sets`): those of the miComplete that found the markers, not the files of
the miComplete installed now, which may weigh its markers otherwise.
`metator.markers.stored_quality` does it all from the map, with the sets of
its entry (`metator.markers.map_sets`; `quality_of` counts the texts of a
group with them):

```python
import metacooler
import metator.markers as mtmk

mc = metacooler.MetaCooler("binning/sample.metacool")
group = ["NODE_522", "NODE_2487"]
print(mtmk.stored_quality(mc, group))

# The same, set by set, from the entry of the map:
sets = mc.metadata["metator"]["markers"]["sets"]
contigs = mc.contigs.set_index("name")
for name, spec in sets.items():
    counts = mtmk.sum_markers(mtmk.decode_markers(text) for text in contigs.loc[group, spec["column"]])
    print(name, mtmk.marker_quality(counts, spec["markers"], spec["weights"]))
```

Over a genome whose genes the run found in single mode (trained on its
contigs, in map order), the numbers are exactly those that miComplete gives
the FASTA file of the genome in map order, as `fasta/mags/` and
`metacooler fasta` write it. The cluster summary holds those of the clusters
as the validation checked them, with their contigs in the order of its
partition, on which pyrodigal's training depends a little: the numbers are
the same but for a rare gene, and for a final cluster that `mge --reassign`
took contigs from, which the summary still counts (sum the markers over the
contigs of the cluster, the round `final`, for those).

Over any other group of contigs (a curated bin: two MAGs merged, a contig
moved from one MAG to another, contigs recruited from `__unbinned__`), the
numbers are an estimate, from the genes found in the genomes of the map, and
not those that miComplete gives the FASTA file of the group: miComplete
trains pyrodigal on the group itself, which finds other genes, and counts
the E-values of the hits over the proteins of the group. On the mock
community MM133 (20,566 contigs), of 28 curated groups (the first half of
the contigs of 12 MAGs, 8 merges of two MAGs, 8 MAGs with the 30 longest
contigs of `__unbinned__`), 23 came out exact, and 5 differed (the four
numbers, in the set that metaTOR keeps, Bact105 for all five):

| Curated group | From the stored markers | From miComplete |
|---|---|---|
| Two MAGs of two organisms merged (GC 66.5% and 39.0%) | 1.0, 1.6667, 1.0, 1.5427 | 0.9905, 1.1827, 0.9915, 1.1059 |
| Two other MAGs merged | 0.9905, 1.375, 0.9918, 1.4454 | 0.9905, 1.3654, 0.9918, 1.4454 |
| The first half of the contigs of a MAG | 0.981, 1.0291, 0.9767, 1.0228 | 0.981, 1.0194, 0.9767, 1.0228 |
| A MAG with 30 contigs of `__unbinned__` | 0.9905, 1.4135, 0.9814, 1.3415 | 0.9905, 1.4038, 0.9814, 1.3214 |
| Another MAG with 30 contigs of `__unbinned__` | 0.9905, 1.3654, 0.9918, 1.299 | 1.0, 1.3619, 1.0, 1.2966 |

The merge of two organisms parts the most: the markers stored on its
contigs count 175 copies of the 105 markers of Bact105, those of the two
genomes, where miComplete, trained on the two organisms together, finds 123
copies of 104 of them. Its redundancy says less of the contamination than
the markers do.

## The report of a run

[`metator report`](cli.md#report), `export` and `pipeline` write the report
of a run from the map of its clusters, in the output directory:
`<name>.report.txt`, `<name>.report_mqc.yaml` and `<name>.report.html`.

`<name>` is the run's map without `.metacool`: `sample` for
`sample.metacool`, `NAME` for an entry `file.mmetacool::entries/NAME`. A map
made from pairs or BAM files is `<prefix>.metacool` (`pipeline --prefix`) or
`<assembly>.metacool`. MultiQC names its rows by `<name>`, so name the maps
after their samples (`-p` of `pipeline` from pairs or BAM files, `-P` of
`hicstuff pipeline`).

### `<name>.report.txt` { #reporttxt }

A text file of sections, then tables, which a script or a MultiQC module
parses line by line (`metator.report.read_txt` reads it back):

- A line that starts with `#` is a comment. The first is
  `# metaTOR report of <name>`.
- `[name]` starts a section: each line is a key, spaces, then its value (the
  rest of the line).
- `[table: name]` starts a table: a header line, then one row per line, its
  fields separated by spaces (no field holds a space).
- `NA` is a missing value. Integers are written as such; percents (keys or
  columns ending in `_pct`) with 2 decimals; ratios and fractions (`_ratio`,
  `_fraction`, from 0 to 1) with 4; lengths in bp end in `_bp`; other numbers
  in full.
- `format_version` (section `report`) is 1. It changes when a key, column or
  unit changes or goes, not when one is added: a parser should skip what it
  does not know.
- Sections, keys and tables that do not apply are left out: `mges` and
  `hosts`, and their tables, before `metator mge` and `metator host`.

| Section | Keys |
|---|---|
| `report` | `format_version`; `sample` (the name of the run); `map`; `metator_version` (that wrote the report); `written`; `map_written` (when the map was written). |
| `run` | `steps`, the steps the map records, in order, comma-separated (`network`, `partition`, `validation`, `markers`, `mge`, `host`); `input` (the map, or the pairs or BAM files, comma-separated), `assembly`, `enzyme`, `depth` (`NA`: the map gave them); `normalization` (network); `algorithm`, `iterations`, `overlap_pct`, `res_param` (partition); `size_bp`, `rec_iter`, `rec_overlap_pct` (validation); `seed`, `rec_seed`, the `--seed` of the partition and of the validation (`NA`: none); `mge_threshold_bin` (`mge`); `host_min_share_pct`, `host_min_contigs` (`host`); `markers_predate_genomes` (`markers`): `true` when the marker genes of the contigs were found in another binning than that of the map (see [Marker genes of the contigs](#marker-genes-of-the-contigs)), `false` when in its own, `NA` when the map does not record it. |
| `contigs` | `contigs`, `contigs_bp` (the contigs of the map and their length); `contig_n50_bp`, `contig_l50`, `contig_n90_bp`, `contig_max_bp`, `contig_median_bp`; `gc_pct` (weighted by length); `coverage_median` (shotgun). Where the contigs went, as a number, length and share of the total length: in a final cluster (`in_final_clusters_contigs`, `_bp`, `_pct`), in a MAG (`in_mags_`...), in an MGE cluster (`in_mges_contigs`, `_bp`), in no genome (`unbinned_contigs`, `_bp`). |
| `contacts` | The contacts of the reads by where they land, from metacooler's `qc` block: `contacts`; `qc_source` (`pairs`: counted on the reads, `pixels`: on the map); `qc_distance_threshold_bp`; within a contig, closer than that (`intra_contig_short`) or not (`intra_contig_long`); between two contigs of a genome (`inter_contig_intra_genome`); between two genomes (`inter_genome`); `trans_contig_ratio`, the share between two contigs (the "3D ratio"); `intra_genome_fraction`, the share within a genome. The network: `network_contigs`, `network_edges`, `normalization`. |
| `clusters` | `core_clusters`, `overlapping_clusters`, `recursive_steps`, `recursive_clusters`, `final_clusters`, `final_clusters_bp`. The final clusters by quality class, each a number and a length (`_bp`): `hq_mags` (completeness ≥ 90%, contamination ≤ 5%), `mq_mags` (70%, 10%), `lq_mags` (50%, 10%), `contaminated_mags` (completeness ≥ 50%, contamination > 10%), `other_clusters` (completeness < 50%). `mags`: the genomes of the map but the MGE genomes. |
| `mges` | `mge_clusters`, `mge_contigs`, `mge_bp` (the MGE clusters, their contigs and length); `mge_genomes` (`--reassign`); `mags_losing_contigs`, the final clusters that held contigs of an MGE genome; `mge_random_clusters` (`--random`). |
| `hosts` | `min_share_pct`, `min_contigs` (the thresholds of `metator host`); `mge_clusters_with_host`, `mge_clusters_with_one_host`, `mge_clusters_with_several_hosts`, `mge_clusters_without_host`; `mags_hosting_mges`. |
| `commands` | The command line of each step, by step. |

| Table | Columns |
|---|---|
| `rounds` | One row per clustering round of the map, in order: `round`, `step` (the command that stored it), `clusters`, `contigs`, `length_bp`. |
| `mags` | One row per MAG, by quality class, then longest first: `name`, `quality` (`HQ`, `MQ`, `LQ`, `Contaminated`, `Other`), `length_bp`, `contigs`, `completeness_pct`, `contamination_pct` (weighted completeness, and weighted redundancy minus 1, in percent), `gc_pct`, `n50_bp`, `hic_abundance`, `sg_abundance`, `intra_genome_ratio` (the share of its contacts within it), `n_mges`, `mges` (comma-separated). |
| `mge_clusters` | One row per MGE cluster: `name`, `contigs`, `length_bp`, `gc_pct`, `n_hosts`, `hosts`, `host_shares_pct` (the share of its signal with each host; lists comma-separated). |
| `hosts` | One row per MGE cluster and host: `mge`, `mag`, `share_pct`. |

```text
# metaTOR report of sample
# format 1: [section] then a key and its value per line; [table: name] then a header and a row
# per line, on spaces; NA: missing. See https://koszullab.github.io/metator/formats/#the-report-of-a-run

[report]
format_version   1
sample           sample
...

[clusters]
core_clusters         5
...
hq_mags               4
hq_mags_bp            468000
mq_mags               1
...

[table: mags]
name                 quality  length_bp  contigs  completeness_pct  contamination_pct  gc_pct  n50_bp  ...
metator_00001_00000  HQ       140000     11       100.00            0.00               49.87   22978   ...
```

### `<name>.report_mqc.yaml` { #report_mqcyaml }

[MultiQC](https://multiqc.info) custom content (`id: metator`): one row of
its General Statistics table, with the numbers that compare across runs,
named as in `<name>.report.txt`: contig length and N50,
`trans_contig_ratio`, `in_mags_pct`, the HQ, MQ, LQ and contaminated MAGs,
and, after `metator mge` and `metator host`, the MGE clusters and those with
a host. `multiqc` finds the files ending in `_mqc.yaml` and gathers the runs
in one table.

```yaml
id: metator
section_name: metaTOR
plot_type: generalstats
pconfig:
- contigs_bp:
    title: Assembly
    ...
data:
  sample:
    contigs_bp: 578000
    contig_n50_bp: 13074
    trans_contig_ratio: 0.526217
    in_mags_pct: 96.88581314878893
    hq_mags: 4
    mq_mags: 1
    ...
```

### `<name>.report.html` { #reporthtml }

A page to open in a browser: the numbers and tables of the run, with links
to the FASTA and cool files of the MAGs, the tables, the plots and the other
two report files, when they are beside it (as `export` and `pipeline` write
them). Its plots are drawn by [plotly.js](https://plotly.com/javascript/),
loaded from its CDN (a pinned version): offline, the page shows only its
numbers and tables.

## Logs

Each command writes its log, `<name>.metator_<date>.log` (date as
`%Y%m%d%H%M%S`, e.g. `sample.metator_20261123225838.log`), in its output
directory; `mge` writes it in the current directory. `<name>` is that of the
[report](#the-report-of-a-run). The log starts with the version of metaTOR,
the date and the command line, then follows what the run did.
