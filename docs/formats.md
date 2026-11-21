# File formats

The files metaTOR reads and writes. Every table is tab-separated.

## Inputs

### Assembly

A FASTA file of the contigs, compressed with gzip or not (`--assembly`): the
assembly that the contact map was made from.

### Depth

The output of `jgi_summarize_bam_contig_depths` of MetaBAT2 (`--depth`), a
header line, then one line per contig:

```text
contigName	contigLen	totalAvgDepth	sample.bam	sample.bam-var
NODE_522	22786	4.76595	4.76595	15.5294
```

metaTOR reads the length and the coverage of each contig, by name, from the
first three columns. It must hold every contig of the assembly, with its
length, as `jgi_summarize_bam_contig_depths` writes it from reads aligned on
that assembly: a run with the depth file of another assembly stops with an
error.

### Pairs

The pairs files that `metator mge` reads, those that hicstuff writes with the
map (`<prefix>.valid_idx.pairs.gz`), follow the
[pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md)
of the 4D Nucleome, with 7 columns: read name, contig and position of the
first mate, contig and position of the second, and the strand of each.

```text
## pairs format v1.0
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#chromsize: NODE_522 22786
...
H9:1:HGJYYBBXY:8:1102:3792:13693	NODE_100063	161	NODE_100063	402	-	-
```

hicstuff writes them sorted by position and compressed with bgzip, with two
more columns (the restriction fragment of each mate); `pairix -p pairs` indexes
them (`.px2`), which `metator mge` needs to read a compressed pairs file.
`metacooler create pairs` makes a contact map from them, or from the pairs
files of another aligner.

### Lists of contigs

The MGE contigs of `metator mge --mges`: one contig name per line.

The clusters of `metator mge --binning`: a contig and its cluster per line,
without header, as `binning.txt` or the collection exported by anvi'o.

## The contig data tables

`contig_data_partition.txt` and `contig_data_final.txt` hold one row per
contig of the map, in its order, with a header. Each step adds its columns to
those of the step before; the first come from the network, which
`metator network` stores in the contigs table of the map:

| Column | Written by | What it holds |
|---|---|---|
| `ID` | network | The ID of the contig, from 1: its position in the map, the ID of its node in the network |
| `Name` | network | Its name |
| `Size` | network | Its length, in bp |
| `GC_content` | network | Its GC content, in percent |
| `Hit` | network | Its hits: the ends of the contacts of the map that lie on it (a contact within the contig counts twice) |
| `Shotgun_coverage` | network | Its coverage in the depth file, or in the `Shotgun_coverage` column of the map, when it gives every contig a value; `-` otherwise |
| `Restriction_site` | network | Its number of restriction sites plus one (its restriction fragments), from `--enzyme`, or in the `Restriction_site` column of the map, when it gives every contig a value; `-` otherwise |
| `Core_cluster_ID`, `Core_cluster_contigs`, `Core_cluster_size` | partition | Its core cluster: the ID, number of contigs and size in bp of the cluster |
| `Overlapping_cluster_ID`, `Overlapping_cluster_contigs`, `Overlapping_cluster_size` | partition | Its overlapping cluster, likewise |
| `Recursive_cluster_ID`, `Recursive_cluster_contigs`, `Recursive_cluster_size` | validation | Its recursive cluster, likewise: `00000`, `-`, `-` when its cluster was not partitioned again |
| `Final_cluster` | validation | Its final cluster, or `ND` (not determined) |

```text
	ID	Name	Size	GC_content	Hit	Shotgun_coverage	Restriction_site	Core_cluster_ID	Core_cluster_contigs	Core_cluster_size	Overlapping_cluster_ID	Overlapping_cluster_contigs	Overlapping_cluster_size	Recursive_cluster_ID	Recursive_cluster_contigs	Recursive_cluster_size	Final_cluster
0	1	NODE_522	22786	62.0775915035548	1683	4.76595	162	00001	2	23301	00001	15	49021	00000	-	-	metator_00001_00000
1	2	NODE_1404	15015	64.2823842823842	952	4.00585	118	00002	1	15015	00002	5	19711	00000	-	-	ND
```

`contig_data_final.txt` starts with a column without header, the row number
from 0. A contig with no edge in the network (no hits, or no contact with
another contig) has `-` in the columns of its core and overlapping clusters. As
inputs, the tables are read by the headers of their columns: other columns
may be added, in any order, such as the bins of another binning.

## The network

The network of contigs is in the map, which `metator network` stores it in
(see [Building the network](guides/network.md#building-the-network)): metaTOR
2.0 writes no `network.txt`. metacooler exports it, one line per pair of
contigs in contact, by their names, with the normalised weight of their edge:

```console
metacooler network sample.metacool --level contig --balance metator_empirical_hit
```

## The clusters

### Cluster summary

`cluster_summary.txt` holds one row per final cluster, named in its first
column. Its first columns are those of miComplete:

| Column | What it holds |
|---|---|
| `Length` | The size of the cluster, in bp |
| `GC-content` | Its GC content, in percent |
| `Present Markers` | The single-copy marker genes found in it |
| `Completeness`, `Redundancy` | The share of the markers found, and their mean number of copies (1 for no duplicate) |
| `Weighted completeness`, `Weighted redundancy` | The same, with miComplete's weights of the markers of the set |
| `Contigs`, `N50`, `L50`, `N90`, `L90` | Its number of contigs, and the N50, L50, N90 and L90 of their sizes |
| `CDs` | Its number of coding sequences, found by prodigal |
| `recursive` | Whether the cluster comes from the recursive partition |
| `step` | The recursive step that made the cluster, 0 for an overlapping cluster |
| `parent` | The cluster it was split from (an overlapping cluster at step 1, a recursive cluster after), empty otherwise |
| `HiC_abundance` | The share of the hits of all the contigs that are on the cluster, in percent |
| `SG_abundance` | The share of the shotgun reads (size × coverage) of all the contigs that are on the cluster, in percent; with `--depth` only |

```text
	Length	GC-content	Present Markers	Completeness	Redundancy	Weighted completeness	Weighted redundancy	Contigs	N50	L50	N90	L90	CDs	recursive	step	parent	HiC_abundance	SG_abundance
metator_00001_00000	49021.00	62.03	1.00	0.01	1.00	0.02	1.00	15.00	9019.00	2.00	1390.00	9.00	70.00	True	0		2.4736	2.6314
```

Each marker set, bacterial (105 markers) and archaeal (131), is tried on
every cluster, which keeps the one that gives it the higher weighted
completeness. The quality classes of the log take the weighted redundancy
minus one as the contamination of a cluster; the recursive
partition divides it by the weighted completeness (see
[Cluster quality](concepts.md#cluster-quality)).

### binning.txt

Each contig of a final cluster and its cluster, without header: the collection
that anvi'o imports (see [Curating MAGs with anvi'o](guides/curation.md)).

```text
NODE_522	metator_00001_00000
NODE_2487	metator_00007_00000
```

### miComplete results

`overlapping_micomplete_results.txt` and `recursive_micomplete_step_<n>.txt`
hold miComplete's results on the overlapping clusters and on the recursive
clusters of each step: three comment lines (`##`), then the columns of the
cluster summary up to `CDs`, and the marker set kept, `Bacteria` or `Archaea`.

### Clustering rounds

`clusters.txt` holds every cluster of every round of the run, in the order
they were made: `overlapping` (the overlapping clusters that miComplete
checks), `recursive_<k>` (those that step `k` of the validation makes by
partitioning a contaminated cluster again) and `final`. Tab-separated, with a
header: the round, the name, the round and name of its parent (the cluster it
was split from, or for a final cluster the cluster it is; `NA` for none), then
miComplete's columns, `NA` where the cluster has none.

```text
round	name	parent_round	parent	Length	GC-content	...	Completeness	...
overlapping	metator_00001_00000	NA	NA	2817353.0	52.30	...	0.9655	...
recursive_1	metator_00002_00001	overlapping	metator_00002_00000	1504321.0	41.80	...	0.7103	...
final	metator_00001_00000	overlapping	metator_00001_00000	2817353.0	52.30	...	0.9655	...
```

`cluster_members.txt` holds the contigs of each: one row per contig of a
cluster, with its round, cluster and contig, under a header. The map of the
clusters stores both (see
[Clustering rounds](guides/metacool.md#clustering-rounds)).

### FASTA files of the clusters

`overlapping_clusters/`, `recursive_clusters/step_<n>/` (the recursive
clusters of step `n`) and `final_clusters/` hold one FASTA file per cluster,
`<cluster>.fa`, with the records of its contigs as they are in the assembly.
The final clusters also hold the pyfastx index of each, `<cluster>.fa.fxi`,
which the validation leaves there.

## The contact map of the clusters

The map that a run starts from, and writes the map of its clusters over, is
a `.metacool` file of [metacooler](https://github.com/js2264/metacooler): a
cooler whose contigs are grouped into genomes, one per final cluster in the
map of the clusters, which keeps the bins of the map. See
[The contact map of the clusters](guides/metacool.md).

## MGE cluster summary

`mges_cluster_summary.tsv` holds one row per MGE cluster, with a header:

| Column | What it holds |
|---|---|
| `BinName` | `MetaTOR_MGE_00001` and on |
| `BinLength`, `GC`, `Hit` | The size, GC content and hits of the cluster |
| `BinningScore` | The lowest score at which its contigs were merged, `NA` for a single contig |
| `ContigsNumber`, `Contigs` | Its number of contigs, and their names, separated by commas |
| `MetagenomicBin`, `AssociationScore`, `AssociatedBins` | For its host MAGs: `None`, `NA` and `NA`, as `metator mge` does not look for them |

```text
BinName	BinLength	GC	Hit	BinningScore	ContigsNumber	Contigs	MetagenomicBin	AssociationScore	AssociatedBins
MetaTOR_MGE_00064	30610.00	68.93	18383.00	1.20	2.00	NODE_8199_len_17842,NODE_10116_len_12768	None	NA	NA
MetaTOR_MGE_00105	23587.00	37.36	2662.00	NA	1.00	NODE_10949_len_23587	None	NA	NA
```

With `--map`, the map of the run stores these rows as the round `mge` of its
clustering rounds (see
[Clustering mobile genetic elements](guides/mge.md#in-the-contact-map-of-the-run)):
one cluster per row, named `BinName`, with the contigs of `Contigs` and the
other columns as written, but `BinLength` and `ContigsNumber`, which are the
`length` and `n_contigs` that metacooler computes from the contigs.

## Logs

Each command writes its log, `metator_<command>_<date>.log`: the version of
metaTOR, the date and the command line, then what the run did.
