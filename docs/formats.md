# File formats

The files metaTOR reads and writes. Every table is tab-separated.

## Inputs

### Assembly

A FASTA file of the contigs, compressed with gzip or not (`--assembly`): the
assembly that the contact map was made from. A contig is read by its name,
the first word of its header: an assembly that names a contig by several
records is refused. The restriction sites of `--enzyme` are counted whatever
the case of the bases: those of a soft-masked assembly, whose repeats are in
lower case, too.

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
error. `jgi_summarize_bam_contig_depths` writes the length of a contig of a
million bases or more in scientific notation, to 6 significant digits
(`1.21754e+06` for 1,217,543 bases): metaTOR reads it as the length of the
contig rounded so.

### Pairs

The pairs files that the map is made from, those that hicstuff writes with
the map (`<prefix>.valid_idx.pairs.gz`), follow the
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
more columns (the restriction fragment of each mate). `metacooler create
pairs` makes a contact map from them, or from the pairs files of another
aligner, with the contact statistics of the reads of each contig that
`metator mge` reads (see
[Clustering mobile genetic elements](guides/mge.md#running-it)); and so do
`metator network` and `metator pipeline`, given them as their arguments, one
bin per contig (see [The inputs of a run](cli.md#inputs)). They read the
contig and the position of each mate (columns 2 to 5), in any order, and the
`#chromsize` lines of the header, which must give the contigs of the assembly
their lengths.

From BAM files, a run writes the pairs of each library in its temporary
directory, `library_<i>.pairs`, as they come: the read name `.` (the map does
not keep it), each mate at the 5' end of its alignment (1-based: its first
aligned base on the forward strand, its last on the reverse strand), the
contig that comes first in the assembly first, and no `#chromsize` line (the
contigs of the map are those of the assembly):

```text
## pairs format v1.0
#shape: upper triangle
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2
.	NODE_100063	161	NODE_100063	402	-	-
```

### Lists of contigs

The MGE contigs of `metator mge --mges`: one contig name per line (the first
word of each line).

## The contig data tables

`contig_data_final.txt` holds one row per contig of the map, in its order,
with a header. [`metator dump`](cli.md#dump) writes it from the map of the
clusters of a run; the validation writes its own in its temporary directory,
which `-N` keeps. Each step adds its columns to those of the step before; the
first come from the network, which `metator network` stores in the contigs
table of the map, then the core and overlapping clusters, which
`metator partition` stores as the clustering rounds of the map (the
validation writes the table of the partition, `contig_data_partition.txt`, in
its temporary directory):

| Column | Written by | What it holds |
|---|---|---|
| `ID` | network | The ID of the contig, from 1: its position in the map as the run read it, the ID of its node in the network. A run that regroups the map (`validation`, `pipeline`) moves its contigs: join a table to the network of another map by contig name (`metator.network.relabel_network`), not by ID. The tables of the validation only: `metator dump` leaves it out |
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

The validation's `contig_data_final.txt` starts with a column without
header, the row number from 0, then `ID`; that of `metator dump` starts with
`Name`, its rows in the order of the map of the clusters, which the
validation regroups by cluster. A contig with no edge in the network (no
hits, or no contact with another contig) has `-` in the columns of its core
and overlapping clusters. `Core_cluster_ID` is a number when every contig is
in a core cluster, and has five digits, as the other IDs, next to the `-` of
a contig in none. As inputs, the tables are read by the headers of their
columns: other columns may be added, in any order, such as the bins of
another binning.

## The network

The network of contigs is in the map, which `metator network` stores it in
(see [Building the network](guides/network.md#building-the-network)): the
steps write no `network.txt`. [`metator dump`](cli.md#dump) writes it from
the map of the clusters of a run: one line per pair of contigs in contact,
without a header, with three columns, `contig1` and `contig2`, the names of
the two contigs, and the weight of their edge (the normalization that
`metator network` stored; the count of their contacts for `None`). The lines
come in the order of the map: each contig, from the first, with the contigs
after it. The contacts within a contig, and the edges of weight 0, are no
edge of the network. 1.x wrote the IDs of the contigs, in the order of their
names; 2.0 writes their names.

metacooler exports the network of any map that holds it, the map of the
contigs once `metator network` ran included. It also writes a line per
contig for the contacts within it (`node1` equal to `node2`), which are no
edge of the network, and would write the edges of weight 0, which metaTOR
leaves out: drop them, as `awk` does here.

```console
metacooler network sample.metacool --level contig --balance metator_empirical_hit --header \
    | awk -F '\t' '$1 != $2 && $3 != 0' > network.tsv
```

`network.tsv` has three columns: `node1` and `node2`, the names of the two
contigs of an edge, and `value`, its weight (the raw count of their contacts
without `--balance`, for the normalisation `None`). Without `--header`,
metacooler writes no header line.

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
| `CDs` | Its number of coding sequences, found by Prodigal (pyrodigal) |
| `recursive` | Whether the cluster comes from the recursive partition |
| `step` | The recursive step that made the cluster, 0 for an overlapping cluster |
| `parent` | The cluster it was split from (an overlapping cluster at step 1, a recursive cluster after), empty otherwise |
| `HiC_abundance` | The share of the hits of all the contigs that are on the cluster, in percent |
| `SG_abundance` | The share of the shotgun reads (size × coverage) of all the contigs of the map that are on the cluster, in percent; with `--depth` only. The clusters of a map made with a contig filter add up to 100%, where 1.x counted every contig of the assembly |

```text
	Length	GC-content	Present Markers	Completeness	Redundancy	Weighted completeness	Weighted redundancy	Contigs	N50	L50	N90	L90	CDs	recursive	step	parent	HiC_abundance	SG_abundance
metator_00001_00000	49021.00	62.03	1.00	0.01	1.00	0.02	1.00	15.00	9019.00	2.00	1390.00	9.00	70.00	False	0		2.4736	2.6314
```

Each marker set, bacterial (105 markers) and archaeal (131), is tried on
every cluster, which keeps the one that gives it the higher weighted
completeness. The quality classes of the log take the weighted redundancy
minus one as the contamination of a cluster; the recursive
partition divides it by the weighted completeness (see
[Cluster quality](concepts.md#cluster-quality)).

miComplete's numbers have two decimals, whole ones too (`Length`), the
abundances four, and `step` none. The map of the clusters holds the same
values, in its genomes table,
and records how the table is written, under the `validation` entry of the
`metator` key of its metadata: `summary_rows`, the clusters in the order of
the lines; `summary_columns`, the columns in order; and `summary_formats`,
the format of each column as printf writes it (`%.2f`, `%.4f`, `%d`, or
`%s` for text and booleans), a missing value being an empty field. From
them, the map alone writes the table again, byte for byte, as
[`metator dump`](cli.md#dump) does; after `metator host`, with the number of
MGE clusters that each cluster hosts and their names last (`n_mges`, `mges`:
see [Hosts of the MGE clusters](#hosts-of-the-mge-clusters)).

`metator dump` writes the genomes of the map: after
[`metator mge --reassign`](guides/mge.md#mge-genomes), the final clusters,
without those whose contigs all went to MGE genomes, then the MGE genomes,
with four more columns:

| Column | What it holds |
|---|---|
| `GC` | The GC content of an MGE genome, in percent, empty for a final cluster |
| `Hit` | Its hits, the sum of those of its contigs |
| `BinningScore` | Its clustering score: the lowest score at which its contigs were merged, empty for a single contig |
| `mge` | `True` for an MGE genome, `False` for a final cluster, as in the column `mge` of the genomes table |

The final clusters keep the values of their validation, which describe them
with the MGE contigs that they lost.

### binning.txt

Each contig of a genome of the map and its genome, without header: the
collection that anvi'o imports (see
[Curating MAGs with anvi'o](guides/curation.md)). Its genomes are the final
clusters, and the MGE genomes after `metator mge --reassign`.

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
they were made: `core` and `overlapping` (every core and overlapping cluster
of the partition, those that miComplete checks with its columns),
`recursive_<k>` (those that step `k` of the validation makes by partitioning
a contaminated cluster again) and `final`; that of `metator dump` holds the
round `mge` that `metator mge` stores too (and `mge_random`, when
`metator mge --random` stores one). Tab-separated, with a
header: the round, the name, the round and name of its parent (the cluster it
was split from, or for a final cluster the cluster it is; `NA` for none), its
number of contigs and their length (`n_contigs`, `length`, as the map
computes them from its contigs), then miComplete's columns, `NA` where the
cluster has none (and those of an MGE cluster: see
[MGE clusters](#mge-clusters)).

```text
round	name	parent_round	parent	n_contigs	length	Length	GC-content	...	Completeness	...
overlapping	metator_00001_00000	NA	NA	31	2817353	2817353.0	52.30	...	0.9655	...
recursive_1	metator_00002_00001	overlapping	metator_00002_00000	18	1504321	1504321.0	41.80	...	0.7103	...
final	metator_00001_00000	overlapping	metator_00001_00000	31	2817353	2817353.0	52.30	...	0.9655	...
```

`cluster_members.txt` holds the contigs of each: one row per contig of a
cluster, with its round, cluster and contig, under a header. The map of the
clusters stores both (see
[Clustering rounds](guides/metacool.md#clustering-rounds)).

### FASTA files of the clusters

`fasta/mags_with-mges/`, in the output directory, holds one FASTA file per
final cluster, `<cluster>.fa`, with the records of its contigs as they are in
the assembly, the MGE contigs that the validation clustered with it included.
`fasta/mags/` (`pipeline`) holds each MAG, a genome of the map of the run but
the MGE genomes, twice: `<MAG>.fa`, the records of its contigs as they are in
the assembly, in the order of the map, and `<MAG>.stitched.fa`, one record
named after it, its contigs, the longest first, joined by 200 `N`, in lines of
80 bases. `fasta/mges/` (`mge`, `pipeline --mges`) holds each MGE cluster the
same two ways, its contigs the longest first, ties by name. metacooler writes the files of the clusters that the validation
checks (`metacooler.write_fasta`), in its temporary directory (`-N` keeps
it): the overlapping clusters in `overlapping_clusters/`, in the order of
the partition (by core cluster, then by ID), and the recursive clusters of
step `n` in `recursive_clusters/step_<n>/`. A final cluster is a copy of the
file that miComplete checked. A gzip-compressed assembly is decompressed
once per run, in the temporary directory, and every file is written from
that copy. No pyfastx index is made, beside the assembly or the clusters.

### Cool files of the MAGs

`cool/` (`pipeline`) holds one [cooler](https://github.com/open2c/cooler) file
per MAG, `<MAG>.cool`, which `metacooler extract` writes from the map of the
run: the contacts of its contigs, and with `--mges` of those of the MGE
clusters that it hosts, at the bins of the map, each contig a chromosome of
the file. `cooler balance` balances each file, the whole of it at once, with
its defaults: the `weight` column of its bins, which cooler, cooltools and
HiGlass read (`clr.matrix(balance=True)`); a file whose every bin is masked,
too small for `cooler balance`, keeps weights that are all `NaN`, with a
warning in the log.

## The contact map of the clusters

The map that a run starts from, and writes the map of its clusters over, is
a `.metacool` file of [metacooler](https://github.com/js2264/metacooler): a
cooler whose contigs are grouped into genomes, one per final cluster in the
map of the clusters, which keeps the bins of the map. See
[The contact map of the clusters](guides/metacool.md).

## MGE clusters

`metator mge` stores its MGE clusters in the map, as the round `mge` of its
clustering rounds (`mge_random` for `--random`), one cluster per MGE cluster
(see
[Clustering mobile genetic elements](guides/mge.md#in-the-contact-map-of-the-run)),
and so does `metator pipeline --mges`, in its step `mge`:

| Column | What it holds |
|---|---|
| `name` | `MetaTOR_MGE_00001` and on, the name of its FASTA files in `fasta/mges/` |
| `length`, `n_contigs` | The size of the cluster and its number of contigs, which metacooler computes from its contigs |
| `GC`, `Hit` | The GC content and the hits of the cluster |
| `BinningScore` | The lowest score at which its contigs were merged, missing for a single contig |
| `n_hosts`, `hosts`, `host_shares` | After `metator host`: the number of hosts of the MGE cluster, their names by decreasing share, and their shares in percent, with two decimals, each list separated by commas (see [Hosts of the MGE clusters](#hosts-of-the-mge-clusters)) |

Its contigs are the members of the round (`mc.cluster_members("mge")`), in
the order of the map, as metacooler stores the members of every round:
`fasta/mges/` gives them the longest first, ties by name. `metator dump`
writes the round with the others, in `clusters.txt` (with `n_contigs` and
`length`, the `ContigsNumber` and `BinLength` of 1.x's
`mges_bin_summary.tsv`) and `cluster_members.txt`. `metator mge` looks for
no host: `metator host` does, and adds its columns to the round.

## Hosts of the MGE clusters

`metator host` writes `tables/mge_hosts.txt` in its output directory, as
`metator pipeline --mges` does in its step `host`: one row per
MGE cluster and MAG in contact, the MGE clusters in the order of the round
`mge`, the MAGs of each by decreasing signal (see
[Clustering mobile genetic elements](guides/mge.md#the-hosts)).
Tab-separated, with a header, every number as it is:

| Column | What it holds |
|---|---|
| `mge` | The MGE cluster (`MetaTOR_MGE_00001`...) |
| `mag` | The MAG: a genome of the map |
| `signal` | The sum of the weights of the edges of the network between their contigs, with the normalization that `metator network` stored (the count of their contacts for `None`) |
| `contigs` | The number of contigs of the MAG in contact with the MGE cluster |
| `share` | The signal, in percent of the signal of the MGE cluster with all the MAGs |
| `share_of_all` | The signal, in percent of the signal of the MGE cluster with every other contig: those of the MAGs, those in no MAG, and those of the other MGE clusters |
| `host` | Whether the MAG hosts the MGE cluster (`True`, `False`): a `share` of at least `--min-share` and at least `--min-contigs` contigs |

```text
mge	mag	signal	contigs	share	share_of_all	host
MetaTOR_MGE_00013	metator_00005_00000	0.7090148667776763	37	93.8246167981684	76.06906868895449	True
MetaTOR_MGE_00013	metator_00020_00000	0.005037475316584434	4	0.6666139362589284	0.5404626529450921	False
```

`plots/mge_hosts.pdf` plots the distribution of the signal of the hosts, on a log
scale, and of their share. The map holds the hosts of each MGE cluster in
its round `mge` (see [MGE clusters](#mge-clusters)), and the MGE clusters of
each MAG in its genomes table:

| Column | What it holds |
|---|---|
| `n_mges` | The number of MGE clusters that the MAG hosts, missing for a genome that is no MAG |
| `mges` | Their names, in the order of the round `mge`, separated by commas |

`metator dump` writes them last in `cluster_summary.txt`, after the columns
of the validation (and those of the MGE genomes, after `metator mge
--reassign`, whose `n_mges` and `mges` are empty).

## The report of a run

`metator report` writes the report of a run in its output directory, from the
map of its clusters alone, and so do `metator export` and `metator pipeline`,
in its step `export` (see [`report`](cli.md#report)): `report.txt`,
`report_mqc.yaml` and `report.html`. The run is named after its map, its
`sample`: `sample` for `sample.metacool`, `NAME` for the entry
`file.mmetacool::entries/NAME` of a collection. Give the maps of the runs of a
project the names of their samples (`-p` of `pipeline` from pairs or BAM
files, `-P` of `hicstuff pipeline`): MultiQC gathers the reports of the runs
by these names.

### report.txt

A text file, in sections, then tables, which a script, or a MultiQC module,
parses line by line (`metator.report.read_txt` reads it back):

- a line that starts with `#` is a comment; the first one,
  `# metaTOR report of <sample>`, tells the file;
- `[name]` starts a section: each of its lines is a key, spaces, and its
  value, the rest of the line;
- `[table: name]` starts a table: its first line is its header, then a row
  per line, each field of a row in a column of the header, separated by
  spaces (no field holds a space);
- `NA` is a missing value; an integer is written as such; a value in percent
  (a key or a column that ends in `_pct`) with 2 decimals, a ratio or a
  fraction (`_ratio`, `_fraction`), from 0 to 1, with 4; a length in bp ends
  in `_bp`; other numbers are written in full.

The keys, the columns and their units are those of the version of the
format, the `format_version` of the section `report`, 1 here: a key, a column
or a unit that changes, or goes, changes it, and one that comes does not, so
a parser skips the keys and the columns that it does not know. A section, a
key or a table that does not apply to a run is left out: the sections `mges`
and `hosts` and their tables before `metator mge` and `metator host`.

| Section | Keys |
|---|---|
| `report` | `format_version`; `sample`, the name of the run; `map`, its map; `metator_version`, the version of metaTOR that wrote the report; `written`, when; `map_written`, when the map was written. |
| `run` | `steps`, the steps that the map records, in order, separated by commas (`network`, `partition`, `validation`, `mge`, `host`); `input`, the input of the run (its map, or its pairs or BAM files, separated by commas), `assembly`, `enzyme` and `depth`, as the run recorded them (`NA`: the map gave them); `normalization` of the network; `algorithm`, `iterations`, `overlap_pct` and `res_param` of the partition; `size_bp`, `rec_iter` and `rec_overlap_pct` of the validation; `seed` and `rec_seed`, the `--seed` of the partition and of the validation (both that of `pipeline --seed`; `NA`: no seed); `mge_threshold_bin` of `mge`; `host_min_share_pct` and `host_min_contigs` of `host`. |
| `contigs` | `contigs` and `contigs_bp`, the contigs of the map and their length; `contig_n50_bp`, `contig_l50`, `contig_n90_bp`, `contig_max_bp`, `contig_median_bp`; `gc_pct`, their GC content, weighted by their length; `coverage_median`, their median shotgun coverage; where they went, their number, length and share of the length of the contigs: in a final cluster (`in_final_clusters_contigs`, `_bp`, `_pct`), in a MAG (`in_mags_`...), in an MGE cluster (`in_mges_contigs`, `_bp`), in no genome of the map (`unbinned_contigs`, `_bp`). |
| `contacts` | The contacts of the reads by where they land, as metacooler counted them, with the genomes of the map (its `qc` block): `contacts`; `qc_source` (`pairs`: counted on the reads, `pixels`: on the map); `qc_distance_threshold_bp`; within a contig, less than that apart (`intra_contig_short`) or more (`intra_contig_long`); between two contigs of a genome (`inter_contig_intra_genome`); between two genomes (`inter_genome`); `trans_contig_ratio`, the share of the contacts between two contigs, the "3D ratio"; `intra_genome_fraction`, the share within a genome. Then the network of the contigs: `network_contigs`, `network_edges` and its `normalization`. |
| `clusters` | The clusters of the rounds: `core_clusters`, `overlapping_clusters`, `recursive_steps`, `recursive_clusters`, `final_clusters` and `final_clusters_bp`. The final clusters by quality class, each its number and length (`_bp`): `hq_mags` (completeness of at least 90% and contamination of at most 5%), `mq_mags` (70%, 10%), `lq_mags` (50%, 10%), `contaminated_mags` (completeness of at least 50%, contamination above 10%) and `other_clusters` (completeness below 50%), classes of the values of miComplete, as the log of the run counts them. `mags`, the MAGs: the genomes of the map but the MGE genomes. |
| `mges` | `mge_clusters`, `mge_contigs` and `mge_bp`, the MGE clusters of `metator mge`, their contigs and length; `mge_genomes`, the MGE genomes of the map (`--reassign`); `mags_losing_contigs`, the final clusters that held contigs of an MGE genome; `mge_random_clusters`, those of `--random`. |
| `hosts` | `min_share_pct` and `min_contigs`, the thresholds of `metator host`; `mge_clusters_with_host`, `mge_clusters_with_one_host`, `mge_clusters_with_several_hosts`, `mge_clusters_without_host`; `mags_hosting_mges`. |
| `commands` | The command line of each step, by step. |

| Table | Columns |
|---|---|
| `rounds` | One row per clustering round of the map, in order: `round`, `step` (the command that stored it), `clusters`, `contigs`, `length_bp`. |
| `mags` | One row per MAG, by quality class, then the longest first: `name`, `quality` (`HQ`, `MQ`, `LQ`, `Contaminated`, `Other`), `length_bp`, `contigs`, `completeness_pct`, `contamination_pct` (miComplete's weighted completeness, and its weighted redundancy minus 1, in percent), `gc_pct`, `n50_bp`, `hic_abundance`, `sg_abundance`, `intra_genome_ratio` (the share of the contacts of the MAG within it), `n_mges` and `mges`, the MGE clusters that it hosts, separated by commas. |
| `mge_clusters` | One row per MGE cluster: `name`, `contigs`, `length_bp`, `gc_pct`, `n_hosts`, `hosts` and `host_shares_pct`, its hosts and the share of its signal with each, separated by commas. |
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

### report_mqc.yaml

The custom content of [MultiQC](https://multiqc.info): a row of its General
Statistics table for the run, `id: metator`, with the numbers of the run that
compare across runs: the length of its contigs and their N50, its 3D ratio
(`trans_contig_ratio`), the share of the contigs in a MAG (`in_mags_pct`), its
HQ, MQ, LQ and contaminated MAGs, and, after `metator mge` and `metator host`,
its MGE clusters and those with a host, named as in `report.txt`. MultiQC
finds the files whose names end in `_mqc.yaml`: `multiqc` on a directory
that holds the output directories of the runs gathers their reports in one
table, a row per run.

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

### report.html

A page to read the run with, which a browser opens from the output directory:
its numbers and tables are in the page, and its links to the FASTA files and
the cool files of the MAGs, the tables and the plots of the run, which it
links to when they are beside it, as `export` and `pipeline` write them. Its
interactive plots are drawn by [plotly.js](https://plotly.com/javascript/),
which the page loads from its CDN when it is opened, a pinned version: without
an internet connection, the page shows its numbers and tables, and says that
its plots need it.

## Logs

Each command writes its log, `metator_<command>_<date>.log`: the version of
metaTOR, the date and the command line, then what the run did.
