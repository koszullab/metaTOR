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
map from them (see [The inputs of a run](cli.md#inputs) and
[Clustering mobile genetic elements](guides/mge.md#running-it)). They read
columns 2 to 5, in any order, and the `#chromsize` lines, which must give
the contigs of the assembly their lengths.

From BAM files, a run writes the pairs of each library in its temporary
directory, `library_<i>.pairs`: read name `.`, each mate at the 5' end of
its alignment (1-based), the contig that comes first in the assembly first,
and no `#chromsize` line.

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

- Both marker sets, bacterial (105 markers) and archaeal (131), are tried on
  every cluster, which keeps the one with the higher weighted completeness.
- Contamination is the weighted redundancy minus one; the recursive
  partition divides it by the weighted completeness (see
  [Cluster quality](concepts.md#cluster-quality)).
- miComplete's numbers have two decimals (whole ones too), the abundances
  four, and `step` none.
- The genomes table of the map holds the same values. The `validation`
  entry of the map's `metator` metadata records how to write the table:
  `summary_rows` (the clusters in order), `summary_columns` (the columns in
  order) and `summary_formats` (each column's printf format: `%.2f`,
  `%.4f`, `%d`, or `%s`; a missing value is an empty field). From them,
  [`metator dump`](cli.md#dump) writes the table again, byte for byte.
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
file per MAG, `<MAG>.cool`, written by `metacooler extract` from the map:
the contacts of its contigs (and, with `--mges`, of the MGE clusters it
hosts), at the bins of the map, each contig a chromosome. `cooler balance`
balances each file with its defaults, in the `weight` column of its bins
(`clr.matrix(balance=True)`). A file too small for `cooler balance` keeps
`NaN` weights, with a warning in the log.

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
map order. `metator dump` writes the round in `clusters.txt` (where
`n_contigs` and `length` are 1.x's `ContigsNumber` and `BinLength` of
`mges_bin_summary.tsv`) and `cluster_members.txt`.

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

## The report of a run

`metator report` writes the report of a run from the map of its clusters
alone, and so do `metator export` and `metator pipeline` (see
[`report`](cli.md#report)): `<name>.report.txt`, `<name>.report_mqc.yaml`
and `<name>.report.html`, in the output directory.

`<name>` is the name of the run's map without `.metacool`: `sample` for
`sample.metacool`, `NAME` for an entry `file.mmetacool::entries/NAME`. A map
made from pairs or BAM files is `<prefix>.metacool` (`pipeline --prefix`) or
`<assembly>.metacool`. MultiQC names the rows of the runs by these names: name
the maps after their samples (`-p` of `pipeline` from pairs or BAM files,
`-P` of `hicstuff pipeline`).

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
| `report` | `format_version`; `sample`, the name of the run; `map`, its map; `metator_version`, the version of metaTOR that wrote the report; `written`, when; `map_written`, when the map was written. |
| `run` | `steps`, the steps that the map records, in order, comma-separated (`network`, `partition`, `validation`, `mge`, `host`); `input` (the map, or the pairs or BAM files, comma-separated), `assembly`, `enzyme` and `depth`, as the run recorded them (`NA`: the map gave them); `normalization` of the network; `algorithm`, `iterations`, `overlap_pct` and `res_param` of the partition; `size_bp`, `rec_iter` and `rec_overlap_pct` of the validation; `seed` and `rec_seed`, the `--seed` of the partition and of the validation (both from `pipeline --seed`; `NA`: no seed); `mge_threshold_bin` of `mge`; `host_min_share_pct` and `host_min_contigs` of `host`. |
| `contigs` | `contigs` and `contigs_bp`, the contigs of the map and their length; `contig_n50_bp`, `contig_l50`, `contig_n90_bp`, `contig_max_bp`, `contig_median_bp`; `gc_pct`, their GC content, weighted by length; `coverage_median`, their median shotgun coverage. Where they went, as a number, length and share of the length of the contigs: in a final cluster (`in_final_clusters_contigs`, `_bp`, `_pct`), in a MAG (`in_mags_`...), in an MGE cluster (`in_mges_contigs`, `_bp`), in no genome (`unbinned_contigs`, `_bp`). |
| `contacts` | The contacts of the reads by where they land, as metacooler counted them with the genomes of the map (its `qc` block): `contacts`; `qc_source` (`pairs`: counted on the reads, `pixels`: on the map); `qc_distance_threshold_bp`; within a contig, less than that apart (`intra_contig_short`) or more (`intra_contig_long`); between two contigs of a genome (`inter_contig_intra_genome`); between two genomes (`inter_genome`); `trans_contig_ratio`, the share of contacts between two contigs (the "3D ratio"); `intra_genome_fraction`, the share within a genome. Then the network: `network_contigs`, `network_edges` and its `normalization`. |
| `clusters` | `core_clusters`, `overlapping_clusters`, `recursive_steps`, `recursive_clusters`, `final_clusters` and `final_clusters_bp`. The final clusters by quality class (from miComplete's values, as the log counts them), each a number and a length (`_bp`): `hq_mags` (completeness ≥ 90%, contamination ≤ 5%), `mq_mags` (70%, 10%), `lq_mags` (50%, 10%), `contaminated_mags` (completeness ≥ 50%, contamination > 10%) and `other_clusters` (completeness < 50%). `mags`, the MAGs: the genomes of the map but the MGE genomes. |
| `mges` | `mge_clusters`, `mge_contigs` and `mge_bp`, the MGE clusters of `metator mge`, their contigs and length; `mge_genomes`, the MGE genomes of the map (`--reassign`); `mags_losing_contigs`, the final clusters that held contigs of an MGE genome; `mge_random_clusters`, those of `--random`. |
| `hosts` | `min_share_pct` and `min_contigs`, the thresholds of `metator host`; `mge_clusters_with_host`, `mge_clusters_with_one_host`, `mge_clusters_with_several_hosts`, `mge_clusters_without_host`; `mags_hosting_mges`. |
| `commands` | The command line of each step, by step. |

| Table | Columns |
|---|---|
| `rounds` | One row per clustering round of the map, in order: `round`, `step` (the command that stored it), `clusters`, `contigs`, `length_bp`. |
| `mags` | One row per MAG, by quality class, then longest first: `name`, `quality` (`HQ`, `MQ`, `LQ`, `Contaminated`, `Other`), `length_bp`, `contigs`, `completeness_pct`, `contamination_pct` (miComplete's weighted completeness, and its weighted redundancy minus 1, in percent), `gc_pct`, `n50_bp`, `hic_abundance`, `sg_abundance`, `intra_genome_ratio` (the share of the MAG's contacts within it), `n_mges` and `mges`, the MGE clusters it hosts, comma-separated. |
| `mge_clusters` | One row per MGE cluster: `name`, `contigs`, `length_bp`, `gc_pct`, `n_hosts`, `hosts` and `host_shares_pct`, its hosts and the share of its signal with each, comma-separated. |
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

[MultiQC](https://multiqc.info) custom content: a row of its General
Statistics table for the run (`id: metator`), with the numbers that compare
across runs: the length of the contigs and their N50, the 3D ratio
(`trans_contig_ratio`), the share of the contigs in a MAG (`in_mags_pct`),
the HQ, MQ, LQ and contaminated MAGs, and, after `metator mge` and
`metator host`, the MGE clusters and those with a host, named as in
`<name>.report.txt`. MultiQC finds files ending in `_mqc.yaml`: `multiqc` on
a directory of output directories gathers the runs in one table.

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

A page to read the run in a browser, from the output directory: its numbers
and tables, and links to the FASTA and cool files of the MAGs, the tables,
the plots and the other two report files, when they are beside it (as
`export` and `pipeline` write them). Its interactive plots are drawn by
[plotly.js](https://plotly.com/javascript/), a pinned version loaded from its
CDN: offline, the page shows its numbers and tables, and says its plots need
a connection.

## Logs

Each command writes its log, `<name>.metator_<date>.log` (date as
`%Y%m%d%H%M%S`, e.g. `sample.metator_20261123225838.log`), in its output
directory; `mge` writes it in the current directory. `<name>` is that of the
[report](#the-report-of-a-run). The log starts with the version of metaTOR,
the date and the command line, then follows what the run did.
