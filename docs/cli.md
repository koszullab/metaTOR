# CLI reference

The `metator` command has one subcommand per task. This page documents every
command and option; `-h/--help` prints the same in a terminal.

```console
metator --help
metator <command> --help
metator --version
```

## Conventions

| Convention | Meaning |
|---|---|
| Several libraries | One map, which sums them: make the map of each library, then merge them with `metacooler merge` (see [Other tools](other-tools.md)). |
| Several enzymes | Separated by commas: `DpnII,HinfI`. |
| The map of a run | A run writes its map to its output directory (`-o`), under the name of `MAP` (`NAME.metacool` for an entry of a collection), and leaves `MAP` as it was: metacooler copies `MAP` next to that map, the run writes what it stores into the copy, and metacooler renames it into place once it is validated. `--inplace` updates `MAP` itself instead, and so does a run whose output directory holds `MAP` already: the steps of a run, one by one, in one output directory (`network MAP -o DIR`, then `partition DIR/MAP -o DIR`...). |
| `-F, --force` | Replace what the command writes itself, when an earlier run left it: the map of the output directory, the final clusters (`final_clusters/`) and their plots (in `plot/`, whose other files are left as they are), the tables of a dump. Without it, a command refuses to overwrite them. The other results of an earlier run in the output directory (the tables of a development build of 2.0, the outputs of 1.x, an earlier summary of `mge`) are refused, with or without it, and never removed. A run that updates its map in place warns, before any work, of the clusters and clustering rounds that it replaces there. |
| Writing a map | Through metacooler: a run holds the lock of the map it writes, `<map>.lock` (that of the output directory, or `MAP` in place), from its checks to its end, and writes the map in an edit session, one staged copy next to it (`.<map>.<uuid>.<n>.d/`), validated, then renamed onto it, so a run that stops leaves the map as it was. `MAP` updated in place through a symbolic link is updated where the link points; a symbolic link in the output directory, where the map of the run goes, is replaced (`-F`), not its target, and a hard link of `MAP` there is refused. A map that another process holds, a map updated in place without write permission, or one whose directory is not writable (with room for a copy of the map and for the new map), is refused before any work. A run killed by SIGTERM (a job scheduler's time limit) or SIGKILL leaves its lock and its staged copy, which the next run that writes the map on the same host takes over and removes; on another host, delete the files that the error names. |
| `-T, --tmpdir` | Each run works in a directory of its own inside it, removed at the end, unless `--no-clean-up` or the run fails: the log says which. |
| `-N, --no-clean-up` | Keep the temporary directory of the run: that of a `validation` or `pipeline` run holds its tables. |
| Log files | Each command writes its log, `metator_<command>_<date>.log`, into its output directory; `mge`, into the current directory. The log starts, as hicstuff's does, with the version, the date, the command line (`cmd`) and the value of each option, then a line of dashes; `pipeline` marks its three steps (`[1/3] network`...), the partition counts its iterations as they go, the validation each recursive step (`Recursive step 1/9`); and the log ends with what the run wrote and how long it took. The bookkeeping of a run (the files it writes on the way, metacooler's accounts of its work) is logged at the DEBUG level, which the log does not show. |

New to metaTOR? Start with the [Quickstart](quickstart.md); for hits, core,
overlapping and recursive clusters, see [Concepts](concepts.md). metaTOR 2.0
removed `qc`, `contactmap`, `metacool`, `scaffold` and `pairs`, and the heatmap
and the GC and coverage plots of the clusters: [Other tools](other-tools.md)
says what replaces each one.

## Clustering

### `pipeline` { #pipeline }

Run the three steps of the clustering, `network`, `partition` and `validation`,
in one command, and write the contact map of the clusters. **Guide:**
[Running the pipeline](guides/pipeline.md).

```console
metator pipeline [OPTIONS] MAP
```

The run starts from a `.metacool` contact map of the contigs, `MAP`, which
hicstuff makes from the reads (`hicstuff pipeline -M metacool`, from
hicstuff 4.0.0), or metacooler from pairs files
(`metacooler create pairs`): see
[Other tools](other-tools.md#from-the-reads-to-the-map).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | A `.metacool` contact map of the contigs, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`): the run reads its network from it, and writes the map of its clusters, with the network of its contigs, as `network` stores it, to the output directory, under the name of `MAP` (see [The map of a run](#conventions)). See [Starting from a contact map](guides/metacool.md#starting-from-a-contact-map). |
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip; a compressed one is decompressed once, into the temporary directory of the run, for the FASTA files of the clusters. |
| `-d, --depth` | none | The coverage of the contigs, from `jgi_summarize_bam_contig_depths` of MetaBAT2. The `Shotgun_coverage` column of the map gives it without `--depth`, which replaces it: see [Starting from a contact map](guides/metacool.md#starting-from-a-contact-map). |
| `-e, --enzyme` | none | The restriction enzymes of the library. The `Restriction_site` column of the map gives the restriction sites without `--enzyme`, which replaces it. |
| `-m, --normalization` | `empirical_hit` | `None`, `abundance`, `length`, `RS`, `empirical_hit` or `theoritical_hit`. `RS` and `theoritical_hit` need `--enzyme`; `abundance` and `theoritical_hit` need `--depth`, or the columns of the map that give them. See [Normalising the network](guides/network.md#normalising-the-network). |
| `-A, --algorithm` | `louvain` | The partition algorithm: `louvain` or `leiden`. |
| `-i, --iterations` | `100` | The iterations of the partition. |
| `-O, --overlap` | `80` | The percentage of iterations in which two core clusters must fall together to be merged into one overlapping cluster. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, that an overlapping cluster must reach to be checked with miComplete, and that a cluster of the recursive partition must exceed. |
| `-j, --rec-iter` | `10` | The iterations of the recursive partition of the contaminated clusters. |
| `-P, --rec-overlap` | `90` | The overlap threshold of the recursive partition, in percent. |
| `-o, --outdir` | the current directory | The output directory: the map of the clusters, `final_clusters/`, `plot/` (see [`validation`](#validation)) and the log. |
| `-p, --prefix` | none | Name the clusters `<prefix>_metator_…`. Letters, digits, `_` and `-` only. |
| `-F, --force` | off | Replace the map, `final_clusters/` and the plots (in `plot/`, whose other files are left as they are) of an earlier run in the output directory. The other results of an earlier run in the output directory (the tables and `recursive_clusters/` of a development build of 2.0, the outputs of 1.x) are refused, with or without `-F`, and never removed. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run works in a directory of its own: the FASTA files of the clusters that miComplete checks (at most the size of the assembly), miComplete's results, the tables of the validation, and, from a compressed assembly, its decompressed copy. Not in `final_clusters/`, which the run makes anew, nor in a directory of the results of an earlier run (`recursive_clusters/`, `final_bin/`...), which the run refuses. |
| `-t, --threads` | `1` | Threads for the partition and the validation (miComplete). |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run, which the log names: the tables of the validation (`cluster_summary.txt`, `contig_data_final.txt`, `binning.txt`, `clusters.txt`, `cluster_members.txt`), miComplete's results and the FASTA files of the clusters that it checked. |
| `--inplace` | off | Update `MAP` in place with the map of the clusters, rather than write it to the output directory: refused for an entry of a collection. A run whose output directory holds `MAP` already updates it in place without `--inplace`. |

The validation is the step that needs the most memory: miComplete runs
Prodigal and HMMER (pyrodigal and pyhmmer) on every cluster.

```bash
# The map of the contigs, from the reads, with hicstuff
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz

# The clusters: binning/sample.metacool, the map of the clusters; hicstuff/sample.metacool is left as it was
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -d depth.txt -t 16 -o binning

# The same, updating hicstuff/sample.metacool in place
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -o binning --inplace
```

### `network` { #network }

Build the network of contigs from the map, and store it in the map: the first
step of the pipeline, whose network `partition` and `validation` read from the
map. **Guide:** [Building the network](guides/network.md).

```console
metator network [OPTIONS] MAP
```

The run stores in the map what it measured on each contig, in its contigs
table (`GC_content`, `Hit`, and `Shotgun_coverage` and `Restriction_site` when
it has them), and the normalization of the network, as the weights
`metator_<normalization>` of its bins (none for `None`), which replace those of
an earlier network; and its command line, under the `network` entry of the
`metator` key of the metadata of the map. The clusters of the map are kept,
its partition too: when the run stores another network than the one the
partition was made on, it says so with a warning, and `metator validation
MAP` refuses the partition until `metator partition MAP` partitions the new
network ([`metator dump MAP`](#dump) refuses the map too). The run writes the
map to its output directory, under the name of `MAP`, or updates `MAP` in
place with `--inplace` (see [The map of a run](#conventions)). It writes no
table in the output directory: `network.txt` is gone, and
`contig_data_network.txt`, the contig
data that it stores in the map, is written in its temporary directory, which
`-N` keeps. [`metator dump`](#dump) writes the network, `network.txt`, from
the map of the clusters of a run; metacooler writes it from any map that holds
it, one that only `network` ran on included (see
[The network](formats.md#the-network)).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | As in `pipeline`: the run writes it, with the network of its contigs, to the output directory. |
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip. |
| `-d, --depth` | none | As in `pipeline`. |
| `-e, --enzyme` | none | As in `pipeline`. |
| `-n, --normalization` | `empirical_hit` | As `--normalization` (`-m`) in `pipeline`. |
| `-o, --outdir` | the current directory | The output directory, of the map of the run and of its log. |
| `-F, --force` | off | Replace the map of an earlier run in the output directory. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run writes, in a directory of its own, the contig data that it stores in the map (`contig_data_network.txt`). |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run, which the log names, and its `contig_data_network.txt`. |
| `--inplace` | off | Store the network in `MAP` itself, in place, rather than in the map of the output directory: refused for an entry of a collection. |

```bash
# binning/sample.metacool, with the network of the contigs
metator network sample.metacool -a assembly.fa -e DpnII -o binning
```

### `partition` { #partition }

Partition the network that `network` stored in the map, many times, into core
and overlapping clusters: the second step. It stores them in the map, and
writes no table and no FASTA file. **Guide:**
[Partition and validation](guides/binning.md).

```console
metator partition [OPTIONS] MAP
```

The run reads the network from the map, which `metator network MAP` stored
there, and refuses a map without it. It stores every core cluster
(`core_00001`...) and every overlapping cluster
(`<prefix>_metator_00001_00000`...), whatever its size, with their contigs, as
the rounds `core` and `overlapping` of the clustering rounds of the map, in
place of those it held, and what the run was as the `partition` entry of the
`metator` key of the metadata of the map, the network that it partitioned
among it (`network_digest`), in the map that it writes to its output
directory (`MAP` itself when it is there already, or with `--inplace`). Its
genomes and contigs tables are left as they were.
`validation` checks the overlapping clusters that reach its `--size`. See
[Step by step](guides/metacool.md#step-by-step).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The `.metacool` contact map of the contigs of the run, which `metator network` stored the network in, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`): the run writes it, with its partition, to the output directory. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `100` | The iterations of the partition. |
| `-O, --overlap` | `80` | The overlap threshold, in percent. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-o, --outdir` | `./` | The output directory, of the map of the run and of its log. |
| `-p, --prefix` | none | Name the overlapping clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace the map of an earlier run in the output directory. |
| `-T, --tmpdir` | `./tmp` | The temporary directory: `partition` writes no temporary file in it. |
| `-t, --threads` | `1` | Threads for the partition. |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run: `partition` leaves nothing in it. |
| `--inplace` | off | Store the partition in `MAP` itself, in place, rather than in the map of the output directory: refused for an entry of a collection. |

```bash
# The map that metator network wrote in binning, updated in place
metator partition binning/sample.metacool -t 16 -o binning
```

### `validation` { #validation }

Check the clusters with miComplete, and partition the contaminated ones
again: the third step. Writes `final_clusters/` and the plots of the final
clusters (`plot/`), and stores the final
clusters, their summary and every clustering round in the map; the tables
of the validation (`cluster_summary.txt`, `contig_data_final.txt`,
`binning.txt`, `clusters.txt`, `cluster_members.txt`, miComplete's results)
and `recursive_clusters/` stay in its temporary directory, which `-N`
keeps, and [`metator dump`](#dump) writes the tables from the map. An
earlier version wrote them in the output directory: there, they are refused,
with or without `-F`, and never removed. **Guide:**
[Partition and validation](guides/binning.md#validation).

```console
metator validation [OPTIONS] MAP
```

The clusters that the run checks are the overlapping clusters of the
partition that `metator partition MAP` stored in the map that reach `--size`:
it writes their FASTA files from the assembly, and the contig data of the
partition, in its temporary directory (`-N` keeps them), which needs room for
them (at most the size of the assembly) and, from a gzip-compressed
assembly, for its decompressed copy, made once for the run. A cluster is
partitioned again when its weighted completeness reaches 33% and its
contamination, (weighted redundancy − 1) / weighted completeness, reaches
5%. Its parts replace it when the most complete of them keeps enough of its
completeness.

Once the clusters are final, the run draws two plots of them in `plot/`, in
its output directory, from miComplete's numbers: `bins_distribution.png`, the
completeness and the contamination (weighted redundancy − 1, drawn from 0 to
105%) of each final cluster, the clusters by decreasing completeness; and
`bins_size_distribution.png`, a pie chart of the size of the clusters of each
quality class (HQ, MQ and LQ MAGs, contaminated clusters and the others, the
classes that the log counts: see [Cluster quality](concepts.md#cluster-quality))
and of the contigs in no cluster, of the contigs of the map in all. A plot
that cannot be drawn is a warning, and the run goes on. The plots of an
earlier run are refused without `-F`, which replaces them; the other files of
`plot/` are left as they are, but for the plots of 1.x that 2.0 no longer
draws, refused with or without `-F` as the other outputs of 1.x.

The run partitions again the network that `metator network MAP` stored in the
map, and refuses a map without it, or without the partition of `metator
partition MAP`, or whose partition was made on another network than the one
stored in it (`metator network MAP` ran again since): run `metator partition
MAP` again. It then writes the map of its clusters, with all its
clustering rounds and the normalization of its network, as `metator pipeline
MAP` does at its end: to its output directory (`MAP` itself when it is there
already, or with `--inplace`). See
[Step by step](guides/metacool.md#step-by-step).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The `.metacool` contact map of the contigs of the run, which `metator network` stored the network in, and `metator partition` its partition, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`): the run writes the map of its clusters from it, to the output directory. |
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip; a compressed one is decompressed once, into the temporary directory of the run, for the FASTA files of the clusters. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `10` | The iterations of the recursive partition. |
| `-O, --overlap` | `90` | The overlap threshold of the recursive partition, in percent. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, that an overlapping cluster of the partition must reach to be checked, and that a cluster of the recursive partition must exceed. |
| `-o, --outdir` | `./` | The output directory: the map of the clusters, `final_clusters/`, `plot/` and the log. |
| `-p, --prefix` | none | Name the recursive clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace the map, `final_clusters/` and the plots (in `plot/`, whose other files are left as they are) of an earlier run in the output directory. The other results of an earlier run in the output directory (the tables and `recursive_clusters/` of a development build of 2.0, the outputs of 1.x) are refused, with or without `-F`, and never removed. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run writes, in a directory of its own, the clusters that miComplete checks (at most the size of the assembly), miComplete's results and the tables of the validation, and, from a compressed assembly, its decompressed copy. Not in `final_clusters/`, which the run makes anew, nor in a directory of the results of an earlier run (`recursive_clusters/`, `final_bin/`...), which the run refuses. |
| `-t, --threads` | `1` | Threads for the partition and miComplete. |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run, which the log names: the tables of the validation (`cluster_summary.txt`, `contig_data_final.txt`, `binning.txt`, `clusters.txt`, `cluster_members.txt`), miComplete's results and the FASTA files of the clusters that it checked. |
| `--inplace` | off | Write the map of the clusters over `MAP`, in place, rather than to the output directory: refused for an entry of a collection. |

```bash
# The map of a run step by step, in binning: it becomes the map of the clusters
metator validation binning/sample.metacool -a assembly.fa -t 16 -o binning
```

## After the clustering

### `dump` { #dump }

Write the tables of a run from the map of its clusters: those of its
validation, which a run keeps in its temporary directory only (with `-N`),
and its network. **API:** `metator.dump.dump`.

```console
metator dump [OPTIONS] MAP
```

The map is the one that `metator validation MAP` or `metator pipeline MAP`
wrote, which holds every result of the run: the network that `metator network`
stored, the partition that `metator partition` stored, and the clustering
rounds, the final clusters and their summary that the validation stored. It
is read only, without a lock, and the tables are of one state of it: a map
that another run replaces or writes while it is read is refused, and no
table is written (run the dump again once that run has ended). A map without
the clusters of a validation (the `validation` entry of the `metator` key of
its metadata, and the clustering round `final`) is refused, and so is one
without the network or the partition of its run, each with the step to run,
and one that a development build of metaTOR 2.0 validated, which does not
record how the summary of its clusters is written (`summary_rows`,
`summary_formats`): run `metator validation` on it again. The tables follow
the genomes of the map, which need not be the final clusters of the run: the
MGE genomes that `metator mge --reassign` makes (see [`mge`](#mge)), or the
genomes of a curation by hand with `metacooler reassign`, each with the values
that the genomes table gives it (those of the validation for a final
cluster, which describe it with the contigs that it had). A map whose network
is not the one that its clusters were made from (`metator network` ran on it
after the run, with another normalization: the `network_digest` of its
`partition` entry, see [`partition`](#partition)) is refused, since its
network and its clusters would not go together: run `metator network -n
<normalization>` with the normalization of the run, which stores its network
again, or `metator partition` and `metator validation` again, which cluster
the network stored now. The output directory is made if it is not there;
one that holds a table of an earlier dump is refused, without `-F`, which
replaces them (a symbolic link there is replaced, not its target), and so is
a directory where a table goes. Give the dump a directory of its own, not
the output directory of a run: a `validation` or `pipeline` run there would
take its tables for those that an earlier version wrote there, which it
refuses, with or without `-F`. The tables are written under
temporary names, then renamed onto their own once all are written: a dump
that fails leaves the output directory as it was.

| Table | What it holds |
|---|---|
| `contig_data_final.txt` | One row per contig of the map, in its order: `Name`, `Size`, `GC_content`, `Hit`, `Shotgun_coverage`, `Restriction_site`, then the ID, number of contigs and size of its core, overlapping and recursive clusters, then `Final_cluster`, its genome in the map (`ND` for none; see [The contig data tables](formats.md#the-contig-data-tables)). No `ID` column: the ID of a contig was its position in the map that the run read. |
| `cluster_summary.txt` | The summary of the genomes of the map, one row per genome: its final clusters, in the order of the `validation` entry of its metadata (`summary_rows`), then its other genomes, in the order of its genomes table (the MGE genomes of `metator mge --reassign`). The values are those of the genomes table, written as the `validation` entry records them: its columns in order (`summary_columns`), each in its format (`summary_formats`). A map with MGE genomes adds their columns, `GC`, `Hit` and `BinningScore` (empty for the final clusters), and `mge`, `True` for an MGE genome and `False` for a final cluster. After [`metator host`](#host), the number of MGE clusters that each MAG hosts and their names (`n_mges`, `mges`, empty for a genome that is no MAG), last. |
| `binning.txt` | Each contig of a genome of the map and its genome, without a header: the binning that anvi'o imports, its final clusters, and its MGE genomes after `metator mge --reassign`. |
| `clusters.txt`, `cluster_members.txt` | Every clustering round of the map (`core`, `overlapping`, `recursive_<k>`, `final`, the round `mge` that `metator mge` stores, and `mge_random` when `metator mge --random` stores one), every cluster with its parent, its number of contigs and their length (`n_contigs`, `length`, which the map computes) and the columns its round stored, and the contigs of each cluster. |
| `network.txt` | The network of the contigs, without a header: `contig1`, `contig2` (their names) and the weight of their edge, the normalization that `metator network` stored (the count of their contacts for `None`); one row per pair of contigs in contact, in the order of the map, without the contacts within a contig (see [The network](formats.md#the-network)). |

The values are those of the tables that the run kept in its temporary
directory: `cluster_summary.txt` is the same file until the genomes of the map
change (`metator mge --reassign`, `metacooler reassign`) or `metator host`
stores the hosts of the MGE clusters in it, and so is
`clusters.txt` until a later command stores a round in the map, whose rounds
keep the final clusters of the run, whatever its genomes. After `metator mge`,
`clusters.txt` and `cluster_members.txt` list the rounds `mge` and
`mge_random` after those of the run, and every row of `clusters.txt` gains
their columns, `GC`, `Hit` and `BinningScore` (`NA` for the rounds of the
run); after `metator host`, the hosts of each MGE cluster too (`n_hosts`,
`hosts`, `host_shares`, `NA` for the clusters of the other rounds).
The rows of `contig_data_final.txt`, `binning.txt` and `cluster_members.txt`
come in the order of the map, which a validation regroups by cluster: the
run's tables followed the map as it read it.

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The map of the clusters of a run, which `metator validation` or `metator pipeline` wrote, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`). Read only. |
| `-o, --outdir` | the current directory | The output directory of the tables, made if it is not there. |
| `-F, --force` | off | Replace the tables of an earlier dump that the output directory holds. |

```bash
# The tables of a run, from the map of its clusters
metator dump binning/sample.metacool -o tables
```

### `mge` { #mge }

Cluster the contigs of mobile genetic elements from the contacts between them,
and store the MGE clusters in the contact map.
**API:** `metator.mge.mge_binning`, `metator.mge.write_mges_fasta`,
`metator.metacool.write_map`.
**Guide:** [Clustering mobile genetic elements](guides/mge.md).

```console
metator mge [OPTIONS] MAP
```

`MAP` is a contact map that `metator network` stored the network of its
contigs in: the map that `network`, `pipeline` or `validation` wrote in its
output directory (or updated in place), or an entry of a `.mmetacool`
collection (`file.mmetacool::entries/<name>`). The contacts of the
MGE contigs come from it: between two contigs, its contacts summed by pair of
contigs, whatever its bin size; within a contig, those between reads more
than 1 kb apart, which its contact statistics count (metacooler's QC block,
made from the pairs at 1001 bp: by `hicstuff pipeline -M metacool` or
`metacooler create pairs`, with metacooler 0.4.0 or later). The GC content
and hits of each contig are those that `metator network` stored in it; its
clusters, if it has any, are not read. A map without those contact
statistics, with them measured from its pixels or at another distance, or
without the columns of `network`, is refused before anything is written,
with how to make one; so is a map that holds none of the MGE contigs. An
MGE contig that the map lacks is left out, with a warning. The MGE contigs
are taken the longest first, ties by name, whatever the order of the map:
that order breaks the ties of the scores, and numbers the clusters of one
contig (see [Clustering mobile genetic elements](guides/mge.md#how-it-works)).

The run writes `mges_binned.fa` and stores its MGE clusters in the map: as
the round `mge` of its clustering rounds, with their contigs, length, GC
content, hits and clustering score, in place of the round `mge` of an earlier
run, after the other rounds, which are kept (`metator dump` writes it with
them: see [`dump`](#dump)); the entry `mge` of the `metator` key of its
metadata records the run, and what it read of the map (`read`). Its genomes
and contigs tables are left as they were, but with `--reassign`, which makes
each MGE cluster a genome of the map, named after it (`MetaTOR_MGE_<ID>`),
with its GC content, hits and clustering score (`GC`, `Hit`,
`BinningScore`): its contigs leave the genomes that held them, final clusters
or none. The other genomes are kept, with their other contigs and the values
of their validation, which describe them with the contigs that they lose (the
log names them), and a final cluster whose contigs are all MGE contigs is no
genome any more, with a warning; the column `mge` of the genomes table is
`True` for the MGE genomes, and `False` for the others. `metator dump`
follows them (see [`dump`](#dump)). A run that replaces the round `mge`
undoes the MGE genomes of the map first, those that the column `mge` names
and those named after a cluster of its round `mge` (as `metacooler reassign
--round mge` makes them): their contigs go back to their final cluster, or to
none. So a run with `--reassign` replaces them with its own, and one without
leaves the map with its final clusters alone (a final cluster that had lost
all its contigs comes back without its values, with a warning: run
`validation` again to measure them). `--random` leaves them as they are, and
is refused with `--reassign`. The run writes the map to its
output directory, as `metator partition MAP` writes it (`MAP` itself when it
is there already, or with `--inplace`: see
[Step by step](guides/metacool.md#step-by-step)), and refused as it refuses
it, before any work; so is an assembly that is not a FASTA file, plain or
compressed with gzip, as `validation` refuses it, and the summary of the MGE
clusters that an earlier version wrote in the output directory
(`mges_bin_summary.tsv`, or `mges_cluster_summary.tsv` of a development
build), with or without `-F`: it is never removed. `--random` stores the round `mge_random` and writes
`mges_binned_random.fa`, and leaves the round `mge` as it was. A run that
replaces the clusters of the map (`partition`, `validation`, `pipeline`)
drops the rounds `mge` and `mge_random`, and their entries, with a warning:
run `mge` again on its map (`validation` and `pipeline` make their final
clusters the genomes of the map, without MGE genomes). A run that replaces
the round `mge` (not `--random`) drops the hosts that [`metator host`](#host)
stored in the map, with a warning: run `host` again.

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip. |
| `-m, --mges` | *required* | The MGE contigs, one name per line. |
| `-s, --threshold-bin` | `0.8` | The score above which two contigs, or clusters, are merged. |
| `-r, --random` | off | Shuffle the contigs between the clusters: a random clustering of the same shape, stored as the round `mge_random` and written to `mges_binned_random.fa`. |
| `--reassign` | off | Make each MGE cluster a genome of the map, with its GC content, hits and clustering score: its contigs leave the final clusters that held them, which keep their values. The column `mge` of the genomes table says which genomes are MGE genomes. Refused with `--random`. |
| `-o, --outdir` | the current directory | The output directory: the map of the run and `mges_binned.fa`. |
| `-F, --force` | off | Replace the map of an earlier run in the output directory. A summary of the MGE clusters that an earlier version wrote in the output directory (`mges_bin_summary.tsv` or `mges_cluster_summary.tsv`) is refused, with or without `-F`, and never removed. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where metacooler copies a compressed assembly, plain, while it reads it. `mge` writes no temporary file of its own. |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run: `mge` leaves nothing in it. |
| `--inplace` | off | Store the MGE clusters in `MAP` itself, in place, rather than in the map of the output directory: refused for an entry of a collection. |

```bash
# The map of the clusters of a run, in binning: updated in place
metator mge binning/sample.metacool -a assembly.fa -m mges.txt -o binning
# The same, the MGE clusters made genomes of the map, then the tables:
# binning.txt lists the MAGs and the MGE genomes
metator mge binning/sample.metacool -a assembly.fa -m mges.txt -o binning --reassign
metator dump binning/sample.metacool -o tables
```

### `host` { #host }

Find the MAGs that host the MGE clusters of the map, from the contacts
between them, and store them in the map.
**API:** `metator.host.find_hosts`, `metator.host.write_hosts`,
`metator.metacool.write_map`.
**Guide:** [Clustering mobile genetic elements](guides/mge.md#the-hosts).

```console
metator host [OPTIONS] MAP
```

`MAP` is the map that `metator mge` stored its MGE clusters in (the round
`mge`), which holds the network that `metator network` stored, and the MAGs
as its genomes: the final clusters of `validation` or `pipeline`, or the bins
of another binner, which `metacooler reassign` made the genomes of the map;
or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`).
`__unbinned__` is no MAG, and neither is a genome named after an MGE cluster
(an MGE genome of `metator mge --reassign`). The contigs of an MGE cluster
count for it only: a MAG that holds some of them (a validation clusters the
MGE contigs with the other contigs) is taken without them. A map without the
network of `metator network`, without MGE clusters or without MAGs is refused
before anything is written, with what to run.

The signal of an MGE cluster with a MAG is the sum of the weights of the
edges of the network between their contigs, with the normalization that
`metator network` stored. A MAG hosts the MGE cluster when its signal is at
least `--min-share` percent of the signal of the MGE cluster with all the
MAGs, and the MGE cluster is in contact with at least `--min-contigs` of its
contigs: an MGE cluster has one host, several, or none. Its contacts with the
contigs in no MAG and with the other MGE clusters do not count in that share,
so an MGE cluster with few contacts with the MAGs (a phage that the library
caught in its free particles) can have a host from a small part of its
contacts: `share_of_all`, in `mge_hosts.txt`, says how small.

The run writes `mge_hosts.txt`, every MGE cluster and MAG in contact, and
`mge_hosts.png`, the signal and the share of the hosts (see
[Hosts of the MGE clusters](formats.md#hosts-of-the-mge-clusters)), in its
output directory, in place of those of an earlier run, and stores the hosts
in the map: the number of hosts of each MGE cluster, their names and their
shares, as columns of the round `mge` (`n_hosts`, `hosts`, `host_shares`);
the number of MGE clusters that each MAG hosts and their names, as columns
of the genomes table (`n_mges`, `mges`); and the run, under the entry `host`
of the `metator` key of its metadata, with the digest of the network that it
read (`network_digest`). They replace those of an earlier run; the rounds,
their clusters and contigs, and the rest of the map are left as they were.
The run writes the map to its output directory, as `metator mge MAP` writes
it (`MAP` itself when it is there already, or with `--inplace`: see
[Step by step](guides/metacool.md#step-by-step)), and refused as it refuses
it, before any work. The hosts describe the round `mge`, the genomes and the
network of the map: a run that replaces one of them drops them, with a
warning once its map is written (run `host` again then): `mge`, but with
`--random`; `partition`, `validation` and `pipeline`; and `network`, when it
stores another network than the one the hosts were found on (its
`network_digest`).

| Option | Default | Description |
|---|---|---|
| `-s, --min-share` | `10` | The share of the signal of an MGE cluster with the MAGs, in percent, from which a MAG is a host. |
| `-c, --min-contigs` | `5` | The number of contigs of a MAG in contact with an MGE cluster from which the MAG is a host. |
| `-o, --outdir` | the current directory | The output directory: the map of the run, `mge_hosts.txt`, `mge_hosts.png` and the log. |
| `-F, --force` | off | Replace the map of an earlier run in the output directory. `mge_hosts.txt` and `mge_hosts.png` are replaced without it. |
| `--inplace` | off | Store the hosts in `MAP` itself, in place, rather than in the map of the output directory: refused for an entry of a collection. |

```bash
# The hosts of the MGE clusters that mge stored in the map of a run, in
# binning: updated in place
metator host binning/sample.metacool -o binning
```
