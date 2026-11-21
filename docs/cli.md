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
| `-F, --force` | Replace the outputs of an earlier run, in the output directory and at `--output-map`. Without it, a command refuses to overwrite them. A run from a map updates that map in place without it, and warns, before any work, of the clusters and clustering rounds that it replaces there: `--output-map` writes the map of the run to another file instead, and leaves the map as it was. |
| Writing a map | Through metacooler: a run holds the lock of the map it writes, `<map>.lock` (FILE's with `--output-map`), from its checks to its end, and writes the map in an edit session, one staged copy next to it (`.<map>.<uuid>.<n>.d/`), validated, then renamed over it, so a run that stops leaves the map as it was. `MAP` given as a symbolic link is updated where the link points; a symbolic link at `--output-map` is replaced (`-F`), not its target, and a hard link of `MAP` there is refused. A map that another process holds, a map without write permission, or one whose directory is not writable (with room for a copy of the map and for the new map), is refused before any work: write a copy with `--output-map`. A run killed by SIGTERM (a job scheduler's time limit) or SIGKILL leaves its lock and its staged copy, which the next run that writes the map on the same host takes over and removes; on another host, delete the files that the error names. |
| `-T, --tmpdir` | Each run works in a directory of its own inside it, removed at the end, unless `--no-clean-up` or the run fails: the log says which. |
| `-N, --no-clean-up` | Keep the temporary directory of the run: that of a `validation` or `pipeline` run holds its tables. |
| Log files | Each command writes its log, `metator_<command>_<date>.log`, into its output directory; `mge`, into the current directory. |

New to metaTOR? Start with the [Quickstart](quickstart.md); for hits, core,
overlapping and recursive clusters, see [Concepts](concepts.md). metaTOR 2.0
removed `qc`, `contactmap`, `metacool`, `scaffold` and `pairs`, and the plots:
[Other tools](other-tools.md) says what replaces each one.

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
hicstuff's next release), or metacooler from pairs files
(`metacooler create pairs`): see
[Other tools](other-tools.md#from-the-reads-to-the-map).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | A `.metacool` contact map of the contigs, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`): the run reads its network from it, and writes the map of its clusters over it, with the network of its contigs, as `network` stores it. See [Starting from a contact map](guides/metacool.md#starting-from-a-contact-map). |
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
| `-o, --outdir` | the current directory | The output directory. |
| `-p, --prefix` | none | Name the clusters `<prefix>_metator_…`. Letters, digits, `_` and `-` only. |
| `-F, --force` | off | Replace `final_clusters/` of an earlier run, remove the tables and `recursive_clusters/` that an earlier version wrote in the output directory, and replace an existing `--output-map` file. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run works in a directory of its own: the FASTA files of the clusters that miComplete checks (at most the size of the assembly), miComplete's results, the tables of the validation, and, from a compressed assembly, its decompressed copy. Not in `final_clusters/`, which the run makes anew, nor in `recursive_clusters/`, which `-F` removes. |
| `-t, --threads` | `1` | Threads for the partition and the validation (miComplete). |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run, which the log names: the tables of the validation (`cluster_summary.txt`, `contig_data_final.txt`, `binning.txt`, `clusters.txt`, `cluster_members.txt`), miComplete's results and the FASTA files of the clusters that it checked. |
| `--output-map` | none | Write the contact map of the clusters to this file (replaced only with `-F`), and leave `MAP` as it was. By default, the run updates `MAP` in place. |

The validation is the step that needs the most memory: miComplete runs
prodigal and hmmer on every cluster.

```bash
# The map of the contigs, from the reads, with hicstuff
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz

# The clusters: hicstuff/sample.metacool becomes the map of the clusters
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -d depth.txt -t 16 -o binning

# The same, leaving the map of hicstuff as it was
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -o binning --output-map binning/sample.metacool
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
network ([`metator dump MAP`](#dump) refuses the map too). The map is updated
in place, or written to `--output-map`. It writes no table: `network.txt` and
`contig_data_network.txt` are gone. [`metator dump`](#dump) writes the
network, `network.txt`, from the map of the clusters of a run; metacooler
writes it from any map that holds it (see
[The network](formats.md#the-network)).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | As in `pipeline`: the network is stored in it, in place, or in `--output-map`, which an entry of a collection needs. |
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip. |
| `-d, --depth` | none | As in `pipeline`. |
| `-e, --enzyme` | none | As in `pipeline`. |
| `-n, --normalization` | `empirical_hit` | As `--normalization` (`-m`) in `pipeline`. |
| `-o, --outdir` | the current directory | The output directory, of the log of the run. |
| `-F, --force` | off | Replace an existing `--output-map` file. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-N, --no-clean-up` | off | Keep the temporary files: the contig data that the run stores in the map, `contig_data_network.txt`. |
| `--output-map` | none | Store the network in a copy of `MAP` written to this file (replaced only with `-F`), and leave `MAP` as it was. By default, the run updates `MAP` in place. |

```bash
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
among it (`network_digest`): the map is updated in place, or written to
`--output-map`. Its genomes and contigs tables are left as they were.
`validation` checks the overlapping clusters that reach its `--size`. See
[Step by step](guides/metacool.md#step-by-step).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The `.metacool` contact map of the contigs of the run, which `metator network` stored the network in, and which the run stores its partition in; or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`), with `--output-map`, the map it writes it to. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `100` | The iterations of the partition. |
| `-O, --overlap` | `80` | The overlap threshold, in percent. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-o, --outdir` | `./` | The output directory, of the log of the run. |
| `-p, --prefix` | none | Name the overlapping clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace an existing `--output-map` file. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-t, --threads` | `1` | Threads for the partition. |
| `-N, --no-clean-up` | off | Keep the temporary files. |
| `--output-map` | none | Write the map with the partition to this file (replaced only with `-F`), and leave `MAP` as it was. By default, the run updates `MAP` in place. |

```bash
metator partition sample.metacool -t 16 -o binning
```

### `validation` { #validation }

Check the clusters with miComplete, and partition the contaminated ones
again: the third step. Writes `final_clusters/`, and stores the final
clusters, their summary and every clustering round in the map; the tables
of the validation (`cluster_summary.txt`, `contig_data_final.txt`,
`binning.txt`, `clusters.txt`, `cluster_members.txt`, miComplete's results)
and `recursive_clusters/` stay in its temporary directory, which `-N`
keeps, and [`metator dump`](#dump) writes the tables from the map. An
earlier version wrote them in the output directory: there, they are refused
without `-F`, which removes them. **Guide:**
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

The run partitions again the network that `metator network MAP` stored in the
map, and refuses a map without it, or without the partition of `metator
partition MAP`, or whose partition was made on another network than the one
stored in it (`metator network MAP` ran again since): run `metator partition
MAP` again. It then writes the map of its clusters, with all its
clustering rounds and the normalization of its network, as `metator pipeline
MAP` does at its end: in place, or to `--output-map`. See
[Step by step](guides/metacool.md#step-by-step).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The `.metacool` contact map of the contigs of the run, which `metator network` stored the network in, and `metator partition` its partition, and which becomes the map of its clusters; or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`), with `--output-map`, the map it writes them to. |
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip; a compressed one is decompressed once, into the temporary directory of the run, for the FASTA files of the clusters. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `10` | The iterations of the recursive partition. |
| `-O, --overlap` | `90` | The overlap threshold of the recursive partition, in percent. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, that an overlapping cluster of the partition must reach to be checked, and that a cluster of the recursive partition must exceed. |
| `-o, --outdir` | `./` | The output directory. |
| `-p, --prefix` | none | Name the recursive clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace `final_clusters/`, remove the tables and `recursive_clusters/` that an earlier version wrote in the output directory, and replace an existing `--output-map` file. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run writes, in a directory of its own, the clusters that miComplete checks (at most the size of the assembly), miComplete's results and the tables of the validation, and, from a compressed assembly, its decompressed copy. Not in `final_clusters/`, which the run makes anew, nor in `recursive_clusters/`, which `-F` removes. |
| `-t, --threads` | `1` | Threads for the partition and miComplete. |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run, which the log names: the tables of the validation (`cluster_summary.txt`, `contig_data_final.txt`, `binning.txt`, `clusters.txt`, `cluster_members.txt`), miComplete's results and the FASTA files of the clusters that it checked. |
| `--output-map` | none | Write the map of the clusters to this file (replaced only with `-F`), and leave `MAP` as it was. By default, the run updates `MAP` in place. |

```bash
# The map of a run step by step: it becomes the map of the clusters
metator validation sample.metacool -a assembly.fa -t 16 -o binning
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
`summary_formats`): run `metator validation` on it again. A map whose
clusters were edited after the run (by `metacooler reassign`: a final
cluster that is no longer a genome of the map, or a genome that is no final
cluster) is refused too, naming them: its summary and its rounds no longer
describe its genomes. Run `metator validation` on it again, or write the
edited assignment with `metacooler dump -t contigs MAP`. A map whose network
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
refuses without `-F`, and removes with it. The tables are written under
temporary names, then renamed onto their own once all are written: a dump
that fails leaves the output directory as it was.

| Table | What it holds |
|---|---|
| `contig_data_final.txt` | One row per contig of the map, in its order: `Name`, `Size`, `GC_content`, `Hit`, `Shotgun_coverage`, `Restriction_site`, then the ID, number of contigs and size of its core, overlapping and recursive clusters, then `Final_cluster` (see [The contig data tables](formats.md#the-contig-data-tables)). No `ID` column: the ID of a contig was its position in the map that the run read. |
| `cluster_summary.txt` | The summary of the final clusters, one row per cluster: its values are those of the genomes table of the map, written as the `validation` entry of its metadata records it, its rows in order (`summary_rows`), its columns in order (`summary_columns`) and the format of each (`summary_formats`). |
| `binning.txt` | Each contig of a final cluster and its cluster, without a header: the binning that anvi'o imports. |
| `clusters.txt`, `cluster_members.txt` | Every clustering round of the map (`core`, `overlapping`, `recursive_<k>`, `final`, the round `mge` that `metator mge` stores, and `mge_random` when `metator mge --random` stores one), every cluster with its parent and the columns its round stored, and the contigs of each cluster. |
| `network.txt` | The network of the contigs, without a header: `contig1`, `contig2` (their names) and the weight of their edge, the normalization that `metator network` stored (the count of their contacts for `None`); one row per pair of contigs in contact, in the order of the map, without the contacts within a contig (see [The network](formats.md#the-network)). |

The values are those of the tables that the run kept in its temporary
directory: `cluster_summary.txt` is the same file, and so is `clusters.txt`,
but for the rounds that a later command stored in the map (`metator mge
--map`), which it lists after those of the run, with `cluster_members.txt`.
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
metator dump sample.metacool -o tables
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
contigs in: that of `network MAP`, `pipeline MAP` or `validation MAP` (their
`--output-map`, if they wrote one), or an entry of a `.mmetacool` collection
(`file.mmetacool::entries/<name>`, with `--output-map`). The contacts of the
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
and contigs tables are left as they were. The map is updated in place, or
written to `--output-map`, as `metator partition MAP` writes it (see
[Step by step](guides/metacool.md#step-by-step)), and refused as it refuses
it, before any work; so is an assembly that is not a FASTA file, plain or
compressed with gzip, as `validation` refuses it, and the summary of the MGE
clusters that an earlier version wrote in the output directory
(`mges_bin_summary.tsv`, or `mges_cluster_summary.tsv` of a development
build), without `-F`, which removes it. `--random` stores the round `mge_random` and writes
`mges_binned_random.fa`, and leaves the round `mge` as it was. A run that
replaces the clusters of the map (`partition`, `validation`, `pipeline`)
drops the rounds `mge` and `mge_random`, and their entries, with a warning:
run `mge` again on its map.

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip. |
| `-m, --mges` | *required* | The MGE contigs, one name per line. |
| `-s, --threshold-bin` | `0.8` | The score above which two contigs, or clusters, are merged. |
| `-r, --random` | off | Shuffle the contigs between the clusters: a random clustering of the same shape, stored as the round `mge_random` and written to `mges_binned_random.fa`. |
| `-o, --outdir` | the current directory | The output directory. |
| `-F, --force` | off | Replace an existing `--output-map` file, and remove the summary of the MGE clusters that an earlier version wrote in the output directory (`mges_bin_summary.tsv` or `mges_cluster_summary.tsv`), which the round `mge` of the map replaces. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where metacooler copies a compressed assembly, plain, while it reads it. `mge` writes no temporary file of its own. |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run: `mge` leaves nothing in it. |
| `--output-map` | none | Write the map with the MGE clusters to this file (replaced only with `-F`), and leave `MAP` as it was. |

```bash
# The map of the run, after network, pipeline or validation
metator mge sample.metacool -a assembly.fa -m mges.txt -o mge
```
