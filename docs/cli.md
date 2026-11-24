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
| Several libraries | One map, which sums them: give the pairs files or the BAM files of every library to `network` or `pipeline`, which sum them in the map that they make (see [The inputs of a run](#inputs)); or make the map of each library, then merge them with `metacooler merge` (see [Other tools](other-tools.md)). |
| Several enzymes | Separated by commas: `DpnII,HinfI`. |
| The map of a run | A run writes its map to its output directory (`-o`), under the name of `MAP` (`NAME.metacool` for an entry of a collection), and leaves `MAP` as it was: metacooler copies `MAP` next to that map, the run writes what it stores into the copy, and metacooler renames it into place once it is validated. `--inplace` updates `MAP` itself instead, and so does a run whose output directory holds `MAP` already: the steps of a run, one by one, in one output directory (`network MAP -o DIR`, then `partition DIR/MAP -o DIR`...). |
| Outputs | Beside its map and its log, a run writes its files in folders of its output directory: `fasta/mags_with-mges/`, the FASTA file of each final cluster, as the validation made it, with the MGE contigs that it holds (`validation`, `pipeline`); `fasta/mags/`, each MAG, a genome of the map but the MGE genomes, as its contigs (`<MAG>.fa`) and stitched (`<MAG>.stitched.fa`: one record, its contigs joined by 200 `N`) (`pipeline`); `fasta/mges/`, each MGE cluster, the same two ways (`mge`, `pipeline --mges`; `fasta/mges_random/` for `mge --random`); `cool/`, a cool file per MAG, with the MGE clusters that it hosts, balanced (`pipeline`); `tables/`, the tables of `dump` (`pipeline`) and `mge_hosts.txt` (`host`, `pipeline --mges`); `plots/`, the plots, as PDF files (`validation`, `host`, `pipeline`). [`export`](#export) writes them all again from the map of a run. |
| `-F, --force` | Replace what the command writes itself, when an earlier run left it: the map of the output directory, the folders that it makes anew (`fasta/mags_with-mges/`, `fasta/mags/`, `fasta/mges/`, `cool/`), the plots of the clusters (in `plots/`, whose other files are left as they are), the tables of a dump. Without it, a command refuses to overwrite them. The other results of an earlier run in the output directory (the tables, `final_clusters/`, `plot/`, `mges_binned.fa` and `mge_hosts.txt` of a development build of 2.0, the outputs of 1.x, an earlier summary of `mge`) are refused, with or without it, and never removed. A run that updates its map in place warns, before any work, of the clusters and clustering rounds that it replaces there. |
| Writing a map | Through metacooler: a run holds the lock of the map it writes, `<map>.lock` (that of the output directory, or `MAP` in place), from its checks to its end, and writes the map in an edit session, one staged copy next to it (`.<map>.<uuid>.<n>.d/`), validated, then renamed onto it, so a run that stops leaves the map as it was. `MAP` updated in place through a symbolic link is updated where the link points; a symbolic link in the output directory, where the map of the run goes, is replaced (`-F`), not its target, and a hard link of `MAP` there is refused. A map that another process holds, a map updated in place without write permission, or one whose directory is not writable (with room for a copy of the map and for the new map), is refused before any work. A run killed by SIGTERM (a job scheduler's time limit) or SIGKILL leaves its lock and its staged copy, which the next run that writes the map on the same host takes over and removes; on another host, delete the files that the error names. |
| `-T, --tmpdir` | Each run works in a directory of its own inside it, removed at the end, unless `--no-clean-up` or the run fails: the log says which. |
| `-N, --no-clean-up` | Keep the temporary directory of the run: that of a `validation` or `pipeline` run holds its tables. |
| Log files | Each command writes its log, `metator_<command>_<date>.log`, into its output directory; `mge`, into the current directory. The log starts, as hicstuff's does, with the version, the date, the command line (`cmd`) and the value of each option, then a line of dashes; `pipeline` marks its steps: three (`[1/3] network`...), one more first from pairs or BAM files (`[1/4] map`...), and two more last with `--mges` (... `[5/5] host`, or `[6/6] host` from pairs or BAM files); the partition counts its iterations as they go, the validation each recursive step (`Recursive step 1/9`); and the log ends with what the run wrote and how long it took. The bookkeeping of a run (the files it writes on the way, metacooler's accounts of its work) is logged at the DEBUG level, which the log does not show. |
| Random numbers | The partition (Louvain or Leiden, in `partition`, `pipeline` and the recursive partition of `validation`) and the random clustering of `mge --random` draw random numbers, other ones at each run: two runs on the same map give close clusters, but not the same. `--seed S` makes a run give the same clusters each time, from the same map and options (and versions of igraph and metaTOR): a step is seeded once, before its iterations, which still draw other numbers one after the other (each iteration of Louvain or Leiden gives another membership, as the core and overlapping clusters need), so that the step as a whole, all its iterations in their order, is the same. `pipeline --seed S` seeds its partition and its validation alike, and gives the clusters of `partition --seed S`, then `validation --seed S`. The metadata of the map store the seed of the run, with its other parameters. |

New to metaTOR? Start with the [Quickstart](quickstart.md); for hits, core,
overlapping and recursive clusters, see [Concepts](concepts.md). metaTOR 2.0
removed `qc`, `contactmap`, `metacool`, `scaffold` and `pairs`, and the heatmap
and the GC and coverage plots of the clusters: [Other tools](other-tools.md)
says what replaces each one.

## The inputs of a run { #inputs }

`network` and `pipeline` start from the contact map of the contigs, or make it
from the pairs files or the BAM files of the run: their arguments,
`INPUT...`.

| Input | What the run does |
|---|---|
| A `.metacool` map, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`) | Reads its network from it, and writes it to the output directory, under its name (see [The map of a run](#conventions)). |
| Pairs files, one per library: plain, or compressed with gzip, bzip2, zstd or LZ4 | Makes the map of the contigs from them first, their contacts summed, as `metacooler create pairs` makes it: one bin per contig of the assembly, every contig in `__unbinned__`, with the contact statistics of the reads, which `mge` reads. |
| BAM files (or SAM files) of both mates of each read (paired-end alignments, as `bwa mem -5SP` writes them), one per library; or of one mate each (as bowtie2 aligns them, one mate at a time), an R1 file and an R2 file per library | Takes their pairs first: the primary alignments of the two mates of each read, both aligned with a mapping quality of at least `--min-quality` (30) on contigs of the assembly, neither flagged as a PCR or optical duplicate (0x400, as Picard MarkDuplicates and samtools markdup flag them) nor as failing QC (0x200), each at the 5' end of its alignment (its first aligned base on the forward strand, its last on the reverse strand), as pairtools and hicstuff's mapping mode 5spm record it. Then makes the map from them, as from pairs files. |

- **Each argument** is a path, a comma-separated list or a glob pattern, which
  stands for the files it matches, sorted by name, as hicstuff 4.0 reads its
  inputs: `a.pairs b.pairs`, `a.pairs,b.pairs` and `'*.pairs'` (quoted, or
  expanded by the shell) are the same two libraries. A path that exists is
  that file, whatever its name holds.
- **BAM files of one mate each** come in two halves: two arguments, the R1
  files then the R2 files (`'*_R1.bam' '*_R2.bam'`, or two comma-separated
  lists), or as many R1 arguments as R2 arguments, as the shell expands
  `*_R1.bam *_R2.bam`. The i-th R1 file and the i-th R2 file are the mates of
  one library, matched by read name (a mate suffix, `/1` or `/2`, left out):
  in the order the aligner wrote both, or sorted by name (`samtools sort -n`
  or `-N`); in another order, the run sorts them by name first, with the
  samtools of pysam, in its temporary directory.
- **BAM files of both mates** are read as the aligner writes them, the two
  mates of a read in a row, or sorted by name; sorted by coordinate (as their
  header says, or as the run finds when a read comes back after others), the
  run groups their reads by name first (`samtools collate`, in its temporary
  directory).
- **What each file is** comes from its first bytes, whatever its name: an
  HDF5 file is a map, a BAM or SAM file is alignments, and any other file is a
  pairs file. Reads (FASTQ), an assembly (FASTA) and CRAM files are refused;
  so are a map given with other inputs, pairs files with BAM files, BAM files
  of both mates with BAM files of one mate each, and a file given twice.
- **The map that the run makes** goes to its output directory, as
  `<prefix>.metacool` (`pipeline --prefix`), or `<assembly>.metacool`, named
  after the FASTA file of the assembly without its extensions
  (`assembly.fa.gz` makes `assembly.metacool`). The run then updates it in
  place, as the map of its output directory: a `pipeline` run makes it the map
  of its clusters. A map there already is replaced only with `-F`, once the
  new map is written; `--inplace` is refused. A run that stops later leaves
  the map of the contigs, which a new run can start from without reading the
  pairs again.
- **The sequences** of the BAM files (their `@SQ` lines) and of the pairs
  files (their `#chromsize` lines) must be contigs of the assembly, with their
  lengths: a sequence of another length, or files none of whose sequences is
  a contig of the assembly, are refused before the pairs are read. A sequence
  that is not a contig of the assembly (a spike-in) is left out, with a
  warning, and so are the pairs with a read on it. Inputs that give no contact
  on the contigs of the assembly are refused: the run writes no map, and
  leaves one there as it was.
- **The pairs of the BAM files** are written in the temporary directory of the
  run, a file per library (`library_<i>.pairs`, about 25 bytes per pair, plus
  the names of its two contigs, which `-N` keeps), then read once by
  metacooler, with those of the other
  libraries. `--threads` reads the BAM files and bins the pairs on several
  threads.

## Clustering

### `pipeline` { #pipeline }

Run the three steps of the clustering, `network`, `partition` and `validation`,
in one command, and write the contact map of the clusters; with `--mges`, go
on with `mge` and `host` on that map. **Guide:**
[Running the pipeline](guides/pipeline.md).

```console
metator pipeline [OPTIONS] INPUT...
```

The run starts from a `.metacool` contact map of the contigs, `MAP`, which
hicstuff makes from the reads (`hicstuff pipeline -M metacool`, from
hicstuff 4.0.0), or metacooler from pairs files
(`metacooler create pairs`): see
[Other tools](other-tools.md#from-the-reads-to-the-map). Or it makes that map
first, from the pairs files or the BAM files of the run, in its output
directory (see [The inputs of a run](#inputs)): a step more, which the log
marks `[1/4] map`.

The run writes the map of its clusters at the end of the validation, which
writes the FASTA files of the final clusters to `fasta/mags_with-mges/` and
their plots to `plots/`. Without `--mges`, the run ends with the validation,
and says so in its log. With `--mges`, the list of the MGE contigs of the
assembly, it goes on with two more steps on the map of the clusters, which it
updates in place, as the commands do one by one in its output directory:
`mge`, as [`mge --reassign`](#mge), with `--mge-threshold-bin`, clusters the
MGE contigs, makes each MGE cluster a genome of the map and writes
`fasta/mges/`; `host`, as [`host`](#host), with `--mge-min-share` and
`--mge-min-contigs`, finds the MAGs that host the MGE clusters and writes
`tables/mge_hosts.txt` and `plots/mge_hosts.pdf`. They are steps 4 and 5 of a
run from a map, 5 and 6 of a run from pairs or BAM files. A run that stops in
one of them leaves the map of the clusters of the validation (with its MGE
clusters, in `host`), and its log says what to run on it to complete the map:
`mge --reassign` then `host`, or `host`. What they would refuse is refused
before any work, the reason logged: a list that cannot be read, a map without
the contact statistics of the reads of its contigs at 1001 bp (see
[`mge`](#mge); the map that the run makes from pairs or BAM files has them), a
list none of whose contigs is in the map (from pairs or BAM files, in the
assembly, before the run makes its map), and the MGE files that an earlier
version wrote in the output directory.

At its end, the run writes its outputs from its map (see
[Outputs](#conventions)): `fasta/mags/`, each MAG (its genomes but the MGE
genomes, which took their contigs out of the MAGs with `--mges`) as its
contigs and stitched, with `metacooler fasta`; `cool/`, a cool file per MAG,
with the MGE clusters that it hosts, extracted from the map by
`metacooler extract` and balanced as `cooler balance` balances it, the whole of
it at once, with its defaults (a cool file whose every bin is masked, too
small for `cooler balance`, is a warning); and `tables/`, the tables of
[`dump`](#dump).

| Option | Default | Description |
|---|---|---|
| `INPUT...` | *required* | A `.metacool` contact map of the contigs (`MAP`), or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`): the run reads its network from it, and writes the map of its clusters, with the network of its contigs, as `network` stores it, to the output directory, under the name of `MAP` (see [The map of a run](#conventions)). See [Starting from a contact map](guides/metacool.md#starting-from-a-contact-map). Or the pairs files or the BAM files of the run, which it makes the map of the contigs from, in its output directory, then updates in place: see [The inputs of a run](#inputs). |
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip; a compressed one is decompressed once, into the temporary directory of the run, for the FASTA files of the clusters. |
| `-d, --depth` | none | The coverage of the contigs, from `jgi_summarize_bam_contig_depths` of MetaBAT2. The `Shotgun_coverage` column of the map gives it without `--depth`, which replaces it: see [Starting from a contact map](guides/metacool.md#starting-from-a-contact-map). |
| `-e, --enzyme` | none | The restriction enzymes of the library. The `Restriction_site` column of the map gives the restriction sites without `--enzyme`, which replaces it. |
| `-m, --normalization` | `empirical_hit` | `None`, `abundance`, `length`, `RS`, `empirical_hit` or `theoritical_hit`. `RS` and `theoritical_hit` need `--enzyme`; `abundance` and `theoritical_hit` need `--depth`, or the columns of the map that give them. See [Normalising the network](guides/network.md#normalising-the-network). |
| `-q, --min-quality` | `30` | BAM files: the mapping quality that both mates of a read must reach for its pair to count. Refused with a map or pairs files, whose contacts are read as they are. |
| `-A, --algorithm` | `louvain` | The partition algorithm: `louvain` or `leiden`. |
| `-i, --iterations` | `100` | The iterations of the partition. |
| `-O, --overlap` | `80` | The percentage of iterations in which two core clusters must fall together to be merged into one overlapping cluster. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, that an overlapping cluster must reach to be checked with miComplete, and that a cluster of the recursive partition must exceed. |
| `-j, --rec-iter` | `10` | The iterations of the recursive partition of the contaminated clusters. |
| `-P, --rec-overlap` | `90` | The overlap threshold of the recursive partition, in percent. |
| `--seed` | none | The seed of the random numbers of Louvain or Leiden, set once before the partition and again before the recursive partition of the validation, their iterations drawing other numbers one after the other: from the same map, the same seed gives the same clusters, those of `partition --seed S` then `validation --seed S`. Without it, each run draws other random numbers (see the random numbers, in [Conventions](#conventions)). |
| `--mges` | none | The MGE contigs, one name per line, as `mge --mges` reads them: after the validation, the run goes on with two more steps, `mge` (as `mge --reassign`) and `host`. Without it, the run ends with the validation. |
| `--mge-threshold-bin` | `0.8` | The step `mge`: the score above which two MGE contigs, or two MGE clusters, are merged, as `--threshold-bin` of `mge`. |
| `--mge-min-share` | `10` | The step `host`: the share of the signal of an MGE cluster with the MAGs, in percent, from which a MAG is a host, as `--min-share` of `host`. |
| `--mge-min-contigs` | `5` | The step `host`: the number of contigs of a MAG in contact with an MGE cluster from which the MAG is a host, as `--min-contigs` of `host`. |
| `-o, --outdir` | the current directory | The output directory: the map of the clusters, the log, and the outputs of the run: `fasta/` (`mags_with-mges/`, `mags/`, and `mges/` with `--mges`), `cool/`, `tables/` and `plots/` (see [Outputs](#conventions)). |
| `-p, --prefix` | none | Name the clusters `<prefix>_metator_…`, and the map that the run makes from pairs or BAM files `<prefix>.metacool` (by default, named after the assembly). Letters, digits, `_` and `-` only. |
| `-F, --force` | off | Replace the map, the folders that the run makes anew (`fasta/mags_with-mges/`, `fasta/mags/`, `cool/`, and `fasta/mges/` with `--mges`), the tables of `tables/` and the plots of the clusters (in `plots/`, whose other files are left as they are) of an earlier run in the output directory. `tables/mge_hosts.txt` and `plots/mge_hosts.pdf` are replaced without it, as `host` replaces them. The other results of an earlier run in the output directory (the tables, `final_clusters/`, `plot/`, `mges_binned.fa` and `mge_hosts.txt` of a development build of 2.0, the outputs of 1.x) are refused, with or without `-F`, and never removed. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run works in a directory of its own: the FASTA files of the clusters that miComplete checks (at most the size of the assembly), miComplete's results, the tables of the validation, and, from a compressed assembly, its decompressed copy; from BAM files, their pairs too (about 25 bytes per pair, plus the names of its two contigs). Not in a folder that the run makes anew (`fasta/mags_with-mges/`...), nor in a directory of the results of an earlier run (`recursive_clusters/`, `final_clusters/`...), which the run refuses. |
| `-t, --threads` | `1` | Threads for the partition and the validation (miComplete), and to read the BAM files and bin the pairs that the run makes its map from. |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run, which the log names: the tables of the validation (`cluster_summary.txt`, `contig_data_final.txt`, `binning.txt`, `clusters.txt`, `cluster_members.txt`), miComplete's results and the FASTA files of the clusters that it checked. |
| `--inplace` | off | Update `MAP` in place with the map of the clusters, rather than write it to the output directory: refused for an entry of a collection, and for pairs or BAM files, whose map the run makes in its output directory and updates in place. A run whose output directory holds `MAP` already updates it in place without `--inplace`. |

The validation is the step that needs the most memory: miComplete runs
Prodigal and HMMER (pyrodigal and pyhmmer) on every cluster.

```bash
# The map of the contigs, from the reads, with hicstuff
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz

# The clusters: binning/sample.metacool, the map of the clusters; hicstuff/sample.metacool is left as it was
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -d depth.txt -t 16 -o binning

# The same, updating hicstuff/sample.metacool in place
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -o binning --inplace

# From the pairs files of two libraries: binning/sample.metacool, made from them, becomes the map of the clusters
metator pipeline 'hicstuff_*/*.valid_idx.pairs.gz' -a assembly.fa -e DpnII -d depth.txt -p sample -o binning

# From BAM files: of both mates (bwa mem -5SP), one per library; or of one mate each, the R1 files then the R2 files
metator pipeline L1.bam,L2.bam -a assembly.fa -e DpnII -o binning
metator pipeline 'L*_R1.bam' 'L*_R2.bam' -a assembly.fa -e DpnII -o binning

# The clusters, then the MGE clusters and their hosts: binning/fasta/mges/,
# binning/tables/mge_hosts.txt, and the map of the clusters with both
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -t 16 -o binning --mges mges.txt
```

### `network` { #network }

Build the network of contigs from the map, and store it in the map: the first
step of the pipeline, whose network `partition` and `validation` read from the
map. **Guide:** [Building the network](guides/network.md).

```console
metator network [OPTIONS] INPUT...
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
place with `--inplace` (see [The map of a run](#conventions)); from pairs or
BAM files, it makes the map of the contigs first, in its output directory,
named after the assembly, then stores the network in it (see
[The inputs of a run](#inputs)). It writes no
table in the output directory: `network.txt` is gone, and
`contig_data_network.txt`, the contig
data that it stores in the map, is written in its temporary directory, which
`-N` keeps. [`metator dump`](#dump) writes the network, `network.txt`, from
the map of the clusters of a run; metacooler writes it from any map that holds
it, one that only `network` ran on included (see
[The network](formats.md#the-network)).

| Option | Default | Description |
|---|---|---|
| `INPUT...` | *required* | As in `pipeline`: a map, which the run writes, with the network of its contigs, to the output directory; or the pairs files or the BAM files of the run, which it makes the map of the contigs from, `<assembly>.metacool`, in its output directory. |
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip. |
| `-d, --depth` | none | As in `pipeline`. |
| `-e, --enzyme` | none | As in `pipeline`. |
| `-n, --normalization` | `empirical_hit` | As `--normalization` (`-m`) in `pipeline`. |
| `-q, --min-quality` | `30` | As in `pipeline`: BAM files only. |
| `-o, --outdir` | the current directory | The output directory, of the map of the run and of its log. |
| `-F, --force` | off | Replace the map of an earlier run in the output directory: that of `MAP`, or the map that the run makes from pairs or BAM files. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run writes, in a directory of its own, the contig data that it stores in the map (`contig_data_network.txt`), and the pairs of the BAM files that it makes its map from. |
| `-t, --threads` | `1` | Threads to read the BAM files and bin the pairs that the run makes its map from. |
| `-N, --no-clean-up` | off | Keep the temporary directory of the run, which the log names, and its `contig_data_network.txt`. |
| `--inplace` | off | Store the network in `MAP` itself, in place, rather than in the map of the output directory: refused for an entry of a collection, and for pairs or BAM files. |

```bash
# binning/sample.metacool, with the network of the contigs
metator network sample.metacool -a assembly.fa -e DpnII -o binning

# binning/assembly.metacool, the map of the contigs made from the pairs files, with their network
metator network L1.pairs.gz L2.pairs.gz -a assembly.fa -e DpnII -o binning
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
| `--seed` | none | The seed of the random numbers of Louvain or Leiden, set once before the iterations, which draw other numbers one after the other: from the same map, the same seed gives the same partition. |
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
again: the third step. Writes the FASTA files of the final clusters
(`fasta/mags_with-mges/`: each final cluster, the MGE contigs that it holds
included) and their plots (`plots/`, as PDF files), and stores the final
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

Once the clusters are final, the run draws two plots of them in `plots/`, in
its output directory, as PDF files, from miComplete's numbers:
`bins_distribution.pdf`, the
completeness and the contamination (weighted redundancy − 1, drawn from 0 to
105%) of each final cluster, the clusters by decreasing completeness; and
`bins_size_distribution.pdf`, a pie chart of the size of the clusters of each
quality class (HQ, MQ and LQ MAGs, contaminated clusters and the others, the
classes that the log counts: see [Cluster quality](concepts.md#cluster-quality))
and of the contigs in no cluster, of the contigs of the map in all. A plot
that cannot be drawn is a warning, and the run goes on. The plots of an
earlier run are refused without `-F`, which replaces them; the other files of
`plots/` are left as they are. The `plot/` of an earlier version (its PNG
files, those of 1.x that 2.0 no longer draws) is refused, with or without
`-F`, as the other outputs of an earlier version.

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
| `--seed` | none | The seed of the random numbers of the recursive partition, set once before it, its iterations drawing other numbers one after the other: from the same map, the same seed gives the same clusters, and `partition --seed S` then `validation --seed S` give those of `pipeline --seed S`. |
| `-o, --outdir` | `./` | The output directory: the map of the clusters, `fasta/mags_with-mges/`, `plots/` and the log. |
| `-p, --prefix` | none | Name the recursive clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace the map, `fasta/mags_with-mges/` and the plots (in `plots/`, whose other files are left as they are) of an earlier run in the output directory. The other results of an earlier run in the output directory (the tables, `recursive_clusters/`, `final_clusters/` and `plot/` of a development build of 2.0, the outputs of 1.x) are refused, with or without `-F`, and never removed. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run writes, in a directory of its own, the clusters that miComplete checks (at most the size of the assembly), miComplete's results and the tables of the validation, and, from a compressed assembly, its decompressed copy. Not in a folder that the run makes anew (`fasta/mags_with-mges/`...), nor in a directory of the results of an earlier run (`recursive_clusters/`, `final_clusters/`...), which the run refuses. |
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
and store the MGE clusters in the contact map. [`pipeline --mges`](#pipeline)
runs it, with `--reassign`, as its step `mge`, on the map of its clusters.
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

The run writes the FASTA files of its MGE clusters to `fasta/mges/`, two for
each, named after it: `MetaTOR_MGE_<ID>.fa`, its contigs as records of their
own, and `MetaTOR_MGE_<ID>.stitched.fa`, one record, its contigs joined by 200
`N`, both the longest first, ties by name; and it stores its MGE clusters in
the map: as
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
clusters, or their FASTA file, that an earlier version wrote in the output
directory (`mges_bin_summary.tsv`, `mges_binned.fa`, or `mges_cluster_summary.tsv`
of a development build), with or without `-F`: it is never removed. The
`fasta/mges/` of an earlier run is refused without `-F`, which replaces it.
`--random` stores the round `mge_random` and writes `fasta/mges_random/`, and
leaves the round `mge` and `fasta/mges/` as they were. A run that
replaces the clusters of the map (`partition`, `validation`, `pipeline`)
drops the rounds `mge` and `mge_random`, and their entries, with a warning:
run `mge` again on its map (`validation` and `pipeline` make their final
clusters the genomes of the map, without MGE genomes); `pipeline --mges`
stores its own round `mge` there in its step `mge`, and warns only of the round
`mge_random`. A run that replaces the round `mge` (not `--random`) drops the
hosts that [`metator host`](#host) stored in the map, with a warning: run
`host` again.

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip. |
| `-m, --mges` | *required* | The MGE contigs, one name per line. |
| `-s, --threshold-bin` | `0.8` | The score above which two contigs, or clusters, are merged. |
| `-r, --random` | off | Shuffle the contigs between the clusters: a random clustering of the same shape, stored as the round `mge_random` and written to `fasta/mges_random/`. |
| `--reassign` | off | Make each MGE cluster a genome of the map, with its GC content, hits and clustering score: its contigs leave the final clusters that held them, which keep their values. The column `mge` of the genomes table says which genomes are MGE genomes. Refused with `--random`. |
| `--seed` | none | The seed of the random clustering of `--random`: from the same map, the same seed gives the same random clustering. |
| `-o, --outdir` | the current directory | The output directory: the map of the run and `fasta/mges/` (`fasta/mges_random/` with `--random`). |
| `-F, --force` | off | Replace the map and `fasta/mges/` (`fasta/mges_random/` with `--random`) of an earlier run in the output directory. A summary or a FASTA file of the MGE clusters that an earlier version wrote in the output directory (`mges_bin_summary.tsv`, `mges_cluster_summary.tsv`, `mges_binned.fa`, `mges_binned_random.fa`) is refused, with or without `-F`, and never removed. |
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
between them, and store them in the map. [`pipeline --mges`](#pipeline) runs
it as its step `host`, after `mge --reassign`.
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
contacts: `share_of_all`, in `tables/mge_hosts.txt`, says how small.

The run writes `tables/mge_hosts.txt`, every MGE cluster and MAG in contact,
and `plots/mge_hosts.pdf`, the signal and the share of the hosts (see
[Hosts of the MGE clusters](formats.md#hosts-of-the-mge-clusters)), in its
output directory, in place of those of an earlier run (the `mge_hosts.txt`
and `mge_hosts.png` that an earlier version wrote at the top of the output
directory are refused, with or without `-F`), and stores the hosts
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
`--random`; `partition`, `validation` and `pipeline` (but `pipeline --mges`,
which finds them again in its step `host`, without the warning); and `network`,
when it stores another network than the one the hosts were found on (its
`network_digest`).

| Option | Default | Description |
|---|---|---|
| `-s, --min-share` | `10` | The share of the signal of an MGE cluster with the MAGs, in percent, from which a MAG is a host. |
| `-c, --min-contigs` | `5` | The number of contigs of a MAG in contact with an MGE cluster from which the MAG is a host. |
| `-o, --outdir` | the current directory | The output directory: the map of the run, `tables/mge_hosts.txt`, `plots/mge_hosts.pdf` and the log. |
| `-F, --force` | off | Replace the map of an earlier run in the output directory. `tables/mge_hosts.txt` and `plots/mge_hosts.pdf` are replaced without it. The `mge_hosts.txt` and `mge_hosts.png` of an earlier version are refused, with or without `-F`, and never removed. |
| `--inplace` | off | Store the hosts in `MAP` itself, in place, rather than in the map of the output directory: refused for an entry of a collection. |

```bash
# The hosts of the MGE clusters that mge stored in the map of a run, in
# binning: updated in place
metator host binning/sample.metacool -o binning
```

### `export` { #export }

Write the outputs of a run again from the map of its clusters: the folders
that `pipeline` writes in its output directory (see
[Outputs](#conventions)), from the map alone, and the assembly for the FASTA
files. **API:** `metator.outputs`.

```console
metator export [OPTIONS] MAP
```

The map is the one that `pipeline` wrote, or `validation`, then `mge` and
`host`: it holds every result of the run, so `export` writes its outputs
again, in another output directory or in place of those of the run, without
running it again. It writes:

| Folder | From the map | Needs |
|---|---|---|
| `fasta/mags_with-mges/` | The FASTA file of each final cluster (the round `final`), with the MGE contigs that it holds: the records of the validation's files, in the order of the map, where the validation wrote them in the order of its partition. | `--assembly` |
| `fasta/mags/` | Each MAG (a genome of the map but the MGE genomes), as its contigs and stitched, as `pipeline` writes it. | `--assembly` |
| `fasta/mges/`, `fasta/mges_random/` | Each MGE cluster of the round `mge` (`mge_random`), its contigs the longest first, as `mge` writes it. | `--assembly`, and the round in the map |
| `cool/` | A cool file per MAG, with the MGE clusters that it hosts, balanced, as `pipeline` writes it. | |
| `tables/` | The tables of [`dump`](#dump); and `mge_hosts.txt`, the hosts of the MGE clusters, found again from the network of the map with the thresholds that the host run recorded (`--min-share` and `--min-contigs` of `host`, or `--mge-min-share` and `--mge-min-contigs` of `pipeline --mges`). | `mge_hosts.txt`: hosts in the map |
| `plots/` | The plots of the final clusters, from the summary that the map holds, rounded to 2 decimals as in `cluster_summary.txt` (the validation drew the values of miComplete: a cluster within 0.005 of a threshold of quality can fall in the next class); and `mge_hosts.pdf`. | `mge_hosts.pdf`: hosts in the map |

Without `--assembly`, the FASTA files are left out, and the log names the
assembly that the run recorded. The map is read only, without a lock: a map
that [`dump`](#dump) refuses is refused, before anything is written, and so
is a map that another run writes while `export` reads it (run it again
then). What cannot be rebuilt from the map is not written: the log of the
run, and the temporary files of its validation (miComplete's results, the
FASTA files of the clusters that it checked: `metacooler fasta MAP --round
overlapping` writes those of a round). The folders that the run makes anew,
the tables of `tables/` and the plots of the clusters of an earlier run in the
output directory are refused without `-F`, and replaced with it;
`tables/mge_hosts.txt` and `plots/mge_hosts.pdf` are replaced without it, as
`host` replaces them; the outputs of an earlier version are refused, with or
without `-F`. The map and the assembly may not be in a folder that the run
makes anew.

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The map of the clusters of a run, which `pipeline` or `validation` wrote, and `mge` and `host` updated, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`). Read only. |
| `-a, --assembly` | none | The assembly that the map was made from: a FASTA file, plain or compressed with gzip. With it, the run writes the FASTA files too (`fasta/`). |
| `-o, --outdir` | the current directory | The output directory: `fasta/` (with `--assembly`), `cool/`, `tables/`, `plots/` and the log. |
| `-F, --force` | off | Replace the folders that the run makes anew (`fasta/mags_with-mges/`, `fasta/mags/`, `fasta/mges/`, `fasta/mges_random/`, `cool/`), the tables of `tables/` and the plots of the clusters (in `plots/`, whose other files are left as they are) of an earlier run in the output directory. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the run writes, in a directory of its own, the decompressed copy of a compressed assembly. Not in a folder that the run makes anew, nor in a directory of the results of an earlier version, which the run refuses. |

```bash
# Every output of a run, from the map of its clusters and its assembly
metator export binning/sample.metacool -a assembly.fa -o outputs

# The same, but the FASTA files: cool/, tables/ and plots/
metator export binning/sample.metacool -o outputs

# The outputs of the run in its own output directory, written again
metator export binning/sample.metacool -a assembly.fa -o binning -F
```
