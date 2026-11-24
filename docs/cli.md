# CLI reference

The `metator` command has one subcommand per task. This page documents every
command and option; `-h/--help` prints the same in a terminal.

```console
metator --help
metator <command> --help
metator --version
```

As in `metator --help`, the commands are grouped in
[Pipeline](#pipeline-commands) (`pipeline`), [Steps](#steps) (`network`,
`partition`, `validation`, `markers`, `mge`, `host`) and [Outputs](#outputs)
(`dump`, `export`, `report`).

New to metaTOR? Start with the [Quickstart](quickstart.md) or
[The end-to-end pipeline](end-to-end.md); for the clusters, see
[Concepts](concepts.md). [Other tools](other-tools.md) replaces the commands
that 2.0 removed (`qc`, `contactmap`, `metacool`, `scaffold`, `pairs`) and the
heatmap, GC and coverage plots.

## Conventions

| Convention | Meaning |
|---|---|
| `<name>` | The name of the run's map without `.metacool`: `sample` for `sample.metacool`, `NAME` for `file.mmetacool::entries/NAME`. A map that `network` or `pipeline` makes from pairs or BAM files is `<prefix>.metacool` (`pipeline --prefix`) or `<assembly>.metacool`. `--prefix` of `partition` and `validation` only names the clusters. |
| The map of a run | A run writes its map to its output directory (`-o`) as `<name>.metacool` and leaves `MAP` as it was. `--inplace` updates `MAP` itself, and so does a run whose output directory already holds `MAP`: the steps can run one by one in one directory (`network MAP -o DIR`, then `partition DIR/MAP -o DIR`…). |
| Writing a map | Through metacooler: the run locks the map (`<map>.lock`) and writes a staged copy beside it, renamed into place once validated, so a run that stops leaves the map as it was. A map held by another process, or that cannot be written, is refused before any work. A killed run leaves its lock and copy: the next run on the same host removes them; on another host, delete the files the error names. |
| Several libraries | Give the pairs or BAM files of every library to `network` or `pipeline`, which sum them in one map (see [The inputs of a run](#inputs)), or merge the maps with `metacooler merge` (see [Other tools](other-tools.md)). |
| Several enzymes | Comma-separated: `DpnII,HinfI`. |
| `-F, --force` | Replace what an earlier run of the command wrote in the output directory (map, folders, cluster plots, tables, report); without it, the command refuses. Outputs of earlier versions (1.x, development builds of 2.0: `final_clusters/`, `plot/`, `mges_binned.fa`…) are always refused and never removed. A run that updates its map in place first warns of the clusters it replaces. |
| `-T, --tmpdir`, `-N, --no-clean-up` | Each run works in its own directory inside `-T`, removed at the end unless `-N` is given or the run fails; the log names it. |
| Log files | `<name>.metator_<date>.log`, the date as `%Y%m%d%H%M%S` (`sample.metator_20261123225838.log`), in the output directory (`mge`: in the current directory). It opens with the version, command line and options, and ends with the files written and the run time. `pipeline` numbers its steps: `[1/5] network` … `[4/5] markers`, `[5/5] export`, plus `map` first from pairs or BAM files, and `mge` and `host` before `export` with `--mges`. |
| Random numbers | Louvain and Leiden (`partition`, `pipeline`, the recursive partition of `validation`) and `mge --random` draw new random numbers at each run, so two runs give close but not identical clusters. `--seed S` makes a run reproducible for the same map, options, igraph and metaTOR; `pipeline --seed S` gives the clusters of `partition --seed S` then `validation --seed S`. The map records the seed. |

### Files of a run { #files }

Besides its map and its log, a run writes:

| Path | Written by | Contents |
|---|---|---|
| `fasta/mags_with-mges/` | `validation`, `pipeline` | One FASTA file per final cluster, MGE contigs included. |
| `fasta/mags/` | `pipeline` | Each MAG (the genomes of the map but the MGE genomes) as its contigs, `<MAG>.fa`, and stitched, `<MAG>.stitched.fa` (one record, contigs joined by 200 `N`). |
| `fasta/mges/` | `mge`, `pipeline --mges` | Each MGE cluster, the same two ways (`fasta/mges_random/` for `mge --random`). |
| `cool/` | `pipeline` | One balanced cool file per MAG, with the MGE clusters it hosts. |
| `tables/` | `pipeline`; `host` | The tables of [`dump`](#dump), and `mge_hosts.txt` (`host`, `pipeline --mges`). |
| `plots/` | `validation`, `host`, `pipeline` | PDF plots of the clusters and of the hosts. |
| `<name>.report.txt`, `<name>.report.html`, `<name>.report_mqc.yaml` | `pipeline` | The report of the run. |

[`export`](#export) writes them all again from the map of a run;
[`report`](#report) writes the report alone.

## The inputs of a run { #inputs }

`network` and `pipeline` take `INPUT...`: a contact map of the contigs, or the
pairs or BAM files to make it from.

| Input | What the run does |
|---|---|
| A `.metacool` map, or an entry `file.mmetacool::entries/NAME` | Reads the network from it and writes the map of the run to the output directory (see [Conventions](#conventions)). |
| Pairs files, one per library: plain, or compressed with gzip, bzip2, zstd or LZ4 | Makes the map first, the libraries summed, as `metacooler create pairs` does: one bin per contig, with the read contact statistics that `mge` needs. |
| BAM or SAM files: both mates in each (as `bwa mem -5SP` writes them), one per library; or one mate each (as bowtie2 aligns them), an R1 and an R2 file per library | Takes their pairs, then makes the map as from pairs files. A pair counts when both primary alignments reach `--min-quality` (30) on contigs of the assembly and neither is flagged as a duplicate (0x400) or as failing QC (0x200). Each read is placed at its 5' end, as pairtools and hicstuff's 5spm mode do. |

- **Arguments**: each is a path, a comma-separated list or a glob pattern
  (its files sorted by name). `a.pairs b.pairs`, `a.pairs,b.pairs` and
  `'*.pairs'` are the same two libraries.
- **File type** comes from the first bytes, not the name: HDF5 is a map,
  BAM/SAM is alignments, anything else is pairs. FASTQ, FASTA and CRAM files
  are refused, as are mixed inputs (a map with other files, pairs with BAM,
  both-mate with one-mate BAM files) and a file given twice.
- **One-mate BAM files**: give the R1 files, then the R2 files
  (`'*_R1.bam' '*_R2.bam'`, or `*_R1.bam *_R2.bam` expanded by the shell). The
  i-th R1 and R2 files are one library, matched by read name. Files that are
  neither in the aligner's order nor sorted by name are sorted by name first,
  in the temporary directory.
- **Both-mate BAM files** sorted by coordinate are grouped by read name first
  (`samtools collate`).
- **Sequences** (`@SQ` lines of BAM files, `#chromsize` lines of pairs files)
  must match the contigs of the assembly: a sequence of another length is
  refused, and one that is no contig (a spike-in) is dropped with a warning,
  with its pairs. Inputs that give no contact on the assembly are refused.
- **The map the run makes** is `<prefix>.metacool` (`pipeline --prefix`) or
  `<assembly>.metacool` (`assembly.fa.gz` gives `assembly.metacool`) in the
  output directory; the run then updates it in place. An existing map there is
  replaced only with `-F`; `--inplace` is refused. If the run stops later, a
  new run can start from that map without reading the pairs again.
- **The pairs of BAM files** go to the temporary directory, one file per
  library (about 25 bytes per pair, plus the contig names; `-N` keeps them).
  `--threads` reads the BAM files and bins the pairs in parallel.

## Pipeline { #pipeline-commands }

### `pipeline` { #pipeline }

Run `network`, `partition`, `validation` and `markers` in one command; with
`--mges`, then `mge` and `host`; last, `export`. **Guide:**
[Running the pipeline](guides/pipeline.md).

```console
metator pipeline [OPTIONS] INPUT...
```

- **Input**: a `.metacool` map of the contigs, made by
  `hicstuff pipeline -M metacool` or `metacooler create pairs` (see
  [Other tools](other-tools.md#from-the-reads-to-the-map)), or pairs or BAM
  files, from which the run makes that map first (see
  [The inputs of a run](#inputs)).
- **Validation**: writes the map of the clusters, `fasta/mags_with-mges/` and
  `plots/`. It is the step that needs the most memory: miComplete runs
  Prodigal and HMMER (pyrodigal, pyhmmer) on every cluster.
- **Markers**: as [`markers`](#markers), finds miComplete's marker genes on
  each contig, in its final cluster, and stores them in the map of the
  clusters, with miComplete's default thresholds.
- **With `--mges`**: the step `mge` (as [`mge --reassign`](#mge)) clusters the
  MGE contigs, makes each MGE cluster a genome of the map and writes
  `fasta/mges/`; the step `host` (as [`host`](#host)) finds their host MAGs
  and writes `tables/mge_hosts.txt` and `plots/mge_hosts.pdf`. What they would
  refuse (an unreadable list, a list with no contig in the map, a map without
  the read contact statistics that `mge` needs) is refused before any work. If
  the run stops in one of them, the log says what to run to finish the map.
- **Export**: as [`export`](#export), writes all the outputs from the map,
  replacing those of the steps, and the report last (see
  [Files of a run](#files)). The cool files are balanced with the defaults of
  `cooler balance`; one too small to balance is a warning.

| Option | Default | Description |
|---|---|---|
| `INPUT...` | *required* | A map, an entry of a collection, or pairs or BAM files (see [The inputs of a run](#inputs) and [Starting from a contact map](guides/metacool.md#starting-from-a-contact-map)). |
| `-a, --assembly` | *required* | The assembly FASTA file, plain or gzip-compressed. |
| `-d, --depth` | none | Contig coverage from MetaBAT2's `jgi_summarize_bam_contig_depths`; replaces the map's `Shotgun_coverage` column. |
| `-e, --enzyme` | none | The restriction enzymes; replaces the map's `Restriction_site` column. |
| `-m, --normalization` | `empirical_hit` | `None`, `abundance`, `length`, `RS`, `empirical_hit` or `theoritical_hit` (see [Normalising the network](guides/network.md#normalising-the-network)); `RS` and `theoritical_hit` need restriction sites, `abundance` and `theoritical_hit` coverage. |
| `-q, --min-quality` | `30` | BAM files only: the mapping quality both mates must reach. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `100` | Iterations of the partition. |
| `-O, --overlap` | `80` | Percentage of iterations in which two core clusters must fall together to merge into one overlapping cluster. |
| `-r, --res-param` | `1.0` | Resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | Size (bp) from which an overlapping cluster is checked, and that a recursive cluster must exceed. |
| `-j, --rec-iter` | `10` | Iterations of the recursive partition. |
| `-P, --rec-overlap` | `90` | Overlap threshold (%) of the recursive partition. |
| `--seed` | none | Seed of the partition and the recursive partition (see [Conventions](#conventions)). |
| `--mges` | none | The MGE contigs, one name per line: adds the steps `mge` and `host`. |
| `--mge-threshold-bin` | `0.8` | `--threshold-bin` of the step `mge`. |
| `--mge-min-share` | `10` | `--min-share` of the step `host`. |
| `--mge-min-contigs` | `5` | `--min-contigs` of the step `host`. |
| `-o, --outdir` | the current directory | The output directory. |
| `-p, --prefix` | none | Name the clusters `<prefix>_metator_…`, and a map made from pairs or BAM files `<prefix>.metacool`; letters, digits, `_` and `-` only. |
| `-F, --force` | off | Replace the map, outputs and report of an earlier run (`tables/mge_hosts.txt` and `plots/mge_hosts.pdf` are replaced without it). |
| `-T, --tmpdir` | `./tmp` | Temporary directory: the clusters that miComplete checks (up to the size of the assembly), its results, the validation tables, the table of the markers, a decompressed assembly and the pairs of BAM files. |
| `-t, --threads` | `1` | Threads for the partition, miComplete, the markers and reading BAM files. |
| `-N, --no-clean-up` | off | Keep the temporary directory: miComplete's results, the validation tables and the table of the markers. |
| `--inplace` | off | Update `MAP` in place; refused for a collection entry and for pairs or BAM files. |

```bash
# The map of the contigs, from the reads, with hicstuff
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz

# The clusters in binning/sample.metacool; hicstuff/sample.metacool is left as it was
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -d depth.txt -t 16 -o binning

# The same, updating hicstuff/sample.metacool in place
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -o binning --inplace

# From the pairs files of two libraries: makes binning/sample.metacool
metator pipeline 'hicstuff_*/*.valid_idx.pairs.gz' -a assembly.fa -e DpnII -d depth.txt -p sample -o binning

# From BAM files: both mates in each; or one mate each, the R1 files then the R2 files
metator pipeline L1.bam,L2.bam -a assembly.fa -e DpnII -o binning
metator pipeline 'L*_R1.bam' 'L*_R2.bam' -a assembly.fa -e DpnII -o binning

# The clusters, then the MGE clusters and their hosts
metator pipeline hicstuff/sample.metacool -a assembly.fa -e DpnII -t 16 -o binning --mges mges.txt
```

## Steps

### `network` { #network }

Build the network of the contigs and store it in the map: the first step.
**Guide:** [Building the network](guides/network.md).

```console
metator network [OPTIONS] INPUT...
```

The run stores in the map:

- the contig data, in the contigs table: `GC_content`, `Hit`, and
  `Shotgun_coverage` and `Restriction_site` when known;
- the network, as the bin weights `metator_<normalization>` (none for
  `None`), replacing an earlier network;
- its command line, under `metator` → `network` in the metadata.

The clusters of the map are kept. If the new network differs from the one
they were made on, the run warns: `validation` refuses the map until
`partition` runs again, and [`dump`](#dump) refuses it too. The run writes no
table: `dump` writes `network.txt` from the map of the clusters, and
metacooler from any map that holds a network (see
[The network](formats.md#the-network)).

| Option | Default | Description |
|---|---|---|
| `INPUT...` | *required* | As in `pipeline`: a map, or pairs or BAM files (see [The inputs of a run](#inputs)). |
| `-a, --assembly` | *required* | The assembly FASTA file, plain or gzip-compressed. |
| `-d, --depth` | none | As in `pipeline`. |
| `-e, --enzyme` | none | As in `pipeline`. |
| `-n, --normalization` | `empirical_hit` | As `--normalization` (`-m`) in `pipeline`. |
| `-q, --min-quality` | `30` | As in `pipeline`: BAM files only. |
| `-o, --outdir` | the current directory | Output directory of the map and the log. |
| `-F, --force` | off | Replace the map of an earlier run. |
| `-T, --tmpdir` | `./tmp` | Temporary directory: the contig data (`contig_data_network.txt`) and the pairs of BAM files. |
| `-t, --threads` | `1` | Threads to read BAM files and bin their pairs. |
| `-N, --no-clean-up` | off | Keep the temporary directory. |
| `--inplace` | off | Store the network in `MAP` itself; refused for a collection entry and for pairs or BAM files. |

```bash
# binning/sample.metacool, with the network of the contigs
metator network sample.metacool -a assembly.fa -e DpnII -o binning

# binning/assembly.metacool, made from the pairs files, with its network
metator network L1.pairs.gz L2.pairs.gz -a assembly.fa -e DpnII -o binning
```

### `partition` { #partition }

Partition the network many times into core and overlapping clusters: the
second step. **Guide:** [Partition and validation](guides/binning.md).

```console
metator partition [OPTIONS] MAP
```

The run reads the network that `network` stored (a map without one is
refused) and stores in the map:

- every core cluster (`core_00001`…) and overlapping cluster
  (`<prefix>_metator_00001_00000`…), whatever its size, as the clustering
  rounds `core` and `overlapping`, replacing earlier ones;
- its parameters under `metator` → `partition` in the metadata, with the
  digest of the network it partitioned (`network_digest`).

It writes no table or FASTA file, and leaves the genomes and contigs tables as
they were. `validation` checks the overlapping clusters that reach its
`--size`. See [Step by step](guides/metacool.md#step-by-step).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | A map with the network of `network`, or an entry of a collection. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `100` | Iterations of the partition. |
| `-O, --overlap` | `80` | Overlap threshold (%). |
| `-r, --res-param` | `1.0` | Resolution of the Leiden algorithm. |
| `--seed` | none | Random seed: the same map and seed give the same partition. |
| `-o, --outdir` | `./` | Output directory of the map and the log. |
| `-p, --prefix` | none | Name the overlapping clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace the map of an earlier run. |
| `-T, --tmpdir` | `./tmp` | Temporary directory; `partition` writes nothing there. |
| `-t, --threads` | `1` | Threads for the partition. |
| `-N, --no-clean-up` | off | Keep the temporary directory; `partition` leaves nothing there. |
| `--inplace` | off | Store the partition in `MAP` itself; refused for a collection entry. |

```bash
# The map that metator network wrote in binning, updated in place
metator partition binning/sample.metacool -t 16 -o binning
```

### `validation` { #validation }

Check the clusters with miComplete and partition the contaminated ones again:
the third step. **Guide:** [Partition and validation](guides/binning.md#validation).

```console
metator validation [OPTIONS] MAP
```

- **Checked clusters**: the overlapping clusters of the partition that reach
  `--size`. A cluster is partitioned again when its weighted completeness
  reaches 33% and its contamination, (weighted redundancy − 1) / weighted
  completeness, reaches 5%. Its parts replace it when the most complete part
  keeps enough of its completeness.
- **Writes** `fasta/mags_with-mges/`, one FASTA file per final cluster (MGE
  contigs included), and two plots in `plots/`: `bins_distribution.pdf`, the
  completeness and contamination of each final cluster, and
  `bins_size_distribution.pdf`, the size of the clusters by quality class (see
  [Cluster quality](concepts.md#cluster-quality)). A plot that fails is a
  warning.
- **Stores** in the map the final clusters, their summary, every clustering
  round and the normalization of the network, as `pipeline` does.
- **Keeps** the validation tables (`cluster_summary.txt`,
  `contig_data_final.txt`, `binning.txt`, `clusters.txt`,
  `cluster_members.txt`) and `recursive_clusters/` in its temporary directory
  only; [`dump`](#dump) writes the tables from the map.
- **Refuses** a map without the network or the partition, or whose partition
  was made on another network (run `partition` again).

See [Step by step](guides/metacool.md#step-by-step).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | A map with the network and partition of `network` and `partition`, or an entry of a collection. |
| `-a, --assembly` | *required* | The assembly FASTA file, plain or gzip-compressed. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `10` | Iterations of the recursive partition. |
| `-O, --overlap` | `90` | Overlap threshold (%) of the recursive partition. |
| `-r, --res-param` | `1.0` | Resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | Size (bp) from which an overlapping cluster is checked, and that a recursive cluster must exceed. |
| `--seed` | none | Seed of the recursive partition: `partition --seed S` then `validation --seed S` give the clusters of `pipeline --seed S`. |
| `-o, --outdir` | `./` | Output directory: the map, `fasta/mags_with-mges/`, `plots/` and the log. |
| `-p, --prefix` | none | Name the recursive clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace the map, `fasta/mags_with-mges/` and the cluster plots of an earlier run. |
| `-T, --tmpdir` | `./tmp` | Temporary directory: the clusters that miComplete checks (up to the size of the assembly), its results, the validation tables and a decompressed assembly. |
| `-t, --threads` | `1` | Threads for the partition and miComplete. |
| `-N, --no-clean-up` | off | Keep the temporary directory: miComplete's results and the validation tables. |
| `--inplace` | off | Write the map of the clusters over `MAP`; refused for a collection entry. |

```bash
# The map of a step-by-step run in binning becomes the map of the clusters
metator validation binning/sample.metacool -a assembly.fa -t 16 -o binning
```

### `markers` { #markers }

Find miComplete's marker genes on each contig of the map, as miComplete finds
them in a bin, and store them in the map: summed over any group of contigs,
they give the completeness, redundancy, weighted completeness and weighted
redundancy that miComplete gives a bin of those contigs, with no gene to find
again. [`pipeline`](#pipeline) runs it after the validation. **API:**
`metator.markers.contig_markers`, `metator.markers.marker_quality`,
`metator.metacool.write_map`.

```console
metator markers [OPTIONS] MAP
```

- **Genomes as bins.** Each genome of the map (the final clusters of
  `validation` or `pipeline`, or another binner's bins set with
  `metacooler reassign`; not `__unbinned__`) is a bin: pyrodigal is trained on
  its contigs in map order, then finds their genes, and pyhmmer searches the
  marker HMMs of both sets (`Bact105`, 105 bacterial markers; `Arch131`, 131
  archaeal ones) in their proteins, E-values counted over the proteins of the
  genome. The hits are read as miComplete reads them, with its thresholds.
  Summed over a genome, the markers give exactly what miComplete gives its
  FASTA file in map order (`fasta/mags/`, `metacooler fasta`).
- **Other contigs.** miComplete finds no gene in a bin under 20 kb (with the
  12-bp linkers between its contigs), too little to train on. The contigs of
  such a genome, and those of `__unbinned__`, get their genes in metagenome
  mode, contig by contig, and their markers counted over the median number
  of proteins of the trained genomes: as in a typical genome of the map
  (on a map without one, over their own proteins, as one bin).
- **Stores** the markers of each contig in the contigs table, one column per
  set (`markers_Bact105`, `markers_Arch131`: text, `Ribosomal_L22,SecY*2`;
  `*0` for a marker that miComplete sees on the contig but drops, which counts
  in the weighted completeness only; empty for none), and the run under
  `metator` → `markers` in the metadata, with how the genes were found and
  the markers and weights of both sets (see
  [Marker genes of the contigs](formats.md#marker-genes-of-the-contigs)).
  They replace those of an earlier run, and so does `pipeline`, whose step
  `markers` finds them again in its clusters.
- **Every other run keeps them**, those that regroup the contigs too
  (`validation`, `mge --reassign`, `metacooler reassign`): they are facts
  about the contigs. Their genes were found in the genomes as the run read
  them: run `markers` again for the numbers that miComplete gives the genomes
  as they are now.
- **Refuses** an assembly that lacks a contig of the map, gives one another
  length, or names a contig by several records, before any work.
- **On a map of many contigs**, the two columns cost time and memory each
  time cooler opens the map, or metacooler reads its contigs table: about
  0.3 s and 150 MB on a map of 20,566 contigs, in proportion to their number,
  and `export` reads the map once per MAG (see
  [Marker genes of the contigs](formats.md#marker-genes-of-the-contigs)).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | A map whose genomes are bins, or an entry of a collection. |
| `-a, --assembly` | *required* | The assembly FASTA file, plain or gzip-compressed. |
| `--evalue` | `4e-10` | E-value threshold of a hit, counted over the proteins of the genome: miComplete's `--evalue`. |
| `--bias` | `0.3` | A hit whose score times this is no more than its bias is dropped: miComplete's `--bias`. |
| `--domain-cutoff` | `1e-05` | A hit whose best domain's E-value exceeds the protein's by more than this is dropped: miComplete's `--domain-cutoff`. |
| `-o, --outdir` | the current directory | Output directory of the map and the log. |
| `-F, --force` | off | Replace the map of an earlier run. |
| `-T, --tmpdir` | `./tmp` | Temporary directory: the table of the markers (`contig_markers.txt`). |
| `-t, --threads` | `1` | Threads: genomes, and groups of the other contigs, searched side by side. |
| `-N, --no-clean-up` | off | Keep the temporary directory, and `contig_markers.txt`. |
| `--inplace` | off | Store the markers in `MAP` itself; refused for a collection entry. |

```bash
# The marker genes of the contigs of the map of a run, in binning, updated in place
metator markers binning/sample.metacool -a assembly.fa -t 16 -o binning
```

### `mge` { #mge }

Cluster the contigs of mobile genetic elements (MGEs) from their contacts and
store the MGE clusters in the map. [`pipeline --mges`](#pipeline) runs it with
`--reassign`. **API:** `metator.mge.mge_binning`,
`metator.mge.write_mges_fasta`, `metator.metacool.write_map`. **Guide:**
[Clustering mobile genetic elements](guides/mge.md).

```console
metator mge [OPTIONS] MAP
```

- **The map** needs the network of `network`, and metacooler's contact
  statistics of the reads, which give the contacts within each contig (its QC
  block, made from the pairs at 1001 bp by `hicstuff pipeline -M metacool` or
  `metacooler create pairs`, metacooler 0.4.0 or later). A map without them,
  or with none of the MGE contigs, is refused before any work; an MGE contig
  missing from the map is skipped with a warning. The MGE contigs are taken
  longest first, ties by name (see [How it works](guides/mge.md#how-it-works)).
- **Writes** `fasta/mges/`, two files per MGE cluster:
  `MetaTOR_MGE_<ID>.fa`, its contigs, and `MetaTOR_MGE_<ID>.stitched.fa`, one
  record, its contigs joined by 200 `N`.
- **Stores** the round `mge` (contigs, length, GC content, hits and
  clustering score of each cluster), replacing an earlier one, and the run
  under `metator` → `mge` in the metadata.
- **`--reassign`** makes each MGE cluster a genome of the map
  (`MetaTOR_MGE_<ID>`, with `GC`, `Hit` and `BinningScore`). Its contigs leave
  the final clusters that held them, which keep their validation values (the
  log names them); a final cluster left without contigs is no longer a genome.
  The column `mge` of the genomes table is `True` for the MGE genomes.
- **Replacing the round `mge`** first undoes the earlier MGE genomes (their
  contigs go back to their final cluster) and drops the hosts that `host`
  stored, with a warning.
- **`--random`** stores the round `mge_random` and writes
  `fasta/mges_random/`, leaving the round `mge` alone; refused with
  `--reassign`.

`partition`, `validation` and `pipeline` drop the rounds `mge` and
`mge_random` when they replace the clusters of the map, with a warning: run
`mge` again (`pipeline --mges` makes its own round `mge`).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | A map with the network of `network`, or an entry of a collection. |
| `-a, --assembly` | *required* | The assembly FASTA file, plain or gzip-compressed. |
| `-m, --mges` | *required* | The MGE contigs, one name per line. |
| `-s, --threshold-bin` | `0.8` | Score above which two contigs, or clusters, are merged. |
| `-r, --random` | off | Shuffle the contigs between the clusters: a random clustering of the same shape. |
| `--reassign` | off | Make each MGE cluster a genome of the map. |
| `--seed` | none | Seed of the random clustering of `--random`. |
| `-o, --outdir` | the current directory | Output directory: the map and `fasta/mges/` (`fasta/mges_random/`). |
| `-F, --force` | off | Replace the map and `fasta/mges/` (`fasta/mges_random/`) of an earlier run. |
| `-T, --tmpdir` | `./tmp` | Temporary directory, where metacooler decompresses a compressed assembly. |
| `-N, --no-clean-up` | off | Keep the temporary directory; `mge` leaves nothing there. |
| `--inplace` | off | Store the MGE clusters in `MAP` itself; refused for a collection entry. |

```bash
# The map of the clusters of a run, in binning, updated in place
metator mge binning/sample.metacool -a assembly.fa -m mges.txt -o binning
# The same, with the MGE clusters as genomes: binning.txt lists the MAGs and the MGE genomes
metator mge binning/sample.metacool -a assembly.fa -m mges.txt -o binning --reassign
metator dump binning/sample.metacool -o tables
```

### `host` { #host }

Find the MAGs that host the MGE clusters of the map from their contacts, and
store them in the map. [`pipeline --mges`](#pipeline) runs it after
`mge --reassign`. **API:** `metator.host.find_hosts`,
`metator.host.write_hosts`, `metator.metacool.write_map`. **Guide:**
[Clustering mobile genetic elements](guides/mge.md#the-hosts).

```console
metator host [OPTIONS] MAP
```

- **The map** needs the network of `network`, the round `mge` of `mge`, and
  MAGs as its genomes: the final clusters of `validation` or `pipeline`, or
  the bins of another binner set with `metacooler reassign`. `__unbinned__`
  and the MGE genomes are not MAGs, and a MAG counts without the MGE contigs
  it holds. A map without any of these is refused before any work.
- **The rule**: the signal of an MGE cluster with a MAG is the sum of the
  network weights between their contigs. A MAG hosts the MGE cluster when it
  has at least `--min-share` % of the cluster's signal with all the MAGs, and
  at least `--min-contigs` of its contigs touch the cluster. An MGE cluster
  can have no host, one or several. Contacts with contigs in no MAG and with
  other MGE clusters do not count, so a host can come from a small part of
  the contacts: `share_of_all` in `mge_hosts.txt` says how small.
- **Writes** `tables/mge_hosts.txt` and `plots/mge_hosts.pdf` (see
  [Hosts of the MGE clusters](formats.md#hosts-of-the-mge-clusters)),
  replacing earlier ones.
- **Stores** the hosts as columns of the round `mge` (`n_hosts`, `hosts`,
  `host_shares`) and of the genomes table (`n_mges`, `mges`), and the run
  under `metator` → `host` in the metadata, with the digest of the network
  (`network_digest`).

A run that later replaces the round `mge` (`mge`, but not `--random`), the
genomes (`partition`, `validation`, `pipeline` without `--mges`) or the
network (`network`) drops the hosts, with a warning: run `host` again.

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | A map with MGE clusters from `mge`, or an entry of a collection. |
| `-s, --min-share` | `10` | Share (%) of an MGE cluster's signal with the MAGs from which a MAG is a host. |
| `-c, --min-contigs` | `5` | Contigs of a MAG in contact with the MGE cluster from which the MAG is a host. |
| `-o, --outdir` | the current directory | Output directory: the map, `tables/mge_hosts.txt`, `plots/mge_hosts.pdf` and the log. |
| `-F, --force` | off | Replace the map of an earlier run (`mge_hosts.txt` and `mge_hosts.pdf` are replaced without it). |
| `--inplace` | off | Store the hosts in `MAP` itself; refused for a collection entry. |

```bash
# The hosts of the MGE clusters of the map of a run, in binning, updated in place
metator host binning/sample.metacool -o binning
```

## Outputs

`dump`, `export` and `report` read the map of the clusters of a run (from
`pipeline`, or from `validation`, then `mge` and `host` if run), or an entry
of a collection. They only read it: a map that another run writes meanwhile is
refused (run the command again). `export` and `report` refuse any map that
`dump` refuses.

### `dump` { #dump }

Write the tables of a run from the map of its clusters: those of its
validation, which a run keeps only in its temporary directory, and its
network. **API:** `metator.dump.dump`.

```console
metator dump [OPTIONS] MAP
```

| Table | Contents |
|---|---|
| `contig_data_final.txt` | One row per contig, in the order of the map: `Name`, `Size`, `GC_content`, `Hit`, `Shotgun_coverage`, `Restriction_site`, the ID, contigs and size of its core, overlapping and recursive clusters, then `Final_cluster`, its genome (`ND` for none). No `ID` column. See [The contig data tables](formats.md#the-contig-data-tables). |
| `cluster_summary.txt` | One row per genome, the final clusters first, then the MGE genomes, with the columns and formats that the validation recorded. MGE genomes add `GC`, `Hit`, `BinningScore` and `mge` (`True` or `False`); after `host`, `n_mges` and `mges` come last. |
| `binning.txt` | Each contig of a genome and its genome, without a header: the binning that anvi'o imports. |
| `clusters.txt`, `cluster_members.txt` | Every clustering round (`core`, `overlapping`, `recursive_<k>`, `final`, `mge`, `mge_random`), each cluster with its parent, `n_contigs`, `length` and the columns of its round; and the contigs of each cluster. After `mge`, `GC`, `Hit` and `BinningScore` (`NA` for the other rounds); after `host`, `n_hosts`, `hosts` and `host_shares`. |
| `network.txt` | The network, without a header: `contig1`, `contig2` and the weight of their edge (the contact count for `None`), one row per pair of contigs in contact, without contacts within a contig. See [The network](formats.md#the-network). |

- **Genomes**: the tables follow the genomes of the map, which need not be the
  final clusters: the MGE genomes of `mge --reassign`, or a curation by hand
  with `metacooler reassign`. Until a later command changes the genomes or
  the rounds, the tables match those that the run kept in its temporary
  directory.
- **Refused**: a map without the clusters of a validation, or without the
  network or the partition of its run (the error names the step to run); a
  map validated by a development build of 2.0 (run `validation` again); a map
  whose network was replaced after the run (run `network -n <normalization>`
  with the run's normalization, or `partition` and `validation` again).
- **Output directory**: give the dump a directory of its own; a `validation`
  or `pipeline` run would refuse its tables as those of an earlier version.
  The tables are written all at once, so a dump that fails leaves the
  directory as it was.

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The map of the clusters of a run, or an entry of a collection. |
| `-o, --outdir` | the current directory | Output directory of the tables, made if needed. |
| `-F, --force` | off | Replace the tables of an earlier dump. |

```bash
# The tables of a run, from the map of its clusters
metator dump binning/sample.metacool -o tables
```

### `export` { #export }

Write the outputs of a run again from the map of its clusters, without
running it again: in another directory, or in place of the run's own.
**API:** `metator.outputs`.

```console
metator export [OPTIONS] MAP
```

| Output | Contents | Needs |
|---|---|---|
| `fasta/mags_with-mges/` | The FASTA file of each final cluster, MGE contigs included, its records in the order of the map. | `--assembly` |
| `fasta/mags/` | Each MAG as its contigs and stitched, as `pipeline` writes it. | `--assembly` |
| `fasta/mges/`, `fasta/mges_random/` | Each MGE cluster of the round `mge` (`mge_random`), as `mge` writes it. | `--assembly`, the round |
| `cool/` | One balanced cool file per MAG, with the MGE clusters it hosts. | |
| `tables/` | The tables of [`dump`](#dump), and `mge_hosts.txt`, found again with the thresholds that the host run recorded. | `mge_hosts.txt`: hosts in the map |
| `plots/` | The plots of the final clusters, from the miComplete values in the map, and `mge_hosts.pdf`. | `mge_hosts.pdf`: hosts in the map |
| `<name>.report.txt`, `<name>.report.html`, `<name>.report_mqc.yaml` | The report, as [`report`](#report) writes it, last. | |

Without `--assembly`, the FASTA files are left out and the log names the
assembly that the run recorded. The log of the run and its temporary files
cannot be rebuilt (`metacooler fasta MAP --round overlapping` writes the FASTA
files of a round). Earlier outputs are replaced only with `-F`, except
`tables/mge_hosts.txt` and `plots/mge_hosts.pdf`, always replaced.

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The map of the clusters of a run, or an entry of a collection. |
| `-a, --assembly` | none | The assembly FASTA file, plain or gzip-compressed; needed for `fasta/`. |
| `-o, --outdir` | the current directory | Output directory, made if needed. |
| `-F, --force` | off | Replace the folders, tables, cluster plots and report of an earlier run. |
| `-T, --tmpdir` | `./tmp` | Temporary directory, for a decompressed assembly. |

```bash
# Every output of a run, from the map of its clusters and its assembly
metator export binning/sample.metacool -a assembly.fa -o outputs

# The same without the FASTA files: cool/, tables/, plots/ and the report
metator export binning/sample.metacool -o outputs

# The outputs of the run written again in its own output directory
metator export binning/sample.metacool -a assembly.fa -o binning -F
```

### `report` { #report }

Write the report of a run from the map of its clusters. `pipeline` and
`export` write it too. **API:** `metator.report`.

```console
metator report [OPTIONS] MAP
```

| Part | Contents |
|---|---|
| The run | Its steps, input and assembly, the parameters of each step and their command lines. |
| The contigs | Number, length, N50, L50, N90, longest and median contig, GC content, median coverage; and how many are in a final cluster, a MAG, an MGE cluster or no cluster. |
| The contacts | The read contacts that metacooler counted (within a contig, within a genome, between genomes, and their ratios), and the network. |
| The clusters | The clusters of each round, and the final clusters by quality class (HQ, MQ, LQ, Contaminated, Other) with their length, from the miComplete values at full precision. |
| The MAGs | One row per MAG: quality class, length, contigs, completeness, contamination, GC content, N50, abundances, share of its contacts within it, and the MGE clusters it hosts. |
| The MGE clusters | After `mge` and `host`: one row per MGE cluster, with its contigs, length, GC content and hosts, and the share of its signal with each. |

It writes three files in the output directory (see
[The report of a run](formats.md#the-report-of-a-run)):

- `<name>.report.txt`: text in a fixed format (key–value sections, then
  tables) for a script or a MultiQC module.
- `<name>.report_mqc.yaml`: [MultiQC](https://multiqc.info) custom content,
  one row of its General Statistics table named `<name>`: `multiqc` on a
  directory of runs gathers them in one table.
- `<name>.report.html`: a page to read the run: its numbers; a table of the
  MAGs that sorts, filters and downloads as CSV, with links to their FASTA
  and cool files when they are beside it; the MGE clusters and their hosts;
  the command lines; and interactive plots, drawn by plotly.js from its CDN.
  Offline, the page shows its numbers and tables without the plots.

An earlier report is replaced only with `-F`.

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The map of the clusters of a run, or an entry of a collection. |
| `-o, --outdir` | the current directory | Output directory of the report and the log, made if needed. |
| `-F, --force` | off | Replace the report of an earlier run. |

```bash
# The report of a run, from the map of its clusters
metator report binning/sample.metacool -o report

# The reports of the runs of a project, gathered in one table by MultiQC
metator pipeline sample_1.metacool -a sample_1.fa -o runs/sample_1
metator pipeline sample_2.metacool -a sample_2.fa -o runs/sample_2
multiqc runs
```
