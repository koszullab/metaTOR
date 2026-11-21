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
| Pairs files | `mge` takes the pairs files of a run as arguments: those that hicstuff writes with the map, `<prefix>.valid_idx.pairs.gz`, once indexed with `pairix -p pairs` (see [Clustering mobile genetic elements](guides/mge.md#running-it)); it stores its MGE clusters in the map given with `--map`. |
| `-F, --force` | Replace the outputs of an earlier run, in the output directory and at `--output-map`. Without it, a command refuses to overwrite them. A run from a map updates that map in place without it, and warns, before any work, of the clusters and clustering rounds that it replaces there: `--output-map` writes the map of the run to another file instead, and leaves the map as it was. |
| Writing a map | Through metacooler: a run holds the lock of the map it writes, `<map>.lock` (FILE's with `--output-map`), from its checks to its end, and writes the map in an edit session, one staged copy next to it (`.<map>.<uuid>.<n>.d/`), validated, then renamed over it, so a run that stops leaves the map as it was. `MAP` given as a symbolic link is updated where the link points; a symbolic link at `--output-map` is replaced (`-F`), not its target, and a hard link of `MAP` there is refused. A map that another process holds, a map without write permission, or one whose directory is not writable (with room for a copy of the map and for the new map), is refused before any work: write a copy with `--output-map`. A run killed by SIGTERM (a job scheduler's time limit) or SIGKILL leaves its lock and its staged copy, which the next run that writes the map on the same host takes over and removes; on another host, delete the files that the error names. |
| `-T, --tmpdir` | Each run works in a directory of its own inside it, removed at the end, unless `--no-clean-up` or the run fails. |
| `-N, --no-clean-up` | Keep the temporary files. |
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
| `-a, --assembly` | *required* | The assembly that the map was made from: a FASTA file, plain or compressed with gzip. |
| `-d, --depth` | none | The coverage of the contigs, from `jgi_summarize_bam_contig_depths` of MetaBAT2. The `Shotgun_coverage` column of the map gives it without `--depth`, which replaces it: see [Starting from a contact map](guides/metacool.md#starting-from-a-contact-map). |
| `-e, --enzyme` | none | The restriction enzymes of the library. The `Restriction_site` column of the map gives the restriction sites without `--enzyme`, which replaces it. |
| `-m, --normalization` | `empirical_hit` | `None`, `abundance`, `length`, `RS`, `empirical_hit` or `theoritical_hit`. `RS` and `theoritical_hit` need `--enzyme`; `abundance` and `theoritical_hit` need `--depth`, or the columns of the map that give them. See [Normalising the network](guides/network.md#normalising-the-network). |
| `-A, --algorithm` | `louvain` | The partition algorithm: `louvain` or `leiden`. |
| `-i, --iterations` | `100` | The iterations of the partition. |
| `-O, --overlap` | `80` | The percentage of iterations in which two core clusters must fall together to be merged into one overlapping cluster. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, a cluster must reach to be kept. |
| `-j, --rec-iter` | `10` | The iterations of the recursive partition of the contaminated clusters. |
| `-P, --rec-overlap` | `90` | The overlap threshold of the recursive partition, in percent. |
| `-o, --outdir` | the current directory | The output directory. |
| `-p, --prefix` | none | Name the clusters `<prefix>_metator_…`. Letters, digits, `_` and `-` only. |
| `-F, --force` | off | Replace the directories of clusters of an earlier run, and an existing `--output-map` file. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-t, --threads` | `1` | Threads for the partition and the validation (miComplete). |
| `-N, --no-clean-up` | off | Keep the temporary files. |
| `--no-metacool` | off | Do not write the contact map of the clusters: the map of the run is left as it was. |
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
`metator` key of the metadata of the map. The map is updated in place, or
written to `--output-map`. It writes no table: `network.txt` and
`contig_data_network.txt` are gone. To export the network, run
`metacooler network MAP --level contig --balance metator_<normalization>`
(the raw counts without `--balance`).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | As in `pipeline`: the network is stored in it, in place, or in `--output-map`, which an entry of a collection needs. |
| `-a, --assembly` | *required* | As in `pipeline`. |
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
and overlapping clusters: the second step. Writes `overlapping_clusters/` and
`contig_data_partition.txt`. **Guide:**
[Partition and validation](guides/binning.md).

```console
metator partition [OPTIONS] MAP
```

The run reads the network and the data of the contigs from the map, which
`metator network MAP` stored there, and refuses a map without them. It stores
its overlapping clusters that reach `--size`, those of
`overlapping_clusters/`, with their contigs, in the clustering rounds of the
map, in place of those it held: the map is updated in place, or written to
`--output-map`. Its genomes and contigs tables are left as they were. The IDs
of the contigs in `contig_data_partition.txt` are their positions in the map.
See [Step by step](guides/metacool.md#step-by-step).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The `.metacool` contact map of the contigs of the run, which `metator network` stored the network in, and which the run stores its overlapping clusters in; or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`), with `--output-map`, the map it writes them to. |
| `-a, --assembly` | *required* | The assembly, as a FASTA file. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `100` | The iterations of the partition. |
| `-O, --overlap` | `80` | The overlap threshold, in percent. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, an overlapping cluster must reach to be written. |
| `-o, --outdir` | `./` | The output directory. |
| `-p, --prefix` | none | Name the clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace `overlapping_clusters/`, and an existing `--output-map` file. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-t, --threads` | `1` | Threads for the partition. |
| `-N, --no-clean-up` | off | Keep the temporary files. |
| `--output-map` | none | Write the map with the overlapping clusters to this file (replaced only with `-F`), and leave `MAP` as it was. By default, the run updates `MAP` in place. |

```bash
metator partition sample.metacool -a assembly.fa -t 16 -o binning
```

### `validation` { #validation }

Check the clusters with miComplete, and partition the contaminated ones
again: the third step. Writes `final_clusters/`, `recursive_clusters/`,
`cluster_summary.txt`, `contig_data_final.txt` and `binning.txt`. **Guide:**
[Partition and validation](guides/binning.md#validation).

```console
metator validation [OPTIONS] MAP
```

A cluster is partitioned again when its weighted completeness reaches 33% and
its contamination, (weighted redundancy − 1) / weighted completeness, reaches
5%. Its parts replace it when the most complete of them keeps enough of its
completeness.

The run partitions the network that `metator network MAP` stored in the map
(refused without it), its contigs named in `contig_data_partition.txt` by
their IDs, then writes the map of its clusters, with all its clustering
rounds and the normalization of its network, as `metator pipeline MAP` does
at its end: in place, or to `--output-map`. See
[Step by step](guides/metacool.md#step-by-step).

| Option | Default | Description |
|---|---|---|
| `MAP` | *required* | The `.metacool` contact map of the contigs of the run, which `metator network` stored the network in, and which becomes the map of its clusters; or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`), with `--output-map`, the map it writes them to. |
| `-a, --assembly` | *required* | The assembly, as a FASTA file. |
| `-c, --contigs` | *required* | The contig data table of the partition, `contig_data_partition.txt`. |
| `-f, --fasta` | *required* | The directory of the clusters to check, `overlapping_clusters/`. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `10` | The iterations of the recursive partition. |
| `-O, --overlap` | `90` | The overlap threshold of the recursive partition, in percent. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, a cluster must reach to be kept. |
| `-o, --outdir` | `./` | The output directory. |
| `-p, --prefix` | none | Name the clusters `<prefix>_metator_…`. |
| `-F, --force` | off | Replace `recursive_clusters/` and `final_clusters/`, and an existing `--output-map` file. Refused when `--fasta` is in one of these directories, or in the target of one that is a symbolic link. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-t, --threads` | `1` | Threads for the partition and miComplete. |
| `-N, --no-clean-up` | off | Keep the temporary files. |
| `--output-map` | none | Write the map of the clusters to this file (replaced only with `-F`), and leave `MAP` as it was. By default, the run updates `MAP` in place. |

```bash
# The map of a run step by step: it becomes the map of the clusters
metator validation sample.metacool -a assembly.fa -c binning/contig_data_partition.txt \
    -f binning/overlapping_clusters -t 16 -o binning
```

## After the clustering

### `mge` { #mge }

Cluster the contigs of mobile genetic elements from the contacts between them.
**API:** `metator.mge.mge_binning`, `metator.metacool.write_round`.
**Guide:** [Clustering mobile genetic elements](guides/mge.md).

```console
metator mge [OPTIONS] PAIRSFILE...
```

`PAIRSFILE` is the pairs files of the run: the contacts within an MGE contig
are counted between reads more than 1 kb apart, which a contact map does not
tell, so a map given in their place is refused. With `--map MAP`, the run
also stores its MGE clusters in the contact map of the run: as the round
`mge` of its clustering rounds, with the columns of `mges_cluster_summary.tsv`
and the contigs of each cluster, in place of the round `mge` of an earlier
run, after the rounds of the clustering run, which are kept; its genomes and
contigs tables, the final clusters, are left as they were. The map is updated
in place, or written to `--output-map`, as `metator partition MAP` writes it
(see [Step by step](guides/metacool.md#step-by-step)), and refused as it
refuses it, before any work. `--random` is refused with `--map`. A run that
replaces the clusters of the map drops the round `mge`, with a warning.

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | The assembly, as a FASTA file. |
| `-c, --contigs` | *required* | The contig data table of the run. |
| `-b, --binning` | *required* | The clusters of the contigs: `binning.txt`, or a collection exported by anvi'o. |
| `-m, --mges` | *required* | The MGE contigs, one name per line. |
| `--map` | none | A `.metacool` contact map of the run, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`), that the run stores its MGE clusters in. |
| `-s, --threshold-bin` | `0.8` | The score above which two contigs, or clusters, are merged. |
| `-r, --random` | off | Shuffle the contigs between the clusters: a random clustering of the same shape. Refused with `--map`. |
| `-o, --outdir` | the current directory | The output directory. |
| `-F, --force` | off | Replace an existing `--output-map` file. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-N, --no-clean-up` | off | Keep the temporary files, and the index of the assembly. |
| `--output-map` | none | With `--map`: write the map with the MGE clusters to this file (replaced only with `-F`), and leave `MAP` as it was. |

```bash
# Index the pairs that hicstuff wrote with the map (once): mge needs the index
# of a compressed pairs file
python -c 'import pypairix; pypairix.build_index("hicstuff/sample.valid_idx.pairs.gz", "pairs")'

metator mge -a assembly.fa -c binning/contig_data_final.txt -b binning/binning.txt \
    -m mges.txt -o mge hicstuff/sample.valid_idx.pairs.gz

# Also store the MGE clusters in the map of the run
metator mge -a assembly.fa -c binning/contig_data_final.txt -b binning/binning.txt \
    -m mges.txt -o mge --map hicstuff/sample.metacool hicstuff/sample.valid_idx.pairs.gz
```
