# Running the pipeline

`metator pipeline` clusters the contigs of an assembly from the contact map of
a metaHiC library, in one command: it builds the network of contigs from the
map, partitions it, validates the clusters, and writes the contact map of the
clusters to its output directory. This page lists its inputs, the files it writes, and
how to run its steps one by one. The [CLI reference](../cli.md#pipeline) lists
every option, and [Concepts](../concepts.md) explains the steps.

```bash
metator pipeline sample.metacool --enzyme DpnII --depth depth.txt --threads 16 \
    --assembly assembly.fa --outdir binning
```

## Before the run: the assembly

metaTOR clusters the contigs of an assembly; it does not assemble. Build the
assembly from a shotgun library of the same sample, with
[MEGAHIT](https://github.com/voutcn/megahit) or
[metaSPAdes](https://github.com/ablab/spades) for instance, and drop the
contigs shorter than 500 bp, which carry too few contacts to cluster:

```bash
megahit -1 shotgun_R1.fq.gz -2 shotgun_R2.fq.gz --min-contig-len 500 -o assembly
```

The shotgun coverage of the contigs, which `--depth` reads, comes from
`jgi_summarize_bam_contig_depths` of MetaBAT2, run on the shotgun reads
aligned back to the assembly. It is optional: the default normalisation does
not use it.

## Inputs

A run starts from a `.metacool` contact map of the contigs, its first
argument, `MAP`: a file, or an entry of a `.mmetacool` collection
(`file.mmetacool::entries/<name>`). The map is made before any clustering:

| From | The map |
|---|---|
| The reads of a library | `hicstuff pipeline -g assembly.fa -e DpnII -M metacool -P sample R1.fq.gz R2.fq.gz` aligns them and writes `sample.metacool` (from hicstuff 4.0.0: see [Other tools](../other-tools.md#from-the-reads-to-the-map)) |
| Pairs files | `metacooler create pairs --pairs sample.pairs.gz --bins 1000 -o sample.metacool` |
| Several libraries | The map of each, summed with `metacooler merge L1.metacool L2.metacool -o sample.metacool` |

See [Other tools](../other-tools.md#from-the-reads-to-the-map) for these
commands, and [Starting from a contact map](metacool.md#starting-from-a-contact-map)
for what the run reads from the map.

`--assembly` (`-a`) is the assembly that the map was made from, as a FASTA
file, plain or compressed with gzip. A compressed one is decompressed once,
into the temporary directory of the run, for the FASTA files of the
clusters: see `--tmpdir` below.

## Outputs

The pipeline writes into `--outdir` (`-o`, the current directory by default):

| File | What it holds |
|---|---|
| `<MAP>` | The contact map of the clusters, named as `MAP` (`sample.metacool` from `sample.metacool`) |
| `final_clusters/` | One FASTA file per final cluster |
| `metator_pipeline_<date>.log` | The log of the run |

What the run finds is in the map. At its end, the run writes the contact map
of all the contigs, grouped by cluster, to the output directory, under the
name of `MAP` (see [The contact map of the clusters](metacool.md)): the
summary of the clusters in its genomes table, every clustering round (core,
overlapping, recursive, final), and the network of the contigs, as `metator
network` stores it (see
[Building the network](network.md#building-the-network)). `MAP` is left as it
was: metacooler copies it next to the map of the run, which the run edits, and
renames into place once it is validated. `--inplace` writes the map of the
clusters over `MAP` instead, and so does a run whose output directory holds
`MAP` already. A run that stops before leaves `MAP` as it was, which a run
restarts from (`-F` to replace the map and the final clusters it began).

The tables of the run are written in its temporary directory (`--tmpdir`, in
a directory of its own, which the log names), which `--no-clean-up` (`-N`)
keeps:

| File | What it holds |
|---|---|
| `cluster_summary.txt` | One row per final cluster: size, GC content, miComplete's completeness and redundancy, N50, abundance… |
| `contig_data_final.txt` | One row per contig: its measures, and the core, overlapping, recursive and final clusters it belongs to |
| `binning.txt` | The contigs of the final clusters and their cluster |
| `clusters.txt`, `cluster_members.txt` | Every cluster of every round, with its parent and miComplete's columns, and the contigs of each: see [Clustering rounds](metacool.md#clustering-rounds) |
| `overlapping_clusters/`, `recursive_clusters/` | The FASTA files of the clusters that miComplete checked: the overlapping clusters, and the recursive clusters, in `step_<n>/` for step `n` |
| `overlapping_micomplete_results.txt`, `recursive_micomplete_step_<n>.txt` | miComplete's results |
| `contig_data_partition.txt` | The contig data of the partition |

The tables list the contigs of the map, in its order, by their position in it
(`ID`). The map holds what they hold: `metator dump MAP` writes them from the
map of the clusters (see [`dump`](../cli.md#dump)), in the order of that map,
without the `ID`. [File formats](../formats.md) describes each table.

A run refuses to overwrite the final clusters of an earlier one unless given
`--force` (`-F`). The other results of an earlier run in the output directory
(the tables and `recursive_clusters/` of a development build of 2.0, the
outputs of 1.x) are refused, with or without `-F`, and never removed: move
them away, or give the run another `--outdir`.

## Choosing the options

| Option | Default | When to change it |
|---|---|---|
| `--enzyme` | none | Give the enzyme of the library: it counts the restriction sites of the contigs, for the `RS` and `theoritical_hit` normalisations. The `Restriction_site` column of the map gives them without it |
| `--depth` | none | Give it for the `abundance` and `theoritical_hit` normalisations, and the shotgun abundance of the clusters. The `Shotgun_coverage` column of the map gives the coverage without it |
| `--normalization` | `empirical_hit` | See [Normalising the network](network.md#normalising-the-network) |
| `--iterations`, `--overlap` | 100, 80 | See [Partition](binning.md#partition) |
| `--rec-iter`, `--rec-overlap` | 10, 90 | See [Validation](binning.md#validation) |
| `--size` | 500000 | The size, in base pairs, an overlapping cluster must reach to be checked with miComplete, and a recursive cluster must exceed |
| `--prefix` | none | Names the clusters `<prefix>_metator_…` |
| `--threads` | 1 | Threads for the partition and the validation |
| `--tmpdir` | `./tmp` | Where the temporary files go; deleted at the end, unless `--no-clean-up`. It holds the tables of the run, the FASTA files of the clusters that miComplete checks (at most the size of the assembly) and, from a compressed assembly, its decompressed copy: give it room for both. Not in `final_clusters/` of `--outdir`, which the run makes anew, nor in a directory of the results of an earlier run there (`recursive_clusters/`, `final_bin/`...), which the run refuses |

The validation is the step that needs the most memory and time: miComplete
runs Prodigal and HMMER (pyrodigal and pyhmmer) on every cluster.

## Running the steps one by one

The three steps of the pipeline are also commands of their own, to run one
after the other, or to rerun the partition with other parameters without
building the network again:

```bash
# 1. Build the network from the map, and store it in binning/sample.metacool
metator network sample.metacool -a assembly.fa -e DpnII -d depth.txt -o binning

# 2. Partition it: the map stores its core and overlapping clusters, in place
metator partition binning/sample.metacool -t 16 -o binning

# 3. Validate the clusters: final_clusters/, and the map of the clusters
metator validation binning/sample.metacool -a assembly.fa -t 16 -o binning
```

Run step by step, the final clusters go to `final_clusters/`. The map is
what the steps share: `metator network` stores the network in the map of its
output directory, which `partition` and `validation` read and update in
place (it is the map of their output directory too), `partition` stores its
partition there, which `validation` reads, and they write the map that
`metator pipeline MAP` writes from the same clusters: see
[Step by step](metacool.md#step-by-step).

To cluster without validating, run the first two steps alone: the clusters
are then the overlapping clusters of the map, its round `overlapping`, which
`metacooler dump -t members` lists with their contigs, and
`metacooler fasta MAP --fasta assembly.fa --round overlapping -o DIR` writes
as FASTA files.

## The log

The log starts, as hicstuff's does, with the version of metaTOR, the date and
the command line (`cmd`), then the value of each option, and a line of
dashes. It then follows the three steps of the run, each from a line of its
own (`[1/3] network`, `[2/3] partition`, `[3/3] validation`): the contacts of
the map and those between contigs, the iterations of the partition as they
go (a line every tenth of them, at least 10 seconds apart, on a large
network), the core and overlapping clusters found, the clusters that
miComplete checked, each recursive step of the validation (`Recursive step
2/9: ...`, and when no further step is needed), and the number of clusters in
each quality class (see [Cluster quality](../concepts.md#cluster-quality)). It
ends with the map that the run wrote and how long it took:

```console
INFO :: [1/3] network: the contigs of sample.metacool and their contacts
INFO :: 41563 contacts in sample.metacool, 18632 of them between contigs (3D ratio: 0.4483).
...
INFO :: Final contact map written after 0h 0m 5s: binning/sample.metacool
```

The **3D ratio** is the share of the contacts that join two different
contigs: the signal that the clustering uses.
