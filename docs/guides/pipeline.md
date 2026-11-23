# Running the pipeline

`metator pipeline` clusters the contigs of an assembly from the contact map of
a metaHiC library, in one command: it builds the network of contigs from the
map, partitions it, validates the clusters, and writes the contact map of the
clusters over the map. This page lists its inputs, the files it writes, and
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
| The reads of a library | `hicstuff pipeline -g assembly.fa -e DpnII -M metacool -P sample R1.fq.gz R2.fq.gz` aligns them and writes `sample.metacool` (from hicstuff's next release: see [Other tools](../other-tools.md#from-the-reads-to-the-map)) |
| Pairs files | `metacooler create pairs --pairs sample.pairs.gz --bins 1000 -o sample.metacool` |
| Several libraries | The map of each, summed with `metacooler merge L1.metacool L2.metacool -o sample.metacool` |

See [Other tools](../other-tools.md#from-the-reads-to-the-map) for these
commands, and [Starting from a contact map](metacool.md#starting-from-a-contact-map)
for what the run reads from the map.

`--assembly` (`-a`) is the assembly that the map was made from, as a FASTA
file, plain or compressed with gzip.

## Outputs

The pipeline writes into `--outdir` (`-o`, the current directory by default):

| File | What it holds |
|---|---|
| `final_clusters/` | One FASTA file per final cluster |
| `cluster_summary.txt` | One row per final cluster: size, GC content, miComplete's completeness and redundancy, N50, abundance… |
| `contig_data_final.txt` | One row per contig: its measures, and the core, overlapping, recursive and final clusters it belongs to |
| `binning.txt` | The contigs of the final clusters and their cluster |
| `network.txt` | The normalised network of contigs |
| `contig_data_network.txt` | One row per contig: its measures, as the network step wrote them |
| `overlapping_clusters/` | The FASTA files of the overlapping clusters, before validation |
| `recursive_clusters/` | The FASTA files of the recursive clusters, in `step_<n>/` for step `n` |
| `overlapping_micomplete_results.txt` | miComplete's results on the overlapping clusters |
| `recursive_micomplete_step_<n>.txt` | miComplete's results on the recursive clusters of step `n` |
| `clusters.txt`, `cluster_members.txt` | Every cluster of every round (overlapping, recursive, final), with its parent and miComplete's columns, and the contigs of each: see [Clustering rounds](metacool.md#clustering-rounds) |
| `metator_pipeline_<date>.log` | The log of the run |

At its end, the run writes the contact map of all the contigs, grouped by
cluster, over `MAP` (see [The contact map of the clusters](metacool.md)):
`--output-map FILE` writes it to FILE and leaves `MAP` as it was, and
`--no-metacool` writes none. A run that stops before leaves `MAP` as it was,
which a run restarts from (`-F` to replace the directories of clusters it
began). [File formats](../formats.md) describes each table.

A run refuses to overwrite the clusters of an earlier one unless given
`--force` (`-F`).

## Choosing the options

| Option | Default | When to change it |
|---|---|---|
| `--enzyme` | none | Give the enzyme of the library: it counts the restriction sites of the contigs, for the `RS` and `theoritical_hit` normalisations. The `Restriction_site` column of the map gives them without it |
| `--depth` | none | Give it for the `abundance` and `theoritical_hit` normalisations, and the shotgun abundance of the clusters. The `Shotgun_coverage` column of the map gives the coverage without it |
| `--normalization` | `empirical_hit` | See [Normalising the network](network.md#normalising-the-network) |
| `--iterations`, `--overlap` | 100, 80 | See [Partition](binning.md#partition) |
| `--rec-iter`, `--rec-overlap` | 10, 90 | See [Validation](binning.md#validation) |
| `--size` | 500000 | The size, in base pairs, a cluster must reach to be kept |
| `--prefix` | none | Names the clusters `<prefix>_metator_…` |
| `--threads` | 1 | Threads for the partition and the validation |
| `--tmpdir` | `./tmp` | Where the temporary files go; deleted at the end, unless `--no-clean-up` |

The validation is the step that needs the most memory and time: miComplete
runs Prodigal and HMMER (pyrodigal and pyhmmer) on every cluster.

## Running the steps one by one

The three steps of the pipeline are also commands of their own, to run one
after the other, or to rerun the partition with other parameters without
building the network again:

```bash
# 1. Build the network from the map: network.txt, contig_data_network.txt
metator network sample.metacool -a assembly.fa -e DpnII -d depth.txt -t 16 -o binning

# 2. Partition it: overlapping_clusters/, contig_data_partition.txt
metator partition -a assembly.fa -c binning/contig_data_network.txt \
    -n binning/network.txt -t 16 -o binning

# 3. Validate the clusters: final_clusters/, cluster_summary.txt, contig_data_final.txt
metator validation -a assembly.fa -c binning/contig_data_partition.txt \
    -n binning/network.txt -f binning/overlapping_clusters -t 16 -o binning
```

Run step by step, the final clusters go to `final_clusters/`, and the
intermediate tables stay. `metator network` only reads the map; given the map
as their first argument, `partition` and `validation` write the map that
`metator pipeline MAP` writes from the same clusters: see
[Step by step](metacool.md#step-by-step).

To cluster without validating, run the first two steps alone: the clusters
are then the overlapping clusters of `overlapping_clusters/`, and
`contig_data_partition.txt` gives the cluster of each contig.

## The log

The log starts with the version of metaTOR, the date and the command, then
the parameters, and records each step: the contacts of the map and those
between contigs, the core and overlapping clusters found, the clusters
validated, and at the end the number of clusters in each quality class (see
[Cluster quality](../concepts.md#cluster-quality)).

```console
INFO :: 41563 contacts in the library.
INFO :: 18632 contacts inter-contigs in the library.
INFO :: 3D ratio : 0.44828332892235884
```

The **3D ratio** is the share of the contacts that join two different
contigs: the signal that the clustering uses.
