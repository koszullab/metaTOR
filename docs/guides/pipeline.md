# Running the pipeline

`metator pipeline` clusters the contigs of an assembly from a metaHiC library,
in one command: it aligns the reads, builds the network of contigs, partitions
it, and validates the clusters. This page lists its inputs, the files it writes,
and how to run its steps one by one. The [CLI reference](../cli.md#pipeline)
lists every option, and [Concepts](../concepts.md) explains the steps.

```bash
metator pipeline --aligner bwa --enzyme DpnII --depth depth.txt --threads 16 \
    --forward R1.fq.gz --reverse R2.fq.gz --assembly assembly.fa --outdir binning
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

A run starts from the reads, their alignments or their pairs, given to
`--forward` (`-1`) and `--reverse` (`-2`), as `--start` (`-S`) says:

| `--start` | `--forward`, `--reverse` | The run starts with |
|---|---|---|
| `fastq` (default) | The R1 and R2 reads, in FASTQ files | Aligning the reads |
| `bam` | The alignments: with bowtie2, one BAM file per mate; with bwa (`--aligner bwa`), one BAM file of both mates, given to `--forward` | Pairing the mates |
| `pair` | One pairs file per library, given to `--forward` | Building the network |

Or it starts from a `.metacool` contact map of the contigs, its first
argument, in place of `--forward`, `--reverse` and `--start`
(`metator pipeline sample.metacool -a assembly.fa`): it starts with building
the network.

Several libraries are given as comma-separated lists, in the same order in
`--forward` and `--reverse`: `-1 L1_R1.fq.gz,L2_R1.fq.gz -2 L1_R2.fq.gz,L2_R2.fq.gz`.
Each library is aligned on its own, and their contacts are summed in the
network.

The BAM files of `--start bam` must be sorted by read name, both files of a
library in the same order. A pairs file of `--start pair` has the 7 columns of
the [pairs format](../formats.md#pairs): read name, then the contig, position
and strand of each mate. A run from a contact map is described in
[The contact map of the clusters](metacool.md#starting-from-a-contact-map).

`--assembly` (`-a`) is the assembly, as a FASTA file, or the basename of its
bowtie2 or bwa index (with bwa, `<basename>.fa` or `<basename>.fasta` must be
next to it). From a FASTA file, an index found next to it (as
`bwa index assembly.fa` or `bowtie2-build assembly.fa assembly.fa` writes it)
is used, and otherwise one is built in the temporary directory, at every run.

## Outputs

The pipeline writes into `--outdir` (`-o`, the current directory by default):

| File | What it holds |
|---|---|
| `final_clusters/` | One FASTA file per final cluster |
| `cluster_summary.txt` | One row per final cluster: size, GC content, miComplete's completeness and redundancy, N50, abundance… |
| `contig_data_final.txt` | One row per contig: its measures, and the core, overlapping, recursive and final clusters it belongs to |
| `binning.txt` | The contigs of the final clusters and their cluster |
| `metator.metacool` | The contact map of all the contigs, grouped by cluster (`<prefix>.metacool` with `--prefix`): from reads, alignments or pairs; a run from a map updates that map instead, or writes `--output-map` |
| `metator.mmetacool` | With several libraries, the contact map of each library, in no cluster (`<prefix>.mmetacool` with `--prefix`) |
| `network.txt` | The normalised network of contigs, all libraries summed |
| `network_<i>.txt` | The network of library `i`, from 0 |
| `alignment_<i>.bam` | With bwa: the alignments of library `i`, both mates, sorted by read name |
| `alignment_<i>_for.bam`, `alignment_<i>_rev.bam` | With bowtie2, or in iterative mode: the alignments of each mate |
| `overlapping_clusters/` | The FASTA files of the overlapping clusters, before validation |
| `recursive_clusters/` | The FASTA files of the recursive clusters, in `step_<n>/` for step `n` |
| `overlapping_micomplete_results.txt` | miComplete's results on the overlapping clusters |
| `recursive_micomplete_step_<n>.txt` | miComplete's results on the recursive clusters of step `n` |
| `clusters.txt`, `cluster_members.txt` | Every cluster of every round (overlapping, recursive, final), with its parent and miComplete's columns, and the contigs of each: see [Clustering rounds](metacool.md#clustering-rounds) |
| `clustering_matrix*.npz` | With `--cluster-matrix`: how often each two contigs were clustered together |
| `metator_pipeline_<date>.log` | The log of the run |

The pairs files of each library are written to the temporary directory as
they are made, read once for the contact map of the contigs and the
network, then removed. With `--no-clean-up` (`-N`), the run sorts,
compresses and indexes them in the output directory instead,
`alignment_<i>_sorted.pairs.gz` and their index, and keeps them, for
`metator mge` later. Without `-N`, a run that stops
before its network is built leaves them in its temporary directory; one
that stops later leaves none, but leaves the contact map of its contigs,
which a run restarts from (`metator pipeline binning/metator.metacool … -F`,
`-F` to replace the directories of clusters it began).
[File formats](../formats.md) describes each table.

A run refuses to overwrite the clusters of an earlier one unless given
`--force` (`-F`).

## Choosing the options

| Option | Default | When to change it |
|---|---|---|
| `--aligner`, `--aligner-mode` | `bowtie2`, `normal` | `--aligner bwa` is recommended; see [Aligning the reads](network.md#aligning-the-reads) |
| `--enzyme` | none | Give the enzyme of the library: it counts the restriction sites of the contigs, for the `RS` and `theoritical_hit` normalisations and `--aligner-mode cutsite`. From a map, its `Restriction_site` column gives them without it |
| `--depth` | none | Give it for the `abundance` and `theoritical_hit` normalisations, and the shotgun abundance of the clusters. From a map, its `Shotgun_coverage` column gives the coverage without it |
| `--normalization` | `empirical_hit` | See [Normalising the network](network.md#normalising-the-network) |
| `--iterations`, `--overlap` | 100, 80 | See [Partition](binning.md#partition) |
| `--rec-iter`, `--rec-overlap` | 10, 90 | See [Validation](binning.md#validation) |
| `--size` | 500000 | The size, in base pairs, a cluster must reach to be kept |
| `--prefix` | none | Names the clusters `<prefix>_metator_…`, the map `<prefix>.metacool` of a run from reads, alignments or pairs and, with several libraries, the maps of the libraries `<prefix>.mmetacool` |
| `--threads` | 1 | Threads for the alignment and the partition |
| `--tmpdir` | `./tmp` | Where the index, the pairs files (without `--no-clean-up`) and the temporary files go; deleted at the end, unless `--no-clean-up` |

The validation is the step that needs the most memory and time: miComplete
runs Prodigal and HMMER (pyrodigal and pyhmmer) on every cluster.

## Running the steps one by one

The three steps of the pipeline are also commands of their own, to run one
after the other, or to rerun the partition with other parameters without
aligning the reads again:

```bash
# 1. Align the reads and build the network: network.txt, contig_data_network.txt
metator network -1 R1.fq.gz -2 R2.fq.gz -a assembly.fa -b bwa -e DpnII -d depth.txt \
    -t 16 -o binning

# 2. Partition it: overlapping_clusters/, contig_data_partition.txt
metator partition -a assembly.fa -c binning/contig_data_network.txt \
    -n binning/network.txt -t 16 -o binning

# 3. Validate the clusters: final_clusters/, cluster_summary.txt, contig_data_final.txt
metator validation -a assembly.fa -c binning/contig_data_partition.txt \
    -n binning/network.txt -f binning/overlapping_clusters -t 16 -o binning
```

Run step by step, the final clusters go to `final_clusters/`, and the
intermediate tables stay. `metator network` keeps the pairs files. From a
contact map of the contigs, each step takes the map as its first argument,
and the steps write the map that `metator pipeline MAP` writes from the same
clusters: see [Step by step](metacool.md#step-by-step).

To cluster without validating, run the first two steps alone: the clusters
are then the overlapping clusters of `overlapping_clusters/`, and
`contig_data_partition.txt` gives the cluster of each contig.

## The log

The log starts with the version of metaTOR, the date and the command, then
the parameters, and records each step: the pairs aligned, the contacts between
contigs, the core and overlapping clusters found, the clusters validated, and
at the end the number of clusters in each quality class (see
[Cluster quality](../concepts.md#cluster-quality)).

```console
INFO :: 61974 contacts in the library.
INFO :: 28510 contacts inter-contigs in the library.
INFO :: 3D ratio : 0.46003162616581145
```

The **3D ratio** is the share of the contacts that join two different
contigs: the signal that the clustering uses.
