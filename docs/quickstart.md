# Quickstart

This page clusters the contigs of an assembly from the reads of a metaHiC
library, and looks at the clusters and at their contact map. It uses the test
data of metaTOR's repository: 70,805 read pairs on an assembly of 1,219
contigs, 2.2 Mb in all, so the run takes about a minute.

## Input files

| File | What it is |
|---|---|
| `assembly.fa` | The metagenome assembly, in FASTA format |
| `for_paired.fq.gz` | The first mate of each metaHiC read pair (R1) |
| `rev_paired.fq.gz` | The second mate, in the same order (R2) |
| `depth.txt` | *(optional)* The shotgun coverage of each contig, from `jgi_summarize_bam_contig_depths` |

```bash
base=https://raw.githubusercontent.com/koszullab/metator/master/tests_data
curl -LO "$base/assembly.fa"
curl -LO "$base/for_paired.fq.gz"
curl -LO "$base/rev_paired.fq.gz"
curl -LO "$base/depth.txt"
```

## Step 1: run the pipeline

```bash
metator pipeline --forward for_paired.fq.gz --reverse rev_paired.fq.gz \
    --assembly assembly.fa --aligner bwa --enzyme DpnII --depth depth.txt \
    --size 20000 --threads 4 --outdir binning
```

`--aligner bwa` aligns the read pairs with `bwa mem -5SP`, the recommended
aligner: see [Aligning reads](guides/network.md#aligning-the-reads).
`--enzyme` and `--depth` add the restriction sites and the shotgun coverage of
each contig to its data; the default normalisation needs neither.

!!! note "`--size 20000` is for this test library only"

    A cluster is kept when it reaches `--size` base pairs: 500 kb by default,
    the size of a small genome. This library is far too small for clusters that
    large, and at the default size the run stops before validation, with no
    cluster to check. On a real sample, leave `--size` out.

The run logs each step (excerpt):

```console
INFO :: Alignment of for_paired.fq.gz and rev_paired.fq.gz:
INFO :: 61974 pairs aligned.
INFO :: Information of binning/metator.metacool:
INFO :: 61974 contacts in the library.
INFO :: 28510 contacts inter-contigs in the library.
INFO :: 3D ratio : 0.46003162616581145
INFO :: Detect core clusters:
INFO :: 550 core clusters were found.
INFO :: Detect overlapping clusters:
INFO :: 290 overlapping clusters were found.
INFO :: Extract clusters:
INFO :: 13 clusters have been extracted
INFO :: Lauch miComplete quality check.
INFO :: No more contaminated cluster has been found after 1 steps.
INFO :: 13 clusters have been kept after the recursive iterations.
INFO :: HQ MAGs: 0	Total Size: 0
INFO :: MQ MAGs: 0	Total Size: 0
INFO :: LQ MAGs: 0	Total Size: 0
INFO :: Contaminated potential MAGs: 0	Total Size: 0
INFO :: Other clusters: 13	Total Size: 380921
INFO :: Writing binning/metator.metacool from the contact map binning/metator.metacool: 1219 contigs, 197 of them in 13 clusters.
```

The Louvain algorithm is randomised: another run finds about as many clusters,
not always the same.

## Step 2: look at what it wrote

```console
$ ls binning
alignment_0.bam
binning.txt
cluster_members.txt
cluster_summary.txt
clusters.txt
contig_data_final.txt
final_clusters
metator.metacool
metator_pipeline_20261118125820.log
network.txt
network_0.txt
overlapping_clusters
overlapping_micomplete_results.txt
recursive_clusters
```

| File | What it holds |
|---|---|
| `final_clusters/` | One FASTA file per final cluster: the result |
| `cluster_summary.txt` | One row per final cluster: size, GC content, completeness, redundancy, N50, Hi-C abundance… |
| `clusters.txt`, `cluster_members.txt` | Every cluster of every round (overlapping, recursive, final), with its parent and miComplete's columns, and the contigs of each |
| `contig_data_final.txt` | One row per contig: its size, GC content, hits, and the clusters it went through |
| `binning.txt` | The contigs of the final clusters and their cluster, as anvi'o imports a binning |
| `metator.metacool` | The contact map of all the contigs, grouped by cluster |
| `network.txt` | The normalised network of contigs |
| `metator_pipeline_<date>.log` | The log of the run |

[Running the pipeline](guides/pipeline.md#outputs) lists every output, and
[File formats](formats.md) describes the tables.

## Step 3: read the cluster summary

```console
$ cut -f 1-6,9,10 binning/cluster_summary.txt | head -4
	Length	GC-content	Present Markers	Completeness	Redundancy	Contigs	N50
metator_00001_00000	49021.00	62.03	1.00	0.01	1.00	15.00	9019.00
metator_00007_00000	24565.00	62.08	1.00	0.01	1.00	12.00	1777.00
metator_00009_00000	24472.00	60.11	0.00	0.00	0.00	9.00	4643.00
```

Completeness and redundancy are miComplete's, as fractions: a complete
genome without contamination has a completeness and a redundancy of 1. The
clusters of this test library hold a few marker genes at most (a completeness
under 10%), so they are classed as *Other clusters*; on a real sample, the log
counts the high-, medium- and low-quality MAGs.

## Step 4: open the contact map

`metator.metacool` is a [metacooler](https://github.com/js2264/metacooler)
file: a [cooler](https://github.com/open2c/cooler) with one bin per contig,
whose contigs are grouped into one genome per final cluster. metacooler, which
comes with metaTOR, reads its genomes:

```bash
metacooler info binning/metator.metacool
metacooler dump -t genomes -H binning/metator.metacool
```

The genomes table holds the columns of `cluster_summary.txt`, and the contigs
that are in no cluster are in the `__unbinned__` genome.
[The contact map of the clusters](guides/metacool.md) explains the file.

## The same from Python

```python
import metacooler
import metator.io as mio
import pandas as pd

clusters = pd.read_csv("binning/cluster_summary.txt", sep="\t", index_col=0)
contigs = mio.read_contig_data("binning/contig_data_final.txt")

mc = metacooler.MetaCooler("binning/metator.metacool")
mc.genomes   # one row per cluster, and __unbinned__, with the columns of cluster_summary.txt
mc.contigs   # one row per contig: GC content, hits, coverage, restriction sites
```

## Next steps

- [Concepts](concepts.md): the network, the partition, the validation.
- [Running the pipeline](guides/pipeline.md): every input, output and option,
  and the steps run one by one.
- [Other tools](other-tools.md): the contact map of one cluster, its
  scaffold, and the quality of a library, with metacooler, hicstuff and
  instaGRAAL.
- [CLI reference](cli.md): every command.
