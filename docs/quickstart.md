# Quickstart

This page makes the contact map of a metaHiC library, clusters the contigs of
the assembly from it, and looks at the clusters and at their contact map. It
uses the test data of metaTOR's repository: the 41,563 pairs of a library on an
assembly of 1,219 contigs, 2.2 Mb in all, so the run takes under a minute.

## Input files

| File | What it is |
|---|---|
| `assembly.fa` | The metagenome assembly, in FASTA format |
| `alignment.pairs` | The pairs of the library: its read pairs aligned on the assembly, in the [pairs format](formats.md#pairs) |
| `depth.txt` | *(optional)* The shotgun coverage of each contig, from `jgi_summarize_bam_contig_depths` |

```bash
base=https://raw.githubusercontent.com/koszullab/metator/master/tests_data
curl -LO "$base/assembly.fa"
curl -LO "$base/alignment.pairs"
curl -LO "$base/depth.txt"
```

## Step 1: make the contact map

metaTOR starts from a contact map of the contigs, a `.metacool` file, made
before any clustering. metacooler, which comes with metaTOR, makes it from the
pairs:

```bash
metacooler create pairs --pairs alignment.pairs --bins 1000 --no-ps \
    --assembly-name assembly -o sample.metacool
```

The map cuts the contigs into 1 kb bins (`--bins 1000`) and puts them all in
one genome, `__unbinned__`. metacooler warns that it clamped 719 pairs: those
of this test file at position 0, which it counts in the first bin of their
contig.

!!! tip "From the reads"

    hicstuff aligns the reads of a library and makes the map in one command:

    ```bash
    hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz
    ```

    `-M metacool` comes with hicstuff's next release: until then, install
    hicstuff from its master branch, or make the map from the pairs files of
    hicstuff 3.3.0 with `metacooler create pairs`, as above.

    It writes the map `hicstuff/sample.metacool`, with one bin per restriction
    fragment, and the pairs `hicstuff/sample.valid_idx.pairs.gz`, which
    `metator mge` reads once indexed with pairix. `metacooler merge` sums the
    maps of several libraries: see
    [Other tools](other-tools.md#from-the-reads-to-the-map).

## Step 2: run the pipeline

```bash
metator pipeline sample.metacool --assembly assembly.fa --enzyme DpnII \
    --depth depth.txt --size 20000 --threads 4 --outdir binning
```

The run reads the network of the contigs from the map. `--enzyme` and
`--depth` add the restriction sites and the shotgun coverage of each contig to
its data; the default normalisation needs neither. At its end, the run writes
the contact map of the clusters over the map: `--output-map` writes it to
another file, and leaves the map as it was.

!!! note "`--size 20000` is for this test library only"

    A cluster is kept when it reaches `--size` base pairs: 500 kb by default,
    the size of a small genome. This library is far too small for clusters that
    large, and at the default size the run stops before validation, with no
    cluster to check. On a real sample, leave `--size` out.

The run logs each step (excerpt):

```console
INFO :: Information of sample.metacool:
INFO :: 41563 contacts in the library.
INFO :: 18632 contacts inter-contigs in the library.
INFO :: 3D ratio : 0.44828332892235884
INFO :: Detect core clusters:
INFO :: 429 core clusters were found.
INFO :: Detect overlapping clusters:
INFO :: 230 overlapping clusters were found.
INFO :: Extract clusters:
INFO :: 24 clusters have been extracted
INFO :: Lauch miComplete quality check.
INFO :: No more contaminated cluster has been found after 1 steps.
INFO :: 24 clusters have been kept after the recursive iterations.
INFO :: HQ MAGs: 0	Total Size: 0
INFO :: MQ MAGs: 0	Total Size: 0
INFO :: LQ MAGs: 0	Total Size: 0
INFO :: Contaminated potential MAGs: 0	Total Size: 0
INFO :: Other clusters: 24	Total Size: 679917
INFO :: Writing sample.metacool from the contact map sample.metacool: 1219 contigs, 301 of them in 24 clusters.
INFO :: sample.metacool is written: the contact map of the clusters.
```

The Louvain algorithm is randomised: another run finds about as many clusters,
not always the same.

## Step 3: look at what it wrote

```console
$ ls binning
binning.txt
cluster_members.txt
cluster_summary.txt
clusters.txt
contig_data_final.txt
final_clusters
metator_pipeline_20261121023231.log
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
| `metator_pipeline_<date>.log` | The log of the run |

`sample.metacool`, next to them, is now the contact map of all the contigs,
grouped by cluster, with the network of the run: what metaTOR measured on
each contig, and the weights of the normalization of the network. [Running the pipeline](guides/pipeline.md#outputs) lists
every output, and [File formats](formats.md) describes the tables.

## Step 4: read the cluster summary

```console
$ cut -f 1-6,9,10 binning/cluster_summary.txt | head -4
	Length	GC-content	Present Markers	Completeness	Redundancy	Contigs	N50
metator_00001_00000	36399.00	62.14	1.00	0.01	1.00	9.00	22786.00
metator_00006_00000	30476.00	61.77	1.00	0.01	1.00	9.00	7448.00
metator_00005_00000	20637.00	63.40	0.00	0.00	0.00	6.00	11969.00
```

Completeness and redundancy are miComplete's, as fractions: a complete
genome without contamination has a completeness and a redundancy of 1. The
clusters of this test library hold a few marker genes at most (a completeness
under 10%), so they are classed as *Other clusters*; on a real sample, the log
counts the high-, medium- and low-quality MAGs.

## Step 5: open the contact map

`sample.metacool` is a [metacooler](https://github.com/js2264/metacooler)
file: a [cooler](https://github.com/open2c/cooler) whose chromosomes are the
contigs, grouped into one genome per final cluster. metacooler reads its
genomes:

```bash
metacooler info sample.metacool
metacooler dump -t genomes -H sample.metacool
```

The genomes table holds the columns of `cluster_summary.txt`, and the contigs
that are in no cluster are in the `__unbinned__` genome.
[The contact map of the clusters](guides/metacool.md) explains the file.

## The same from Python

```python
import metacooler
import pandas as pd

clusters = pd.read_csv("binning/cluster_summary.txt", sep="\t", index_col=0)
contigs = pd.read_csv("binning/contig_data_final.txt", sep="\t", index_col=0).set_index("Name", drop=False)

mc = metacooler.MetaCooler("sample.metacool")
mc.genomes   # one row per cluster, and __unbinned__, with the columns of cluster_summary.txt
mc.contigs   # one row per contig: GC content, hits, coverage, restriction sites
```

## Next steps

- [Concepts](concepts.md): the network, the partition, the validation.
- [Running the pipeline](guides/pipeline.md): every input, output and option,
  and the steps run one by one.
- [Other tools](other-tools.md): the map from the reads, the contact map of one
  cluster, its scaffold, and the quality of a library, with hicstuff,
  metacooler and instaGRAAL.
- [CLI reference](cli.md): every command.
