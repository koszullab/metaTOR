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

    `-M metacool` needs hicstuff 4.0.0 or later: with an earlier hicstuff,
    make the map from its pairs files with `metacooler create pairs`, as
    above.

    It writes the map `hicstuff/sample.metacool`, with one bin per restriction
    fragment, and the pairs `hicstuff/sample.valid_idx.pairs.gz`.
    `metacooler merge` sums the maps of several libraries: see
    [Other tools](other-tools.md#from-the-reads-to-the-map).

!!! tip "Or let the pipeline make the map"

    `metator pipeline` and `metator network` make the map themselves, one bin
    per contig, from the pairs files or the BAM files of the run, several
    libraries summed:

    ```bash
    metator pipeline alignment.pairs --assembly assembly.fa --enzyme DpnII \
        --depth depth.txt --size 20000 --threads 4 --outdir binning --no-clean-up
    ```

    makes `binning/assembly.metacool` from `alignment.pairs` (named after the
    assembly; `--prefix sample` names it `sample.metacool`), then runs the
    three steps below on it, and makes it the map of the clusters. See
    [The inputs of a run](cli.md#inputs).

## Step 2: run the pipeline

```bash
metator pipeline sample.metacool --assembly assembly.fa --enzyme DpnII \
    --depth depth.txt --size 20000 --threads 4 --outdir binning --no-clean-up
```

The run reads the network of the contigs from the map. `--enzyme` and
`--depth` add the restriction sites and the shotgun coverage of each contig to
its data; the default normalisation needs neither. At its end, the run writes
the contact map of the clusters to its output directory,
`binning/sample.metacool`, and leaves `sample.metacool` as it was
(`--inplace` would update it in place). `--no-clean-up` keeps its temporary
directory, in `tmp/`, with the tables of the run.

!!! note "`--size 20000` is for this test library only"

    A cluster is kept when it reaches `--size` base pairs: 500 kb by default,
    the size of a small genome. This library is far too small for clusters that
    large, and at the default size the run stops before validation, with no
    cluster to check. On a real sample, leave `--size` out.

The log starts with the version, the command line and the value of each
option, then follows the four steps of the run, says that it has no steps
`mge` and `host`, without the list of the MGE contigs that `--mges` gives it
(see [Running the pipeline](guides/pipeline.md#mobile-genetic-elements-and-their-hosts)),
then where its last step, `export`, wrote its outputs, and ends with how long
it took (excerpt):

```console
INFO :: ## MetaTOR: v1.4.0 log file
INFO :: ## date: 2026-11-23 15:23:35
INFO :: ## cmd: metator pipeline sample.metacool --assembly assembly.fa --enzyme DpnII --depth depth.txt --size 20000 --threads 4 --outdir binning --no-clean-up
INFO :: ## map: sample.metacool
INFO :: ## assembly: assembly.fa
...
INFO :: ## inplace: False
INFO :: ---
INFO :: The run writes its map to binning/sample.metacool, and leaves sample.metacool as it is (--inplace updates it).
INFO :: [1/4] network: the contigs of sample.metacool and their contacts
INFO :: 41563 contacts in sample.metacool, 18632 of them between contigs (3D ratio: 0.4483).
INFO :: [2/4] partition: the clusters of the network
INFO :: 100 iterations of louvain on the network of 1058 contigs and 13415 edges.
INFO :: louvain: 100/100 iterations, done in 0h 0m 1s.
INFO :: 420 core clusters, merged into 229 overlapping clusters.
INFO :: [3/4] validation: miComplete checks the clusters, and the contaminated ones are partitioned again
INFO :: 26 of the 229 overlapping clusters reach 20000 bp, 0.764 Mb in all: miComplete checks them.
INFO :: miComplete checked the 26 clusters in 0h 0m 3s, and found 0 of them contaminated, which the recursive steps partition again.
INFO :: Recursive step 1/9: no cluster is contaminated, so no recursive step is needed.
INFO :: 26 clusters kept after the recursive steps, 764377 bp in all:
INFO :: HQ MAGs: 0	Total Size: 0
INFO :: MQ MAGs: 0	Total Size: 0
INFO :: LQ MAGs: 0	Total Size: 0
INFO :: Contaminated potential MAGs: 0	Total Size: 0
INFO :: Other clusters: 26	Total Size: 764377
INFO :: Writing binning/sample.metacool from sample.metacool: 1219 contigs, 347 of them in 26 clusters.
INFO :: No list of MGE contigs (--mges): the run has no steps mge and host, which cluster the MGE contigs and find their hosts: give --mges to run them.
INFO :: [4/4] export: the outputs of the run, from its map
INFO :: The FASTA files of the 26 final cluster(s), with the MGE contigs that they hold: binning/fasta/mags_with-mges
INFO :: The FASTA files of the 26 MAG(s) of the run, each as its contigs and stitched: binning/fasta/mags
INFO :: The cool files of the 26 MAG(s) of the run, 0 of them with the MGE clusters that they host, balanced: binning/cool
INFO :: The tables of binning/sample.metacool are written in binning/tables: 1219 contigs, 27 genomes, 3 clustering rounds and 13415 edges of the network.
INFO :: The temporary directory of the run, ./tmp/7891e04f, and its files are kept (-N).
INFO :: Final contact map written after 0h 0m 5s: binning/sample.metacool
```

On a large sample the partition takes minutes: the log counts its iterations
as they go, and each recursive step of the validation (`Recursive step 2/9:
...`), and says when no further step is needed.

The Louvain algorithm is randomised: another run finds about as many clusters,
not always the same.

## Step 3: look at what it wrote

```console
$ ls binning
cool
fasta
metator_pipeline_20261123152335.log
plots
sample.metacool
tables
$ ls binning/fasta binning/tables
binning/fasta:
mags  mags_with-mges

binning/tables:
binning.txt  cluster_members.txt  cluster_summary.txt  clusters.txt
contig_data_final.txt  network.txt
```

The run stored everything else in the map of the clusters,
`binning/sample.metacool`: its network, its clusters and their summary. It
wrote its tables from the map, as `metator dump` does, in `binning/tables`;
the tables of its validation, the same, are in its temporary directory, which
`--no-clean-up` kept (`ls tmp/*/`; the log names it).

| File | What it holds |
|---|---|
| `sample.metacool` | The contact map of the clusters (see below) |
| `fasta/mags/` | Each MAG, two FASTA files: its contigs (`<MAG>.fa`) and stitched (`<MAG>.stitched.fa`, one record, its contigs joined by 200 `N`): the result |
| `fasta/mags_with-mges/` | One FASTA file per final cluster, as the validation made it (the same contigs as `fasta/mags/`, without `--mges`) |
| `cool/` | A cool file per MAG, balanced: its contacts, at the bins of the map |
| `plots/` | The plots of the clusters, as PDF files: their completeness and contamination (`bins_distribution.pdf`), and the size of each quality class (`bins_size_distribution.pdf`) |
| `metator_pipeline_<date>.log` | The log of the run |
| `tables/cluster_summary.txt` | One row per final cluster: size, GC content, completeness, redundancy, N50, Hi-C abundance… |
| `tables/clusters.txt`, `tables/cluster_members.txt` | Every cluster of every round (core, overlapping, recursive, final), with its parent, its number of contigs and their length, and miComplete's columns, and the contigs of each |
| `tables/contig_data_final.txt` | One row per contig: its size, GC content, hits, and the clusters it went through |
| `tables/binning.txt` | The contigs of the final clusters and their cluster, as anvi'o imports a binning |
| `tables/network.txt` | The network of the contigs, by their names |

`binning/sample.metacool` is the contact map of all the contigs, grouped by
cluster, with the network of the run: what metaTOR measured on each contig,
and the weights of the normalization of the network. `sample.metacool`, which
the run started from, is left as it was. [Running the pipeline](guides/pipeline.md#outputs) lists
every output, and [File formats](formats.md) describes the tables.

## Step 4: read the cluster summary

```console
$ cut -f 1-6,9,10 binning/tables/cluster_summary.txt | head -4
	Length	GC-content	Present Markers	Completeness	Redundancy	Contigs	N50
metator_00001_00000	36399.00	62.14	1.00	0.01	1.00	9.00	22786.00
metator_00006_00000	30476.00	61.77	1.00	0.01	1.00	9.00	7448.00
metator_00005_00000	20637.00	63.40	0.00	0.00	0.00	6.00	11969.00
```

Completeness and redundancy are miComplete's, as fractions: a complete
genome without contamination has a completeness and a redundancy of 1. The
clusters of this test library hold a few marker genes at most (a completeness
under 10%), so they are classed as *Other clusters*; on a real sample, the log
counts the high-, medium- and low-quality MAGs, and
`binning/plots/bins_size_distribution.pdf` shows their share of the assembly.

## Step 5: open the contact map

`binning/sample.metacool` is a [metacooler](https://github.com/js2264/metacooler)
file: a [cooler](https://github.com/open2c/cooler) whose chromosomes are the
contigs, grouped into one genome per final cluster. metacooler reads its
genomes:

```bash
metacooler info binning/sample.metacool
metacooler dump -t genomes -H binning/sample.metacool
```

The genomes table holds the columns of `cluster_summary.txt`, and the contigs
that are in no cluster are in the `__unbinned__` genome.
[The contact map of the clusters](guides/metacool.md) explains the file.

## The same from Python

```python
import glob

import metacooler
import pandas as pd

# The tables of the run, in its temporary directory, which --no-clean-up kept
(tables,) = glob.glob("tmp/*/")
clusters = pd.read_csv(f"{tables}/cluster_summary.txt", sep="\t", index_col=0)
contigs = pd.read_csv(f"{tables}/contig_data_final.txt", sep="\t", index_col=0, dtype={"Name": str})
contigs = contigs.set_index("Name", drop=False)

mc = metacooler.MetaCooler("binning/sample.metacool")
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
