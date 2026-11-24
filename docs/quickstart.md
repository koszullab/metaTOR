# Quickstart

This page clusters the contigs of the test data of metaTOR's repository:
41,563 read pairs on an assembly of 1,219 contigs (2.2 Mb). The run takes
under a minute. For a real sample, see
[The end-to-end pipeline](end-to-end.md).

## Input files

| File | What it is |
|---|---|
| `assembly.fa` | The metagenome assembly (FASTA) |
| `alignment.pairs` | The Hi-C read pairs aligned on the assembly ([pairs format](formats.md#pairs)) |
| `depth.txt` | *(optional)* The shotgun coverage of each contig, from `jgi_summarize_bam_contig_depths` |

```bash
base=https://raw.githubusercontent.com/koszullab/metator/master/tests_data
curl -LO "$base/assembly.fa"
curl -LO "$base/alignment.pairs"
curl -LO "$base/depth.txt"
```

## Step 1: make the contact map

metaTOR starts from a contact map of the contigs, a `.metacool` file.
metacooler, installed with metaTOR, makes it from the pairs:

```bash
metacooler create pairs --pairs alignment.pairs --bins 1000 --no-ps \
    --assembly-name assembly -o sample.metacool
```

This cuts the contigs into 1 kb bins, all in the genome `__unbinned__`.
metacooler warns that it clamped 719 pairs: reads at position 0, counted in
the first bin of their contig.

!!! tip "From the reads"

    hicstuff 4.0.0 or later makes the map from the reads, one bin per
    restriction fragment (see [Other tools](other-tools.md#from-the-reads-to-the-map)):

    ```bash
    hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff -P sample R1.fq.gz R2.fq.gz
    ```

!!! tip "Or let the pipeline make the map"

    `metator pipeline` and `metator network` also take pairs or BAM files,
    and make the map themselves, one bin per contig
    (see [The inputs of a run](cli.md#inputs)):

    ```bash
    metator pipeline alignment.pairs --assembly assembly.fa --enzyme DpnII \
        --depth depth.txt --size 20000 --threads 4 --outdir binning --no-clean-up
    ```

    The map is then `binning/assembly.metacool`, or `binning/sample.metacool`
    with `--prefix sample`.

## Step 2: run the pipeline

```bash
metator pipeline sample.metacool --assembly assembly.fa --enzyme DpnII \
    --depth depth.txt --size 20000 --threads 4 --outdir binning --no-clean-up
```

- `--enzyme` and `--depth` add restriction fragments and shotgun coverage to
  the contig data. The default normalisation needs neither.
- The contact map of the clusters goes to `binning/sample.metacool`;
  `sample.metacool` is left unchanged (`--inplace` updates it instead).
- `--no-clean-up` keeps the temporary directory, in `tmp/`.

!!! note "`--size 20000` is for this test library only"

    miComplete checks only the clusters of at least `--size` bp (500 kb by
    default). At the default, this small library has none, and the run stops
    before validation. On a real sample, leave `--size` out.

The log follows the steps of the run (excerpt):

```console
INFO :: ## MetaTOR: v1.4.0 log file
INFO :: ## cmd: metator pipeline sample.metacool --assembly assembly.fa --enzyme DpnII --depth depth.txt --size 20000 --threads 4 --outdir binning --no-clean-up
...
INFO :: [1/5] network: the contigs of sample.metacool and their contacts
INFO :: [2/5] partition: the clusters of the network
INFO :: 420 core clusters, merged into 231 overlapping clusters.
INFO :: [3/5] validation: miComplete checks the clusters, and the contaminated ones are partitioned again
INFO :: 24 of the 231 overlapping clusters reach 20000 bp, 0.72 Mb in all: miComplete checks them.
...
INFO :: Other clusters: 24	Total Size: 720373
INFO :: [4/5] markers: the marker genes of the contigs, as miComplete finds them in their clusters
...
INFO :: [5/5] export: the outputs of the run, from its map
...
INFO :: The report of the run, from its map: binning/sample.report.txt, binning/sample.report.html, binning/sample.report_mqc.yaml
INFO :: Final contact map written after 0h 0m 10s: binning/sample.metacool
```

Without `--mges`, the run skips the `mge` and `host` steps (see
[Running the pipeline](guides/pipeline.md#mobile-genetic-elements-and-their-hosts)).
On this tiny library, `export` also warns that most cool files cannot be
balanced. Louvain is randomised: another run finds about as many clusters,
not always the same ones.

## Step 3: look at the outputs

```console
$ ls binning
cool
fasta
plots
sample.metacool
sample.metator_20261123152335.log
sample.report.html
sample.report.txt
sample.report_mqc.yaml
tables
```

| File | What it holds |
|---|---|
| `sample.metacool` | The contact map of the clusters (see [step 5](#step-5-open-the-contact-map)) |
| `fasta/mags/` | The MAGs: each as its contigs (`<MAG>.fa`) and stitched (`<MAG>.stitched.fa`, contigs joined by 200 `N`) |
| `fasta/mags_with-mges/` | One FASTA file per final cluster (the same contigs as `fasta/mags/` without `--mges`) |
| `cool/` | One balanced cool file per MAG |
| `plots/` | Completeness and contamination (`bins_distribution.pdf`), size of each quality class (`bins_size_distribution.pdf`) |
| `sample.report.html` | The report, for a browser: summary, a table of the MAGs, interactive plots |
| `sample.report.txt`, `sample.report_mqc.yaml` | The same report as text, and for [MultiQC](https://multiqc.info) |
| `sample.metator_<date>.log` | The log |
| `tables/cluster_summary.txt` | One row per final cluster: size, GC content, completeness, redundancy, N50, Hi-C abundance… |
| `tables/clusters.txt`, `tables/cluster_members.txt` | Every cluster of every round (core, overlapping, recursive, final), and its contigs |
| `tables/contig_data_final.txt` | One row per contig: size, GC content, hits, clusters |
| `tables/binning.txt` | Each contig of a final cluster and its cluster, for anvi'o |
| `tables/network.txt` | The network of the contigs |

`metator dump` writes the same tables from the map.
[Running the pipeline](guides/pipeline.md#outputs) lists every output, and
[File formats](formats.md) describes them.

## Step 4: read the cluster summary

```console
$ cut -f 1-6,9,10 binning/tables/cluster_summary.txt | head -4
	Length	GC-content	Present Markers	Completeness	Redundancy	Contigs	N50
metator_00005_00000	20637.00	63.40	0.00	0.00	0.00	6.00	11969.00
metator_00007_00000	40709.00	62.29	4.00	0.04	1.00	17.00	3159.00
metator_00006_00000	30476.00	61.77	1.00	0.01	1.00	9.00	7448.00
```

Completeness and redundancy come from miComplete, as fractions: 1 and 1 for
a complete, clean genome. These test clusters hold few marker genes, so they
are all *Other clusters* (see [Cluster quality](concepts.md#cluster-quality)).
`binning/sample.report.html` shows the run at a glance (see
[The report of a run](formats.md#the-report-of-a-run)).

## Step 5: open the contact map

`binning/sample.metacool` is a [metacooler](https://github.com/js2264/metacooler)
file: a [cooler](https://github.com/open2c/cooler) whose contigs are grouped
into one genome per final cluster, plus `__unbinned__` for the rest.

```bash
metacooler info binning/sample.metacool
metacooler dump -t genomes -H binning/sample.metacool
```

The genomes table holds the columns of `cluster_summary.txt`. See
[The contact map of the clusters](guides/metacool.md).

## The same from Python

```python
import metacooler
import pandas as pd

clusters = pd.read_csv("binning/tables/cluster_summary.txt", sep="\t", index_col=0)
contigs = pd.read_csv("binning/tables/contig_data_final.txt", sep="\t", index_col="Name", dtype={"Name": str})

mc = metacooler.MetaCooler("binning/sample.metacool")
mc.genomes   # one row per cluster, and __unbinned__, with the columns of cluster_summary.txt
mc.contigs   # one row per contig: GC content, hits, coverage, restriction sites
```

## Next steps

- [The end-to-end pipeline](end-to-end.md): from the reads to MAGs, MGE
  clusters and their hosts.
- [Concepts](concepts.md): how metaTOR clusters contigs.
- [Running the pipeline](guides/pipeline.md): every input, output and option.
- [Other tools](other-tools.md): what hicstuff, metacooler and instaGRAAL do
  for metaTOR.
- [CLI reference](cli.md): every command.
