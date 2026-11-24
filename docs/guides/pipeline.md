# Running the pipeline

`metator pipeline` bins the contigs of an assembly from a Hi-C contact map in
one command: it builds the network of contigs, partitions it, validates the
clusters, writes the contact map of the clusters and stores in it the marker
genes of each contig. With `--mges`, it also clusters the mobile genetic
elements and finds their hosts. The
[CLI reference](../cli.md#pipeline) lists every option, and
[Concepts](../concepts.md) explains the steps.

```bash
metator pipeline sample.metacool --enzyme DpnII --depth depth.txt --threads 16 \
    --assembly assembly.fa --outdir binning --mges mges.txt
```

## Before the run: the assembly

metaTOR does not assemble. [The end-to-end pipeline](../end-to-end.md) shows
how to get the assembly (MEGAHIT or metaSPAdes, on a shotgun library of the
same sample), the contact map and the list of MGE contigs. Drop contigs
shorter than 500 bp: they carry too few contacts to cluster.

`--depth` is optional: it is the output of MetaBAT2's
`jgi_summarize_bam_contig_depths` on the shotgun reads aligned back to the
assembly. The default normalisation does not use it.

## Inputs

The first argument is a `.metacool` contact map of the contigs (`MAP`): a
file, or an entry of a `.mmetacool` collection
(`file.mmetacool::entries/<name>`).

| Starting from | How to get the map |
|---|---|
| Reads | `hicstuff pipeline -g assembly.fa -e DpnII -M metacool -P sample R1.fq.gz R2.fq.gz` writes `sample.metacool` (hicstuff 4.0.0 or later) |
| Pairs files | `metacooler create pairs --pairs sample.pairs.gz --bins 1000 -o sample.metacool`, or give the pairs files to `metator pipeline` (below) |
| BAM files | Give them to `metator pipeline` (below) |
| Several libraries | `metacooler merge L1.metacool L2.metacool -o sample.metacool` sums their maps; or give the pairs or BAM files of every library to `metator pipeline` |

See [Other tools](../other-tools.md#from-the-reads-to-the-map) for these
commands, and [Starting from a contact map](metacool.md#starting-from-a-contact-map)
for what the run reads from the map.

`--assembly` (`-a`) is the FASTA file the map was made from, plain or
gzipped.

### From pairs or BAM files

Given pairs or BAM files instead of a map, the run first makes the map in its
output directory (one bin per contig, libraries summed), then clusters it.
Each argument is a file, a comma-separated list or a glob pattern.

```bash
# Pairs files of hicstuff, one per library: binning/sample.metacool
metator pipeline 'hicstuff_*/*.valid_idx.pairs.gz' -p sample --enzyme DpnII --depth depth.txt \
    --threads 16 --assembly assembly.fa --outdir binning

# BAM files with both mates of each read (bwa mem -5SP), one per library: binning/assembly.metacool
metator pipeline L1.bam,L2.bam --enzyme DpnII --threads 16 --assembly assembly.fa --outdir binning

# BAM files with one mate each (bowtie2): the R1 files, then the R2 files
metator pipeline 'L*_R1.bam' 'L*_R2.bam' --enzyme DpnII --threads 16 --assembly assembly.fa --outdir binning
```

- The map is named after `--prefix` (`sample.metacool`), or else after the
  assembly (`assembly.metacool` from `assembly.fa.gz`).
- From BAM files, a pair counts when both mates are primary alignments with a
  mapping quality of at least `--min-quality` (`-q`, 30). Reads flagged as
  duplicates or as failing QC are left out; metaTOR flags none itself. See
  [The inputs of a run](../cli.md#inputs) for the details.
- The pairs of BAM files are written to `--tmpdir` first (about 25 bytes per
  pair): give it room.
- A run that stops after making the map leaves it in the output directory.
  Restart from it: `metator pipeline binning/sample.metacool ... --outdir binning`.

## Outputs

For `metator pipeline sample.metacool --outdir binning`, the run writes into
`binning/`:

| File | What it holds |
|---|---|
| `sample.metacool` | The contact map of the clusters (from pairs or BAM files, the map the run made) |
| `fasta/mags_with-mges/` | One FASTA file per final cluster, as the validation made it, MGE contigs included |
| `fasta/mags/` | Each MAG twice: its contigs (`<MAG>.fa`) and stitched into one record, contigs joined by 200 `N` (`<MAG>.stitched.fa`). With `--mges`, without the MGE contigs |
| `fasta/mges/` | With `--mges`: each MGE cluster, the same two ways (see [Outputs](mge.md#outputs) of `mge`) |
| `cool/` | One balanced cool file per MAG, `<MAG>.cool`, with the MGE clusters it hosts |
| `tables/` | The tables of [`metator dump`](../cli.md#dump): `cluster_summary.txt`, `contig_data_final.txt`, `binning.txt`... With `--mges`, `mge_hosts.txt` (see [What host writes](mge.md#what-host-writes)) |
| `plots/` | `bins_distribution.pdf` and `bins_size_distribution.pdf` (see [Viewing the clusters](binning.md#viewing-the-clusters)); with `--mges`, `mge_hosts.pdf` |
| `sample.report.html` | The report of the run, to open in a browser: its numbers, a table of the MAGs with links to their files, the MGE clusters and their hosts, and interactive plots (see [`report`](../cli.md#report)) |
| `sample.report.txt`, `sample.report_mqc.yaml` | The same report as text for scripts, and as a row of the General Statistics table of [MultiQC](https://multiqc.info) (see [The report of a run](../formats.md#the-report-of-a-run)) |
| `sample.metator_<date>.log` | The log of the run |

`sample` is the name of the map without `.metacool`.
[File formats](../formats.md) describes each file.

The map holds everything the run found: the network, every clustering round
and the summary of the final clusters (see
[The contact map of the clusters](metacool.md)). `MAP` itself is left as it
was, unless you give `--inplace` or the output directory already holds it.
The last step, `export`, writes the folders and the report from the map;
[`metator export`](../cli.md#export) writes them again at any time.

The validation also writes its working tables (miComplete's results, the
FASTA files of the clusters it checked) in `--tmpdir`, deleted at the end
unless you give `--no-clean-up` (`-N`).

A run does not overwrite the outputs of an earlier run without `--force`
(`-F`), and refuses the outputs of older metaTOR versions even with it: see
[`pipeline`](../cli.md#pipeline).

## Choosing the options

| Option | Default | When to change it |
|---|---|---|
| `--enzyme` | none | The enzymes of the library, for the `RS` and `theoritical_hit` normalisations. Not needed if the map has a `Restriction_site` column |
| `--depth` | none | For the `abundance` and `theoritical_hit` normalisations, and the abundance of the clusters. Not needed if the map has a `Shotgun_coverage` column |
| `--normalization` | `empirical_hit` | See [Normalising the network](network.md#normalising-the-network) |
| `--iterations`, `--overlap` | 100, 80 | See [Partition](binning.md#partition) |
| `--rec-iter`, `--rec-overlap` | 10, 90 | See [Validation](binning.md#validation) |
| `--size` | 500000 | The size, in bp, an overlapping cluster must reach to be checked with miComplete, and a recursive cluster must exceed |
| `--mges` | none | The MGE contigs, one name per line: adds the steps `mge` and `host` (see [below](#mobile-genetic-elements-and-their-hosts)) |
| `--mge-threshold-bin` | 0.8 | Step `mge`: the score above which MGE contigs are merged (see [How it works](mge.md#how-it-works)) |
| `--mge-min-share`, `--mge-min-contigs` | 10, 5 | Step `host`: the share of signal (%) and the number of contacting contigs a MAG needs to be a host (see [The hosts](mge.md#the-hosts)) |
| `--prefix` | none | Names the clusters `<prefix>_metator_…` |
| `--threads` | 1 | Threads for the partition, the validation and the markers |
| `--tmpdir` | `./tmp` | Temporary files. Give it room for the FASTA files of the checked clusters (at most the size of the assembly) and the decompressed copy of a gzipped assembly |

The validation needs the most memory and time: miComplete runs Prodigal and
HMMER on every cluster.

## The marker genes of the contigs

After the validation, the step `markers` (as
[`metator markers`](../cli.md#markers)) finds miComplete's marker genes on
each contig, in its final cluster, as miComplete finds them, and stores them
in the map of the clusters. Summed over any group of contigs, they give the
completeness and redundancy that miComplete would give it: a viewer, or a
script, recomputes them as a curation moves contigs between MAGs, without
finding a gene again (see
[Marker genes of the contigs](../formats.md#marker-genes-of-the-contigs)).

## Mobile genetic elements and their hosts

`--mges` takes the contigs of plasmids and phages, found beforehand with an
annotation tool such as geNomad, one name per line. After the validation, the
run adds two steps on the map of the clusters:

| Step | Same as | What it does |
|---|---|---|
| `mge` | `metator mge MAP --reassign` | Clusters the MGE contigs and makes each MGE cluster a genome of the map, taking its contigs out of the MAGs; writes `fasta/mges/` |
| `host` | `metator host MAP` | Finds the MAGs that host each MGE cluster; writes `tables/mge_hosts.txt` and `plots/mge_hosts.pdf` |

So `fasta/mags/` holds the MAGs without their MGE contigs (see
[MGE genomes](mge.md#mge-genomes)), `fasta/mags_with-mges/` the final
clusters as the validation made them, and each file of `cool/` includes the
MGE clusters its MAG hosts.

The run checks the MGE inputs before any work: the list must be readable, at
least one of its contigs must be in the map, and the map must hold the
contact statistics that `mge` needs (`hicstuff pipeline -M metacool`,
`metacooler create pairs` from metacooler 0.4.0, and a run from pairs or BAM
files write them).

If the run stops in `markers`, `mge` or `host`, its log says what to run to
finish without redoing the validation. The same commands add the MGE steps to a run
made without `--mges`:

```bash
metator mge binning/sample.metacool -a assembly.fa -m mges.txt --reassign -o binning
metator host binning/sample.metacool -o binning
metator export binning/sample.metacool -a assembly.fa -o binning -F
```

## Running the steps one by one

Each step of the pipeline is a command of its own. Run them one by one to
rerun the partition with other parameters without rebuilding the network:

```bash
# 1. Build the network and store it in binning/sample.metacool
metator network sample.metacool -a assembly.fa -e DpnII -d depth.txt -o binning

# 2. Partition it: core and overlapping clusters, stored in the map
metator partition binning/sample.metacool -t 16 -o binning

# 3. Validate the clusters: fasta/mags_with-mges/, plots/, and the map of the clusters
metator validation binning/sample.metacool -a assembly.fa -t 16 -o binning

# 4. Find the marker genes of the contigs, stored in the map
metator markers binning/sample.metacool -a assembly.fa -t 16 -o binning

# 5. (MGEs) Cluster the MGE contigs: fasta/mges/
metator mge binning/sample.metacool -a assembly.fa -m mges.txt --reassign -o binning

# 6. (MGEs) Find their hosts: tables/mge_hosts.txt, plots/mge_hosts.pdf
metator host binning/sample.metacool -o binning

# 7. Write the other outputs: fasta/mags/, cool/, tables/ and the report
metator export binning/sample.metacool -a assembly.fa -o binning -F
```

The map is what the steps share: each one reads and updates
`binning/sample.metacool` in place, and together they write the same map as
`metator pipeline` (see [Step by step](metacool.md#step-by-step)).

To cluster without validating, stop after step 2. The overlapping clusters
are then the round `overlapping` of the map;
`metacooler fasta MAP --fasta assembly.fa --round overlapping -o DIR` writes
them as FASTA files.

## The log

The log, `sample.metator_<date>.log` in the output directory, starts with
the version, the date, the command line and the options. It then marks each
step: `[1/5] network`, `[2/5] partition`, `[3/5] validation`, `[4/5] markers`,
`[5/5] export` (`[1/6] map` first from pairs or BAM files; `[5/7] mge` and
`[6/7] host` with `--mges`). It reports the contacts, the progress of the partition, each
recursive step of the validation, the number of clusters in each
[quality class](../concepts.md#cluster-quality), and with `--mges` the MGE
clusters and their hosts.

```console
INFO :: [1/5] network: the contigs of sample.metacool and their contacts
INFO :: 41563 contacts in sample.metacool, 18632 of them between contigs (3D ratio: 0.4483).
...
INFO :: [5/5] export: the outputs of the run, from its map
...
INFO :: Final contact map written after 0h 0m 9s: binning/sample.metacool
```

The **3D ratio** is the share of contacts that join two different contigs:
the signal that the clustering uses.
