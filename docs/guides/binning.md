# Partition and validation

The second and third steps of metaTOR turn the network into clusters.
`metator partition` finds groups of contigs in the network stored in the map,
and `metator validation` checks them with miComplete and splits the
contaminated ones. `metator pipeline` runs both. [Concepts](../concepts.md)
explains the principle; this page covers the parameters and the outputs.

![The steps and parameters of metaTOR's clustering](../images/metator_figures_param.png)

## Partition

```bash
# The map that metator network sample.metacool -o binning wrote, updated in place
metator partition binning/sample.metacool -t 16 -o binning
```

| Option | Default | Effect |
|---|---|---|
| `--algorithm` (`-A`) | `louvain` | `louvain` or `leiden`. Louvain is faster, and gave better clusters in our benchmarks on gut and seawater samples |
| `--iterations` (`-i`) | 100 | How many times the network is partitioned. More iterations give more, smaller core clusters; 100 found every cluster in our benchmarks |
| `--overlap` (`-O`) | 80 | Two core clusters merge into one overlapping cluster when they fall together in at least this percentage of the iterations. Higher gives less contaminated but less complete clusters |
| `--res-param` (`-r`) | 1.0 | The resolution of Leiden (1 is plain modularity, as in Louvain): higher gives more, smaller communities |
| `--prefix` (`-p`) | none | Names the overlapping clusters `<prefix>_metator_…` |

`partition` stores the core clusters (`core_00001`…) and the overlapping
clusters (`metator_00001_00000`…) in the map, as the rounds `core` and
`overlapping`. It writes no table and no FASTA file.

## Validation

```bash
metator validation binning/sample.metacool -a assembly.fa -t 16 -o binning
```

The validation checks the overlapping clusters of at least `--size` bp
(500 kb). It writes them as FASTA files in its temporary directory: give it
room for up to the size of the assembly, plus the decompressed copy of a
gzipped assembly. If no cluster reaches `--size`, the run stops with an
error. If the network was rebuilt since the partition, run `partition` again
first.

[miComplete](https://pypi.org/project/micomplete/) checks each cluster: it
predicts genes with Prodigal and finds the single-copy marker genes of
bacteria (105) and archaea (131) with HMMER (through pyrodigal and pyhmmer).
Each cluster keeps the marker set that gives it the higher weighted
completeness.

A cluster is **contaminated** when its weighted completeness is at least 33%
and its contamination, (weighted redundancy − 1) / weighted completeness, is
at least 5%. Each contaminated cluster is partitioned again, on the network
of its own contigs:

| Option | Default | Effect |
|---|---|---|
| `--iterations` (`-i`) of `validation`, `--rec-iter` (`-j`) of `pipeline` | 10 | How many times each contaminated cluster is partitioned |
| `--overlap` (`-O`) of `validation`, `--rec-overlap` (`-P`) of `pipeline` | 90 | The overlap threshold of the recursive partition |
| `--algorithm`, `--res-param`, `--prefix` | | As in the partition |
| `--size` (`-s`) | 500000 | The size, in bp, an overlapping cluster must reach to be checked, and a recursive cluster must exceed |

The parts larger than `--size`, the **recursive clusters**, are checked in
turn. They replace the cluster when the most complete of them keeps enough
of its completeness: more than 0.4 + (c − 0.4) / 1.6 for a cluster of
completeness c. Otherwise the cluster is kept whole. This repeats on the
recursive clusters that are still contaminated, for at most 9 steps. When a
cluster is replaced, the contigs of its smaller parts are left out of the
final clusters (`ND`).

## Outputs

`validation` and `pipeline` write into the output directory:

- `fasta/mags_with-mges/`: one FASTA file per final cluster, MGE contigs
  included;
- `plots/`: the plots of the final clusters (see
  [Viewing the clusters](#viewing-the-clusters));
- the map of the clusters (see [The contact map of the clusters](metacool.md)),
  and the log, `sample.metator_<date>.log`.

Their other files go to the temporary directory of the run (`-T`), deleted
at the end unless `-N`: `recursive_clusters/` (FASTA files, in `step_<n>/`),
`cluster_summary.txt`, `contig_data_final.txt`, `binning.txt`, and
miComplete's results (`overlapping_micomplete_results.txt`,
`recursive_micomplete_step_<n>.txt`). [`metator dump`](../cli.md#dump)
writes the tables from the map at any time; [File formats](../formats.md)
describes them.

A final cluster is named after its overlapping and recursive clusters:
`metator_00012_00003` is recursive cluster 3 of overlapping cluster 12, and
`metator_00001_00000` an overlapping cluster kept whole.

### Viewing the clusters

The log counts the final clusters in each quality class (HQ, MQ and LQ MAGs,
contaminated clusters and others: see
[Cluster quality](../concepts.md#cluster-quality)), and the run draws two
plots in `plots/`:

| Plot | Shows |
|---|---|
| `bins_distribution.pdf` | The completeness (red) and contamination (black) of each final cluster, by decreasing completeness, with dashed lines at the class thresholds. Contamination is capped at 105% |
| `bins_size_distribution.pdf` | The share of the contigs of the run in each quality class and in no cluster, with the number and size of the clusters of each class |

`metacooler show` draws the contacts between the final clusters (see
[Other tools](../other-tools.md)):

```bash
metacooler show binning/sample.metacool -o clusters.png
```

## What comes next

- [Curating MAGs with anvi'o](curation.md): remove the contamination that is
  left, by hand.
- [Other tools](../other-tools.md): the contact map of a cluster, and its
  scaffold, with metacooler, hicstuff and instaGRAAL.
- [Clustering mobile genetic elements](mge.md): the plasmids and phages, and
  their hosts.
