# Partition and validation

The second and third steps of metaTOR turn the network into clusters:
`metator partition` finds groups of contigs in the network, and
`metator validation` checks them with miComplete and splits the contaminated
ones. `metator pipeline` runs both. [Concepts](../concepts.md) explains the
principle; this page covers the parameters and the outputs.

![The steps and parameters of metaTOR's clustering](../images/metator_figures_param.png)

## Partition

```bash
metator partition -a assembly.fa -c binning/contig_data_network.txt \
    -n binning/network.txt -t 16 -o binning
```

| Option | Default | Effect |
|---|---|---|
| `--algorithm` (`-A`) | `louvain` | The community detection algorithm: `louvain` or `leiden`. Louvain is faster, and gave better clusters in our benchmarks on gut and seawater samples |
| `--iterations` (`-i`) | 100 | How many times the network is partitioned. More iterations split the contigs into more, smaller core clusters; 100 is enough to find every cluster in our benchmarks |
| `--overlap` (`-O`) | 80 | Two core clusters are merged into one overlapping cluster when they fall in the same community in at least this percentage of the iterations. A higher threshold gives less contaminated, but less complete, clusters |
| `--res-param` (`-r`) | 1.0 | The resolution of the modularity that the Leiden algorithm optimises (1 is plain modularity, as in Louvain): the higher, the more, smaller communities |
| `--size` (`-s`) | 500000 | The size, in base pairs, an overlapping cluster must reach to be written |
| `--cluster-matrix` (`-C`) | off | Also write `clustering_matrix_partition.npz`: how often each two contigs fell in the same community |
| `--prefix` (`-p`) | none | Names the clusters `<prefix>_metator_…` |

`partition` writes the overlapping clusters of at least `--size` base pairs to
`overlapping_clusters/`, one FASTA file each, and the contig data table with the
core and overlapping cluster of each contig to `contig_data_partition.txt`. When
no cluster reaches `--size`, there is nothing to validate: the pipeline stops
there with an error.

## Validation

```bash
metator validation -a assembly.fa -c binning/contig_data_partition.txt \
    -n binning/network.txt -f binning/overlapping_clusters -t 16 -o binning
```

Validation runs [miComplete](https://pypi.org/project/micomplete/) on every
cluster: Prodigal predicts its genes, and HMMER finds among them the
single-copy marker genes of bacteria (105 markers) and of archaea
(131 markers), through pyrodigal and pyhmmer, their Python bindings. Each
cluster keeps the set that gives it the higher weighted completeness.

A cluster is **contaminated** when its weighted completeness reaches 33% and its
contamination, (weighted redundancy − 1) / weighted completeness, reaches 5%.
Each contaminated cluster is partitioned again, on the network of its own
contigs:

| Option | Default | Effect |
|---|---|---|
| `--iterations` (`-i`) of `validation`, `--rec-iter` (`-j`) of `pipeline` | 10 | How many times each contaminated cluster is partitioned |
| `--overlap` (`-O`) of `validation`, `--rec-overlap` (`-P`) of `pipeline` | 90 | The overlap threshold of the recursive partition |
| `--algorithm`, `--res-param`, `--size`, `--prefix` | | As in the partition |
| `--cluster-matrix` (`-C`) | off | Also write `clustering_matrix_recursive.npz`: how often each two contigs of a cluster partitioned again fell together, in its latest partition (and, in `pipeline`, `clustering_matrix.npz`, the mean of both matrices) |

The parts of a cluster, its **recursive clusters**, are checked with
miComplete in turn. They replace the cluster when the most complete of them
keeps enough of its completeness: more than 0.4 + (c − 0.4) / 1.6, for a
cluster of completeness c. Otherwise the cluster is kept whole. The step
repeats on the recursive clusters that are still contaminated, until none is
replaced, for 9 steps at most.

A part is a recursive cluster only when it is larger than `--size`: when a
cluster is replaced by its recursive clusters, the contigs of its smaller
parts are left out of the final clusters (`ND`).

## Outputs

| File | What it holds |
|---|---|
| `final_clusters/` | One FASTA file per final cluster |
| `recursive_clusters/` | One FASTA file per recursive cluster, in `step_<n>/` for step `n` |
| `cluster_summary.txt` | One row per final cluster, with miComplete's measures and the cluster's abundance |
| `contig_data_final.txt` | The contig data table, with the recursive and final cluster of each contig |
| `binning.txt` | The contigs of the final clusters and their cluster, without header |
| `overlapping_micomplete_results.txt` | miComplete's results on the overlapping clusters |
| `recursive_micomplete_step_<n>.txt` | miComplete's results on the recursive clusters of step `n` |

A final cluster is named after its overlapping cluster and its recursive
cluster: `metator_00012_00003` is recursive cluster 3 of overlapping
cluster 12, and `metator_00001_00000` an overlapping cluster kept whole.
`cluster_summary.txt` and `contig_data_final.txt` are described in
[File formats](../formats.md).

### Viewing the clusters

metaTOR draws no plot since 2.0. The log of the validation counts the
clusters in each quality class (HQ, MQ and LQ MAGs, contaminated clusters and
the others), and `cluster_summary.txt` gives the completeness, contamination
and size of each one. `metacooler show` draws the contacts between the final
clusters from the contact map of the run (see
[Other tools](../other-tools.md)):

```bash
metacooler show binning/metator.metacool -o clusters.png
```

## What comes next

- [Curating MAGs with anvi'o](curation.md): remove what contamination is left,
  by hand.
- [Other tools](../other-tools.md): the contact map of a cluster, and its
  scaffold, with metacooler, hicstuff and instaGRAAL.
- [Clustering mobile genetic elements](mge.md): the plasmids and phages, and
  their hosts.
