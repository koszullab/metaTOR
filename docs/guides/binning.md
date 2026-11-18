# Partition and validation

The second and third steps of metaTOR turn the network into bins:
`metator partition` finds groups of contigs in the network, and
`metator validation` checks them with miComplete and splits the contaminated
ones. `metator pipeline` runs both. [Concepts](../concepts.md) explains the
principle; this page covers the parameters and the outputs.

![The steps and parameters of metaTOR's binning](../images/metator_figures_param.png)

## Partition

```bash
metator partition -a assembly.fa -c binning/contig_data_network.txt \
    -n binning/network.txt -t 16 -o binning
```

| Option | Default | Effect |
|---|---|---|
| `--algorithm` (`-A`) | `louvain` | The community detection algorithm: `louvain` or `leiden`. Louvain is faster, and gave better bins in our benchmarks on gut and seawater samples |
| `--iterations` (`-i`) | 100 | How many times the network is partitioned. More iterations split the contigs into more, smaller core bins; 100 is enough to find every bin in our benchmarks |
| `--overlap` (`-O`) | 80 | Two core bins are merged into one overlapping bin when they fall in the same community in at least this percentage of the iterations. A higher threshold gives less contaminated, but less complete, bins |
| `--res-param` (`-r`) | 1.0 | The resolution of the Leiden algorithm (Constant Potts Model): the higher, the more, smaller communities |
| `--size` (`-s`) | 500000 | The size, in base pairs, an overlapping bin must reach to be written |
| `--cluster-matrix` (`-C`) | off | Also write `clustering_matrix_partition.npz`: how often each two contigs fell in the same community |
| `--prefix` (`-p`) | none | Names the bins `<prefix>_metator_…` |

`partition` writes the overlapping bins of at least `--size` base pairs to
`overlapping_bin/`, one FASTA file each, and the contig data table with the
core and overlapping bin of each contig to `contig_data_partition.txt`. When
no bin reaches `--size`, there is nothing to validate: the pipeline stops
there with an error.

## Validation

```bash
metator validation -a assembly.fa -c binning/contig_data_partition.txt \
    -n binning/network.txt -f binning/overlapping_bin -t 16 -o binning
```

Validation runs [miComplete](https://pypi.org/project/micomplete/) on every
bin: prodigal predicts its genes, and hmmer finds among them the single-copy
marker genes of bacteria (105 markers) and of archaea (131 markers). Each bin
keeps the set that gives it the higher weighted completeness.

A bin is **contaminated** when its weighted completeness reaches 33% and its
contamination, (weighted redundancy − 1) / weighted completeness, reaches 5%.
Each contaminated bin is partitioned again, on the network of its own
contigs:

| Option | Default | Effect |
|---|---|---|
| `--iterations` (`-i`) of `validation`, `--rec-iter` (`-j`) of `pipeline` | 10 | How many times each contaminated bin is partitioned |
| `--overlap` (`-O`) of `validation`, `--rec-overlap` (`-P`) of `pipeline` | 90 | The overlap threshold of the recursive partition |
| `--algorithm`, `--res-param`, `--size`, `--prefix` | | As in the partition |
| `--cluster-matrix` (`-C`) | off | Also write `clustering_matrix_recursive.npz` (and, in `pipeline`, `clustering_matrix.npz`, the mean of both) |

The parts of a bin, its **recursive bins**, are checked with miComplete in
turn. They replace the bin when the most complete of them keeps enough of its
completeness: more than 0.4 + (c − 0.4) / 1.6, for a bin of completeness c.
Otherwise the bin is kept whole. The step repeats on the recursive bins that
are still contaminated, until none is replaced.

## Outputs

| File | What it holds |
|---|---|
| `final_bin/` (`final_bin_unscaffold/` in `pipeline`) | One FASTA file per final bin |
| `recursive_bin/` | One FASTA file per recursive bin |
| `bin_summary.txt` | One row per final bin, with miComplete's measures and the bin's abundance |
| `contig_data_final.txt` | The contig data table, with the recursive and final bin of each contig |
| `binning.txt` | Each binned contig and its final bin, without header |
| `overlapping_micomplete_results.txt` | miComplete's results on the overlapping bins |
| `recursive_micomplete_step_<n>.txt` | miComplete's results on the recursive bins of step `n` |
| `plot/` | The plots of the bins (below) |

A final bin is named after its overlapping bin and its recursive bin:
`metator_00012_00003` is recursive bin 3 of overlapping bin 12, and
`metator_00001_00000` an overlapping bin kept whole. `bin_summary.txt` and
`contig_data_final.txt` are described in [File formats](../formats.md).

### Plots

| Plot | Shows |
|---|---|
| `bins_distribution.png` | The completeness and contamination of the bins |
| `bins_size_distribution.png` | The share of the assembly in each quality class, and unbinned |
| `MAGs_GC_distribution.png`, `MAGs_HiC_cov_distribution.png` | The GC content and Hi-C coverage of the contigs of each MAG |
| `MAGs_SG_cov_distribution.png` | The shotgun coverage of the contigs of each MAG, with `--depth` |
| `network_heatmap.png` | The contacts between the bins, from the pairs files of the output directory |

|Completeness and contamination|Quality classes|
|:-:|:-:|
|![bins_distribution](../images/bins_distribution.png)|![bins_size_distribution](../images/bins_size_distribution.png)|

The heatmap shows the contacts within each final bin (on the diagonal) and
between them, at 50 kb resolution, the bins in the order of
`bin_summary.txt`: a good binning has a bright diagonal of blocks on a dark
background. It is drawn from the `.pairs.gz` files of the output directory,
and skipped, with a warning, when there are none.

![network_heatmap](../images/network_heatmap.png)

## What comes next

- [Curating MAGs with anvi'o](curation.md): remove what contamination is left,
  by hand.
- [Viewing and scaffolding MAGs](visualisation.md): the contact map of a bin,
  and its scaffold.
- [Binning mobile genetic elements](mge.md): the plasmids and phages, and
  their hosts.
