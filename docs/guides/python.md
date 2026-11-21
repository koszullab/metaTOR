# Using the Python API

Each command of metaTOR runs functions of its library, which can be called
from Python to run a step, read the outputs, or check other clusters.
This page shows the pieces to put together; the
[Python API reference](../api/network.md) lists them all.

## The modules

| Module | Usual name | What it holds |
|---|---|---|
| `metator.network` | `mtn` | The contig data table, and the network built from a contact map |
| `metator.partition` | `mtp` | Partitioning the network into core and overlapping clusters |
| `metator.validation` | `mtv` | miComplete, and the recursive partition of contaminated clusters |
| `metator.metacool` | `mtmc` | Writing the `.metacool` contact map of the clusters |
| `metator.mge`, `metator.host` | `mtm` | Clustering MGE contigs, and finding their hosts |
| `metator.io` | `mio` | Reading and writing tables and FASTA files |
| `metator.log` | `mtl` | The logger, `metator.log.logger` |

## Reading the outputs

pandas reads `cluster_summary.txt`, indexed by cluster name, and the contig
data table, indexed here by contig name:

```python
import pandas as pd

summary = pd.read_csv("binning/cluster_summary.txt", sep="\t", index_col=0)
contigs = pd.read_csv("binning/contig_data_final.txt", sep="\t", index_col=0).set_index("Name", drop=False)

# The high-quality MAGs, and their contigs.
hq = summary[(summary["Weighted completeness"] >= 0.9) & (summary["Weighted redundancy"] <= 1.05)]
hq_contigs = contigs[contigs["Final_cluster"].isin(hq.index)]
```

[File formats](../formats.md) describes their columns.

## The contact map of the clusters

The `.metacool` file is read with
[metacooler](https://github.com/js2264/metacooler), or as any cool file with
cooler (see [The contact map of the clusters](metacool.md)):

```python
import metacooler

mc = metacooler.MetaCooler("sample.metacool")
mc.genomes  # one row per cluster, and __unbinned__, with the columns of cluster_summary.txt
mc.contigs  # one row per contig
```

## Running a step

The functions behind the commands take their options as arguments. The
partition of a network, as `metator partition` runs it:

```python
import os

import metator.partition as mtp

os.makedirs("binning/overlapping_clusters", exist_ok=True)
contig_data_file = mtp.partition(
    algorithm="louvain",
    assembly="assembly.fa",
    contig_data_file="binning/contig_data_network.txt",
    iterations=100,
    network_file="binning/network.txt",
    outdir="binning",
    fasta_dir="binning/overlapping_clusters",
    overlapping_parameter=0.8,  # --overlap 80
    resolution_parameter=1.0,
    size=500000,
    threads=4,
    prefix="metator",  # "<prefix>_metator" for --prefix <prefix>
)
```

| Command | Function |
|---|---|
| `metator network` | `metator.network.create_contig_data` and `contig_data_from_map`, then `metator.network.metacool_to_contacts` |
| `metator partition` | `metator.partition.partition` |
| `metator validation` | `metator.validation.recursive_decontamination` |
| `metator mge` | `metator.mge.mge_binning` |

The commands also check their inputs, write their log file and clean up
their temporary files, which the functions leave to the caller.

## Checking other clusters

`metator.validation.micomplete_quality` runs miComplete on a directory of
FASTA files, one per cluster, as the validation does, and
`metator.validation.bin_quality` gives the quality class of each from its
results: see [Curating MAGs with anvi'o](curation.md#5-the-curated-clusters),
which checks clusters curated by hand this way.

## Logging

metaTOR logs through the `metator_logger` logger of the `logging` module, at
level `INFO`, to the standard error. `metator.log.set_file_handler` adds a log
file, or replaces the one it added before, as each command does:

```python
import metator.log as mtl

mtl.set_file_handler("my_run.log")
```
