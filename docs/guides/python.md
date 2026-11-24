# Using the Python API

Each command of metaTOR runs functions of its library, which you can call
from Python to run a step, read the outputs, or check other clusters. The
[Python API reference](../api/network.md) lists them all.

## The modules

| Module | Usual name | What it holds |
|---|---|---|
| `metator.inputs` | `mti` | The inputs of a run: the map, or the pairs or BAM files it makes the map from |
| `metator.network` | `mtn` | The contig data table, and the network built from a contact map |
| `metator.partition` | `mtp` | Partitioning the network into core and overlapping clusters |
| `metator.validation` | `mtv` | miComplete, and the recursive partition of contaminated clusters |
| `metator.figures` | `mtf` | The quality plots of the final clusters |
| `metator.metacool` | `mtmc` | Writing the `.metacool` contact map of the clusters |
| `metator.mge` | `mtm` | Clustering MGE contigs |
| `metator.host` | `mth` | Finding the MAGs that host the MGE clusters |
| `metator.io` | `mio` | Reading and writing tables and FASTA files |
| `metator.log` | `mtl` | The logger, `metator.log.logger` |

## Reading the outputs

pandas reads the tables that `metator dump MAP -o tables` writes from the map
of the clusters ([File formats](../formats.md) describes their columns). The
validation's own `contig_data_final.txt`, in its temporary directory (`-N`),
starts with a row number and the `ID` column: read it with `index_col=0` too.

```python
import pandas as pd

summary = pd.read_csv("tables/cluster_summary.txt", sep="\t", index_col=0)
contigs = pd.read_csv("tables/contig_data_final.txt", sep="\t", dtype={"Name": str})
contigs = contigs.set_index("Name", drop=False)

# The high-quality MAGs, and their contigs.
hq = summary[(summary["Weighted completeness"] >= 0.9) & (summary["Weighted redundancy"] <= 1.05)]
hq_contigs = contigs[contigs["Final_cluster"].isin(hq.index)]
```

## The contact map of the clusters

Read the `.metacool` file with
[metacooler](https://github.com/js2264/metacooler), or as any cool file with
cooler (see [The contact map of the clusters](metacool.md)):

```python
import metacooler

# The map of the clusters that metator pipeline sample.metacool -o binning wrote
mc = metacooler.MetaCooler("binning/sample.metacool")
mc.genomes  # one row per cluster, and __unbinned__, with the columns of cluster_summary.txt
mc.contigs  # one row per contig
```

## Running a step

The functions behind the commands take their options as arguments. Here,
the partition of the network that `metator network` stored in a map, as
`metator partition` runs it:

```python
import metator.network as mtn
import metator.partition as mtp

clusters = mtp.partition(
    algorithm="louvain",
    # The network of the map: one edge per pair of contigs in contact, by
    # their positions in the map, weighed by the stored normalization.
    network=mtn.map_network("binning/sample.metacool"),
    iterations=100,
    overlapping_parameter=0.8,  # --overlap 80
    resolution_parameter=1.0,
    threads=4,
    prefix="metator",  # "<prefix>_metator" for --prefix <prefix>
)
# The core and the overlapping clusters, by number, each the IDs of its
# contigs, and the names of the overlapping clusters.
clusters.core, clusters.overlapping, clusters.names
```

| Command | Function |
|---|---|
| `metator network` | `metator.network.create_contig_data` and `contig_data_from_map`, `contig_contacts` and `contig_data_table`, then `metator.metacool.write_map` (`network=`) |
| `metator partition` | `metator.partition.partition`, `metator.clusters.partition_rounds`, then `metator.metacool.write_map` (`partition=`, `rounds=`) |
| `metator validation` | `metator.metacool.read_partition`, `metator.validation.write_overlapping_fasta` and `recursive_decontamination`, then `metator.metacool.write_map` (`clusters=`, a `FinalClusters`, and `network=`) |
| `metator mge` | `metator.mge.mge_binning` and `write_mges_fasta`, then `metator.metacool.write_map` (`mge=`) |
| `metator host` | `metator.network.stored_network`, `metator.host.find_hosts`, `write_hosts` and `plot_hosts`, then `metator.metacool.write_map` (`hosts=`) |

The commands also check their inputs, write their log and clean up their
temporary files; the functions leave that to you.

## Checking other clusters

`metator.validation.micomplete_quality` runs miComplete on a directory of
FASTA files, one per cluster, as the validation does.
`metator.validation.bin_quality` gives the quality class of each from its
results, and `metator.figures` draws the validation's plots. See
[Curating MAGs with anvi'o](curation.md#5-the-curated-clusters) for an
example.

## Logging

metaTOR logs through the `metator_logger` logger of `logging`, at level
`INFO`, to standard error. `metator.log.set_file_handler` adds a log file, or
replaces the one it added before:

```python
import metator.log as mtl

mtl.set_file_handler("my_run.log")
```
