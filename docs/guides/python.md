# Using the Python API

Each command of metaTOR runs functions of its library, which can be called
from Python to run a step, read the outputs, or redo a plot on other bins.
This page shows the pieces to put together; the
[Python API reference](../api/network.md) lists them all.

## The modules

| Module | Usual name | What it holds |
|---|---|---|
| `metator.align` | `mta` | Aligning the reads and pairing the mates |
| `metator.network` | `mtn` | The contig data table, and the network built from the pairs or a contact map |
| `metator.partition` | `mtp` | Partitioning the network into core and overlapping bins |
| `metator.validation` | `mtv` | miComplete, and the recursive partition of contaminated bins |
| `metator.metacool` | `mtmc` | Writing the `.metacool` contact map of the bins |
| `metator.contact_map` | `mtc` | The contact map of one contig or bin |
| `metator.scaffold`, `metator.regions` | `mts` | Scaffolding a bin |
| `metator.quality_check` | `mtq` | The quality of a library |
| `metator.mge`, `metator.host` | `mtm` | Binning MGE contigs, and finding their hosts |
| `metator.figures` | `mtf` | The plots |
| `metator.io` | `mio` | Reading and writing tables, pairs and FASTA files |
| `metator.log` | `mtl` | The logger, `metator.log.logger` |

## Reading the outputs

`read_bin_summary` and `read_contig_data` read the two main tables as pandas
data frames, indexed by bin name and by contig name:

```python
import metator.io as mio

summary = mio.read_bin_summary("binning/bin_summary.txt")
contigs = mio.read_contig_data("binning/contig_data_final.txt")

# The high-quality MAGs, and their contigs.
hq = summary[(summary["Weighted completeness"] >= 0.9) & (summary["Weighted redundancy"] <= 1.05)]
hq_contigs = contigs[contigs["Final_bin"].isin(hq.index)]
```

[File formats](../formats.md) describes their columns.

## The contact map of the bins

The `.metacool` file is read with
[metacooler](https://github.com/js2264/metacooler), or as any cool file with
cooler (see [The contact map of the bins](metacool.md)):

```python
import metacooler

mc = metacooler.MetaCooler("binning/metator.metacool")
mc.genomes  # one row per bin, with the columns of bin_summary.txt
mc.contigs  # one row per contig
```

## The clustering matrices

With `--cluster-matrix`, the partition writes how often each two contigs fell
in the same community, as a sparse matrix in a `.npz` file. Row and column
*i* are the contig of ID *i* of the contig data table (from 1: row and column
0 are empty), and only the upper triangle is stored:

```python
import metator.io as mio
from scipy import sparse

matrix = sparse.load_npz("binning/clustering_matrix_partition.npz").tocsr()
matrix[1, 1196]  # how often contigs 1 and 1196 were binned together, from 0 to 1

# As a text file: a header line (rows, columns, non-zero values), then row, column, value
mio.save_sparse_matrix(matrix.tocoo(), "binning/clustering_matrix_partition.tsv")
```

## Running a step

The functions behind the commands take their options as arguments. The
partition of a network, as `metator partition` runs it:

```python
import os

import metator.partition as mtp

os.makedirs("binning/overlapping_bin", exist_ok=True)
clustering_matrix_file, contig_data_file = mtp.partition(
    algorithm="louvain",
    assembly="assembly.fa",
    cluster_matrix=True,
    contig_data_file="binning/contig_data_network.txt",
    iterations=100,
    network_file="binning/network.txt",
    outdir="binning",
    fasta_dir="binning/overlapping_bin",
    overlapping_parameter=0.8,  # --overlap 80
    resolution_parameter=1.0,
    size=500000,
    threads=4,
    prefix="metator",  # "<prefix>_metator" for --prefix <prefix>
)
```

| Command | Function |
|---|---|
| `metator network` | `metator.align.get_contact_pairs`, then `metator.network.alignment_to_contacts` (or `metacool_to_contacts`, from a contact map) |
| `metator partition` | `metator.partition.partition` |
| `metator validation` | `metator.validation.recursive_decontamination` |
| `metator metacool` | `metator.metacool.write_metacool` |
| `metator contactmap` | `metator.contact_map.generate_contact_map` |
| `metator scaffold` | `metator.scaffold.get_scaffolds` |
| `metator qc` | `metator.quality_check.quality_check` |
| `metator mge` | `metator.mge.mge_binning` |
| `metator pairs` | `metator.io.sort_and_index_pairs` |

The commands also check their inputs, write their log file and clean up
their temporary files, which the functions leave to the caller.

## Checking other bins

`metator.validation.micomplete_quality` runs miComplete on a directory of
FASTA files, one per bin, as the validation does, and the functions of
`metator.figures` draw the plots of the validation from its results: see
[Curating MAGs with anvi'o](curation.md#5-the-curated-bins), which checks bins
curated by hand this way.

## Logging

metaTOR logs through the `metator_logger` logger of the `logging` module, at
level `INFO`, to the standard error. `metator.log.set_file_handler` adds a log
file, or replaces the one it added before, as each command does:

```python
import metator.log as mtl

mtl.set_file_handler("my_run.log")
```
