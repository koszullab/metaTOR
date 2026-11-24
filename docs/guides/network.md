# Building the network

`metator network`, the first step of the pipeline, builds the network of
contigs from a `.metacool` contact map and stores it in the map, where
`metator partition` and `metator validation` read it.

## The contact map

The map comes from hicstuff (`hicstuff pipeline -M metacool`) or metacooler
(see [Other tools](../other-tools.md#from-the-reads-to-the-map)).
`metator network` can also make it from pairs or BAM files, in its output
directory (see [The inputs of a run](../cli.md#inputs)):

```bash
# binning/assembly.metacool, from the BAM files of two libraries, with its network
metator network L1.bam L2.bam -a assembly.fa -e DpnII -o binning
```

With hicstuff, the alignment options are its own (`--aligner`, `--mapping`,
`--quality-min`). With long reads (2×100 bp or more) and an enzyme that cuts
often, `hicstuff pipeline -m cutsite` recovers the reads that span a
ligation junction.

## Contig data and hits

`metator network` first measures every contig of the map, the **contig
data**, and stores it in the contigs table of the map:

| Column | What it is |
|---|---|
| `GC_content` | Its GC content |
| `Hit` | Its **hits**: the ends of contacts that lie on it (a contact within the contig counts twice) |
| `Restriction_site` | Its restriction fragments (sites + 1), with `--enzyme` |
| `Shotgun_coverage` | Its shotgun coverage, with `--depth` |

Without `--enzyme` or `--depth`, the `Restriction_site` and
`Shotgun_coverage` columns already in the map are used (see
[Starting from a contact map](metacool.md#starting-from-a-contact-map)).
Contigs of the assembly missing from the map have no contacts and appear in
no table.

## Building the network

The edge between two contigs is the sum of their contacts in the map. There
is no edge from a contig to itself.

`metator network` stores the network in the map: its contacts, and the
normalisation weights (below) as the bin weights `metator_<normalization>`.
It writes no table. Its nodes are the contigs, by their position in the map
(the `ID` column of `contig_data_final.txt`, from 1). `validation` and
`pipeline` reorder the map, so join tables from different maps by contig
name, not by ID.

To get the network as a table, `metator dump` writes `network.txt` from the
map of the clusters, and metacooler exports it from any map that holds it
(see [The network](../formats.md#the-network)):

```bash
metator dump binning/sample.metacool -o tables
metacooler network binning/sample.metacool --level contig --balance metator_empirical_hit --header \
    | awk -F '\t' '$1 != $2 && $3 != 0' > network.tsv
```

## Normalising the network

The contacts between two contigs grow with their size, abundance and number
of restriction sites. `--normalization` (`-n` in `network`, `-m` in
`pipeline`) divides them by a factor of the two contigs: the product of one
weight per contig, which the map stores, as cooler applies balancing
weights. A contig with no hits (or no coverage, when the normalisation uses
it) gets no weight, and so no edges.

| `--normalization` | The contacts of contigs A and B are divided by | Needs |
|---|---|---|
| `empirical_hit` (default) | √(hits<sub>A</sub> × hits<sub>B</sub>) | |
| `None` | Nothing: the raw contacts | |
| `length` | √(hits<sub>A</sub>/length<sub>A</sub> × hits<sub>B</sub>/length<sub>B</sub>) | |
| `abundance` | √(coverage<sub>A</sub> × coverage<sub>B</sub>) | `--depth` |
| `RS` | 0.01 × √(hits<sub>A</sub>/RS<sub>A</sub> × hits<sub>B</sub>/RS<sub>B</sub>), RS being the restriction fragments | `--enzyme` |
| `theoritical_hit` | The geometric mean of the hits expected from coverage, length and restriction sites | `--depth`, `--enzyme` |

The `Shotgun_coverage` and `Restriction_site` columns of the map can stand
for `--depth` and `--enzyme`. `empirical_hit` needs no model of the library
and works best in our benchmarks.
