# Building the network

The first step of metaTOR, `metator network` (run by `metator pipeline`),
builds the network of contigs from a `.metacool` contact map of the contigs,
its first argument, which it only reads. This page covers the contig data, the
network and its normalisation.

## The contact map

The map is made before any clustering: hicstuff makes it from the reads of a
library (`hicstuff pipeline -M metacool`), and metacooler from pairs files
(`metacooler create pairs`), or from the maps of several libraries, which
`metacooler merge` sums (see
[Other tools](../other-tools.md#from-the-reads-to-the-map)). The aligner, the
mapping mode and the minimum mapping quality are those of hicstuff
(`--aligner`, `--mapping` and `--quality-min` of `hicstuff pipeline`), which
metaTOR 1.x ran with its own options. Cutting the reads at their ligation
sites recovers the reads that span a junction, which neither part of aligns
well: it is worth it with long reads (2×100 bp or more) and an enzyme that
cuts often (`hicstuff pipeline -m cutsite`).

## Contig data and hits

Before building the network, metaTOR measures every contig of the assembly,
in the **contig data** table (`contig_data_network.txt`): its size, its GC
content, its **hits** (the ends of the contacts of the map that lie on it: a
contact within the contig counts twice), its number of restriction sites plus
one, that is of restriction fragments (with `--enzyme`, `-` otherwise), and its shotgun
coverage (with `--depth`, `-` otherwise). The restriction sites and the
coverage come from the `Restriction_site` and `Shotgun_coverage` columns of
the map where the options give none (see
[Starting from a contact map](metacool.md#starting-from-a-contact-map)).

## Building the network

The contacts of the map between two different contigs, summed over their
bins, are the edge between them: every pair counts (metaTOR 2.0 removed
`--edge`, which left out the pairs with a mate near the end of its contig).
`--self-contacts` (`-s`) adds the contacts of each contig with itself.

`network.txt` holds the network, one line per pair of contigs in contact: the
IDs of the two contigs (as in the contig data table, from 1) and the
normalised weight of the edge. The contacts of several libraries are those
of their map, summed.

## Normalising the network

Contigs differ in size, abundance and number of restriction sites, and the
number of contacts between two contigs grows with all three. `--normalization`
(`-n` in `network`, `-m` in `pipeline`) divides the contacts between two
contigs by a factor of the two:

| `--normalization` | The contacts of contigs A and B are divided by | Needs |
|---|---|---|
| `empirical_hit` (default) | √(hits<sub>A</sub> × hits<sub>B</sub>) | |
| `None` | Nothing: the raw contacts | |
| `length` | √(hits<sub>A</sub>/length<sub>A</sub> × hits<sub>B</sub>/length<sub>B</sub>) | |
| `abundance` | √(coverage<sub>A</sub> × coverage<sub>B</sub>) | `--depth` |
| `RS` | 0.01 × √(hits<sub>A</sub>/RS<sub>A</sub> × hits<sub>B</sub>/RS<sub>B</sub>), RS being the restriction fragments (sites + 1) | `--enzyme` |
| `theoritical_hit` | The geometric mean of the hits expected from coverage, length and restriction sites | `--depth`, `--enzyme` |

The `Shotgun_coverage` and `Restriction_site` columns of the map stand for
`--depth` and `--enzyme`.

`empirical_hit` needs no model of the library and works best in our
benchmarks.

### Your own normalisation

To weigh the edges otherwise (with a model such as
[HiCzin](https://github.com/dyxstat/HiCzin)), build the raw network, write
your weights in its third column, and partition it:

```bash
metator network sample.metacool -a assembly.fa -n None -o binning
# ... write binning/my_network.txt from binning/network.txt ...
metator partition -a assembly.fa -c binning/contig_data_network.txt \
    -n binning/my_network.txt -o binning
metator validation -a assembly.fa -c binning/contig_data_partition.txt \
    -n binning/my_network.txt -f binning/overlapping_clusters -o binning
```
