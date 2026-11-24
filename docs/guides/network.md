# Building the network

The first step of metaTOR, `metator network` (run by `metator pipeline`),
builds the network of contigs from a `.metacool` contact map of the contigs,
its first argument, and stores it in the map, where `metator partition` and
`metator validation` read it. This page covers the contig data, the network
and its normalisation.

## The contact map

The map is made before any clustering: hicstuff makes it from the reads of a
library (`hicstuff pipeline -M metacool`), and metacooler from pairs files
(`metacooler create pairs`), or from the maps of several libraries, which
`metacooler merge` sums (see
[Other tools](../other-tools.md#from-the-reads-to-the-map)). `metator network`
and `metator pipeline` make it too, from the pairs files or the BAM files of
the run, their arguments, one bin per contig, several libraries summed, in
their output directory (see [The inputs of a run](../cli.md#inputs)):

```bash
# binning/assembly.metacool, made from the BAM files of both mates of two libraries, with its network
metator network L1.bam L2.bam -a assembly.fa -e DpnII -o binning
```

The aligner, the mapping mode and the minimum mapping quality are those of
hicstuff (`--aligner`, `--mapping` and `--quality-min` of `hicstuff
pipeline`), which metaTOR 1.x ran with its own options; from BAM files,
`--min-quality` (`-q`, 30 by default) is the mapping quality that both mates
of a read must reach for its pair to count. Cutting the reads at their ligation
sites recovers the reads that span a junction, which neither part of aligns
well: it is worth it with long reads (2×100 bp or more) and an enzyme that
cuts often (`hicstuff pipeline -m cutsite`).

## Contig data and hits

Before building the network, metaTOR measures every contig of the map, the
**contig data**: its size, its GC content, its **hits** (the ends of the
contacts of the map that lie on it: a contact within the contig counts
twice), its number of restriction sites plus one, that is of restriction
fragments (with `--enzyme`), and its shotgun coverage (with `--depth`). The
restriction sites and the coverage come from the `Restriction_site` and
`Shotgun_coverage` columns of the map where the options give none (see
[Starting from a contact map](metacool.md#starting-from-a-contact-map)).
`metator network` stores them in the contigs table of the map: `GC_content`,
`Hit`, and `Shotgun_coverage` and `Restriction_site` when it has them. The
contigs of the assembly that the map lacks (a map made with a contig filter)
have no contact: they are in no table of the run.

## Building the network

The contacts of the map between two different contigs, summed over their
bins, are the edge between them: every pair counts (metaTOR 2.0 removed
`--edge`, which left out the pairs with a mate near the end of its contig).
The contacts of a contig with itself are no edge (metaTOR 2.0 removed
`--self-contacts`, which added them).

The network is the map itself: its contacts, summed by pair of contigs, and
the weights of its normalisation (below), which `metator network` stores in
the map as the weights `metator_<normalization>` of its bins. The contacts of
several libraries are those of their map, summed. The nodes of the network are
the contigs of the map, by their position in it (the `ID` column of
`contig_data_final.txt`, from 1): a run that regroups the map
(`validation`, `pipeline`) moves them, so join a table to the network of
another map by contig name, not by ID. `metator network` writes no table of
it. `metator dump` writes it from the map of the clusters of a run, as
`network.txt`; metacooler exports it from any map that holds it, by contig
names, once the contacts within each contig, which metacooler writes too and
are no edge, are dropped (see [The network](../formats.md#the-network)):

```bash
# From the map of the clusters of a run (metator pipeline sample.metacool -o binning):
# network.txt, among its tables
metator dump binning/sample.metacool -o tables
# From any map that holds it: node1, node2, value (the raw counts without --balance)
metacooler network binning/sample.metacool --level contig --balance metator_empirical_hit --header \
    | awk -F '\t' '$1 != $2 && $3 != 0' > network.tsv
```

## Normalising the network

Contigs differ in size, abundance and number of restriction sites, and the
number of contacts between two contigs grows with all three. `--normalization`
(`-n` in `network`, `-m` in `pipeline`) divides the contacts between two
contigs by a factor of the two, the product of one weight per contig: the map
stores the weights, and an edge weighs its contacts times the weights of its
two contigs, as cooler applies balancing weights. A contig whose weight is
infinite (no hit, or no coverage) has none, and its contacts are no edge.

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
