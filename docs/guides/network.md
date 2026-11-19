# Aligning reads and building the network

The first step of metaTOR, `metator network` (run by `metator pipeline`),
aligns the metaHiC reads on the assembly, pairs the mates, and builds the
network of contigs from the pairs. This page covers its choices: the aligner,
the mapping mode, and the normalisation of the network.

## Aligning the reads

`--aligner` (`-b`) chooses the aligner. **We recommend `--aligner bwa`**:

| `--aligner` | How the reads are aligned | Writes |
|---|---|---|
| `bwa` (recommended) | `bwa mem -5SP` on both mates at once: of a read cut in two by a ligation junction, the 5'-most part is the primary alignment, as in the 4DN, Juicer and Arima pipelines | `alignment_<i>.bam` |
| `bowtie2` (default) | Each mate on its own, with `bowtie2 --very-sensitive-local` | `alignment_<i>_for.bam`, `alignment_<i>_rev.bam` |

A pair is kept when both mates align with a mapping quality of at least
`--min-quality` (`-q`, 30 by default). The alignments are written to the
output directory, sorted by read name; the pairs go to
`alignment_<i>_sorted.pairs.gz`, sorted, compressed and indexed with pairix.

`--aligner-mode` (`-B`) chooses how the reads are aligned, with hicstuff's
functions:

- **`normal`** (default): as above.
- **`iterative`**: each mate is truncated to its first 20 bases and aligned;
  those that do not align with a mapping quality of `--min-quality` are
  extended by 20 bases and aligned again, up to their full length. Each mate
  is aligned on its own, with either aligner.
- **`cutsite`**: each read is cut at the ligation sites of `--enzyme`, the
  pieces of the two mates are paired again, and the pairs are aligned in
  `normal` mode. Without `--enzyme`, the reads are aligned in `iterative`
  mode instead.

Cutting the reads at their ligation sites recovers the reads that span a
junction, which neither part of aligns well: it is worth it with long reads
(2×100 bp or more) and an enzyme that cuts often. hicstuff can cut the reads
beforehand, for any aligner:

```bash
hicstuff cutsite -F R1.fq.gz -R R2.fq.gz -p cut/sample -e DpnII -t 16
```

## Contig data and hits

Before building the network, metaTOR measures every contig of the assembly,
in the **contig data** table (`contig_data_network.txt`): its size, its GC
content, its **hits** (the ends of the pairs that lie on it, in all libraries:
a pair within the contig counts twice), its number of restriction sites plus
one, that is of restriction fragments (with `--enzyme`, `-` otherwise), and its shotgun
coverage (with `--depth`, `-` otherwise). From a map, the restriction sites
and the coverage come from its `Restriction_site` and `Shotgun_coverage`
columns where the options give none (see
[Starting from a contact map](metacool.md#starting-from-a-contact-map)). With
several libraries, the hits of each contig in each library go to
`hit_data_alignment.txt`.

## Building the network

Each pair with its mates on two different contigs adds a contact between
them. `--edge` (`-E`, 0 by default) leaves out the pairs with a mate within
that many base pairs of the end of its contig (in `metator network` only:
`pipeline` refuses any other value than 0, as it counts every pair, as its
contact map does); `--self-contacts` (`-s`) adds
the contacts of each contig with itself, which `--edge` does not concern.

`network.txt` holds the network, one line per pair of contigs in contact: the
IDs of the two contigs (as in the contig data table, from 1) and the
normalised weight of the edge. It normalises the contacts of all the
libraries, summed; each library also gets the network of its own contacts,
`network_<i>.txt`.

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

From a map, its `Shotgun_coverage` and `Restriction_site` columns stand for
`--depth` and `--enzyme`.

`empirical_hit` needs no model of the library and works best in our
benchmarks.

### Your own normalisation

To weigh the edges otherwise (with a model such as
[HiCzin](https://github.com/dyxstat/HiCzin)), build the raw network, write
your weights in its third column, and partition it:

```bash
metator network -1 R1.fq.gz -2 R2.fq.gz -a assembly.fa -b bwa -n None -o binning
# ... write binning/my_network.txt from binning/network.txt ...
metator partition -a assembly.fa -c binning/contig_data_network.txt \
    -n binning/my_network.txt -o binning
metator validation -a assembly.fa -c binning/contig_data_partition.txt \
    -n binning/my_network.txt -f binning/overlapping_bin -o binning
```

## From a contact map

`--start metacool` builds the network from a `.metacool` contact map of the
contigs instead of from reads: see
[Starting from a contact map](metacool.md#starting-from-a-contact-map).
