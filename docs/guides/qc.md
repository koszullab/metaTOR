# Checking a library

`metator qc` measures the quality of a metaHiC library on the genomes that
the binning recovered best. Within the contigs of good MAGs, each pair can be
told apart as a real contact or an artefact of the protocol, and a pair that
joins two of these MAGs can only be noise.

## What it measures

1. **The MAGs.** It takes the bins of `bin_summary.txt` with a weighted
   completeness above 70% and a weighted redundancy under 1.15, and their
   contigs longer than 100 kb.
2. **The pairs.** It reads, from the pairs files, the pairs between these
   contigs, and finds the restriction fragment of each mate, from the
   restriction sites of `--enzyme`.
3. **The events.** Each pair within a contig is one of hicstuff's
   [3C events](https://koszullab.github.io/hicstuff/guides/filter/):
   - an **uncut** (*religated* in the log): the two mates face each other
     (`+-`) across a few restriction sites, a fragment that was not cut, or
     was cut and ligated again;
   - a **loop**: the mates face away from each other (`-+`) across a few
     restriction sites, the two ends of a fragment ligated together;
   - a **weird**: the mates are on the same fragment and the same strand;
   - **informative** otherwise.

   A pair between two contigs of the same MAG is informative, and one between
   two MAGs is noise.

"A few" is a threshold for the uncuts and one for the loops, in restriction
sites: hicstuff sets them from the library (as `hicstuff filter` does), or
`--threshold` (`-t`) gives them, as `uncuts,loops`: `-t 3,3`.

## Running it

```bash
metator qc -a assembly.fa -e DpnII -O binning -p sample -P binning/alignment_0_sorted.pairs.gz
```

It reads `bin_summary.txt` and `contig_data_final.txt` in `--metator-dir`
(`-O`, the current directory by default), or where `--bin-summary` (`-b`) and
`--contig-data` (`-c`) say, and the pairs files of the run (kept with
`metator pipeline --no-clean-up`).

| Option | Default | Effect |
|---|---|---|
| `--enzyme` (`-e`) | *required* | The enzymes of the library, separated by commas |
| `--threshold` (`-t`) | set from the library | The uncut and loop thresholds |
| `--prefix` (`-p`) | none | The name of the sample, on the plots and their files |
| `--plot` (`-P`) | off | Write the plots |
| `--outdir` (`-o`) | `metator_qc`, or `metator_qc_<prefix>` | Where the plots and the log go |
| `--no-clean-up` (`-n`) | off | Keep the temporary files: the pairs extracted, with their fragments |

With `--plot`, it writes `<prefix>_event.pdf`, the distribution of the events
from which the thresholds are set (not drawn when `--threshold` gives them),
and `<prefix>_camembert_plot.pdf`, the share of each kind of pair.

## The results

The results are in the log, `metator_qc_<date>.log`. On a simulated
metagenome of two genomes, in contigs larger than 100 kb:

```console
INFO :: 30000 pairs have been extracted.
INFO :: Filtering with thresholds: uncuts=2 loops=2
INFO :: Number of MAGs: 2
INFO :: Contacts within contigs: 27188.
INFO :: Contacts within MAGs: 29736.
INFO :: Contacts between MAGs: 264.
INFO :: Religated ratio: 4.77%.
INFO :: Loop ratio: 4.73%.
INFO :: Weirds ratio: 3.77%.
INFO :: Informative contacts estimation: 85.97%.
INFO :: Ratio inter/intra contigs: 9.88%.
INFO :: Noise contact ratio: 0.88%
INFO :: Noise score: 1.74E-02
```

| Line | What it is |
|---|---|
| `Number of MAGs` | The MAGs selected |
| `Contacts within contigs`, `Contacts within MAGs`, `Contacts between MAGs` | The pairs extracted within a contig; within the MAGs (within a contig, or between two contigs of the same MAG); and between two MAGs |
| `Religated ratio`, `Loop ratio`, `Weirds ratio` | The uncuts, loops and weirds, in percent of the pairs within the MAGs (within a contig, or between two contigs of the same MAG) |
| `Informative contacts estimation` | The informative pairs, in percent of all the pairs extracted |
| `Ratio inter/intra contigs` | The share of the informative pairs that join two contigs |
| `Noise contact ratio` | The pairs between two MAGs, in percent of all the pairs extracted |
| `Noise score` | The contacts per pair of MAGs, against those per MAG: *x* / (*y* + *x*), with *x* the pairs between MAGs divided by the number of pairs of MAGs, and *y* the pairs within MAGs divided by the number of MAGs. 0 for a library with no noise |

When no MAG passes, or none has a contig longer than 100 kb, as on the small
test library of the [Quickstart](../quickstart.md), no pair is extracted, and
the log says why, such as `The 2 high-quality MAGs have no contig longer than
100 kb: no quality check.`

## From a contact map

`qc` also takes the `.metacool` contact map of the run (or an entry of a
`.mmetacool` collection) in place of its pairs files, and measures from it
the noise of the library, and only it:

```bash
metator qc -a assembly.fa -e DpnII -O binning -p sample binning/metator.metacool
```

The map holds the number of pairs between each two of its bins: those
between the contigs of the MAGs, summed by pair of contigs, are the pairs
that `qc` extracts from the pairs files the map was made of (each pair whose
two reads lie on the contigs, within their length), when their pairs are in
the upper triangle, as metaTOR's alignment writes them. So the contacts
within contigs, within MAGs and between MAGs, the noise contact ratio and
the noise score are those of these pairs files, line for line. From pairs
files with the reads of some pairs in the other order, `qc` reads the pairs
of two contigs in one order only and leaves out the others, which the map
counts. On the simulated metagenome above:

```console
INFO :: 30000 pairs have been extracted from binning/metator.metacool.
INFO :: Number of MAGs: 2
INFO :: Contacts within contigs: 27188.
INFO :: Contacts within MAGs: 29736.
INFO :: Contacts between MAGs: 264.
INFO :: Noise contact ratio: 0.88%
INFO :: Noise score: 1.74E-02
INFO :: The religated, loop and weirds ratios, the informative contacts estimation and the ratio inter/intra contigs are not measured from a contact map: they need the positions and strands of the reads, which only pairs files hold.
```

The events of the pairs within a contig (uncuts, loops, weirds and
informative pairs) need the positions and strands of their reads, which a
map does not hold: they are not measured, and `--threshold` and `--plot`,
which set and draw them, are refused. `--assembly` and `--enzyme` are not
used. Refused too, with the reason in the log: a map with pairs files, or
with another map; a map that metacooler cannot open (a `.mmetacool`
collection rather than one of its entries, for example); a map that gives a
contig of the MAGs another length than the contig data table (a map of
another assembly); and a map that stores both triangles of its matrix. A
contig of the MAGs that the map lacks (a map made with a contig filter) has
no contacts there, with a warning that names it.
